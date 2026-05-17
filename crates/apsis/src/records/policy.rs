//! Snapshot cadence policies for the record writer.

/// Snapshot cadence — controls how often `RecordHook` writes a `Snapshot`
/// frame between the mandatory initial and final bookends.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum RecordPolicy {
    /// Initial + final bookend snapshots; events always written.
    #[default]
    BookendsAndEvents,
    /// Bookends + a `Snapshot` whenever `steps % N == 0`.
    EveryNSteps(u32),
    /// Bookends + a `Snapshot` whenever sim time crosses a multiple of `dt`.
    /// Caller is expected to pass a strictly positive `dt`; values `<= 0`
    /// collapse to "fire every step", silently equivalent to [`Self::Dense`].
    EveryTime(f64),
    /// A `Snapshot` every step. Debug mode.
    Dense,
}

impl RecordPolicy {
    /// Decide whether the writer should emit a Snapshot for this post-step
    /// fire. `t_last_snapshot` tracks the last time a Snapshot was emitted
    /// (None before the initial bookend).
    pub fn should_snapshot(&self, t: f64, steps: u64, t_last_snapshot: Option<f64>) -> bool {
        match *self {
            Self::BookendsAndEvents => false,
            Self::EveryNSteps(n) => steps > 0 && steps.is_multiple_of(n as u64),
            Self::EveryTime(dt) => match t_last_snapshot {
                None => true,
                Some(t_last) => t >= t_last + dt,
            },
            Self::Dense => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_bookends_and_events() {
        assert_eq!(RecordPolicy::default(), RecordPolicy::BookendsAndEvents);
    }

    #[test]
    fn bookends_only_never_snapshots_between() {
        let p = RecordPolicy::BookendsAndEvents;
        assert!(!p.should_snapshot(0.0, 0, None));
        assert!(!p.should_snapshot(1.5, 100, Some(0.0)));
    }

    #[test]
    fn every_n_steps_fires_on_multiples() {
        let p = RecordPolicy::EveryNSteps(100);
        assert!(!p.should_snapshot(0.0, 0, None));
        assert!(!p.should_snapshot(0.1, 50, None));
        assert!(p.should_snapshot(0.2, 100, None));
        assert!(!p.should_snapshot(0.21, 101, Some(0.2)));
        assert!(p.should_snapshot(0.4, 200, Some(0.2)));
    }

    #[test]
    fn every_time_fires_when_interval_crossed() {
        let p = RecordPolicy::EveryTime(0.1);
        assert!(p.should_snapshot(0.05, 1, None));
        assert!(!p.should_snapshot(0.05, 2, Some(0.05)));
        assert!(p.should_snapshot(0.16, 3, Some(0.05)));
    }

    #[test]
    fn dense_fires_every_step() {
        let p = RecordPolicy::Dense;
        for s in 1..10u64 {
            assert!(p.should_snapshot(s as f64 * 0.01, s, None));
        }
    }
}
