//! Orchestration for **Precision Run** mode — the run-to-completion
//! execution discipline required by
//! [`ExecutionProfile::Precision`](crate::physics::integrator::traits::ExecutionProfile::Precision)
//! integrators (today, IAS15).
//!
//! # Architectural role
//!
//! The [`PrecisionRunController`] is the single source of truth for
//! "a precision run is in progress". Both the physics thread and the
//! UI read its state:
//!
//! * **Physics thread** consults the controller between sub-steps to
//!   decide whether to keep integrating, pause cleanly, or abort.
//!   Pausing observes a substep boundary — never mid-iteration —
//!   so determinism is preserved.
//! * **UI** reads `state()` and [`telemetry()`](Self::telemetry) to
//!   render the Precision Run panel, progress bar, and throughput
//!   readouts. It issues control intent (`request_pause`,
//!   `request_abort`, etc.) and *observes* the state change on the
//!   next physics tick — it does not mutate the controller
//!   directly beyond those intent methods.
//!
//! # State machine
//!
//! ```text
//!                          ┌──────┐
//!    set_integrator(IAS15) │ Idle │
//!    chooses run → start() └──┬───┘
//!                             │
//!                             ▼
//!                        ┌─────────┐    request_pause()        ┌────────┐
//!                   ┌───→│ Running │ ────────────────────────→ │ Paused │
//!                   │    └───┬─────┘                           └───┬────┘
//!                   │        │                                     │
//!                   │        │  t_sim reaches t_target             │
//!                   │        │  (physics reports via             ┌ resume()
//!                   │        │   mark_completed)                  │
//!                   │        ▼                                    │
//!                   │   ┌───────────┐                             │
//!                   │   │ Completed │  ← user clicks Close → Idle │
//!                   │   └───────────┘                             │
//!                   │                                             │
//!                   │           request_abort()                   │
//!                   └─── Aborting ──── (physics observes) ────────┘
//! ```
//!
//! Intent methods (`request_pause`, `request_abort`) are
//! **non-blocking** — they flip the state to a transient signal and
//! return. The physics thread observes on its next substep boundary
//! and calls the confirming method (`mark_paused`, `mark_aborted`)
//! once the substep completes. This separation keeps determinism:
//! no intermediate state is ever half-committed.

use std::time::{Duration, Instant};

pub mod telemetry;

pub use telemetry::{Telemetry, TelemetryBuilder};

/// State a precision run can be in. Transitions are driven by
/// [`PrecisionRunController`] methods with explicit contracts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RunState {
    /// No run in progress. The simulator is in its normal real-time
    /// or idle state.
    Idle,

    /// A run is actively advancing. `t_target` is the absolute
    /// simulation time the run is trying to reach; `t_start` is the
    /// simulation time at which the run began; `wall_start` is the
    /// wall-clock moment at which the physics thread started
    /// stepping. Progress percent is `(t_sim - t_start) / (t_target - t_start)`.
    Running {
        t_target: f64,
        t_start: f64,
        wall_start: Instant,
    },

    /// The UI has requested a pause; the physics thread will
    /// transition to [`Paused`](RunState::Paused) at the next
    /// substep boundary. Distinct from `Paused` so callers can
    /// render "Pausing…" while the request is in flight.
    Pausing {
        t_target: f64,
        t_start: f64,
        wall_elapsed: Duration,
    },

    /// The run is paused at a clean substep boundary. No physics
    /// work is happening. `wall_elapsed` is the cumulative wall
    /// time spent *running* so far (pause time excluded) — used to
    /// keep throughput calculations honest when resuming.
    Paused {
        t_target: f64,
        t_start: f64,
        wall_elapsed: Duration,
    },

    /// The UI has requested an abort; the physics thread will
    /// transition to [`Completed`](RunState::Completed) with
    /// [`RunOutcome::UserAborted`] at the next substep boundary.
    Aborting { t_target: f64, t_start: f64 },

    /// The run has ended. [`RunOutcome`] classifies why.
    Completed { outcome: RunOutcome },
}

/// Why a precision run ended. Drives the Summary panel's framing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RunOutcome {
    /// The run reached `t_target` under its own power.
    Reached,
    /// The user clicked Abort before completion. Partial state is
    /// still available in the controller / system; the UI decides
    /// whether to offer Export of the partial run.
    UserAborted,
    /// The physics thread observed an unrecoverable error (integrator
    /// panic, force-model failure). Today no producer emits this —
    /// included so the state-machine is exhaustive and future error
    /// paths do not require an enum-variant migration.
    Errored,
}

/// Controller for a precision run. Owns the state machine, observed
/// telemetry, and the intent signals flowing from UI to physics.
pub struct PrecisionRunController {
    state: RunState,
    telemetry: Telemetry,
}

impl PrecisionRunController {
    /// Construct a fresh controller in [`RunState::Idle`]. This is
    /// the only state in which `start` is legal.
    pub fn new() -> Self {
        Self {
            state: RunState::Idle,
            telemetry: Telemetry::default(),
        }
    }

    /// Current state. Cheap to read (copies a small enum); the UI
    /// may call every frame.
    pub fn state(&self) -> RunState {
        self.state
    }

    /// Immutable view of the latest telemetry sample. Updated by
    /// [`update_telemetry`](Self::update_telemetry) from the physics
    /// thread on a regular cadence.
    pub fn telemetry(&self) -> &Telemetry {
        &self.telemetry
    }

    /// Progress as a fraction `[0.0, 1.0]` for the UI progress bar.
    /// Returns `0.0` outside `Running`/`Pausing`/`Paused`.
    pub fn progress(&self, t_sim: f64) -> f32 {
        match self.state {
            RunState::Running { t_target, t_start, .. }
            | RunState::Pausing { t_target, t_start, .. }
            | RunState::Paused { t_target, t_start, .. }
            | RunState::Aborting { t_target, t_start, .. } => {
                let span = (t_target - t_start).max(f64::MIN_POSITIVE);
                let done = (t_sim - t_start).max(0.0);
                ((done / span) as f32).clamp(0.0, 1.0)
            }
            RunState::Idle => 0.0,
            RunState::Completed { outcome } => match outcome {
                RunOutcome::Reached => 1.0,
                RunOutcome::UserAborted | RunOutcome::Errored => {
                    // Leave the bar where telemetry last put it by
                    // returning the last recorded fraction. Callers
                    // that want a different rendering can detect
                    // `RunOutcome::UserAborted` explicitly.
                    self.telemetry.last_progress_fraction
                }
            },
        }
    }

    // ── State transitions ────────────────────────────────────────────────────

    /// Transition `Idle → Running`. Panics in debug if called from
    /// another state; callers should gate on `state()`.
    pub fn start(&mut self, t_target: f64, t_start_sim: f64) {
        debug_assert!(
            matches!(self.state, RunState::Idle | RunState::Completed { .. }),
            "start() requires Idle or Completed state; got {:?}",
            self.state
        );
        self.state = RunState::Running {
            t_target,
            t_start: t_start_sim,
            wall_start: Instant::now(),
        };
        self.telemetry = Telemetry::default();
    }

    /// Flip the state to `Pausing` (intent signal). The physics
    /// thread observes this at its next substep boundary and calls
    /// [`mark_paused`](Self::mark_paused) once the substep lands.
    /// No-op if not currently `Running`.
    pub fn request_pause(&mut self) {
        if let RunState::Running { t_target, t_start, wall_start } = self.state {
            self.state = RunState::Pausing {
                t_target,
                t_start,
                wall_elapsed: wall_start.elapsed(),
            };
        }
    }

    /// Called by the physics thread after the substep-in-flight
    /// completes, in response to a prior `request_pause`. Promotes
    /// `Pausing → Paused`. No-op outside `Pausing`.
    pub fn mark_paused(&mut self) {
        if let RunState::Pausing { t_target, t_start, wall_elapsed } = self.state {
            self.state = RunState::Paused { t_target, t_start, wall_elapsed };
        }
    }

    /// Resume from `Paused → Running`. No-op outside `Paused`.
    /// `wall_elapsed` is preserved across the pause so throughput
    /// remains honest.
    pub fn resume(&mut self) {
        if let RunState::Paused { t_target, t_start, wall_elapsed } = self.state {
            // Reconstruct `wall_start` so `wall_start.elapsed()` again
            // yields the correct cumulative value.
            self.state = RunState::Running {
                t_target,
                t_start,
                wall_start: Instant::now() - wall_elapsed,
            };
        }
    }

    /// Flip the state to `Aborting` (intent signal). The physics
    /// thread observes this at its next substep boundary and calls
    /// [`mark_aborted`](Self::mark_aborted).
    ///
    /// Legal from `Running`, `Pausing`, and `Paused` — the user may
    /// abort without resuming first. No-op from other states.
    pub fn request_abort(&mut self) {
        let (t_target, t_start) = match self.state {
            RunState::Running { t_target, t_start, .. }
            | RunState::Pausing { t_target, t_start, .. }
            | RunState::Paused { t_target, t_start, .. } => (t_target, t_start),
            _ => return,
        };
        self.state = RunState::Aborting { t_target, t_start };
    }

    /// Called by the physics thread after the substep-in-flight
    /// completes, in response to a prior `request_abort`. Promotes
    /// `Aborting → Completed { UserAborted }`.
    pub fn mark_aborted(&mut self) {
        if matches!(self.state, RunState::Aborting { .. }) {
            self.state = RunState::Completed { outcome: RunOutcome::UserAborted };
        }
    }

    /// Called by the physics thread when `t_sim >= t_target` lands
    /// cleanly. Promotes `Running → Completed { Reached }`.
    /// No-op outside `Running`.
    pub fn mark_completed(&mut self) {
        if matches!(self.state, RunState::Running { .. }) {
            self.state = RunState::Completed { outcome: RunOutcome::Reached };
        }
    }

    /// Called by the physics thread when an unrecoverable integrator
    /// or force-model error occurs. Promotes any active state to
    /// `Completed { Errored }`. No-op from `Idle` or `Completed`.
    pub fn mark_errored(&mut self) {
        if !matches!(self.state, RunState::Idle | RunState::Completed { .. }) {
            self.state = RunState::Completed { outcome: RunOutcome::Errored };
        }
    }

    /// Acknowledge a finished run and return to `Idle`. Called by
    /// the UI when the user closes the Summary panel. No-op outside
    /// `Completed`.
    pub fn acknowledge(&mut self) {
        if matches!(self.state, RunState::Completed { .. }) {
            self.state = RunState::Idle;
        }
    }

    // ── Telemetry update ─────────────────────────────────────────────────────

    /// Replace the telemetry snapshot. Called by the physics thread
    /// on a regular cadence (typically every substep or every
    /// wall-second, whichever is less frequent).
    pub fn update_telemetry(&mut self, new_telemetry: Telemetry) {
        self.telemetry = new_telemetry;
    }
}

impl Default for PrecisionRunController {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fresh() -> PrecisionRunController {
        PrecisionRunController::new()
    }

    #[test]
    fn fresh_controller_is_idle() {
        assert_eq!(fresh().state(), RunState::Idle);
    }

    #[test]
    fn start_transitions_to_running() {
        let mut c = fresh();
        c.start(10.0, 0.0);
        assert!(matches!(c.state(), RunState::Running { .. }));
    }

    #[test]
    fn pause_flow_preserves_wall_elapsed() {
        let mut c = fresh();
        c.start(10.0, 0.0);
        // Simulate a bit of wall time before the pause request.
        std::thread::sleep(std::time::Duration::from_millis(3));
        c.request_pause();
        assert!(matches!(c.state(), RunState::Pausing { .. }));
        c.mark_paused();
        match c.state() {
            RunState::Paused { wall_elapsed, .. } => {
                assert!(
                    wall_elapsed >= std::time::Duration::from_millis(3),
                    "wall_elapsed should carry across the pause transition"
                );
            }
            other => panic!("expected Paused, got {:?}", other),
        }
    }

    #[test]
    fn resume_from_paused_restores_running() {
        let mut c = fresh();
        c.start(10.0, 0.0);
        c.request_pause();
        c.mark_paused();
        c.resume();
        assert!(matches!(c.state(), RunState::Running { .. }));
    }

    #[test]
    fn abort_from_running_via_signal_then_mark() {
        let mut c = fresh();
        c.start(10.0, 0.0);
        c.request_abort();
        assert!(matches!(c.state(), RunState::Aborting { .. }));
        c.mark_aborted();
        assert_eq!(
            c.state(),
            RunState::Completed { outcome: RunOutcome::UserAborted }
        );
    }

    #[test]
    fn abort_from_paused_is_allowed() {
        let mut c = fresh();
        c.start(10.0, 0.0);
        c.request_pause();
        c.mark_paused();
        c.request_abort();
        assert!(matches!(c.state(), RunState::Aborting { .. }));
    }

    #[test]
    fn mark_completed_only_from_running() {
        let mut c = fresh();
        c.start(10.0, 0.0);
        c.mark_completed();
        assert_eq!(
            c.state(),
            RunState::Completed { outcome: RunOutcome::Reached }
        );
    }

    #[test]
    fn acknowledge_returns_to_idle() {
        let mut c = fresh();
        c.start(10.0, 0.0);
        c.mark_completed();
        c.acknowledge();
        assert_eq!(c.state(), RunState::Idle);
    }

    #[test]
    fn progress_is_zero_when_idle() {
        let c = fresh();
        assert_eq!(c.progress(5.0), 0.0);
    }

    #[test]
    fn progress_reflects_running_fraction() {
        let mut c = fresh();
        c.start(10.0, 0.0);
        assert_eq!(c.progress(5.0), 0.5);
        assert_eq!(c.progress(0.0), 0.0);
        assert_eq!(c.progress(10.0), 1.0);
        assert_eq!(c.progress(15.0), 1.0); // clamped
    }

    #[test]
    fn telemetry_is_reset_on_start() {
        let mut c = fresh();
        c.update_telemetry(Telemetry {
            substeps: 999,
            ..Default::default()
        });
        c.start(10.0, 0.0);
        assert_eq!(c.telemetry().substeps, 0);
    }
}
