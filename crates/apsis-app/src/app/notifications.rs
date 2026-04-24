//! Consumer-side store for [`apsis::core::log`] events — powers
//! the Notification Center panel in the top bar.
//!
//! # Design
//!
//! Subscribes to the event bus once at app startup. Incoming events
//! are buffered into a bounded ring (fixed capacity) so a long run
//! cannot exhaust memory. A simple time-windowed coalesce groups
//! repeated events with the same [`Event::coalesce_key`] — what the
//! user sees is one entry with a multiplier, not N identical lines.
//!
//! # Coalesce semantics
//!
//! An event arriving with `coalesce_key = Some(k)` is merged with
//! the most recent entry that also has `coalesce_key = Some(k)` and
//! the same [`Level`], *if* that entry's `last_at` is within
//! [`COALESCE_WINDOW`] of now. Otherwise it becomes a new entry.
//! Outside the window, a fresh run of the same key starts a new
//! entry — this keeps the display honest when a pathological burst
//! is separated from a later isolated hit.
//!
//! # Threading
//!
//! Shared via `Arc<Mutex<NotificationStore>>`. The bus invokes the
//! subscriber synchronously from whichever thread published, so the
//! lock must be held briefly. All store methods are O(k) over the
//! most recent coalescing window, never O(N) over the full buffer.

use apsis::core::log::{Event, Level};
use std::collections::VecDeque;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// UI filter selector for the Notification Center. Owned by
/// `SimulationApp`; used by the panel to decide which entries to
/// render.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NotificationFilter {
    #[default]
    All,
    Info,
    Warn,
    Error,
}

impl NotificationFilter {
    pub fn matches(self, level: Level) -> bool {
        match self {
            Self::All => true,
            Self::Info => level == Level::Info,
            Self::Warn => level == Level::Warn,
            Self::Error => level == Level::Error,
        }
    }
}

/// Time window within which same-key events merge.
pub const COALESCE_WINDOW: Duration = Duration::from_millis(500);

/// Ring-buffer capacity — the hard upper bound on how many entries
/// the store remembers.
pub const MAX_ENTRIES: usize = 200;

/// How far back we scan when looking for a coalesce target. Limits
/// the per-publish work regardless of buffer size — bounded tail.
const COALESCE_SCAN_DEPTH: usize = 16;

/// One visible entry in the notification center. May represent a
/// single event or a coalesced run.
#[derive(Debug, Clone)]
pub struct NotificationEntry {
    pub event: Event,
    /// Number of events merged into this entry. `1` for a pristine
    /// entry; bumped on each coalesce hit.
    pub count: u32,
    pub first_at: Instant,
    pub last_at: Instant,
}

impl NotificationEntry {
    fn new(event: Event) -> Self {
        let now = Instant::now();
        Self { event, count: 1, first_at: now, last_at: now }
    }
}

pub struct NotificationStore {
    entries: VecDeque<NotificationEntry>,
    unread: usize,
}

impl NotificationStore {
    pub fn new() -> Self {
        Self { entries: VecDeque::with_capacity(MAX_ENTRIES), unread: 0 }
    }

    pub fn ingest(&mut self, event: Event) {
        if let Some(key) = event.coalesce_key {
            let now = Instant::now();
            let scan_limit = COALESCE_SCAN_DEPTH.min(self.entries.len());
            for entry in self.entries.iter_mut().rev().take(scan_limit) {
                if now.duration_since(entry.last_at) > COALESCE_WINDOW {
                    break;
                }
                if entry.event.coalesce_key == Some(key) && entry.event.level == event.level {
                    entry.count = entry.count.saturating_add(1);
                    entry.last_at = now;
                    self.unread = self.unread.saturating_add(1);
                    return;
                }
            }
        }

        if self.entries.len() >= MAX_ENTRIES {
            self.entries.pop_front();
        }
        self.entries.push_back(NotificationEntry::new(event));
        self.unread = self.unread.saturating_add(1);
    }

    pub fn entries(&self) -> impl DoubleEndedIterator<Item = &NotificationEntry> {
        self.entries.iter()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn unread_count(&self) -> usize {
        self.unread
    }

    pub fn mark_all_read(&mut self) {
        self.unread = 0;
    }

    pub fn clear(&mut self) {
        self.entries.clear();
        self.unread = 0;
    }

    /// Count entries matching a level filter. Used by the UI filter
    /// chips to render counts next to each level.
    pub fn count_at_level(&self, level: Level) -> usize {
        self.entries.iter().filter(|e| e.event.level == level).count()
    }
}

impl Default for NotificationStore {
    fn default() -> Self {
        Self::new()
    }
}

/// Hook the shared store into the event bus. Returns the
/// subscription id — the caller is responsible for keeping the
/// `Arc` alive; dropping it detaches the subscriber on the next
/// event (the bus holds an `Arc` on the callback closure, which
/// in turn holds this `Arc`, so the store survives as long as
/// the bus callback does).
pub fn attach_to_bus(store: Arc<Mutex<NotificationStore>>) -> apsis::core::log::SubscriptionId {
    apsis::core::log::subscribe(move |event: &Event| {
        if let Ok(mut s) = store.lock() {
            s.ingest(event.clone());
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use apsis::core::log::Source;

    fn event(level: Level, message: &'static str, coalesce: Option<&'static str>) -> Event {
        let mut e = Event::new(level, Source::System, message);
        if let Some(k) = coalesce {
            e = e.with_coalesce_key(k);
        }
        e
    }

    #[test]
    fn fresh_store_is_empty() {
        let s = NotificationStore::new();
        assert!(s.is_empty());
        assert_eq!(s.unread_count(), 0);
    }

    #[test]
    fn ingest_appends_without_coalesce_key() {
        let mut s = NotificationStore::new();
        s.ingest(event(Level::Info, "a", None));
        s.ingest(event(Level::Info, "a", None));
        assert_eq!(s.len(), 2);
        assert_eq!(s.unread_count(), 2);
    }

    #[test]
    fn same_coalesce_key_within_window_merges() {
        let mut s = NotificationStore::new();
        s.ingest(event(Level::Warn, "x", Some("k")));
        s.ingest(event(Level::Warn, "x", Some("k")));
        s.ingest(event(Level::Warn, "x", Some("k")));
        assert_eq!(s.len(), 1);
        assert_eq!(s.entries().next().unwrap().count, 3);
        assert_eq!(s.unread_count(), 3);
    }

    #[test]
    fn different_coalesce_keys_do_not_merge() {
        let mut s = NotificationStore::new();
        s.ingest(event(Level::Warn, "x", Some("k1")));
        s.ingest(event(Level::Warn, "y", Some("k2")));
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn different_levels_do_not_merge_even_with_same_key() {
        let mut s = NotificationStore::new();
        s.ingest(event(Level::Info, "x", Some("k")));
        s.ingest(event(Level::Warn, "x", Some("k")));
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn ring_evicts_oldest() {
        let mut s = NotificationStore::new();
        for _ in 0..(MAX_ENTRIES + 5) {
            s.ingest(event(Level::Info, "x", None));
        }
        assert_eq!(s.len(), MAX_ENTRIES);
    }

    #[test]
    fn clear_resets_unread() {
        let mut s = NotificationStore::new();
        s.ingest(event(Level::Warn, "x", None));
        s.clear();
        assert!(s.is_empty());
        assert_eq!(s.unread_count(), 0);
    }

    #[test]
    fn mark_all_read_preserves_entries() {
        let mut s = NotificationStore::new();
        s.ingest(event(Level::Warn, "x", None));
        s.mark_all_read();
        assert_eq!(s.unread_count(), 0);
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn count_at_level_filters_correctly() {
        let mut s = NotificationStore::new();
        s.ingest(event(Level::Info, "i", None));
        s.ingest(event(Level::Warn, "w", None));
        s.ingest(event(Level::Warn, "w", None));
        s.ingest(event(Level::Error, "e", None));
        assert_eq!(s.count_at_level(Level::Info), 1);
        assert_eq!(s.count_at_level(Level::Warn), 2);
        assert_eq!(s.count_at_level(Level::Error), 1);
    }
}
