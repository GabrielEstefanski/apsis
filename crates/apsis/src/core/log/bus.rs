//! Process-global event bus — single point of publication for
//! structured diagnostics.
//!
//! # Design
//!
//! A broadcast with synchronous fan-out. Every [`publish`] call
//! walks the subscriber list and invokes each callback with the
//! event. Callbacks must be fast and non-blocking; heavy work
//! (disk I/O, rendering) belongs in whatever the callback hands the
//! event off to.
//!
//! The bus is lazy-initialised via [`OnceLock`]. The first access
//! registers the default **stderr bridge** so existing behaviour
//! (line-oriented WARN / INFO / ERROR on stderr) is preserved for
//! headless runs, benches, and tests without any additional setup.
//!
//! # Why synchronous, not async
//!
//! Events are rare by construction (scenario-stiffness signals,
//! integrator mode switches, user actions). The overhead of a
//! channel + consumer task is unnecessary. Synchronous delivery
//! also gives us deterministic ordering, which matters for the UI's
//! notification feed.
//!
//! # Why `OnceLock`, not `Lazy`
//!
//! `std::sync::OnceLock` is the stable-Rust equivalent of
//! `once_cell::Lazy`. No extra dependency required.

use std::sync::{Arc, Mutex, OnceLock};

use super::event::Event;

/// Identifier returned by [`subscribe`]; pass it to
/// [`unsubscribe`] to deregister the callback cleanly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SubscriptionId(u64);

type Callback = Arc<dyn Fn(&Event) + Send + Sync + 'static>;

struct Subscriber {
    id: SubscriptionId,
    callback: Callback,
}

/// Subscriber registry. Public methods go through the global
/// accessors at module level ([`publish`], [`subscribe`],
/// [`unsubscribe`]); the type itself is opaque.
pub(crate) struct EventBus {
    subscribers: Mutex<Vec<Subscriber>>,
    next_id: Mutex<u64>,
}

impl EventBus {
    fn new_with_default_subscribers() -> Self {
        let bus = Self { subscribers: Mutex::new(Vec::new()), next_id: Mutex::new(0) };
        // Default subscriber: line-oriented stderr bridge. Matches
        // the pre-bus `eprintln!` output format for continuity with
        // bench log scrubbers and dev-loop terminals.
        bus.subscribe_inner(Arc::new(|event: &Event| {
            eprintln!("{}", event.format_single_line());
        }));
        bus
    }

    fn subscribe_inner(&self, callback: Callback) -> SubscriptionId {
        let id = {
            let mut next = self.next_id.lock().unwrap();
            let id = SubscriptionId(*next);
            *next = next.wrapping_add(1);
            id
        };
        self.subscribers.lock().unwrap().push(Subscriber { id, callback });
        id
    }

    fn publish_inner(&self, event: &Event) {
        // Clone the callbacks under the lock so the lock is not held
        // during user-supplied callback execution — a slow subscriber
        // should never block other publishers.
        let callbacks: Vec<Callback> = {
            let subs = self.subscribers.lock().unwrap();
            subs.iter().map(|s| s.callback.clone()).collect()
        };
        for cb in callbacks {
            cb(event);
        }
    }

    fn unsubscribe_inner(&self, id: SubscriptionId) {
        self.subscribers.lock().unwrap().retain(|s| s.id != id);
    }
}

static BUS: OnceLock<EventBus> = OnceLock::new();

fn bus() -> &'static EventBus {
    BUS.get_or_init(EventBus::new_with_default_subscribers)
}

/// Publish an event. Delivered synchronously to every current
/// subscriber. The call is infallible; subscribers that panic are
/// the subscriber's problem, not the publisher's.
pub fn publish(event: Event) {
    bus().publish_inner(&event);
}

/// Register a subscriber. The callback is invoked once per
/// [`publish`] call for the lifetime of the registration (or until
/// [`unsubscribe`] is called with the returned id).
///
/// Callbacks must be fast. Hand off long work (rendering, I/O) to
/// another thread or queue; holding up the bus delays every other
/// subscriber.
pub fn subscribe<F>(callback: F) -> SubscriptionId
where
    F: Fn(&Event) + Send + Sync + 'static,
{
    bus().subscribe_inner(Arc::new(callback))
}

/// Deregister a subscriber previously registered with
/// [`subscribe`]. Calling with an unknown id is a no-op.
pub fn unsubscribe(id: SubscriptionId) {
    bus().unsubscribe_inner(id);
}

#[cfg(test)]
mod tests {
    use super::super::event::{Event, Level, Source};
    use super::*;

    // The event bus is a process-global singleton; unit tests run in
    // parallel by default and all share it. To keep test-to-test
    // isolation, each test filters inbound events by a unique marker
    // string baked into `event.message`. Unrelated traffic from other
    // tests is ignored by the `.starts_with` guard.

    #[test]
    fn subscribe_receives_published_events() {
        const MARKER: &str = "bus_test_sub::";
        let received: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
        let sink = received.clone();
        let id = subscribe(move |event: &Event| {
            if event.message.starts_with(MARKER) {
                sink.lock().unwrap().push(event.message.to_string());
            }
        });

        publish(Event::new(Level::Info, Source::System, "bus_test_sub::first"));
        publish(Event::new(Level::Warn, Source::System, "bus_test_sub::second"));

        let seen = received.lock().unwrap().clone();
        assert!(seen.contains(&"bus_test_sub::first".to_string()));
        assert!(seen.contains(&"bus_test_sub::second".to_string()));

        unsubscribe(id);
    }

    #[test]
    fn unsubscribe_stops_delivery() {
        const MARKER: &str = "bus_test_unsub::";
        let received: Arc<Mutex<u32>> = Arc::new(Mutex::new(0));
        let sink = received.clone();
        let id = subscribe(move |event: &Event| {
            if event.message.starts_with(MARKER) {
                *sink.lock().unwrap() += 1;
            }
        });

        publish(Event::new(Level::Info, Source::System, "bus_test_unsub::a"));
        unsubscribe(id);
        publish(Event::new(Level::Info, Source::System, "bus_test_unsub::b"));

        assert_eq!(
            *received.lock().unwrap(),
            1,
            "only the first event should have been counted (tagged by marker)"
        );
    }
}
