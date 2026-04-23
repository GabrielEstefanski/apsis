//! Structured diagnostic events + process-global event bus.
//!
//! The entry points for producers are the macros
//! [`warn_diag!`](crate::warn_diag), [`info_diag!`](crate::info_diag)
//! and [`error_diag!`](crate::error_diag). Each constructs an
//! [`Event`] and publishes it through [`publish`].
//!
//! Consumers register via [`subscribe`], which returns a
//! [`SubscriptionId`] for later deregistration. A default subscriber
//! writes every event to stderr in the pre-existing
//! `[gravity-sim WARN] …` format, so headless runs and the
//! benchmark harness behave identically to before the bus landed.
//!
//! # Usage
//!
//! ```ignore
//! use crate::core::log::Source;
//!
//! crate::warn_diag!(
//!     Source::Integrator,
//!     "IAS15 dt floor reached; controller accepted degraded step",
//!     dt = 1e-12_f64,
//!     floor_hit_count = 42_u64,
//! );
//! ```
//!
//! The macro captures `stringify!(dt)` / `stringify!(floor_hit_count)`
//! as field names automatically; values are formatted with `Debug`.
//!
//! # Why a bus, not direct `eprintln!`
//!
//! The UI will subscribe to render a notification center (Precision
//! Run mode), so producers cannot assume stderr is the only consumer.
//! The event-bus layer keeps production sites agnostic of who reads.

pub mod bus;
pub mod event;

pub use bus::{publish, subscribe, unsubscribe, SubscriptionId};
pub use event::{Event, Level, Source};

// ── Producer macros ──────────────────────────────────────────────────────────
//
// Each macro is a thin wrapper over `Event::new(...)` + `publish(event)`. The
// ergonomic win is field capture (`key = value` → `(stringify!(key),
// format!("{:?}", value))`) and the fixed `Level` per macro name.
//
// Macros are exported at the crate root via `#[macro_export]`. Call sites
// write `crate::warn_diag!(...)` regardless of where the macro definition
// physically lives — the module path is not part of the public signature.

/// Publish a `Level::Warn` event. See module-level docs for usage.
///
/// Two forms:
///
/// * `warn_diag!(source, "literal message")` — no fields.
/// * `warn_diag!(source, "literal message", key = value, ...)` — one or
///   more structured fields. Values use `Debug` formatting.
#[macro_export]
macro_rules! warn_diag {
    ($source:expr, $msg:literal $(,)?) => {
        $crate::core::log::publish(
            $crate::core::log::Event::new(
                $crate::core::log::Level::Warn,
                $source,
                $msg,
            )
        );
    };
    ($source:expr, $msg:literal, $($key:ident = $val:expr),+ $(,)?) => {
        $crate::core::log::publish({
            let mut event = $crate::core::log::Event::new(
                $crate::core::log::Level::Warn,
                $source,
                $msg,
            );
            $(
                event = event.with_field(stringify!($key), format!("{:?}", &$val));
            )+
            event
        });
    };
}

/// Publish a `Level::Info` event. Same signature as
/// [`warn_diag!`](crate::warn_diag).
#[macro_export]
macro_rules! info_diag {
    ($source:expr, $msg:literal $(,)?) => {
        $crate::core::log::publish(
            $crate::core::log::Event::new(
                $crate::core::log::Level::Info,
                $source,
                $msg,
            )
        );
    };
    ($source:expr, $msg:literal, $($key:ident = $val:expr),+ $(,)?) => {
        $crate::core::log::publish({
            let mut event = $crate::core::log::Event::new(
                $crate::core::log::Level::Info,
                $source,
                $msg,
            );
            $(
                event = event.with_field(stringify!($key), format!("{:?}", &$val));
            )+
            event
        });
    };
}

/// Publish a `Level::Error` event. Same signature as
/// [`warn_diag!`](crate::warn_diag).
#[macro_export]
macro_rules! error_diag {
    ($source:expr, $msg:literal $(,)?) => {
        $crate::core::log::publish(
            $crate::core::log::Event::new(
                $crate::core::log::Level::Error,
                $source,
                $msg,
            )
        );
    };
    ($source:expr, $msg:literal, $($key:ident = $val:expr),+ $(,)?) => {
        $crate::core::log::publish({
            let mut event = $crate::core::log::Event::new(
                $crate::core::log::Level::Error,
                $source,
                $msg,
            );
            $(
                event = event.with_field(stringify!($key), format!("{:?}", &$val));
            )+
            event
        });
    };
}

#[cfg(test)]
mod macro_tests {
    use super::*;
    use std::sync::{Arc, Mutex};

    #[test]
    fn warn_diag_with_fields_round_trips_through_bus() {
        let captured: Arc<Mutex<Option<Event>>> = Arc::new(Mutex::new(None));
        let sink = captured.clone();
        let id = subscribe(move |event: &Event| {
            *sink.lock().unwrap() = Some(event.clone());
        });

        let dt: f64 = 1.5e-12;
        let count: u64 = 7;
        crate::warn_diag!(
            Source::Integrator,
            "unit test — floor reached",
            dt = dt,
            count = count,
        );

        let got = captured.lock().unwrap().clone().expect("event should have been captured");
        assert_eq!(got.level, Level::Warn);
        assert_eq!(got.source, Source::Integrator);
        assert_eq!(got.message, "unit test — floor reached");
        assert_eq!(got.fields.len(), 2);
        assert_eq!(got.fields[0].0, "dt");
        assert_eq!(got.fields[1].0, "count");
        assert!(got.fields[0].1.contains("1.5e-12"), "dt debug formatting preserved");

        unsubscribe(id);
    }
}
