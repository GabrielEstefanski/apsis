//! Lightweight structured diagnostics emitter — a thin stand-in for
//! `tracing::warn!` until the project adopts the `tracing` crate
//! properly.
//!
//! Call sites use the [`warn_diag!`](crate::warn_diag) macro with a
//! literal message and optional `key = value` fields:
//!
//! ```ignore
//! use crate::warn_diag;
//!
//! warn_diag!(
//!     "integrator requires deterministic force",
//!     integrator = "IAS15",
//!     exact_threshold_before = 64_usize,
//!     exact_threshold_after = 10_000_usize,
//! );
//! ```
//!
//! The output goes to stderr as a single line with the shape:
//!
//! ```text
//! [gravity-sim WARN] integrator requires deterministic force
//!                    { integrator="IAS15" exact_threshold_before=64 exact_threshold_after=10000 }
//! ```
//!
//! # Evolution path
//!
//! When the project adopts `tracing`, rewrite the macro body to call
//! `tracing::warn!` with field syntax. Call sites do not change.
//! The single call-site convention already mirrors `tracing`'s
//! key-value field model so the migration is mechanical.
//!
//! # Why not `eprintln!` everywhere
//!
//! Two reasons:
//!
//! * **Consistency.** One prefix, one format, one place to change
//!   when we do switch to `tracing`.
//! * **Grepability.** `[gravity-sim WARN]` is a stable, unique
//!   anchor for downstream consumers (bench log scrubbers, UI toast
//!   channels once they subscribe).

/// Emit a structured warning to stderr with a message and optional
/// `key = value` fields.
///
/// See the module-level docs for format and migration notes.
#[macro_export]
macro_rules! warn_diag {
    ($msg:literal $(,)?) => {
        eprintln!("[gravity-sim WARN] {}", $msg);
    };
    ($msg:literal, $($key:ident = $val:expr),+ $(,)?) => {{
        use std::fmt::Write as _;
        let mut buf = String::with_capacity(128);
        let _ = write!(&mut buf, "[gravity-sim WARN] {} {{", $msg);
        let mut first = true;
        $(
            if !first { let _ = write!(&mut buf, " "); }
            let _ = write!(&mut buf, "{}={:?}", stringify!($key), &$val);
            #[allow(unused_assignments)] { first = false; }
        )+
        let _ = write!(&mut buf, " }}");
        eprintln!("{}", buf);
    }};
}
