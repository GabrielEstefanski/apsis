//! `RecordHook` — writer for apsis records, implemented as a `SimHook`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::core::hooks::{Command, HookContext, SimHook};
use crate::records::format::{FORMAT_VER, MAGIC};
use crate::records::frame::{BodyState, Diagnostic, Frame, ResumeState, Snapshot, Trailer};
use crate::records::header::Header;
use crate::records::policy::{DiagnosticCadence, RecordPolicy};

pub struct RecordHook {
    writer: BufWriter<File>,
    hasher: blake3::Hasher,
    header: Header,
    header_written: bool,
    policy: RecordPolicy,
    diagnostics: DiagnosticCadence,
    capture_resume: bool,
    t_last_snapshot: Option<f64>,
    t_last_diagnostic: Option<f64>,
    frame_count: u64,
    /// Set once `on_finish` runs — `Drop` then skips its safety-net
    /// flush, since the writer is already closed properly.
    closed: bool,
}

impl RecordHook {
    /// Construct with a fully-built [`Header`]. For the common case of
    /// gathering the header from a [`System`], pair this with
    /// [`crate::records::provenance::header_from_system`].
    pub fn with_header(
        path: impl AsRef<Path>,
        header: Header,
        policy: RecordPolicy,
    ) -> std::io::Result<Self> {
        let file = File::create(path.as_ref())?;
        Ok(Self {
            writer: BufWriter::new(file),
            hasher: blake3::Hasher::new(),
            header,
            header_written: false,
            policy,
            diagnostics: DiagnosticCadence::Off,
            capture_resume: false,
            t_last_snapshot: None,
            t_last_diagnostic: None,
            frame_count: 0,
            closed: false,
        })
    }

    /// Enable periodic emission of `Diagnostic` frames (ΔE/E, ΔLz/Lz)
    /// at the given cadence. Default is [`DiagnosticCadence::Off`].
    pub fn with_diagnostics(mut self, cadence: DiagnosticCadence) -> Self {
        self.diagnostics = cadence;
        self
    }

    /// Emit a `ResumeState` frame alongside every `Snapshot`. Required
    /// for mid-run resume via [`crate::records::Record::resume_from`];
    /// off by default because IAS15's serialised scratch is several KB
    /// per snapshot.
    pub fn with_resume_capture(mut self, enabled: bool) -> Self {
        self.capture_resume = enabled;
        self
    }

    fn write_file_header(&mut self) -> std::io::Result<()> {
        let toml = self
            .header
            .to_toml()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let toml_bytes = toml.as_bytes();
        let header_len = toml_bytes.len() as u64;

        let mut prefix = Vec::with_capacity(16);
        prefix.extend_from_slice(MAGIC);
        prefix.extend_from_slice(&FORMAT_VER.to_le_bytes());
        prefix.extend_from_slice(&0u16.to_le_bytes());
        prefix.extend_from_slice(&header_len.to_le_bytes());

        // The trailer's BLAKE3 covers the frame stream only — the header
        // is plaintext + re-parseable, and `created_utc` is per-run
        // wall-clock metadata. Hashing the header would couple the
        // content-addressable trailer to the timestamp, breaking the
        // "same {seed, config} → byte-equal frame stream + trailer"
        // contract the reproducibility gate (records::tests::reproducibility)
        // exercises.
        self.writer.write_all(&prefix)?;
        self.writer.write_all(toml_bytes)?;
        Ok(())
    }

    fn write_frame(&mut self, frame: &Frame) -> std::io::Result<()> {
        let mut buf = Vec::new();
        frame.write(&mut buf)?;
        self.hasher.update(&buf);
        self.writer.write_all(&buf)?;
        self.frame_count += 1;
        Ok(())
    }

    fn snapshot_from_ctx(ctx: &HookContext<'_>) -> Snapshot {
        let bodies = ctx
            .bodies
            .iter()
            .map(|b| BodyState {
                pos: [b.pos_x, b.pos_y, b.pos_z],
                vel: [b.vel_x, b.vel_y, b.vel_z],
            })
            .collect();
        Snapshot { t: ctx.t, bodies }
    }

    fn diagnostic_from_ctx(ctx: &HookContext<'_>) -> Diagnostic {
        // Precision-limited regime serialises as NaN — readers detect
        // with `.is_nan()`. Matches the CSV recorder convention.
        Diagnostic {
            t: ctx.t,
            d_energy_rel: ctx.rel_energy_error.unwrap_or(f64::NAN),
            d_lz_rel: ctx.rel_angular_momentum_error.unwrap_or(f64::NAN),
        }
    }

    fn maybe_write_resume(&mut self, ctx: &HookContext<'_>) {
        if !self.capture_resume {
            return;
        }
        let bytes = ctx.resume_state.clone().unwrap_or_default();
        self.write_frame(&Frame::ResumeState(ResumeState {
            t: ctx.t,
            step_count: ctx.steps,
            bytes,
        }))
        .expect("RecordHook: write resume state");
    }
}

// Writes panic on I/O failure: the hook has no Result channel into
// the integrator loop, and a silent partial record verifies as valid.
impl SimHook for RecordHook {
    fn name(&self) -> &'static str {
        "RecordHook"
    }

    fn wants_resume_state(&self) -> bool {
        self.capture_resume
    }

    fn pre_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
        if !self.header_written {
            self.write_file_header().expect("RecordHook: write header");
            let snap = Self::snapshot_from_ctx(ctx);
            self.write_frame(&Frame::Snapshot(snap)).expect("RecordHook: write initial bookend");
            self.header_written = true;
            self.t_last_snapshot = Some(ctx.t);
            self.maybe_write_resume(ctx);
            if self.diagnostics != DiagnosticCadence::Off {
                let d = Self::diagnostic_from_ctx(ctx);
                self.write_frame(&Frame::Diagnostic(d))
                    .expect("RecordHook: write initial diagnostic");
                self.t_last_diagnostic = Some(ctx.t);
            }
        }
        Vec::new()
    }

    fn post_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
        if self.policy.should_snapshot(ctx.t, ctx.steps, self.t_last_snapshot) {
            let snap = Self::snapshot_from_ctx(ctx);
            self.write_frame(&Frame::Snapshot(snap)).expect("RecordHook: write snapshot");
            self.t_last_snapshot = Some(ctx.t);
            self.maybe_write_resume(ctx);
        }
        if self.diagnostics.should_emit(ctx.t, ctx.steps, self.t_last_diagnostic) {
            let d = Self::diagnostic_from_ctx(ctx);
            self.write_frame(&Frame::Diagnostic(d)).expect("RecordHook: write diagnostic");
            self.t_last_diagnostic = Some(ctx.t);
        }
        Vec::new()
    }

    fn on_finish(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
        if !self.header_written || self.closed {
            return Vec::new();
        }
        // Final bookend (skip when already emitted at this t).
        if self.t_last_snapshot != Some(ctx.t) {
            let snap = Self::snapshot_from_ctx(ctx);
            self.write_frame(&Frame::Snapshot(snap)).expect("RecordHook: write final bookend");
        }
        // Snapshot the hasher before writing the trailer so it isn't covered.
        let frames_blake3 = self.hasher.clone().finalize();
        let trailer = Trailer {
            t: ctx.t,
            step_count: ctx.steps,
            frame_count: self.frame_count,
            blake3: *frames_blake3.as_bytes(),
        };
        self.write_frame(&Frame::Trailer(trailer)).expect("RecordHook: write trailer");
        self.writer.flush().expect("RecordHook: flush writer");
        self.closed = true;
        Vec::new()
    }
}

impl Drop for RecordHook {
    fn drop(&mut self) {
        // Semantic close runs in `on_finish`. This is a flush-only
        // safety net for paths that bypass it; the resulting record
        // has no trailer and `Record::open` rejects it as malformed.
        let _ = self.writer.flush();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::hooks::{HookPhase, HookPhaseKind};
    use crate::domain::body::Body;
    use crate::records::format::MAGIC;
    use crate::records::header::{
        Apsis, BodiesMeta, IntegratorMeta, KernelMeta, Reproducibility, UnitSystemMeta,
    };

    fn fake_header() -> Header {
        Header {
            apsis: Apsis {
                version: "0.1.0".into(),
                git_sha: "test".into(),
                created_utc: "2026-05-16T00:00:00Z".into(),
                rustc_version: "".into(),
                generated_by: "apsis-test".into(),
            },
            reproducibility: Reproducibility { cargo_lock_blake3: "00".repeat(32), seed: 0 },
            unit_system: UnitSystemMeta {
                g: 1.0,
                length: "AU".into(),
                mass: "M_sun".into(),
                time: "yr/2pi".into(),
            },
            integrator: IntegratorMeta {
                kind: "Ias15".into(),
                dt_mode: "Fixed".into(),
                initial_dt: 1e-3,
                params: Default::default(),
            },
            kernel: KernelMeta {
                variant: "Newton".into(),
                softening: None,
                exactness: None,
                continuity: None,
            },
            operators: vec![],
            bodies: BodiesMeta { count: 0, list: vec![] },
        }
    }

    fn make_ctx<'a>(bodies: &'a [Body], t: f64, steps: u64) -> HookContext<'a> {
        HookContext {
            bodies,
            t,
            dt: 1e-3,
            steps,
            rel_energy_error: None,
            rel_angular_momentum_error: None,
            phase: HookPhase(HookPhaseKind::PreStep),
            resume_state: None,
        }
    }

    #[test]
    fn hook_writes_magic_and_initial_bookend() {
        let tmp = std::env::temp_dir().join("apsis-record-test-hook-1.apsis");
        let _ = std::fs::remove_file(&tmp);
        {
            let mut hook =
                RecordHook::with_header(&tmp, fake_header(), RecordPolicy::BookendsAndEvents)
                    .unwrap();
            let bodies = vec![Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0)];
            let ctx = make_ctx(&bodies, 0.0, 0);
            hook.pre_step(&ctx);
            // flushed on drop
        }
        let bytes = std::fs::read(&tmp).unwrap();
        assert_eq!(&bytes[..4], MAGIC, "MAGIC mismatch");
        let _ = std::fs::remove_file(&tmp);
    }
}
