//! `RecordHook` — writer for apsis records, implemented as a `SimHook`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::core::hooks::{CollisionEvent, Command, EscapeEvent, HookContext, SimHook};
use crate::records::frame::{BodyState, Event, Frame, Snapshot, Trailer};
use crate::records::header::Header;
use crate::records::policy::RecordPolicy;

/// File-format version embedded after the magic. Bumping requires the
/// `tests::schema_version` pin + an ADR update.
pub const FORMAT_VER: u16 = 1;
pub const MAGIC: &[u8; 4] = b"APSR";

pub struct RecordHook {
    path: PathBuf,
    writer: BufWriter<File>,
    hasher: blake3::Hasher,
    header: Header,
    header_written: bool,
    policy: RecordPolicy,
    t_last_snapshot: Option<f64>,
    /// Cached pre-step / post-step snapshot used to write the final
    /// bookend on `Drop`. Updated on every `post_step`; the most recent
    /// value is the "final" state by definition.
    last_state: Option<Snapshot>,
    last_steps: u64,
    frame_count: u64,
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
            path: path.as_ref().to_path_buf(),
            writer: BufWriter::new(file),
            hasher: blake3::Hasher::new(),
            header,
            header_written: false,
            policy,
            t_last_snapshot: None,
            last_state: None,
            last_steps: 0,
            frame_count: 0,
        })
    }

    /// The file the record is being written to.
    pub fn path(&self) -> &Path {
        &self.path
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

        self.hasher.update(&prefix);
        self.writer.write_all(&prefix)?;
        self.hasher.update(toml_bytes);
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
}

impl SimHook for RecordHook {
    fn name(&self) -> &'static str {
        "RecordHook"
    }

    fn pre_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
        if !self.header_written {
            self.write_file_header().expect("RecordHook: write header");
            let snap = Self::snapshot_from_ctx(ctx);
            self.write_frame(&Frame::Snapshot(snap.clone()))
                .expect("RecordHook: write initial bookend");
            self.header_written = true;
            self.t_last_snapshot = Some(ctx.t);
            self.last_state = Some(snap);
            self.last_steps = ctx.steps;
        }
        Vec::new()
    }

    fn post_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
        let snap = Self::snapshot_from_ctx(ctx);
        if self.policy.should_snapshot(ctx.t, ctx.steps, self.t_last_snapshot) {
            self.write_frame(&Frame::Snapshot(snap.clone())).expect("RecordHook: write snapshot");
            self.t_last_snapshot = Some(ctx.t);
        }
        // Cache the most recent state for the final bookend on Drop.
        self.last_state = Some(snap);
        self.last_steps = ctx.steps;
        Vec::new()
    }

    fn on_collision(&mut self, ev: &CollisionEvent, _ctx: &HookContext<'_>) -> Vec<Command> {
        let f = Frame::Event(Event::Collision {
            t: ev.t,
            body_a: ev.i as u32,
            body_b: ev.j as u32,
            distance: ev.separation,
        });
        self.write_frame(&f).expect("RecordHook: write collision");
        Vec::new()
    }

    fn on_escape(&mut self, ev: &EscapeEvent, _ctx: &HookContext<'_>) -> Vec<Command> {
        let f = Frame::Event(Event::Escape { t: ev.t, body: ev.body as u32, radius: ev.radius });
        self.write_frame(&f).expect("RecordHook: write escape");
        Vec::new()
    }
}

impl Drop for RecordHook {
    fn drop(&mut self) {
        if !self.header_written {
            return;
        }
        // Write the final bookend Snapshot from the cached last state.
        // If the last cached state is at the same time as t_last_snapshot
        // (i.e. the writer already emitted it under a non-bookend policy),
        // skip the duplicate; otherwise emit.
        let cached = self.last_state.take();
        let final_t = if let Some(snap) = cached {
            let must_emit_bookend = match self.t_last_snapshot {
                None => true,
                Some(t) => (snap.t - t).abs() > f64::EPSILON,
            };
            let snap_t = snap.t;
            if must_emit_bookend {
                let _ = self.write_frame(&Frame::Snapshot(snap));
            }
            snap_t
        } else {
            0.0
        };

        // Snapshot the hasher BEFORE writing the trailer. The trailer's
        // payload carries this digest; the bytes we write next are the
        // trailer itself, which is excluded from its own hash by
        // construction.
        let pre_trailer_hash = self.hasher.clone().finalize();
        let trailer = Trailer {
            t: final_t,
            step_count: self.last_steps,
            frame_count: self.frame_count + 1,
            blake3: *pre_trailer_hash.as_bytes(),
        };
        let mut out = Vec::new();
        if Frame::Trailer(trailer).write(&mut out).is_ok() {
            let _ = self.writer.write_all(&out);
        }
        let _ = self.writer.flush();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::hooks::{HookPhase, HookPhaseKind};
    use crate::domain::body::Body;
    use crate::records::header::{
        Apsis, BodiesMeta, IntegratorMeta, KernelMeta, Reproducibility, UnitSystemMeta,
    };

    fn fake_header() -> Header {
        Header {
            apsis: Apsis {
                version: "0.1.0".into(),
                git_sha: "test".into(),
                created_utc: "2026-05-16T00:00:00Z".into(),
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
            kernel: KernelMeta { variant: "Newton".into(), softening: None },
            operators: vec![],
            bodies: BodiesMeta { count: 0, list: vec![] },
        }
    }

    fn make_ctx<'a>(
        bodies: &'a [Body],
        names: &'a [String],
        t: f64,
        steps: u64,
    ) -> HookContext<'a> {
        HookContext {
            bodies,
            names,
            t,
            dt: 1e-3,
            steps,
            rel_energy_error: 0.0,
            rel_angular_momentum_error: 0.0,
            phase: HookPhase(HookPhaseKind::PreStep),
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
            let names = vec!["sun".to_string()];
            let ctx = make_ctx(&bodies, &names, 0.0, 0);
            hook.pre_step(&ctx);
            // Drop writes trailer.
        }
        let bytes = std::fs::read(&tmp).unwrap();
        assert_eq!(&bytes[..4], MAGIC, "MAGIC mismatch");
        let _ = std::fs::remove_file(&tmp);
    }
}
