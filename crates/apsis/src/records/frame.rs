//! Binary frame codec for Apsis Records.
//!
//! See `docs/adr/011-apsis-record.md` §"File format" for byte-level layout.

use std::io::{self, Read, Write};

pub const KIND_SNAPSHOT: u8 = 0x00;
pub const KIND_DIAGNOSTIC: u8 = 0x03;
pub const KIND_RESUME_STATE: u8 = 0x04;
pub const KIND_TRAILER: u8 = 0xFF;

/// Per-body dynamic state.
#[derive(Debug, Clone, PartialEq)]
pub struct BodyState {
    pub pos: [f64; 3],
    pub vel: [f64; 3],
}

#[derive(Debug, Clone, PartialEq)]
pub struct Snapshot {
    pub t: f64,
    pub bodies: Vec<BodyState>,
}

/// Relative drift `(X − X₀) / |X₀|` of conserved scalars at time `t`.
#[derive(Debug, Clone, PartialEq)]
pub struct Diagnostic {
    pub t: f64,
    pub d_energy_rel: f64,
    pub d_lz_rel: f64,
}

/// Per-integrator scratch state captured at the moment of a Snapshot.
/// Format is integrator-internal; consumers route by the record header's
/// `integrator.kind` to the matching [`Integrator::restore_resume_state`]
/// implementation. `step_count` carries the System-level step counter
/// at capture so periodic schedules (e.g. COM recentering every 97
/// steps) resume on the same cadence as the original run.
#[derive(Debug, Clone, PartialEq)]
pub struct ResumeState {
    pub t: f64,
    pub step_count: u64,
    pub bytes: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Trailer {
    pub t: f64,
    pub step_count: u64,
    pub frame_count: u64,
    pub blake3: [u8; 32],
}

#[derive(Debug, Clone, PartialEq)]
pub enum Frame {
    Snapshot(Snapshot),
    Diagnostic(Diagnostic),
    ResumeState(ResumeState),
    Trailer(Trailer),
}

impl Frame {
    pub fn write<W: Write>(&self, w: &mut W) -> io::Result<()> {
        match self {
            Frame::Snapshot(s) => {
                let n = s.bodies.len() as u32;
                let payload_len = 4 + n as usize * 48;
                write_header(w, KIND_SNAPSHOT, s.t, payload_len as u32)?;
                w.write_all(&n.to_le_bytes())?;
                for b in &s.bodies {
                    for c in b.pos.iter().chain(b.vel.iter()) {
                        w.write_all(&c.to_le_bytes())?;
                    }
                }
            },
            Frame::Diagnostic(d) => {
                let payload_len = 2 * 8;
                write_header(w, KIND_DIAGNOSTIC, d.t, payload_len as u32)?;
                w.write_all(&d.d_energy_rel.to_le_bytes())?;
                w.write_all(&d.d_lz_rel.to_le_bytes())?;
            },
            Frame::ResumeState(rs) => {
                let n = rs.bytes.len() as u32;
                let payload_len = 8 + 4 + rs.bytes.len();
                write_header(w, KIND_RESUME_STATE, rs.t, payload_len as u32)?;
                w.write_all(&rs.step_count.to_le_bytes())?;
                w.write_all(&n.to_le_bytes())?;
                w.write_all(&rs.bytes)?;
            },
            Frame::Trailer(tr) => {
                let payload_len = 8 + 8 + 32;
                write_header(w, KIND_TRAILER, tr.t, payload_len as u32)?;
                w.write_all(&tr.step_count.to_le_bytes())?;
                w.write_all(&tr.frame_count.to_le_bytes())?;
                w.write_all(&tr.blake3)?;
            },
        }
        Ok(())
    }

    pub fn read<R: Read>(r: &mut R) -> io::Result<Option<Self>> {
        let mut kind = [0u8; 1];
        match r.read(&mut kind)? {
            0 => return Ok(None),
            1 => {},
            _ => unreachable!("read of 1 byte returned >1"),
        }
        let mut t_buf = [0u8; 8];
        r.read_exact(&mut t_buf)?;
        let t = f64::from_le_bytes(t_buf);
        let mut len_buf = [0u8; 4];
        r.read_exact(&mut len_buf)?;
        let payload_len = u32::from_le_bytes(len_buf) as usize;
        let mut payload = vec![0u8; payload_len];
        r.read_exact(&mut payload)?;

        let frame = match kind[0] {
            KIND_SNAPSHOT => {
                let n = u32::from_le_bytes(payload[..4].try_into().unwrap()) as usize;
                let mut bodies = Vec::with_capacity(n);
                let mut off = 4;
                for _ in 0..n {
                    let pos = [
                        f64::from_le_bytes(payload[off..off + 8].try_into().unwrap()),
                        f64::from_le_bytes(payload[off + 8..off + 16].try_into().unwrap()),
                        f64::from_le_bytes(payload[off + 16..off + 24].try_into().unwrap()),
                    ];
                    let vel = [
                        f64::from_le_bytes(payload[off + 24..off + 32].try_into().unwrap()),
                        f64::from_le_bytes(payload[off + 32..off + 40].try_into().unwrap()),
                        f64::from_le_bytes(payload[off + 40..off + 48].try_into().unwrap()),
                    ];
                    bodies.push(BodyState { pos, vel });
                    off += 48;
                }
                Frame::Snapshot(Snapshot { t, bodies })
            },
            KIND_DIAGNOSTIC => {
                let d_energy_rel = f64::from_le_bytes(payload[..8].try_into().unwrap());
                let d_lz_rel = f64::from_le_bytes(payload[8..16].try_into().unwrap());
                Frame::Diagnostic(Diagnostic { t, d_energy_rel, d_lz_rel })
            },
            KIND_RESUME_STATE => {
                let step_count = u64::from_le_bytes(payload[..8].try_into().unwrap());
                let n = u32::from_le_bytes(payload[8..12].try_into().unwrap()) as usize;
                let bytes = payload[12..12 + n].to_vec();
                Frame::ResumeState(ResumeState { t, step_count, bytes })
            },
            KIND_TRAILER => {
                let step_count = u64::from_le_bytes(payload[..8].try_into().unwrap());
                let frame_count = u64::from_le_bytes(payload[8..16].try_into().unwrap());
                let mut blake3 = [0u8; 32];
                blake3.copy_from_slice(&payload[16..48]);
                Frame::Trailer(Trailer { t, step_count, frame_count, blake3 })
            },
            other => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("unknown frame kind: 0x{other:02X}"),
                ));
            },
        };
        Ok(Some(frame))
    }
}

fn write_header<W: Write>(w: &mut W, kind: u8, t: f64, payload_len: u32) -> io::Result<()> {
    w.write_all(&[kind])?;
    w.write_all(&t.to_le_bytes())?;
    w.write_all(&payload_len.to_le_bytes())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn round_trip(f: &Frame) -> Frame {
        let mut buf = Vec::new();
        f.write(&mut buf).unwrap();
        Frame::read(&mut Cursor::new(buf)).unwrap().unwrap()
    }

    #[test]
    fn snapshot_round_trip() {
        let snap = Frame::Snapshot(Snapshot {
            t: 1.5,
            bodies: vec![
                BodyState { pos: [1.0, 2.0, 3.0], vel: [-0.1, 0.2, -0.3] },
                BodyState { pos: [0.0, 0.0, 0.0], vel: [0.0, 1.0, 0.0] },
            ],
        });
        assert_eq!(round_trip(&snap), snap);
    }

    #[test]
    fn diagnostic_round_trip() {
        let d = Frame::Diagnostic(Diagnostic { t: 5.0, d_energy_rel: -3.1e-13, d_lz_rel: 1.2e-14 });
        assert_eq!(round_trip(&d), d);
    }

    #[test]
    fn resume_state_round_trip() {
        let rs = Frame::ResumeState(ResumeState {
            t: 7.25,
            step_count: 12345,
            bytes: vec![0xDE, 0xAD, 0xBE, 0xEF, 0x01, 0x02, 0x03],
        });
        assert_eq!(round_trip(&rs), rs);
    }

    #[test]
    fn resume_state_empty_payload_round_trips() {
        let rs = Frame::ResumeState(ResumeState { t: 0.0, step_count: 0, bytes: vec![] });
        assert_eq!(round_trip(&rs), rs);
    }

    #[test]
    fn trailer_round_trip() {
        let tr = Frame::Trailer(Trailer {
            t: 1000.0,
            step_count: 999_999,
            frame_count: 42,
            blake3: [0xAB; 32],
        });
        assert_eq!(round_trip(&tr), tr);
    }

    #[test]
    fn read_returns_none_at_eof() {
        let mut empty = Cursor::new(Vec::<u8>::new());
        assert!(Frame::read(&mut empty).unwrap().is_none());
    }
}
