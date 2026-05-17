//! Binary frame codec for Apsis Records.
//!
//! See `docs/adr/011-apsis-record.md` §"File format" for byte-level layout.

use std::io::{self, Read, Write};

pub const KIND_SNAPSHOT: u8 = 0x00;
pub const KIND_COLLISION: u8 = 0x01;
pub const KIND_ESCAPE: u8 = 0x02;
pub const KIND_DIAGNOSTIC: u8 = 0x03;
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

#[derive(Debug, Clone, PartialEq)]
pub enum Event {
    Collision { t: f64, body_a: u32, body_b: u32, distance: f64 },
    Escape { t: f64, body: u32, radius: f64 },
}

/// Relative drift `(X − X₀) / |X₀|` of conserved scalars at time `t`.
#[derive(Debug, Clone, PartialEq)]
pub struct Diagnostic {
    pub t: f64,
    pub d_energy_rel: f64,
    pub d_lz_rel: f64,
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
    Event(Event),
    Diagnostic(Diagnostic),
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
            Frame::Event(Event::Collision { t, body_a, body_b, distance }) => {
                let payload_len = 1 + 4 + 4 + 8;
                write_header(w, KIND_COLLISION, *t, payload_len as u32)?;
                // Sub-kind byte reserved for forward compat with operator-emitted
                // events (0x10–0xFE range, post-v0.1). Currently always 0.
                w.write_all(&[0u8])?;
                w.write_all(&body_a.to_le_bytes())?;
                w.write_all(&body_b.to_le_bytes())?;
                w.write_all(&distance.to_le_bytes())?;
            },
            Frame::Event(Event::Escape { t, body, radius }) => {
                let payload_len = 1 + 4 + 8;
                write_header(w, KIND_ESCAPE, *t, payload_len as u32)?;
                w.write_all(&[0u8])?;
                w.write_all(&body.to_le_bytes())?;
                w.write_all(&radius.to_le_bytes())?;
            },
            Frame::Diagnostic(d) => {
                let payload_len = 2 * 8;
                write_header(w, KIND_DIAGNOSTIC, d.t, payload_len as u32)?;
                w.write_all(&d.d_energy_rel.to_le_bytes())?;
                w.write_all(&d.d_lz_rel.to_le_bytes())?;
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
            KIND_COLLISION => {
                // payload[0] = sub-kind (currently always 0; reserved for future variants)
                let body_a = u32::from_le_bytes(payload[1..5].try_into().unwrap());
                let body_b = u32::from_le_bytes(payload[5..9].try_into().unwrap());
                let distance = f64::from_le_bytes(payload[9..17].try_into().unwrap());
                Frame::Event(Event::Collision { t, body_a, body_b, distance })
            },
            KIND_ESCAPE => {
                let body = u32::from_le_bytes(payload[1..5].try_into().unwrap());
                let radius = f64::from_le_bytes(payload[5..13].try_into().unwrap());
                Frame::Event(Event::Escape { t, body, radius })
            },
            KIND_DIAGNOSTIC => {
                let d_energy_rel = f64::from_le_bytes(payload[..8].try_into().unwrap());
                let d_lz_rel = f64::from_le_bytes(payload[8..16].try_into().unwrap());
                Frame::Diagnostic(Diagnostic { t, d_energy_rel, d_lz_rel })
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
    fn collision_round_trip() {
        let ev = Frame::Event(Event::Collision { t: 12.4, body_a: 0, body_b: 3, distance: 1e-5 });
        assert_eq!(round_trip(&ev), ev);
    }

    #[test]
    fn escape_round_trip() {
        let ev = Frame::Event(Event::Escape { t: 200.0, body: 7, radius: 1e6 });
        assert_eq!(round_trip(&ev), ev);
    }

    #[test]
    fn diagnostic_round_trip() {
        let d = Frame::Diagnostic(Diagnostic {
            t: 5.0,
            d_energy_rel: -3.1e-13,
            d_lz_rel: 1.2e-14,
        });
        assert_eq!(round_trip(&d), d);
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
