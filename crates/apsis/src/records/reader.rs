//! Reader for `.apsis` records.

use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use crate::records::format::{FORMAT_VER, MAGIC};
use crate::records::frame::{Frame, Snapshot, Trailer};
use crate::records::header::Header;

#[derive(Debug)]
pub struct Record {
    path: PathBuf,
    header: Header,
    /// File-byte offset where the first frame begins.
    frames_start: u64,
    /// Cached trailer (validated on open).
    trailer: Trailer,
}

impl Record {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, RecordError> {
        let path = path.as_ref().to_path_buf();
        let mut file = BufReader::new(File::open(&path)?);

        let mut prefix = [0u8; 16];
        file.read_exact(&mut prefix)?;
        if &prefix[..4] != MAGIC {
            return Err(RecordError::BadMagic);
        }
        let fmt_ver = u16::from_le_bytes(prefix[4..6].try_into().unwrap());
        if fmt_ver != FORMAT_VER {
            return Err(RecordError::UnsupportedFormatVersion(fmt_ver));
        }
        let header_len = u64::from_le_bytes(prefix[8..16].try_into().unwrap()) as usize;

        let mut header_buf = vec![0u8; header_len];
        file.read_exact(&mut header_buf)?;
        let header_str =
            std::str::from_utf8(&header_buf).map_err(|_| RecordError::BadHeaderUtf8)?;
        let header: Header =
            Header::from_toml(header_str).map_err(|e| RecordError::HeaderParse(e.to_string()))?;

        let frames_start = (16 + header_len) as u64;

        let mut last_trailer: Option<Trailer> = None;
        file.seek(SeekFrom::Start(frames_start))?;
        while let Some(frame) = Frame::read(&mut file)? {
            if let Frame::Trailer(t) = frame {
                last_trailer = Some(t);
                break;
            }
        }
        let trailer = last_trailer.ok_or(RecordError::MissingTrailer)?;

        Ok(Self { path, header, frames_start, trailer })
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn trailer(&self) -> &Trailer {
        &self.trailer
    }

    /// Initial + final bookend Snapshots. The format contract guarantees
    /// both frames are Snapshots (initial bookend is written before any
    /// other frame, final bookend is written immediately before the
    /// trailer); an open record that omits either is malformed and
    /// returns [`RecordError::MissingBookend`].
    pub fn bookends(&self) -> Result<(Snapshot, Snapshot), RecordError> {
        let mut iter = self.frames()?;
        let first = match iter.next() {
            Some(Ok(Frame::Snapshot(s))) => s,
            Some(Ok(_)) | None => return Err(RecordError::MissingBookend),
            Some(Err(e)) => return Err(e),
        };
        let mut last: Option<Snapshot> = None;
        for f in iter {
            if let Frame::Snapshot(s) = f? {
                last = Some(s);
            }
        }
        let last = last.ok_or(RecordError::MissingBookend)?;
        Ok((first, last))
    }

    pub fn events(
        &self,
    ) -> Result<impl Iterator<Item = Result<crate::records::frame::Event, RecordError>>, RecordError>
    {
        let frames = self.frames()?;
        Ok(frames.filter_map(|f| match f {
            Ok(Frame::Event(e)) => Some(Ok(e)),
            Ok(_) => None,
            Err(e) => Some(Err(e)),
        }))
    }

    pub fn dense(
        &self,
    ) -> Result<
        impl Iterator<Item = Result<crate::records::frame::Snapshot, RecordError>>,
        RecordError,
    > {
        let frames = self.frames()?;
        Ok(frames.filter_map(|f| match f {
            Ok(Frame::Snapshot(s)) => Some(Ok(s)),
            Ok(_) => None,
            Err(e) => Some(Err(e)),
        }))
    }

    fn frames(&self) -> Result<FrameIter, RecordError> {
        let mut file = BufReader::new(File::open(&self.path)?);
        file.seek(SeekFrom::Start(self.frames_start))?;
        Ok(FrameIter { file })
    }
}

struct FrameIter {
    file: BufReader<File>,
}

impl Iterator for FrameIter {
    type Item = Result<Frame, RecordError>;
    fn next(&mut self) -> Option<Self::Item> {
        match Frame::read(&mut self.file) {
            Ok(None) => None,
            Ok(Some(Frame::Trailer(_))) => None,
            Ok(Some(f)) => Some(Ok(f)),
            Err(e) => Some(Err(RecordError::Io(e))),
        }
    }
}

#[derive(Debug)]
pub enum RecordError {
    Io(std::io::Error),
    BadMagic,
    UnsupportedFormatVersion(u16),
    BadHeaderUtf8,
    HeaderParse(String),
    MissingTrailer,
    MissingBookend,
}

impl std::fmt::Display for RecordError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "io error: {e}"),
            Self::BadMagic => write!(f, "file does not start with APSR magic"),
            Self::UnsupportedFormatVersion(v) => {
                write!(f, "unsupported format version: {v} (expected {FORMAT_VER})")
            },
            Self::BadHeaderUtf8 => write!(f, "header section is not valid UTF-8"),
            Self::HeaderParse(msg) => write!(f, "header TOML parse: {msg}"),
            Self::MissingTrailer => write!(f, "record has no trailer (truncated or partial)"),
            Self::MissingBookend => {
                write!(f, "record is missing one or both bookend Snapshots (malformed)")
            },
        }
    }
}

impl std::error::Error for RecordError {}

impl From<std::io::Error> for RecordError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}
