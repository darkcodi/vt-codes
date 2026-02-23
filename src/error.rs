use core::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    InvalidInputLength,
    BufferTooSmall,
    CorruptedData,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidInputLength => write!(f, "Invalid input length"),
            Error::BufferTooSmall => write!(f, "Buffer too small"),
            Error::CorruptedData => write!(f, "Corrupted data, cannot recover"),
        }
    }
}