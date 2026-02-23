//! Error types for VT code encoding and decoding.
//!
//! This module defines the error conditions that can occur during
//! VT code operations.

use core::fmt;

/// Errors that can occur during VT code encoding or decoding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// The input length is invalid for the operation.
    ///
    /// For encoding: occurs when the message length is 0 or exceeds 240 bytes.
    /// For decoding: occurs when the received length is 0 or exceeds 256 bytes.
    InvalidInputLength,

    /// The provided buffer is too small for the operation.
    ///
    /// For encoding: the buffer must be at least n bytes where n is the
    /// codeword length (determined by the message length).
    /// For decoding: the buffer must be at least as large as the received length.
    BufferTooSmall,

    /// The received data is corrupted and cannot be recovered.
    ///
    /// This occurs when the decoder cannot find a valid VT codeword,
    /// typically due to:
    /// - More than one insertion or deletion in the data
    /// - Bit flips or other corruption not covered by VT codes
    /// - Invalid or malicious data
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