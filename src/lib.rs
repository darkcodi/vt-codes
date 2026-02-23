#![no_std]
#![forbid(unsafe_code)]

//! # VT-Codes - Varshamov-Tenengolts Codes
//!
//! This library provides in-place encoding and decoding using q-ary VT (Varshamov-Tenengolts)
//! codes for `no_std` environments. VT codes can correct a single insertion or deletion
//! in transmitted data — unlike traditional error correction codes that only handle bit flips.
//!
//! ## Overview
//!
//! VT codes are a class of error-correcting codes designed to handle synchronization errors.
//! While most codes (like Hamming or Reed-Solomon) correct substitutions (bit flips),
//! VT codes correct *indels* — insertions and deletions of symbols.
//!
//! This library implements q-ary VT codes with q=255, meaning each symbol is a byte.
//! The construction adds redundancy to enable correction of exactly one insertion
//! or deletion in the encoded message.
//!
//! ## Features
//!
//! - `no_std` compatible — works in embedded and bare-metal environments
//! - In-place encoding and decoding — minimal memory overhead
//! - Zero allocations — no heap required
//! - Corrects single insertions or deletions (not bit flips)
//!
//! ## Functions
//!
//! - [`vt_encode_in_place`] — Encode data in-place, returning the codeword length
//! - [`vt_decode_in_place`] — Decode and correct indels in-place, returning message length
//!
//! ## Example
//!
//! ```rust
//! use vt_codes::{vt_encode_in_place, vt_decode_in_place};
//!
//! let mut buf = [0u8; 256];
//! let data = b"Hello world";
//! buf[..data.len()].copy_from_slice(data);
//!
//! // Encode
//! let enc_len = vt_encode_in_place(&mut buf, data.len()).unwrap();
//!
//! // Simulate a deletion
//! buf.copy_within(5.., 4); // Delete byte at position 4
//!
//! // Decode (corrects the deletion)
//! let dec_len = vt_decode_in_place(&mut buf, enc_len - 1).unwrap();
//!
//! assert_eq!(&buf[..dec_len], data);
//! ```

use crate::buffer::Scratch;
use crate::vt::{decode_internal_q255, encode_q_ary_into, find_k_q255, find_smallest_n_checked_q255};

mod error;
mod buffer;
mod vt;
mod utils;

pub use error::Error;

/// Encodes data in-place using q-ary VT codes with q=255 (byte symbols).
///
/// # Algorithm Overview
///
/// VT (Varshamov-Tenengolts) codes add redundancy to enable correction of
/// single insertions or deletions. This implementation uses q=255 (byte symbols)
/// and constructs a codeword of length n from a message of length k, where
/// n ≥ k and the specific relationship allows detection/correction of indels.
///
/// The encoding computes a checksum value α and constructs a codeword where
/// the weighted sum of all symbols (with position weights) equals α modulo q.
///
/// # How it Works
///
/// 1. Find the smallest codeword length n that can accommodate the k-byte message
/// 2. Compute the checksum α from the message symbols and their positions
/// 3. Construct the codeword by embedding the message and adding redundancy
/// 4. The resulting codeword satisfies: Σ(i * c\[i\]) ≡ α (mod q) for all i
///
/// This mathematical property allows the decoder to detect where an insertion
/// or deletion occurred and recover the original data.
///
/// # Buffer Layout
///
/// - **Input:** `buf[0..len]` contains the message to encode
/// - **Output:** `buf[0..n]` contains the encoded codeword
///
/// The buffer is modified in-place. The codeword length n is always ≥ k,
/// with the exact value determined by the message length.
///
/// # Constraints
///
/// - Message length: 1–240 bytes (k)
/// - Codeword length: 6–256 bytes (n)
/// - Buffer must be at least n bytes
/// - Returns `Error::InvalidInputLength` if len is 0 or > 240
/// - Returns `Error::BufferTooSmall` if buffer is too small for the codeword
///
/// # In-Place Operation
///
/// This function modifies the buffer in-place:
/// - Reads the message from `buf[0..len]`
/// - Writes the codeword to `buf[0..n]`
/// - Uses internal scratch space to avoid clobbering input while reading
/// - Returns the codeword length n
///
/// # Example
///
/// ```rust
/// let mut buf = [0u8; 256];
/// buf[..5].copy_from_slice(b"Hello");
///
/// let enc_len = vt_codes::vt_encode_in_place(&mut buf, 5)?;
/// // buf[0..enc_len] now contains the VT codeword
/// # Ok::<(), vt_codes::Error>(())
/// ```
pub fn vt_encode_in_place(
    buf: &mut [u8],
    len: usize,
) -> Result<usize, Error> {
    let mut scratch = Scratch::new();
    if len == 0 {
        return Err(Error::InvalidInputLength);
    }

    // With this construction, maximum message length occurs at n=255.
    let max_k = find_k_q255(255) as usize; // = 240 for this construction
    if len > max_k {
        return Err(Error::InvalidInputLength);
    }

    // For this VT construction (q=255), the smallest n such that k fits ends up with find_k(n)==k.
    let k = len as u8;
    let n = find_smallest_n_checked_q255(k).ok_or(Error::InvalidInputLength)?;
    let n_usize = n as usize;

    if buf.len() < n_usize {
        return Err(Error::BufferTooSmall);
    }

    // Encode into scratch, then copy out (avoid clobbering input while reading it).
    encode_q_ary_into(
        n,
        &buf[..len],
        &mut scratch.cw[..n_usize],
        &mut scratch.alpha,
    );
    buf[..n_usize].copy_from_slice(&scratch.cw[..n_usize]);

    Ok(n_usize)
}

/// Decodes data in-place using q-ary VT codes with q=255 (byte symbols).
///
/// Corrects a single deletion or insertion in the received data.
///
/// # Algorithm Overview
///
/// VT decoding works by testing all possible message lengths and attempting
/// to reconstruct the original message. The decoder tries:
///
/// 1. The received length as-is (no error)
/// 2. The received length minus one (insertion occurred)
/// 3. The received length plus one (deletion occurred)
///
/// For each candidate length n, the decoder attempts to find a valid VT
/// codeword that matches the received data. If successful, it returns the
/// decoded message.
///
/// # How it Works
///
/// - **Deletion correction:** When a symbol was deleted, the received data
///   has length n-1. The decoder finds which symbol is missing and where.
///
/// - **Insertion correction:** When a symbol was inserted, the received data
///   has length n+1. The decoder identifies and removes the extra symbol.
///
/// The decoder validates candidates by checking if the weighted sum
/// property holds: Σ(i * c\[i\]) ≡ α (mod q)
///
/// # Error Handling
///
/// Returns `Error::CorruptedData` if:
/// - The data has more than one insertion/deletion
/// - The data is too corrupted to recover
/// - No valid VT codeword can be found
///
/// # Constraints
///
/// - Received length: 1–256 bytes
/// - Buffer must be at least as large as the received length
/// - Returns `Error::InvalidInputLength` if len is 0 or > 256
///
/// # In-Place Operation
///
/// This function modifies the buffer in-place:
/// - Reads the received data from `buf[0..len]`
/// - Writes the decoded message to `buf[0..k]`
/// - Returns the original message length k
///
/// # Example
///
/// ```rust
/// let mut buf = [0u8; 256];
/// buf[..5].copy_from_slice(b"Hello");
///
/// // Encode
/// let enc_len = vt_codes::vt_encode_in_place(&mut buf, 5)?;
///
/// // Simulate a deletion
/// buf.copy_within(2.., 1); // Delete byte at position 1
///
/// // Decode (corrects the deletion)
/// let dec_len = vt_codes::vt_decode_in_place(&mut buf, enc_len - 1)?;
/// assert_eq!(&buf[..dec_len], b"Hello");
/// # Ok::<(), vt_codes::Error>(())
/// ```
pub fn vt_decode_in_place(
    buf: &mut [u8],
    len: usize,
) -> Result<usize, Error> {
    let mut scratch = Scratch::new();
    if len == 0 || len > 256 || len > buf.len() {
        return Err(Error::InvalidInputLength);
    }

    scratch.rx[..len].copy_from_slice(&buf[..len]);

    let (cands, cand_count): ([u8; 3], usize) = if len == 256 {
        ([255, 0, 0], 1)
    } else {
        let l = len as u8;
        ([l, l.saturating_add(1), l.saturating_sub(1)], 3)
    };

    for i in 0..cand_count {
        let n = cands[i];
        if n < 6 {
            continue;
        }

        let k = find_k_q255(n);
        if k == 0 {
            continue;
        }

        // Split borrows: rx immutably, others mutably (disjoint fields => OK).
        let received = &scratch.rx[..len];
        let cw = &mut scratch.cw;
        let tmp = &mut scratch.tmp;
        let alpha = &mut scratch.alpha;

        if decode_internal_q255(n, received, buf, cw, tmp, alpha) {
            return Ok(k as usize);
        }
    }

    Err(Error::CorruptedData)
}

#[cfg(test)]
mod tests {
    extern crate std;

    use super::*;

    #[test]
    fn round_trip_no_error() {
        let data = b"Hello world";
        let len = data.len();

        let mut buf = [0u8; 256];
        buf[..len].copy_from_slice(data);

        let n = vt_encode_in_place(&mut buf, len).unwrap();
        let k = vt_decode_in_place(&mut buf, n).unwrap();
        assert_eq!(&buf[..k], data);
    }

    #[test]
    fn test_insertions() {
        let data = b"Hello world";
        let len = data.len();

        let mut buf = [0u8; 256];
        buf[..len].copy_from_slice(data);

        let n = vt_encode_in_place(&mut buf, len).unwrap();

        for ins_pos in 0..=n {
            let mut with_insertion = std::vec::Vec::from(&buf[..n]);
            with_insertion.insert(ins_pos, 0x2A);

            let mut decode_buf = [0u8; 256];
            decode_buf[..with_insertion.len()].copy_from_slice(&with_insertion);

            let decoded_len = vt_decode_in_place(&mut decode_buf, with_insertion.len()).unwrap();
            assert_eq!(
                &decode_buf[..decoded_len],
                data,
                "failed at insertion pos {ins_pos}"
            );
        }
    }

    #[test]
    fn test_deletions() {
        let data = b"Hello world";
        let len = data.len();

        let mut buf = [0u8; 256];
        buf[..len].copy_from_slice(data);

        let n = vt_encode_in_place(&mut buf, len).unwrap();

        for del_pos in 0..n {
            let mut with_deletion = std::vec::Vec::from(&buf[..n]);
            with_deletion.remove(del_pos);

            let mut decode_buf = [0u8; 256];
            decode_buf[..with_deletion.len()].copy_from_slice(&with_deletion);

            let decoded_len = vt_decode_in_place(&mut decode_buf, with_deletion.len()).unwrap();
            assert_eq!(
                &decode_buf[..decoded_len],
                data,
                "failed at deletion pos {del_pos}"
            );
        }
    }
}
