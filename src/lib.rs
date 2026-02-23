#![no_std]
#![forbid(unsafe_code)]

use crate::buffer::Scratch;
use crate::error::Error;
use crate::vt::{decode_internal_q255, encode_q_ary_into, find_k_q255, find_smallest_n_checked_q255};

mod error;
mod buffer;
mod vt;
mod utils;

/// Encodes data in-place using q-ary VT codes with q=255 (byte symbols).
///
/// Buffer layout:
/// - input message: `buf[0..len]`
/// - output codeword: `buf[0..n]`
///
/// No heap allocation. Uses `scratch` only.
pub fn vt_encode_in_place(
    buf: &mut [u8],
    len: usize,
    scratch: &mut Scratch,
) -> Result<usize, Error> {
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
/// - `buf[0..len]`: received (length may be n-1, n, or n+1)
/// - on success, writes decoded message into `buf[0..k]`
///
/// No heap allocation. Uses `scratch` only.
pub fn vt_decode_in_place(
    buf: &mut [u8],
    len: usize,
    scratch: &mut Scratch,
) -> Result<usize, Error> {
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

        let mut scratch = Scratch::new();

        let n = vt_encode_in_place(&mut buf, len, &mut scratch).unwrap();
        let k = vt_decode_in_place(&mut buf, n, &mut scratch).unwrap();
        assert_eq!(&buf[..k], data);
    }

    #[test]
    fn test_insertions() {
        let data = b"Hello world";
        let len = data.len();

        let mut buf = [0u8; 256];
        buf[..len].copy_from_slice(data);

        let mut scratch = Scratch::new();
        let n = vt_encode_in_place(&mut buf, len, &mut scratch).unwrap();

        for ins_pos in 0..=n {
            let mut with_insertion = std::vec::Vec::from(&buf[..n]);
            with_insertion.insert(ins_pos, 0x2A);

            let mut decode_buf = [0u8; 256];
            decode_buf[..with_insertion.len()].copy_from_slice(&with_insertion);

            let mut scratch2 = Scratch::new();
            let decoded_len =
                vt_decode_in_place(&mut decode_buf, with_insertion.len(), &mut scratch2).unwrap();
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

        let mut scratch = Scratch::new();
        let n = vt_encode_in_place(&mut buf, len, &mut scratch).unwrap();

        for del_pos in 0..n {
            let mut with_deletion = std::vec::Vec::from(&buf[..n]);
            with_deletion.remove(del_pos);

            let mut decode_buf = [0u8; 256];
            decode_buf[..with_deletion.len()].copy_from_slice(&with_deletion);

            let mut scratch2 = Scratch::new();
            let decoded_len =
                vt_decode_in_place(&mut decode_buf, with_deletion.len(), &mut scratch2).unwrap();
            assert_eq!(
                &decode_buf[..decoded_len],
                data,
                "failed at deletion pos {del_pos}"
            );
        }
    }
}
