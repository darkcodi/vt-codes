use std::collections::HashMap;
use std::collections::HashSet;

// --- Error type ---

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    InvalidInputLength,
    BufferTooSmall,
    CorruptedData,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidInputLength => write!(f, "Invalid input length"),
            Error::BufferTooSmall => write!(f, "Buffer too small"),
            Error::CorruptedData => write!(f, "Corrupted data, cannot recover"),
        }
    }
}

impl std::error::Error for Error {}

// --- Public API ---

/// Encodes data in-place using Varshamov-Tenengolts codes.
///
/// # Arguments
/// * `buf` - Buffer containing input data at `[0..len]` and output codeword at `[0..n]`
/// * `len` - Length of input message in bytes
///
/// # Returns
/// * `Ok(n)` - Length of the encoded codeword
/// * `Err(Error::BufferTooSmall)` - If buffer is too small to hold the codeword
/// * `Err(Error::InvalidInputLength)` - If input length is invalid
///
/// # Examples
/// ```
/// use vt_ecc::{vt_encode_in_place, vt_decode_in_place, Error};
/// let mut buf = vec![0u8; 32];
/// buf[..11].copy_from_slice(b"Hello world");
/// let n = vt_encode_in_place(&mut buf, 11).unwrap();
/// assert_eq!(n, 20); // codeword length for k=11 with q=255
/// // Decode back
/// let k = vt_decode_in_place(&mut buf, n).unwrap();
/// assert_eq!(&buf[..k], b"Hello world");
/// # Ok::<(), Error>(())
/// ```
pub fn vt_encode_in_place(buf: &mut [u8], len: usize) -> Result<usize, Error> {
    if len == 0 {
        return Err(Error::InvalidInputLength);
    }

    let max_k = find_k(255, 255) as usize;
    if len > max_k {
        return Err(Error::InvalidInputLength);
    }

    let k = len as u8;
    let n = find_smallest_n_checked(k, 255).ok_or(Error::InvalidInputLength)?;
    let n_usize = n as usize;

    if buf.len() < n_usize {
        return Err(Error::BufferTooSmall);
    }

    let t = ceil_log2(n);
    let (systematic_positions, table_1_l, table_1_r, _table_1_rev, table_2, _table_2_rev) = generate_tables(n, t);

    let encoded = encode_q_ary(n, k, &systematic_positions, &table_1_l, &table_1_r, &table_2, &buf[..len]);
    buf[..n_usize].copy_from_slice(&encoded);

    Ok(n_usize)
}

/// Decodes data in-place using Varshamov-Tenengolts codes.
///
/// Corrects a single deletion or insertion in the received data.
///
/// # Arguments
/// * `buf` - Buffer containing received data at `[0..len]`, will contain decoded message at `[0..k]`
/// * `len` - Length of received data (may have 1 deletion/insertion)
///
/// # Returns
/// * `Ok(k)` - Length of the decoded message
/// * `Err(Error::InvalidInputLength)` - If input length is invalid
/// * `Err(Error::CorruptedData)` - If data is too corrupted to recover
///
/// # Examples
/// ```
/// use vt_ecc::{vt_encode_in_place, vt_decode_in_place, Error};
/// let mut buf = vec![0u8; 32];
/// buf[..11].copy_from_slice(b"Hello world");
/// let n = vt_encode_in_place(&mut buf, 11).unwrap();
/// // After transmission with potential errors:
/// let k = vt_decode_in_place(&mut buf, n).unwrap();
/// assert_eq!(&buf[..k], b"Hello world");
/// # Ok::<(), Error>(())
/// ```
pub fn vt_decode_in_place(buf: &mut [u8], len: usize) -> Result<usize, Error> {
    if len == 0 || len > 256 || len > buf.len() {
        return Err(Error::InvalidInputLength);
    }

    // Try different possible n values: len, len+1, len-1
    let possible_n: Vec<u8> = vec![
        len as u8,
        (len as u8).saturating_add(1),
        (len as u8).saturating_sub(1),
    ];

    for n in possible_n {
        if n < 2 {
            continue;
        }

        let k = find_k(n, 255);
        if k == 0 {
            continue;
        }

        let t = ceil_log2(n);
        let (systematic_positions, table_1_l, table_1_r, table_1_rev, table_2, table_2_rev) = generate_tables(n, t);

        let k_usize = k as usize;

        if let Some(decoded) = decode_internal(
            n,
            &buf[..len.min(buf.len())],
            &systematic_positions,
            &table_1_l,
            &table_1_r,
            &table_1_rev,
            &table_2,
            &table_2_rev,
        ) {
            buf[..k_usize].copy_from_slice(&decoded);
            return Ok(k_usize);
        }
    }

    Err(Error::CorruptedData)
}

// --- Alpha & syndrome ---

fn convert_y_to_alpha(y: &[u8]) -> Vec<u8> {
    (0..y.len() - 1)
        .map(|i| if y[i + 1] >= y[i] { 1 } else { 0 })
        .collect()
}

fn compute_syndrome_binary(m: u8, a: u8, y: &[u8]) -> u8 {
    let mut s: i64 = 0;
    for (i, &yi) in y.iter().enumerate() {
        s += (i as i64 + 1) * yi as i64;
    }
    mod_floor(a as i64 - s, m as i64) as u8
}

fn compute_syndrome_q_ary(m: u8, a: u8, b: u8, y: &[u8]) -> (u8, u8) {
    let alpha = convert_y_to_alpha(y);
    let s1 = compute_syndrome_binary(m, a, &alpha);
    let sum: i64 = y.iter().map(|&v| v as i64).sum();
    // q=255, so modulus is 256
    let s2 = mod_floor(b as i64 - sum, 256) as u8;
    (s1, s2)
}

// --- find_k and find_smallest_n ---

fn find_k(n: u8, q: u8) -> u8 {
    if q == 1 {
        n - ceil_log2(n + 1)
    } else {
        let t = ceil_log2(n);
        if q == 2 {
            if n < 7 {
                return 0;
            }
            let alphabet_size = (q as i64 + 1) as f64;
            if power_of_two(n - 1) {
                ((n as i64 - 3 * t as i64 + 3) as f64 * alphabet_size.log2()).floor() as i64 as u8 + 2 * (t - 4) + 1
            } else {
                ((n as i64 - 3 * t as i64 + 3) as f64 * alphabet_size.log2()).floor() as i64 as u8 + 2 * (t - 3)
            }
        } else {
            // For q > 2 (specifically q=255), we work with bytes directly
            if n < 6 {
                return 0;
            }
            // Each position stores 1 byte directly, no conversion needed
            let base = 0i64.max(n as i64 - 3 * t as i64 + 3);
            // Each iteration of step 2 stores 1 byte
            let bytes_per_tuple = 1i64;
            // Step 2b stores 1 byte
            let bytes_single = 1i64;
            if power_of_two(n - 1) {
                (base + bytes_per_tuple * 0i64.max(t as i64 - 4) + 2 * bytes_single) as u8
            } else {
                (base + bytes_per_tuple * 0i64.max(t as i64 - 3) + bytes_single) as u8
            }
        }
    }
}

fn find_smallest_n_checked(k: u8, q: u8) -> Option<u8> {
    let log_q_plus_1 = if q == 255 { 8 } else { ceil_log2(q.wrapping_add(1)) };

    let mut n: u16 = (k as u16) / (log_q_plus_1 as u16).max(1);
    while n <= 255 {
        let nu8 = n as u8;
        if find_k(nu8, q) >= k {
            return Some(nu8);
        }
        n += 1;
    }
    None
}

// --- Error correction ---

fn correct_binary_indel(n: u8, m: u8, a: u8, y: &[u8]) -> Option<Vec<u8>> {
    let nu = n as usize;
    let s = compute_syndrome_binary(m, a, y);
    let w: u8 = y.iter().sum();

    if y.len() == nu - 1 {
        // deletion
        let mut y_decoded = vec![0u8; nu];
        if s == 0 {
            y_decoded[..nu - 1].copy_from_slice(y);
        } else if s <= w {
            let mut num_ones_seen: u8 = 0;
            let mut found = false;
            for i in (0..nu - 1).rev() {
                if y[i] == 1 {
                    num_ones_seen += 1;
                    if num_ones_seen == s {
                        y_decoded[..i].copy_from_slice(&y[..i]);
                        // insert 0 at position i
                        y_decoded[i + 1..].copy_from_slice(&y[i..]);
                        found = true;
                        break;
                    }
                }
            }
            if !found {
                return None;
            }
        } else {
            // 1 deleted, s-w-1 = number of 0s to left
            let target = s - w - 1;
            if target == 0 {
                y_decoded[0] = 1;
                y_decoded[1..].copy_from_slice(y);
            } else {
                let mut num_zeros_seen: u8 = 0;
                let mut found = false;
                for i in 0..nu - 1 {
                    if y[i] == 0 {
                        num_zeros_seen += 1;
                        if num_zeros_seen == target {
                            y_decoded[..=i].copy_from_slice(&y[..=i]);
                            y_decoded[i + 1] = 1;
                            y_decoded[i + 2..].copy_from_slice(&y[i + 1..]);
                            found = true;
                            break;
                        }
                    }
                }
                if !found {
                    return None;
                }
            }
        }
        Some(y_decoded)
    } else {
        // insertion (y.len() == n+1)
        let mut y_decoded = vec![0u8; nu];
        let s_val = mod_floor(m as i64 - n as i64 - 1, m as i64) as u8;
        let m_minus_w = mod_floor(m as i64 - w as i64, m as i64) as u8;
        if s == s_val || s == 0 {
            // last entry inserted (s == m-n-1 mod m, which covers s==0 when m-n-1==0)
            // but Python checks s == m-n-1 or s == 0 separately
            y_decoded.copy_from_slice(&y[..nu]);
        } else if s == m_minus_w {
            // first entry inserted
            y_decoded.copy_from_slice(&y[1..=nu]);
        } else if s > m_minus_w {
            // 0 was inserted, m-s 1's to the right
            let target = mod_floor(m as i64 - s as i64, m as i64) as u8;
            let mut num_ones_seen: u8 = 0;
            let mut found = false;
            for i in (2..=nu).rev() {
                if y[i] == 1 {
                    num_ones_seen += 1;
                    if num_ones_seen == target {
                        if y[i - 1] == 0 {
                            y_decoded[..i - 1].copy_from_slice(&y[..i - 1]);
                            y_decoded[i - 1..].copy_from_slice(&y[i..]);
                            found = true;
                        }
                        break;
                    }
                }
            }
            if !found {
                return None;
            }
        } else {
            // 1 was inserted, m-w-s 0's to the left
            let target = mod_floor(m as i64 - w as i64 - s as i64, m as i64) as u8;
            let mut num_zeros_seen: u8 = 0;
            let mut found = false;
            for i in 0..nu {
                if y[i] == 0 {
                    num_zeros_seen += 1;
                    if num_zeros_seen == target {
                        if y[i + 1] == 1 {
                            y_decoded[..=i].copy_from_slice(&y[..=i]);
                            y_decoded[i + 1..].copy_from_slice(&y[i + 2..]);
                            found = true;
                        }
                        break;
                    }
                }
            }
            if !found {
                return None;
            }
        }
        Some(y_decoded)
    }
}

fn correct_q_ary_indel(n: u8, m: u8, a: u8, b: u8, y: &[u8]) -> Option<Vec<u8>> {
    let alpha = convert_y_to_alpha(y);
    let alpha_corrected = correct_binary_indel(n - 1, m, a, &alpha)?;
    if compute_syndrome_binary(m, a, &alpha_corrected) != 0 {
        return None;
    }

    let nu = n as usize;
    let mut y_decoded = vec![0u8; nu];

    if alpha.len() == nu - 2 {
        // deletion
        let sum: i64 = y.iter().map(|&v| v as i64).sum();
        // q=255, so modulus is 256
        let error_symbol = mod_floor(b as i64 - sum, 256) as u8;

        // find diff_pos
        let diff_pos = if alpha == alpha_corrected[..alpha.len()] {
            nu - 2
        } else {
            let mut dp = 0;
            for i in 0..alpha.len() {
                if alpha[i] != alpha_corrected[i] {
                    dp = i;
                    break;
                }
            }
            dp
        };

        // find del_pos
        let mut del_pos_found = false;
        for del_pos in (0..=diff_pos + 1).rev() {
            let ok = if del_pos == 0 {
                alpha_corrected[0] == (if y[0] >= error_symbol { 1 } else { 0 })
            } else if del_pos == nu - 1 {
                alpha_corrected[nu - 2] == (if error_symbol >= y[nu - 2] { 1 } else { 0 })
            } else {
                (alpha_corrected[del_pos - 1] == (if error_symbol >= y[del_pos - 1] { 1 } else { 0 }))
                    && (alpha_corrected[del_pos] == (if y[del_pos] >= error_symbol { 1 } else { 0 }))
            };
            if ok {
                y_decoded[..del_pos].copy_from_slice(&y[..del_pos]);
                y_decoded[del_pos] = error_symbol;
                y_decoded[del_pos + 1..].copy_from_slice(&y[del_pos..]);
                del_pos_found = true;
                break;
            }
        }
        if !del_pos_found {
            return None;
        }
    } else {
        // insertion
        let sum: i64 = y.iter().map(|&v| v as i64).sum();
        // q=255, so modulus is 256
        let error_symbol = mod_floor(sum - b as i64, 256) as u8;

        // find diff_pos
        let diff_pos = if alpha[..alpha_corrected.len()] == alpha_corrected[..] {
            nu - 1
        } else {
            let mut dp = 0;
            for i in 0..alpha_corrected.len() {
                if alpha[i] != alpha_corrected[i] {
                    dp = i;
                    break;
                }
            }
            dp
        };

        // find ins_pos
        let mut ins_pos_found = false;
        for ins_pos in (0..=diff_pos + 1).rev() {
            let ok = if ins_pos == 0 || ins_pos == nu {
                y[ins_pos] == error_symbol
            } else {
                (y[ins_pos] == error_symbol)
                    && (alpha_corrected[ins_pos - 1]
                    == (if y[ins_pos + 1] >= y[ins_pos - 1] { 1 } else { 0 }))
            };
            if ok {
                y_decoded[..ins_pos].copy_from_slice(&y[..ins_pos]);
                y_decoded[ins_pos..].copy_from_slice(&y[ins_pos + 1..]);
                ins_pos_found = true;
                break;
            }
        }
        if !ins_pos_found {
            return None;
        }
    }

    if compute_syndrome_q_ary(m, a, b, &y_decoded) == (0, 0) {
        Some(y_decoded)
    } else {
        None
    }
}

// --- Encoding/decoding core functions ---

fn encode_q_ary(n: u8, k: u8, systematic_positions: &[u8], table_1_l: &[u8], table_1_r: &[u8], table_2: &[u8], x: &[u8]) -> Vec<u8> {
    let nu = n as usize;
    let t = ceil_log2(n);
    let mut y = vec![0u8; nu];

    // step 1: encode bytes in non-dyadic positions
    let step_1_num_bytes = 0i64.max(n as i64 - 3 * t as i64 + 3) as usize;
    if step_1_num_bytes > 0 {
        for (i, &pos) in systematic_positions.iter().enumerate() {
            if i < step_1_num_bytes {
                y[pos as usize] = x[i];
            }
        }
    }

    // step 2: encode bytes in near-dyadic positions
    let mut bytes_done = step_1_num_bytes;
    for j in 3..t {
        let pj = 1u32 << j;
        if pj as u8 == n - 1 {
            // special case: store in y[2^j - 1]
            y[(pj - 1) as usize] = x[bytes_done].wrapping_add(1);
            bytes_done += 1;
            break;
        }
        let table_1_index = x[bytes_done] as usize;
        y[(pj - 1) as usize] = table_1_r[table_1_index];
        y[(pj + 1) as usize] = table_1_l[table_1_index];
        bytes_done += 1;
    }

    // set y[3] and y[5]
    y[3] = 255;
    let table_2_index = x[bytes_done] as usize;
    y[5] = table_2[table_2_index];
    bytes_done += 1;

    // step 3: set alpha at positions except dyadic
    let mut alpha = convert_y_to_alpha(&y);
    for j in 2..t {
        let pj = 1u32 << j;
        if pj as u8 == n - 1 {
            break;
        }
        alpha[(pj + 1 - 1) as usize] = if y[(pj + 1) as usize] >= y[(pj - 1) as usize] {
            1
        } else {
            0
        };
    }
    alpha[2] = 1;

    // step 4: set alpha at dyadic positions using VT conditions
    for j in 0..t {
        alpha[(1usize << j) - 1] = 0;
    }
    let m = n;
    let mut syndrome = compute_syndrome_binary(m, 0, &alpha);
    if syndrome != 0 {
        for j in (0..t).rev() {
            let pos = 1u32 << j;
            if syndrome >= pos as u8 {
                alpha[(pos - 1) as usize] = 1;
                syndrome -= pos as u8;
                if syndrome == 0 {
                    break;
                }
            }
        }
    }

    // step 5: set symbols of y at dyadic positions except 1 and 2
    for j in 2..t {
        let pos = 1usize << j;
        if alpha[pos - 1] == 0 {
            y[pos] = y[pos - 1].wrapping_sub(1);
        } else {
            y[pos] = y[pos - 1];
        }
    }

    // step 6: set positions 0, 1, 2
    let sum: i64 = y[3..].iter().map(|&v| v as i64).sum();
    let w = mod_floor(0i64 - sum, 256) as u8;
    let (val_x, val_y, val_z) = if w == 1 {
        (0, 2, 255)
    } else if w == 2 {
        (1, 2, 255)
    } else {
        (0, 1, mod_floor(w as i64 - 1, 256) as u8)
    };
    if alpha[0] == 0 && alpha[1] == 0 {
        y[0] = val_z;
        y[1] = val_y;
        y[2] = val_x;
    } else if alpha[0] == 0 && alpha[1] == 1 {
        y[0] = val_z;
        y[1] = val_x;
        y[2] = val_y;
    } else if alpha[0] == 1 && alpha[1] == 0 {
        y[0] = val_x;
        y[1] = val_z;
        y[2] = val_y;
    } else {
        y[0] = val_x;
        y[1] = val_y;
        y[2] = val_z;
    }

    y
}

fn decode_codeword_q_ary(n: u8, k: u8, systematic_positions: &[u8], table_1_rev: &HashMap<(u8, u8), u16>, table_2_rev: &HashMap<u8, u8>, y: &[u8]) -> Option<Vec<u8>> {
    let t = ceil_log2(n);
    let mut x = vec![0u8; k as usize];

    // step 1: decode bytes from non-dyadic positions
    let step_1_num_bytes = 0i64.max(n as i64 - 3 * t as i64 + 3) as usize;
    if step_1_num_bytes > 0 {
        for (i, &pos) in systematic_positions.iter().enumerate() {
            if i < step_1_num_bytes {
                x[i] = y[pos as usize];
            }
        }
    }

    // step 2: decode bytes from near-dyadic positions
    let mut bytes_done = step_1_num_bytes;
    for j in 3..t {
        let pj = 1u32 << j;
        if pj as u8 == n - 1 {
            if y[(pj - 1) as usize] == 0 {
                return None;
            }
            x[bytes_done] = y[(pj - 1) as usize].wrapping_sub(1);
            bytes_done += 1;
            break;
        }
        let r = y[(pj - 1) as usize];
        let l = y[(pj + 1) as usize];
        if let Some(&idx) = table_1_rev.get(&(r, l)) {
            x[bytes_done] = idx as u8;
        } else {
            return None;
        }
        bytes_done += 1;
    }

    // step 2b: y[5]
    if y[3] != 255 {
        return None;
    }
    if let Some(&idx) = table_2_rev.get(&y[5]) {
        x[bytes_done] = idx;
    } else {
        return None;
    }
    bytes_done += 1;

    Some(x)
}

fn generate_tables(n: u8, t: u8) -> (Vec<u8>, Vec<u8>, Vec<u8>, HashMap<(u8, u8), u16>, Vec<u8>, HashMap<u8, u8>) {
    // table 1
    // For q=255 we map a single byte (0..255) -> a pair (r,l), so size must be 256.
    let table_1_size = 256usize;
    let mut table_1_l = vec![0u8; table_1_size];
    let mut table_1_r = vec![0u8; table_1_size];
    let mut pos = 0;
    for r in 0..=255 {
        if pos == table_1_size {
            break;
        }
        if r == 0 {
            continue;
        }
        for l in 0..=255 {
            if pos == table_1_size {
                break;
            }
            if l == r - 1 {
                continue;
            }
            table_1_l[pos] = l;
            table_1_r[pos] = r;
            pos += 1;
        }
    }

    let mut table_1_rev = HashMap::new();
    for i in 0..table_1_size {
        table_1_rev.insert((table_1_r[i], table_1_l[i]), i as u16);
    }

    // table 2
    // Also a direct byte->symbol map: size must be 256 and must be bijective.
    let table_2_size = 256usize;
    let mut table_2 = vec![0u8; table_2_size];
    for i in 0..table_2_size {
        table_2[i] = i as u8;
    }
    let mut table_2_rev = HashMap::new();
    for i in 0..table_2_size {
        table_2_rev.insert(table_2[i], i as u8);
    }

    // systematic positions
    let mut non_sys: Vec<u8> = vec![0, 1, 2, 3, 4, 5];
    for j in 3..t {
        let pj = 1u32 << j;
        non_sys.push((pj - 1) as u8);
        non_sys.push(pj as u8);
        non_sys.push((pj + 1) as u8);
    }
    let non_sys_set: HashSet<u8> = non_sys.into_iter().collect();
    let systematic_positions: Vec<u8> = (0..n).filter(|x| !non_sys_set.contains(x)).collect();

    (systematic_positions, table_1_l, table_1_r, table_1_rev, table_2, table_2_rev)
}

fn is_codeword(n: u8, y: &[u8]) -> bool {
    if y.len() != n as usize {
        return false;
    }
    compute_syndrome_q_ary(n, 0, 0, y) == (0, 0)
}

fn decode_internal(
    n: u8,
    y: &[u8],
    systematic_positions: &[u8],
    table_1_l: &[u8],
    table_1_r: &[u8],
    table_1_rev: &HashMap<(u8, u8), u16>,
    table_2: &[u8],
    table_2_rev: &HashMap<u8, u8>,
) -> Option<Vec<u8>> {
    let n_y = y.len() as i64;
    if n_y < n as i64 - 1 || n_y > n as i64 + 1 {
        return None;
    }

    let corrected = if n_y != n as i64 {
        correct_q_ary_indel(n, n, 0, 0, y)?
    } else {
        y.to_vec()
    };

    if !is_codeword(n, &corrected) {
        return None;
    }

    let k = find_k(n, 255);
    let decoded = decode_codeword_q_ary(n, k, systematic_positions, table_1_rev, table_2_rev, &corrected)?;

    // CRITICAL: Ensure the candidate is a codeword from *our* encoder construction,
    // not merely something that satisfies the VT syndromes for some other n.
    let reencoded = encode_q_ary(
        n,
        k,
        systematic_positions,
        table_1_l,
        table_1_r,
        table_2,
        &decoded,
    );
    if reencoded == corrected { Some(decoded) } else { None }
}

// --- Utils ---

/// Python-style non-negative modulo: ((a % m) + m) % m
fn mod_floor(a: i64, m: i64) -> i64 {
    ((a % m) + m) % m
}

fn ceil_log2(x: u8) -> u8 {
    (x as f64).log2().ceil() as u8
}

fn power_of_two(num: u8) -> bool {
    num > 0 && (num & (num - 1)) == 0
}

#[cfg(test)]
mod tests {
    extern crate alloc;
    extern crate std;

    use super::*;
    use alloc::vec;

    #[test]
    fn test_insertions() {
        // arrange
        let data = b"Hello world";
        let len = data.len();
        let mut buf = vec![0u8; 256];
        buf[..len].copy_from_slice(data);

        // act
        let n = vt_encode_in_place(&mut buf, len).unwrap();
        std::println!("Encoded frame ({} bytes): {:02X?}", n, &buf[..n]);

        // assert: try inserting a byte at every position and verify it decodes correctly (removing the inserted byte).
        for ins_pos in 0..=n {
            let mut with_insertion = buf[..n].to_vec();
            with_insertion.insert(ins_pos, 0x2A);
            std::println!("With insertion at position {}: {:02X?}", ins_pos, &with_insertion);

            let mut decode_buf = vec![0u8; 256];
            decode_buf[..with_insertion.len()].copy_from_slice(&with_insertion);
            let decoded_len = vt_decode_in_place(&mut decode_buf, with_insertion.len()).unwrap();

            assert_eq!(
                &decode_buf[..decoded_len],
                data,
                "Failed to decode with insertion at position {ins_pos}"
            );
        }
    }

    #[test]
    fn test_deletions() {
        // arrange
        let data = b"Hello world";
        let len = data.len();
        let mut buf = vec![0u8; 256];
        buf[..len].copy_from_slice(data);

        // act
        let n = vt_encode_in_place(&mut buf, len).unwrap();
        std::println!("Encoded frame ({} bytes): {:02X?}", n, &buf[..n]);

        // assert: try deleting a byte at every position and verify it reports a valid deletion position.
        for del_pos in 0..n {
            let mut with_deletion = buf[..n].to_vec();
            with_deletion.remove(del_pos);
            std::println!("With deletion at position {}: {:02X?}", del_pos, &with_deletion);

            let mut decode_buf = vec![0u8; 256];
            decode_buf[..with_deletion.len()].copy_from_slice(&with_deletion);
            let decoded_len = vt_decode_in_place(&mut decode_buf, with_deletion.len()).unwrap();

            assert_eq!(
                &decode_buf[..decoded_len],
                data,
                "Failed to decode with deletion at position {del_pos}"
            );
        }
    }
}
