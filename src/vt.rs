use crate::utils::{ceil_log2, mod_floor, power_of_two};

// ===== Core (q=255) =====

pub(crate) fn decode_internal_q255(
    n: u8,
    received: &[u8],
    out_buf: &mut [u8],
    cw: &mut [u8; 256],
    tmp: &mut [u8; 256],
    alpha: &mut [u8; 256],
) -> bool {
    let nu = n as usize;

    if received.len() + 1 < nu || received.len() > nu + 1 {
        return false;
    }

    // Build corrected codeword in cw (length n).
    if received.len() != nu {
        if correct_q_ary_indel_into(
            n,
            n,
            0,
            0,
            received,
            &mut cw[..nu],
            &mut alpha[..],
            &mut tmp[..],
        )
            .is_none()
        {
            return false;
        }
    } else {
        cw[..nu].copy_from_slice(received);
    }

    if !is_codeword_q255(n, &cw[..nu], &mut alpha[..]) {
        return false;
    }

    let k = find_k_q255(n) as usize;
    if out_buf.len() < k {
        return false;
    }

    if decode_codeword_q255_into(n, &cw[..nu], &mut out_buf[..k]).is_none() {
        return false;
    }

    // Re-encode and compare.
    encode_q_ary_into(n, &out_buf[..k], &mut tmp[..nu], &mut alpha[..]);
    tmp[..nu] == cw[..nu]
}

// ===== Alpha & syndromes =====

fn convert_y_to_alpha(y: &[u8], out_alpha: &mut [u8]) -> usize {
    if y.len() < 2 {
        return 0;
    }
    let alen = y.len() - 1;
    for i in 0..alen {
        out_alpha[i] = if y[i + 1] >= y[i] { 1 } else { 0 };
    }
    alen
}

fn compute_syndrome_binary(m: u8, a: u8, bits: &[u8]) -> u8 {
    let mut s: i64 = 0;
    for (i, &bi) in bits.iter().enumerate() {
        s += (i as i64 + 1) * (bi as i64);
    }
    mod_floor(a as i64 - s, m as i64) as u8
}

fn compute_syndrome_q_ary(m: u8, a: u8, b: u8, y: &[u8], alpha_buf: &mut [u8]) -> (u8, u8) {
    let alen = convert_y_to_alpha(y, alpha_buf);
    let s1 = compute_syndrome_binary(m, a, &alpha_buf[..alen]);
    let sum: i64 = y.iter().map(|&v| v as i64).sum();
    let s2 = mod_floor(b as i64 - sum, 256) as u8; // q=255 => modulus 256
    (s1, s2)
}

fn is_codeword_q255(n: u8, y: &[u8], alpha_buf: &mut [u8]) -> bool {
    y.len() == n as usize && compute_syndrome_q_ary(n, 0, 0, y, alpha_buf) == (0, 0)
}

// ===== find_k / find_smallest_n (q=255 only) =====

pub(crate) fn find_k_q255(n: u8) -> u8 {
    if n < 6 {
        return 0;
    }
    let t = ceil_log2(n);
    let base = core::cmp::max(0i64, n as i64 - 3 * t as i64 + 3);
    if power_of_two(n - 1) {
        (base + core::cmp::max(0i64, t as i64 - 4) + 2) as u8
    } else {
        (base + core::cmp::max(0i64, t as i64 - 3) + 1) as u8
    }
}

pub(crate) fn find_smallest_n_checked_q255(k: u8) -> Option<u8> {
    // Minimum n for this construction to be meaningful is 6.
    let mut n = core::cmp::max(6u16, k as u16);
    while n <= 255 {
        let nu8 = n as u8;
        if find_k_q255(nu8) >= k {
            return Some(nu8);
        }
        n += 1;
    }
    None
}

// ===== Table 1 (no storage) =====
//
// The original generate_tables() built a 256-entry mapping idx -> (r,l)
// but due to the early stop at 256 entries, it degenerates to:
//   idx 0..254: (r=1, l=idx+1)
//   idx 255:    (r=2, l=0)
// and the reverse mapping is the inverse of that.

fn table1_pair(idx: u8) -> (u8, u8) {
    if idx < 255 {
        (1, idx.wrapping_add(1)) // (r,l)
    } else {
        (2, 0)
    }
}

fn table1_rev(r: u8, l: u8) -> Option<u8> {
    if r == 1 && l != 0 {
        Some(l - 1)
    } else if r == 2 && l == 0 {
        Some(255)
    } else {
        None
    }
}

// ===== Systematic position predicate (no HashSet / no list) =====

fn is_non_systematic(pos: u8, t: u8) -> bool {
    if pos <= 5 {
        return true;
    }
    // near-dyadic triples for j in 3..t: (2^j-1, 2^j, 2^j+1)
    for j in 3..t {
        let pj = (1u16 << (j as u16)) as u16;
        let pjm1 = pj.wrapping_sub(1) as u8;
        let pju8 = pj as u8;
        let pjp1 = pj.wrapping_add(1) as u8;
        if pos == pjm1 || pos == pju8 || pos == pjp1 {
            return true;
        }
    }
    false
}

// ===== Encoding/decoding (q=255) =====

pub(crate) fn encode_q_ary_into(n: u8, x: &[u8], out_y: &mut [u8], alpha_buf: &mut [u8]) {
    let nu = n as usize;
    debug_assert!(out_y.len() >= nu);

    for v in &mut out_y[..nu] {
        *v = 0;
    }

    let t = ceil_log2(n);

    // step 1: encode bytes in systematic positions (all except non-systematic)
    let step_1_num_bytes = core::cmp::max(0i64, n as i64 - 3 * t as i64 + 3) as usize;
    if step_1_num_bytes > 0 {
        let mut xi = 0usize;
        for pos in 0..n {
            if !is_non_systematic(pos, t) {
                out_y[pos as usize] = x[xi];
                xi += 1;
                if xi == step_1_num_bytes {
                    break;
                }
            }
        }
    }

    // step 2: encode bytes in near-dyadic positions
    let mut bytes_done = step_1_num_bytes;
    for j in 3..t {
        let pj = (1u16 << (j as u16)) as u8;

        if pj == n - 1 {
            out_y[(pj as usize).wrapping_sub(1)] = x[bytes_done].wrapping_add(1);
            bytes_done += 1;
            break;
        }

        let (r, l) = table1_pair(x[bytes_done]);
        out_y[(pj as usize).wrapping_sub(1)] = r; // y[2^j - 1]
        out_y[(pj as usize).wrapping_add(1)] = l; // y[2^j + 1]
        bytes_done += 1;
    }

    // set y[3] and y[5]
    out_y[3] = 255;
    out_y[5] = x[bytes_done]; // table_2 is identity

    // step 3: set alpha at positions except dyadic
    let alen = convert_y_to_alpha(&out_y[..nu], alpha_buf);
    let alpha = &mut alpha_buf[..alen];

    for j in 2..t {
        let pj = (1u16 << (j as u16)) as u8;
        if pj == n - 1 {
            break;
        }
        alpha[pj as usize] = if out_y[(pj as usize) + 1] >= out_y[(pj as usize) - 1] {
            1
        } else {
            0
        };
    }
    alpha[2] = 1;

    // step 4: set alpha at dyadic positions using VT condition
    for j in 0..t {
        let p = (1usize << (j as usize)) - 1;
        alpha[p] = 0;
    }
    let m = n;
    let mut syndrome = compute_syndrome_binary(m, 0, alpha);
    if syndrome != 0 {
        for j in (0..t).rev() {
            let pos = (1u16 << (j as u16)) as u8;
            if syndrome >= pos {
                alpha[(pos as usize) - 1] = 1;
                syndrome -= pos;
                if syndrome == 0 {
                    break;
                }
            }
        }
    }

    // step 5: set symbols of y at dyadic positions except 1 and 2
    for j in 2..t {
        let pos = 1usize << (j as usize);
        if alpha[pos - 1] == 0 {
            out_y[pos] = out_y[pos - 1].wrapping_sub(1);
        } else {
            out_y[pos] = out_y[pos - 1];
        }
    }

    // step 6: set positions 0, 1, 2
    let sum: i64 = out_y[3..nu].iter().map(|&v| v as i64).sum();
    let w = mod_floor(0i64 - sum, 256) as u8;
    let (val_x, val_y, val_z) = if w == 1 {
        (0, 2, 255)
    } else if w == 2 {
        (1, 2, 255)
    } else {
        (0, 1, mod_floor(w as i64 - 1, 256) as u8)
    };

    match (alpha[0], alpha[1]) {
        (0, 0) => {
            out_y[0] = val_z;
            out_y[1] = val_y;
            out_y[2] = val_x;
        }
        (0, 1) => {
            out_y[0] = val_z;
            out_y[1] = val_x;
            out_y[2] = val_y;
        }
        (1, 0) => {
            out_y[0] = val_x;
            out_y[1] = val_z;
            out_y[2] = val_y;
        }
        _ => {
            out_y[0] = val_x;
            out_y[1] = val_y;
            out_y[2] = val_z;
        }
    }
}

fn decode_codeword_q255_into(n: u8, y: &[u8], out_x: &mut [u8]) -> Option<()> {
    let t = ceil_log2(n);
    let k = find_k_q255(n) as usize;
    if out_x.len() < k {
        return None;
    }

    // step 1: decode bytes from systematic positions
    let step_1_num_bytes = core::cmp::max(0i64, n as i64 - 3 * t as i64 + 3) as usize;
    if step_1_num_bytes > 0 {
        let mut xi = 0usize;
        for pos in 0..n {
            if !is_non_systematic(pos, t) {
                out_x[xi] = y[pos as usize];
                xi += 1;
                if xi == step_1_num_bytes {
                    break;
                }
            }
        }
    }

    // step 2: decode bytes from near-dyadic positions
    let mut bytes_done = step_1_num_bytes;
    for j in 3..t {
        let pj = (1u16 << (j as u16)) as u8;

        if pj == n - 1 {
            let v = y[(pj as usize) - 1];
            if v == 0 {
                return None;
            }
            out_x[bytes_done] = v.wrapping_sub(1);
            bytes_done += 1;
            break;
        }

        let r = y[(pj as usize) - 1];
        let l = y[(pj as usize) + 1];
        out_x[bytes_done] = table1_rev(r, l)?;
        bytes_done += 1;
    }

    // step 2b: y[5]
    if y[3] != 255 {
        return None;
    }
    out_x[bytes_done] = y[5]; // table_2_rev is identity

    Some(())
}

// ===== Error correction (binary + q-ary) =====

fn correct_binary_indel_into(n: u8, m: u8, a: u8, y: &[u8], out: &mut [u8]) -> Option<()> {
    let nu = n as usize;
    if out.len() < nu {
        return None;
    }

    let s = compute_syndrome_binary(m, a, y);
    let w: u8 = y.iter().copied().fold(0u8, |acc, v| acc.wrapping_add(v));

    if y.len() == nu - 1 {
        // deletion
        if s == 0 {
            out[..nu - 1].copy_from_slice(y);
            return Some(());
        }

        if s <= w {
            // insert 0; s = number of 1s to the right
            let mut ones_seen: u8 = 0;
            for i in (0..nu - 1).rev() {
                if y[i] == 1 {
                    ones_seen = ones_seen.wrapping_add(1);
                    if ones_seen == s {
                        out[..i].copy_from_slice(&y[..i]);
                        out[i] = 0;
                        out[i + 1..nu].copy_from_slice(&y[i..nu - 1]);
                        return Some(());
                    }
                }
            }
            None
        } else {
            // insert 1; s-w-1 = number of 0s to the left
            let target = s.wrapping_sub(w).wrapping_sub(1);
            if target == 0 {
                out[0] = 1;
                out[1..nu].copy_from_slice(y);
                return Some(());
            }

            let mut zeros_seen: u8 = 0;
            for i in 0..nu - 1 {
                if y[i] == 0 {
                    zeros_seen = zeros_seen.wrapping_add(1);
                    if zeros_seen == target {
                        out[..=i].copy_from_slice(&y[..=i]);
                        out[i + 1] = 1;
                        out[i + 2..nu].copy_from_slice(&y[i + 1..nu - 1]);
                        return Some(());
                    }
                }
            }
            None
        }
    } else {
        // insertion (y.len() == n+1)
        let s_val = mod_floor(m as i64 - n as i64 - 1, m as i64) as u8;
        let m_minus_w = mod_floor(m as i64 - w as i64, m as i64) as u8;

        if s == s_val || s == 0 {
            // last entry inserted
            out.copy_from_slice(&y[..nu]);
            return Some(());
        }
        if s == m_minus_w {
            // first entry inserted
            out.copy_from_slice(&y[1..=nu]);
            return Some(());
        }

        if s > m_minus_w {
            // inserted 0; m-s 1's to the right
            let target = mod_floor(m as i64 - s as i64, m as i64) as u8;
            let mut ones_seen: u8 = 0;
            for i in (2..=nu).rev() {
                if y[i] == 1 {
                    ones_seen = ones_seen.wrapping_add(1);
                    if ones_seen == target {
                        if y[i - 1] == 0 {
                            out[..i - 1].copy_from_slice(&y[..i - 1]);
                            out[i - 1..nu].copy_from_slice(&y[i..i + (nu - (i - 1))]);
                            return Some(());
                        }
                        break;
                    }
                }
            }
            None
        } else {
            // inserted 1; m-w-s 0's to the left
            let target = mod_floor(m as i64 - w as i64 - s as i64, m as i64) as u8;
            let mut zeros_seen: u8 = 0;
            for i in 0..nu {
                if y[i] == 0 {
                    zeros_seen = zeros_seen.wrapping_add(1);
                    if zeros_seen == target {
                        if y[i + 1] == 1 {
                            out[..=i].copy_from_slice(&y[..=i]);
                            out[i + 1..nu].copy_from_slice(&y[i + 2..(i + 2 + (nu - (i + 1)))]);
                            return Some(());
                        }
                        break;
                    }
                }
            }
            None
        }
    }
}

fn correct_q_ary_indel_into(
    n: u8,
    m: u8,
    a: u8,
    b: u8,
    y: &[u8],
    out_y: &mut [u8],      // length n
    alpha_buf: &mut [u8],  // work
    alpha_corr: &mut [u8], // work
) -> Option<()> {
    let nu = n as usize;
    if out_y.len() < nu {
        return None;
    }

    // alpha from received y
    let alen = convert_y_to_alpha(y, alpha_buf);

    // correct alpha using binary VT on length n-1
    let alpha_expected_len = (n - 1) as usize;
    if correct_binary_indel_into(
        n - 1,
        m,
        a,
        &alpha_buf[..alen],
        &mut alpha_corr[..alpha_expected_len],
    )
        .is_none()
    {
        return None;
    }
    if compute_syndrome_binary(m, a, &alpha_corr[..alpha_expected_len]) != 0 {
        return None;
    }

    if alen == nu - 2 {
        // deletion (y.len() == n-1)
        let sum: i64 = y.iter().map(|&v| v as i64).sum();
        let error_symbol = mod_floor(b as i64 - sum, 256) as u8;

        // diff_pos
        let diff_pos = if &alpha_buf[..alen] == &alpha_corr[..alen] {
            nu - 2
        } else {
            let mut dp = 0usize;
            for i in 0..alen {
                if alpha_buf[i] != alpha_corr[i] {
                    dp = i;
                    break;
                }
            }
            dp
        };

        // find del_pos
        for del_pos in (0..=diff_pos + 1).rev() {
            let ok = if del_pos == 0 {
                alpha_corr[0] == if y[0] >= error_symbol { 1 } else { 0 }
            } else if del_pos == nu - 1 {
                alpha_corr[nu - 2] == if error_symbol >= y[nu - 2] { 1 } else { 0 }
            } else {
                (alpha_corr[del_pos - 1] == if error_symbol >= y[del_pos - 1] { 1 } else { 0 })
                    && (alpha_corr[del_pos] == if y[del_pos] >= error_symbol { 1 } else { 0 })
            };

            if ok {
                out_y[..del_pos].copy_from_slice(&y[..del_pos]);
                out_y[del_pos] = error_symbol;
                out_y[del_pos + 1..nu].copy_from_slice(&y[del_pos..nu - 1]);
                if compute_syndrome_q_ary(m, a, b, &out_y[..nu], alpha_buf) == (0, 0) {
                    return Some(());
                } else {
                    return None;
                }
            }
        }
        None
    } else {
        // insertion (y.len() == n+1)
        let sum: i64 = y.iter().map(|&v| v as i64).sum();
        let error_symbol = mod_floor(sum - b as i64, 256) as u8;

        // diff_pos
        let diff_pos = if &alpha_buf[..alpha_expected_len] == &alpha_corr[..alpha_expected_len] {
            nu - 1
        } else {
            let mut dp = 0usize;
            for i in 0..alpha_expected_len {
                if alpha_buf[i] != alpha_corr[i] {
                    dp = i;
                    break;
                }
            }
            dp
        };

        // find ins_pos
        for ins_pos in (0..=diff_pos + 1).rev() {
            let ok = if ins_pos == 0 || ins_pos == nu {
                y[ins_pos] == error_symbol
            } else {
                (y[ins_pos] == error_symbol)
                    && (alpha_corr[ins_pos - 1]
                    == if y[ins_pos + 1] >= y[ins_pos - 1] {
                    1
                } else {
                    0
                })
            };

            if ok {
                out_y[..ins_pos].copy_from_slice(&y[..ins_pos]);
                out_y[ins_pos..nu].copy_from_slice(&y[ins_pos + 1..nu + 1]);
                if compute_syndrome_q_ary(m, a, b, &out_y[..nu], alpha_buf) == (0, 0) {
                    return Some(());
                } else {
                    return None;
                }
            }
        }
        None
    }
}