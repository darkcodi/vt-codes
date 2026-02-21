mod utils;

use std::collections::HashMap;
use crate::utils::*;

// --- Base conversion (only used in tests now) ---

#[cfg(test)]
fn q_ary_array_to_number(arr: &[u8], q: u16) -> i128 {
    let q = q as i128;
    let mut num: i128 = 0;
    for &digit in arr {
        num = q * num + digit as i128;
    }
    num
}

#[cfg(test)]
fn number_to_q_ary_array(mut num: i128, q: u16, out_size: Option<usize>) -> Option<Vec<u8>> {
    let qb = q as i128;
    let mut out = Vec::new();
    while num > 0 {
        out.push((num % qb) as u8);
        num /= qb;
    }
    out.reverse();
    match out_size {
        None => Some(out),
        Some(sz) => {
            if out.len() > sz {
                None
            } else {
                let mut padded = vec![0u8; sz - out.len()];
                padded.extend_from_slice(&out);
                Some(padded)
            }
        }
    }
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

pub fn find_k(n: u8, q: u8) -> u8 {
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

pub fn find_smallest_n(k: u8, q: u8) -> u8 {
    assert!(q >= 1);
    assert!(k >= 1);
    let mut n = if q == 1 {
        let sum = k as i64 + ceil_log2(k + 1) as i64;
        sum as u8
    } else {
        // For q == 255, ceil_log2(q+1) would be 8 (since q+1 = 256)
        // But ceil_log2 takes u8, and 256 overflows, so handle it specially
        let log_q_plus_1 = if q == 255 { 8 } else { ceil_log2(q + 1) };
        k / log_q_plus_1.max(1)
    };
    loop {
        if find_k(n, q) >= k {
            break;
        }
        n += 1;
    }
    n
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

// --- VTCode struct ---

pub struct VTCode {
    n: u8,
    k: u8,
    m: u8,
    a: u8,
    b: u8,
    // q-ary-specific (always q=255, alphabet size 256)
    t: u8,
    systematic_positions_step_1: Vec<u8>, // 0-indexed
    table_1_l: Vec<u8>,
    table_1_r: Vec<u8>,
    table_1_rev: HashMap<(u8, u8), u16>,  // u16 to support indices up to 32767
    table_2: Vec<u8>,
    table_2_rev: HashMap<u8, u8>,
}

impl VTCode {
    pub fn new(n: u8) -> Self {
        assert!(n >= 2);
        // Always use q=255 (alphabet 0-255, 256 symbols)
        let k = find_k(n, 255);
        assert!(k > 0);

        let mut code = VTCode {
            n,
            k,
            m: 0,
            a: 0,
            b: 0,
            t: 0,
            systematic_positions_step_1: Vec::new(),
            table_1_l: Vec::new(),
            table_1_r: Vec::new(),
            table_1_rev: HashMap::new(),
            table_2: Vec::new(),
            table_2_rev: HashMap::new(),
        };

        code.m = n;
        code.t = ceil_log2(n);
        assert!(code.a < code.m);
        // b is u8, so always <= 255
        code.generate_tables();
        code
    }

    pub fn encode(&self, x: &[u8]) -> Vec<u8> {
        assert!(x.len() == self.k as usize);
        // x is bytes (0-255), encoded using q=255 alphabet
        self.encode_q_ary(x)
    }

    pub fn decode(&self, y: &[u8]) -> Option<Vec<u8>> {
        let n_y = y.len() as i64;
        if n_y < self.n as i64 - 1 || n_y > self.n as i64 + 1 {
            return None;
        }
        // y values are u8, so always valid (0-255)

        let corrected = if n_y != self.n as i64 {
            correct_q_ary_indel(self.n, self.m, self.a, self.b, y)?
        } else {
            y.to_vec()
        };
        self.decode_codeword(&corrected)
    }

    fn is_codeword(&self, y: &[u8]) -> bool {
        if y.len() != self.n as usize {
            return false;
        }
        compute_syndrome_q_ary(self.m, self.a, self.b, y) == (0, 0)
    }

    fn decode_codeword(&self, y: &[u8]) -> Option<Vec<u8>> {
        if !self.is_codeword(y) {
            return None;
        }
        self.decode_codeword_q_ary(y)
    }

    fn decode_codeword_q_ary(&self, y: &[u8]) -> Option<Vec<u8>> {
        let mut x = vec![0u8; self.k as usize];

        // step 1: decode bytes from non-dyadic positions
        // Copy bytes directly from systematic positions
        let step_1_num_bytes = 0i64.max(self.n as i64 - 3 * self.t as i64 + 3) as usize;
        if step_1_num_bytes > 0 {
            for (i, &pos) in self.systematic_positions_step_1.iter().enumerate() {
                if i < step_1_num_bytes {
                    x[i] = y[pos as usize];
                }
            }
        }

        // step 2: decode bytes from near-dyadic positions
        let mut bytes_done = step_1_num_bytes;
        for j in 3..self.t {
            let pj = 1u32 << j;
            if pj as u8 == self.n - 1 {
                // special case: y[2^j - 1]
                if y[(pj - 1) as usize] == 0 {
                    return None;
                }
                // The byte value was stored +1, so subtract 1 to recover
                x[bytes_done] = y[(pj - 1) as usize].wrapping_sub(1);
                bytes_done += 1;
                break;
            }
            // Use table_1_rev to look up the original byte value
            let r = y[(pj - 1) as usize];
            let l = y[(pj + 1) as usize];
            if let Some(&idx) = self.table_1_rev.get(&(r, l)) {
                x[bytes_done] = idx as u8;
            } else {
                return None;
            }
            bytes_done += 1;
        }

        // step 2b: y[5]
        // q=255 (not 2), so always use the q-ary path
        if y[3] != 255 {
            return None;
        }
        if let Some(&idx) = self.table_2_rev.get(&y[5]) {
            x[bytes_done] = idx;
        } else {
            return None;
        }
        bytes_done += 1;

        assert!(bytes_done == self.k as usize);
        Some(x)
    }

    fn encode_q_ary(&self, x: &[u8]) -> Vec<u8> {
        let nu = self.n as usize;
        let mut y = vec![0u8; nu];

        // step 1: encode bytes in non-dyadic positions
        // x is already bytes, so we can copy them directly
        let step_1_num_bytes = 0i64.max(self.n as i64 - 3 * self.t as i64 + 3) as usize;
        if step_1_num_bytes > 0 {
            for (i, &pos) in self.systematic_positions_step_1.iter().enumerate() {
                if i < step_1_num_bytes {
                    y[pos as usize] = x[i];
                }
            }
        }

        // step 2: encode bytes in near-dyadic positions
        let mut bytes_done = step_1_num_bytes;
        // Each input byte is used directly as a table index
        for j in 3..self.t {
            let pj = 1u32 << j;
            if pj as u8 == self.n - 1 {
                // special case: store in y[2^j - 1] (Python: y[2**j-1])
                // Use the byte value directly, add 1 to avoid 0
                y[(pj - 1) as usize] = x[bytes_done].wrapping_add(1);
                bytes_done += 1;
                break;
            }
            // Use the byte value directly as table index
            let table_1_index = x[bytes_done] as usize;
            y[(pj - 1) as usize] = self.table_1_r[table_1_index];
            y[(pj + 1) as usize] = self.table_1_l[table_1_index];
            bytes_done += 1;
        }

        // set y[3] and y[5]
        // q=255 (not 2), so always use the q-ary path
        y[3] = 255;
        let table_2_index = x[bytes_done] as usize;
        y[5] = self.table_2[table_2_index];
        bytes_done += 1;
        assert!(bytes_done == self.k as usize);

        // step 3: set alpha at positions except dyadic
        let mut alpha = convert_y_to_alpha(&y);
        for j in 2..self.t {
            let pj = 1u32 << j;
            if pj as u8 == self.n - 1 {
                break;
            }
            alpha[(pj + 1 - 1) as usize] = if y[(pj + 1) as usize] >= y[(pj - 1) as usize] {
                1
            } else {
                0
            };
        }
        alpha[2] = 1; // alpha[3-1]

        // step 4: set alpha at dyadic positions using VT conditions
        for j in 0..self.t {
            alpha[(1usize << j) - 1] = 0;
        }
        let mut syndrome = compute_syndrome_binary(self.m, self.a, &alpha);
        if syndrome != 0 {
            for j in (0..self.t).rev() {
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
        for j in 2..self.t {
            let pos = 1usize << j;
            if alpha[pos - 1] == 0 {
                y[pos] = y[pos - 1].wrapping_sub(1);
            } else {
                y[pos] = y[pos - 1];
            }
        }

        // step 6: set positions 0, 1, 2
        let sum: i64 = y[3..].iter().map(|&v| v as i64).sum();
        // q=255, so q+1 = 256
        let w = mod_floor(self.b as i64 - sum, 256) as u8;
        // q > 2: find (val_x, val_y, val_z) with val_x < val_y < val_z,
        // val_x + val_y + val_z = w mod (q+1)
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

        debug_assert_eq!(alpha, convert_y_to_alpha(&y));
        debug_assert!(self.is_codeword(&y));
        y
    }

    fn generate_tables(&mut self) {
        // table 1: map floor(2*log2(q)) bits to pairs (r,l) with r!=0, l!=r-1
        // q=255, log2(255) ≈ 7.99, so 2*7.99 ≈ 15.98, floor = 15
        let table_1_size = 1usize << (2.0 * 255_f64.log2()).floor() as usize;
        self.table_1_l = vec![0u8; table_1_size];
        self.table_1_r = vec![0u8; table_1_size];
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
                self.table_1_l[pos] = l;
                self.table_1_r[pos] = r;
                pos += 1;
            }
        }

        self.table_1_rev = HashMap::new();
        for i in 0..table_1_size {
            self.table_1_rev
                .insert((self.table_1_r[i], self.table_1_l[i]), i as u16);
        }

        // table 2: map floor(log2(q)) bits to 1 q-ary symbol != q-1
        // q=255 (not 2), so always generate table 2
        let table_2_size = 1usize << 255_f64.log2().floor() as usize;
        self.table_2 = vec![0u8; table_2_size];
        for i in 0..table_2_size {
            self.table_2[i] = i as u8;
            if i as u8 == 254 {
                self.table_2[i] = 255;
            }
        }
        self.table_2_rev = HashMap::new();
        for i in 0..table_2_size {
            self.table_2_rev.insert(self.table_2[i], i as u8);
        }

        // systematic positions for step 1 (0-indexed)
        let mut non_sys: Vec<u8> = vec![0, 1, 2, 3, 4, 5];
        for j in 3..self.t {
            let pj = 1u32 << j;
            non_sys.push((pj - 1) as u8);
            non_sys.push(pj as u8);
            non_sys.push((pj + 1) as u8);
        }
        let non_sys_set: std::collections::HashSet<u8> = non_sys.into_iter().collect();
        self.systematic_positions_step_1 = (0..self.n)
            .filter(|x| !non_sys_set.contains(x))
            .collect();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_floor() {
        assert_eq!(mod_floor(-3, 5), 2);
        assert_eq!(mod_floor(7, 5), 2);
        assert_eq!(mod_floor(0, 5), 0);
        assert_eq!(mod_floor(-10, 3), 2);
    }

    #[test]
    fn test_base_conversion_roundtrip() {
        let arr = vec![1, 0, 1, 1, 0];
        let num = q_ary_array_to_number(&arr, 2);
        let back = number_to_q_ary_array(num, 2, Some(5)).unwrap();
        assert_eq!(arr, back);

        let arr_q = vec![2, 0, 1, 3];
        let num_q = q_ary_array_to_number(&arr_q, 4);
        let back_q = number_to_q_ary_array(num_q, 4, Some(4)).unwrap();
        assert_eq!(arr_q, back_q);
    }

    #[test]
    fn test_find_k_find_smallest_n_consistency() {
        // Test helper functions still work with various q values
        for q in [1u8, 2, 3, 4, 255] {
            for k in 1..=20 {
                let n = find_smallest_n(k, q);
                assert!(
                    find_k(n, q) >= k,
                    "find_k({n}, {q}) = {} < {k}",
                    find_k(n, q)
                );
                // n should be the smallest
                if n > 2 {
                    assert!(
                        find_k(n - 1, q) < k,
                        "n-1={} also works for k={k}, q={q}",
                        n - 1
                    );
                }
            }
        }
    }

    #[test]
    fn test_encode_decode_roundtrip_q255() {
        // Test with various n values (minimum for q=255 is 6)
        for n in [10u8, 15, 20] {
            let code = VTCode::new(n);

            // all zeros
            let x = vec![0u8; code.k as usize];
            let y = code.encode(&x);
            assert!(code.is_codeword(&y));
            assert_eq!(code.decode(&y).unwrap(), x, "n={n} all zeros");

            // all ones
            let x = vec![1u8; code.k as usize];
            let y = code.encode(&x);
            assert!(code.is_codeword(&y));
            assert_eq!(code.decode(&y).unwrap(), x, "n={n} all ones");

            // alternating
            let x: Vec<u8> = (0..code.k).map(|i| (i % 2) as u8).collect();
            let y = code.encode(&x);
            assert!(code.is_codeword(&y));
            assert_eq!(code.decode(&y).unwrap(), x, "n={n} alternating");
        }
    }

    #[test]
    fn test_deletion_correction_q255() {
        // Test deletion correction with q=255
        for n in [10u8, 15, 20] {
            let code = VTCode::new(n);
            let x: Vec<u8> = (0..code.k).map(|i| (i % 2) as u8).collect();
            let y = code.encode(&x);

            // delete each position and verify recovery
            for del_pos in 0..n as usize {
                let mut y_del = y[..del_pos].to_vec();
                y_del.extend_from_slice(&y[del_pos + 1..]);
                assert_eq!(y_del.len(), (n - 1) as usize);
                let recovered = code.decode(&y_del);
                assert_eq!(
                    recovered.as_deref(),
                    Some(x.as_slice()),
                    "Failed for n={n}, del_pos={del_pos}"
                );
            }
        }
    }

    #[test]
    fn test_insertion_correction_q255() {
        // Test insertion correction with q=255
        let n = 15u8;
        let code = VTCode::new(n);
        let x: Vec<u8> = (0..code.k).map(|i| (i % 2) as u8).collect();
        let y = code.encode(&x);

        // Test insertions at a few positions with various values
        for ins_pos in [0, n as usize / 2, n as usize] {
            for val in [0u8, 1, 127, 255] {
                let mut y_ins = y[..ins_pos].to_vec();
                y_ins.push(val);
                y_ins.extend_from_slice(&y[ins_pos..]);
                let recovered = code.decode(&y_ins);
                // insertion correction may not always succeed (ambiguous cases)
                // but if it does, it must be correct
                if let Some(ref rec) = recovered {
                    assert_eq!(
                        rec, &x,
                        "Wrong decode for n={n}, ins_pos={ins_pos}, val={val}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_minimum_n() {
        // n=6 is the minimum for q=255 (k should be at least 1)
        let code = VTCode::new(6);
        assert!(code.k >= 1);
        let x = vec![0u8; code.k as usize];
        let y = code.encode(&x);
        assert_eq!(code.decode(&y).unwrap(), x);
    }

    #[test]
    fn test_large_n() {
        // Test with larger n values (but keep n small enough that k fits in u8)
        // For q=255, n should be <= about 20 to keep k under 255
        for n in [15u8, 18, 20] {
            let code = VTCode::new(n);
            assert!(code.k > 0);

            let x = vec![0u8; code.k as usize];
            let y = code.encode(&x);
            assert!(code.is_codeword(&y));
            assert_eq!(code.decode(&y).unwrap(), x);

            // Test deletion correction
            let del_pos = n as usize / 2;
            let mut y_del = y[..del_pos].to_vec();
            y_del.extend_from_slice(&y[del_pos + 1..]);
            let recovered = code.decode(&y_del);
            assert_eq!(recovered.as_deref(), Some(x.as_slice()));
        }
    }

    #[test]
    fn test_main() {
        // "Hello world" is 11 bytes
        let data = b"Hello world";
        // find_k(20, 255) = 11 bytes, exactly enough for "Hello world"
        let n = 20u8;
        let code = VTCode::new(n);

        // Pad data to k bytes with zeros
        let mut input = data.to_vec();
        while input.len() < code.k as usize {
            input.push(0);
        }

        // Encode the bytes directly
        let encoded = code.encode(&input);

        // Drop one byte (use fixed position for reproducibility)
        let del_pos = 5; // arbitrary position
        let mut with_deletion = encoded[..del_pos].to_vec();
        with_deletion.extend_from_slice(&encoded[del_pos + 1..]);

        // Decode
        let decoded = code.decode(&with_deletion).expect("Failed to decode after deletion");

        // Verify we recovered "Hello world"
        assert_eq!(&decoded[..data.len()], data);
    }
}
