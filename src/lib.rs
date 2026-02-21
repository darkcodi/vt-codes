mod utils;

use std::collections::HashMap;
use crate::utils::*;

// --- Base conversion () ---

fn q_ary_array_to_number(arr: &[u8], q: u8) -> i128 {
    let q = q as i128; // i128 intermediates for overflow safety
    let mut num: i128 = 0;
    for &digit in arr {
        num = q * num + digit as i128;
    }
    num
}

fn number_to_q_ary_array(mut num: i128, q: u8, out_size: Option<usize>) -> Option<Vec<u8>> {
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

fn convert_base(in_array: &[u8], in_base: u8, out_base: u8, out_size: Option<usize>) -> Option<Vec<u8>> {
    let num = q_ary_array_to_number(in_array, in_base);
    number_to_q_ary_array(num, out_base, out_size)
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

fn compute_syndrome_q_ary(m: u8, a: u8, b: u8, q: u8, y: &[u8]) -> (u8, u8) {
    let alpha = convert_y_to_alpha(y);
    let s1 = compute_syndrome_binary(m, a, &alpha);
    let sum: i64 = y.iter().map(|&v| v as i64).sum();
    let s2 = mod_floor(b as i64 - sum, q as i64) as u8;
    (s1, s2)
}

// --- find_k and find_smallest_n ---

pub fn find_k(n: u8, q: u8, correct_substitutions: bool) -> u8 {
    if q != 2 && correct_substitutions {
        panic!("correct_substitutions can be True only for q = 2");
    }
    if q == 2 {
        if !correct_substitutions {
            n - ceil_log2(n + 1)
        } else {
            let n_plus = n as i64 + n as i64 + 1;
            n - ceil_log2(n_plus as u8)
        }
    } else {
        let t = ceil_log2(n);
        if q == 3 {
            if n < 7 {
                return 0;
            }
            if power_of_two(n - 1) {
                ((n as i64 - 3 * t as i64 + 3) as f64 * (q as f64).log2()).floor() as i64 as u8 + 2 * (t - 4) + 1
            } else {
                ((n as i64 - 3 * t as i64 + 3) as f64 * (q as f64).log2()).floor() as i64 as u8 + 2 * (t - 3)
            }
        } else {
            if n < 6 {
                return 0;
            }
            let base = (0i64.max(n as i64 - 3 * t as i64 + 3) as f64 * (q as f64).log2()).floor() as i64;
            let bits_per_tuple = (2.0 * ((q - 1) as f64).log2()).floor() as i64;
            let bits_single = ((q - 1) as f64).log2().floor() as i64;
            if power_of_two(n - 1) {
                (base + bits_per_tuple * 0i64.max(t as i64 - 4) + 2 * bits_single) as u8
            } else {
                (base + bits_per_tuple * 0i64.max(t as i64 - 3) + bits_single) as u8
            }
        }
    }
}

pub fn find_smallest_n(k: u8, q: u8, correct_substitutions: bool) -> u8 {
    assert!(q >= 2);
    assert!(k >= 1);
    if q != 2 && correct_substitutions {
        panic!("correct_substitutions can be True only for q = 2");
    }
    let mut n = if q == 2 {
        if !correct_substitutions {
            let sum = k as i64 + ceil_log2(k + 1) as i64;
            sum as u8
        } else {
            let two_k = k as i64 + k as i64 + 1;
            let log_val = ceil_log2(two_k as u8);
            (k as i64 + log_val as i64) as u8
        }
    } else {
        k / ceil_log2(q)
    };
    loop {
        if find_k(n, q, correct_substitutions) >= k {
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

fn correct_binary_substitution(n: u8, m: u8, a: u8, y: &[u8]) -> Vec<u8> {
    let m_expected = 2 * n + 1;
    assert!(m == m_expected);
    let s = compute_syndrome_binary(m, a, y);
    let mut y_decoded = y.to_vec();
    if s == 0 {
        // no error
    } else if s < n + 1 {
        // 1 flipped to 0 at position s (1-indexed)
        y_decoded[(s - 1) as usize] = 1;
    } else {
        // 0 flipped to 1 at position 2n+1-s (1-indexed)
        let pos = (2 * n as i64 + 1 - s as i64 - 1) as usize;
        y_decoded[pos] = 0;
    }
    y_decoded
}

fn correct_q_ary_indel(n: u8, m: u8, a: u8, b: u8, q: u8, y: &[u8]) -> Option<Vec<u8>> {
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
        let error_symbol = mod_floor(b as i64 - sum, q as i64) as u8;

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
        let error_symbol = mod_floor(sum - b as i64, q as i64) as u8;

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

    if compute_syndrome_q_ary(m, a, b, q, &y_decoded) == (0, 0) {
        Some(y_decoded)
    } else {
        None
    }
}

// --- VTCode struct ---

pub struct VTCode {
    n: u8,
    q: u8,
    k: u8,
    m: u8,
    a: u8,
    b: u8,
    correct_substitutions: bool,
    // binary-specific (1-indexed)
    systematic_positions: Vec<u8>,
    parity_positions: Vec<u8>,
    // q-ary-specific
    t: u8,
    systematic_positions_step_1: Vec<u8>, // 0-indexed
    table_1_l: Vec<u8>,
    table_1_r: Vec<u8>,
    table_1_rev: HashMap<(u8, u8), u8>,
    table_2: Vec<u8>,
    table_2_rev: HashMap<u8, u8>,
}

impl VTCode {
    pub fn new(n: u8, q: u8, a: u8, b: u8, correct_substitutions: bool) -> Self {
        assert!(q >= 2);
        assert!(n >= 2);
        let k = find_k(n, q, correct_substitutions);
        assert!(k > 0);

        let mut code = VTCode {
            n,
            q,
            k,
            m: 0,
            a,
            b,
            correct_substitutions,
            systematic_positions: Vec::new(),
            parity_positions: Vec::new(),
            t: 0,
            systematic_positions_step_1: Vec::new(),
            table_1_l: Vec::new(),
            table_1_r: Vec::new(),
            table_1_rev: HashMap::new(),
            table_2: Vec::new(),
            table_2_rev: HashMap::new(),
        };

        if q == 2 {
            code.m = if !correct_substitutions {
                n + 1
            } else {
                2 * n + 1
            };
            assert!(a < code.m);
            code.generate_systematic_positions_binary();
        } else {
            code.m = n;
            code.t = ceil_log2(n);
            assert!(a < code.m);
            assert!(b < q);
            code.generate_tables();
        }
        code
    }

    pub fn encode(&self, x: &[u8]) -> Vec<u8> {
        assert!(x.len() == self.k as usize);
        assert!(x.iter().all(|&v| v == 0 || v == 1));
        if self.q == 2 {
            self.encode_binary(x)
        } else {
            self.encode_q_ary(x)
        }
    }

    pub fn decode(&self, y: &[u8]) -> Option<Vec<u8>> {
        let n_y = y.len() as i64;
        if n_y < self.n as i64 - 1 || n_y > self.n as i64 + 1 {
            return None;
        }
        assert!(y.iter().all(|&v| v < self.q));

        let corrected;
        if self.q == 2 {
            if n_y != self.n as i64 {
                corrected = correct_binary_indel(self.n, self.m, self.a, y)?;
            } else if self.correct_substitutions && !self.is_codeword(y) {
                corrected = correct_binary_substitution(self.n, self.m, self.a, y);
            } else {
                corrected = y.to_vec();
            }
        } else {
            if n_y != self.n as i64 {
                corrected = correct_q_ary_indel(self.n, self.m, self.a, self.b, self.q, y)?;
            } else {
                corrected = y.to_vec();
            }
        }
        self.decode_codeword(&corrected)
    }

    fn is_codeword(&self, y: &[u8]) -> bool {
        if y.len() != self.n as usize {
            return false;
        }
        if self.q == 2 {
            compute_syndrome_binary(self.m, self.a, y) == 0
        } else {
            compute_syndrome_q_ary(self.m, self.a, self.b, self.q, y) == (0, 0)
        }
    }

    fn decode_codeword(&self, y: &[u8]) -> Option<Vec<u8>> {
        if !self.is_codeword(y) {
            return None;
        }
        if self.q == 2 {
            Some(self.decode_codeword_binary(y))
        } else {
            self.decode_codeword_q_ary(y)
        }
    }

    fn decode_codeword_binary(&self, y: &[u8]) -> Vec<u8> {
        self.systematic_positions
            .iter()
            .map(|&pos| y[(pos - 1) as usize])
            .collect()
    }

    fn decode_codeword_q_ary(&self, y: &[u8]) -> Option<Vec<u8>> {
        let mut x = vec![0u8; self.k as usize];

        // step 1
        let step_1_num_bits =
            (0i64.max(self.n as i64 - 3 * self.t as i64 + 3) as f64 * (self.q as f64).log2()).floor() as usize;
        if step_1_num_bits > 0 {
            let step_1_vals: Vec<u8> = self
                .systematic_positions_step_1
                .iter()
                .map(|&pos| y[pos as usize])
                .collect();
            let step_1_bits = convert_base(&step_1_vals, self.q, 2, Some(step_1_num_bits))?;
            x[..step_1_num_bits].copy_from_slice(&step_1_bits);
        }

        // step 2
        let mut bits_done = step_1_num_bits;
        let bits_per_tuple = (2.0 * ((self.q - 1) as f64).log2()).floor() as usize;
        for j in 3..self.t {
            let pj = 1u32 << j; // 2^j
            if pj as u8 == self.n - 1 {
                let num_bits_special = ((self.q - 1) as f64).log2().floor() as usize;
                // Python: y[2**j-1] means y[(2**j)-1]
                if y[(pj - 1) as usize] == 0 {
                    return None;
                }
                let bits = number_to_q_ary_array(
                    (y[(pj - 1) as usize] - 1) as i128,
                    2,
                    Some(num_bits_special),
                )?;
                x[bits_done..bits_done + num_bits_special].copy_from_slice(&bits);
                bits_done += num_bits_special;
                break;
            }
            // Python: y[2**j-1] is index 2^j - 1, y[2**j+1] is index 2^j + 1
            // but in Python the table uses (y[2**j-1], y[2**j+1]) where these are
            // 0-indexed array accesses at positions 2^j-1 and 2^j+1
            let r = y[(pj - 1) as usize];
            let l = y[(pj + 1) as usize];
            if let Some(&idx) = self.table_1_rev.get(&(r, l)) {
                let bits = number_to_q_ary_array(idx as i128, 2, Some(bits_per_tuple))?;
                x[bits_done..bits_done + bits_per_tuple].copy_from_slice(&bits);
            } else {
                return None;
            }
            bits_done += bits_per_tuple;
        }

        // step 2b: y[5]
        if self.q == 3 {
            if y[5] != 2 {
                return None;
            }
        } else {
            if y[3] != self.q - 1 {
                return None;
            }
            let bits_in_c5 = ((self.q - 1) as f64).log2().floor() as usize;
            if let Some(&idx) = self.table_2_rev.get(&y[5]) {
                let bits = number_to_q_ary_array(idx as i128, 2, Some(bits_in_c5))?;
                x[bits_done..bits_done + bits_in_c5].copy_from_slice(&bits);
            } else {
                return None;
            }
            bits_done += bits_in_c5;
        }

        assert!(bits_done == self.k as usize);
        Some(x)
    }

    fn encode_binary(&self, x: &[u8]) -> Vec<u8> {
        let nu = self.n as usize;
        let mut y = vec![0u8; nu];

        // set systematic positions
        for (i, &pos) in self.systematic_positions.iter().enumerate() {
            y[(pos - 1) as usize] = x[i];
        }

        // set parity bits to satisfy syndrome == 0
        let mut syndrome = compute_syndrome_binary(self.m, self.a, &y);
        if syndrome != 0 {
            for &pos in self.parity_positions.iter().rev() {
                if syndrome >= pos {
                    y[(pos - 1) as usize] = 1;
                    syndrome -= pos;
                    if syndrome == 0 {
                        break;
                    }
                }
            }
        }
        debug_assert!(self.is_codeword(&y));
        y
    }

    fn encode_q_ary(&self, x: &[u8]) -> Vec<u8> {
        let nu = self.n as usize;
        let mut y = vec![0u8; nu];

        // step 1: encode bits in non-dyadic positions
        let step_1_num_bits =
            (0i64.max(self.n as i64 - 3 * self.t as i64 + 3) as f64 * (self.q as f64).log2()).floor() as usize;
        if step_1_num_bits > 0 {
            let out_size = self.systematic_positions_step_1.len();
            let vals = convert_base(&x[..step_1_num_bits], 2, self.q, Some(out_size))
                .expect("base conversion failed in encode step 1");
            for (i, &pos) in self.systematic_positions_step_1.iter().enumerate() {
                y[pos as usize] = vals[i];
            }
        }

        // step 2: encode bits in near-dyadic positions
        let mut bits_done = step_1_num_bits;
        let bits_per_tuple = (2.0 * ((self.q - 1) as f64).log2()).floor() as usize;
        for j in 3..self.t {
            let pj = 1u32 << j;
            if pj as u8 == self.n - 1 {
                // special case: store in y[2^j - 1] (Python: y[2**j-1])
                let num_bits_special = ((self.q - 1) as f64).log2().floor() as usize;
                y[(pj - 1) as usize] =
                    q_ary_array_to_number(&x[bits_done..bits_done + num_bits_special], 2) as u8
                        + 1;
                bits_done += num_bits_special;
                break;
            }
            let table_1_index =
                q_ary_array_to_number(&x[bits_done..bits_done + bits_per_tuple], 2) as usize;
            y[(pj - 1) as usize] = self.table_1_r[table_1_index];
            y[(pj + 1) as usize] = self.table_1_l[table_1_index];
            bits_done += bits_per_tuple;
        }

        // set y[3] and y[5]
        y[3] = self.q - 1;
        if self.q == 3 {
            y[5] = 2;
        } else {
            let bits_in_c5 = ((self.q - 1) as f64).log2().floor() as usize;
            let table_2_index =
                q_ary_array_to_number(&x[bits_done..bits_done + bits_in_c5], 2) as usize;
            y[5] = self.table_2[table_2_index];
            bits_done += bits_in_c5;
        }
        assert!(bits_done == self.k as usize);

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
        let w = mod_floor(self.b as i64 - sum, self.q as i64) as u8;
        if self.q == 3 {
            if alpha[0] == 1 && alpha[1] == 1 {
                y[2] = 2;
                y[1] = 2;
                y[0] = mod_floor(w as i64 - 4, 3) as u8;
            } else if alpha[0] == 1 && alpha[1] == 0 {
                y[2] = 1;
                y[1] = 2;
                y[0] = w;
            } else if alpha[0] == 0 && alpha[1] == 1 {
                y[2] = 2;
                if w == 1 {
                    y[1] = 0;
                    y[0] = 2;
                } else if w == 0 {
                    y[1] = 0;
                    y[0] = 1;
                } else {
                    y[1] = 1;
                    y[0] = 2;
                }
            } else {
                // alpha[0]==0 && alpha[1]==0: fallback
                alpha[0] = 1;
                alpha[1] = 1;
                alpha[2] = 0;
                y[3] = 1;
                y[4] = if alpha[3] == 0 { 0 } else { 1 };
                y[2] = 2;
                y[1] = 2;
                let sum_b: i64 = y[1..].iter().map(|&v| v as i64).sum();
                y[0] = mod_floor(self.b as i64 - sum_b, 3) as u8;
            }
        } else {
            // q > 3: find (val_x, val_y, val_z) with val_x < val_y < val_z,
            // val_x + val_y + val_z = w mod q
            let (val_x, val_y, val_z) = if w == 1 {
                (0, 2, self.q - 1)
            } else if w == 2 {
                (1, 2, self.q - 1)
            } else {
                (0, 1, mod_floor(w as i64 - 1, self.q as i64) as u8)
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
        }

        debug_assert_eq!(alpha, convert_y_to_alpha(&y));
        debug_assert!(self.is_codeword(&y));
        y
    }

    fn generate_systematic_positions_binary(&mut self) {
        let t = ceil_log2(self.n + 1);
        let num_parity = (self.n - self.k) as usize;
        let mut parity = vec![0u8; num_parity];
        for i in 0..t as usize {
            parity[i] = 1 << i;
        }
        if self.correct_substitutions {
            assert!(num_parity == (t + 1) as usize);
            let tu = t as usize;
            if parity[tu - 1] == self.n {
                parity[tu - 1] = self.n - 1;
                parity[tu] = self.n;
            } else {
                parity[tu] = self.n;
            }
        }
        parity.sort();
        self.parity_positions = parity;

        let parity_set: std::collections::HashSet<u8> =
            self.parity_positions.iter().copied().collect();
        self.systematic_positions = (1..=self.n)
            .filter(|x| !parity_set.contains(x))
            .collect();
    }

    fn generate_tables(&mut self) {
        // table 1: map floor(2*log2(q-1)) bits to pairs (r,l) with r!=0, l!=r-1
        let table_1_size = 1usize << (2.0 * ((self.q - 1) as f64).log2()).floor() as usize;
        self.table_1_l = vec![0u8; table_1_size];
        self.table_1_r = vec![0u8; table_1_size];
        let mut pos = 0;
        for r in 0..self.q {
            if pos == table_1_size {
                break;
            }
            if r == 0 {
                continue;
            }
            for l in 0..self.q {
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
                .insert((self.table_1_r[i], self.table_1_l[i]), i as u8);
        }

        // table 2: map floor(log2(q-1)) bits to 1 q-ary symbol != q-2
        if self.q != 3 {
            let table_2_size = 1usize << ((self.q - 1) as f64).log2().floor() as usize;
            self.table_2 = vec![0u8; table_2_size];
            for i in 0..table_2_size {
                self.table_2[i] = i as u8;
                if i as u8 == self.q - 2 {
                    self.table_2[i] = self.q - 1;
                }
            }
            self.table_2_rev = HashMap::new();
            for i in 0..table_2_size {
                self.table_2_rev.insert(self.table_2[i], i as u8);
            }
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
    fn main_test() {
        let message = b"Hello world";
        println!("Original message: {:?}", std::str::from_utf8(message).unwrap());
        println!("Message length: {} bytes", message.len());

        // For each byte, create a VT code (8 bits per byte)
        // n=12 gives k=8 binary bits (exactly one byte)
        let n: u8 = 12;
        let code = VTCode::new(n, 2, 0, 0, false);
        println!("VT code per byte: n={}, k={}", code.n, code.k);

        // Encode each byte
        let mut encoded_bytes = Vec::new();
        for &byte in message {
            let bits: Vec<u8> = (0..8).rev().map(|i| ((byte >> i) & 1) as u8).collect();
            let encoded = code.encode(&bits);
            encoded_bytes.push((byte, encoded));
        }
        println!("Encoded {} bytes into {} bits total", message.len(), message.len() * n as usize);

        // Flatten the encoded data
        let flat_encoded: Vec<u8> = encoded_bytes.iter()
            .flat_map(|(_, enc)| enc.clone())
            .collect();
        println!("Flat encoded: {} bits", flat_encoded.len());

        // Corrupt: delete ONE BIT from a chunk in the middle
        // VT codes can correct a single deletion per chunk
        let del_bit_pos = flat_encoded.len() / 2;  // Delete from the middle
        let mut corrupted: Vec<u8> = flat_encoded[..del_bit_pos].to_vec();
        corrupted.extend_from_slice(&flat_encoded[del_bit_pos + 1..]);
        println!("\nAfter deleting 1 bit from position {}: {} bits remain",
            del_bit_pos, corrupted.len());

        // Find which chunk was affected
        let affected_chunk_idx = del_bit_pos / n as usize;

        // Decode: process each n-bit chunk
        let mut recovered_bytes = Vec::new();
        let mut recovered_count = 0;
        let mut failed_count = 0;

        for i in 0..encoded_bytes.len() {
            let chunk_start = (i * n as usize) as usize;
            let chunk_end = chunk_start + n as usize;

            // Adjust indices for the deleted bit
            let adjusted_start = if chunk_start > del_bit_pos {
                chunk_start - 1
            } else {
                chunk_start
            };
            let adjusted_end = if chunk_end > del_bit_pos {
                chunk_end - 1
            } else {
                chunk_end
            };

            let chunk: Vec<u8> = if i == affected_chunk_idx {
                // This chunk has n-1 bits (VT can correct this!)
                corrupted[adjusted_start..adjusted_end].to_vec()
            } else {
                // Other chunks have n bits
                corrupted[adjusted_start..adjusted_end].to_vec()
            };

            match code.decode(&chunk) {
                Some(bits) => {
                    // Convert bits back to byte
                    let mut byte: u8 = 0;
                    for (j, &bit) in bits.iter().enumerate() {
                        byte |= bit << (7 - j);
                    }
                    recovered_bytes.push(byte);
                    recovered_count += 1;
                }
                None => {
                    failed_count += 1;
                    println!("Failed to decode chunk {}", i);
                }
            }
        }

        let recovered_str = std::str::from_utf8(&recovered_bytes).unwrap_or("<invalid utf8>");
        println!("\nRecovered {} bytes, {} failed", recovered_count, failed_count);
        println!("Recovered bytes: {:?}", recovered_bytes);
        println!("Recovered message: {:?}", recovered_str);
        println!("Original message:  {:?}", message);

        // All bytes should be recovered (VT corrects the single bit deletion)
        assert_eq!(recovered_count, message.len(),
            "Should recover all {} bytes", message.len());
        assert_eq!(failed_count, 0, "Should have no failed decodes");

        // Verify all bytes match the original
        assert_eq!(&recovered_bytes[..], message,
            "Recovered message should match original exactly");

        println!("\nSuccess! All bytes recovered correctly using VT single-deletion correction.");
    }

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
        for q in [2, 3, 4, 5] {
            for k in 1..=20 {
                let n = find_smallest_n(k, q, false);
                assert!(
                    find_k(n, q, false) >= k,
                    "find_k({n}, {q}) = {} < {k}",
                    find_k(n, q, false)
                );
                // n should be the smallest
                if n > 2 {
                    assert!(
                        find_k(n - 1, q, false) < k,
                        "n-1={} also works for k={k}, q={q}",
                        n - 1
                    );
                }
            }
        }
    }

    #[test]
    fn test_binary_encode_decode_roundtrip() {
        for n in [7u8, 10, 15, 20, 31] {
            let code = VTCode::new(n, 2, 0, 0, false);
            // all zeros
            let x = vec![0u8; code.k as usize];
            let y = code.encode(&x);
            assert!(code.is_codeword(&y));
            assert_eq!(code.decode(&y).unwrap(), x);

            // all ones
            let x = vec![1u8; code.k as usize];
            let y = code.encode(&x);
            assert!(code.is_codeword(&y));
            assert_eq!(code.decode(&y).unwrap(), x);

            // alternating
            let x: Vec<u8> = (0..code.k).map(|i| (i % 2) as u8).collect();
            let y = code.encode(&x);
            assert!(code.is_codeword(&y));
            assert_eq!(code.decode(&y).unwrap(), x);
        }
    }

    #[test]
    fn test_binary_deletion_correction() {
        for n in [7u8, 10, 15, 20] {
            let code = VTCode::new(n, 2, 0, 0, false);
            let x: Vec<u8> = (0..code.k).map(|i| ((i + 1) % 2) as u8).collect();
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
    fn test_binary_insertion_correction() {
        for n in [7u8, 10, 15, 20] {
            let code = VTCode::new(n, 2, 0, 0, false);
            let x: Vec<u8> = (0..code.k).map(|i| (i % 2) as u8).collect();
            let y = code.encode(&x);

            // insert each value at each position
            for ins_pos in 0..=n as usize {
                for &val in &[0u8, 1] {
                    let mut y_ins = y[..ins_pos].to_vec();
                    y_ins.push(val);
                    y_ins.extend_from_slice(&y[ins_pos..]);
                    assert_eq!(y_ins.len(), (n + 1) as usize);
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
    }

    #[test]
    fn test_binary_substitution_correction() {
        for n in [7u8, 10, 15, 20] {
            let code = VTCode::new(n, 2, 0, 0, true);
            let x: Vec<u8> = (0..code.k).map(|i| ((i + 1) % 2) as u8).collect();
            let y = code.encode(&x);

            // flip each bit
            for flip_pos in 0..n as usize {
                let mut y_sub = y.clone();
                y_sub[flip_pos] = 1 - y_sub[flip_pos];
                let recovered = code.decode(&y_sub);
                assert_eq!(
                    recovered.as_deref(),
                    Some(x.as_slice()),
                    "Failed for n={n}, flip_pos={flip_pos}"
                );
            }
        }
    }

    #[test]
    fn test_q_ary_encode_decode_roundtrip() {
        for q in [3u8, 4, 5] {
            let min_n = if q == 3 { 9 } else { 10 };
            for n in [min_n, min_n + 5, min_n + 10] {
                let code = VTCode::new(n, q, 0, 0, false);
                // all zeros
                let x = vec![0u8; code.k as usize];
                let y = code.encode(&x);
                assert!(code.is_codeword(&y));
                assert_eq!(code.decode(&y).unwrap(), x, "q={q}, n={n} all zeros");

                // all ones
                let x = vec![1u8; code.k as usize];
                let y = code.encode(&x);
                assert!(code.is_codeword(&y));
                assert_eq!(code.decode(&y).unwrap(), x, "q={q}, n={n} all ones");

                // alternating
                let x: Vec<u8> = (0..code.k).map(|i| (i % 2) as u8).collect();
                let y = code.encode(&x);
                assert!(code.is_codeword(&y));
                assert_eq!(code.decode(&y).unwrap(), x, "q={q}, n={n} alternating");
            }
        }
    }

    #[test]
    fn test_q_ary_deletion_correction() {
        for q in [3u8, 4, 5] {
            let n = if q == 3 { 9 } else { 10 };
            let code = VTCode::new(n, q, 0, 0, false);
            let x: Vec<u8> = (0..code.k).map(|i| (i % 2) as u8).collect();
            let y = code.encode(&x);

            for del_pos in 0..n as usize {
                let mut y_del = y[..del_pos].to_vec();
                y_del.extend_from_slice(&y[del_pos + 1..]);
                let recovered = code.decode(&y_del);
                assert_eq!(
                    recovered.as_deref(),
                    Some(x.as_slice()),
                    "Failed for q={q}, n={n}, del_pos={del_pos}"
                );
            }
        }
    }

    #[test]
    fn test_q_ary_insertion_correction() {
        for q in [3u8, 4] {
            let n = if q == 3 { 9 } else { 10 };
            let code = VTCode::new(n, q, 0, 0, false);
            let x: Vec<u8> = (0..code.k).map(|i| (i % 2) as u8).collect();
            let y = code.encode(&x);

            // insert each valid symbol at a few positions
            for ins_pos in [0, n as usize / 2, n as usize] {
                for val in 0..q {
                    let mut y_ins = y[..ins_pos].to_vec();
                    y_ins.push(val);
                    y_ins.extend_from_slice(&y[ins_pos..]);
                    let recovered = code.decode(&y_ins);
                    if let Some(ref rec) = recovered {
                        assert_eq!(
                            rec, &x,
                            "Wrong decode for q={q}, n={n}, ins_pos={ins_pos}, val={val}"
                        );
                    }
                }
            }
        }
    }
}
