/// Python-style non-negative modulo: ((a % m) + m) % m
pub(crate) fn mod_floor(a: i64, m: i64) -> i64 {
    ((a % m) + m) % m
}

pub(crate) fn ceil_log2(x: u8) -> u8 {
    if x <= 1 {
        return 0;
    }
    let lz = (x - 1).leading_zeros() as u8; // 0..=8
    8 - lz
}

pub(crate) fn power_of_two(num: u8) -> bool {
    num > 0 && (num & (num - 1)) == 0
}
