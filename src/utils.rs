/// Python-style non-negative modulo: ((a % m) + m) % m
pub fn mod_floor(a: i64, m: i64) -> i64 {
    ((a % m) + m) % m
}

pub fn ceil_log2(x: i64) -> i64 {
    (x as f64).log2().ceil() as i64
}

pub fn power_of_two(num: i64) -> bool {
    num > 0 && (num & (num - 1)) == 0
}
