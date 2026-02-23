# 🔧 vt-codes

[![crates.io](https://img.shields.io/crates/v/vt-codes)](https://crates.io/crates/vt-codes)
[![docs.rs](https://img.shields.io/docsrs/vt-codes)](https://docs.rs/vt-codes)
[![MIT License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

> Varshamov-Tenengolts (VT) codes for single insertion/deletion correction. `no_std`, zero allocations.

**VT codes** are forward error correction codes that can correct a single insertion or deletion in transmitted data — unlike most error correction codes that only handle bit flips.

## ✨ Why?

- 🔧 **In-place** — encode/decode without extra buffers
- 📦 **`no_std`** — works on bare-metal & embedded
- 🚀 **Zero alloc** — no heap, ever
- 📭 **Zero deps** — no external dependencies
- 🛡️ **Indel protection** — corrects single insertions/deletions

## 🚀 Quick look

```rust
let mut buf = [0u8; 256];
buf[..11].copy_from_slice(b"Hello world");

// Encode
let enc_len = vt_codes::vt_encode_in_place(&mut buf, 11)?;

// Simulate a deletion
buf.copy_within(3.., 2); // Remove byte at position 2

// Decode (corrects the deletion)
let dec_len = vt_codes::vt_decode_in_place(&mut buf, enc_len - 1)?;

assert_eq!(&buf[..dec_len], b"Hello world");
```

## 📦 Install

```toml
[dependencies]
vt-codes = "0.1"
```

## Limitations

- **Max message length:** 240 bytes
- **Stack usage:** ~1 KB for internal scratch buffers (4 × 256 bytes)

## License

**MIT** — do whatever.
