/// Reusable scratch buffers.
///
/// RAM cost: 3 * 256 = 768 bytes.
#[derive(Clone)]
pub(crate) struct Scratch {
    pub(crate) rx: [u8; 256],    // immutable source copy of received bytes
    pub(crate) cw: [u8; 256],    // corrected / codeword
    pub(crate) tmp: [u8; 256],   // re-encode / alpha_corr work
    pub(crate) alpha: [u8; 256], // alpha work
}

impl Scratch {
    pub(crate) const fn new() -> Self {
        Self {
            rx: [0; 256],
            cw: [0; 256],
            tmp: [0; 256],
            alpha: [0; 256],
        }
    }
}