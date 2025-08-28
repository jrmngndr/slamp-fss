pub mod config;
pub mod error;
pub mod field_arithmetic;
pub mod linalg_util;
pub mod prg;
pub mod util;
pub mod key_gen;
pub mod eval;

pub use config::*;
pub use error::*;
pub use field_arithmetic::*;
pub use key_gen::*;
pub use eval::*;

pub mod slamp_fss {
    use super::*;

    /// Generate SLAMP-FSS keys for given interpolation points
    ///
    /// # Arguments
    /// * `points` - Vector of (input, output) pairs where inputs are in 64-bit domain
    ///
    /// # Returns
    /// * `Ok((sk1, sk2))` - Two secret key shares on success
    /// * `Err(SlampFssError)` - Error if key generation fails
    pub fn generate_keys(points: &[(u128, u128)]) -> Result<(SecretKeyG, SecretKeyG)> {
        gen_key(points)
    }

    /// Evaluate the SLAMP-FSS function using two key shares
    ///
    /// # Arguments
    /// * `x` - Input value (should be in 64-bit domain, positioned in MSB)
    /// * `sk1` - First secret key share
    /// * `sk2` - Second secret key share
    ///
    /// # Returns
    /// * Function output value
    pub fn evaluate(x: u128, sk1: &SecretKeyG, sk2: &SecretKeyG) -> u128 {
        eval_two_party(x, sk1, sk2)
    }

    /// Evaluate using a single key share
    ///
    /// # Arguments
    /// * `x` - Input value
    /// * `sk` - Secret key share
    ///
    /// # Returns
    /// * Function output value
    pub fn evaluate_single(x: u128, sk: &SecretKeyG) -> u128 {
        eval_function(x, sk)
    }
}

 