//! # Unirand Crate
//!
//! This crate implements Marsaglia's Universal Random Number Generator.
//! [More details on the original paper](https://www.sciencedirect.com/science/article/abs/pii/016771529090092L).
//! 
//! ## Overview
//!
//! The RNG uses a sequence of operations to generate uniformly distributed
//! random numbers between 0 and 1. It has been designed with simplicity and 
//! reproducibility in mind.
//!
//! ## Usage Instructions
//!
//! To use this crate, add the following dependency to your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! unirand = "0.2.0"
//! ```
//!
//! Then, you can initialise and use the RNG in your project as follows:
//!
//! ```rust
//!     use unirand::MarsagliaUniRng;
//!
//!     let mut rng = MarsagliaUniRng::new();
//!     rng.rinit(170);
//!     println!("Random number: {}", rng.uni());
//! ```
//!
//! ## Further Information
//!
//! See the documentation for individual functions and methods below for more details.

const LEN_U: usize = 98; // Length of the random values array.

/// A struct representing Marsaglia's Universal Random Number Generator.
pub struct MarsagliaUniRng {
    uni_u: [f32; LEN_U], // Array holding the recent random numbers.
    uni_c: f32,          // Correction to avoid periodicity.
    uni_cd: f32,         // Correction delta value.
    uni_cm: f32,         // Correction modulus.
    uni_ui: usize,       // Current position in the random values array.
    uni_uj: usize,       // Second index used for generating new numbers.
    initialised: bool,   // Set to true once rstart has been called via rinit.
}

impl Default for MarsagliaUniRng {
    fn default() -> Self {
        Self::new()
    }
}

impl MarsagliaUniRng {
    /// Creates a new instance of the RNG.
    ///
    /// # Example
    ///
    /// ```rust
    /// use unirand::MarsagliaUniRng;
    ///
    /// let rng = MarsagliaUniRng::new();
    /// ```
    pub fn new() -> Self {
        Self {
            uni_u: [0.0; LEN_U],
            uni_c: 0.0,
            uni_cd: 0.0,
            uni_cm: 0.0,
            uni_ui: 0,
            uni_uj: 0,
            initialised: false,
        }
    }

    /// Generates a new random float value in [0, 1).
    ///
    /// # Panics
    ///
    /// Panics if called before `rinit`.
    ///
    /// # Example
    ///
    /// ```rust
    /// use unirand::MarsagliaUniRng;
    ///
    /// let mut rng = MarsagliaUniRng::new();
    /// rng.rinit(170);
    /// let number = rng.uni();
    /// println!("Random number: {}", number);
    /// ```
    pub fn uni(&mut self) -> f32 {
        if !self.initialised {
            panic!("uni: called before rinit -- generator not initialised");
        }
        let mut luni = self.uni_u[self.uni_ui] - self.uni_u[self.uni_uj];
        if luni < 0.0 {
            luni += 1.0;
        }
        self.uni_u[self.uni_ui] = luni;
        
        // Adjust indices for the next random number generation.
        // Wraps at 1 -> 97, mirroring the Fortran: I97 = I97 - 1; IF (I97 .EQ. 0) I97 = 97.
        if self.uni_ui == 1 {
            self.uni_ui = 97;
        } else {
            self.uni_ui -= 1;
        }
        if self.uni_uj == 1 {
            self.uni_uj = 97;
        } else {
            self.uni_uj -= 1;
        }
    
        self.uni_c -= self.uni_cd;
        if self.uni_c < 0.0 {
            self.uni_c += self.uni_cm;
        }
    
        luni -= self.uni_c;
        if luni < 0.0 {
            luni += 1.0;
        }
        luni
    }

    /// Initialises the random values array using four seeds.
    /// Called internally by `rinit` after seed validation.
    ///
    /// # Parameters
    ///
    /// - `i`, `j`, `k`, `l`: Pre-validated seed values.
    fn rstart(&mut self, mut i: i32, mut j: i32, mut k: i32, mut l: i32) {
        for ii in 1..=97 {
            let mut s = 0.0;
            let mut t = 0.5;
            for _ in 1..=24 {
                let m = ((i * j % 179) * k) % 179;
                i = j;
                j = k;
                k = m;
                l = (53 * l + 1) % 169;
                if l * m % 64 >= 32 {
                    s += t;
                }
                t *= 0.5;
            }
            self.uni_u[ii] = s;
        }
        // Set fixed correction values.
        self.uni_c = 362436.0 / 16777216.0;
        self.uni_cd = 7654321.0 / 16777216.0;
        self.uni_cm = 16777213.0 / 16777216.0;
        self.uni_ui = 97;
        self.uni_uj = 33;
        self.initialised = true;
    }

    /// Validates and decomposes a single seed into four seeds, then initialises the random values array.
    ///
    /// # Panics
    ///
    /// Panics if `ijkl` is outside the valid range `0..=900_000_000`.
    pub fn rinit(&mut self, ijkl: i32) {
        if !(0..=900_000_000).contains(&ijkl) {
            panic!("rinit: ijkl = {ijkl} -- out of range");
        }

        let ij = ijkl / 30082;
        let kl = ijkl - (30082 * ij);
        let i = ((ij / 177) % 177) + 2;
        let j = (ij % 177) + 2;
        let k = ((kl / 169) % 178) + 1;
        let l = kl % 169;

        // Ranges are guaranteed by the decomposition arithmetic above;
        // verified here in debug builds only.
        debug_assert!((2..=178).contains(&i), "i = {i} out of range");
        debug_assert!((2..=178).contains(&j), "j = {j} out of range");
        debug_assert!((1..=178).contains(&k), "k = {k} out of range");
        debug_assert!((0..=168).contains(&l), "l = {l} out of range");

        self.rstart(i, j, k, l);
    }
}

/// Implements `RngCore` from the `rand_core` crate, enabling use with the broader
/// Rust random number ecosystem (distributions, shuffling, sampling, etc.).
///
/// Note: the underlying generator produces f32 values with 24-bit mantissa precision.
/// `next_u32` therefore has full entropy in the upper 24 bits only; the lower 8 bits
/// are always zero. `next_u64` combines two `next_u32` calls and carries the same
/// limitation. The generator must be initialised with `rinit` before use.
impl rand_core::RngCore for MarsagliaUniRng {
    fn next_u32(&mut self) -> u32 {
        // Scale [0.0, 1.0) to [0, 2^32) and truncate; lower 8 bits are always
        // zero due to f32's 24-bit mantissa.
        (self.uni() * 4_294_967_296.0_f32) as u32
    }

    fn next_u64(&mut self) -> u64 {
        rand_core::impls::next_u64_via_u32(self)
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        rand_core::impls::fill_bytes_via_next(self, dest);
    }
}

/// Produces an infinite sequence of uniform random values in [0, 1).
/// The generator must be initialised with `rinit` before iteration begins;
/// calling `next` on an uninitialised generator will panic.
impl Iterator for MarsagliaUniRng {
    type Item = f32;

    fn next(&mut self) -> Option<f32> {
        Some(self.uni())
    }
}

#[cfg(test)]
mod tests {
    use super::MarsagliaUniRng;

    /// This test checks that a known valid seed produces the expected output.
    #[test]
    fn test_rng_output() {
        let mut rng = MarsagliaUniRng::new();
        rng.rinit(170);
        let random_value = rng.uni();
        let expected = 0.68753344;
        let tolerance = 1e-6;
        assert!(
            (random_value - expected).abs() < tolerance,
            "Expected {}, got {}",
            expected,
            random_value
        );
    }

    /// This test verifies that out-of-range seeds cause expected panics.
    #[test]
    #[should_panic(expected = "rinit: ijkl = -1 -- out of range")]
    fn test_rinit_panics_low_seed() {
        let mut rng = MarsagliaUniRng::new();
        rng.rinit(-1); // Below valid range
    }

    #[test]
    #[should_panic(expected = "rinit: ijkl = 900000001 -- out of range")]
    fn test_rinit_panics_high_seed() {
        let mut rng = MarsagliaUniRng::new();
        rng.rinit(900_000_001); // Above valid range
    }

    /// This is a Statistical Quality Test (SQT) to ensure the RNG produces a uniform distribution.
    #[test]
    fn test_rng_statistics() {
        let mut rng = MarsagliaUniRng::new();
        rng.rinit(170);
        let n = 10_000;
        let sum: f32 = (0..n).map(|_| rng.uni()).sum();
        let mean = sum / n as f32;
        assert!(
            (mean - 0.5).abs() < 0.01,
            "Mean out of expected range: {}",
            mean
        );
    }

    /// This test checks for reproducibility with repeated initialisations using the same seed.
    #[test]
    fn test_rng_reproducibility() {
        let mut rng1 = MarsagliaUniRng::new();
        let mut rng2 = MarsagliaUniRng::new();
        rng1.rinit(42);
        rng2.rinit(42);
        for _ in 0..100 {
            assert!((rng1.uni() - rng2.uni()).abs() < 1e-7);
        }
    }

    /// Verifies that next_u32 produces values consistent with the underlying uni() output.
    #[test]
    fn test_next_u32_matches_uni() {
        use rand_core::RngCore;
        let mut rng_u32 = MarsagliaUniRng::new();
        let mut rng_uni = MarsagliaUniRng::new();
        rng_u32.rinit(42);
        rng_uni.rinit(42);
        for _ in 0..100 {
            let u32_val = rng_u32.next_u32();
            let uni_val = (rng_uni.uni() * 4_294_967_296.0_f32) as u32;
            assert_eq!(u32_val, uni_val);
        }
    }

    /// Verifies that next_u64 produces a valid u64 from two successive next_u32 calls.
    #[test]
    fn test_next_u64_consistent() {
        use rand_core::RngCore;
        let mut rng1 = MarsagliaUniRng::new();
        let mut rng2 = MarsagliaUniRng::new();
        rng1.rinit(42);
        rng2.rinit(42);
        for _ in 0..50 {
            let u64_val = rng1.next_u64();
            // next_u64_via_u32 places the first call in the low bits, second in high bits.
            let lo = rng2.next_u32() as u64;
            let hi = rng2.next_u32() as u64;
            assert_eq!(u64_val, (hi << 32) | lo);
        }
    }

    /// Verifies that fill_bytes fills a buffer without panicking and produces
    /// consistent output across identical seeds.
    #[test]
    fn test_fill_bytes_reproducible() {
        use rand_core::RngCore;
        let mut rng1 = MarsagliaUniRng::new();
        let mut rng2 = MarsagliaUniRng::new();
        rng1.rinit(42);
        rng2.rinit(42);
        let mut buf1 = [0u8; 32];
        let mut buf2 = [0u8; 32];
        rng1.fill_bytes(&mut buf1);
        rng2.fill_bytes(&mut buf2);
        assert_eq!(buf1, buf2);
    }

    /// Verifies that the iterator produces the same sequence as direct uni() calls.
    #[test]
    fn test_iterator_matches_uni() {
        let mut rng_direct = MarsagliaUniRng::new();
        let mut rng_iter = MarsagliaUniRng::new();
        rng_direct.rinit(42);
        rng_iter.rinit(42);
        for (iter_val, direct_val) in rng_iter.take(100).zip((0..100).map(|_| rng_direct.uni())) {
            assert!((iter_val - direct_val).abs() < 1e-7);
        }
    }

    /// Verifies that values produced via the iterator lie within [0.0, 1.0).
    #[test]
    fn test_iterator_output_range() {
        let mut rng = MarsagliaUniRng::new();
        rng.rinit(170);
        for (i, v) in rng.take(1_000).enumerate() {
            assert!(v >= 0.0 && v < 1.0, "Iterator value out of range at step {}: {}", i, v);
        }
    }

    /// Verifies that iterating an uninitialised generator panics with a clear message.
    #[test]
    #[should_panic(expected = "uni: called before rinit -- generator not initialised")]
    fn test_iterator_before_rinit() {
        let mut rng = MarsagliaUniRng::new();
        rng.next();
    }

    /// Verifies that calling uni() before rinit panics with a clear message.
    #[test]
    #[should_panic(expected = "uni: called before rinit -- generator not initialised")]
    fn test_uni_before_rinit() {
        let mut rng = MarsagliaUniRng::new();
        rng.uni();
    }

    /// Verifies that all generated values lie within [0.0, 1.0).
    #[test]
    fn test_rng_output_range() {
        let mut rng = MarsagliaUniRng::new();
        rng.rinit(170);
        for _ in 0..10_000 {
            let v = rng.uni();
            assert!(v >= 0.0 && v < 1.0, "Value out of range: {}", v);
        }
    }

    /// Regression test for the index wrap-around bug: generates enough values to exercise
    /// both the uni_ui cycle (wraps every 97 calls) and the uni_uj first wrap (after 33 calls),
    /// confirming all outputs remain in range across wrap boundaries.
    #[test]
    fn test_rng_wrap_around() {
        let mut rng = MarsagliaUniRng::new();
        rng.rinit(170);
        // 200 values covers two full uni_ui cycles and multiple uni_uj wraps.
        for i in 0..200 {
            let v = rng.uni();
            assert!(v >= 0.0 && v < 1.0, "Value out of range at step {}: {}", i, v);
        }
    }

    /// Verifies that the minimum valid seed (0) and maximum valid seed (900_000_000)
    /// are accepted and produce in-range output.
    #[test]
    fn test_rng_boundary_seeds() {
        for &seed in &[0, 900_000_000] {
            let mut rng = MarsagliaUniRng::new();
            rng.rinit(seed);
            let v = rng.uni();
            assert!(v >= 0.0 && v < 1.0, "Seed {} produced out-of-range value: {}", seed, v);
        }
    }

    /// Verifies that re-initialising with the same seed resets state and yields an
    /// identical sequence to a freshly constructed instance.
    #[test]
    fn test_rng_reinitialisation() {
        let mut rng = MarsagliaUniRng::new();
        rng.rinit(42);
        let first_run: Vec<f32> = (0..50).map(|_| rng.uni()).collect();

        // Re-seed the same instance and confirm the sequence is identical.
        rng.rinit(42);
        let second_run: Vec<f32> = (0..50).map(|_| rng.uni()).collect();

        for (i, (a, b)) in first_run.iter().zip(second_run.iter()).enumerate() {
            assert!((a - b).abs() < 1e-7, "Mismatch at step {}: {} vs {}", i, a, b);
        }
    }

    /// Verifies that distinct seeds produce distinct sequences.
    #[test]
    fn test_rng_distinct_seeds_diverge() {
        let mut rng1 = MarsagliaUniRng::new();
        let mut rng2 = MarsagliaUniRng::new();
        rng1.rinit(1);
        rng2.rinit(2);
        let any_differ = (0..50).any(|_| (rng1.uni() - rng2.uni()).abs() > 1e-7);
        assert!(any_differ, "Different seeds produced identical sequences");
    }

    /// Checks that the variance of the output is close to 1/12, as expected for a
    /// uniform distribution on [0, 1).
    #[test]
    fn test_rng_variance() {
        let mut rng = MarsagliaUniRng::new();
        rng.rinit(170);
        let n = 10_000;
        let values: Vec<f32> = (0..n).map(|_| rng.uni()).collect();
        let mean = values.iter().sum::<f32>() / n as f32;
        let variance = values.iter().map(|x| (x - mean).powi(2)).sum::<f32>() / n as f32;
        let expected_variance: f32 = 1.0 / 12.0;
        assert!(
            (variance - expected_variance).abs() < 0.005,
            "Variance out of expected range: {}",
            variance
        );
    }
}
