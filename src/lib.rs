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
//! unirand = "0.1.2"
//! ```
//!
//! Then, you can initialise and use the RNG in your project as follows:
//!
//! ```rust
//!     use unirand::MarsagliaUniRng;

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
        }
    }

    /// Generates a new random float value between 0 and 1.
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
        let mut luni = self.uni_u[self.uni_ui] - self.uni_u[self.uni_uj];
        if luni < 0.0 {
            luni += 1.0;
        }
        self.uni_u[self.uni_ui] = luni;
        
        // Adjust indices for the next random number generation.
        if self.uni_ui == 0 {
            self.uni_ui = 97;
        } else {
            self.uni_ui -= 1;
        }
        if self.uni_uj == 0 {
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
    ///
    /// # Parameters
    ///
    /// - `i`, `j`, `k`, `l`: The seed values used for initialisation.
    pub fn rstart(&mut self, mut i: i32, mut j: i32, mut k: i32, mut l: i32) {
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
    }

    /// Validates and decomposes a single seed into four seeds, then initialises the random values array.
    ///
    /// # Panics
    ///
    /// Panics if the seed (`ijkl`) is out of range or if the generated seeds are invalid.
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

        if !(1..=178).contains(&i) {
            panic!("rinit: i = {i} -- out of range");
        }
        if !(2..=178).contains(&j) {
            panic!("rinit: j = {j} -- out of range");
        }
        if !(1..=178).contains(&k) {
            panic!("rinit: k = {k} -- out of range");
        }
        if !(0..=168).contains(&l) {
            panic!("rinit: l = {l} -- out of range");
        }
        if i == 1 && j == 1 && k == 1 {
            panic!("rinit: 1 1 1 not allowed for 1st 3 seeds");
        }

        self.rstart(i, j, k, l);
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
}
