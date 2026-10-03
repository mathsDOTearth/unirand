//! # Unirand Crate
//!
//! This crate implements Marsaglia's Universal Random Number Generator, as defined in
//! G. Marsaglia, A. Zaman and W. W. Tsang, "Toward a universal random number generator",
//! *Statistics & Probability Letters* 9(1), 1990, pp. 35-39
//! ([article](https://www.sciencedirect.com/science/article/abs/pii/016771529090092L)).
//!
//! ## Overview
//!
//! The generator combines a lagged Fibonacci generator, `x(n) = x(n-97) - x(n-33) mod 1`,
//! with an arithmetic sequence `c(n) = c(n-1) - 7654321/2^24 mod 16777213/2^24`. Every
//! value it produces is an exact multiple of 2^-24, and the period is approximately 2^144.
//! Because all quantities are exact 24-bit fractions, the state is held as 24-bit integers;
//! the output is bit-for-bit identical to the single-precision reference implementation on
//! any platform.
//!
//! The generator is not suitable for cryptographic or security-sensitive applications.
//!
//! ## Usage
//!
//! ```toml
//! [dependencies]
//! unirand = "0.3.1"
//! ```
//!
//! ```rust
//! use unirand::MarsagliaUniRng;
//!
//! let mut rng = MarsagliaUniRng::new(170);
//! println!("Random number: {}", rng.uni());
//! ```
//!
//! ## Seeding
//!
//! Two seeding schemes are provided:
//!
//! - [`MarsagliaUniRng::new`] and [`MarsagliaUniRng::try_new`] accept a single seed in
//!   `0..=900_000_000`, decomposed into four seeds as in F. James's RANMAR
//!   (*Computer Physics Communications* 60, 1990).
//! - [`MarsagliaUniRng::from_seeds`] accepts the four seeds `i, j, k, l` of the paper
//!   directly.
//!
//! The single seed `1802 * 30082 + 9373` is equivalent to the four seeds `(12, 34, 56, 78)`.
//!
//! ## Floating-point output
//!
//! [`MarsagliaUniRng::uni`] returns each value of the paper as an `f32`, and
//! [`MarsagliaUniRng::uni_f64`] returns the same value converted exactly to `f64`. Every
//! value lies in [0, 1) and is an exact multiple of 2^-24, so an `f64` obtained this way
//! has 24-bit resolution. `0.0` occurs with probability 2^-24 per value and `1.0` never
//! occurs; where an open lower bound is required, as for `ln(u)`, `1.0 - u` lies in (0, 1].
//!
//! ```rust
//! use unirand::MarsagliaUniRng;
//!
//! let mut rng = MarsagliaUniRng::new(170);
//! let single: f32 = rng.uni();
//! let double: f64 = rng.uni_f64();
//! let many: Vec<f64> = rng.by_ref().take(1_000).map(f64::from).collect();
//! assert!((0.0..1.0).contains(&single) && (0.0..1.0).contains(&double));
//! assert!(many.iter().all(|u| (0.0..1.0).contains(u)));
//! ```
//!
//! The floating-point methods of `rand`, such as `random::<f32>()` and `random::<f64>()`,
//! draw on the `rand_core` interface, which consumes two or three values of the generator
//! per result. Their output is uniform but is not the sequence defined in the paper; use
//! `uni` or `uni_f64` where that sequence is required.

#![no_std]

use core::fmt;
use core::iter::FusedIterator;

/// Length of the lagged Fibonacci state, the larger lag of the generator.
const LEN_U: usize = 97;
/// Initial index of the larger lag (97), expressed 0-based.
const INITIAL_I97: usize = 96;
/// Initial index of the smaller lag (33), expressed 0-based.
const INITIAL_J97: usize = 32;
/// Mask reducing a value modulo 2^24.
const MASK_24: u32 = (1 << 24) - 1;
/// Scale factor converting a 24-bit integer to a fraction in [0, 1).
const SCALE_24: f32 = 1.0 / 16_777_216.0;
/// Initial value of the arithmetic sequence, 362436 / 2^24.
const C_INITIAL: u32 = 362_436;
/// Decrement of the arithmetic sequence, 7654321 / 2^24.
const C_DELTA: u32 = 7_654_321;
/// Modulus of the arithmetic sequence, 16777213 / 2^24.
const C_MODULUS: u32 = 16_777_213;

/// Error returned when a seed lies outside its valid range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeedError {
    /// The single seed exceeds [`MarsagliaUniRng::MAX_SEED`].
    SeedOutOfRange(u32),
    /// The four seeds violate the constraints of the paper: `i`, `j`, `k` must lie in
    /// `1..=178` and not all equal 1, and `l` must lie in `0..=168`.
    InvalidSeeds { i: u32, j: u32, k: u32, l: u32 },
}

impl fmt::Display for SeedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            SeedError::SeedOutOfRange(seed) => write!(
                f,
                "seed {seed} out of range 0..={}",
                MarsagliaUniRng::MAX_SEED
            ),
            SeedError::InvalidSeeds { i, j, k, l } => write!(
                f,
                "seeds (i, j, k, l) = ({i}, {j}, {k}, {l}) invalid: i, j, k must lie in \
                 1..=178 and not all equal 1, l must lie in 0..=168"
            ),
        }
    }
}

impl core::error::Error for SeedError {}

/// Decomposes a single seed into the four seeds of the paper, following James's RANMAR:
/// `ij = seed / 30082` determines `i` and `j`, and `kl = seed mod 30082` determines `k`
/// and `l`.
fn decompose_seed(seed: u32) -> (u32, u32, u32, u32) {
    let ij = seed / 30082;
    let kl = seed % 30082;
    let i = (ij / 177) % 177 + 2;
    let j = ij % 177 + 2;
    let k = (kl / 169) % 178 + 1;
    let l = kl % 169;
    (i, j, k, l)
}

/// Returns whether four seeds satisfy the constraints of the paper: `i`, `j`, `k` in
/// `1..=178` and not all equal 1, and `l` in `0..=168`.
fn seeds_are_valid(i: u32, j: u32, k: u32, l: u32) -> bool {
    [i, j, k].iter().all(|seed| (1..=178).contains(seed))
        && !(i == 1 && j == 1 && k == 1)
        && l <= 168
}

/// Marsaglia's Universal Random Number Generator.
///
/// The state is held as 24-bit integers `x` representing the fractions `x / 2^24` of the
/// paper, so that all arithmetic is exact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MarsagliaUniRng {
    u: [u32; LEN_U], // Lagged Fibonacci state, U(1..=97) of the paper, stored 0-based.
    c: u32,          // Current term of the arithmetic sequence.
    i97: usize,      // Index of the larger lag.
    j97: usize,      // Index of the smaller lag.
}

impl MarsagliaUniRng {
    /// Largest valid single seed.
    pub const MAX_SEED: u32 = 900_000_000;

    /// Creates a generator from a single seed in `0..=MAX_SEED`.
    ///
    /// # Panics
    ///
    /// Panics if `seed` exceeds [`MAX_SEED`](Self::MAX_SEED). See
    /// [`try_new`](Self::try_new) for a non-panicking alternative.
    ///
    /// # Example
    ///
    /// ```rust
    /// use unirand::MarsagliaUniRng;
    ///
    /// let mut rng = MarsagliaUniRng::new(170);
    /// let number = rng.uni();
    /// assert!((0.0..1.0).contains(&number));
    /// ```
    #[must_use]
    pub fn new(seed: u32) -> Self {
        match Self::try_new(seed) {
            Ok(rng) => rng,
            Err(error) => panic!("MarsagliaUniRng::new: {error}"),
        }
    }

    /// Creates a generator from a single seed, returning an error if `seed` exceeds
    /// [`MAX_SEED`](Self::MAX_SEED).
    ///
    /// The seed is decomposed into `ij = seed / 30082` and `kl = seed mod 30082`, and
    /// thence into the four seeds of the paper, following James's RANMAR.
    ///
    /// # Errors
    ///
    /// Returns [`SeedError::SeedOutOfRange`] if `seed` exceeds [`MAX_SEED`](Self::MAX_SEED).
    pub fn try_new(seed: u32) -> Result<Self, SeedError> {
        if seed > Self::MAX_SEED {
            return Err(SeedError::SeedOutOfRange(seed));
        }
        let (i, j, k, l) = decompose_seed(seed);
        Self::from_seeds(i, j, k, l)
    }

    /// Creates a generator from the four seeds of the paper.
    ///
    /// `i`, `j` and `k` must lie in `1..=178` and must not all equal 1; `l` must lie in
    /// `0..=168`.
    ///
    /// # Errors
    ///
    /// Returns [`SeedError::InvalidSeeds`] if any of these constraints is violated.
    ///
    /// # Example
    ///
    /// ```rust
    /// use unirand::MarsagliaUniRng;
    ///
    /// let rng = MarsagliaUniRng::from_seeds(12, 34, 56, 78).unwrap();
    /// assert_eq!(rng, MarsagliaUniRng::new(1802 * 30082 + 9373));
    /// ```
    pub fn from_seeds(i: u32, j: u32, k: u32, l: u32) -> Result<Self, SeedError> {
        if !seeds_are_valid(i, j, k, l) {
            return Err(SeedError::InvalidSeeds { i, j, k, l });
        }
        Ok(Self::rstart(i, j, k, l))
    }

    /// Re-seeds the generator in place from a single seed, discarding the current state.
    ///
    /// # Panics
    ///
    /// Panics if `seed` exceeds [`MAX_SEED`](Self::MAX_SEED).
    pub fn rinit(&mut self, seed: u32) {
        *self = Self::new(seed);
    }

    /// Fills the lagged Fibonacci state from pre-validated seeds. Each of the 97 entries is
    /// a 24-bit fraction whose bits, most significant first, are produced by a combination
    /// of a 3-lag multiplicative generator modulo 179 and a linear congruential generator
    /// modulo 169.
    #[must_use]
    fn rstart(mut i: u32, mut j: u32, mut k: u32, mut l: u32) -> Self {
        let mut u = [0; LEN_U];
        for entry in &mut u {
            let mut bits = 0;
            for _ in 0..24 {
                let m = ((i * j % 179) * k) % 179;
                i = j;
                j = k;
                k = m;
                l = (53 * l + 1) % 169;
                bits = (bits << 1) | u32::from(l * m % 64 >= 32);
            }
            *entry = bits;
        }
        Self {
            u,
            c: C_INITIAL,
            i97: INITIAL_I97,
            j97: INITIAL_J97,
        }
    }

    /// Advances the generator and returns the next value as a 24-bit integer in
    /// `0..2^24`, equal to `uni() * 2^24`.
    #[inline]
    pub fn uni_u24(&mut self) -> u32 {
        // x(n) = x(n-97) - x(n-33) mod 1, computed modulo 2^24.
        let x = self.u[self.i97].wrapping_sub(self.u[self.j97]) & MASK_24;
        self.u[self.i97] = x;
        self.i97 = if self.i97 == 0 {
            LEN_U - 1
        } else {
            self.i97 - 1
        };
        self.j97 = if self.j97 == 0 {
            LEN_U - 1
        } else {
            self.j97 - 1
        };

        // c(n) = c(n-1) - cd mod cm; the modulus is not a power of two, so it is
        // reduced explicitly.
        self.c = if self.c >= C_DELTA {
            self.c - C_DELTA
        } else {
            self.c + (C_MODULUS - C_DELTA)
        };

        x.wrapping_sub(self.c) & MASK_24
    }

    /// Advances the generator and returns the next value in [0, 1). The result is an exact
    /// multiple of 2^-24; `0.0` may occur, `1.0` cannot.
    ///
    /// # Example
    ///
    /// ```rust
    /// use unirand::MarsagliaUniRng;
    ///
    /// let mut rng = MarsagliaUniRng::new(170);
    /// println!("Random number: {}", rng.uni());
    /// ```
    #[inline]
    pub fn uni(&mut self) -> f32 {
        // Exact: a 24-bit integer is representable in f32, and scaling by a power of two
        // introduces no rounding.
        self.uni_u24() as f32 * SCALE_24
    }

    /// Advances the generator and returns the next value in [0, 1) as an `f64`.
    ///
    /// The value is identical to that returned by [`uni`](Self::uni), converted exactly to
    /// `f64`; it remains an exact multiple of 2^-24 and therefore has 24-bit resolution, not
    /// the 53-bit resolution of the `f64` format. `0.0` may occur, `1.0` cannot.
    ///
    /// # Example
    ///
    /// ```rust
    /// use unirand::MarsagliaUniRng;
    ///
    /// let mut first = MarsagliaUniRng::new(170);
    /// let mut second = first.clone();
    /// assert_eq!(first.uni_f64(), f64::from(second.uni()));
    /// ```
    #[inline]
    pub fn uni_f64(&mut self) -> f64 {
        f64::from(self.uni())
    }
}

/// Implements `TryRng` from the `rand_core` crate with an infallible error type, which
/// provides `rand_core::Rng` through its blanket implementation and so enables use with the
/// broader Rust random number ecosystem (distributions, shuffling, sampling, etc.).
///
/// Each call to the underlying generator yields 24 random bits. `next_u32` concatenates the
/// bits of one value with the 8 most significant bits of the next; `next_u64` concatenates
/// 24, 24 and 16 bits from three values; `fill_bytes` takes 3 bytes from each value. All
/// output bits are therefore fully random.
impl rand_core::TryRng for MarsagliaUniRng {
    type Error = rand_core::Infallible;

    #[inline]
    fn try_next_u32(&mut self) -> Result<u32, Self::Error> {
        let high = self.uni_u24();
        let low = self.uni_u24();
        Ok((high << 8) | (low >> 16))
    }

    #[inline]
    fn try_next_u64(&mut self) -> Result<u64, Self::Error> {
        let high = u64::from(self.uni_u24());
        let middle = u64::from(self.uni_u24());
        let low = u64::from(self.uni_u24());
        Ok((high << 40) | (middle << 16) | (low >> 8))
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), Self::Error> {
        let mut chunks = dest.chunks_exact_mut(3);
        for chunk in &mut chunks {
            chunk.copy_from_slice(&self.uni_u24().to_be_bytes()[1..]);
        }
        let remainder = chunks.into_remainder();
        if !remainder.is_empty() {
            let bytes = self.uni_u24().to_be_bytes();
            remainder.copy_from_slice(&bytes[1..=remainder.len()]);
        }
        Ok(())
    }
}

/// Implements `SeedableRng`, so that the generator may be constructed with
/// `seed_from_u64`, `from_rng` and related methods.
///
/// The 4-byte seed is read as a little-endian `u32` and reduced modulo `MAX_SEED + 1`.
impl rand_core::SeedableRng for MarsagliaUniRng {
    type Seed = [u8; 4];

    fn from_seed(seed: Self::Seed) -> Self {
        Self::new(u32::from_le_bytes(seed) % (Self::MAX_SEED + 1))
    }
}

/// Produces an infinite sequence of uniform random values in [0, 1).
///
/// Adaptors such as `take` consume the generator; use `by_ref()` to retain it, e.g.
/// `rng.by_ref().take(5)`.
impl Iterator for MarsagliaUniRng {
    type Item = f32;

    #[inline]
    fn next(&mut self) -> Option<f32> {
        Some(self.uni())
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (usize::MAX, None)
    }
}

impl FusedIterator for MarsagliaUniRng {}

/// Compiles and runs the code examples in the README as doctests.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;

#[cfg(test)]
mod tests {
    extern crate std;

    use super::{MarsagliaUniRng, SeedError};
    use rand_core::{Rng, SeedableRng};
    use std::vec::Vec;

    /// Single seed equivalent to the four seeds (12, 34, 56, 78) of the published test.
    const REFERENCE_SEED: u32 = 1802 * 30082 + 9373;

    /// Published reference values: after 20000 calls from seeds (12, 34, 56, 78), the
    /// next six outputs multiplied by 2^24.
    const REFERENCE_VALUES: [u32; 6] = [
        6_533_892, 14_220_222, 7_275_067, 6_172_232, 8_354_498, 10_633_180,
    ];

    #[test]
    fn test_published_reference_values() {
        let mut rng = MarsagliaUniRng::from_seeds(12, 34, 56, 78).unwrap();
        for _ in 0..20_000 {
            rng.uni();
        }
        for &expected in &REFERENCE_VALUES {
            assert_eq!(rng.uni() * 16_777_216.0, expected as f32);
        }
    }

    #[test]
    fn test_single_seed_decomposition() {
        assert_eq!(
            MarsagliaUniRng::new(REFERENCE_SEED),
            MarsagliaUniRng::from_seeds(12, 34, 56, 78).unwrap()
        );
    }

    #[test]
    fn test_uni_u24_matches_uni() {
        let mut rng_int = MarsagliaUniRng::new(REFERENCE_SEED);
        let mut rng_float = MarsagliaUniRng::new(REFERENCE_SEED);
        for _ in 0..1_000 {
            assert_eq!(rng_int.uni_u24() as f32, rng_float.uni() * 16_777_216.0);
        }
    }

    /// `uni_f64` must return the value of `uni` converted exactly, so that multiplying by
    /// 2^24 recovers the integer form.
    #[test]
    fn test_uni_f64_matches_uni() {
        let mut rng_f64 = MarsagliaUniRng::new(REFERENCE_SEED);
        let mut rng_f32 = MarsagliaUniRng::new(REFERENCE_SEED);
        let mut rng_int = MarsagliaUniRng::new(REFERENCE_SEED);
        for _ in 0..1_000 {
            let value = rng_f64.uni_f64();
            assert_eq!(value, f64::from(rng_f32.uni()));
            assert_eq!(value * 16_777_216.0, f64::from(rng_int.uni_u24()));
        }
    }

    /// Regression value for seed 170 from previous releases.
    #[test]
    fn test_rng_output() {
        let mut rng = MarsagliaUniRng::new(170);
        assert_eq!(rng.uni(), 0.687_533_44);
    }

    #[test]
    #[should_panic(expected = "seed 900000001 out of range")]
    fn test_new_panics_high_seed() {
        let _ = MarsagliaUniRng::new(900_000_001);
    }

    #[test]
    fn test_try_new_rejects_high_seed() {
        assert_eq!(
            MarsagliaUniRng::try_new(900_000_001),
            Err(SeedError::SeedOutOfRange(900_000_001))
        );
    }

    #[test]
    fn test_from_seeds_validation() {
        assert!(MarsagliaUniRng::from_seeds(1, 1, 1, 0).is_err());
        assert!(MarsagliaUniRng::from_seeds(0, 2, 3, 4).is_err());
        assert!(MarsagliaUniRng::from_seeds(2, 179, 3, 4).is_err());
        assert!(MarsagliaUniRng::from_seeds(2, 3, 4, 169).is_err());
        assert!(MarsagliaUniRng::from_seeds(1, 1, 2, 0).is_ok());
        assert!(MarsagliaUniRng::from_seeds(178, 178, 178, 168).is_ok());
    }

    #[test]
    fn test_rng_boundary_seeds() {
        for seed in [0, MarsagliaUniRng::MAX_SEED] {
            let v = MarsagliaUniRng::new(seed).uni();
            assert!((0.0..1.0).contains(&v), "Seed {seed} produced {v}");
        }
    }

    #[test]
    fn test_rng_statistics() {
        let n = 10_000;
        let sum: f64 = MarsagliaUniRng::new(170).take(n).map(f64::from).sum();
        let mean = sum / n as f64;
        assert!(
            (mean - 0.5).abs() < 0.01,
            "Mean out of expected range: {mean}"
        );
    }

    /// Checks that the variance is close to 1/12, as expected for a uniform distribution.
    #[test]
    fn test_rng_variance() {
        let values: Vec<f64> = MarsagliaUniRng::new(170)
            .take(10_000)
            .map(f64::from)
            .collect();
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let variance = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
        assert!(
            (variance - 1.0 / 12.0).abs() < 0.005,
            "Variance out of expected range: {variance}"
        );
    }

    /// Exercises several cycles of both lag indices, including their wrap-around.
    #[test]
    fn test_rng_output_range() {
        for (step, v) in MarsagliaUniRng::new(170).take(10_000).enumerate() {
            assert!(
                (0.0..1.0).contains(&v),
                "Value out of range at step {step}: {v}"
            );
        }
    }

    #[test]
    fn test_rng_reinitialisation() {
        let mut rng = MarsagliaUniRng::new(42);
        let first_run: Vec<f32> = rng.by_ref().take(50).collect();
        rng.rinit(42);
        let second_run: Vec<f32> = rng.take(50).collect();
        assert_eq!(first_run, second_run);
    }

    #[test]
    fn test_clone_continues_identically() {
        let mut rng = MarsagliaUniRng::new(42);
        rng.by_ref().take(123).for_each(drop);
        let mut copy = rng.clone();
        for _ in 0..1_000 {
            assert_eq!(rng.uni(), copy.uni());
        }
    }

    #[test]
    fn test_rng_distinct_seeds_diverge() {
        let first: Vec<f32> = MarsagliaUniRng::new(1).take(50).collect();
        let second: Vec<f32> = MarsagliaUniRng::new(2).take(50).collect();
        assert_ne!(first, second);
    }

    #[test]
    fn test_next_u32_composition() {
        let mut rng = MarsagliaUniRng::new(42);
        let mut reference = MarsagliaUniRng::new(42);
        for _ in 0..100 {
            let high = reference.uni_u24();
            let low = reference.uni_u24();
            assert_eq!(rng.next_u32(), (high << 8) | (low >> 16));
        }
    }

    #[test]
    fn test_next_u64_composition() {
        let mut rng = MarsagliaUniRng::new(42);
        let mut reference = MarsagliaUniRng::new(42);
        for _ in 0..100 {
            let high = u64::from(reference.uni_u24());
            let middle = u64::from(reference.uni_u24());
            let low = u64::from(reference.uni_u24());
            assert_eq!(rng.next_u64(), (high << 40) | (middle << 16) | (low >> 8));
        }
    }

    /// Every bit position of `next_u32`, `next_u64` and every byte position of
    /// `fill_bytes` must take both values; previous releases left the low 8 bits zero.
    #[test]
    fn test_all_output_bits_vary() {
        let mut rng = MarsagliaUniRng::new(170);
        let (mut or32, mut and32) = (0u32, u32::MAX);
        let (mut or64, mut and64) = (0u64, u64::MAX);
        for _ in 0..1_000 {
            let value32 = rng.next_u32();
            let value64 = rng.next_u64();
            or32 |= value32;
            and32 &= value32;
            or64 |= value64;
            and64 &= value64;
        }
        assert_eq!((or32, and32), (u32::MAX, 0));
        assert_eq!((or64, and64), (u64::MAX, 0));

        let mut buffer = [0u8; 4_000];
        rng.fill_bytes(&mut buffer);
        for position in 0..8 {
            assert!(
                buffer
                    .iter()
                    .skip(position)
                    .step_by(8)
                    .any(|&byte| byte != 0),
                "Byte position {position} is always zero"
            );
        }
    }

    /// The correction sequence must satisfy the closed form of the paper,
    /// `c(n) = (362436 - n * 7654321) mod 16777213`, over one full period. A single period
    /// is required because each boundary case of the modular reduction, such as
    /// `c(n-1) = 7654321`, occurs only once per period.
    #[test]
    fn test_correction_sequence_full_period() {
        use super::{C_DELTA, C_INITIAL, C_MODULUS};

        let mut rng = MarsagliaUniRng::new(0);
        for step in 1..=i64::from(C_MODULUS) {
            rng.uni_u24();
            let expected =
                (i64::from(C_INITIAL) - step * i64::from(C_DELTA)).rem_euclid(i64::from(C_MODULUS));
            assert_eq!(i64::from(rng.c), expected, "step {step}");
        }
    }

    #[test]
    fn test_seedable_rng() {
        let seed = 170u32.to_le_bytes();
        assert_eq!(MarsagliaUniRng::from_seed(seed), MarsagliaUniRng::new(170));
        assert_eq!(
            MarsagliaUniRng::seed_from_u64(7),
            MarsagliaUniRng::seed_from_u64(7)
        );
        // Out-of-range seeds are reduced rather than rejected.
        MarsagliaUniRng::from_seed(u32::MAX.to_le_bytes()).uni();
    }

    #[test]
    fn test_iterator_size_hint() {
        let rng = MarsagliaUniRng::new(1);
        assert_eq!(rng.size_hint(), (usize::MAX, None));
    }

    /// Returns the 64-bit FNV-1a hash of the bit patterns of the first `count` values of
    /// `uni()` for the given seed.
    fn hash_uni_sequence(seed: u32, count: usize) -> u64 {
        MarsagliaUniRng::new(seed)
            .take(count)
            .fold(0xcbf2_9ce4_8422_2325, |hash, value| {
                (hash ^ u64::from(value.to_bits())).wrapping_mul(0x0000_0100_0000_01b3)
            })
    }

    /// Long-run regression test: the first 1,000,000 values of `uni()` for each seed must be
    /// bit-identical to those of the `f32` implementation released as 0.2.0, which was
    /// verified against the published reference values.
    #[test]
    fn test_uni_sequence_regression() {
        let reference_hashes: [(u32, u64); 7] = [
            (0, 0x6c47_c3ed_f4b0_7ab2),
            (1, 0xa49f_6f5b_81a4_7b93),
            (42, 0x334e_f419_6d4c_2eb2),
            (170, 0xd7c6_66bf_3c27_d0f3),
            (REFERENCE_SEED, 0xffda_35af_ee83_4e7d),
            (123_456_789, 0x5a9d_f787_4864_c054),
            (MarsagliaUniRng::MAX_SEED, 0xb7d0_937a_eee2_7e8a),
        ];
        for (seed, expected) in reference_hashes {
            assert_eq!(
                hash_uni_sequence(seed, 1_000_000),
                expected,
                "Sequence for seed {seed} differs from 0.2.0"
            );
        }
    }

    /// Every seed in `0..=MAX_SEED` must decompose into valid seeds. `i` and `j` depend
    /// only on `ij = seed / 30082`, and `k` and `l` only on `kl = seed mod 30082`, so the
    /// range is covered by checking every attainable `ij` and every attainable `kl`.
    /// Requiring `i, j >= 2` also excludes the forbidden case `i = j = k = 1`.
    #[test]
    fn test_every_seed_decomposes_validly() {
        for ij in 0..=MarsagliaUniRng::MAX_SEED / 30082 {
            let (i, j, _, _) = super::decompose_seed(ij * 30082);
            assert!(
                (2..=178).contains(&i) && (2..=178).contains(&j),
                "ij = {ij} gives i = {i}, j = {j}"
            );
        }
        for kl in 0..30082 {
            let (_, _, k, l) = super::decompose_seed(kl);
            assert!(
                (1..=178).contains(&k) && l <= 168,
                "kl = {kl} gives k = {k}, l = {l}"
            );
            assert!(super::seeds_are_valid(2, 2, k, l));
        }
        assert!(MarsagliaUniRng::try_new(MarsagliaUniRng::MAX_SEED).is_ok());
    }

    /// `fill_bytes` must reproduce the documented layout for every remainder length and
    /// consume exactly `ceil(n / 3)` values for a buffer of `n` bytes.
    #[test]
    fn test_fill_bytes_all_lengths() {
        for length in 0..=10 {
            let mut rng = MarsagliaUniRng::new(42);
            let mut reference = MarsagliaUniRng::new(42);
            let mut buffer = [0u8; 10];
            rng.fill_bytes(&mut buffer[..length]);

            let draws = length.div_ceil(3);
            let expected: Vec<u8> = (0..draws)
                .flat_map(|_| reference.uni_u24().to_be_bytes()[1..].to_vec())
                .take(length)
                .collect();
            assert_eq!(&buffer[..length], expected.as_slice(), "length {length}");
            assert!(buffer[length..].iter().all(|&byte| byte == 0));
            assert_eq!(
                rng, reference,
                "length {length} consumed the wrong number of values"
            );
        }
    }

    /// Returns the chi-squared statistic of byte counts against a uniform distribution.
    fn chi_squared_bytes(bytes: &[u8]) -> f64 {
        let mut counts = [0u64; 256];
        for &byte in bytes {
            counts[usize::from(byte)] += 1;
        }
        let expected = bytes.len() as f64 / 256.0;
        counts
            .iter()
            .map(|&count| (count as f64 - expected).powi(2) / expected)
            .sum()
    }

    /// Critical value of the chi-squared distribution with 255 degrees of freedom at the
    /// 0.001 significance level.
    const CHI_SQUARED_255_CRITICAL: f64 = 330.5;

    /// Chi-squared tests of byte uniformity through the `rand` API, which previously
    /// returned 0 for every `u8` and zero for every fourth byte of `fill`. The seed is fixed,
    /// so the test is deterministic.
    #[test]
    fn test_rand_byte_uniformity() {
        use rand::RngExt;

        let sample_count = 256 * 1_000;
        let mut rng = MarsagliaUniRng::new(170);

        let sampled: Vec<u8> = (0..sample_count).map(|_| rng.random::<u8>()).collect();
        let statistic = chi_squared_bytes(&sampled);
        assert!(
            statistic < CHI_SQUARED_255_CRITICAL,
            "random::<u8>() chi-squared = {statistic}"
        );

        let mut filled = std::vec![0u8; sample_count];
        rng.fill(filled.as_mut_slice());
        let statistic = chi_squared_bytes(&filled);
        assert!(
            statistic < CHI_SQUARED_255_CRITICAL,
            "fill() chi-squared = {statistic}"
        );
    }

    #[test]
    fn test_seed_error_messages() {
        assert_eq!(
            std::format!("{}", SeedError::SeedOutOfRange(900_000_001)),
            "seed 900000001 out of range 0..=900000000"
        );
        assert_eq!(
            std::format!("{}", MarsagliaUniRng::from_seeds(1, 1, 1, 0).unwrap_err()),
            "seeds (i, j, k, l) = (1, 1, 1, 0) invalid: i, j, k must lie in 1..=178 and not \
             all equal 1, l must lie in 0..=168"
        );
    }
}
