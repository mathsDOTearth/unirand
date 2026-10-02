# unirand
A Rust implementation of Marsaglia's Universal Random Number Generator

This program is based on "Toward a universal random number generator" 
by George Marsaglia, Arif Zaman, Wai Wan Tsang 
published in Statistics & Probability Letters Volume 9, Issue 1, January 1990, Pages 35-39
https://www.sciencedirect.com/science/article/abs/pii/016771529090092L

The generator combines a lagged Fibonacci generator with an arithmetic sequence
to produce uniformly distributed values in [0, 1), each an exact multiple of
2^-24, with a period of approximately 2^144. The state is held as exact 24-bit
integers, so the output is bit-for-bit identical to the reference implementation
of the paper on every platform. The crate is `no_std`.

# Security Warning
This RNG is not suitable for cryptographic or security-sensitive applications.
It is designed for statistical simulation and makes no guarantees of
unpredictability or resistance to state recovery.

# Usage

`Cargo.toml`  
```toml
[dependencies]
unirand = "0.3.0"
```

## Basic usage

Create a generator from a seed in `0..=900_000_000` and generate values via `uni()`:

```rust
use unirand::MarsagliaUniRng;

fn main() {
    let mut rng = MarsagliaUniRng::new(170);
    for _ in 0..5 {
        println!("Random number: {}", rng.uni());
    }
}
```

`MarsagliaUniRng::try_new` returns a `Result` instead of panicking on an
out-of-range seed. `uni_u24()` returns the same values as 24-bit integers
(`uni() * 2^24`).

## Seeding with the four seeds of the paper

The paper initialises the generator from four seeds `i, j, k` in `1..=178`
(not all 1) and `l` in `0..=168`. These may be supplied directly:

```rust
use unirand::MarsagliaUniRng;

fn main() {
    // Published test: after 20000 values, the next six multiplied by 2^24 are
    // 6533892, 14220222, 7275067, 6172232, 8354498, 10633180.
    let mut rng = MarsagliaUniRng::from_seeds(12, 34, 56, 78).unwrap();
    for _ in 0..20_000 {
        rng.uni();
    }
    assert_eq!(rng.uni_u24(), 6533892);
}
```

The single seed of `new` is decomposed into four seeds following F. James's
RANMAR; `MarsagliaUniRng::new(1802 * 30082 + 9373)` is equivalent to the
four seeds `(12, 34, 56, 78)`.

## Iterator usage

`MarsagliaUniRng` implements `Iterator<Item = f32>`, giving access to the full
Rust iterator adaptor API. Adaptors such as `take` consume the generator, so
use `by_ref()` to keep it:

```rust
use unirand::MarsagliaUniRng;

fn main() {
    let mut rng = MarsagliaUniRng::new(170);

    // Collect 5 values.
    let values: Vec<f32> = rng.by_ref().take(5).collect();
    println!("{:?}", values);

    // Sum 1000 values.
    rng.rinit(170);
    let sum: f32 = rng.by_ref().take(1_000).sum();
    println!("Sum: {}", sum);

    // Filter values above 0.9.
    rng.rinit(170);
    let high: Vec<f32> = rng.take(10_000).filter(|&x| x > 0.9).collect();
    println!("Values above 0.9: {}", high.len());
}
```

## rand_core usage

`MarsagliaUniRng` implements the `rand_core` 0.10 traits `Rng` (through
`TryRng`) and `SeedableRng`, enabling use with the broader Rust random number
ecosystem. Add `rand` 0.10 to your dependencies to access distributions,
shuffling, and sampling:

`Cargo.toml`
```toml
[dependencies]
unirand = "0.3.0"
rand = "0.10"
```

```rust
use unirand::MarsagliaUniRng;
use rand::prelude::*;

fn main() {
    let mut rng = MarsagliaUniRng::new(170);

    // Generate a random integer in a range.
    let n: u32 = rng.random_range(1..=100);
    println!("Random integer 1-100: {}", n);

    // Shuffle a vector.
    let mut data = vec![1, 2, 3, 4, 5];
    data.shuffle(&mut rng);
    println!("Shuffled: {:?}", data);

    // Fill a buffer with random bytes.
    let mut buf = [0u8; 8];
    rng.fill(&mut buf);
    println!("Random bytes: {:?}", buf);
}
```

Each value of the underlying generator carries 24 random bits. `next_u32`
combines 24 bits from one value with 8 from the next, `next_u64` combines
24 + 24 + 16 bits from three values, and `fill_bytes` takes 3 bytes per value,
so every output bit is random. The `rand_core` output stream is therefore not a
direct transcription of `uni()`; use `uni()` or `uni_u24()` where the sequence
defined in the paper is required.

# Change Log

## version 0.3.0
This release contains breaking API changes. The sequence produced by `uni()` is
unchanged and bit-for-bit identical to previous releases and to the paper.  
Fixed `next_u32`, which left the low 8 bits zero: `rand`'s `random::<u8>()` always returned 0, `random::<u16>()` had only 8 random bits, and `fill_bytes` set every fourth byte to 0. `next_u32`, `next_u64` and `fill_bytes` now produce fully random bits, so their output differs from 0.2.0  
Fixed the README iterator example, which did not compile; README examples are now compiled and run as doctests  
Construction now requires a seed: `new(seed)`, `try_new(seed) -> Result<_, SeedError>`, or `from_seeds(i, j, k, l)` for the four seeds of the paper. Seeds are `u32`. The uninitialised state, and its runtime check in `uni()`, are removed, as is `Default`  
`rinit(seed)` remains for re-seeding in place  
Updated to `rand_core` 0.10, for use with `rand` 0.10: `RngCore` is replaced by `TryRng` with an infallible error type, which provides `rand_core::Rng`. The `rand` integration no longer works with `rand` 0.9; the remainder of the API is unaffected  
Implemented `rand_core::SeedableRng`  
Added `uni_u24()`, returning each value as an exact 24-bit integer  
State is held as exact 24-bit integers rather than `f32`  
Derived `Clone`, `Debug`, `PartialEq` and `Eq`; implemented `FusedIterator` and an infinite `size_hint`  
The crate is now `no_std`; edition 2024, minimum Rust 1.85  
Added a test against the published reference values  

## version 0.2.0
Implemented `Iterator<Item = f32>` which enables idiomatic use of Rust iterator adaptors  
Implemented `rand_core::RngCore` which integrates with the Rust random number ecosystem  

## version 0.1.4
`rstart` made private so now `rinit` is the sole validated public entry point  
Dead code removed from `rinit` now the unreachable range checks replaced with `debug_assert!`  
Added initialisation guard so `uni()` panics with a clear message if called before `rinit`  

## version 0.1.3
Fixed index wrap-around bug in `uni()`: indices now cycle over 1..=97, preventing access to `uni_u[0]` which was never seeded and held the constant 0.0  
Added tests: output range, wrap-around regression, boundary seeds, re-initialisation, distinct seed divergence, variance  

## version 0.1.2
Corrected code comments and documentation  
Added more documentation  
Added more tests  

## version 0.1.1
Corrected code to pass `cargo clippy`  
Added out of range tests  

## version 0.1.0
Initial version
