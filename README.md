# unirand
A Rust implementation of Marsaglia's Universal Random Number Generator

This program is based on "Toward a universal random number generator" 
by George Marsaglia, Arif Zaman, Wai Wan Tsang 
published in Statistics & Probability Letters Volume 9, Issue 1, January 1990, Pages 35-39
https://www.sciencedirect.com/science/article/abs/pii/016771529090092L

The RNG uses a sequence of operations to generate uniformly distributed
random numbers between 0 and 1. It has been designed with simplicity and 
reproducibility in mind.

# Security Warning
This RNG is not suitable for cryptographic or security-sensitive applications.
It is designed for statistical simulation and makes no guarantees of
unpredictability or resistance to state recovery.

# Usage

`Cargo.toml`  
```toml
[dependencies]
unirand = "0.2.0"
```

## Basic usage

Initialise with a seed and generate values directly via `uni()`:

```rust
use unirand::MarsagliaUniRng;

fn main() {
    let mut rng = MarsagliaUniRng::new();
    rng.rinit(170);
    for _ in 0..5 {
        println!("Random number: {}", rng.uni());
    }
}
```

## Iterator usage

`MarsagliaUniRng` implements `Iterator<Item = f32>`, giving access to the full
Rust iterator adaptor API:

```rust
use unirand::MarsagliaUniRng;

fn main() {
    let mut rng = MarsagliaUniRng::new();
    rng.rinit(170);

    // Collect 5 values.
    let values: Vec<f32> = rng.take(5).collect();
    println!("{:?}", values);

    // Sum 1000 values.
    rng.rinit(170);
    let sum: f32 = rng.take(1_000).sum();
    println!("Sum: {}", sum);

    // Filter values above 0.9.
    rng.rinit(170);
    let high: Vec<f32> = rng.take(10_000).filter(|&x| x > 0.9).collect();
    println!("Values above 0.9: {}", high.len());
}
```

## rand_core usage

`MarsagliaUniRng` implements `rand_core::RngCore`, enabling use with the broader
Rust random number ecosystem. Add `rand` to your dependencies to access
distributions, shuffling, and sampling:

`Cargo.toml`
```toml
[dependencies]
unirand = "0.2.0"
rand = "0.9"
```

```rust
use unirand::MarsagliaUniRng;
use rand::prelude::*;

fn main() {
    let mut rng = MarsagliaUniRng::new();
    rng.rinit(170);

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

Note: `next_u32` has full entropy in the upper 24 bits only, reflecting the f32
precision of the underlying generator. The lower 8 bits are always zero.

# Change Log

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
