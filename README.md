# unirand
A Rust implementation of Marsaglia's Universal Random Number Generator

This program is based on "Toward a universal random number generator" 
by George Marsaglia, Arif Zaman, Wai Wan Tsang 
published in Statistics & Probability Letters Volume 9, Issue 1, January 1990, Pages 35-39
https://www.sciencedirect.com/science/article/abs/pii/016771529090092L?via%3Dihub

The RNG uses a sequence of operations to generate uniformly distributed
random numbers between 0 and 1. It has been designed with simplicity and 
reproducibility in mind.

# SECURITY WARNING
This random number generator is probably not suitable for use in Cryptography or Security 
application as it has not been rigorously tested with that use in mind.

# Usage Example
`Cargo.toml`  
```toml
[dependencies]
unirand = "0.1.2"
```

Then, you can initialise and use the RNG in your project as follows:

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

# Change Log
## version 0.1.2
Corrected code comments and documentation  
Added more documentation  
Added more tests  

## version 0.1.1
Corrected code to pass `cargo clippy`  
Added out of range tests  

## version 0.1.0
Initial version
