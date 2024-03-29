// Rust version of Marsaglia's UNI random number generator
// Ported from C version provided by Edinburough University
// Ported by R. Neale info@maths.earth 2023 under MIT License

const LEN_U: usize = 98; // Constant defining the length of the random values array.

// A struct representing Marsaglia's Universal Random Number Generator.
struct MarsagliaUniRng {
    uni_u: [f32; LEN_U], // Array holding the recent random numbers.
    uni_c: f32, 		// Correction to avoid periodicity.
    uni_cd: f32, 		// Correction delta value.
    uni_cm: f32,		// Correction modulus.
    uni_ui: usize,		// Current position in the random values array.
    uni_uj: usize,
}

impl MarsagliaUniRng {
// Constructor for the random number generator.
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
// Generate a new random float value between 0 and 1
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

// Initialises the random values array using four seeds.
    pub fn rstart(&mut self, i: i32, j: i32, k: i32, l: i32) {
        let mut i = i;
        let mut j = j;
        let mut k = k;
        let mut l = l;
        for ii in 1..=97 {
            let mut s = 0.0;
            let mut t = 0.5;
            for _jj in 1..=24 {
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
// Set fixed correction values
        self.uni_c = 362436.0 / 16777216.0;
        self.uni_cd = 7654321.0 / 16777216.0;
        self.uni_cm = 16777213.0 / 16777216.0;
        self.uni_ui = 97;
        self.uni_uj = 33;
    }

// Validates and decomposes a single seed into four seeds, then initialises the random values array.
    pub fn rinit(&mut self, ijkl: i32) {
        if ijkl < 0 || ijkl > 900_000_000 {
            panic!("rinit: ijkl = {} -- out of range", ijkl);
        }

        let ij = ijkl / 30082;
        let kl = ijkl - (30082 * ij);
        let i = ((ij / 177) % 177) + 2;
        let j = (ij % 177) + 2;
        let k = ((kl / 169) % 178) + 1;
        let l = kl % 169;

        if i <= 0 || i > 178 {
            panic!("rinit: i = {} -- out of range", i);
        }
        if j <= 0 || j > 178 {
            panic!("rinit: j = {} -- out of range", j);
        }
        if k <= 0 || k > 178 {
            panic!("rinit: k = {} -- out of range", k);
        }
        if l < 0 || l > 168 {
            panic!("rinit: l = {} -- out of range", l);
        }
        if i == 1 && j == 1 && k == 1 {
            panic!("rinit: 1 1 1 not allowed for 1st 3 seeds");
        }

        self.rstart(i, j, k, l);
    }
}

// Example of how to use the Marsaglia's UNI random number generator.
fn main() {
    let mut rng = MarsagliaUniRng::new(); 	// Create a new random number generator instance.
    rng.rinit(170);							// Initialise it with a seed.
    println!("Random Number: {}", rng.uni()); // Generate and print a random number.
}
