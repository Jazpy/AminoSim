use rand::rngs::ThreadRng;
use rand::distributions::{Uniform, Distribution};

#[derive(Clone)]
pub struct Sequence {
    pub nucleotides: Vec<u8>,
    size: usize,
    freq_table: Vec<(u8, f64)>,
    max_freq: f64
}

fn get_cumulative(t: &Vec<(u8, f64)>) -> f64 {
    let mut cumulative_freq: f64 = 0.0;

    // Build with cumulative values
    assert!(t.len() > 0, "Empty frequency table");
    for &(_, f) in t.iter() {
        // Validate values in table
        assert!(f > 0.0, "Can't have nucleotide frequencies <= 0");

        cumulative_freq += f;
    }

    cumulative_freq
}

impl Sequence {
    pub fn new(t: &Vec<(u8, f64)>, l: usize) -> Sequence {
        let cumulative_freq = get_cumulative(t);

        // Build our empty sequence
        let mut ret = Sequence {
            nucleotides: Vec::<u8>::new(),
            size: 0,
            freq_table: t.clone(),
            max_freq: cumulative_freq
        };

        // Append 'l' nucleotides to our sequence
        ret.append(l);
        ret
    }

    pub fn from_vec(s: Vec<u8>, t: &Vec<(u8, f64)>) -> Sequence {
        let cumulative_freq = get_cumulative(t);

        // Attach given vec to our Sequence object
        let len = s.len();
        Sequence {
            nucleotides: s,
            size: len,
            freq_table: t.clone(),
            max_freq: cumulative_freq
        }
    }

    fn sample(&self, generator: Uniform<f64>, mut rng: ThreadRng) -> u8 {
        let mut r: f64 = generator.sample(&mut rng);

        for &(c, f) in self.freq_table.iter() {
            if r < f {
                return c
            }

            r -= f;
        }

        assert!(false, "Something went terribly wrong in Sequence's sampler");
        return 0
    }

    pub fn append(&mut self, l: usize) {
        // Initialize RNG
        let rng = rand::thread_rng();
        let generator = Uniform::from(0.0..self.max_freq);

        for _ in 0..l {
            self.nucleotides.push(self.sample(generator, rng));
        }

        self.size += l;
    }

    #[allow(dead_code)]
    pub fn print(&self) {
        unsafe {
            println!("{}", std::str::from_utf8_unchecked(&self.nucleotides));
        }
    }

    #[allow(dead_code)]
    pub fn to_string(&self) -> &str {
        unsafe {
            return std::str::from_utf8_unchecked(&self.nucleotides);
        }
    }
}
