extern crate test;

mod xxhash3;
const ALPHA1: f64 = 0.7213;
const ALPHA2: f64 = 1.079;

pub struct LogLogBeta {
    precision_bits: u64,
    max_precision_bits: u64,
    max_x: u64,
    register_count: u64,
    alpha: f64,
    registers: Vec<u8>,
}

impl LogLogBeta {
    pub fn new(error_rate: f64) -> LogLogBeta {
        let precision_bits = (1.04 / error_rate).powi(2).log2().ceil() as u64;
        let count = 1u64 << precision_bits;
        LogLogBeta {
            precision_bits,
            max_precision_bits: 64 - precision_bits,
            max_x: std::u64::MAX >> (64 - precision_bits),
            register_count: count,
            alpha: ALPHA1 / (1.0 + (ALPHA2 / (count as f64))),
            registers: vec![0; count as usize],
        }
    }

    pub fn add_hash(&mut self, hash: u64) {
        let tmp = (hash << self.precision_bits) ^ self.max_x;
        let val = (tmp.leading_zeros() + 1) as u8;

        let k = hash >> self.max_precision_bits;
        if self.registers[k as usize] < val {
            self.registers[k as usize] = val as u8;
        }
    }

    pub fn add(&mut self, bytes: &[u8]) {
        self.add_hash(xxhash3::hash_64(bytes))
    }

    fn sum_registers(&self) -> (f64, f64) {
        let mut zero_count = 0.0;
        let mut sum = 0.0;

        for val in self.registers.iter() {
            let x = *val;
            if x == 0u8 {
                zero_count += 1.0;
            }

            sum += 1.0 / 2.0f64.powf(x as f64);
        }

        (sum, zero_count)
    }

    pub fn beta(&self, ez: f64) -> f64 {
        let zl = (ez + 1.0f64).ln();
        -0.370393911 * ez
            + 0.070471823 * zl
            + 0.17393686 * (zl.powi(2))
            + 0.16339839 * (zl.powi(3))
            + -0.09237745 * (zl.powi(4))
            + 0.03738027 * (zl.powi(5))
            + -0.005384159 * (zl.powi(6))
            + 0.00042419 * (zl.powi(7))
    }

    pub fn cardinality(&self) -> u64 {
        let (sum, zero_count) = self.sum_registers();
        let reg_count = self.register_count as f64;

        ((self.alpha * reg_count * ((reg_count - zero_count) as f64))
            / (self.beta(zero_count as f64) + sum)) as u64
    }

    pub fn merge(&mut self, merge_me: &LogLogBeta) {
        for register_ix in 0..merge_me.register_count {
            self.registers[register_ix as usize] = std::cmp::max(
                self.registers[register_ix as usize],
                merge_me.registers[register_ix as usize],
            );
        }
    }
}

#[cfg(test)]
mod tests {
    extern crate rand;
    use super::*;
    use rand::distributions::Alphanumeric;
    use rand::{thread_rng, Rng};

    fn estimate_error(actual: i64, expected: i64) -> f64 {
        let delta = (expected - actual).abs();
        (delta as f64) / (expected as f64)
    }

    fn gen_string() -> String {
        let mut rng = thread_rng();

        let str_length: u64 = (rng.gen::<u64>() % 256) + 1;
        rng.sample_iter(&Alphanumeric)
            .take(str_length as usize)
            .collect()
    }

    #[test]
    fn add_one() {
        let mut ll = LogLogBeta::new(0.01);
        ll.add("hello world".as_bytes());
    }

    #[test]
    fn check_error_rate() {
        let mut count = 100u64;

        while count <= 1000000 {
            let mut ll = LogLogBeta::new(0.01);
            for _ in 0..count {
                ll.add(gen_string().as_bytes());
            }
            let est = ll.cardinality();
            let error_rate = estimate_error(est as i64, count as i64);
            println!("actual {} expected {}: rate {}\n", est, count, error_rate);
            assert_eq!(true, error_rate < 0.02);
            count *= 10;
        }
    }

    #[test]
    fn add_1000() {
        let count = 1000u64;

        let mut ll = LogLogBeta::new(0.01);
        for _ in 0..count {
            ll.add(gen_string().as_bytes());
        }
        let est = ll.cardinality();
        let error_rate = estimate_error(est as i64, count as i64);
        println!("actual {} expected {}: rate {}\n", est, count, error_rate);
        assert_eq!(true, error_rate < 0.02);
    }

    #[test]
    fn merge() {
        let count = 1000u64;

        let mut ll1 = LogLogBeta::new(0.01);
        let mut ll2 = LogLogBeta::new(0.01);
        for _ in 0..count {
            ll1.add(gen_string().as_bytes());
            ll2.add(gen_string().as_bytes());
        }

        ll1.merge(&ll2);
        let est = ll1.cardinality();
        let error_rate = estimate_error(est as i64, 2 * count as i64);
        println!("actual {} expected {}: rate {}\n", est, count, error_rate);
        assert_eq!(true, error_rate < 0.02);
    }

    use test::{black_box, Bencher};

    #[bench]
    fn add(b: &mut Bencher) {
        let mut ll = LogLogBeta::new(0.01);
        b.iter(|| black_box(ll.add_hash(9394790207897)));
    }

    #[bench]
    fn cardinality(b: &mut Bencher) {
        let mut ll = LogLogBeta::new(0.01);
        for _ in 0..1000 {
            ll.add(gen_string().as_bytes());
        }
        b.iter(|| black_box(ll.cardinality()));
    }

    #[bench]
    fn perf_merge(b: &mut Bencher) {
        let mut ll1 = LogLogBeta::new(0.01);
        let mut ll2 = LogLogBeta::new(0.01);
        for _ in 0..1000 {
            ll1.add(gen_string().as_bytes());
            ll2.add(gen_string().as_bytes());
        }
        b.iter(|| black_box(ll1.merge(&ll2)));
    }

    #[bench]
    fn full(b: &mut Bencher) {
        let mut ll = LogLogBeta::new(0.01);
        let count = 100000;
        let mut hashes = vec![0u64; count];
        for i in 0..count {
            hashes[i] = xxhash3::hash_64(gen_string().as_bytes())
        }

        b.iter(|| {
            for hash in &hashes {
                ll.add_hash(*hash)
            }
            black_box(ll.cardinality())
        });
    }
}
