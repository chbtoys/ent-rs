//! # ent-rs
//!
//! `ent-rs` is a library for analyzing the entropy and randomness of binary data.
//! It provides byte/bit entropy, chi-square testing, mean, Pi estimation, and serial correlation.
//!
//! ```rust
//! use ent_rs::EntStats;
//! let data = b"example data";
//! let stats = EntStats::from_data(data, false);
//! println!("Entropy: {}", stats.entropy);
//! ```

use statrs::function::erf::erfc;
use std::f64::consts::SQRT_2;

/// Result of statistical analysis on binary data.
#[derive(Debug, Clone)]
pub struct EntStats {
    /// Shannon entropy in bits per byte (or bit).
    pub entropy: f64,
    /// Ideal compression percentage based on entropy.
    pub compression_percent: f64,
    /// Chi-square test value.
    pub chisquare: f64,
    /// p-value of chi-square test.
    pub p_value: f64,
    /// Arithmetic mean of all data bytes.
    pub mean: f64,
    /// Estimated value of Pi from Monte Carlo method.
    pub pi_estimate: f64,
    /// Serial correlation coefficient between adjacent values.
    pub serial_correlation: f64,
    /// Byte frequency table: (value, count, fraction).
    pub byte_frequencies: Option<Vec<(u8, usize, f64)>>,
    /// Bit frequency table: [(count, fraction) for 0, 1].
    pub bit_frequencies: Option<[(usize, f64); 2]>,
}

impl EntStats {
    /// Compute entropy statistics from byte slice, using bit mode or byte mode.
    pub fn from_data(data: &[u8], bit_mode: bool) -> Self {
        let entropy = calculate_entropy(data, bit_mode);
        let compression_percent = if bit_mode {
            100.0 * (1.0 - entropy)
        } else {
            100.0 * (1.0 - entropy / 8.0)
        };
        let (chisquare, p_value) = calculate_chisquare(data, bit_mode);
        let mean = calculate_mean(data);
        let pi_estimate = estimate_pi(data);
        let serial_correlation = serial_correlation(data);

        let (byte_frequencies, bit_frequencies) = if bit_mode {
            (None, Some(bit_occurrences(data)))
        } else {
            (Some(byte_occurrences(data)), None)
        };

        EntStats {
            entropy,
            compression_percent,
            chisquare,
            p_value,
            mean,
            pi_estimate,
            serial_correlation,
            byte_frequencies,
            bit_frequencies,
        }
    }
}

// Internal computation functions

fn calculate_entropy(data: &[u8], bit_mode: bool) -> f64 {
    let mut freq = if bit_mode {
        vec![0f64; 2]
    } else {
        vec![0f64; 256]
    };

    if bit_mode {
        for &b in data {
            for i in 0..8 {
                freq[(b >> i) as usize & 1] += 1.0;
            }
        }
        let total = 8.0 * data.len() as f64;
        for f in freq.iter_mut() {
            *f /= total;
        }
    } else {
        for &b in data {
            freq[b as usize] += 1.0;
        }
        let total = data.len() as f64;
        for f in freq.iter_mut() {
            *f /= total;
        }
    }

    freq.iter()
        .filter(|&&p| p > 0.0)
        .map(|&p| -p * p.log2())
        .sum()
}

fn calculate_chisquare(data: &[u8], bit_mode: bool) -> (f64, f64) {
    if bit_mode {
        let mut count = [0usize; 2];
        for &b in data {
            for i in 0..8 {
                count[(b >> i) as usize & 1] += 1;
            }
        }
        let total = data.len() * 8;
        let expected = total as f64 / 2.0;
        let chisq = count
            .iter()
            .map(|&obs| {
                let diff = obs as f64 - expected;
                diff * diff / expected
            })
            .sum::<f64>();
        let z = (chisq - 1.0).sqrt();
        (chisq, 1.0 - 0.5 * erfc(-z / SQRT_2))
    } else {
        let mut count = [0usize; 256];
        for &b in data {
            count[b as usize] += 1;
        }
        let total = data.len();
        let expected = total as f64 / 256.0;
        let chisq = count
            .iter()
            .map(|&obs| {
                let diff = obs as f64 - expected;
                diff * diff / expected
            })
            .sum::<f64>();
        let z = (chisq - 255.0).sqrt();
        (chisq, 1.0 - 0.5 * erfc(-z / SQRT_2))
    }
}

fn calculate_mean(data: &[u8]) -> f64 {
    data.iter().map(|&b| b as f64).sum::<f64>() / data.len() as f64
}

fn estimate_pi(data: &[u8]) -> f64 {
    let mut hits = 0;
    let mut total = 0;
    let r_sq = 1u64 << 48;

    for chunk in data.chunks_exact(6) {
        let x = ((chunk[0] as u64) << 16) | ((chunk[1] as u64) << 8) | chunk[2] as u64;
        let y = ((chunk[3] as u64) << 16) | ((chunk[4] as u64) << 8) | chunk[5] as u64;
        let dist_sq = x * x + y * y;
        if dist_sq < r_sq {
            hits += 1;
        }
        total += 1;
    }

    if total > 0 {
        4.0 * hits as f64 / total as f64
    } else {
        0.0
    }
}

fn serial_correlation(data: &[u8]) -> f64 {
    if data.len() < 2 {
        return -99999.0;
    }

    let mut sum_x = 0f64;
    let mut sum_y = 0f64;
    let mut sum_xy = 0f64;
    let mut sum_x2 = 0f64;
    let mut sum_y2 = 0f64;

    for i in 1..data.len() {
        let x = data[i - 1] as f64;
        let y = data[i] as f64;
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
        sum_y2 += y * y;
    }

    let n = (data.len() - 1) as f64;
    let num = n * sum_xy - sum_x * sum_y;
    let denom = ((n * sum_x2 - sum_x.powi(2)) * (n * sum_y2 - sum_y.powi(2))).sqrt();

    if denom == 0.0 {
        -99999.0
    } else {
        num / denom
    }
}

fn byte_occurrences(data: &[u8]) -> Vec<(u8, usize, f64)> {
    let mut counts = [0usize; 256];
    for &b in data {
        counts[b as usize] += 1;
    }
    let total = data.len() as f64;
    (0..=255)
        .map(|i| (i as u8, counts[i], counts[i] as f64 / total))
        .collect()
}

fn bit_occurrences(data: &[u8]) -> [(usize, f64); 2] {
    let mut count = [0usize; 2];
    for &b in data {
        for i in 0..8 {
            count[(b >> i) as usize & 1] += 1;
        }
    }
    let total = (data.len() * 8) as f64;
    [
        (count[0], count[0] as f64 / total),
        (count[1], count[1] as f64 / total),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_entropy_on_uniform_data() {
        let data = vec![0xAA; 4096]; // Constant pattern
        let stats = EntStats::from_data(&data, false);
        assert!(stats.entropy < 1.0, "Expected low entropy for uniform data");
        assert!(stats.compression_percent > 85.0);
    }

    #[test]
    fn test_entropy_on_random_data() {
        let data: Vec<u8> = (0..=255).cycle().take(4096).collect(); // Pseudo-random-like distribution
        let stats = EntStats::from_data(&data, false);
        assert!(
            stats.entropy > 7.5,
            "Expected high entropy for diverse data"
        );
        assert!(stats.compression_percent < 10.0);
    }

    #[test]
    fn test_chisquare_and_pvalue_validity() {
        // 75% zeros, 25% uniform noise
        let mut data = vec![0u8; 3072];
        data.extend((0..=255).cycle().take(1024));
        let stats = EntStats::from_data(&data, false);
        assert!(
            stats.chisquare > 0.0,
            "Chi-square should be > 0 for biased input"
        );
        assert!(
            (0.0..=1.0).contains(&stats.p_value),
            "p-value should be in [0, 1]"
        );
    }

    #[test]
    fn test_mean_value_byte_mode() {
        let data = vec![0x00, 0xFF];
        let stats = EntStats::from_data(&data, false);
        assert!(
            (stats.mean - 127.5).abs() < 1.0,
            "Mean should be close to 127.5"
        );
    }

    #[test]
    fn test_pi_estimation_sanity() {
        let data: Vec<u8> = (0..=255).cycle().take(8192).collect(); // Somewhat randomized
        let stats = EntStats::from_data(&data, false);
        assert!(
            (2.5..=3.8).contains(&stats.pi_estimate),
            "Pi estimate should be within realistic bounds"
        );
    }

    #[test]
    fn test_serial_correlation_constant() {
        let data = vec![0x33; 2048];
        let stats = EntStats::from_data(&data, false);
        assert_eq!(
            stats.serial_correlation, -99999.0,
            "All values are equal, correlation should be undefined"
        );
    }

    #[test]
    fn test_bit_mode_entropy_and_frequencies() {
        let data = vec![0b10101010u8; 1024]; // Equal 1s and 0s
        let stats = EntStats::from_data(&data, true);
        assert!(
            (stats.entropy - 1.0).abs() < 0.01,
            "Expected full bit entropy"
        );
        assert!(stats.bit_frequencies.is_some());
        let freqs = stats.bit_frequencies.unwrap();
        assert!((freqs[0].1 - 0.5).abs() < 0.01);
        assert!((freqs[1].1 - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_byte_frequency_distribution_length() {
        let data: Vec<u8> = (0..=255).cycle().take(4096).collect();
        let stats = EntStats::from_data(&data, false);
        let freqs = stats.byte_frequencies.as_ref().unwrap();
        assert_eq!(freqs.len(), 256);
    }
}
