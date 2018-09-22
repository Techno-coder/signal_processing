use super::Hertz;
use super::SampleRate;

/// Calculates the width of a bin for a real valued Fourier transform
/// Hence `bin_count` is half of what it is normally (complex valued Fourier transform)
pub fn bin_width(sample_rate: SampleRate, bin_count: usize) -> f64 {
	sample_rate as f64 / (bin_count * 2) as f64
}

/// Calculates the central frequency of a bin
pub fn bin_center_frequency(sample_rate: SampleRate, bin_count: usize, bin_index: usize) -> Hertz {
	assert!(bin_index < bin_count);
	let width = bin_width(sample_rate, bin_count);
	width * bin_index as f64
}

pub fn bin_frequency_range(sample_rate: SampleRate, bin_count: usize, bin_index: usize) -> (Hertz, Hertz) {
	let half_width = bin_width(sample_rate, bin_count) / 2.0;
	let center_frequency = bin_center_frequency(sample_rate, bin_count, bin_index);
	(center_frequency - half_width, center_frequency + half_width)
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_bin_width() {
		assert_eq!(bin_width(40000, 32), 625.0);
	}

	#[test]
	fn test_bin_center_frequency() {
		assert_eq!(bin_center_frequency(40000, 32, 31), 19375.0);
	}

	#[test]
	fn test_bin_frequency_range() {
		assert_eq!(bin_frequency_range(40000, 32, 0), (-312.5, 312.5));
		assert_eq!(bin_frequency_range(40000, 32, 2), (937.5, 1562.5));
	}
}
