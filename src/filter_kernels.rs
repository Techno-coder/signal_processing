use super::Sample;

/// Impulse response representing the slope between samples
pub fn first_difference() -> Vec<Sample> {
	vec![1.0, -1.0]
}

pub fn low_pass_square_pulse(amplitude: f64, width: usize) -> Vec<Sample> {
	vec![amplitude; width]
}

pub fn high_pass_square_pulse(delta_amplitude: f64, width: usize) -> Vec<Sample> {
	assert!(width > 1);
	let dampen_amplitude = delta_amplitude / (width - 1) as f64;
	let mut filter = vec![dampen_amplitude; width];
	filter[width / 2] = delta_amplitude;
	filter
}

#[cfg(test)]
mod tests {
	use crate::convolution;
	use super::*;

	#[test]
	fn test_first_difference() {
		let signal = [0.0, 1.0, 3.0, 8.0, 13.0, -1.0];
		assert_eq!(convolution::convolve_signal(&signal, &first_difference()),
		           vec![0.0, 1.0, 2.0, 5.0, 5.0, -14.0, 1.0])
	}
}
