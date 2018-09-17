use super::Sample;

/// Calculates a convolution of a signal
pub fn convolve_signal(signal: &[Sample], impulse_response: &[Sample]) -> Vec<Sample> {
	let mut convolution = vec![0.0; signal.len() + impulse_response.len() - 1];
	for (sample_index, sample) in signal.iter().enumerate() {
		for (impulse_index, impulse) in impulse_response.iter().enumerate() {
			convolution[sample_index + impulse_index] += sample * impulse;
		}
	}
	convolution
}

/// Calculates a single sample of the output convolution
pub fn convolve_single(signal: &[Sample], impulse_response: &[Sample], index: usize) -> Sample {
	assert!(index < signal.len() + impulse_response.len() - 1);
	let mut output_sample = 0.0;
	for (impulse_index, impulse) in impulse_response.iter().enumerate() {
		if impulse_index <= index && (index - impulse_index) < signal.len() {
			output_sample += signal[index - impulse_index] * impulse;
		}
	}
	output_sample
}

/// Convolution representing the slope between samples
pub fn first_difference(signal: &[Sample]) -> Vec<Sample> {
	convolve_signal(signal, &[1.0, -1.0])
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_convolve_signal() {
		let signal = [0.0, 1.0, 2.0, 3.0, 2.0, 0.0];
		let impulse_response = [1.0, 2.0];
		let convolution = convolve_signal(&signal, &impulse_response);
		assert_eq!(convolution, vec![0.0, 1.0, 2.0 + 2.0, 4.0 + 3.0, 6.0 + 2.0, 4.0, 0.0])
	}

	#[test]
	fn test_convolve_single() {
		let signal = [0.0, 1.0, 2.0, 3.0, 2.0, 0.0];
		let impulse_response = [1.0, 2.0];
		assert_eq!(convolve_single(&signal, &impulse_response, 4), 8.0);
		assert_eq!(convolve_single(&signal, &impulse_response, 1), 1.0);
		assert_eq!(convolve_single(&signal, &impulse_response, 6), 0.0);
	}

	#[test]
	fn test_first_difference() {
		let signal = [0.0, 1.0, 3.0, 8.0, 13.0, -1.0];
		assert_eq!(first_difference(&signal), vec![0.0, 1.0, 2.0, 5.0, 5.0, -14.0, 1.0])
	}
}
