use crate::fourier_transform::FourierTransform;
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

/// Calculates a convolution with discrete fourier transforms
pub fn convolve_fourier<T>(signal: &[Sample], impulse_response: &[Sample]) -> Vec<Sample> where T: FourierTransform {
	let convolution_length = signal.len() + impulse_response.len() - 1;
	let signal_frequencies = T::analysis_extend(&signal, convolution_length);
	let kernel_frequencies = T::analysis_extend(&impulse_response, convolution_length);

	let frequency_length = (convolution_length + 1) / 2;
	let output_frequencies: Vec<_> = (0..frequency_length)
		.map(|index| signal_frequencies[index] * kernel_frequencies[index]).collect();
	T::synthesis(&output_frequencies, convolution_length)
}

#[cfg(test)]
mod tests {
	use crate::fourier_transform::CorrelationFourier;
	use crate::math;
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
	fn test_convolve_fourier() {
		let signal = [0.0, 1.0, 2.0, 3.0, 2.0, 0.0];
		let impulse_response = [1.0, 2.0];
		let convolution: Vec<f64> = convolve_fourier::<CorrelationFourier>(&signal, &impulse_response)
			.into_iter().map(math::approximate).collect();
		assert_eq!(convolution, vec![0.0, 1.0, 2.0 + 2.0, 4.0 + 3.0, 6.0 + 2.0, 4.0, 0.0])
	}
}
