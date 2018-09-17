use std::f64;
use std::f64::consts;
use super::Sample;

pub fn cosine_basis_single(frequency: f64, signal_length: usize, index: usize) -> Sample {
	(2.0 * consts::PI * frequency * (index as f64 / signal_length as f64)).cos()
}

pub fn cosine_basis(frequency: f64, signal_length: usize) -> Vec<Sample> {
	(0..signal_length).map(|i| cosine_basis_single(frequency, signal_length, i)).collect()
}

pub fn sine_basis_single(frequency: f64, signal_length: usize, index: usize) -> Sample {
	(2.0 * consts::PI * frequency * (index as f64 / signal_length as f64)).sin()
}

pub fn sine_basis(frequency: f64, signal_length: usize) -> Vec<Sample> {
	(0..signal_length).map(|i| sine_basis_single(frequency, signal_length, i)).collect()
}

pub fn normalize_cosine_amplitude(frequency: f64, amplitude: Sample, signal_length: usize) -> Sample {
	if frequency == 0.0 || frequency == ((signal_length + 1) / 2) as f64 {
		amplitude / signal_length as f64
	} else {
		amplitude / (signal_length as f64 / 2.0)
	}
}

pub fn normalize_sine_amplitude(amplitude: Sample, signal_length: usize) -> Sample {
	-amplitude / (signal_length as f64 / 2.0)
}

pub fn synthesis(cosine_amplitudes: &[Sample], sine_amplitudes: &[Sample], signal_length: usize) -> Vec<Sample> {
	let upper_bound = (signal_length + 1) / 2;
	assert!(cosine_amplitudes.len() >= upper_bound);
	assert_eq!(cosine_amplitudes.len(), sine_amplitudes.len());
	(0..signal_length)
		.map(|i| {
			let cosines: f64 = (0..upper_bound)
				.map(|k| normalize_cosine_amplitude(k as f64, cosine_amplitudes[k], signal_length) *
					cosine_basis_single(k as f64, signal_length, i)).sum();
			let sines: f64 = (0..upper_bound)
				.map(|k| normalize_sine_amplitude(sine_amplitudes[k], signal_length) *
					sine_basis_single(k as f64, signal_length, i)).sum();
			cosines + sines
		})
		.collect()
}

pub fn fourier_transform(signal: &[Sample]) -> (Vec<Sample>, Vec<Sample>) {
	let upper_bound = (signal.len() + 1) / 2;
	((0..upper_bound).map(|k| (0..signal.len()).map(|i| signal[i] *
		cosine_basis_single(k as f64, signal.len(), i)).sum()).collect(),
	 (0..upper_bound).map(|k| (0..signal.len()).map(|i| -signal[i] *
		 sine_basis_single(k as f64, signal.len(), i)).sum()).collect())
}

#[cfg(test)]
mod tests {
	use crate::math;
	use super::*;

	#[test]
	fn test_cosine_basis() {
		assert!(cosine_basis(0.0, 10).iter().all(|x| x == &1.0));
		assert!(cosine_basis(16.0, 32).iter().step_by(2).all(|x| x == &1.0));
		assert!(cosine_basis(16.0, 32).iter().skip(1).step_by(2).all(|x| x == &-1.0));
	}

	#[test]
	fn test_sine_basis() {
		const ZERO_EPSILON: f64 = 0.0000000000001;
		assert!(sine_basis(0.0, 10).iter().all(|x| x == &0.0));
		assert!(sine_basis(16.0, 32).iter().all(|x| &-ZERO_EPSILON < x && x < &ZERO_EPSILON));
	}

	#[test]
	fn test_synthesis() {
		let cosine_amplitudes = [15.0, -2.5, -2.5, -2.5, -2.5];
		let sine_amplitudes = [0.0, 3.4409548, 0.81229924, -0.81229924, -3.4409548];
		let synthesis: Vec<_> = synthesis(&cosine_amplitudes, &sine_amplitudes, 5)
			.into_iter().map(math::approximate).collect();
		assert_eq!(synthesis, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
	}

	#[test]
	fn test_fourier_transform() {
		let signal = [1.0, 2.0, 3.0, 4.0, 5.0];
		let (cosine_amplitudes, sine_amplitudes) = fourier_transform(&signal);
		let synthesis: Vec<_> = synthesis(&cosine_amplitudes, &sine_amplitudes, 5)
			.into_iter().map(math::approximate).collect();
		assert_eq!(synthesis, signal);
	}
}