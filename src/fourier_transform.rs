use crate::frequency::Frequency;
use crate::rectangular::Rectangular;
use std::f64;
use std::f64::consts;
use super::Sample;

pub trait FourierTransform {
	fn analysis_extend(signal: &[Sample], signal_length: usize) -> Vec<Frequency<Rectangular>>;
	fn synthesis(frequencies: &[Frequency<Rectangular>], signal_length: usize) -> Vec<Sample>;

	fn analysis(signal: &[Sample]) -> Vec<Frequency<Rectangular>> {
		Self::analysis_extend(signal, signal.len())
	}
}

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

pub struct CorrelationFourier();

impl FourierTransform for CorrelationFourier {
	fn analysis_extend(signal: &[Sample], signal_length: usize) -> Vec<Frequency<Rectangular>> {
		let upper_bound = (signal_length + 1) / 2;
		(0..upper_bound).map(|k| {
			let cosine = (0..signal_length).map(|i| signal.get(i).unwrap_or(&0.0) *
				cosine_basis_single(k as f64, signal_length, i)).sum();
			let sine = (0..signal_length).map(|i| -signal.get(i).unwrap_or(&0.0) *
				sine_basis_single(k as f64, signal_length, i)).sum();
			Rectangular { cosine, sine }.into()
		}).collect()
	}

	fn synthesis(frequencies: &[Frequency<Rectangular>], signal_length: usize) -> Vec<Sample> {
		let upper_bound = (signal_length + 1) / 2;
		assert!(frequencies.len() >= upper_bound);
		(0..signal_length).map(|i| {
			(0..upper_bound).map(|k| {
				let cosine = normalize_cosine_amplitude(k as f64, frequencies[k].cosine, signal_length) *
					cosine_basis_single(k as f64, signal_length, i);
				let sine = normalize_sine_amplitude(frequencies[k].sine, signal_length) *
					sine_basis_single(k as f64, signal_length, i);
				cosine + sine
			}).sum()
		}).collect()
	}
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
		assert!(sine_basis(0.0, 10).iter().all(|x| x == &0.0));
		assert!(sine_basis(16.0, 32).into_iter().map(math::approximate).all(|x| x == 0.0));
	}

	#[test]
	fn test_synthesis() {
		let frequencies = [
			Frequency(Rectangular { cosine: 15.0, sine: 0.0 }),
			Frequency(Rectangular { cosine: -2.5, sine: 3.4409548 }),
			Frequency(Rectangular { cosine: -2.5, sine: 0.81229924 }),
		];
		let synthesis: Vec<_> = CorrelationFourier::synthesis(&frequencies, 5)
			.into_iter().map(math::approximate).collect();
		assert_eq!(synthesis, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
	}

	#[test]
	fn test_analysis() {
		let signal = [1.0, 2.0, 3.0, 4.0, 5.0];
		let frequencies = CorrelationFourier::analysis(&signal);
		let synthesis: Vec<_> = CorrelationFourier::synthesis(&frequencies, 5)
			.into_iter().map(math::approximate).collect();
		assert_eq!(synthesis, signal);
	}
}