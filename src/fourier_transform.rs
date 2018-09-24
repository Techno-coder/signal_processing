use crate::bin::Bin;
use crate::rectangular::Rectangular;
use std::f64;
use std::f64::consts;
use super::Sample;

pub trait FourierTransform: Send + Sync {
	fn analysis_extend(signal: &[Sample], signal_length: usize) -> Vec<Bin<Rectangular>>;
	fn synthesis(bins: &[Bin<Rectangular>], signal_length: usize) -> Vec<Sample>;

	fn analysis(signal: &[Sample]) -> Vec<Bin<Rectangular>> {
		Self::analysis_extend(signal, signal.len())
	}
}

pub fn cosine_basis_single(bin_index: usize, signal_length: usize, index: usize) -> Sample {
	(2.0 * consts::PI * bin_index as f64 * (index as f64 / signal_length as f64)).cos()
}

pub fn cosine_basis(bin_index: usize, signal_length: usize) -> Vec<Sample> {
	(0..signal_length).map(|i| cosine_basis_single(bin_index, signal_length, i)).collect()
}

pub fn sine_basis_single(bin_index: usize, signal_length: usize, index: usize) -> Sample {
	(2.0 * consts::PI * bin_index as f64 * (index as f64 / signal_length as f64)).sin()
}

pub fn sine_basis(bin_index: usize, signal_length: usize) -> Vec<Sample> {
	(0..signal_length).map(|i| sine_basis_single(bin_index, signal_length, i)).collect()
}

pub fn normalize_cosine_amplitude(bin_index: usize, amplitude: Sample, signal_length: usize) -> Sample {
	if bin_index == 0 || bin_index == ((signal_length + 1) / 2) {
		amplitude / signal_length as f64
	} else {
		amplitude / (signal_length as f64 / 2.0)
	}
}

pub fn normalize_sine_amplitude(amplitude: Sample, signal_length: usize) -> Sample {
	-amplitude / (signal_length as f64 / 2.0)
}

/// Calculates the number of bins generated for a real valued Fourier transform
pub fn bin_count(signal_length: usize) -> usize {
	(signal_length / 2) + 1
}

pub struct CorrelationFourier();

impl FourierTransform for CorrelationFourier {
	fn analysis_extend(signal: &[Sample], signal_length: usize) -> Vec<Bin<Rectangular>> {
		let bin_count = bin_count(signal_length);
		(0..bin_count).map(|k| {
			let cosine = (0..signal_length).map(|i| signal.get(i).unwrap_or(&0.0) *
				cosine_basis_single(k, signal_length, i)).sum();
			let sine = (0..signal_length).map(|i| -signal.get(i).unwrap_or(&0.0) *
				sine_basis_single(k, signal_length, i)).sum();
			Rectangular { cosine, sine }.into()
		}).collect()
	}

	fn synthesis(bins: &[Bin<Rectangular>], signal_length: usize) -> Vec<Sample> {
		let bin_count = bin_count(signal_length);
		assert!(bins.len() >= bin_count);
		(0..signal_length).map(|i| {
			(0..bin_count).map(|k| {
				let cosine = normalize_cosine_amplitude(k, bins[k].cosine, signal_length) *
					cosine_basis_single(k, signal_length, i);
				let sine = normalize_sine_amplitude(bins[k].sine, signal_length) *
					sine_basis_single(k, signal_length, i);
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
		assert!(cosine_basis(0, 10).iter().all(|x| x == &1.0));
		assert!(cosine_basis(16, 32).iter().step_by(2).all(|x| x == &1.0));
		assert!(cosine_basis(16, 32).iter().skip(1).step_by(2).all(|x| x == &-1.0));
	}

	#[test]
	fn test_sine_basis() {
		assert!(sine_basis(0, 10).iter().all(|x| x == &0.0));
		assert!(sine_basis(16, 32).into_iter().map(math::approximate).all(|x| x == 0.0));
	}

	#[test]
	fn test_synthesis() {
		let bins = [
			Bin(Rectangular { cosine: 15.0, sine: 0.0 }),
			Bin(Rectangular { cosine: -2.5, sine: 3.4409548 }),
			Bin(Rectangular { cosine: -2.5, sine: 0.81229924 })
		];
		let synthesis: Vec<_> = CorrelationFourier::synthesis(&bins, 5)
			.into_iter().map(math::approximate).collect();
		assert_eq!(synthesis, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
	}

	#[test]
	fn test_synthesis_even_length() {
		let bins = [
			Bin(Rectangular { cosine: 10.0, sine: 0.0 }),
			Bin(Rectangular { cosine: -2.0, sine: 2.0 }),
			Bin(Rectangular { cosine: -2.0, sine: 0.0 }),
		];
		let synthesis: Vec<_> = CorrelationFourier::synthesis(&bins, 4);
		assert_eq!(synthesis, vec![1.0, 2.0, 3.0, 4.0]);
	}

	#[test]
	fn test_analysis() {
		let signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
		let bins = CorrelationFourier::analysis(&signal);
		let synthesis: Vec<_> = CorrelationFourier::synthesis(&bins, 6)
			.into_iter().map(math::approximate).collect();
		assert_eq!(synthesis, signal);
	}
}