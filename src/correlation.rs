use super::fourier_transform::FourierTransform;
use super::fourier_transform;
use super::polar::Polar;
use super::rectangular::Rectangular;
use super::Sample;

pub fn correlate_signal(signal: &[Sample], target: &[Sample]) -> Vec<Sample> {
	(0..signal.len()).map(|index| correlate_single(signal, target, index)).collect()
}

pub fn correlate_single(signal: &[Sample], target: &[Sample], index: usize) -> Sample {
	assert!(index < signal.len());
	let mut output_sample = 0.0;
	for (target_index, target) in target.iter().enumerate() {
		if index + target_index < signal.len() {
			output_sample += signal[index + target_index] * target;
		}
	}
	output_sample
}

pub fn correlation<T>(signal: &[Sample], target: &[Sample]) -> f64 where T: FourierTransform {
	correlate_fourier::<T>(signal, target).iter().sum()
}

pub fn correlate_fourier<T>(signal: &[Sample], target: &[Sample]) -> Vec<Sample> where T: FourierTransform {
	let signal_bins = T::analysis(&signal);
	let target_bins = T::analysis_extend(&target, signal.len());

	let mut output_bins = Vec::new();
	for index in 0..(fourier_transform::bin_count(signal.len())) {
		let target_bin: Polar = target_bins[index].take().into();
		let target_bin: Rectangular = target_bin.complex_conjugate().into();
		output_bins.push(signal_bins[index] * target_bin.into());
	}
	T::synthesis(&output_bins, signal.len())
}

#[cfg(test)]
mod tests {
	use crate::fourier_transform::CorrelationFourier;
	use super::*;

	#[test]
	fn test_correlate_single() {
		let signal = [0.0, 1.0, 2.0];
		let target = [1.0];
		assert_eq!(correlate_single(&signal, &target, 0), 0.0);
		assert_eq!(correlate_single(&signal, &target, 1), 1.0);
		assert_eq!(correlate_single(&signal, &target, 2), 2.0);
	}

	#[test]
	fn test_correlate_signal() {
		let signal = [0.0, 1.0, 2.0, 5.0, 2.0, 1.0];
		let target = [1.0, 2.0];
		assert_eq!(correlate_signal(&signal, &target), vec![2.0, 5.0, 12.0, 9.0, 4.0, 1.0]);
	}

	#[test]
	fn test_correlation() {
		let signal = [0.0, 1.0, 2.0, 5.0, 2.0, 1.0];
		let target = [1.0, 2.0];
		assert_eq!(correlation::<CorrelationFourier>(&signal, &target), vec![2.0, 5.0, 12.0, 9.0, 4.0, 1.0].into_iter().sum());
	}
}