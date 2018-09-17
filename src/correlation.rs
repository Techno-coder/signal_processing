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

#[cfg(test)]
mod tests {
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
}