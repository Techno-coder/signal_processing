use super::Sample;

pub fn mean(samples: &[Sample]) -> f64 {
	let sum: f64 = samples.iter().sum();
	sum / samples.len() as f64
}

pub fn deviation(sample: Sample, mean: f64) -> f64 {
	sample - mean
}

pub fn variance(samples: &[Sample]) -> f64 {
	let mean = mean(samples);
	let square_sum: f64 = samples.iter()
	                              .map(|x| deviation(*x, mean))
	                              .map(|x| x * x)
	                              .sum();
	square_sum / (samples.len() - 1) as f64
}

pub fn standard_deviation(samples: &[Sample]) -> f64 {
	variance(samples).sqrt()
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_mean() {
		assert_eq!(mean(&[1.0, 2.0, 3.0, 4.0]), 2.5);
	}

	#[test]
	fn test_deviation() {
		assert_eq!(deviation(5.0, -2.5), 7.5);
	}

	#[test]
	fn test_variance() {
		assert!(variance(&[0.0]).is_nan());
		assert_eq!(variance(&[-1.0, 0.0, 1.0, 2.0]), (6.0 - 1.0) / 3.0);
	}

	#[test]
	fn test_standard_deviation() {
		assert!(standard_deviation(&[0.0]).is_nan());
		assert_eq!(standard_deviation(&[-1.0, 0.0, 1.0, 2.0]), (5.0 as f64 / 3.0).sqrt());
	}
}