use std::collections::BTreeMap;
use super::IntegralSample;

pub struct Histogram {
	samples: BTreeMap<IntegralSample, usize>,
	sample_count: usize,
}

impl Histogram {
	pub fn new() -> Histogram {
		Histogram {
			samples: BTreeMap::new(),
			sample_count: 0,
		}
	}

	pub fn append(&mut self, sample: IntegralSample) {
		use std::ops::AddAssign;
		self.samples.entry(sample).or_insert(0).add_assign(1);
		self.sample_count += 1;
	}

	pub fn mean(&self) -> f64 {
		let sum: i64 = self.samples.iter()
		                   .map(|(x, count)| *x * *count as i64)
		                   .sum();
		sum as f64 / self.sample_count as f64
	}

	pub fn variance(&self) -> f64 {
		use crate::statistics;
		let mean = self.mean();
		let square_sum: f64 = self.samples.iter()
		                          .map(|(x, count)| (statistics::deviation(*x as f64, mean), count))
		                          .map(|(x, count)| (x * x) * *count as f64)
		                          .sum();
		square_sum / (self.sample_count - 1) as f64
	}

	pub fn standard_deviation(&self) -> f64 {
		self.variance().sqrt()
	}
}

#[cfg(test)]
mod tests {
	use crate::statistics;
	use super::*;

	#[test]
	fn test_histogram() {
		let mut histogram = Histogram::new();
		histogram.append(0);
		assert_eq!(histogram.mean(), 0.0);
		assert!(histogram.variance().is_nan());
		assert!(histogram.standard_deviation().is_nan());

		histogram.append(1);
		assert_eq!(histogram.mean(), 0.5);
		assert_eq!(histogram.variance(), statistics::variance(&[0.0, 1.0]));
		assert_eq!(histogram.standard_deviation(), statistics::standard_deviation(&[0.0, 1.0]));

		histogram.append(1);
		assert_eq!(histogram.mean(), 2.0 / 3.0);
		assert_eq!(histogram.variance(), statistics::variance(&[0.0, 1.0, 1.0]));
		assert_eq!(histogram.standard_deviation(), statistics::standard_deviation(&[0.0, 1.0, 1.0]));
	}
}