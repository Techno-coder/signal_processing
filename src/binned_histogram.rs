use crate::histogram::Histogram;
use std::ops::Deref;
use super::Sample;

pub struct BinnedHistogram {
	histogram: Histogram,
	bin_interval: u64,
	bin_interval_center: u64,
}

impl BinnedHistogram {
	/// Bin width must be a whole number
	pub fn new(bin_interval: u64) -> BinnedHistogram {
		BinnedHistogram {
			histogram: Histogram::new(),
			bin_interval,
			bin_interval_center: bin_interval / 2,
		}
	}

	pub fn append(&mut self, sample: Sample) {
		let bin_index = (sample / self.bin_interval as f64).round() as i64;
		let bin_sample = bin_index * self.bin_interval as i64;
		let normalized_sample = bin_sample + self.bin_interval_center as i64;
		self.histogram.append(normalized_sample);
	}
}

impl Deref for BinnedHistogram {
	type Target = Histogram;

	fn deref(&self) -> &<Self as Deref>::Target {
		&self.histogram
	}
}

#[cfg(test)]
mod tests {
	use crate::statistics;
	use super::*;

	#[test]
	fn test_binned_histogram() {
		let mut histogram = BinnedHistogram::new(1000);
		histogram.append(0.0);
		assert_eq!(histogram.mean(), 500.0);
		assert!(histogram.variance().is_nan());

		histogram.append(1000.0);
		assert_eq!(histogram.mean(), statistics::mean(&[500.0, 1500.0]));
		assert_eq!(histogram.variance(), statistics::variance(&[500.0, 1500.0]));

		histogram.append(-1000.0);
		assert_eq!(histogram.mean(), statistics::mean(&[-500.0, 500.0, 1500.0]));
		assert_eq!(histogram.variance(), statistics::variance(&[-500.0, 500.0, 1500.0]));
	}
}