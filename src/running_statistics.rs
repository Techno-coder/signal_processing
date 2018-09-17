use std::iter::Iterator;
use super::Sample;

pub struct RunningVariance<I> where I: IntoIterator<Item=Sample> {
	samples_processed: u32,
	sample_sum: f64,
	square_sum: f64,
	samples: I::IntoIter,
}

impl<I> RunningVariance<I> where I: IntoIterator<Item=Sample> {
	pub fn new(samples: I) -> Option<RunningVariance<I>> {
		let mut samples = samples.into_iter();
		let sample = samples.next()?;
		Some(RunningVariance {
			samples_processed: 1,
			sample_sum: sample,
			square_sum: sample * sample,
			samples,
		})
	}
}

impl<I> Iterator for RunningVariance<I> where I: IntoIterator<Item=Sample> {
	type Item = f64;

	fn next(&mut self) -> Option<<Self as Iterator>::Item> {
		let sample = self.samples.next()?;
		self.samples_processed += 1;
		self.sample_sum += sample;
		self.square_sum += sample * sample;

		let sample_sum_mean = (self.sample_sum * self.sample_sum) / self.samples_processed as f64;
		Some((self.square_sum - sample_sum_mean) / (self.samples_processed - 1) as f64)
	}
}

pub struct RunningStandardDeviation<I> where I: IntoIterator<Item=Sample> {
	running_variance: RunningVariance<I>,
}

impl<I> RunningStandardDeviation<I> where I: IntoIterator<Item=Sample> {
	pub fn new(samples: I) -> Option<RunningStandardDeviation<I>> {
		Some(RunningStandardDeviation {
			running_variance: RunningVariance::new(samples)?,
		})
	}
}

impl<I> Iterator for RunningStandardDeviation<I> where I: IntoIterator<Item=Sample> {
	type Item = f64;

	fn next(&mut self) -> Option<<Self as Iterator>::Item> {
		Some(self.running_variance.next()?.sqrt())
	}
}

#[cfg(test)]
mod tests {
	use crate::statistics;
	use super::*;

	#[test]
	fn test_running_variance() {
		assert!(RunningVariance::new(vec![]).is_none());
		let variance = RunningVariance::new(vec![0.0, 1.0, 2.0, 3.0])
			.unwrap().last().unwrap();
		assert_eq!(variance, statistics::variance(&[0.0, 1.0, 2.0, 3.0]));
	}

	#[test]
	fn test_running_standard_deviation() {
		assert!(RunningStandardDeviation::new(vec![]).is_none());
		let standard_deviation = RunningStandardDeviation::new(vec![0.0, 1.0, 2.0, 3.0])
			.unwrap().last().unwrap();
		assert_eq!(standard_deviation, statistics::standard_deviation(&[0.0, 1.0, 2.0, 3.0]));
	}
}