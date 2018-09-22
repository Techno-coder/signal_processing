use std::f64::consts;
use super::Sample;

#[derive(Debug, Clone)]
pub struct Window {
	window: Vec<Sample>,
}

impl Window {
	pub fn generate<F>(length: usize) -> Window where F: WindowFunction {
		assert!(length > 0);
		Window {
			window: F::generate(length),
		}
	}

	pub fn apply(&self, signal: &[Sample]) -> Vec<Sample> {
		signal.iter().enumerate()
		      .map(|(index, sample)| self.apply_single(sample, index))
		      .collect()
	}

	pub fn apply_single(&self, sample: &Sample, index: usize) -> Sample {
		self.window.get(index).unwrap_or(&0.0) * sample
	}

	pub fn normalize_amplitude(&self, overlap_factor: f64) -> Window {
		assert!(overlap_factor > 0.0);
		Window {
			window: self.window.iter().map(|x| x / overlap_factor).collect()
		}
	}

	pub fn width(&self) -> usize {
		self.window.len()
	}
}

pub trait WindowFunction {
	fn generate(length: usize) -> Vec<Sample>;
}

pub struct Sine();

impl WindowFunction for Sine {
	fn generate(length: usize) -> Vec<Sample> {
		let denominator = (length - 1) as f64;
		(0..length)
			.map(|n| ((consts::PI * n as f64) / denominator).sin())
			.collect()
	}
}

pub struct Hann();

impl WindowFunction for Hann {
	fn generate(length: usize) -> Vec<Sample> {
		let mut window = Sine::generate(length);
		window.iter_mut().for_each(|x| *x = *x * *x);
		window
	}
}

pub struct Dirichlet();

impl WindowFunction for Dirichlet {
	fn generate(length: usize) -> Vec<Sample> {
		vec![1.0; length]
	}
}

#[cfg(test)]
mod tests {
	use crate::math;
	use super::*;

	#[test]
	fn test_dirichlet_window() {
		let signal = [1.0, 2.0, 3.0, 4.0, 5.0];
		let window = Window::generate::<Dirichlet>(3);
		assert_eq!(window.apply(&signal), &[1.0, 2.0, 3.0, 0.0, 0.0]);
	}

	#[test]
	fn test_hann_window() {
		let signal = [1.0; 5];
		let window = Window::generate::<Hann>(5);
		let output: Vec<_> = window.apply(&signal).into_iter().map(math::approximate).collect();
		assert_eq!(&output, &[0.0, 0.5, 1.0, 0.5, 0.0]);

		let signal = [3.0, 4.0, 5.0, 6.0];
		let window = Window::generate::<Hann>(4);
		let output: Vec<_> = window.apply(&signal).into_iter().map(math::approximate).collect();
		assert_eq!(&output, &[0.0, 3.0, 3.75, 0.0]);
	}

	#[test]
	fn test_hann_window_overlap() {
		let window = Window::generate::<Hann>(512);
		let mut signal = [0.0; 768];
		for (index, sample) in window.apply(&[1.0; 512]).into_iter().enumerate() {
			signal[index] += sample;
			signal[index + 256] += sample;
		}

		signal.iter_mut().for_each(|x| *x = x.round());
		assert_eq!(&signal[256..512], &[1.0 as f64; 256][..]);
	}
}