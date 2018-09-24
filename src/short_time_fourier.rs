use crate::bin::Bin;
use crate::fourier_transform::FourierTransform;
use crate::rectangular::Rectangular;
use crate::window::Window;
use rayon::prelude::*;
use std::collections::VecDeque;
use std::iter;
use std::iter::FromIterator;
use std::marker::PhantomData;
use super::Sample;

pub struct ShortTimeAnalyser<'a, T> {
	signal: &'a [Sample],
	window: &'a Window,

	frame_spacing: usize,
	loop_frame_count: usize,
	_transform: PhantomData<T>,
}

impl<'a, T> ShortTimeAnalyser<'a, T> where T: FourierTransform {
	pub fn new(signal: &'a [Sample], overlap: usize, window: &'a Window) -> Self {
		assert!(overlap < window.width());
		let frame_spacing = window.width() - overlap;
		ShortTimeAnalyser {
			signal,
			window,
			frame_spacing,
			loop_frame_count: (signal.len() - window.width()) / frame_spacing,
			_transform: Default::default(),
		}
	}

	pub fn calculate_frame(&self, frame_index: usize) -> Vec<Bin<Rectangular>> {
		assert!(frame_index < self.total_frames());
		let frame_start = frame_index * self.frame_spacing;
		if frame_index <= self.loop_frame_count {
			debug_assert!(frame_start + self.window.width() <= self.signal.len());
			let frame = &self.signal[frame_start..frame_start + self.window.width()];
			T::analysis(&self.window.apply(frame))
		} else {
			debug_assert!(frame_start + self.window.width() > self.signal.len());
			let final_frame = self.window.apply(&self.signal[frame_start..]);
			T::analysis_extend(&final_frame, self.window.width())
		}
	}

	pub fn total_frames(&self) -> usize {
		self.loop_frame_count + 1
	}

	pub fn calculate_all(&self) -> Vec<Vec<Bin<Rectangular>>> {
		(0..self.total_frames()).into_par_iter().map(|frame_index| self.calculate_frame(frame_index)).collect()
	}
}

pub struct ShortTimeSynthesiser<'a, T> {
	samples: VecDeque<Sample>,
	window: &'a Window,
	overlap: usize,
	frame_spacing: usize,
	frame_complete_length: usize,
	overlapping_frames_count: usize,
	_transform: PhantomData<T>,
}

impl<'a, T> ShortTimeSynthesiser<'a, T> where T: FourierTransform {
	pub fn new(overlap: usize, window: &'a Window) -> Self {
		assert!(overlap < window.width());
		let frame_spacing = window.width() - overlap;
		ShortTimeSynthesiser {
			samples: VecDeque::from_iter(iter::repeat(0.0).take(overlap)),
			window,
			overlap,
			frame_spacing,
			frame_complete_length: window.width() + overlap,
			overlapping_frames_count: (window.width() as f64 / frame_spacing as f64).ceil() as usize,
			_transform: Default::default(),
		}
	}

	pub fn push_frames(&mut self, frames: &Vec<Vec<Bin<Rectangular>>>) {
		let complete_end = self.samples.len() - self.overlap;
		(0..(self.frame_spacing * frames.len())).for_each(|_| self.samples.push_back(0.0));
		for level in 0..self.overlapping_frames_count {
			let parallel_frame_count = (frames.len() + self.overlapping_frames_count)
				.saturating_sub(level + 1) / self.overlapping_frames_count;
			(0..parallel_frame_count).into_par_iter().for_each(|index| {
				let frame_index = level + (index * self.overlapping_frames_count);
				let frame = T::synthesis(&frames[frame_index], self.window.width());
				let frame_start = complete_end + (frame_index * self.frame_spacing);
				for (index, sample) in frame.iter().enumerate() {
					let sample_index = frame_start + index;
					let pointer = &self.samples[sample_index] as *const _ as *mut Sample;
					unsafe { *pointer += self.window.apply_single(sample, index); }
				}
			});
		}
	}

	pub fn flush_ready(&mut self) -> Vec<Sample> {
		let mut samples = Vec::new();
		while self.samples.len() > self.frame_complete_length {
			samples.push(self.samples.pop_front().unwrap());
		}
		samples
	}

	pub fn flush_all(self) -> Vec<Sample> {
		self.samples.into_iter().collect()
	}
}

#[cfg(test)]
mod tests {
	use crate::fourier_transform::CorrelationFourier;
	use crate::math;
	use crate::utility;
	use crate::window;
	use super::*;
	use test::Bencher;

	#[test]
	fn test_analysis() {
		let signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0, 6.0];
		let window = Window::generate::<window::Hann>(4);
		let matrix = ShortTimeAnalyser::<CorrelationFourier>::new(&signal, 2, &window).calculate_all();
		let matrix: Vec<Vec<_>> = matrix.into_iter().map(|frame| frame.into_iter().map(|bin| {
			Bin(Rectangular {
				cosine: math::approximate(bin.cosine),
				sine: math::approximate(bin.sine),
			})
		}).collect()).collect();
		assert_eq!(&matrix, &[
			[
				Rectangular { cosine: 3.75, sine: 0.0 }.into(),
				Rectangular { cosine: -2.25, sine: -1.5 }.into(),
				Rectangular { cosine: 0.75, sine: 0.0 }.into(),
			],
			[
				Rectangular { cosine: 6.75, sine: 0.0 }.into(),
				Rectangular { cosine: -3.75, sine: -3.0 }.into(),
				Rectangular { cosine: 0.75, sine: 0.0 }.into(),
			],
			[
				Rectangular { cosine: 9.0, sine: 0.0 }.into(),
				Rectangular { cosine: -4.5, sine: -4.5 }.into(),
				Rectangular { cosine: 0.0, sine: 0.0 }.into(),
			]
		]);
	}

	#[test]
	fn test_synthesis() {
		let signal: Vec<_> = (0..200).map(|x| x as f64).collect();
		let analysis_window = Window::generate::<window::Hann>(100);
		let matrix = ShortTimeAnalyser::<CorrelationFourier>::new(&signal, 50, &analysis_window).calculate_all();

		let synthesis_window = Window::generate::<window::Dirichlet>(100);
		let mut synthesiser = ShortTimeSynthesiser::<CorrelationFourier>::new(50, &synthesis_window);
		synthesiser.push_frames(&matrix);
		let signal: Vec<_> = synthesiser.flush_all().iter().map(|sample| sample.round()).collect();
		assert_eq!(utility::find_peak(&signal), Some(&151.0));
	}

	#[test]
	fn test_analysis_total_frames() {
		let signal: Vec<_> = (0..11025).map(|x| x as f64).collect();
		let analysis_window = Window::generate::<window::Hann>(128);
		let analyser = ShortTimeAnalyser::<CorrelationFourier>::new(&signal, 120, &analysis_window);
		assert_eq!(analyser.total_frames(), analyser.calculate_all().len());
	}

	#[bench]
	#[cfg(feature = "fast_fourier")]
	fn bench_analysis_and_synthesis(bench: &mut Bencher) {
		use crate::fast_fourier::FastFourier;
		let signal: Vec<_> = (0..8192).map(|x| x as f64).collect();
		let window = Window::generate::<window::Sine>(256);
		let overlap = 128;
		bench.iter(|| {
			let analysis = ShortTimeAnalyser::<FastFourier>::new(&signal, overlap, &window).calculate_all();
			let mut synthesiser = ShortTimeSynthesiser::<FastFourier>::new(overlap, &window);
			synthesiser.push_frames(&analysis);
			synthesiser.flush_all()
		});
	}
}