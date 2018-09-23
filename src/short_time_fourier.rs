use crate::bin::Bin;
use crate::fourier_transform::FourierTransform;
use crate::rectangular::Rectangular;
use crate::window::Window;
use rayon::prelude::*;
use super::Sample;

pub fn analysis<T>(signal: &[Sample], overlap: usize, window: &Window)
                   -> Vec<Vec<Bin<Rectangular>>> where T: FourierTransform {
	assert!(overlap < window.width());
	let frame_spacing = window.width() - overlap;
	let final_frame_start = signal.len() - window.width();
	let loop_frame_count = final_frame_start / frame_spacing;

	let mut matrix: Vec<_> = (0..loop_frame_count).into_par_iter().map(|frame_index| {
		let frame_start = frame_index * frame_spacing;
		let frame = &signal[frame_start..frame_start + window.width()];
		T::analysis(&window.apply(frame))
	}).collect();

	let final_frame = window.apply(&signal[final_frame_start..]);
	matrix.push(T::analysis_extend(&final_frame, window.width()));
	matrix
}

pub fn synthesis<T>(matrix: &Vec<Vec<Bin<Rectangular>>>, overlap: usize, window: &Window)
                    -> Vec<Sample> where T: FourierTransform {
	assert!(overlap < window.width());
	let frames: Vec<_> = matrix.par_iter().map(|frame| {
		T::synthesis(&frame, window.width())
	}).collect();
	let frame_spacing = window.width() - overlap;
	let signal_length = (frames.len() * frame_spacing) + overlap;
	let mut signal = vec![Sample::default(); signal_length];

	let mut frame_start = 0;
	for frame in frames.iter() {
		for (index, sample) in frame.iter().enumerate() {
			let signal_index = frame_start + index;
			signal[signal_index] += window.apply_single(sample, index);
		}
		frame_start += frame_spacing;
	}
	signal
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
		let matrix = analysis::<CorrelationFourier>(&signal, 2, &window);
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
		let matrix = analysis::<CorrelationFourier>(&signal, 50, &analysis_window);

		let synthesis_window = Window::generate::<window::Dirichlet>(100);
		let mut signal = synthesis::<CorrelationFourier>(&matrix, 50, &synthesis_window);
		signal.iter_mut().for_each(|sample| *sample = sample.round());
		assert_eq!(utility::find_peak(&signal), Some(&151.0));
	}

	#[bench]
	#[cfg(feature = "fast_fourier")]
	fn bench_analysis_and_synthesis(bench: &mut Bencher) {
		use crate::fast_fourier::FastFourier;
		let signal: Vec<_> = (0..8192).map(|x| x as f64).collect();
		let window = Window::generate::<window::Sine>(256);
		let overlap = 128;
		bench.iter(|| synthesis::<FastFourier>(&analysis::<FastFourier>(&signal, overlap, &window), overlap, &window));
	}
}