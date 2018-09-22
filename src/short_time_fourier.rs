use crate::fourier_transform::FourierTransform;
use crate::frequency::Frequency;
use crate::rectangular::Rectangular;
use crate::window::Window;
use super::Sample;

pub fn analysis<T>(signal: &[Sample], overlap: usize, window: &Window)
                   -> Vec<Vec<Frequency<Rectangular>>> where T: FourierTransform {
	assert!(overlap < window.width());
	let chunk_spacing = window.width() - overlap;
	let mut matrix = Vec::new();

	let mut chunk_start = 0;
	while chunk_start + window.width() < signal.len() {
		let chunk = &signal[chunk_start..chunk_start + window.width()];
		matrix.push(T::analysis(&window.apply(chunk)));
		chunk_start += chunk_spacing;
	}

	let final_chunk = window.apply(&signal[chunk_start..]);
	matrix.push(T::analysis_extend(&final_chunk, window.width()));
	matrix
}

pub fn synthesis<T>(matrix: &Vec<Vec<Frequency<Rectangular>>>, overlap: usize, window: &Window)
                    -> Vec<Sample> where T: FourierTransform {
	assert!(overlap < window.width());
	let chunks: Vec<_> = matrix.iter().map(|chunk| {
		T::synthesis(&chunk, window.width())
	}).collect();
	let chunk_spacing = window.width() - overlap;
	let signal_length = (chunks.len() * chunk_spacing) + overlap;
	let mut signal = vec![Sample::default(); signal_length];

	let mut chunk_start = 0;
	for chunk in chunks.iter() {
		for (index, sample) in chunk.iter().enumerate() {
			let signal_index = chunk_start + index;
			signal[signal_index] += *sample;
		}
		chunk_start += chunk_spacing;
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

	#[test]
	fn test_analysis() {
		let signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0, 6.0];
		let window = Window::generate::<window::Hann>(4);
		let matrix = analysis::<CorrelationFourier>(&signal, 2, &window);
		let matrix: Vec<Vec<_>> = matrix.into_iter().map(|chunk| chunk.into_iter().map(|frequency| {
			Frequency(Rectangular {
				cosine: math::approximate(frequency.cosine),
				sine: math::approximate(frequency.sine),
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
		let window = Window::generate::<window::Hann>(100);
		let matrix = analysis::<CorrelationFourier>(&signal, 50, &window);
		let mut signal = synthesis::<CorrelationFourier>(&matrix, 50, &window);
		signal.iter_mut().for_each(|sample| *sample = sample.round());
		assert_eq!(utility::find_peak(&signal), Some(&151.0));
	}
}