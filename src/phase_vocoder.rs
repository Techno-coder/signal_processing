use crate::bin::Bin;
///! A phase vocoder allows pitching shifting or time scaling without
///! changing the other domain.
///! Adapted from http://blogs.zynaptiq.com/bernsee/pitch-shifting-using-the-ft/

use crate::bin_frequency;
use crate::fourier_transform;
use crate::fourier_transform::FourierTransform;
use crate::polar::Polar;
use crate::rectangular::Rectangular;
use crate::short_time_fourier;
use crate::window::Window;
use rayon::prelude::*;
use std::f64::consts;
use super::Hertz;
use super::Sample;
use super::SampleRate;

#[derive(Debug, Default, Copy, Clone)]
pub struct PhaseVocoderBin {
	pub magnitude: Sample,
	pub frequency: Hertz,
}

pub trait PhaseVocoderProcessor: Send + Sync {
	fn process(&self, bins: Vec<PhaseVocoderBin>, output_bins: &mut Vec<PhaseVocoderBin>);
}

#[derive(Debug, Clone)]
pub struct PitchShift {
	pub pitch_shift_ratio: f64,
}

impl PhaseVocoderProcessor for PitchShift {
	fn process(&self, bins: Vec<PhaseVocoderBin>, output_bins: &mut Vec<PhaseVocoderBin>) {
		for bin_index in 0..bins.len() {
			let pitch_shifted_index = (bin_index as f64 * self.pitch_shift_ratio) as i64;
			if 0 <= pitch_shifted_index && pitch_shifted_index < bins.len() as i64 {
				let pitch_shifted_index = pitch_shifted_index as usize;
				output_bins[pitch_shifted_index].magnitude += bins[bin_index].magnitude;
				output_bins[pitch_shifted_index].frequency = bins[bin_index].frequency * self.pitch_shift_ratio;
			}
		}
	}
}

pub struct IdentityProcessor();

impl PhaseVocoderProcessor for IdentityProcessor {
	fn process(&self, bins: Vec<PhaseVocoderBin>, output_bins: &mut Vec<PhaseVocoderBin>) {
		use std::mem;
		mem::replace(output_bins, bins);
	}
}

pub fn process_signal<T, P>(signal: &[Sample], sample_rate: SampleRate, overlap: usize, window: &Window, processor: P)
                            -> Vec<Sample> where T: FourierTransform, P: PhaseVocoderProcessor {
	let chunk_size = ::num_cpus::get() * 64;
	let frame_step_size = window.width() - overlap;
	let phase_step = phase_step(frame_step_size, window);
	let bin_count = fourier_transform::bin_count(window.width());
	let bin_width = bin_frequency::bin_width(sample_rate, bin_count);
	let overlap_factor = overlap_factor(overlap, window);
	let analyser = short_time_fourier::ShortTimeAnalyser::<T>::new(signal, overlap, window);

	let mut bin_phase_accumulate = vec![0.0; bin_count];
	let mut previous_chunk_frame = vec![Bin(Polar::default()); bin_count];
	let window = &window.normalize_amplitude(overlap_factor);
	let mut synthesiser = short_time_fourier::ShortTimeSynthesiser::<T>::new(overlap, window);
	let mut samples = Vec::new();

	let mut chunk_frame_index = 0;
	while chunk_frame_index < analyser.total_frames() {
		let chunk_end = usize::min(analyser.total_frames(), chunk_frame_index + chunk_size);
		let index_end = chunk_end - chunk_frame_index;
		debug_assert!(index_end <= chunk_size);

		let polar_frames: Vec<Vec<_>> = (chunk_frame_index..chunk_end).into_par_iter().map(|frame_index| {
			let frame = analyser.calculate_frame(frame_index);
			frame.iter().map(|bin| Bin(Polar::from(bin.take()))).collect()
		}).collect();

		let processed_frames: Vec<_> = (0..index_end).into_par_iter().map(|frame_index| {
			let frame = (0..bin_count).map(|bin_index| {
				let bin = polar_frames[frame_index][bin_index];
				let previous_frame_bin = if frame_index > 0 {
					polar_frames[frame_index - 1][bin_index]
				} else {
					previous_chunk_frame[bin_index]
				};
				let frequency = true_frequency(phase_step, bin_width, bin_index, &bin, &previous_frame_bin);
				PhaseVocoderBin { magnitude: bin.magnitude, frequency }
			}).collect();

			let mut output_bins = vec![PhaseVocoderBin::default(); bin_count];
			processor.process(frame, &mut output_bins);
			output_bins
		}).collect();

		let frames = vec![vec![Bin(Rectangular::default()); bin_count]; index_end];
		bin_phase_accumulate.par_iter_mut().enumerate().for_each(|(bin_index, accumulate)| {
			for (frame_index, frame) in processed_frames.iter().enumerate() {
				let bin = &frame[bin_index];
				let phase = original_phase(overlap_factor, phase_step, bin_width, bin_index,
				                           bin.frequency, *accumulate);
				let pointer = &frames[frame_index][bin_index] as *const _ as *mut Bin<Rectangular>;
				unsafe { *pointer = Bin(Polar { magnitude: bin.magnitude, phase }.into()); }
				*accumulate = phase;
			}
		});
		synthesiser.push_frames(&frames);

		samples.append(&mut synthesiser.flush_ready());
		previous_chunk_frame = polar_frames.last().unwrap().to_vec();
		chunk_frame_index += chunk_size;
	}

	samples.append(&mut synthesiser.flush_all());
	samples
}

pub fn overlap_factor(overlap: usize, window: &Window) -> f64 {
	window.width() as f64 / (window.width() - overlap) as f64
}

/// Constant factor for converting phase angles to phase shift in samples
pub fn phase_step(frame_step_size: usize, window: &Window) -> f64 {
	2.0 * consts::PI * (frame_step_size as f64 / window.width() as f64)
}

/// Estimates the true frequency of a bin by analysing the change in phase
/// of a bin across each frame
pub fn true_frequency(phase_step: f64, bin_width: f64, bin_index: usize,
                      bin: &Bin<Polar>, previous_frame_bin: &Bin<Polar>) -> Hertz {
	let phase_difference = bin.phase - previous_frame_bin.phase;
	let expected_phase_difference = bin_index as f64 * phase_step;
	let actual_phase_difference = phase_difference - expected_phase_difference;

	let phase_offset = map_into_circular_interval(actual_phase_difference);
	let center_frequency_difference = phase_offset / phase_step;
	(bin_index as f64 + center_frequency_difference) * bin_width
}

/// Calculates the original phase of the bin by working backwards from the
/// true frequency
pub fn original_phase(overlap_factor: f64, phase_step: f64, bin_width: f64, bin_index: usize,
                      true_frequency: f64, bin_phase_accumulate: f64) -> f64 {
	let center_frequency_difference = (true_frequency / bin_width) - bin_index as f64;
	let actual_phase_difference = 2.0 * consts::PI * (center_frequency_difference / overlap_factor);
	let phase_difference = actual_phase_difference + (bin_index as f64 * phase_step);
	bin_phase_accumulate + phase_difference
}

pub fn map_into_circular_interval(value: f64) -> f64 {
	let mut multiplier = (value / consts::PI) as i64;
	if multiplier >= 0 {
		multiplier += multiplier & 1;
	} else {
		multiplier -= multiplier & 1;
	}
	value - (consts::PI * multiplier as f64)
}

#[cfg(test)]
mod tests {
	use crate::window;
	use super::*;
	use test::Bencher;

	#[test]
	fn test_overlap_factor() {
		let window = Window::generate::<window::Dirichlet>(512);
		assert_eq!(overlap_factor(256, &window), 2.0);
		assert_eq!(overlap_factor(384, &window), 4.0);
	}

	#[test]
	fn test_phase_step() {
		let frame_step_size = 256;
		let window = Window::generate::<window::Dirichlet>(512);
		assert_eq!(phase_step(frame_step_size, &window), consts::PI);
	}

	#[test]
	fn test_true_frequency() {
		let sample_rate = 44100;
		let overlap = 256;
		let window = Window::generate::<window::Dirichlet>(512);

		let frame_step_size = window.width() - overlap;
		let bin_count = fourier_transform::bin_count(window.width());
		let bin_width = bin_frequency::bin_width(sample_rate, bin_count);
		let phase_step = phase_step(frame_step_size, &window);

		let previous_bin = Polar { magnitude: 5.0, phase: 0.0 }.into();
		let bin = Polar { magnitude: 5.0, phase: 1.0 }.into();
		assert_eq!(true_frequency(phase_step, bin_width, 0, &bin, &previous_bin), 27.31024509864819);

		let previous_bin = Polar { magnitude: 10.0, phase: 3.0 }.into();
		let bin = Polar { magnitude: 20.0, phase: 4.0 }.into();
		assert_eq!(true_frequency(phase_step, bin_width, 10, &bin, &previous_bin), 885.2868987951463);

		let previous_bin = Polar { magnitude: 100.0, phase: 35.0 }.into();
		let bin = Polar { magnitude: 20.0, phase: 35.0 }.into();
		assert_eq!(true_frequency(phase_step, bin_width, 40, &bin, &previous_bin), 3431.906614785992);
	}

	#[test]
	fn test_original_phase() {
		let sample_rate = 44100;
		let overlap = 256;
		let window = Window::generate::<window::Dirichlet>(512);
		let overlap_factor = overlap_factor(overlap, &window);

		let frame_step_size = window.width() - overlap;
		let bin_count = fourier_transform::bin_count(window.width());
		let bin_width = bin_frequency::bin_width(sample_rate, bin_count);
		let phase_step = phase_step(frame_step_size, &window);

		let accumulate_phase = 200.0;
		let true_frequency = 3000.0;
		assert_eq!(original_phase(overlap_factor, phase_step, bin_width, 40,
		                          true_frequency, accumulate_phase), 309.84888598266355);

		let accumulate_phase = 100.0;
		let true_frequency = 2500.0;
		assert_eq!(original_phase(overlap_factor, phase_step, bin_width, 30,
		                          true_frequency, accumulate_phase), 191.54073831888627);
	}

	#[bench]
	#[cfg(feature = "fast_fourier")]
	fn bench_pitch_shifter(bench: &mut Bencher) {
		use crate::fast_fourier::FastFourier;
		let sample_rate = 44100;
		let overlap = 240;
		let window = Window::generate::<window::Sine>(256);
		let processor = PitchShift { pitch_shift_ratio: 1.5 };
		let signal: Vec<_> = (0..22050).map(|x| x as f64).collect();
		bench.iter(|| process_signal::<FastFourier, _>(&signal, sample_rate, overlap,
		                                               &window, processor.clone()));
	}
}
