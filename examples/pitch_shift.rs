extern crate hound;
extern crate signal_processing;

use signal_processing::fast_fourier::FastFourier;
use signal_processing::window;
use std::i32;

fn main() {
	let mut reader = hound::WavReader::open("input_file.wav").unwrap();
	let spec = reader.spec();
	let input_samples: Vec<i32> = reader.samples::<i32>().map(|x| x.unwrap()).collect();
	let signal: Vec<_> = input_samples.iter().map(|x| *x as f64).collect();

	use signal_processing::phase_vocoder;
	let pitch_shifter = phase_vocoder::PitchShift { pitch_shift_ratio: 1.2 };
	let window = window::Window::generate::<window::Sine>(512);
	let signal = phase_vocoder::process_signal::<FastFourier, _>(&signal, spec.sample_rate.into(), 448, &window, pitch_shifter);
	let output_samples: Vec<_> = signal.iter().map(|x| *x as i16).collect();

	let mut writer = hound::WavWriter::create("pitch_shifted_file.wav", spec).unwrap();
	output_samples.iter().for_each(|sample| writer.write_sample(*sample).unwrap());
	writer.finalize().unwrap();
}