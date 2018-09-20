use crate::fourier_transform::FourierTransform;
use crate::frequency::Frequency;
use crate::rectangular::Rectangular;
use crate::utility;
use num_complex::Complex64;
use rustfft::FFTplanner;
use super::Sample;

pub struct FastFourier();

impl FourierTransform for FastFourier {
	fn analysis_extend(signal: &[Sample], signal_length: usize) -> Vec<Frequency<Rectangular>> {
		let mut signal: Vec<_> = signal.iter().map(|real| Complex64::new(*real, 0.0)).collect();
		utility::pad_default(&mut signal, signal_length);
		let mut output = vec![Complex64::default(); signal_length];

		let planner = FFTplanner::new(false).plan_fft(signal_length);
		planner.process(&mut signal, &mut output);
		output.into_iter().take((signal_length + 1) / 2)
		      .map(|complex| Rectangular { cosine: complex.re, sine: complex.im }.into())
		      .collect()
	}

	fn synthesis(frequencies: &[Frequency<Rectangular>], signal_length: usize) -> Vec<Sample> {
		let mut frequencies: Vec<_> = frequencies.iter().map(|complex| Complex64::new(complex.cosine, complex.sine)).collect();
		(1..frequencies.len()).rev().for_each(|index| frequencies.push(frequencies[index].conj()));
		utility::pad_default(&mut frequencies, signal_length);
		let mut output = vec![Complex64::default(); signal_length];

		let planner = FFTplanner::new(true).plan_fft(signal_length);
		planner.process(&mut frequencies, &mut output);
		output.into_iter()
		      .map(|complex| complex.re)
		      .map(|sample| sample * (1.0 / signal_length as f64))
		      .collect()
	}
}

#[cfg(test)]
mod tests {
	use crate::convolution;
	use crate::correlation;
	use crate::math;
	use super::*;

	#[test]
	fn test_synthesis() {
		let frequencies = [
			Frequency(Rectangular { cosine: 15.0, sine: 0.0 }),
			Frequency(Rectangular { cosine: -2.5, sine: 3.4409548 }),
			Frequency(Rectangular { cosine: -2.5, sine: 0.81229924 }),
		];
		let synthesis: Vec<_> = FastFourier::synthesis(&frequencies, 5)
			.into_iter().map(math::approximate).collect();
		assert_eq!(synthesis, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
	}

	#[test]
	fn test_analysis() {
		let signal = [1.0, 2.0, 3.0, 4.0, 5.0];
		let frequencies = FastFourier::analysis(&signal);
		let synthesis: Vec<_> = FastFourier::synthesis(&frequencies, 5)
			.into_iter().map(math::approximate).collect();
		assert_eq!(synthesis, signal);
	}

	#[test]
	fn test_convolution() {
		let signal = [0.0, 1.0, 2.0, 3.0, 2.0, 0.0];
		let impulse_response = [1.0, 2.0];
		let convolution: Vec<f64> = convolution::convolve_fourier::<FastFourier>(&signal, &impulse_response)
			.into_iter().map(math::approximate).collect();
		assert_eq!(convolution, vec![0.0, 1.0, 2.0 + 2.0, 4.0 + 3.0, 6.0 + 2.0, 4.0, 0.0])
	}

	#[test]
	fn test_correlation() {
		let signal = [0.0, 1.0, 2.0, 5.0, 2.0, 1.0];
		let target = [1.0, 2.0];
		let correlation = correlation::correlation::<FastFourier>(&signal, &target);
		assert_eq!(math::approximate(correlation), vec![2.0, 5.0, 12.0, 9.0, 4.0, 1.0].into_iter().sum());
	}
}