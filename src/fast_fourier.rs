use crate::bin::Bin;
use crate::fourier_transform;
use crate::fourier_transform::FourierTransform;
use crate::rectangular::Rectangular;
use crate::utility;
use num_complex::Complex64;
use rustfft::FFTplanner;
use super::Sample;

pub struct FastFourier();

impl FourierTransform for FastFourier {
	fn analysis_extend(signal: &[Sample], signal_length: usize) -> Vec<Bin<Rectangular>> {
		let upper_bound = fourier_transform::bin_count(signal_length);
		let mut signal: Vec<_> = signal.iter().map(|real| Complex64::new(*real, 0.0)).collect();
		utility::pad_default(&mut signal, signal_length);
		let mut output = vec![Complex64::default(); signal_length];

		let planner = FFTplanner::new(false).plan_fft(signal_length);
		planner.process(&mut signal, &mut output);
		output.into_iter().take(upper_bound)
		      .map(|complex| Rectangular { cosine: complex.re, sine: complex.im }.into())
		      .collect()
	}

	fn synthesis(bins: &[Bin<Rectangular>], signal_length: usize) -> Vec<Sample> {
		let mirror_bound = fourier_transform::bin_count(signal_length - 1);
		let mut bins: Vec<_> = bins.iter().map(|complex| Complex64::new(complex.cosine, complex.sine)).collect();
		utility::pad_default(&mut bins, mirror_bound);
		(1..mirror_bound).rev().for_each(|index| bins.push(bins[index].conj()));
		let mut output = vec![Complex64::default(); signal_length];

		let planner = FFTplanner::new(true).plan_fft(signal_length);
		planner.process(&mut bins, &mut output);
		output.into_iter().map(|complex| complex.re / signal_length as f64).collect()
	}
}

#[cfg(test)]
mod tests {
	use crate::convolution;
	use crate::correlation;
	use crate::math;
	use super::*;
	use test::Bencher;

	#[test]
	fn test_synthesis() {
		let bins = [
			Bin(Rectangular { cosine: 15.0, sine: 0.0 }),
			Bin(Rectangular { cosine: -2.5, sine: 3.4409548 }),
			Bin(Rectangular { cosine: -2.5, sine: 0.81229924 }),
		];
		let synthesis: Vec<_> = FastFourier::synthesis(&bins, 5)
			.into_iter().map(math::approximate).collect();
		assert_eq!(synthesis, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
	}

	#[test]
	fn test_synthesis_even_length() {
		let bins = [
			Bin(Rectangular { cosine: 10.0, sine: 0.0 }),
			Bin(Rectangular { cosine: -2.0, sine: 2.0 }),
			Bin(Rectangular { cosine: -2.0, sine: 0.0 }),
		];
		let synthesis: Vec<_> = FastFourier::synthesis(&bins, 4);
		assert_eq!(&synthesis, &[1.0, 2.0, 3.0, 4.0]);
	}

	#[test]
	fn test_analysis() {
		let signal = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
		let bins = FastFourier::analysis(&signal);
		let synthesis: Vec<_> = FastFourier::synthesis(&bins, 6)
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

	#[bench]
	fn bench_analysis_and_synthesis(bench: &mut Bencher) {
		let signal: Vec<_> = (0..8192).map(|x| x as f64).collect();
		bench.iter(|| FastFourier::synthesis(&FastFourier::analysis(&signal), 8192));
	}
}