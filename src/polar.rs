use std::f64::consts;
use super::Sample;

/// Converts rectangular to polar form as (Magnitude, Phase)
pub fn to_polar(real: Sample, imaginary: Sample) -> (Sample, Sample) {
	let magnitude = (real * real + imaginary * imaginary).sqrt();
	let phase = if real == 0.0 {
		if imaginary.is_sign_positive() {
			consts::FRAC_PI_2
		} else {
			-consts::FRAC_PI_2
		}
	} else {
		let phase = (imaginary / real).atan();
		if real.is_sign_negative() {
			if imaginary.is_sign_negative() {
				phase - consts::PI
			} else {
				phase + consts::PI
			}
		} else {
			phase
		}
	};
	(magnitude, phase)
}

/// Allow the magnitude to be negative to resolve cutting of the phase
pub fn unwrap_phase(previous_unwrapped: Sample, phase: Sample) -> Sample {
	let multiplier = (previous_unwrapped - phase) / (2.0 * consts::PI);
	phase + multiplier.trunc() * 2.0 * consts::PI
}

/// Converts polar to rectangular form as (Real, Imaginary)
pub fn to_rectangular(magnitude: Sample, phase: Sample) -> (Sample, Sample) {
	(magnitude * phase.cos(), magnitude * phase.sin())
}

pub fn to_polar_spectrum(cosines: &[Sample], sines: &[Sample]) -> (Vec<Sample>, Vec<Sample>) {
	assert_eq!(cosines.len(), sines.len());
	let mut magnitudes = Vec::new();
	let mut phases = Vec::new();
	for (cosine, sine) in cosines.iter().zip(sines.iter()) {
		let (magnitude, mut phase) = to_polar(*cosine, *sine);
		if !phases.is_empty() {
			phase = unwrap_phase(*phases.last().unwrap(), phase);
		}
		magnitudes.push(magnitude);
		phases.push(phase);
	}
	(magnitudes, phases)
}

#[cfg(test)]
mod tests {
	use crate::math;
	use super::*;

	#[test]
	fn test_to_polar() {
		let (magnitude, phase) = to_polar(3.0, 4.0);
		assert_eq!(magnitude, 5.0);
		assert_eq!(math::approximate(phase), 0.9272952079772949);
		let (_, phase) = to_polar(-3.0, 4.0);
		assert_eq!(math::approximate(phase), 2.2142975330352783);

		let (_, phase) = to_polar(5.0, 0.0);
		assert_eq!(phase, 0.0);
		let (_, phase) = to_polar(-5.0, 0.0);
		assert_eq!(phase, consts::PI);
	}
}
