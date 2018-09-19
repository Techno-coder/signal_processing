use crate::rectangular::Rectangular;
use std::f64::consts;
use super::Sample;

#[derive(Debug, PartialOrd, PartialEq, Copy, Clone)]
pub struct Polar {
	pub magnitude: f64,
	pub phase: f64,
}

impl Polar {
	/// Allow the magnitude to be negative to resolve cutting of the phase
	pub fn unwrap_phase(&mut self, previous_phase: f64) {
		let multiplier = (previous_phase - self.phase) / (2.0 * consts::PI);
		self.phase += multiplier.trunc() * 2.0 * consts::PI;
	}
}

impl From<Rectangular> for Polar {
	fn from(other: Rectangular) -> Self {
		let magnitude = (other.cosine * other.cosine + other.sine * other.sine).sqrt();
		let phase = if other.cosine == 0.0 {
			if other.sine.is_sign_positive() {
				consts::FRAC_PI_2
			} else {
				-consts::FRAC_PI_2
			}
		} else {
			let phase = (other.sine / other.cosine).atan();
			if other.cosine.is_sign_negative() {
				if other.sine.is_sign_negative() {
					phase - consts::PI
				} else {
					phase + consts::PI
				}
			} else {
				phase
			}
		};
		Polar { magnitude, phase }
	}
}

pub fn to_polar_spectrum(cosines: &[Sample], sines: &[Sample]) -> (Vec<Sample>, Vec<Sample>) {
	assert_eq!(cosines.len(), sines.len());
	let mut magnitudes = Vec::new();
	let mut phases = Vec::new();
	for (cosine, sine) in cosines.iter().zip(sines.iter()) {
		let mut polar: Polar = Rectangular { cosine: *cosine, sine: *sine }.into();
		if !phases.is_empty() {
			polar.unwrap_phase(*phases.last().unwrap());
		}
		magnitudes.push(polar.magnitude);
		phases.push(polar.phase);
	}
	(magnitudes, phases)
}

#[cfg(test)]
mod tests {
	use crate::math;
	use super::*;

	#[test]
	fn test_into_polar() {
		let polar: Polar = Rectangular { cosine: 3.0, sine: 4.0 }.into();
		assert_eq!(polar.magnitude, 5.0);
		assert_eq!(math::approximate(polar.phase), 0.9272952079772949);
		let polar: Polar = Rectangular { cosine: -3.0, sine: 4.0 }.into();
		assert_eq!(math::approximate(polar.phase), 2.2142975330352783);

		let polar: Polar = Rectangular { cosine: 5.0, sine: 0.0 }.into();
		assert_eq!(polar.phase, 0.0);
		let polar: Polar = Rectangular { cosine: -5.0, sine: 0.0 }.into();
		assert_eq!(polar.phase, consts::PI);
	}
}
