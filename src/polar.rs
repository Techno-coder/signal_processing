use crate::bin::Bin;
use crate::rectangular::Rectangular;
use std::f64::consts;
use super::Sample;

#[derive(Debug, Default, PartialOrd, PartialEq, Copy, Clone)]
pub struct Polar {
	pub magnitude: Sample,
	pub phase: Sample,
}

impl Polar {
	/// Allow the magnitude to be negative to resolve cutting of the phase
	pub fn unwrap_phase(&mut self, previous_phase: Sample) {
		let multiplier = (previous_phase - self.phase) / (2.0 * consts::PI);
		self.phase += multiplier.trunc() * 2.0 * consts::PI;
	}

	pub fn complex_conjugate(&self) -> Polar {
		Polar {
			magnitude: self.magnitude,
			phase: -self.phase,
		}
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

pub fn to_polar_spectrum(bins: &[Bin<Rectangular>]) -> Vec<Bin<Polar>> {
	let mut previous_phase = None;
	bins.iter().map(|bin| {
		let mut polar: Bin<Polar> = (*bin).into();
		if let Some(previous_phase) = previous_phase {
			polar.unwrap_phase(previous_phase);
		}
		previous_phase = Some(polar.phase);
		polar
	}).collect()
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
