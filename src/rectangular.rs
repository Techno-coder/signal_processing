use crate::polar::Polar;
use super::Sample;

#[derive(Debug, Default, PartialOrd, PartialEq, Copy, Clone)]
pub struct Rectangular {
	pub cosine: Sample,
	pub sine: Sample,
}

impl From<Polar> for Rectangular {
	fn from(other: Polar) -> Self {
		Rectangular {
			cosine: other.magnitude * other.phase.cos(),
			sine: other.magnitude * other.phase.sin(),
		}
	}
}
