use crate::polar::Polar;

#[derive(Debug, PartialOrd, PartialEq, Copy, Clone)]
pub struct Rectangular {
	pub cosine: f64,
	pub sine: f64,
}

impl From<Polar> for Rectangular {
	fn from(other: Polar) -> Self {
		Rectangular {
			cosine: other.magnitude * other.phase.cos(),
			sine: other.magnitude * other.phase.sin(),
		}
	}
}
