use std::f64::consts;
use super::Sample;

pub fn approximate(sample: Sample) -> Sample {
	const ZERO_EPSILON: f64 = 0.0000000000001;
	if -ZERO_EPSILON < sample && sample < ZERO_EPSILON {
		0.0
	} else {
		sample as f32 as f64
	}
}

pub fn sinc(x: f64) -> f64 {
	if x != 0.0 {
		let x = x * consts::PI;
		x.sin() / x
	} else {
		1.0
	}
}
