use super::Sample;

pub fn approximate(sample: Sample) -> Sample {
	const ZERO_EPSILON: f64 = 0.0000000000001;
	if -ZERO_EPSILON < sample && sample < ZERO_EPSILON {
		0.0
	} else {
		sample as f32 as f64
	}
}
