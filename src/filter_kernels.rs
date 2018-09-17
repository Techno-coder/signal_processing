use super::Sample;

pub fn low_pass_square_pulse(amplitude: f64, width: usize) -> Vec<Sample> {
	vec![amplitude; width]
}

pub fn high_pass_square_pulse(delta_amplitude: f64, width: usize) -> Vec<Sample> {
	assert!(width > 1);
	let dampen_amplitude = delta_amplitude / (width - 1) as f64;
	let mut filter = vec![dampen_amplitude; width];
	filter[width / 2] = delta_amplitude;
	filter
}
