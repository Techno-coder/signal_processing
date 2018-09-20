use std::f64::consts;
use std::iter;
use super::Sample;

/// Computes the inverse discrete fourier transform of a rectangular pulse in the frequency domain
pub fn rectangular_pulse(cutoff_frequency: f64, signal_length: usize) -> Vec<Sample> {
	iter::once(1.0).chain((1..signal_length)
		.map(|i| (2.0 * cutoff_frequency * i as f64 * consts::PI).sin() / (i as f64 * consts::PI))
	).collect()
}