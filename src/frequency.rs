use super::Sample;

/// Multiplies two frequency domain components together
pub fn multiply_rectangular((a_real, a_imaginary): (Sample, Sample),
                            (b_real, b_imaginary): (Sample, Sample)) -> (Sample, Sample) {
	(a_real * b_real - a_imaginary * b_imaginary,
	 a_imaginary * b_real + a_real * b_imaginary)
}

/// Calculates `a / b` where `a` and `b` are part of the frequency domain
pub fn divide_rectangular((a_real, a_imaginary): (Sample, Sample),
                          (b_real, b_imaginary): (Sample, Sample)) -> (Sample, Sample) {
	let denominator = b_real * b_real + b_imaginary * b_imaginary;
	let real = (a_real * b_real + a_imaginary * b_imaginary) / denominator;
	let imaginary = (a_imaginary * b_real - a_real * b_imaginary) / denominator;
	(real, imaginary)
}

pub fn multiply_polar((a_magnitude, a_phase): (Sample, Sample),
                      (b_magnitude, b_phase): (Sample, Sample)) -> (Sample, Sample) {
	(a_magnitude * b_magnitude, a_phase + b_phase)
}

/// Calculates `a / b` where `a` and `b` are part of the frequency domain in polar form
pub fn divide_polar((a_magnitude, a_phase): (Sample, Sample),
                    (b_magnitude, b_phase): (Sample, Sample)) -> (Sample, Sample) {
	(a_magnitude / b_magnitude, a_phase - b_phase)
}

#[cfg(test)]
mod tests {
	use crate::polar;
	use super::*;

	#[test]
	fn test_multiply_rectangular() {
		assert_eq!(multiply_rectangular((1.0, 2.0), (3.0, 4.0)), (-5.0, 10.0));
	}

	#[test]
	fn test_divide_rectangular() {
		let a = (1.0, 2.0);
		let b = (3.0, 4.0);
		let c = multiply_rectangular(a, b);
		assert_eq!(divide_rectangular(c, a), b);
		assert_eq!(divide_rectangular(c, b), a);
	}

	#[test]
	fn test_multiply_polar() {
		let a = (1.0, 2.0);
		let b = (3.0, 4.0);
		let polar_a = polar::to_polar(a.0, a.1);
		let polar_b = polar::to_polar(b.0, b.1);
		let (magnitude, phase) = multiply_polar(polar_a, polar_b);
		let polar_rectangular = polar::to_rectangular(magnitude, phase);
		assert_eq!(polar_rectangular, multiply_rectangular(a, b));
	}

	#[test]
	fn test_divide_polar() {
		let a = (1.0, 2.0);
		let b = (3.0, 4.0);
		let c = multiply_polar(a, b);
		assert_eq!(divide_polar(c, a), b);
		assert_eq!(divide_polar(c, b), a);
	}
}
