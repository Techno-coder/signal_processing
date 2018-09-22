pub fn pad_default<T>(vector: &mut Vec<T>, length: usize) where T: Default {
	while vector.len() < length {
		vector.push(T::default());
	}
}

/// Returns the highest element in the array that forms a pyramid
pub fn find_peak<T>(data: &[T]) -> Option<&T> where T: PartialOrd {
	let mut increasing = true;
	let mut highest = data.first()?;
	for element in &data[1..] {
		if increasing {
			if element >= highest {
				highest = element;
			} else {
				increasing = false;
			}
		} else {
			if element >= highest {
				return None;
			}
		}
	}
	Some(highest)
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_pad_zeros() {
		let mut vector = Vec::new();
		pad_default(&mut vector, 10);
		assert_eq!(vector.len(), 10);
		assert!(vector.iter().all(|x: &f64| x == &0.0));
	}

	#[test]
	fn test_find_peak() {
		let data = [1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0];
		assert_eq!(find_peak(&data), Some(&4.0));

		let data = [1.0, 2.0, 3.0, 2.0, 1.0, 2.0, 3.0, 2.0, 1.0];
		assert_eq!(find_peak(&data), None);
	}
}