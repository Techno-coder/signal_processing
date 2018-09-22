pub fn pad_default<T>(vector: &mut Vec<T>, length: usize) where T: Default {
	while vector.len() < length {
		vector.push(T::default());
	}
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
}