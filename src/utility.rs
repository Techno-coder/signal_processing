use super::*;

pub fn pad_zeros(vector: &mut Vec<Sample>, length: usize) {
	assert!(vector.len() <= length);
	while vector.len() != length {
		vector.push(0.0);
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_pad_zeros() {
		let mut vector = Vec::new();
		pad_zeros(&mut vector, 10);
		assert_eq!(vector.len(), 10);
		assert!(vector.iter().all(|x| x == &0.0));
	}
}