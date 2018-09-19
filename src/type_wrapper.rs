macro_rules! wrapper {
    ($type: ident) => {
		impl<T> $type<T> {
			pub fn take(self) -> T {
				self.0
			}
		}

		impl<T> ops::Deref for $type<T> {
			type Target = T;

			fn deref(&self) -> &<Self as ops::Deref>::Target {
				&self.0
			}
		}

		impl<T> ops::DerefMut for $type<T> {
			fn deref_mut(&mut self) -> &mut <Self as ops::Deref>::Target {
				&mut self.0
			}
		}

		impl<T> From<T> for $type<T> {
			fn from(rhs: T) -> Self {
				$type(rhs)
			}
		}
    };
}

macro_rules! cross_cast {
	($type: ident, $a: ident, $b: ident) => {
		impl From<$type<$b>> for $type<$a> {
			fn from(other: $type<$b>) -> Self {
				$type(other.take().into())
			}
		}

		impl From<$type<$a>> for $type<$b> {
			fn from(other: $type<$a>) -> Self {
				$type(other.take().into())
			}
		}
	}
}
