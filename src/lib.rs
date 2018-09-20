#![feature(self_struct_ctor)]

#[cfg(feature = "fast_fourier")]
extern crate rustfft;
#[cfg(feature = "fast_fourier")]
extern crate num_complex;

#[macro_use]
pub mod type_wrapper;
pub mod binned_histogram;
pub mod convolution;
pub mod correlation;
pub mod filter_kernels;
pub mod histogram;
pub mod running_statistics;
pub mod statistics;
pub mod fourier_transform;
pub mod math;
pub mod polar;
pub mod frequency;
pub mod rectangular;
pub mod transform_pairs;
pub mod utility;

#[cfg(feature = "fast_fourier")]
pub mod fast_fourier;

pub type Sample = f64;
pub type IntegralSample = i64;
