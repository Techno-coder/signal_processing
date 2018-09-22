#![feature(self_struct_ctor)]

#[cfg(test)]
#[macro_use]
extern crate pretty_assertions;

#[cfg(feature = "fast_fourier")]
extern crate num_complex;
#[cfg(feature = "fast_fourier")]
extern crate rustfft;

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
pub mod bin;
pub mod rectangular;
pub mod transform_pairs;
pub mod utility;
pub mod window;
pub mod short_time_fourier;
pub mod phase_vocoder;
pub mod bin_frequency;

#[cfg(feature = "fast_fourier")]
pub mod fast_fourier;

pub type Sample = f64;
pub type IntegralSample = i64;
pub type Hertz = f64;
pub type SampleRate = u64;
