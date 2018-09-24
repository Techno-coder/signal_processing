#![feature(self_struct_ctor)]
#![feature(test)]

#[cfg(feature = "fast_fourier")]
extern crate num_complex;
#[cfg(test)]
#[macro_use]
extern crate pretty_assertions;
extern crate rayon;
#[cfg(feature = "fast_fourier")]
extern crate rustfft;
extern crate test;
extern crate num_cpus;

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
