#![feature(self_struct_ctor)]

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

pub type Sample = f64;
pub type IntegralSample = i64;
