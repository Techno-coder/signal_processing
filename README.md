# signal_processing
Rust library for various signal processing algorithms and structures

This is a work in progress library that adapts many algorithms from "The Scientist and Engineer's Guide to
Digital Signal Processing" book. I created this library with the ultimate goal of making a pitch shifter and
to document my learning in digital signal processing along the way.

## Algorithms and Data Structures
- Phase Vocoder
  - Pitch shifter
- Convolution
- Correlation
- Fourier transform

## Why 64 bit floats?
I've tested my library with both 32 bit floats and 64 bit floats and I've found that 32 bit floats have a
significant amount of rounding errors to the point that audio processed with 32 bit floats have a very noticeable
drop in quality. I'm sure that many of the rounding errors can be fixed by changing some floats to be 64 bit but
I believe that this is not worth the effort. Processors can do floating point arithmetic for both 32 bit and 64 bit
floats in around one clock cycle now so the only performance benefit gained would be through cache and memory usage.

## Usages of `unsafe`
There are several segments of `unsafe` code. Most notably is in the Phase Vocoder and Short Time Fourier Transform files.
This unsafe code mutates array elements that cannot be done normally as it is in a parallelised section. Without the
parallelisation the code runs 6 - 10 times slower so I believe this is definitely worth the "unsafeness" (`cargo bench`). 
Theoretically this should not cause any panics or have any race conditions.
