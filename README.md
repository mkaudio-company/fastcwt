# fastcwt
[![](https://img.shields.io/crates/v/fastcwt.svg)](https://crates.io/crates/fastcwt)
[![](https://img.shields.io/crates/l/fastcwt.svg)](https://crates.io/crates/fastcwt)
[![](https://docs.rs/fastcwt/badge.svg)](https://docs.rs/fastcwt/)

Rust-lang Continuous Wavelet Transform(CWT) library inspired by fCWT.

This crate is a direct translation of fCWT Library written in C++ by Arts, L.P.A. and van den Broek, E.L. (https://github.com/fastlib/fCWT)

I changed certain functions that I cannot translate, it seems like it could be eliminated, or the difference between fftw Library and rustfft crate. (fCWT used fftw, and fastcwt used rustfft.)

# Usage



# Changelog
0.1.5 - Parallelized FFT using rayon crate.

0.1.4 - Added error messages in assert!() and #![forbid(unsafe_code)] macro.

0.1.3 - Get rid of unsafe codes and fn find2power().

0.1.1, 0.1.2 - Minor fixes.

0.1.0 - Initial release.

# Citation
Arts, L.P.A., van den Broek, E.L. The fast continuous wavelet transformation (fCWT) for real-time, high-quality, noise-resistant time–frequency analysis. Nat Comput Sci 2, 47–58 (2022). https://doi.org/10.1038/s43588-021-00183-z