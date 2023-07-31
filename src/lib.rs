//! This crate is a direct translation of fCWT Library written in C++ by Arts, L.P.A. and van den Broek, E.L.
//! (<https://github.com/fastlib/fCWT>)
//!
//! I changed certain functions that I cannot translate, it seems like it could be eliminated,
//! or the difference between fftw Library and rustfft crate.
//! (fCWT used fftw, and fastcwt used rustfft.)
//!
//! ### Usage
//! ```
//! use fastcwt::*;
//! use rand::prelude::*;
//!
//! let wavelet = Wavelet::create(1.0); //Create a Morlet wavelet.
//! let scale = Scales::create(ScaleTypes::LinFreq, 48000, 20.0, 20000.0, 1000); //Create a scale factor.
//!
//! let mut transform = FastCWT::create(wavelet, true); // Create a fCWT instance.
//!
//! let mut input = vec![];
//! for _ in 0 .. 48000
//! {
//!     input.push(thread_rng().gen_range(-1.0 .. 1.0))
//! };
//!
//! let result = transform.cwt(1000, input.as_slice(), scale); //Store the result.
//! ```
//!
//! Changelog
//!
//! 0.1.7 - Used boxed slice instead of vec in Scales struct.
//! 
//! 0.1.6 - Transfered owndership to company account.
//! 
//! 0.1.5 - Parallelized FFT using rayon crate.
//!
//! 0.1.4 - Added error messages in assert!() and #![forbid(unsafe_code)] macro.
//!
//! 0.1.3 - Get rid of unsafe codes and fn find2power().
//!
//! 0.1.1, 0.1.2 - Minor fixes.
//!
//! 0.1.0 - Initial release.
//!
//! ### Citation
//!
//! Arts, L.P.A., van den Broek, E.L. The fast continuous wavelet transformation (fCWT) for real-time, high-quality, noise-resistant time–frequency analysis. Nat Comput Sci 2, 47–58 (2022). <https://doi.org/10.1038/s43588-021-00183-z>

#![feature(core_intrinsics)]
#![forbid(unsafe_code)]

use rustfft;
use rayon::prelude::*;

/// Scale types selection for Scale object.
#[derive(PartialEq)]
pub enum ScaleTypes
{
    /// Linear scale.
    Linear,
    /// Logarithmic scale.
    Log,
    /// Linear for frequency.
    LinFreq
}
/// Morlet wavelet object.
pub struct Wavelet
{
    width : usize,
    imag_freq : bool,
    double_sided : bool,
    mother : Vec<f64>,
    fb : f64
}
impl Wavelet
{
    /// Create a wavelet object.
    ///
    /// bandwidth           - bandwidth of the Morlet wavelet
    pub fn create(bandwidth : f64) -> Wavelet
    {
        return Wavelet
        {
            width : 0,
            imag_freq : false,
            double_sided : false,
            mother : vec![],
            fb : bandwidth
        }
    }
    fn generate(& mut self, size : usize)
    {
        //Frequency domain, because we only need size. Default scale is always 2;
        self.width = size;

        let toradians = 2.0 * std::f64::consts::PI / size as f64;
        let norm = (2.0 * std::f64::consts::PI).sqrt() * (1.0 / std::f64::consts::PI).powf(0.25);

        //calculate array
        for w in 0 .. self.width
        {
            let mut tmp1 = 2.0 * (w as f64 * toradians) * self.fb - 2.0 * std::f64::consts::PI * self.fb;
            tmp1 = - tmp1.powf(2.0) / 2.0;
            self.mother.push(norm * (tmp1).exp());
        }
    }
}

/// Scale factor for the wavelet transform.
pub struct Scales
{
    scales : Box<[f64]>,
    fs : usize,
    num_scales : usize
}
impl Scales
{
    /// Create the scale factor for the transform.
    ///
    /// st                  - Log | Linear for logarithmic or linear distribution of scales across frequency range
    ///
    /// afs                 - Sample frequency
    ///
    /// af0                 - Beginning of the frequency range
    ///
    /// af1                 - End of the frequency range
    ///
    /// af_num              - Number of wavelets to generate across frequency range
    pub fn create(st : ScaleTypes, afs : usize, af0 : f64, af1 : f64, af_num : usize) -> Scales
    {
        let mut scales = Scales
        {
            scales: vec![0.0;af_num].into_boxed_slice(),
            fs : afs,
            num_scales: af_num,
        };
        match st
        {
            ScaleTypes::Linear => { scales.calculate_linscale_array(afs, af0, af1, af_num); }
            ScaleTypes::Log => { scales.calculate_logscale_array(2.0, afs, af0, af1, af_num); }
            ScaleTypes::LinFreq => { scales.calculate_linfreq_array(afs, af0, af1, af_num); }
        }
        return scales;
    }
    pub fn get_scales(& self) -> Box<[f64]> { return self.scales.clone(); }
    pub fn get_frequencies(& self, p_freqs : & mut Vec<f64>) -> Vec<f64>
    {
        let mut frequencies = vec![];
        for i in 0..p_freqs.len() { frequencies.push(self.fs as f64 / self.scales[i]); }
        return frequencies;
    }
    fn calculate_logscale_array(& mut self, base : f64, fs : usize, f0 : f64, f1 : f64, f_num : usize)
    {
        let nf0 = f0;
        let nf1 = f1;
        let s0 = fs as f64 / nf1;
        let s1 = fs as f64 / nf0;

        //Cannot pass the nyquist frequency
        assert!(f1 as usize <= fs / 2, "Max frequency cannot be higher than the Nyquist frequency.");

        let power0 = s0.log(std::f64::consts::E) / base.log(std::f64::consts::E);
        let power1 = s1.log(std::f64::consts::E) / base.log(std::f64::consts::E);
        let dpower = power1 - power0;

        for i in 0 .. f_num
        {
            let power = power0 + dpower / ((f_num - 1) * i) as f64;
            self.scales[i] = base.powf(power);
        }
    }
    fn calculate_linscale_array(& mut self, fs : usize, f0 : f64, f1 : f64, f_num : usize)
    {
        //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;

        //Cannot pass the nyquist frequency
        assert!(f1 <= (fs / 2) as f64, "Max frequency cannot be higher than the Nyquist frequency.");
        let df = f1 - f0;

        for i in 0 .. f_num { self.scales[f_num - i - 1] = fs as f64 / f0 + (df / f_num as f64) * i as f64; }
    }
    fn calculate_linfreq_array(& mut self, fs : usize, f0 : f64, f1 : f64, f_num : usize)
    {
        //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
        let s0 = fs as f64 / f1;
        let s1 = fs as f64 / f0;

        //Cannot pass the nyquist frequency
        assert!(f1 <= fs as f64 / 2.0, "Max frequency cannot be higher than the Nyquist frequency.");
        let ds = s1 - s0;

        for i in 0 .. f_num { self.scales[i] = s0 + (ds / f_num as f64) * i as f64; }
    }
}

/// Actual continuous wavelet transform.
pub struct FastCWT
{
    wavelet : Wavelet,
    use_normalization : bool
}
impl FastCWT
{
    /// # Arguments
    /// wavelet             - Wavelet object.
    ///
    /// optplan             - Use FFT optimization plans if true.
    pub fn create(wavelet : Wavelet, optplan : bool) -> FastCWT { return FastCWT { wavelet, use_normalization : optplan, } }
    /// # Arguments
    /// input     - Input data in vector format
    ///
    /// scales    - Scales object
    pub fn cwt(& mut self, num : usize, input : & [f64], scales : Scales) -> Vec<rustfft::num_complex::Complex<f64>>
    {
        //Find nearest power of 2
        let newsize = num.next_power_of_two();
        let mut buffer = vec![];

        //Copy input to new input buffer
        for data in input { buffer.push(rustfft::num_complex::Complex::new(* data, 0.0)); }

        let mut planner = rustfft::FftPlanner::new();

        // //Perform forward FFT on input signal
        planner.plan_fft_forward(buffer.len()).process(& mut buffer);

        //Generate mother wavelet function
        self.wavelet.generate(newsize);
        for i in 1 .. newsize >> 1 { buffer[newsize - i] = buffer[i]; }

        let buffer = std::sync::Arc::new(std::sync::Mutex::new(Some(buffer)));
        let buffer_clone = buffer.clone();
        
        (0 .. scales.num_scales).into_par_iter().for_each(|i|
        {
            if let Ok(mut guard) = buffer_clone.try_lock()
            {
                if let Some(buffer_clone) = guard.as_mut()
                {
                    let mut planner = rustfft::FftPlanner::new();
                    //FFT-base convolution in the frequency domain
                    self.daughter_wavelet_multiplication(buffer_clone, self.wavelet.mother.as_slice(), scales.scales[i],num, self.wavelet.imag_freq, self.wavelet.double_sided);

                    planner.plan_fft_forward(buffer_clone.len()).process(buffer_clone);
                    if self.use_normalization { for n in 0 .. newsize { buffer_clone[n] = buffer_clone[n] / newsize as f64; } }
                }
            }
        });
        return buffer.lock().unwrap().take().unwrap();
    }
    fn daughter_wavelet_multiplication(& self, buffer : & mut [rustfft::num_complex::Complex<f64>], mother : &[f64], scale : f64, i_size : usize, imaginary : bool, doublesided : bool)
    {
        let endpoint = std::cmp::min((i_size as f64 / 2.0) as usize, (i_size as f64 * 2.0 / scale) as usize);
        let step = scale / 2.0;

        let maximum = i_size as f64 - 1.0;
        let s1 = i_size - 1;

        for n in 0 .. endpoint
        {
            let tmp = std::cmp::min(maximum as usize, step as usize * n);

            if doublesided
            {
                buffer[s1 - n].re = if imaginary { buffer[s1 - n].re * mother[tmp as usize] } else { buffer[s1 - n].re * mother[tmp as usize] * -1.0 };
                buffer[s1 - n].im = buffer[s1 - n].im * mother[tmp as usize];
            } else { buffer[n].re = buffer[n].re * mother[tmp as usize]; buffer[n].im = buffer[n].im * mother[tmp as usize]; }
        }
    }
}