//! This crate is a direct translation of fCWT Library written in C++ by Arts, L.P.A. and van den Broek, E.L.
//! (<https://github.com/fastlib/fCWT>)
//!
//! I changed certain functions that I cannot translate, it seems like it could be eliminated,
//! or the difference between fftw Library and rustfft crate.
//! (fCWT used fftw, and fastcwt used rustfft.)
//!
//! Citation
//!
//! Arts, L.P.A., van den Broek, E.L. The fast continuous wavelet transformation (fCWT) for real-time, high-quality, noise-resistant time–frequency analysis. Nat Comput Sci 2, 47–58 (2022). <https://doi.org/10.1038/s43588-021-00183-z>

#![feature(core_intrinsics)]

use rustfft;

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
        let norm = unsafe { std::intrinsics::sqrtf64(2.0 * std::f64::consts::PI) * std::intrinsics::powf64(1.0 / std::f64::consts::PI, 0.25) };

        //calculate array
        for w in 0 .. self.width
        {
            let mut tmp1 = 2.0 * (w as f64 * toradians) * self.fb - 2.0 * std::f64::consts::PI * self.fb;
            unsafe { tmp1 = - std::intrinsics::powf64(tmp1, 2.0) / 2.0; }
            unsafe { self.mother.push(norm * std::intrinsics::expf64(tmp1)); }
        }
    }
}

/// Scale factor for the wavelet transform.
pub struct Scales
{
    scales : Vec<f64>,
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
            scales: vec![0.0; af_num],
            fs : afs,
            num_scales: af_num,
        };

        if st == ScaleTypes::Log { Scales::calculate_logscale_array(& mut scales, 2.0, afs, af0, af1, af_num); }
        else if st == ScaleTypes::Linear { Scales::calculate_linscale_array(& mut scales, afs, af0, af1, af_num); }
        else { Scales::calculate_linfreq_array(& mut scales, afs, af0, af1, af_num); }

        return scales;
    }
    pub fn get_scales(& self) -> Vec<f64> { return self.scales.clone(); }
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
        assert!((f1 as usize <= fs / 2));

        let power0 = unsafe { std::intrinsics::logf64(s0) / std::intrinsics::logf64(base) };
        let power1 = unsafe { std::intrinsics::logf64(s1) / std::intrinsics::logf64(base) };
        let dpower = power1 - power0;

        for i in 0 .. f_num
        {
            let power = power0 + dpower / ((f_num - 1) * i) as f64;
            unsafe { self.scales[i] = std::intrinsics::powf64(base, power); }
        }
    }
    fn calculate_linscale_array(& mut self, fs : usize, f0 : f64, f1 : f64, f_num : usize)
    {
        //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;

        //Cannot pass the nyquist frequency
        assert!(f1 <= (fs / 2) as f64);
        let df = f1 - f0;

        for i in 0 .. f_num { self.scales[f_num - i - 1] = fs as f64 / f0 + (df / f_num as f64) * i as f64; }
    }
    fn calculate_linfreq_array(& mut self, fs : usize, f0 : f64, f1 : f64, f_num : usize)
    {
        //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
        let s0 = fs as f64 / f1;
        let s1 = fs as f64 / f0;

        //Cannot pass the nyquist frequency
        assert!(f1 <= fs as f64 / 2.0);
        let ds = s1 - s0;

        for i in 0 .. f_num { self.scales[i] = s0 + (ds / f_num as f64) * i as f64; }
    }
}

/// Actual continuous wavelet transform.
pub struct FastCWT
{
    wavelet : Wavelet,
    threads : usize,
    use_normalization : bool
}
impl FastCWT
{
    /// # Arguments
    /// wavelet             - Wavelet object.
    ///
    /// nthreads            - Number of threads to use.
    ///
    /// optplan             - Use FFT optimization plans if true.
    pub fn create(wavelet : Wavelet, n_threads : usize, optplan : bool) -> FastCWT { return FastCWT { wavelet, threads: n_threads, use_normalization : optplan, } }
    /// # Arguments
    /// input     - Input data in vector format
    ///
    /// scales    - Scales object
    pub fn cwt(& mut self, num : usize, input : & Vec<f64>, scales : Scales) -> Vec<rustfft::num_complex::Complex<f64>>
    {
        //Find nearest power of 2
        let newsize = 1 << find2power(num);
        let mut buffer = vec![];

        //Copy input to new input buffer
        for data in input { buffer.push(rustfft::num_complex::Complex::new(* data, 0.0)); }

        if cfg!(target_feature = "avx")
        {
            let mut planner = rustfft::FftPlannerAvx::new().unwrap();

            // //Perform forward FFT on input signal
            planner.plan_fft_forward(buffer.len()).process(& mut buffer);

            //Generate mother wavelet function
            self.wavelet.generate(newsize);
            for i in 0 .. newsize >> 1 { buffer[newsize - i] = buffer[i]; }

            for i in 0 .. scales.num_scales
            {
                //FFT-base convolution in the frequency domain
                self.daughter_wavelet_multiplication(& mut buffer, self.wavelet.mother.clone(), scales.scales[i],num, self.wavelet.imag_freq, self.wavelet.double_sided);

                planner.plan_fft_forward(buffer.len()).process(& mut buffer);
                if self.use_normalization
                {
                    let batchsize = unsafe { std::intrinsics::ceilf64(newsize as f64 / self.threads as f64) as usize };

                    for m in 0 .. self.threads
                    {
                        let start = batchsize * m;
                        let end = std::cmp::min(newsize, batchsize * ( m + 1));

                        for n in start .. end { buffer[n] = buffer[n] / newsize as f64; }
                    }
                }
            };
        }
        else if cfg!(target_feature = "neon")
        {
            let mut planner = rustfft::FftPlannerNeon::new().unwrap();

            // //Perform forward FFT on input signal
            planner.plan_fft_forward(buffer.len()).process(& mut buffer);

            //Generate mother wavelet function
            self.wavelet.generate(newsize);
            for i in 1 .. newsize >> 1 { buffer[newsize - i] = buffer[i]; }

            for i in 0 .. scales.num_scales
            {
                //FFT-base convolution in the frequency domain
                self.daughter_wavelet_multiplication(& mut buffer, self.wavelet.mother.clone(), scales.scales[i],num, self.wavelet.imag_freq, self.wavelet.double_sided);

                planner.plan_fft_forward(buffer.len()).process(& mut buffer);
                if self.use_normalization
                {
                    let batchsize = unsafe { std::intrinsics::ceilf64(newsize as f64 / self.threads as f64) as usize };

                    for m in 0 .. self.threads
                    {
                        let start = batchsize * m;
                        let end = std::cmp::min(newsize, batchsize * ( m + 1));

                        for n in start .. end { buffer[n] = buffer[n] / newsize as f64; }
                    }
                }
            };
        }
        else
        {
            let mut planner = rustfft::FftPlannerScalar::new();
            // //Perform forward FFT on input signal
            planner.plan_fft_forward(buffer.len()).process(& mut buffer);

            //Generate mother wavelet function
            self.wavelet.generate(newsize);
            for i in 0 .. newsize >> 1 { buffer[newsize - i] = buffer[i]; }

            for i in 0 .. scales.num_scales
            {
                //FFT-base convolution in the frequency domain
                self.daughter_wavelet_multiplication(& mut buffer, self.wavelet.mother.clone(), scales.scales[i],num, self.wavelet.imag_freq, self.wavelet.double_sided);

                planner.plan_fft_forward(buffer.len()).process(& mut buffer);
                if self.use_normalization
                {
                    let batchsize = unsafe { std::intrinsics::ceilf64(newsize as f64 / self.threads as f64) as usize };

                    for m in 0 .. self.threads
                    {
                        let start = batchsize * m;
                        let end = std::cmp::min(newsize, batchsize * ( m + 1));

                        for n in start .. end { buffer[n] = buffer[n] / newsize as f64; }
                    }
                }
            };
        }
        return buffer;
    }
    fn daughter_wavelet_multiplication(& self, buffer : & mut Vec<rustfft::num_complex::Complex<f64>>, mother : Vec<f64>, scale : f64, i_size : usize, imaginary : bool, doublesided : bool)
    {
        let endpoint = std::cmp::min((i_size as f64 / 2.0) as usize, (i_size as f64 * 2.0 / scale) as usize);
        let step = scale / 2.0;

        let athreads = std::cmp::min(self.threads, std::cmp::max(1, endpoint / 16));
        let batchsize = endpoint / athreads;
        let maximum = i_size as f64 - 1.0;
        let s1 = i_size - 1;

        for m in 0 .. athreads
        {
            let start = batchsize * m;
            let end = batchsize * (m + 1);

            for n in start .. end
            {
                let tmp = std::cmp::min(maximum as usize, step as usize * n);

                match doublesided
                {
                    true =>
                        {
                            buffer[s1-n].re = match imaginary { true => { buffer[s1-n].re * mother[tmp as usize] } false => { buffer[s1-n].re * mother[tmp as usize] * - 1.0 } };
                            buffer[s1-n].im = buffer[s1-n].im * mother[tmp as usize];
                        }
                    false => { buffer[n].re = buffer[n].re * mother[tmp as usize]; buffer[n].im = buffer[n].im * mother[tmp as usize]; }
                }
            }
        }
    }
}

fn find2power(n : usize) -> usize
{
    let mut m = 0;
    let mut m2 = 1 << m; /* 2 to the power of m */
    while (m2 as isize - n as isize) < 0
    {
        m = m + 1;
        m2 <<= 1; /* m2 = m2*2 */
    }
    return m;
}