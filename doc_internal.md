Internal Subroutines
====================

**These functions/subroutines should probably not be called by the user. Use at your own risk.**


Instrument Deconvolution (Internal)
-----------------------------------

Apply a instrument change, either removal or additon, in the frequency domain.
This is also called a deconvolution and convolution::

    void ztransfer(float *dat, int npts, double delta, double *sre, double *sim,
                   double *xre, double *xim, int nfreq, int nfft, double delfrq,
                   double *F)

- dat - Input data
- npts - length of input data
- delta - time sampling of data
- sre - Instrument to be Removed, Real component, "FROM transfer function"
        The real component should be inverted (1/z) on input
- sim - Instrument to be Removed, Imaginary component, "FROM transfer function"
        The real component should be inverted (1/z) on input
- xre - Instrument to be Added, Real component, "TO transfer function"
- xim - Instrument to be Added, Imaginary component, "TO transfer function"
- nfreq - nfft / 2 + 1
- nfft - Length of sre, sim, xre, xim (Frequency domain)
- delfrq - Sample interval of points in the frequency domain 
           1.0 / (nfft * delta)
- F - Frequencies over which to apply the action
    - F[0] - Low  Frequency at which action is not performed, but set to 0.0
    - F[1] - Low  Frequency at which action is performed
    - F[2] - High Frequency at which action is performed
    - F[3] - High Frequency at which action is not performed, but set to 0.0

Note: Wow, this is super complicated with lots of pitfalls.
It should be streamlined for typical use cases, like removal and convolution.


Transfer - Read PoleZero files
------------------------------

Filter data using a Infinite Impulse Response Filter.

.. code-block:: c

      void xapiir(float *data, int nsamps, char *aproto,
                  double trbndw, double a,
                  int iord, char *type,
                  double flo, double fhi, double ts, int passes)


- `data` - input time series to be filtered, overwritten on output
- `nsamps` - length of data
- `aproto` - Filter Prototype

   - "BU" - Butterworth filter
   - "BU" - Bessel filter
   - "BU" - Chebyshev Type I filter
   - "BU" - Chebyshev Type II filter
- `trbndw` - Transition Bandwidth, only used in Chebyshev Type I and II Filters
- `a` - Attenuation factor, amplitude reached at stopband edge, only used in Chebyshev Type I and II Filters
- `iord` - Filter Order, not to exceed 10, 4-5 should be sufficient
- `type` - Filter Type

   - "LP" - Lowpass
   - "HP" - Highpass
   - "BP" - Bandpass
   - "BR" - Bandreject
- `flo` - Low frequency cutoff, not used in Lowpass filters
- `fhi` - High frequneyc cutoff, not used in Highpass filters
- `ts` - Sampling rate (seconds)
- `passes` - Number of filter passes

   - 1 - Forward filter
   - 2 - Forward and Reverse filtering, i.e. zero-phase filtering


Cross Correlation (Internal)
----------------------------

Compute the cross-correlation of two signals.  This is called by correlate::

       void
       crscor(float *data1, float *data2, int nsamps, int nwin, int wlen, 
              char *type, float *c, int *nfft, char *err, int err_s)

- data1 - first signal
- data2 - second signal
- nsamps - length of longest signal
- nwin - number of windows
- wlen -number of samples in each window; max 2048
- type - type of window
   - HAM Hamming window
   - HAN Hanning window
   - R - Rectangular window
   - C - 10% cosine taper window
   - T - Triangular window
- c - output cross correlation coefficients, length 2*wlen - 1
- nfft - number of sample in correlation signal
- err - Error messages
- err_s - Length of message err

Each window is processed as follows and summed into C::

       scale_i = rms(data_i)
       Data1 = FFT(window(data1) / scale)
       Data2 = FFT(window(data2) / scale)
       C = C + (conjg(Data1) x Data2)

Then transformed back to the time domain::

       c = iFFT(C)


Find Slope of Data
------------------

Fit a line to the a data series using a linear least squares approach::

    void lifite(double x1, double dx, float *y, int n, float *a, float *b, 
                float *siga, float *sigb, float *sig, float *cc)

- x1 - Begining value, b-value typically
- dx - Time sampling of data
- y - Input data series (Amplitude)
- n - length of y
- a - Output slope of linear fit
- b - Output trend of linear fit
- siga - Standard deviation of slope
- sigb - Standard deviation of trend
- sig - Standard deviation of data
- cc - Correlation coefficent between data and linear fit


Find Slope of Data - Unevenly spaced data
-----------------------------------------

Fit a line to an unevely spaced data using a linear least squares approach ::

     void lifitu(float *x, float *y, int n, float *a, float *b,
                 float *siga, float *sigb, float *sig, float *cc)

- x - Input time values of the data
- y - Input data values (amplitude)
- n - length of x and y
- a - Output slope of the linear fit
- b - Output trend of the linear fit
- siga - Standard deviation of the slope
- sigb - Standard deviation of the trend
- sig - Standard deviation of data
- cc - Correlation coefficent between data and linear fit


FIR Filter
----------

Calculate Hilbert or derivative of a signal with a FIR Filter::

      void firtrn(char *ftype, float *x, int n, float *buffer, float *y)

- ftype - Desired Transform
     - HILBERT
     - DERIVATIVE of Hilbert Transform
- x - input signal
- n - length of input signal
- buffer - temp storage, min size of 4297
- y - output array

Hilbert transform coefficients::

      C_i = (2.0/pi*(2*i-1)) * (0.54 + 0.46 *cos( 2 * pi * ( 2*i - 1 ) / 201

Overlap
-------

Overlap - save routine::

      void overlp(float *input, int npts, float *output, float *c, int nc, int nfft, float *buffer, float *cbuff)

- input - input signal
- npts - length of input signal
- output - filtered output, may be the same as input
- c - coefficent sequence of filter
- nc - length of coeffience sequence
- nfft - lenght of FFT in convolutions
- buffer - temp storage - must be 2 * nfft
- cbuff - temp storage - must be 2 * nfft

Shift and Filling
-----------------

Shift a signal in place with zero-filling::

      void zshft(float *signal, int n, in ishft)

- signal - input signal to be shifted
- n - length of signal
- ishft - number of samples to shift
   - > 0 shift to the right
   - < 0 shift to the left

zero
----

Multiply (or set) a signal to 0.0::

      void zero(float *a, int n)

- a - signal to set to 0.0
- n - length of signal

FFT
---

Compute Fast Fourier transform of a sequence::

       void fft(float *xreal, float *ximag, int n, int idir)

- xreal - Real part of signal
- ximag - Imaginary part of signal
- n - length of signal xreal and ximag
        signal must be a power of 2
- idir - Direction of tranform
    - -1 Forward
    -  1 Inverse transform (normalization performed)

Copy Double Array
-----------------

Copy a double precision floating point array::

      void copydouble(double *source, int length, double *sink)

- source - input array
- length - length of input and output array
- sink - output array

Copy Float Array
----------------

Copy a single precision floating point array::

      void copyfloat(float *src, float *dest, int n)

- src - input array
- dest - output array
- n - length of input and output array

Next Power of Two
-----------------

Find the next power of 2 greater than number::

      int next2(int num)

- num - Number to find the next power of 2 greater than

IIR - Apply Filter
------------------

This is not what you want, please see xapiir()

Apply an IIR (Infinite Impulse Response) Filter to a data sequence.  The filter is assumed to be stored as second
order sections.  The filtering is done in place.  Zero-phase (forward plus reverse filtering) is an option::

   void apply(float *data, int nsamps, int zp, float *sn, float *sd, int nsects)

- data - Array containing data on input and filtered data on output
- nsamps - Length of array  data
- zp - If Zero Phase filtering is requested
   - TRUE Zero phase Two pass filtering (forward + reverse filters)
   - FALSE Single Pass filtering
- sn - Numerator polynomials for 2nd Order Sections
- sd - Denominator polynomials for 2nd Order Sections
- nsects - Number of 2nd Order Sections


The filter is applied as such::

       y_n = b[0] * x[n]   + b[1] * x[n-1] + b[2] * x[n-2] -
                             a[1] * y[n-1] + b[2] * y[n-2]
       where
         - N = nsamps
         - b = sn
         - a = sd


IIR - Design a Filter
---------------------

Design IIR Digital Filters from Analog Prototypes::

    void design(int iord, char *type, char *aproto, double a, double trbndw, double fl,
                double fh, double ts, float *sn, float *sd, int *nsects)

- iord - Filter Order, Maximum of 10
- type - Filter Type
     - 'LP'  Lowpass
     - 'HP'  Highpass
     - 'BP'  Bandpass
     - 'BR'  Bandreject
- aproto - Analog Prototype
     - 'BU'  Butterworth
     - 'BE'  Bessel
     - 'C1'  Cheyshev Type I
     - 'C2'  Cheyshev Type II
- a -  Chebyshev stopband Attenuation Factor
- trbndw -  Chebyshev transition bandwidth, fraction of lowpass prototype passband width
- fl - Low Frequency cutoff
- fh - High Frequency cutoff
- ts - Sampling Rate / Delta
- sn - Array containing numerator coefficients of 2nd Order Sections, Packed Head to Tail
- sd - Array containing denominator coefficients of 2nd Order Sections, Packed Head to Tail
- nsects -  Length of arrays  sn and  sd



IIR - Convert lp to lp
----------------------

Subroutine to generate second order section parameterization from an pole-zero description for lowpass filters::

    void lp(complexfp, complexfz, char *rtype, int rtype_s, double dcvalue,
            int nsects, float *sn, float *sd)


- p - Array of Poles
- z - Array of Zeros
- rtype - Character array containing root type information
     - "SP"  Single real pole
     - "CP"  Complex conjugate pole pair
     - "CPZ" Complex conjugate pole and zero pairs
- rtype_s - Length of  rtype
- dcvalue - Zero-frequency value of prototype filter
- nsects - Number of second-order sections
- sn - Output Numerator polynomials for second order sections.
- sd - Output Denominator polynomials for second order sections.



IIR - Convert lp to bp
----------------------

Subroutine to convert an prototype lowpass filter to a bandpass filter via the analog polynomial transformation.
The lowpass filter is described in terms of its poles and zeros (as input to this routine).  The output consists
of the parameters for second order sections.::


    void lptbp(complexfp, complexfz, char *rtype, int rtype_s, double dcvalue,
               int *nsects, double fl, double fh, float *sn, float *sd)

- p - Array of Poles
- z - Array of Zeros
- rtype - Character array containing root type information
     - "SP"  Single real pole
     - "CP"  Complex conjugate pole pair
     - "CPZ" Complex conjugate pole and zero pairs
- rtype_s - Length of rtype
- dcvalue - Zero-frequency value of prototype filter
- nsects - Number of second-order sections.
    On output this subroutine doubles the number of
    sections.
- fl - Low Frequency cutoff
- fh - High Frequency cutoff
- sn - Output Numerator polynomials for second order sections.
- sd - Output Denominator polynomials for second order sections.



IIR - Convert lp to br
----------------------

 Subroutine to convert a lowpass filter to a band reject filter via an analog polynomial transformation.
 The lowpass filter is described in terms of its poles and zeros (as input to this routine).  The output
 consists of the parameters for second order sections.::

   void  lptbr(complexfp, complexfz, char *rtype, int rtype_s, double dcvalue,
              int *nsects, double fl, double fh, float *sn, float *sd)


- p - Array of Poles
- z - Array of Zeros
- rtype - Character array containing root type information
     - "SP"  Single real pole
     - "CP"  Complex conjugate pole pair
     - "CPZ" Complex conjugate pole and zero pairs
- rtype_s - Length of rtype
- dcvalue - Zero-frequency value of prototype filter
- nsects - Number of second-order sections.
    On output this subroutine doubles the number of
    sections.
- fl - Low Frequency cutoff
- fh - High Frequency cutoff
- sn - Output Numerator polynomials for second order sections.
- sd - Output Denominator polynomials for second order sections.


IIR - Convert lp to hp
----------------------

Subroutine to convert a lowpass filter to a highpass filter via an analog polynomial transformation.  The lowpass filter is
described in terms of its poles and zeroes (as input to this routine).  The output consists of the parameters for
second order sections::

    void lpthp(complexfp, complexfz, char *rtype, int rtype_s, double dcvalue,
               int nsects, float *sn, float *sd)

- p - Array of Poles
- z -  Array of Zeros
- rtype - Character array containing root type information
     - "SP"  Single real pole
     - "CP"  Complex conjugate pole pair
     - "CPZ" Complex conjugate pole and zero pairs
- rtype_s - Length of rtype
- dcvalue -  Zero-frequency value of prototype filter
- nsects - Number of second-order sections.
- sn - Output Numerator polynomials for second order sections.
- sd - Output Denominator polynomials for second order sections.




IIR - Frequency warping
-----------------------

Applies tangent frequency warping to compensate for bilinear analog -> digital transformation::

  double warp(double f, double ts)

- f - Original Design Frequency Specification (Hz)
- ts -  Sampling Internal (seconds)

IIR - Alter Filter Cutoffs
--------------------------

Alter the cutoff of a filter.  Assumed that the filter is structured as 2nd order sections.
Changes the cutoffs of a normalized lowpass or highpass filters through a simple polynomial transformation::

    void cutoffs(float *sn, float *sd, int nsects, double f)

- sn - Numerator polynomials for 2nd order sections
- sd - Denominator polynomials for 2nd order sections
- nsects - Number of 2nd order sections
- f - New cutoff frequency


IIR - Calculates Chebparm
-------------------------

Calculate Chebyshev Type I and II Design Parameters::

   void chebparm(double a, double trbndw, int iord, float *eps, float *ripple)

- a - Desired Stopband Attenuation, i.e. max stopband amplitude is 1/ATTEN
- trbndw - Transition bandwidth between stop and passband as a fraction of the passband width
- iord - Filter Order (number of Poles)
- eps - Output Chebyshev passband parameter
- ripple - Passband ripple

Chebyshev parameters are calculated as::

  omega = 1.0 + trbndw
  alpha = (omega + sqrt(omega^2 - 1.0) )^iord
  g = alpha^2 + 1 / 2alpha
  eps = sqrt(a^2 - 1.0) / g
  ripple = 1 / sqrt(1.0 + eps^2)

IIR - Calculate Butterworth Poles
----------------------------------
Compute the Butterworth Poles for a Normalized Low Pass (LP) Filter::

    void bu roots(complexfp, char *rtype, int rtype_s, float *dcvalue, int *nsects, int iord)

- p - Complex Array containing Poles. Contains only one of from each
    - Complex Conjugate Pole-Zero Pair
    - Complex Conjugate Pole Pair
    - Single Real Pole
- rtype - Character Array indicating 2nd Order Section Types
    - 'CPZ' Complex Conjugate Pole-Zero Pair
    - 'CP'  Complex Conjugate Pole Pair
    - 'SP'  Single Real Pole
- rtype_s - Length of string rtype
- dcvalue - Magnitude of the filter at Zero Frequency
- nsects - Number of 2nd Order Sections
- iord - Desired Filter Order, Must be between 1 and 8


IIR - Calculate Bessel Poles
----------------------------

Compute Bessel Poles For a Normalized Low Pass (LP) Filter::

    void beroots(complexfp, char *rtype, int rtype_s, float *dcvalue, int *nsects,
                 int iord)


- p - Complex Array containing Poles. Contains only one of from each
    - Complex Conjugate Pole-Zero Pair
    - Complex Conjugate Pole Pair
    - Single Real Pole
- rtype - Character Array indicating 2nd Order Section Types
    - 'CPZ' Complex Conjugate Pole-Zero Pair
    - 'CP'  Complex Conjugate Pole Pair
    - 'SP'  Single Real Pole
- rtype_s - Length of string rtype
- dcvalue - Magnitude of the filter at Zero Frequency
- nsects - Number of 2nd Order Sections
- iord - Desired Filter Order, Must be between 1 and 8



IIR - Calculate Chebyshev Type I Roots
--------------------------------------

Compute Chebyshev Type I Poles for a Normalized Low Pass (LP) Filter::

   void c1roots(complexfp, char *rtype, int rtype_s, float *dcvalue, int *nsects,
               int iord, double eps)


- p - Complex Array containing Poles. Contains only one of from each
    - Complex Conjugate Pole-Zero Pair
    - Complex Conjugate Pole Pair
    - Single Real Pole
- rtype - Character Array indicating 2nd Order Section Types
    - 'CPZ' Complex Conjugate Pole-Zero Pair
    - 'CP'  Complex Conjugate Pole Pair
    - 'SP'  Single Real Pole
- rtype_s - Length of string rtype
- dcvalue - Magnitude of the filter at Zero Frequency
- nsects - Number of 2nd Order Sections
- iord -  Desired Filter Order, Must be between 1 and 8
- eps - Output Chebyshev Parameter Related to Passband Ripple

IIR - Calculate Chebyshev Type II Roots
---------------------------------------

Compute root for normailzed Low Pass Chebyshev Type II Filter::

    void c2roots(complexfp, complexfz, char *rtype, int rtype_s, float *dcvalue,
                int *nsects, int iord, double a, double omegar)

- p - Complex Array containing Poles. Contains only one of from each
     - Complex Conjugate Pole-Zero Pair
     - Complex Conjugate Pole Pair
     - Single Real Pole
- z - Complex Array containing Zeros Contains only one of from each
     - Complex Conjugate Pole-Zero Pair
     - Complex Conjugate Pole Pair
     - Single Real Pole
- rtype - Character Array indicating 2nd Order Section Types
     - 'CPZ' Complex Conjugate Pole-Zero Pair
     - 'CP'  Complex Conjugate Pole Pair
     - 'SP'  Single Real Pole
- rtype_s - Length of string rtype
- dcvalue - Magnitude of filter at Zero Frequency
- nsects - Number of 2nd order sections
- iord - Input Desired Filter Order
- a - Input Stopband attenuation factor
- omegar - Input Cutoff frequency of stopband passband cutoff is 1.0 Hz


IIR - Bilinear Transform
------------------------

Transform an analog filter to a digital filter via the bilinear transformation.
Assumes both filters are stored as 2nd Order sections and the transform is done in place::

   void bilin2(float *sn, float *sd, int nsects)

- sn -  Array containing numerator polynomial coefficeients. Packed head to tail and using groups of 3.  Length is 3nsects
- sd -  Array containing demoninator polynomial coefficeients. Packed head to tail and using groups of 3.  Length is 3nsects
- nsects -  Number of 2nd order sections.

The bilinear transform for each 2nd order section is apply as such ::

   scale = a0 + a1 + a2

   a_0 = 1.0
   a_1 = 2.0(a_0 - a_2) / scale
   a_2 = (a_2 - a_1 + a_0) / scale

   b_0 = (b_2 + b_1 + b_0) / scale
   b_1 = 2.0(b_0 - b_2) / scale
   b_2 = (b_2 - b_1 + b_0) / scale
