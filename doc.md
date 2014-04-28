# libsac Documentation


- xapiir - Signal Filtering
- crscor - Compute the Cross Correlation of two signals
- envelope - Compute the envelope of a singal
- rms - Compute the RootMeanSquare of an signal

- Internal Subroutines (Not called by the user)
    - fft - Compute the Fast Fourier Transform
    - copydouble - Copy a singal from one array to another
    - copyfloat - Copy a signal from one array to another
    - window - Window a sequence
    - next2 - Find a power of 2 greater than an integer
    - firtrn - FIR Filter
    - overlap - Overlap Save Routine
    - zshft - Shfit a signal with zero-filling
    - zero - Multiply a singal by 0.0
    - apply - Apply a IIR Digital filter to a signal (in-place)
    - design - Design IIR Digital filter from Analog Prototype
    - bilin2 - Convert Analog Filter to a Digital Filter using Bilinear Transform
    - warp - Tangent frequency warping to account for Bilinear Transformation
    - cutoffs - Alter cutoffs of a Analog Filter
    - Analog Filter Transformations
        - lp
        - lptbp
        - lptbr
        - lpthp
    - Analog Filter Prototypes
        - chebparm
        - buroots
        - beroots
        - c1roots
        - c2roots

## firtrn


Calculate Hilbert or derivative of a signal with a FIR Filter

    void firtrn(char *ftype, float *x, int n, float *buffer, float *y)
    
    ftype - Desired Transform
        - HILBERT
        - DERIVATIVE of Hilbert Transform
    x - input signal
    n - length of input signal
    buffer - temp storage, min size of 4297
    y - output array
    
Hilbert transform coefficients

    C_i = (2.0/pi*(2*i-1)) * (0.54 + 0.46 *cos( 2 * pi * ( 2*i - 1 ) / 201



## overlp
Overlap - save routine

      void overlp(float *input, int npts, float *output, float *c, int nc, int nfft, float *buffer, float *cbuff)
      
      input - input signal
      npts - length of input signal
      output - filtered output, may be the same as input
      c - coefficent sequence of filter
      nc - length of coeffience sequence
      nfft - lenght of FFT in convolutions
      buffer - temp storage - must be 2 * nfft
      cbuff - temp storage - must be 2 * nfft
      
Perform convolution of signals input and c to apply filter / transform

## zshft
Shift a signal in place with zero-filling

    void zshft(float *signal, int n, in ishft) 
    
    signal - input signal to be shifted
    n - length of signal
    ishft - number of samples to shift
      - > 0 shift to the right
      - < 0 shift to the left
      

## zero
Multiply (or set) a signal to 0.0

    void zero(float *a, int n)
    
    a - signal to set to 0.0
    n - length of signal

## fft
Compute Fast Fourier transform of a sequence

    void fft(float *xreal, float *ximag, int n, int idir)
    
    xreal - Real part of signal
    ximag - Imaginary part of signal
    n - length of signal xreal and ximag
        signal must be a power of 2
    idir - Direction of tranform
       - -1 Forward
       -  1 Inverse transform (normalization performed)



## crscor
Compute the cross-correlation of two signals

    void crscor(float *data1, float *data2, int nsamps, int nwin, int wlen, char *type, float *c, int *nfft, char *err, int err_s)
    
    data1 - first signal
    data2 - second signal
    nsamps - length of longest signal
    nwin - number of windows
    wlen -number of samples in each window; max 2048
    type - type of window
       - HAM Hamming window
       - HAN Hanning window
       - R - Rectangular window
       - C - 10% cosine taper window
       - T - Triangular window
    c - output cross correlation coefficients, length 2*wlen - 1
    nfft - number of sample in correlation signal
    err - Error messages
    err_s - Length of message err

Each window is processed as follows and summed into C

    scale_i = rms(data_i)
    Data1 = FFT(window(data1) / scale)
    Data2 = FFT(window(data2) / scale)
    C = C + (conjg(Data1) x Data2)

Then transformed back to the time domain
    
    c = iFFT(C)

## window

Window a sequence 

    void window(float x, int n, char *ftype, int fsamp, int wlen, float *y, char *err, int err_s)
    
    x - input array
    n - length of input array
    ftype - type of window to apply
       - HAM Hamming window
       - HAN Hanning window
       - R - Rectangular window
       - C - 10% cosine taper window
       - T - Triangular window
    fsamp - index of first sample of the window
    wlen - window length in samples
    y - output, windowed sample
    err - error condition message
    err_s - length of message err

For all window types
    
    y_i = 0    if i < fsamp or i > fsamp + wlen - 1

Otherwise

#### Rectangular window

    y = x_i  

#### Triangular window

    y = x_i * (1 - abs(i - center) / extent)
    center = (wlen - 1)/2 + fsamp
    extent = center - fsamp
    
#### Hamming window

    y_i = x_i * (f0 + f1 * cos(omega * (i-fsamp) - pi)
    f0 = 0.54
    f1 = 0.46
    
#### Hanning window

    y_i = x_i * (f0 + f1 * cos(omega * (i-fsamp) - pi)
    f0 = 0.5
    f1 = 0.5
#### Cosine window

    if i < fsamp + wlen/10
       y_i = x_i * 0.5 * (1.0 - cos( pi * (i-fsamp) )
    if i > fsamp - wlen/10
       y_i = x_i * 0.5 * (1.0 - cos( pi * (fsamp + wlen - 1 - i) )
    else
       y_i = x_i
rms
---
Compute rms value of an array

    double rms(float *x, int nsamps)
    
    x - input array to find the rms value 
    nsamps - length of input array
    
RMS value is computed as such

rms = sum (x_i^2)

copydouble
----------
Copy a double precision floating point array

    void copydouble(double *source, int length, double *sink) 
    
    source - input array
    length - length of input and output array
    sink - output array
    
copyfloat
---------
Copy a single precision floating point array

    void copyfloat(float *src, float *dest, int n)
    
    src - input array
    dest - output array
    n - length of input and output array
    
next2
-----
Find the next power of 2 greater than number

    int next2(int num)
    
    num - Number to find the next power of 2 greater than

envelope
--------
Envelope of a time series using the Hilbert transform

    void envelope(int n, float *in, float *out)

    n - Length of input and output time series
    in - Input time series
    out - Output time series with envelope applied

The envelope is applied as such where the H(x) is the Hilbert transform

  out = sqrt( H( in(t) )^2 + in(t)^2 )

xapiir
------

apply
-----

design
------

lp
--

lptbp
-----

lptbr
-----

lpthp
-----

warp
----

cutoffs
-------

chebparm
--------

buroots
-------

beroots
-------

c1roots
-------

c2roots
-------

bilin2
------













