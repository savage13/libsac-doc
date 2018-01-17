====================
libsac Documentation
====================

The libsac library contains C / Fortran callable subroutines that utilize the core
funtionality of the SAC (Seismic Analysis Code) program.  The internal routines here
are wrapped in an interface that should be more streamlined to use than previous
versions. Care was taken to make sure the calling routines generate the same output
as what would be gotten from SAC; most functions here include the equivalent SAC
commands.

For exhustive examples, please see doc/examples contained in the SAC distribution.

.. contents::

User Callable Functions
=======================


Instrument Removal / Deconvolution
----------------------------------

.. code-block:: c

    // Remove Evalresp Response 
    int remove_evalresp_simple(float *data, int n, float dt, double limits[4])

    // Remove Evalresp Response 
    int remove_evalresp(float *data, int n, float dt, double limits[4],
                        char *id, char *when, char *resp_file)

    // Remove SAC_PZ Polezero Response
    int remove_polezero(float *data, int n, float dt, double limits[4],
                        char *id, char *when, char *pzfile);

    // Remove SAC_PZ Polezero Response
    int remove_polezero_simple(float *data, int n, float dt, double limits[4]);

Remove an evalresp or polezero specified instrument response.  Best to use the `simple` versions unless extra information is required.  `remove_evalresp_simple` and `remove_polezero_simple` will attempt to find the appropriate response file within the same directory before removing.  `remove_evalresp` and `remove_polezero` requires the specification of a response file.

**Arguments**

- `data` - Time series data
- `n` - Length of data
- `dt` - Sampling interval (seconds)
- `limits` - Frequency range over which response is removed (Hz)
   
   - `limits[4] = { 0.002, 0.005, 12.0, 20.0 };`
   - limits[0] - Low  Frequency Edge (Response = 0)
   - limits[1] - Low  Frequency Edge (Response = 1)
   - limits[2] - High Frequency Edge (Response = 1)
   - limits[3] - High Frequency Edge (Response = 0)
- `id` - Intrument Identifier, NET.STA.LOC.CHA

   - `*.*.*.*` wildcards are valid, but unexpected results may occur
   - `*` will try to automatically determine the instrument id using the current sac file
- `when` - Reference time, YYYY,DDD,HH:MM:SS

   - `*` will try to automatically determine the reference time using the current sac file
- `resp_file` - Evalresp Response File
- `pzfile` - Polezero Response File

**Note:** Data is modified in place.

**Examples**

.. code-block:: fortran

    implicit none
    include "sacf.h"
    real*4 y(1999), b, dt
    integer n, nerr, nmax
    real*8 limits(4)

    ! Function Prototypes
    integer remove_evalresp_simple

    nmax = 1999

    ! Frequency Limits (Hz)
    limits(1) = 0.002
    limits(2) = 0.005
    limits(3) = 12.0
    limits(4) = 20.0

    ! Read in Raw Sac file
    call rsac1("raw.sac", y, n, b, dt, nmax, nerr)

    ! Remove Response using Evalresp
    if(remove_evalresp_simple(y, n, dt, limits) .ne. 0) then
       write(*,*) "Error removing instrument with evalresp"
       return
    endif

**Effective SAC Commands**

.. code-block:: shell

    SAC> read raw.sac
    SAC> transfer from evalresp to none freq 0.002 0.005 12 20


Remove Mean
-----------

.. code-block:: c

    void  remove_mean(float *data, int n)

Remove the mean of a data series.  The mean of the data series is automatically calculated and removed from the data series. 

**Arguments**

- `data` - Input data series
- `n` - length of data

**Note:** Data is modified in place.

**Examples**

.. code-block:: fortran

    implicit none

    integer,parameter :: nmax = 1776
    integer :: npts, nerr
    real*4 :: data(nmax), beg, dt

    ! Read in the data file
    call rsac1('raw.sac', data, npts, beg, dt, nmax, nerr)

    ! Remove the mean of the data in place
    call remove_mean(data, npts)

**Effective SAC Commands**

.. code-block:: shell

    SAC> read raw.sac
    SAC> rmean


Remove Trend
------------

.. code-block:: c

    void remove_trend(float *data, int n, float delta, float b)

Remove the trend of a data series

**Arguments**

- `data` - Input data series, overwritten on output
- `n` - length of data
- `delta` - Time sampling of the data
- `b` - Initial time value of the data series

**Note:** Data is modified in place.

This calls lifite() and rtrend() internally

Trend is removed as

.. code-block:: c

    y[i] = y[i] - yint - slope * (b + delta * i)

where y is the data

**Examples**

.. code-block:: c

    #define NMAX 1969

    float y[NMAX], b, dt;
    int nmax = NMAX;
    int n, nerr;

    // Read in the data file
    rsac1("raw.sac", y, &n, &b, &dt, &nmax, &nerr, -1);

    // Remove the trend of the data in place
    remove_trend(y, n, dt, b);

**Effective SAC Commands**

.. code-block:: shell

    SAC> read raw.sac
    SAC> rtrend


Remove Trend - Unevenly sampled data
------------------------------------

.. code-block:: c

    void rtrend2(float *data, int n, float yint, float slope, float *t)

Removing the trend from an unevenly sampled data ::

**Arguments**

- `data` - Input data series, overwritten on output
- `n` - Length of data
- `yint` - Y intercept of the trend to remove
- `slope` - Slope of the trend to remove
- `t` - time values for the unevenly sampled data

**Note:** Data is modified in place.

Trend is removed as::

    y[i] = y[i] - yint - slope * t[i]

where y is the data


Rotate
------

.. code-block:: c

    void rotate(float *si1, float *si2, int ns, double angle, 
                int lnpi, int lnpo, float *so1, float *so2)

To perform a clockwise rotation of a pair of signals.

**Arguments**

- `sil` - First input signal
- `si2` - Second input signal
- `ns` - Number of points in input signals, `sil` and `si2`
- `angle` - Angle of rotateion in degrees clockwise from `si1`
- `lnpi` - True (1) if input signals are normal polarity
- `lnpo` - True (1) if output signals are normal polarity
- `so1` - First output signals
- `so2` - Second output signals 

**Note:** Input and output signals may be the same arrays

**Examples**

Rotation of two signal in Fortran

.. code-block:: fortran

    implicit none

    integer,parameter :: nmax = 1954
    integer :: npts, nerr
    real :: signal1(nmax), signal2(nmax)
    real :: rotated_signal1(nmax), rotated_signal2(nmax)
    real :: beg, dt, baz, cmpaz1
    real*8 :: angle
    logical :: lnpi, lnpo

    integer sac_compare

    ! Read in the first signal to be rotated
    call rsac1('signal1.sac', signal1, npts, beg, dt, nmax, nerr)

    call getfhv("cmpaz", cmpaz1, nerr)
    call getfhv("baz", baz, nerr);

    ! Read in the second signal to be rotated
    call rsac1('signal2.sac', signal2, npts, beg, dt, nmax, nerr)

    ! Set up parameters for rotation
    lnpi = .true.     ! input signals have "normal" polarity
    lnpo = .true.     ! output signals have "normal polarity

    ! Compute angle to rotate to
    angle = baz + 180.0 - cmpaz1;

    call rotate(signal1, signal2, npts, angle, lnpi, lnpo,
                rotated_signal1, rotated_signal2)


**Effective SAC Commands**

.. code-block:: shell

    SAC> read signal1.sac signal2.sac
    SAC> rotate


Filtering
---------


.. code-block:: c

    void bandpass(float *data, int n, float dt, float low, float high)
    void lowpass(float *data, int n, float dt, float corner)
    void highpass(float *data, int n, float dt, float corner)

    void filter(int prototype,
                int type,
                float *data, int n, float dt,
                float low, float high, int passes, int order,
                float transition,
                float attenuation)

Filter data using a Butterworh filter with two pass, four pole filter.  Data is filtered using an Infinite Impulse Repsonse Filter.  For more detailed filter types use the generic `filter` function

**Arguments**

- `data` - Input and output data
- `n` - Length of data
- `dt` - Time sampling of the data (seconds)
- `low` - low frequency corner
- `high` - high frequency corner
- `corner` - corner of the filter for `lowpass` or `highpass`

- `passes` - Number of passes

    - 1 - forward pass only
    - 2 - forward and backward pass
- `order` - Filter Order, not to exceed 10, 4-5 should be sufficient
- `transition` - Transition Bandwidth, only used in Chebyshev Type I and II Filters
- `attenuation` - Attenuation factor, amplitude reached at stopband edge, only used in Chebyshev Type I and II Filters
- `prototype` - Filter Prototype

   - 0 - Butterworth filter
   - 1 - Bessel filter
   - 2 - Chebyshev Type I filter
   - 3 - Chebyshev Type II filter
- `type` - Filter Type

   - 0 - Bandpass
   - 1 - Highpass
   - 2 - Lowpass
   - 3 - Bandreject

**Examples**

Bandpass filter in C

.. code-block:: c

    #define NMAX 2015
    float y[NMAX], b, dt;
    int n, nerr, nmax = NMAX;

    // Read in the data file
    rsac1("raw.sac", y, &n, &b, &dt, &nmax, &nerr, -1);

    // bandpass filter from 0.10 Hz to 1.00 Hz
    bandpass(y, n, dt, 0.10, 1.00);

Highpass filter in Fortran

.. code-block:: fortran

    implicit none
    integer nmax, n, nerr, sac_compare
    real*4 :: y(2012), b, dt
    nmax = 2012

    ! Read in the data file
    call rsac1("raw.sac", y, n, b, dt, nmax, nerr)

    ! highpass filter at 10.0 Hz
    call highpass(y, n, dt, 10.0)

**Effective SAC Commands**

.. code-block:: shell

    SAC> read raw.sac
    SAC> bp co 0.10 1.0 p 2 n 4

    SAC> read raw.sac
    SAC> hp co 10.0 p 2 n 4


Cross Correlation
-----------------

.. code-block:: c

    void correlate(float *f, int nf, float *g, int ng, float *c, int nc)

Compute the cross-correlation of two signals

**Arguments**

- `f` - First time series
- `nf` - Length of first time series
- `g` - Second time series
- `ng` - Length of second time series
- `c` - Cross correlation time series
- `nc` - Size of c, must be at least (nf + ng - 1)

**Return:** Cross correlation function, length: nf + ng - 1

If the signals are not the same length, then find the longest
signal, make both signals that length by filling the remainder
with zeros (pad at the end) and then run them through crscor

**Examples**

.. code-block: c

    implicit none
    character(len=*) filea, fileb
    real*4 :: amp_max, t_max
    real*4 :: a(1976), b(1976), ba,bb,dt,bc
    real*4 :: c(10000)
    integer :: nmax, nerr,na, nb, nc, i

    ! Function Prototypes
    real*4 correlate_time_begin
    integer correlate_max

    nmax = 1976

    c(:) = 0.0

    ! Read in files to correlate
    call rsac1("file1.sac", a, na, ba, dt, nmax, nerr)
    call rsac1("file2.saC", b, nb, bb, dt, nmax, nerr)

    ! Compute length of correlation
    nc = na + nb - 1

    ! Correlate
    call correlate(a, na, b, nb, c, nc)

    ! Compute begin time of correlation
    bc = correlate_time_begin(dt, na, nb, ba, bb)

    ! Compute maximum value of correlation
    i = correlate_max(c, nc)

**Effective SAC Commands**

.. code-block:: shell

    SAC> read file1.sac file2.sac
    SAC> correlate


Cross Correlation Extras
------------------------

.. code-block:: c

    int correlate_max(float *c, int nc)

Find the maximum of a correlation

**Arguments**

- `c` - float array (returned from correlate function)
- `nc` - length of c

**Return:** Index of maximum value in array

---------------------------------------------------

.. code-block:: c

    float correlate_time(float dt, float b, int i)

Compute the time of a data point given dt and begin time

**Arguments**

- `dt` - Time sampling
- `b` - Begin time
- `i` - data sample

**Return:** time value (b + i * dt)

---------------------------------------------------

.. code-block:: c

    float * correlate_time_array(float dt, float b, int n)

Compute a time array given dt and begin time

**Arguments**

- `dt` - Time sampling
- `b` - Begin time
- `n` - Length of data array

**Return:** time array

***************************************************

.. code-block:: c

    float correlate_time_begin(float dt, float n1, float _n2, float b1, float b2)

Compute begin time from a corealtion of two time series

**Arguments**

- `dt` - Time sampling
- `n1` - Length of first time series
- `n2` - Length of second time series (unused)
- `b1` - Begin time of first time series
- `b2` - Begin time of second time series

**Return:** `-dt * (n1 - 1) + (b2 - b1)`

This accounts for the possible differences in begin times of two time series


Envelope
--------

.. code-block:: c

      void envelope(int n, float *in, float *out)

Compute the envelope of a time series using the Hilbert transform

**Arguments**

- `n` - Length of input and output time series
- `in` - Input time series
- `out` - Output time series with envelope applied

The envelope is applied as such where the H(x) is the Hilbert transform::

      out = sqrt( H( in(t) )^2 + in(t)^2 )

**Examples**

.. code-block:: c

    #define NMAX 1929
    int nlen, nerr, nmax;
    float yarray[NMAX], yenv[NMAX];
    float beg, delta;

    nmax = NMAX;

    // Read in data file
    rsac1("raw.sac", yarray, &nlen, &beg, &delta, &nmax, &nerr, SAC_STRING_LENGTH);

    // Calculate Envelope of data
    envelope(nlen, yarray, yenv);

**Effective SAC Commands**

.. code-block:: shell

    SAC> read raw.sac
    SAC> envelope


Differentiate 2-pt
------------------


.. code-block:: c

     void dif2(float *array, int n, double step, float *output)

Differentiate a data set using a two point differentiation

**Arguments**

- `array` - Input data to differentiate
- `n` - length of ararry
- `step` - Time sampling of input data 
- `output` - Output differentiated data, length n-1

This is the default scheme in the SAC program.

The output array will be 1 data point less than the input array.

Since this is not a centered differeniation, there is an implied shift
in the independent variable by half the step size::

    b_new = b_old + 0.5 * step

Differntiation is performed as::

    out[i] = (1/step) * (in[i+1] - in[i])

Differentiate 3-pt
------------------

.. code-block:: c

     void dif3(float *array, int n, double step, float *output)

Perform "three-point" (centered two-point) differntiation

**Arguments**

- `array` - Input data to differentiate
- `n` - length of ararry
- `step` - Time sampling of input data 
- `output` - Output differentiated data, length n-2

The output array will be 2 data points less than the input array.

There is an implied shift in the time variable by a full step size::

    b_new = b_old + step

Differntiation is performed as::

    out[i] = 1/(2 * step) * (in[i+1] - in[i-1])

Differentiate 5-pt
------------------

.. code-block:: c

     void dif3(float *array, int n, double step, float *output)

Perform "five-point" (centered four-point) differntiation

**Arguments**

- `array` - Input data to differentiate
- `n` - length of ararry
- `step` - Time sampling of input data 
- `output` - Output differentiated data, length n-2

The output array will be 2 data points less than the input array.

There is an implied shift in the time variable by a full step size::

    b_new = b_old + step

Differntiation is performed in the interior as ::

    out[i] = 2/(3  * step) * (in[i+1] - in[i-1]) -
             1/(12 * step) * (in[i+2] - in[i-2])

Differentiation at the end points is::

    out[1]   = (in[3] -   in[i]) / (2 * step)
    out[n-2] = (in[n] - in[n-2]) / (2 * step)

Integerate - Trapezodial
------------------------


.. code-block:: c

    void int_trap(float *y, int n, double delta)

Integrate a data series using the trapezodial method

**Arguments**

- `y` - Input data series, overwritten on output
- `n` - length of y
- `delta` - time sampling of the data series

Integration is performed as::

    out[i] = out[i-1] + (delta/2) * (in[i] + in[i+1])

where the initial out value is 0.0.

The number of points on output should be reduced by 1 ::

     len(out) = len(in) - 1

and the beging value is shifted by -0.5 delta::

     b_out = b_in - 0.5 * delta

Integerate - Rectangular
------------------------


.. code-block:: c

    void int_rect(float *y, int n, double delta)

Integrate a data series using the rectangular method

**Arguments**

- `y` - Input data series, overwritten on output
- `n` - length of y
- `delta` - time sampling of the data series

Integration is performed as::

    out[i] = (delta * in[i]) + out[i-1]

and the initial value is::

    out[0] = in[0] * delta

The number of points and b value are left unchanged here::

    len(out) = len(in)
    b_out    = b_in


RMS
---

.. code-block:: c

      double rms(float *x, int nsamps)

Compute Root Mean Square value of an array

**Arguments**

- `x` - input array to find the rms value
- `nsamps` - length of input array

RMS value is computed as such::

      rms = sqrt( sum (x_i^2) )

**Note** This routine does not divide by the number of points

Window
------

.. code-block:: c

      void window(float x, int n, char *ftype, int fsamp,
                  int wlen, float *y, char *err, int err_s)

Window a sequence

**Arguments**

- `x` - input array
- `n` - length of input array
- `ftype` - type of window to apply

    - `HAM` Hamming window
    - `HAN` Hanning window
    - `R` - Rectangular window
    - `C` - 10% cosine taper window
    - `T` - Triangular window
- `fsamp` - index of first sample of the window
- `wlen` - window length in samples
- `y` - output, windowed sample
- `err` - error condition message
- `err_s` - length of message err

For all window types::

      y_i = 0    if i < fsamp or i > fsamp + wlen - 1

- Rectangular window::

      y = x_i

- Triangular window::

      y = x_i * (1 - abs(i - center) / extent)
      center = (wlen - 1)/2 + fsamp
      extent = center - fsamp

-  Hamming window::

      y_i = x_i * (f0 + f1 * cos(omega * (i-fsamp) - pi)
      f0 = 0.54
      f1 = 0.46

- Hanning window::

      y_i = x_i * (f0 + f1 * cos(omega * (i-fsamp) - pi)
      f0 = 0.5
      f1 = 0.5

-  Cosine window::

      if i < fsamp + wlen/10
         y_i = x_i * 0.5 * (1.0 - cos( pi * (i-fsamp) )
      if i > fsamp - wlen/10
         y_i = x_i * 0.5 * (1.0 - cos( pi * (fsamp + wlen - 1 - i) )
      else
         y_i = x_i

