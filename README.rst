fbm
===

* Exact methods for simulating fractional Brownian motion (fBm) or fractional
  Gaussian noise (fGn) in python.
* *Approximate* simulation of multifractional Brownian motion (mBm) or
  multifractional Gaussian noise (mGn).

Installation
------------

The fbm package is available on PyPI and can be installed via pip:

.. code-block::

    pip install fbm

fractional Brownian motion
--------------------------

Fractional Brownian motion can be generated via either Hosking's method, the
Cholesky method, or the Davies-Harte method. All three methods are
theoretically exact in generating a discretely sampled fBm/fGn.

Usage:

.. code-block:: python

    from fbm import FBM


    f = FBM(n=1024, hurst=0.75, length=1, method='daviesharte')
    # or
    f = FBM(1024, 0.75)

    # Generate a fBm realization
    fbm_sample = f.fbm()

    # Generate a fGn realization
    fgn_sample = f.fgn()

    # Get the times associated with the fBm
    t_values = f.times()

where ``n`` is the number of equispaced increments desired for a fBm with Hurst
parameter ``hurst`` on the interval [0, ``length``]. Method can be
either ``'hosking'``, ``'cholesky'``, or ``'daviesharte'``. The ``fbm()``
method returns a length ``n+1`` array of discrete values for the fBm (includes
0). The ``fgn()`` method returns a length ``n`` array of fBm
increments, or fGn. The ``times()`` method returns a length ``n+1`` array of
times corresponding to the fBm realizations.

The ``n`` and ``hurst`` parameters are required. The ``length`` parameter
defaults to 1 and ``method`` defaults to ``'daviesharte'``.

For simulating multiple realizations use the FBM class provided as above. Some
intermediate values are cached for repeated simulation.

For one-off samples of fBm or fGn there are separate functions available:

.. code-block:: python

    from fbm import fbm, fgn, times


    # Generate a fBm realization
    fbm_sample = fbm(n=1024, hurst=0.75, length=1, method='daviesharte')

    # Generate a fGn realization
    fgn_sample = fgn(n=1024, hurst=0.75, length=1, method='daviesharte')

    # Get the times associated with the fBm
    t_values = times(n=1024, length=1)

For fastest performance use the Davies and Harte method. Note that the
Davies and Harte method can fail if the Hurst parameter ``hurst`` is close to
1 and there are a small amount of increments ``n``. If this occurs, a warning
is printed to the console and it will fallback to using Hosking's method to
generate the realization. See page 412 of the following paper for a more
detailed explanation:

* Wood, Andrew TA, and Grace Chan. "Simulation of stationary Gaussian processes
  in [0, 1] d." Journal of computational and graphical statistics 3, no. 4
  (1994): 409-432.


**Hosking's method:**

* Hosking, Jonathan RM. "Modeling persistence in hydrological time series
  using fractional differencing." Water resources research 20, no. 12 (1984):
  1898-1908.

**Cholesky method:**

* Asmussen, Søren. Stochastic simulation with a view towards stochastic
  processes. University of Aarhus. Centre for Mathematical Physics and
  Stochastics (MaPhySto)[MPS], 1998.

**Davies Harte method:**

* Davies, Robert B., and D. S. Harte. "Tests for Hurst effect." Biometrika 74,
  no. 1 (1987): 95-101.


multifractional Brownian motion
-------------------------------

This package supports *approximate* generation of multifractional
Brownian motion. The current method uses the Riemann–Liouville fractional
integral representation of mBm.

Usage:

.. code-block:: python

    import math
    from fbm import MBM


    # Example Hurst function with respect to time.
    def h(t):
        return 0.25 * math.sin(20*t) + 0.5

    m = MBM(n=1024, hurst=h, length=1, method='riemannliouville')
    # or
    m = MBM(1024, h)

    # Generate a mBm realization
    mbm_sample = m.mbm()

    # Generate a mGn realization
    mgn_sample = m.mgn()

    # Get the times associated with the mBm
    t_values = m.times()


The ``hurst`` argument here should be a callable that accepts one argument
and returns a float in (0, 1).

For one-off samples of mBm or mGn there are separate functions available:

.. code-block:: python

    from fbm import mbm, mgn, times


    # Define a hurst function
    def h(t):
        return 0.75 - 0.5 * t

    # Generate a mbm realization
    mbm_sample = mbm(n=1024, hurst=h, length=1, method='riemannliouville')

    # Generate a fGn realization
    mgn_sample = mgn(n=1024, hurst=h, length=1, method='riemannliouville')

    # Get the times associated with the mBm
    t_values = times(n=1024, length=1)


**Riemann-Liouville representation method:**

*Approximate* method originally proposed for fBm in

* Rambaldi, Sandro, and Ombretta Pinazza. "An accurate fractional Brownian
  motion generator." Physica A: Statistical Mechanics and its Applications 208,
  no. 1 (1994): 21-30.

Adapted to approximate mBm in

* Muniandy, S. V., and S. C. Lim. "Modeling of locally self-similar processes
  using multifractional Brownian motion of Riemann-Liouville type." Physical
  Review E 63, no. 4 (2001): 046104.
