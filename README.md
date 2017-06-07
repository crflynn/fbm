# fbm
Exact methods for simulating fractional Brownian motion (fBm) or fractional Gaussian noise (fGn) in python.

The three methods are Hosking's method, the Cholesky method, and the Davies Harte method. All three methods are exact in generating a discretely sampled fBm/fGn.

Usage:

```python
from fbm import FBM


f = FBM(n=16, H=0.75, L=1, method='daviesharte')

# Generate a fBm realization
fbm_sample = f.sample()

# Generate a fGn realization
fgn_sample = f.sample_noise()

# Get the times associated with the fBm
times = f.times()
```

where `n` is the number of equispaced increments desired for a fBm with Hurst parameter `H` on the interval [0, `L`]. Method can be either `'hosking'`,`'cholesky'`, or `'daviesharte'`. The `sample` method returns a length `n+1` array of discrete values for the fBm (includes 0). The `sample_noise` method returns a length `n` array of fBm increments, or fGn. The `times` method returns a length `n+1` array of times corresponding to the fBm realizations.

For simulating multiple realizations use the FBM class provided as above. For one-off samples of fBm or fGn there are also functions available which handle the FBM object themselves:

```python
from fbm import fbm, fgn, times


# Generate a fBm realization
fbm_sample = fbm(n=16, H=0.75, L=1, method='daviesharte')

# Generate a fGn realization
fgn_sample = fgn(n=16, H=0.75, L=1, method='daviesharte')

# Get the times associated with the fBm
times = times(n=16, L=1)
```

For fastest performance use the Davies and Harte method. It is much faster than both other methods especially for larger increment quantities. Note that the Davies and Harte method can fail if the Hurst parameter `H` is close to 1 and there are a small amount of increments `n`. If this occurs, python will print a warning to the console and fallback to using Hosking's method to generate the realization. See page 412 of the following paper for a more detailed explanation

* Wood, Andrew TA, and Grace Chan. "Simulation of stationary Gaussian processes in [0, 1] d." Journal of computational and graphical statistics 3, no. 4 (1994): 409-432.


**Hosking's method:**

* Hosking, Jonathan RM. "Modeling persistence in hydrological time series using fractional differencing." Water resources research 20, no. 12 (1984): 1898-1908.

**Cholesky method:**

* Asmussen, Søren. Stochastic simulation with a view towards stochastic processes. University of Aarhus. Centre for Mathematical Physics and Stochastics (MaPhySto)[MPS], 1998.

**Davies Harte method:**

* Davies, Robert B., and D. S. Harte. "Tests for Hurst effect." Biometrika 74, no. 1 (1987): 95-101.
