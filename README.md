# fbm
Exact methods for simulating fractional Brownian motion (fBm) in python.

The three methods are Hosking's method, the Cholesky method, and the Davies Harte method. All three methods are exact in generating a discretely sampled fbm using python.

Usage:

```python
fbm_realization, fbm_increments, times = fbm(n, H=0.5, L=1, method='daviesharte')
```

where `n` is the number of equispaced increments desired for a fBm with Hurst parameter `H` on the interval [0, `L`]. Method can be either 'hosking','cholesky', or 'daviesharte'.

The Hosking and Cholesky methods are mathematically the same. The Cholesky script uses the Cholesky decomposition method from numpy's linear algebra library, while Hosking's method performs the same computations directly, which is slightly faster. For best performance use the Davies and Harte method, which is much faster than both other methods especially for larger increment quantities.


Hosking's method

Hosking, Jonathan RM. "Modeling persistence in hydrological time series using fractional differencing." Water resources research 20, no. 12 (1984): 1898-1908.

Cholesky method:

Asmussen, Søren. Stochastic simulation with a view towards stochastic processes. University of Aarhus. Centre for Mathematical Physics and Stochastics (MaPhySto)[MPS], 1998.

Davies Harte method:

Davies, Robert B., and D. S. Harte. "Tests for Hurst effect." Biometrika 74, no. 1 (1987): 95-101.
