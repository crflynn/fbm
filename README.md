# fbm
Exact methods for simulating fractional Brownian motion (fBm) in python.

The three methods are Hosking's method, the Cholesky method, and the Davies Harte method. All three methods are exact in generating a discretely sampled fbm using python.

Usage:

```python
fbmr, fgnr, times = fbm(n, H=0.5, L=1, method='daviesharte')
```

where `n` is the number of equispaced increments desired for a fBm with Hurst parameter `H` on the interval [0, `L`]. Method can be either 'hosking','cholesky', or 'daviesharte'. The function returns 

* `fbmr`, a list of values of a discretely sampled fbm realization,
* `fgnr`, a list of increments of the fbm realization (fractional Gaussian noise)
* `times`, a list of time values corresponding to the values of `fbmr`

The Hosking and Cholesky methods are mathematically the same. The Cholesky function uses the Cholesky decomposition method from numpy's linear algebra library, while Hosking's method performs the same computations directly, which is slightly faster. For best performance use the Davies and Harte method, which is much faster than both other methods especially for larger increment quantities.

The Davies and Harte method can fail if the Hurst parameter `H` is close to 1 and there are a small amount of increments `n`. If this occurs, python will print a warning to the console and fallback to using Hosking's method to generate the realization. See page 412 of the following paper for a more detailed explanation

* Wood, Andrew TA, and Grace Chan. "Simulation of stationary Gaussian processes in [0, 1] d." Journal of computational and graphical statistics 3, no. 4 (1994): 409-432.


**Hosking's method:**

* Hosking, Jonathan RM. "Modeling persistence in hydrological time series using fractional differencing." Water resources research 20, no. 12 (1984): 1898-1908.

**Cholesky method:**

* Asmussen, Søren. Stochastic simulation with a view towards stochastic processes. University of Aarhus. Centre for Mathematical Physics and Stochastics (MaPhySto)[MPS], 1998.

**Davies Harte method:**

* Davies, Robert B., and D. S. Harte. "Tests for Hurst effect." Biometrika 74, no. 1 (1987): 95-101.
