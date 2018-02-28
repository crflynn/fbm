"""Generate realizations of fractional Brownian motion."""
import warnings

import numpy as np


class FBM(object):
    """The FBM class.

    After instantiating with n = number of increments, hurst parameter, length
    of realization (default = 1) and method of generation
    (default daviesharte), call fbm() for fBm, fgn()
    for fGn, or times() to get corresponding time values.
    """

    def __init__(self, n, hurst, length=1, method="daviesharte"):
        """Instantiate the FBM."""
        self._methods = {"daviesharte": self._daviesharte,
                         "cholesky": self._cholesky,
                         "hosking": self._hosking}
        self.n = n
        self.hurst = hurst
        self.length = length
        self.method = method
        self._fgn = self._methods[self.method]
        # Some reusable values to speed up Monte Carlo.
        self._cov = None
        self._eigenvals = None
        self._C = None
        # Flag if some params get changed
        self._changed = False

    def __str__(self):
        """Str method."""
        return "fBm (" + str(self.method) + ") on [0, " + str(self.length) + \
            "] with Hurst value " + str(self.hurst) + " and " + \
            str(self.n) + " increments"

    def __repr__(self):
        """Repr method."""
        return "FBM(n=" + str(self.n) + ", hurst=" + str(self.hurst) + \
            ", length=" + str(self.length) + ", method=\"" + \
            str(self.method) + "\")"

    @property
    def n(self):
        """Get the number of increments."""
        return self._n

    @n.setter
    def n(self, value):
        if not isinstance(value, int) or value <= 0:
            raise TypeError("Number of increments must be a positive int.")
        self._n = value
        self._changed = True

    @property
    def hurst(self):
        """Hurst parameter."""
        return self._hurst

    @hurst.setter
    def hurst(self, value):
        if not isinstance(value, float) or value <= 0 or value >= 1:
            raise ValueError("Hurst parameter must be in interval (0, 1).")
        self._hurst = value
        self._changed = True

    @property
    def length(self):
        """Get the length of process."""
        return self._length

    @length.setter
    def length(self, value):
        if not isinstance(value, (int, float)) or value <= 0:
            raise ValueError("Length of fbm must be greater than 0.")
        self._length = value
        self._changed = True

    @property
    def method(self):
        """Get the algorithm used to generate."""
        return self._method

    @method.setter
    def method(self, value):
        if value not in self._methods:
            raise ValueError("Method must be 'daviesharte', 'hosking' or \
                             'cholesky'.")
        self._method = value
        self._fgn = self._methods[self.method]
        self._changed = True

    def fbm(self):
        """Sample the fractional Brownian motion."""
        return np.insert(self.fgn().cumsum(), [0], 0)

    def fgn(self):
        """Sample the fractional Gaussian noise."""
        scale = (1.0 * self.length / self.n) ** self.hurst
        gn = np.random.normal(0.0, 1.0, self.n)

        # If hurst == 1/2 then just return Gaussian noise
        if self.hurst == 0.5:
            return gn * scale
        else:
            fgn = self._fgn(gn)

        # Scale to interval [0, L]
        return fgn * scale

    def times(self):
        """Get times associated with the fbm/fgn samples."""
        return np.linspace(0, self.length, self.n + 1)

    def _autocovariance(self, k):
        """Autocovariance for fgn."""
        return 0.5 * (abs(k - 1) ** (2 * self.hurst) -
                      2 * abs(k) ** (2 * self.hurst) +
                      abs(k + 1) ** (2 * self.hurst))

    def _daviesharte(self, gn):
        """Generate a fgn realization using Davies-Harte method.

        Uses Davies and Harte method (exact method) from:
        Davies, Robert B., and D. S. Harte. "Tests for Hurst effect."
        Biometrika 74, no. 1 (1987): 95-101.

        Can fail if n is small and hurst close to 1. Falls back to Hosking
        method in that case. See:

        Wood, Andrew TA, and Grace Chan. "Simulation of stationary Gaussian
        processes in [0, 1] d." Journal of computational and graphical
        statistics 3, no. 4 (1994): 409-432.
        """
        # Monte carlo consideration
        if self._eigenvals is None or self._changed:
            # Generate the first row of the circulant matrix
            row_component = [self._autocovariance(i) for i in range(1, self.n)]
            reverse_component = list(reversed(row_component))
            row = [self._autocovariance(0)] + row_component + \
                  [0] + reverse_component

            # Get the eigenvalues of the circulant matrix
            # Discard the imaginary part (should all be zero in theory so
            # imaginary part will be very small)
            self._eigenvals = np.fft.fft(row).real
            self._changed = False

        # If any of the eigenvalues are negative, then the circulant matrix
        # is not positive definite, meaning we cannot use this method. This
        # occurs for situations where n is low and H is close to 1.
        # Fall back to using the Hosking method. See the following for a more
        # detailed explanation:
        #
        # Wood, Andrew TA, and Grace Chan. "Simulation of stationary Gaussian
        #     processes in [0, 1] d." Journal of computational and graphical
        #     statistics 3, no. 4 (1994): 409-432.
        if np.any([ev < 0 for ev in self._eigenvals]):
            warnings.warn(
                "Combination of increments n and Hurst value H "
                "invalid for Davies-Harte method. Reverting to Hosking method."
                " Occurs when n is small and Hurst is close to 1. ")
            # Set method to hosking for future samples.
            self.method = "hosking"
            # Don"t need to store eigenvals anymore.
            self._eigenvals = None
            return self._hosking(gn)

        # Generate second sequence of i.i.d. standard normals
        gn2 = np.random.normal(0.0, 1.0, self.n)

        # Resulting sequence from matrix multiplication of positive definite
        # sqrt(C) matrix with fgn sample can be simulated in this way.
        w = np.zeros(2 * self.n, dtype=complex)
        for i in range(2 * self.n):
            if i == 0:
                w[i] = np.sqrt(self._eigenvals[i] / (2 * self.n)) * gn[i]
            elif i < self.n:
                w[i] = np.sqrt(self._eigenvals[i] / (4 * self.n)) * \
                    (gn[i] + 1j * gn2[i])
            elif i == self.n:
                w[i] = np.sqrt(self._eigenvals[i] / (2 * self.n)) * gn2[0]
            else:
                w[i] = np.sqrt(self._eigenvals[i] / (4 * self.n)) * \
                    (gn[2 * self.n - i] - 1j * gn2[2 * self.n - i])

        # Resulting z is fft of sequence w. Discard small imaginary part (z
        # should be real in theory).
        z = np.fft.fft(w)
        fgn = z[:self.n].real
        return fgn

    def _cholesky(self, gn):
        """Generate a fgn realization using the Cholesky method.

        Uses Cholesky decomposition method (exact method) from:
        Asmussen, S. (1998). Stochastic simulation with a view towards
        stochastic processes. University of Aarhus. Centre for Mathematical
        Physics and Stochastics (MaPhySto)[MPS].
        """
        # Monte carlo consideration
        if self._C is None or self._changed:
            # Generate covariance matrix
            G = np.matrix(np.zeros([self.n, self.n]))
            for i in range(self.n):
                for j in range(i + 1):
                    G[i, j] = self._autocovariance(i - j)

            # Cholesky decomposition
            self._C = np.linalg.cholesky(G)
            self._changed = False

        # Generate fgn
        fgn = self._C * np.matrix(gn).T
        fgn = np.squeeze(np.asarray(fgn))
        return fgn

    def _hosking(self, gn):
        """Generate a fGn realization using Hosking's method.

        Method of generation is Hosking's method (exact method) from his paper:
        Hosking, J. R. (1984). Modeling persistence in hydrological time series
        using fractional differencing. Water resources research, 20(12),
        1898-1908.
        """
        fgn = np.zeros(self.n)
        phi = np.zeros(self.n)
        psi = np.zeros(self.n)
        # Monte carlo consideration
        if self._cov is None or self._changed:
            self._cov = np.array(
                [self._autocovariance(i) for i in range(self.n)])
            self._changed = False

        # First increment from stationary distribution
        fgn[0] = gn[0]
        v = 1
        phi[0] = 0

        # Generate fgn realization with n increments of size 1
        for i in range(1, self.n):
            phi[i - 1] = self._cov[i]
            for j in range(i - 1):
                psi[j] = phi[j]
                phi[i - 1] -= psi[j] * self._cov[i - j - 1]
            phi[i - 1] /= v
            for j in range(i - 1):
                phi[j] = psi[j] - phi[i - 1] * psi[i - j - 2]
            v *= (1 - phi[i - 1] * phi[i - 1])
            for j in range(i):
                fgn[i] += phi[j] * fgn[i - j - 1]
            fgn[i] += np.sqrt(v) * gn[i]

        return fgn


def fbm(n, hurst, length=1, method="daviesharte"):
    """One off sample of fBm."""
    f = FBM(n, hurst, length, method)
    return f.fbm()


def fgn(n, hurst, length=1, method="daviesharte"):
    """One off sample of fGn."""
    f = FBM(n, hurst, length, method)
    return f.fgn()


def times(n, length=1):
    """Generate the times associated with increments n and length."""
    return np.linspace(0, length, n + 1)
