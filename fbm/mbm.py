"""Generate realizations of multifractional Brownian motion."""
import inspect
from math import gamma

import numpy as np


class MBM(object):
    """The MBM class.

    A class for generating multifractional Brownian motion or
    multifractional Gaussian noise using approximate methods.
    """

    def __init__(self, n, hurst, length=1, method="riemannliouville"):
        """Instantiate an MBM."""
        self._methods = {"riemannliouville": self._riemannliouville}
        self.n = n
        self.length = length
        self.hurst = hurst
        self.method = method
        self._mbm = self._methods[self.method]
        self._dt = 1.0 * self.length / self.n
        self._ts = self.times()
        # Flag if some params get changed
        self._changed = False

    def __str__(self):
        """Str method."""
        return "mBm (" + str(self.method) + ") on [0, " + str(self.length) + \
            "] with Hurst function " + self.hurst.__name__ + " and " + \
            str(self.n) + " increments"

    def __repr__(self):
        """Repr method."""
        return "MBM(n=" + str(self.n) + ", hurst=" + self.hurst.__name__ + \
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
        try:
            num_args = len(inspect.getargspec(value).args)
        except Exception:
            raise ValueError(
                "Hurst parameter must be a function of one argument.")
        if not callable(value) or num_args != 1:
            raise ValueError(
                "Hurst parameter must be a function of one argument.")
        self._check_hurst(value)
        self._hurst = value
        self._changed = True

    def _check_hurst(self, value):
        self._hs = [value(t) for t in self.times()]
        for h in self._hs:
            if h <= 0 or h >= 1:
                raise ValueError("Hurst range must be on interval (0, 1).")

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
            raise ValueError("Method must be ...")
        self._method = value
        self._mgn = self._methods[self.method]
        self._changed = True

    def mbm(self):
        """Generate a realization of multifractional Brownian motion."""
        return self._mbm()

    def mgn(self):
        """Generate a realization of multifractional Gaussian noise."""
        return np.diff(self.mbm())

    def _riemannliouville(self):
        """Generate Riemann-Liouville mBm."""
        gn = np.random.normal(0.0, 1.0, self.n)
        if self._changed:
            self._dt = 1.0 * self.length / self.n
            self._ts = self.times()
            self._check_hurst(self.hurst)
            self._changed = False
        mbm = [0]
        coefs = [(g / np.sqrt(self._dt)) * self._dt for g in gn]
        for k in range(1, self.n+1):
            weights = [self._w(t, self._hs[k]) for t in self._ts[1:k+1]]
            seq = [coefs[i-1] * weights[k - i] for i in range(1, k+1)]
            mbm.append(sum(seq))
        return np.array(mbm)

    def times(self):
        """Get times associated with the fbm/fgn samples."""
        return np.linspace(0, self.length, self.n + 1)

    def _w(self, t, hurst):
        """Get the Riemann-Liouville method weight for time t."""
        w = 1.0 / gamma(hurst + 0.5) * \
            np.sqrt((t ** (2 * hurst) -
                    (t - self._dt) ** (2 * hurst)) / (2 * hurst * self._dt))
        return w


def mbm(n, hurst, length=1, method="riemannliouville"):
    """One off sample of mBm."""
    m = MBM(n, hurst, length, method)
    return m.mbm()


def mgn(n, hurst, length=1, method="riemannliouville"):
    """One off sample of mGn."""
    m = MBM(n, hurst, length, method)
    return m.mgn()
