import numpy as np


def hosking(n, H=0.5, L=1):
    # returns fbm, fgn, times

    # Creates a discretely sampled fractional Brownian motion (fbm) realization with n increments, Hurst parameter H
    # and length L. Returns a 1 x (n+1) length array of values, since B(t=0) = 0, with the fractional Gaussian noise
    # realization as a 1 x n length array of values and a 1 x (n+1) array of time values corresponding to the fbm
    # realization.

    # Method of generation is Hosking's method (exact method) from his paper:
    # Hosking, J. R. (1984). Modeling persistence in hydrological time series using fractional differencing. Water
    # resources research, 20(12), 1898-1908.

    # Hosking's method generates a fractional Gaussian noise (fgn) realization. The cumulative sum of this realization
    # gives a fbm.

    # Input checking
    if not isinstance(n, int) or n <= 0:
        raise TypeError('Number of increments must be a positive integer')
    if H <= 0 or H >= 1:
        raise ValueError('Hurst parameter must be in interval (0, 1).')
    if L <= 0:
        raise ValueError('Length of fbm must be greater than 0.')

    # For scaling to interval [0, L]
    increment = float(L) / n
    scale = increment ** H

    gn = np.random.normal(0.0, 1.0, n)

    # If H = 0.5 then just generate a standard Brownian motion, otherwise proceed with Hosking's method
    if H == 0.5:
        fgn = gn
        del gn
    else:
        # Autocovariance function of fgn
        def gamma(k):
            return 0.5 * (abs(k - 1) ** (2 * H) - 2 * abs(k) ** (2 * H) + abs(k + 1) ** (2 * H))

        # Initializations
        fgn = np.zeros(n)
        phi = np.zeros(n)
        psi = np.zeros(n)
        cov = np.array([gamma(i) for i in range(n)])

        # First increment from stationary distribution
        fgn[0] = gn[0]
        v = 1
        phi[0] = 0

        # Generates fgn realization with n increments of size 1
        for i in range(1, n):
            phi[i - 1] = cov[i]
            for j in range(i - 1):
                psi[j] = phi[j]
                phi[i - 1] -= psi[j] * cov[i - j - 1]
            phi[i - 1] /= v
            for j in range(i - 1):
                phi[j] = psi[j] - phi[i - 1] * psi[i - j - 2]
            v *= (1 - phi[i - 1] * phi[i - 1])
            for j in range(i):
                fgn[i] += phi[j] * fgn[i - j - 1]
            fgn[i] += np.sqrt(v) * gn[i]

    # Scale to interval [0, L]
    fgn *= scale

    # Take cumulative sum, return fbm, fgn, timesteps
    fbm = fgn.cumsum()
    fbm = np.insert(fbm, [0], 0)
    t = np.linspace(0, L, n + 1)

    return fbm, fgn, t
