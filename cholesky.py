import numpy as np


def cholesky(n, H=0.5, L=1):
    # returns fbm, fgn, times

    # Creates a discretely sampled fractional Brownian motion (fbm) realization with n increments, Hurst parameter H
    # and length L. Returns a 1 x (n+1) length array of values, since B(t=0) = 0, with the fractional Gaussian noise
    # realization as a 1 x n length array of values and a 1 x (n+1) array of time values corresponding to the fbm
    # realization.

    # Uses Cholesky decomposition method (exact method) from:
    # Asmussen, S. (1998). Stochastic simulation with a view towards stochastic processes. University of Aarhus. Centre
    # for Mathematical Physics and Stochastics (MaPhySto)[MPS].

    # Hosking's method performs the same operations directly rather than taking a Cholesky decomposition of the
    # covariance matrix.

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

    fgn = np.random.normal(0.0, 1.0, n)

    # If H = 0.5 then just generate a standard Brownian motion, otherwise proceed with the Cholesky method
    if H == 0.5:
        pass
    else:
        # Autocovariance function of fgn
        def gamma(k):
            return 0.5 * (abs(k - 1) ** (2 * H) - 2 * abs(k) ** (2 * H) + abs(k + 1) ** (2 * H))

        # Generate covariance matrix
        G = np.matrix(np.zeros([n, n]))
        for i in range(n):
            for j in range(i + 1):
                G[i, j] = gamma(i - j)

        # Cholesky decomposition
        C = np.linalg.cholesky(G)

        # Generate fgn
        fgn = C * np.matrix(fgn).T
        fgn = np.squeeze(np.asarray(fgn))

    # Scale to interval [0, L]
    fgn *= scale

    # Take cumulative sum, return fbm, fgn, timesteps
    fbm = fgn.cumsum()
    fbm = np.insert(fbm, [0], 0)
    t = np.linspace(0, L, n + 1)

    return fbm, fgn, t
