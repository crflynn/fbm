import numpy as np


def daviesharte(n, H=0.5, L=1):
    # returns fbm, fgn, times

    # Creates a discretely sampled fractional Brownian motion (fbm) realization with n increments, Hurst parameter H
    # and length L. Returns a 1 x (n+1) length array of values, since B(t=0) = 0, with the fractional Gaussian noise
    # realization as a 1 x n length array of values and a 1 x (n+1) array of time values corresponding to the fbm
    # realization.

    # Uses Davies and Harte method (exact method) from:
    # Davies, Robert B., and D. S. Harte. "Tests for Hurst effect." Biometrika 74, no. 1 (1987): 95-101.

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

    # If H = 0.5 then just generate a standard Brownian motion, otherwise proceed with the Davies Harte method
    if H == 0.5:
        pass
    else:
        # Autocovariance function of fgn
        def gamma(k):
            return 0.5 * (abs(k - 1) ** (2 * H) - 2 * abs(k) ** (2 * H) + abs(k + 1) ** (2 * H))

        # Generate first row of circulant matrix
        row_component = [gamma(i) for i in range(1, n)]
        reverse_component = [row_component[-i] for i in range(1, n)]
        row = [gamma(0)] + row_component + [0] + reverse_component

        # Get eigenvalues of circulant matrix
        # Discard imaginary part (should all be zero in theory so imaginary part will be very small)
        eigenvals = np.fft.fft(row).real

        # Generate second sequence of i.d.d. standard normals
        fgn2 = np.random.normal(0.0, 1.0, n)

        # Resulting sequence from matrix multiplication of positive definite sqrt(C) matrix with fgn sample can be
        # simulated in this way.
        w = np.zeros(2 * n, dtype=complex)
        for i in range(2 * n):
            if i == 0:
                w[i] = np.sqrt(eigenvals[i] / (2 * n)) * fgn[i]
            elif i < n:
                w[i] = np.sqrt(eigenvals[i] / (4 * n)) * (fgn[i] + 1j * fgn2[i])
            elif i == n:
                w[i] = np.sqrt(eigenvals[i] / (2 * n)) * fgn2[0]
            else:
                w[i] = np.sqrt(eigenvals[i] / (4 * n)) * (fgn[2 * n - i] - 1j * fgn2[2 * n - i])

        # Resulting z is fft of sequence w. Discard small imaginary part (z should be real in theory).
        z = np.fft.fft(w)
        fgn = z[:n].real

    # Scale to interval [0, L]
    fgn *= scale

    # Take cumulative sum, return fbm, fgn, timesteps
    fbm = fgn.cumsum()
    fbm = np.insert(fbm, [0], 0)
    t = np.linspace(0, L, n + 1)

    return fbm, fgn, t
