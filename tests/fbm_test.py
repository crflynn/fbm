"""Test the FBM class."""
# flake8: noqa
import pytest

import numpy as np

from fbm import FBM
from fbm import fbm
from fbm import fgn
from fbm import times


def test_FBM_init(n_good, hurst_good, length_good, fbm_method_good):
    f = FBM(n_good, hurst_good, length_good, fbm_method_good)
    print(str(f))
    print(repr(f))

def test_FBM_init_n_bad(n_bad, hurst_good, length_good, fbm_method_good):
    with pytest.raises((TypeError, ValueError)):
        f = FBM(n_bad, hurst_good, length_good, fbm_method_good)

def test_FBM_init_hurst_bad(n_good, hurst_bad, length_good, fbm_method_good):
    with pytest.raises((TypeError, ValueError)):
        f = FBM(n_good, hurst_bad, length_good, fbm_method_good)

def test_FBM_init_length_bad(n_good, hurst_good, length_bad, fbm_method_good):
    with pytest.raises((TypeError, ValueError)):
        f = FBM(n_good, hurst_good, length_bad, fbm_method_good)

def test_FBM_init_method_bad(n_good, hurst_good, length_good, method_bad):
    with pytest.raises((TypeError, ValueError)):
        f = FBM(n_good, hurst_good, length_good, method_bad)

def test_FBM_fbm(n_good, hurst_good, length_good, fbm_method_good):
    f = FBM(n_good, hurst_good, length_good, fbm_method_good)
    fbm_sample = f.fbm()
    assert isinstance(fbm_sample, np.ndarray)
    assert len(fbm_sample) == n_good + 1

def test_FBM_fgn(n_good, hurst_good, length_good, fbm_method_good):
    f = FBM(n_good, hurst_good, length_good, fbm_method_good)
    fgn_sample = f.fgn()
    assert isinstance(fgn_sample, np.ndarray)
    assert len(fgn_sample) == n_good

def test_FBM_times(n_good, hurst_good, length_good, fbm_method_good):
    f = FBM(n_good, hurst_good, length_good, fbm_method_good)
    ts = f.times()
    assert isinstance(ts, np.ndarray)
    assert len(ts) == n_good + 1

def test_FBM_method_fallback(n_fallback, hurst_fallback):
    f = FBM(5, 0.99, 1, method='daviesharte')
    with pytest.warns(Warning):
        sample = f.fbm()

def test_fbm(n_good, hurst_good, length_good, fbm_method_good):
    fbm_sample = fbm(n_good, hurst_good, length_good, fbm_method_good)
    assert isinstance(fbm_sample, np.ndarray)
    assert len(fbm_sample) == n_good + 1

def test_fgn(n_good, hurst_good, length_good, fbm_method_good):
    fgn_sample = fgn(n_good, hurst_good, length_good, fbm_method_good)
    assert isinstance(fgn_sample, np.ndarray)
    assert len(fgn_sample) == n_good

def test_times(n_good, length_good):
    ts = times(n_good, length_good)
    assert isinstance(ts, np.ndarray)
    assert len(ts) == n_good + 1
