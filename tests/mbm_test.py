"""Test the MBM class."""
# flake8: noqa
import pytest

import numpy as np

from fbm import MBM
from fbm import mbm
from fbm import mgn


def test_MBM_init(n_good, hurst_func_good, length_good, mbm_method_good):
    m = MBM(n_good, hurst_func_good, length_good, mbm_method_good)
    print(str(m))
    print(repr(m))

def test_MBM_change(n_good, hurst_func_good, length_good, mbm_method_good):
    m = MBM(n_good, hurst_func_good, length_good, mbm_method_good)
    m.n = 42
    assert m._changed == True
    sample = m.mbm()
    assert m._changed == False

def test_MBM_init_n_bad(n_bad, hurst_func_good, length_good, mbm_method_good):
    with pytest.raises((TypeError, ValueError)):
        m = MBM(n_bad, hurst_func_good, length_good, mbm_method_good)

def test_MBM_init_hurst_bad(n_good, hurst_func_bad, length_good, mbm_method_good):
    with pytest.raises((TypeError, ValueError)):
        m = MBM(n_good, hurst_func_bad, length_good, mbm_method_good)

def test_MBM_init_length_bad(n_good, hurst_func_good, length_bad, mbm_method_good):
    with pytest.raises((TypeError, ValueError)):
        m = MBM(n_good, hurst_func_good, length_bad, mbm_method_good)

def test_MBM_init_method_bad(n_good, hurst_func_good, length_good, method_bad):
    with pytest.raises((TypeError, ValueError)):
        m = MBM(n_good, hurst_func_good, length_good, method_bad)

def test_MBM_mbm(n_good, hurst_func_good, length_good, mbm_method_good):
    m = MBM(n_good, hurst_func_good, length_good, mbm_method_good)
    mbm_sample = m.mbm()
    assert isinstance(mbm_sample, np.ndarray)
    assert len(mbm_sample) == n_good + 1

def test_MBM_mgn(n_good, hurst_func_good, length_good, mbm_method_good):
    m = MBM(n_good, hurst_func_good, length_good, mbm_method_good)
    mgn_sample = m.mgn()
    assert isinstance(mgn_sample, np.ndarray)
    assert len(mgn_sample) == n_good

def test_MBM_times(n_good, hurst_func_good, length_good, mbm_method_good):
    m = MBM(n_good, hurst_func_good, length_good, mbm_method_good)
    ts = m.times()
    assert isinstance(ts, np.ndarray)
    assert len(ts) == n_good + 1

def test_mbm(n_good, hurst_func_good, length_good, mbm_method_good):
    mbm_sample = mbm(n_good, hurst_func_good, length_good, mbm_method_good)
    assert isinstance(mbm_sample, np.ndarray)
    assert len(mbm_sample) == n_good + 1

def test_mgn(n_good, hurst_func_good, length_good, mbm_method_good):
    mgn_sample = mgn(n_good, hurst_func_good, length_good, mbm_method_good)
    assert isinstance(mgn_sample, np.ndarray)
    assert len(mgn_sample) == n_good
