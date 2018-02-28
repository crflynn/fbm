"""Pytest fixtures."""
# flake8: noqa
import math

import pytest


@pytest.fixture(params=[-1, 0, 1.5, '1'])
def n_bad(request):
    return request.param

@pytest.fixture(params=[4, 16, 256])
def n_good(request):
    return request.param

@pytest.fixture(params=[-1, 0, '1'])
def length_bad(request):
    return request.param

@pytest.fixture(params=[0.5, 1, 2])
def length_good(request):
    return request.param

@pytest.fixture(params=[-1, 0, 1, 2, '0.5'])
def hurst_bad(request):
    return request.param

@pytest.fixture(params=[0.25, 0.5, 0.75])
def hurst_good(request):
    return request.param

# Bad functions for Hurst in mBm
def hb1(t):
    return -0.5 + t

def hb2(t):
    return 0.5 + t

def hb3(t, b):
    return t + b

def hb4():
    return 0.5

@pytest.fixture(params=[hb1, hb2, hb3, hb4, 0.5, '0.5'])
def hurst_func_bad(request):
    return request.param

# Acceptable functions for Hurst in mBm
def h1(t):
    return 0.3 + 0.6 * (t % 1)

def h2(t):
    if t <= 0.5:
        return 0.25
    else:
        return 0.75

def h3(t):
    return 0.25 * math.sin(t) + 0.5

@pytest.fixture(params=[h1, h2, h3])
def hurst_func_good(request):
    return request.param

@pytest.fixture(params=[4])
def n_fallback(request):
    return request.param

@pytest.fixture(params=[0.99])
def hurst_fallback(request):
    return request.param

@pytest.fixture(params=['hosking', 'cholesky', 'daviesharte'])
def fbm_method_good(request):
    return request.param

@pytest.fixture(params=['riemannliouville'])
def mbm_method_good(request):
    return request.param

@pytest.fixture(params=['bad_method'])
def method_bad(request):
    return request.param
