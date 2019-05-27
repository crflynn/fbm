# flake8: noqa
from fbm import MBM
import matplotlib.pyplot as plt
import time
import math


def h(t):
    # return 0.499*math.sin(t) + 0.5
    # return 0.6 * t + 0.3
    return 0.5 * math.exp(-8 * t ** 2) + 0.35


m = MBM(2 ** 8, h, 1)
t = m.times()
mbm = m.mbm()

plt.plot(t, mbm)
plt.plot(t, [h(tt) for tt in t])
plt.show()
