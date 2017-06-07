from unittest import TestCase

import numpy as np

from fbm import FBM, fbm, fgn, times


class FBMTests(TestCase):

    def setUp(self):
        self.f_daviesharte = FBM(n, H, L, method='daviesharte')
        self.f_hosking = FBM(n, H, L, method='hosking')
        self.f_cholesky = FBM(n, H, L, method='cholesky')

    def fbm_test(self):
        np.random.seed(42)
        realization = fbm(1024, 0.75, 1)
        self.assertEqual(realization, _)

    def fbm_daviesharte_test(self):
        np.random.seed(42)
        realization = fbm(n, H, L, method='daviesharte')
        self.assertEqual(realization, _)

    def fbm_hosking_test(self):
        np.random.seed(42)
        realization = fbm(n, H, L, method='hosking')
        self.assertEqual(realization, _)

    def fbm_cholesky_test(self):
        np.random.seed(42)
        realization = fbm(n, H, L, method='cholesky')
        self.assertEqual(realization, _)

    def fgn_daviesharte_test(self):
        np.random.seed(42)
        realization = fgn(n, H, L, method='daviesharte')
        self.assertEqual(realization, _)

    def fgn_hosking_test(self):
        np.random.seed(42)
        realization = fgn(n, H, L, method='hosking')
        self.assertEqual(realization, _)

    def fgn_cholesky_test(self):
        np.random.seed(42)
        realization = fgn(n, H, L, method='cholesky')
        self.assertEqual(realization, _)

    def times_test(self):
        times = times(n, L)
        self.assertEqual(times, _)

    def fbm_class_daviesharte_sample_test(self):
        np.random.seed(42)
        s = self.f_daviesharte.sample()
        self.assertEqual(s, _)

    def fbm_class_cholesky_sample_test(self):
        np.random.seed(42)
        s = self.f_cholesky.sample()
        self.assertEqual(s, _)

    def fbm_class_hosking_sample_test(self):
        np.random.seed(42)
        s = self.f_hosking.sample()
        self.assertEqual(s, _)

    def fbm_class_daviesharte_noise_test(self):
        np.random.seed(42)
        s = self.f_daviesharte.sample_noise()
        self.assertEqual(s, _)

    def fbm_class_cholesky_noise_test(self):
        np.random.seed(42)
        s = self.f_cholesky.sample_noise()
        self.assertEqual(s, _)

    def fbm_class_hosking_noise_test(self):
        np.random.seed(42)
        s = self.f_hosking.sample_noise()
        self.assertEqual(s, _)

    def fbm_class_daviesharte_times_test(self):
        np.random.seed(42)
        s = self.f_daviesharte.times()
        self.assertEqual(s, _)

    def fbm_class_cholesky_times_test(self):
        np.random.seed(42)
        s = self.f_cholesky.times()
        self.assertEqual(s, _)

    def fbm_class_hosking_times_test(self):
        np.random.seed(42)
        s = self.f_hosking.times()
        self.assertEqual(s, _)
