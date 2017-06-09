from unittest import TestCase

import numpy as np

from fbm import FBM, fbm, fgn, times


class FBMTests(TestCase):

    def setUp(self):
        # For seed = 42, n = 5, H = 0.7, L = 1
        self.realizations = {
            'fbm': {
                'daviesharte': np.array([0.        ,  0.29326393,  0.57918201,
                                         0.69236348,  1.17425979,  1.33485224
                                         ]),
                'cholesky': np.array([0.        ,  0.16100061,  0.16997488,
                                      0.38610856,  0.92328558,  1.02602397]),
                'hosking': np.array([ 0.        ,  0.16100061,  0.16997488,
                                      0.38610856,  0.92328558,  1.02602397]),
            },
            'fgn': {
                'daviesharte': np.array([0.29326393,  0.28591808,  0.11318147,
                                         0.48189631,  0.16059245]),
                'cholesky': np.array([0.16100061,  0.00897426,  0.21613369,
                                      0.53717702,  0.10273839]),
                'hosking': np.array([0.16100061,  0.00897426,  0.21613369,
                                     0.53717702,  0.10273839]),
            }
        }
        self.seed = 42
        self.n = 5
        self.H = 0.7
        self.L = 1
        self.t = np.linspace(0, self.L, self.n + 1)
        self.cls_objects = {
            'daviesharte': FBM(self.n, self.H, self.L, method='daviesharte'),
            'cholesky': FBM(self.n, self.H, self.L, method='cholesky'),
            'hosking': FBM(self.n, self.H, self.L, method='hosking'),
        }

    def compare(self, a, b):
        """Just make sure they are within the display precision."""
        print(a)
        print(b)
        print(abs(a - b) < 0.00000001)
        return (abs(a - b) < 0.00000001).all()

    def realization_compare(self, func, process, method):
        np.random.seed(self.seed)
        realization = func(self.n, self.H, self.L, method=method)
        self.assertTrue(self.compare(realization,
                                     self.realizations[process][method]))

    def fbm_daviesharte_test(self):
        self.realization_compare(fbm, 'fbm', 'daviesharte')

    def fbm_cholesky_test(self):
        self.realization_compare(fbm, 'fbm', 'cholesky')

    def fbm_hosking_test(self):
        self.realization_compare(fbm, 'fbm', 'hosking')

    def fgn_daviesharte_test(self):
        self.realization_compare(fgn, 'fgn', 'daviesharte')

    def fgn_cholesky_test(self):
        self.realization_compare(fgn, 'fgn', 'cholesky')

    def fgn_hosking_test(self):
        self.realization_compare(fgn, 'fgn', 'hosking')

    def times_test(self):
        t = times(self.n, self.L)
        self.assertTrue((t == self.t).all())

    def fbm_class_compare(self, process, method):
        np.random.seed(self.seed)
        obj = self.cls_objects[method]
        realization = obj.fbm()
        self.assertTrue(self.compare(realization,
                                     self.realizations[process][method]))

    def fgn_class_compare(self, process, method):
        np.random.seed(self.seed)
        obj = self.cls_objects[method]
        realization = obj.fgn()
        self.assertTrue(self.compare(realization,
                                     self.realizations[process][method]))

    def fbm_class_daviesharte_sample_test(self):
        self.fbm_class_compare('fbm', 'daviesharte')

    def fbm_class_cholesky_sample_test(self):
        self.fbm_class_compare('fbm', 'cholesky')

    def fbm_class_hosking_sample_test(self):
        self.fbm_class_compare('fbm', 'hosking')

    def fgn_class_daviesharte_sample_test(self):
        self.fgn_class_compare('fgn', 'daviesharte')

    def fgn_class_cholesky_sample_test(self):
        self.fgn_class_compare('fgn', 'cholesky')

    def fgn_class_hosking_sample_test(self):
        self.fgn_class_compare('fgn', 'hosking')

    def fbm_class_times_test(self):
        obj = self.cls_objects['daviesharte']
        self.assertTrue((obj.times() == self.t).all())
