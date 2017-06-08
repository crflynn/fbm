#!/usr/bin/env python

from distutils.core import setup


def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='fbm',
      version='0.1.0',
      description='Fractional Brownian Motion',
      long_description=readme(),
      license='MIT',
      author='Christopher Flynn',
      url='https://github.com/crflynn/fbm',
      py_modules=['fbm'],
      zip_safe=False,
      install_requires=['numpy']
      classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: MIT License',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.6',
      ],
      include_package_data=True)
