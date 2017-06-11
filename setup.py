#!/usr/bin/env python

from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='fbm',
      version='0.1.1',
      description='Fractional Brownian Motion',
      long_description=readme(),
      license='MIT',
      author='Christopher Flynn',
      author_email='crf204@gmail.com',
      url='https://github.com/crflynn/fbm',
      packages=['fbm'],
      zip_safe=False,
      install_requires=['numpy'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: MIT License',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: Implementation :: CPython',
      ],
      include_package_data=True)
