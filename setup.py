#!/usr/bin/env python3

from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='itsne',
      version='1.0.0',
      description='Compute interactive t-SNE plots from RNA-seq count data.',
      author='Clay McLeod',
      author_email='clay.mcleod@STJUDE.org',
      scripts=['scripts/itsne-main', 'scripts/itsne-normalize-matrix.R'],
      packages=find_packages(),
      install_requires=requirements,
      python_requires='>=3.0, <3.8')
