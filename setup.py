#! /usr/bin/env python

import os
from setuptools import setup


def read_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


setup(name='tefingerprint',
      version='0.1.4',
      author='Tim Millar',
      author_email='tim.millar@plantandfood.co.nz',
      url='https://github.com/PlantandFoodResearch/TEFingerprint',
      download_url='https://github.com/PlantandFoodResearch/TEFingerprint/archive/v0.1.4.tar.gz',
      description='Toolkit for identifying transposon movement',
      long_description=read_file('README.rst'),
      scripts=['applications/tefingerprint',
               'applications/tef-extract-informative',
               'applications/tef-filter-gff'],
      packages=['tefingerprint'],
      keywords=['biology',
                'bioinformatics',
                'transposon',
                'transposable element',
                'cluster'],
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Science/Research',
                   'Programming Language :: Python :: 3.4',
                   'Topic :: Scientific/Engineering',
                   'Operating System :: Unix']
      )
