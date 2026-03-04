#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Setup script for pyTomo
Compiles Cython extensions
"""

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import sys
import os

# Define Cython extensions
# For conda-forge builds, we rely on the conda-build compiler flags
# For local conda builds, we use the conda compilers (gcc_linux-64) which provide
# the proper glibc with libmvec support for vectorized math functions
extensions = [
    Extension(
        name="geom_mat_gen_cython",
        sources=["geom_mat_gen_cython.pyx"],
        include_dirs=[np.get_include()],
        # Let the compiler optimize - conda's gcc includes proper libmvec support
        libraries=["m"] if sys.platform != "win32" else [],
        language="c",
    )
]

# Cythonize the extensions
ext_modules = cythonize(
    extensions,
    compiler_directives={
        'language_level': "3",  # Use Python 3 syntax
        'boundscheck': False,   # Disable bounds checking for speed
        'wraparound': False,    # Disable negative indexing for speed
        'cdivision': True,      # Use C division semantics
        'infer_types': True,    # Infer types automatically
    },
    annotate=False,  # Set to True to generate HTML annotation files for debugging
)

setup(
    name="pyTomo",
    version="1.0.0",
    description="Tomography reconstruction code for fusion plasmas",
    author="Tomas Odstrcil",
    ext_modules=ext_modules,
    include_dirs=[np.get_include()],
    zip_safe=False,
)