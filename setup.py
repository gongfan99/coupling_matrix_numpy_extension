#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# from distutils.core import setup, Extension
# import numpy as np

# ext_modules = [ Extension('mylib', sources = ['mylib.c']) ]

# setup(
        # name = 'mylib',
        # version = '1.0',
        # include_dirs = [np.get_include()], #Add Include path of numpy
        # ext_modules = ext_modules
      # )
      
from setuptools import setup, Extension, find_packages
import numpy as np

ext_modules = [ Extension('couplingmatrix',
include_dirs = [np.get_include()],
libraries=['npymath'],
library_dirs = [np.get_include() + '\..\lib'],
sources = ['src\couplingmatrix.cc']
) ]

setup(
    name = 'couplingmatrix',
    version = '1.0.0',
    zip_safe = True,
    ext_modules = ext_modules,
    
    # metadata for upload to PyPI
    author="Fan Gong",
    author_email="gongfan99@hotmail.com",
    description="This is a native package about coupling matrix",
    license="MIT",
    keywords="coupling matrix filter",
    url="http://example.com/HelloWorld/"
)