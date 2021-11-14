#!/usr/bin/env python3
# coding: utf-8
"""
@author: Ping Qiu  qiuping1@genomics.cn
@last modified by: Ping Qiu
@file:setup.py
@time:2021/03/02
"""
from setuptools import Extension, setup, find_packages
import sys
from pathlib import Path

# Newer packaging standards may recommend removing the current dir from the
# path, add it back if needed.
if '' not in sys.path:
    sys.path.insert(0, '')
import setup_build

if sys.version_info < (3, 7):
    sys.exit('stereopy requires Python >= 3.7')

setup(
    name='gefpy',
    version='0.2.1',
    setup_requires=['pkgconfig', 'Cython', 'setuptools_scm'],
    description='Spatial transcriptomic analysis in python.',
    long_description=Path('README.md').read_text('utf-8'),
    long_description_content_type="text/markdown",
    url='https://github.com/BGIResearch/gefpy',
    author='BGIResearch',
    author_email='huangzhibo@genomics.cn',
    python_requires='>=3.7',
    install_requires=[
        "pandas >= 1.3.3",
        "h5py >= 3.2.1",
        "setuptools >= 41.0.0",
        "opencv-python >= 4.5.4.58",
        "tifffile"
    ],
    extras_require=dict(
        doc=['sphinx>=3.2'],
        test=['pytest>=4.4', 'pytest-nunit'],
    ),
    packages=find_packages(),
    include_package_data=True,
    ext_modules=[Extension('gefpy.x', ['x.cpp'])],
    classifiers=[
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
    cmdclass = {'build_ext': setup_build.gefpy_build_ext}
)
