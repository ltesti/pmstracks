#! /usr/bin/env python
# -*- coding: utf-8 -*-
# from numpy.distutils.core import setup, Extension
from setuptools import setup, find_packages
import os

# set the version number
with open('pmstracks/_version.py') as f:
    exec(f.read())

setup(
    name="pmstracks",
    version=__version__,
    author="Leonardo Testi",
    author_email="ltesti120a@gmail.com",
    #packages=find_packages(),
    package_dir={'pmstracks': 'pmstracks'},
    packages=['pmstracks'],
    url="https://github.com/ltesti/pmstracks/",
    license="GPL",
    description="Pre main sequence tracks utilities",
    long_description=open("README.md").read() + "\n\n",
    #package_data={"": ["LICENSE", "AUTHORS.rst"], 
    #              "pmstracks/tracks": ['*','*/*','*/*/*']},
    setup_requires=['setuptools_scm'],
    include_package_data=True,
    install_requires=["matplotlib","numpy", "scipy"],
    #data_files=[('pmstracks/tracks', ['*/*','*/*/*'])],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        'Programming Language :: Python :: 2',
    ],
    # ext_modules = ["DiscEvolution"]
)
