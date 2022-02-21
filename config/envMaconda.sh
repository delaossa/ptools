#!/bin/bash

# Adding Anaconda environment (ROOT and HDF5) should be installed there....
export CONDA_BUILD_SYSROOT=/Users/delaossa/local/MacOSX10.9.sdk
export ROOTSYS=$CONDA_PREFIX

# HDF5: C++ interface
export HDF5CPP=$CONDA_PREFIX

# ptools environment:
export PTOOLS=/Users/delaossa/local/ptools
export DYLD_LIBRARY_PATH=$PTOOLS/lib
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PTOOLS/pymacros:$PATH
export PYTHONPATH=$PTOOLS/lib:$PYTHONPATH

