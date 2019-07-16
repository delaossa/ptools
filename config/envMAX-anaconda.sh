#!/bin/bash

# Adding Anaconda environment (ROOT and HDF5) should be installed there....
export ROOTSYS=$CONDA_PREFIX

# HDF5: C++ interface
export HDF5CPP=$CONDA_PREFIX

# ptools environment:
export PTOOLS=/data/netapp/fla/plasma/software-max/ptools-ana
export DYLD_LIBRARY_PATH=$PTOOLS/lib:$DYLD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PTOOLS/pymacros:$PATH
export PYTHONPATH=$PTOOLS/pymacros:$PYTHONPATH
