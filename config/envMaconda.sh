#!/bin/bash

# Adding Anaconda environment (ROOT, HDF5, SZIP) should be installed there....
export ROOTSYS=$CONDA_PREFIX
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

# HDF5: C++ interface
export HDF5CPP=$CONDA_PREFIX
export LD_LIBRARY_PATH=$HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
export PTOOLS=/Users/delaossa/plasma/software/ptools
export LD_LIBRARY_PATH=$PTOOLS/lib:$LD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PTOOLS/pymacros:$PATH
export PYTHONPATH=$PTOOLS/pymacros:$PYTHONPATH

