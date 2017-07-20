#!/bin/bash

# Setting up Root environment....
export ROOTSYS=/Users/delaossa/local/root
export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

# HDF5: C++ interface
export HDF5CPP=/usr/local/opt/hdf5
export DYLD_LIBRARY_PATH=$HDF5CPP/lib:$DYLD_LIBRARY_PATH
export SZIP=/usr/local/opt/szip
export DYLD_LIBRARY_PATH=$SZIP/lib:$DYLD_LIBRARY_PATH

# ptools environment:
export PTOOLS=/Users/delaossa/plasma/software/ptools
export DYLD_LIBRARY_PATH=$PTOOLS/lib:$DYLD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PTOOLS/pymacros:$PATH
export PYTHONPATH=$PTOOLS/lib:$PYTHONPATH

