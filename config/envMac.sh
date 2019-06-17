#!/bin/bash

# Setting up Root environment....
export ROOTSYS=/usr/local/opt/root/
export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

# HDF5: C++ interface
export HDF5CPP=/usr/local/opt/hdf5@1.8
export DYLD_LIBRARY_PATH=$HDF5CPP/lib:$DYLD_LIBRARY_PATH
export SZIP=/usr/local/opt/szip
export DYLD_LIBRARY_PATH=$SZIP/lib:$DYLD_LIBRARY_PATH

# ptools environment:
export PTOOLS=/Users/delaossa/plasma/software/ptools
export DYLD_LIBRARY_PATH=$PTOOLS/lib:$DYLD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PATH
export PYTHONPATH=$PTOOLS/lib:$PTOOLS/pymacros:$PYTHONPATH

