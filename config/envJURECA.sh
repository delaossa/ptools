#!/bin/bash

# Setting up Root environment....
export ROOTSYS=/homec/hhh23/hhh231/plasma/software/root-6.06.04
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH

# HDF5: C++ interface
export HDF5CPP=/homec/hhh23/hhh231/plasma/software/hdf5-1.8.14
export LD_LIBRARY_PATH=$HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
export PTOOLS=/homec/hhh23/hhh231/plasma/software/ptools
export LD_LIBRARY_PATH=$PTOOLS/lib:$LD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PTOOLS/pymacros:$PATH

