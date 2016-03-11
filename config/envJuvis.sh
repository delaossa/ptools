#!/bin/bash

# Setting up Root environment....
export ROOTSYS=/home/hhh23/hhh231/plasma/software/root
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH

# HDF5: C++ interface
export HDF5CPP=/home/hhh23/hhh231/plasma/software/hdf5-1.8.14
export LD_LIBRARY_PATH=$HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
export PTOOLS=/home/hhh23/hhh231/plasma/software/ptools
export LD_LIBRARY_PATH=$PTOOLS/lib:$LD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PATH

