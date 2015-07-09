#!/bin/bash

# Setting up Root environment....
export ROOTSYS=/data/netapp/fla/plasma/software/root
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH

# HDF5: C++ interface
export HDF5CPP=/data/fla/plasma/software/hdf5-1.8.11-noparallel
export LD_LIBRARY_PATH=$HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
export PTOOLS=/data/netapp/fla/plasma/software/ptools
export LD_LIBRARY_PATH=$PTOOLS/lib:$LD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PATH

