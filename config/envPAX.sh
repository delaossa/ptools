#!/bin/bash

# Setting up Root environment....
export ROOTSYS=$SOFT/root-6.06.08
#export ROOTSYS=/afs/desy.de/products/root/amd64_rhel60/6.02.00/
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH

# HDF5: C++ interface
export HDF5CPP=/afs/ifh.de/group/pitz/data/software/hdf5-1.8.7-noparallel
export LD_LIBRARY_PATH=$HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
export PTOOLS=/afs/ifh.de/group/pitz/data/software/ptools
export LD_LIBRARY_PATH=$PTOOLS/lib:$LD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PATH

