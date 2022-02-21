#!/bin/bash

# Global path for the Plasma sowftware
export SOFT=/data/netapp/fla/plasma/software-max

# Setting up Root environment....
export ROOTSYS=$SOFT/root-6.06.08
# ROOT setup
alias thisroot='source $ROOTSYS/bin/thisroot.sh'
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH

# HDF5: C++ interface
export HDF5CPP=$SOFT/hdf5-1.8.14-c++
export LD_LIBRARY_PATH=$HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
export PTOOLS=$SOFT/ptools
export LD_LIBRARY_PATH=$PTOOLS/lib:$LD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PATH
