#!/bin/bash

# Setting up Root environment....
export ROOTSYS=/homea/hhh45/hhh451/software/root-6.06.08
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH

# HDF5: C++ interface
export HDF5CPP=/usr/local/software/jureca/Stages/2016a/software/HDF5/1.8.16-iimpi-8.2.5-GCC-4.9.3-2.25
export LD_LIBRARY_PATH=$HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
export PTOOLS=/homea/hhh45/hhh451/software/ptools
export LD_LIBRARY_PATH=$PTOOLS/lib:$LD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PTOOLS/pymacros:$PATH

