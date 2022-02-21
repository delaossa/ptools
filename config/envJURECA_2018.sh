#!/bin/bash

SOFT=/homea/hhh45/hhh451/software-jureca

# Setting up Root environment....
export ROOTSYS=$SOFT/root-6.14.04
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

# HDF5: C++ interface
export HDF5CPP=/usr/local/software/jureca/Stages/2018a/software/HDF5/1.8.20-iimpi-2018a
export LD_LIBRARY_PATH=$HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
export PTOOLS=$SOFT/ptools
export LD_LIBRARY_PATH=$PTOOLS/lib:$LD_LIBRARY_PATH
export PATH=$PTOOLS/bin:$PTOOLS/scripts:$PTOOLS/pymacros:$PATH

