#!/bin/tcsh

# Global path for the Plasma sowftware
setenv SOFT /afs/desy.de/group/fla/plasma/data002/SOFTWARE    
    
# Setting up Root environment....
setenv ROOTSYS $SOFT/root
setenv LD_LIBRARY_PATH $ROOTSYS/lib:$LD_LIBRARY_PATH
setenv PATH $ROOTSYS/bin:$PATH

# HDF5: C++ interface
setenv HDF5CPP $SOFT/hdf5-1.8.14-c++
setenv LD_LIBRARY_PATH $HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
setenv PTOOLS $SOFT/ptools
setenv LD_LIBRARY_PATH $PTOOLS/lib:$LD_LIBRARY_PATH
setenv PATH $PTOOLS/bin:$PTOOLS/scripts:$PATH

