#!/bin/tcsh

# Setting up Root environment....
setenv ROOTSYS /home/delaossa/plasma/software/root-6.06.08
setenv LD_LIBRARY_PATH $ROOTSYS/lib:$LD_LIBRARY_PATH
setenv PATH $ROOTSYS/bin:$PATH
setenv PYTHONPATH $ROOTSYS/lib:$PYTHONPATH

# HDF5: C++ interface
setenv HDF5CPP /home/delaossa/plasma/software/hdf5-1.8.16-c++
setenv LD_LIBRARY_PATH $HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
setenv PTOOLS /home/delaossa/plasma/software/ptools
setenv LD_LIBRARY_PATH $PTOOLS/lib:$LD_LIBRARY_PATH
setenv PATH $PTOOLS/bin:$PTOOLS/scripts:$PTOOLS/pymacros:$PATH
setenv PYTHONPATH $PTOOLS/lib:$PYTHONPATH

