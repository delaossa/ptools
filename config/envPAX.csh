#!/bin/tcsh

# Setting up Root environment....
setenv ROOTSYS /data/netapp/fla/plasma/software/root
setenv LD_LIBRARY_PATH $ROOTSYS/lib:$LD_LIBRARY_PATH
setenv PATH $ROOTSYS/bin:$PATH

# HDF5: C++ interface
setenv HDF5CPP /afs/ifh.de/group/pitz/data/software/hdf5-1.8.7-noparallel
setenv LD_LIBRARY_PATH $HDF5CPP/lib:$LD_LIBRARY_PATH

# ptools environment:
setenv PTOOLS /afs/ifh.de/group/pitz/data/software/ptools
setenv LD_LIBRARY_PATH $PTOOLS/lib:$LD_LIBRARY_PATH
setenv PATH $PTOOLS/bin:$PTOOLS/scripts:$PATH

