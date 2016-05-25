# PTOOLS

## Overview

PTOOLS is a software package written in C++ that allows for a direct handling, analysis and plotting of the output data from Osiris and HiPACE simulations.
In addition to the standard C++ libraries, PTOOLS uses [HDF5](http://www.hdfgroup.org/HDF5) to read data files and [ROOT](https://root.cern.ch) for data analysis and plotting.

## Installation

The PTOOLS source code is stored and maintained through this git repository :
```https://flauser1.desy.de/gitlab/delaossa/ptools.git
```

One can fork the project to a local copy of the source code by doing:
```git clone git@flauser1.desy.de:delaossa/ptools.git
```

## Usage

The routines of openPMD viewer can be used in two ways :

- Use the **Python API**, in order to write a script that loads the
  data and produces a set of pre-defined plots.

- Use the **interactive GUI inside the IPython Notebook**, in order to interactively
visualize the data.
