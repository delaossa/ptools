# PTOOLS

## Overview

PTOOLS is a software package written in C++ that allows for a direct handling, analysis and plotting of the output data from Osiris and HiPACE simulations.
In addition to the standard C++ libraries, PTOOLS uses [HDF5](http://www.hdfgroup.org/HDF5) to read data files and [ROOT](https://root.cern.ch) for data analysis and plotting.

## Installation

### Basic installation

To install this package :

- Clone this repository using `git`
---

git clone git@flauser1.desy.de:delaossa/ptools.git
```

- `cd` into the directory `ptools` and run
```
python setup.py install
```

## Usage

The routines of openPMD viewer can be used in two ways :

- Use the **Python API**, in order to write a script that loads the
  data and produces a set of pre-defined plots.

- Use the **interactive GUI inside the IPython Notebook**, in order to interactively
visualize the data.

#### Tutorials

The notebooks in the folder `tutorials/` demonstrate how to use both
the API and the interactive GUI. You can view these notebooks online
[here](https://github.com/openPMD/openPMD-viewer/tree/master/tutorials),
or, alternatively, you can run them on your local computer by typing:

`ipython notebook tutorials/`

NB: For [NERSC](http://www.nersc.gov/) users, you can run the tutorials on a
remote machine by logging in at
[https://ipython.nersc.gov](https://ipython.nersc.gov), and by
navigating to your personal copy of the directory `openPMD-viewer/tutorials`.

#### Notebook quick-starter

If you wish to use the **interactive GUI**, the installation of `openPMD-viewer` provides
a convenient executable which automatically
**creates a new pre-filled notebook** and **opens it in a
browser**. To use this executable, simply type in a regular terminal:

`openPMD_notebook`

(This executable is installed by default when running `python setup.py install`.)

## Contributing to the openPMD-viewer

We welcome contributions to the code! Please read [this page](https://github.com/openPMD/openPMD-viewer/blob/master/CONTRIBUTING.md) for
guidelines on how to contribute.

![travis badge](https://travis-ci.org/openPMD/openPMD-viewer.svg?branch=master)
