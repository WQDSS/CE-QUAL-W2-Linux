# This is a port of the CE-QUAL-W2 model to linux
## Changes from vanilla version
1. Added a preprocessor flag to remove any dependency on windows UI contents
1. Added a Makefile to compile on Linux

## Instructions for compiling windows CLI only:
1. Create a new release type in the VS project, and enable the preprocessor flag `CLI_ONLY=1`

## Instructions for compiling in Linux:
1. install intel fortran compiler `ifort`
1. install gnu Make
1. Run `make renames` - this renames some of the files which have space in their names
1. Run `make w2_exe_linux` - this will build the linux executable

## Known issues:
1. Compiling with gfortran doesn't work due to syntax used for some of the printouts
1. Compiling and linking with openmp and MKL is causing some issues right now. These flags are in the makefile but are currently disabled

