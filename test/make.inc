#  make include file
#
SHELL = /bin/sh
#
#  The machine (platform) identifier to append to the library names
#
PLAT =
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
FORTRAN  = gfortran
CC       = gcc
OPTS     = -fPIC -O3 -funroll-all-loops
DRVOPTS  = $(OPTS)
NOOPT    =
LOADER   = gfortran
LOADOPTS =
INCLUDE = include/
INCLUDE_UTILS = ../../../../include/api/utils/
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib
#
#  The location and name of the Reference library.
#
VARLIB      = libnumx_mcvar$(PLAT).a
