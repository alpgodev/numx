#  make include file
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
OPTS     = -fPIC -O3 -funroll-all-loops
DRVOPTS  = $(OPTS)
NOOPT    =
INCLUDE = include/
INCLUDE_UTILS = ../../../api/include/utils/
#
#  The location and name of the Reference library.
#
LIB_NAME = libnumx_mcvar.a
MCVARLIB  = ../../../../lib/math/$(PLAT)/$(LIB_NAME)
