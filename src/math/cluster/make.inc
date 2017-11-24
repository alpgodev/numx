#  make include file
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
OPTS     = -fPIC -x f77-cpp-input -O3 -funroll-all-loops
# -Wall -Wno-unused-function -Wno-unused-dummy-argument -Wno-uninitialized -Wno-maybe-uninitialized -Wno-unused-label -Wno-conversion
DRVOPTS  = $(OPTS)
#
#  The location and name of the Reference library.
#
LIB_NAME = libnumx_cluster.a
CLUSTERLIB = ../../../lib/math/$(PLAT)/$(LIB_NAME)