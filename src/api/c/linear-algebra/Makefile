include ../../../../make.inc
include ../make.inc
#/*================================================================
#   NumX - C/C++ Linux/MacOsx build             
#
#   Module - Linear Algebra
#
#===================================================================*/

# Library name
LIB_NAME = libnumx-linear-algebra-c-1.0

# output libraries directory
SOLIB_NAME = $(LIB_PATH)/$(LIB_NAME).so
ALIB_NAME  = $(LIB_PATH)/$(LIB_NAME).a

all: $(SOLIB_NAME) $(ALIB_NAME)

# any extra libraries that your librarie may depend on.
ifeq ($(PLAT),macosx)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_rnd -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(MATH_LIBS) $(GCC_LIBS)
    LIB_SRCS := ndmc.c	nm1c.c  nmc.c   nminfc.c    pmc.c   pmcc.c	\
                pmmtc.c	pmtmc.c pmtvc.c	pmvc.c  pmxc.c	pvxc.c     \
                rchoc.c rmc.c	schursc.c   smc.c   tmc.c	xmc.c   \
                xvc.c   testlinalgc.c   \
                $(UTILS_PATH)/utils.c   $(PROTECT_PATH)/sysinfo-macosx.c
endif
ifeq ($(PLAT),linux)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_rnd -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(MATH_LIBS) $(GCC_LIBS)
    LIB_SRCS := ndmc.c	nm1c.c  nmc.c   nminfc.c    pmc.c   pmcc.c	\
                pmmtc.c	pmtmc.c pmtvc.c	pmvc.c  pmxc.c	pvxc.c     \
                rchoc.c rmc.c	schursc.c   smc.c   tmc.c	xmc.c   \
                xvc.c   testlinalgc.c   \
#                $(PROTECT_PATH)/decrypt-linux.c	$(PROTECT_PATH)/sysinfo-linux.c	\
                $(UTILS_PATH)/utils.c
endif

LIB_OBJS := $(LIB_SRCS:.c=.o)

# Making Libraries (this line is for static Lib)
$(SOLIB_NAME):  $(LIB_OBJS)
	$(SOLIBTOOL) $(SOLIB_LDFLAGS) -o $(SOLIB_NAME) $(LIB_OBJS) $(SOLIB_LDLIBS)

$(ALIB_NAME): $(LIB_OBJS)
	ar r $(ALIB_NAME) $(LIB_OBJS)

.c.o:
	$(CC)  -c -fPIC $(LIB_CFLAGS) $(CINCL) $(CFLAGS) -o $@ $<

clean:
	-rm -f $(LIB_OBJS)
