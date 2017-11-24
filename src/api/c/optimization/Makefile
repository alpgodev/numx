include ../../../../make.inc
include ../make.inc
#================================================================
#   NumX - C/C++ Linux/MacOsx build             
#
#   Module - Optimization
#
#================================================================

# Library name
LIB_NAME = libnumx-optimization-c-1.0

# output libraries directory
SOLIB_NAME = $(LIB_PATH)/$(LIB_NAME).so
ALIB_NAME  = $(LIB_PATH)/$(LIB_NAME).a

all: $(SOLIB_NAME) $(ALIB_NAME)

# any extra libraries that your librarie may depend on.
ifeq ($(PLAT),macosx)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_stat -lnumx_optim -lnumx_optim_portfolio -lnumx_rnd -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB)  $(MATH_LIBS) $(GCC_LIBS)
    LIB_SRCS := allocitc.c  allocitrfrc.c   allocmvc.c  allocmvrfrc.c   \
                allocmvtcc.c    allocrbc.c  allocrbrfrc.c   allocsrc.c  \
                allocsrrfrc.c   bfgsboxc.c  bfgsboxsc.c calepsvarc.c    \
                covfiltering3c.c    covfilteringc.c cstsecc.c   mvec.c  \
                nlsbc.c nlsc.c  nlsxc.c qpc.c rnlsbc.c rnlsc.c \
                riskbudgetitc.c sdlscic.c   sdlscorc.c  sdlsgenc.c  \
                sdlstracec.c    testOptimizationc.c \
                $(PROTECT_PATH)/sysinfo-macosx.c    $(UTILS_PATH)/utils.c
endif
ifeq ($(PLAT),linux)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_optim_portfolio -lnumx_stat -lnumx_optim -lnumx_rnd -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(MATH_LIBS) $(GCC_LIBS)
    LIB_SRCS := allocitc.c  allocitrfrc.c   allocmvc.c  allocmvrfrc.c   \
                allocmvtcc.c    allocrbc.c  allocrbrfrc.c   allocsrc.c  \
                allocsrrfrc.c   bfgsboxc.c  bfgsboxsc.c calepsvarc.c    \
                covfiltering3c.c    covfilteringc.c cstsecc.c   mvec.c  \
                nlsbc.c nlsc.c  nlsxc.c qpc.c rnlsbc.c rnlsc.c \
                riskbudgetitc.c sdlscic.c   sdlscorc.c  sdlsgenc.c  \
                sdlstracec.c    testOptimizationc.c \
#                $(PROTECT_PATH)/decrypt-linux.c	$(PROTECT_PATH)/sysinfo-linux.c \
                $(UTILS_PATH)/utils.c
endif

LIB_OBJS  = $(LIB_SRCS:.c=.o)

# Making Libraries (his line is for static Lib)
$(SOLIB_NAME):  $(LIB_OBJS)
	$(SOLIBTOOL) $(SOLIB_LDFLAGS) -o $(SOLIB_NAME) $(LIB_OBJS) $(SOLIB_LDLIBS)

$(ALIB_NAME): $(LIB_OBJS)
	ar r $(ALIB_NAME) $(LIB_OBJS)

.c.o:
	$(CC)  -c -fPIC $(LIB_CFLAGS) $(CINCL) $(CFLAGS) -o $@ $<

clean:
	-rm -f $(LIB_OBJS)
