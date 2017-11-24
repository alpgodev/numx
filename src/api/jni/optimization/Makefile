include ../../../../make.inc
include ../make.inc
#/*================================================================
#   NumX - Linux/MacOsx build             
#
#===================================================================*/

# library name
LIB_NAME   = libnumx-optimization-java-1.0

# output libraries directory
ifeq ($(PLAT),macosx)
    SOLIB_NAME = $(LIB_PATH)/$(LIB_NAME).dylib
endif
ifeq ($(PLAT),linux)
    SOLIB_NAME = $(LIB_PATH)/$(LIB_NAME).so
endif
ALIB_NAME  = $(LIB_PATH)/$(LIB_NAME).a

all: $(SOLIB_NAME) $(ALIB_NAME)

# any extra libraries that your librarie may depend on.
ifeq ($(PLAT),macosx)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_optim_portfolio -lnumx_stat -lnumx_optim -lnumx_rnd -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(MATH_LIBS) $(GCC_LIBS)
    LIB_SRCS := allocitjni.c    allocitrfrjni.c allocmvjni.c    allocmvrfrjni.c \
                allocmvtcjni.c  allocrbjni.c    allocrbrfrjni.c allocsrjni.c    \
                allocsrrfrjni.c calepsvarjni.c  covfiltering3jni.c  covfilteringjni.c   \
                cstsecjni.c     mvejni.c    qpjni.c riskbudgetitjni.c   sdlscijni.c     \
                sdlscorjni.c    sdlsgenjni.c    sdlstracejni.c  testOptimizationjni.c   \
                $(PROTECT_PATH)/sysinfo-macosx.c    $(UTILS_PATH)/utils.c
endif
ifeq ($(PLAT),linux)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_optim_portfolio -lnumx_stat -lnumx_optim -lnumx_rnd -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) -lgfortran
    LIB_SRCS := allocitjni.c    allocitrfrjni.c allocmvjni.c    allocmvrfrjni.c \
                allocmvtcjni.c  allocrbjni.c    allocrbrfrjni.c allocsrjni.c    \
                allocsrrfrjni.c calepsvarjni.c  covfiltering3jni.c  covfilteringjni.c   \
                cstsecjni.c     mvejni.c    qpjni.c riskbudgetitjni.c   sdlscijni.c     \
                sdlscorjni.c    sdlsgenjni.c    sdlstracejni.c  testOptimizationjni.c   \
                $(UTILS_PATH)/utils.c
#                $(PROTECT_PATH)/decrypt-linux.c	$(PROTECT_PATH)/sysinfo-linux.c	
                
endif

LIB_OBJS  = $(LIB_SRCS:.c=.o)

# Making Libraries (this line is for static Lib)
$(SOLIB_NAME):  $(LIB_OBJS)
	$(SOLIBTOOL) $(SOLIB_LDFLAGS) -o $(SOLIB_NAME) $(LIB_OBJS) $(SOLIB_LDLIBS)

$(ALIB_NAME): $(LIB_OBJS)
	ar r $(ALIB_NAME) $(LIB_OBJS)

.c.o:
	$(CC)  -c -fPIC $(LIB_CFLAGS) $(CINCL) $(CFLAGS) -o $@ $<

clean:
	-rm -f $(LIB_OBJS)
