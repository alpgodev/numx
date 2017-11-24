include ../../../../make.inc
include ../make.inc
#/*================================================================
#   NumX - Linux/MacOsx build             
#
#   Module - Linear Algebra
#
#===================================================================*/

# library name
LIB_NAME   = libnumx-linear-algebra-java-1.0

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
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_alglin -lnumx_calelm -lnumx_rnd $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(MATH_LIBS) $(GCC_LIBS)
    LIB_SRCS := ndmjni.c    nm1jni.c    nminfjni.c  nmjni.c \
                pmcjni.c    pmjni.c pmmtjni.c   pmtmjni.c   \
                pmtvjni.c   pmvjni.c    pmxjni.c    pvxjni.c    \
                rchojni.c   rmjni.c schursjni.c   smjni.c \
                testlinalgjni.c tmjni.c xmjni.c xvjni.c \
                $(PROTECT_PATH)/sysinfo-macosx.c    $(UTILS_PATH)/utils.c
endif
ifeq ($(PLAT),linux)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_alglin -lnumx_calelm -lnumx_rnd $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(MATH_LIBS) $(GCC_LIBS)
    LIB_SRCS := ndmjni.c    nm1jni.c    nminfjni.c  nmjni.c \
                pmcjni.c    pmjni.c pmmtjni.c   pmtmjni.c   \
                pmtvjni.c   pmvjni.c    pmxjni.c    pvxjni.c    \
                rchojni.c   rmjni.c schursjni.c   smjni.c \
                testlinalgjni.c tmjni.c xmjni.c xvjni.c \
                $(UTILS_PATH)/utils.c
#                $(PROTECT_PATH)/decrypt-linux.c	$(PROTECT_PATH)/sysinfo-linux.c	\
  
endif

LIB_OBJS  = $(LIB_SRCS:.c=.o)

# Making Libraries (this line is for static Lib)
$(SOLIB_NAME):  $(LIB_OBJS)
	$(SOLIBTOOL) $(SOLIB_LDFLAGS) -o $(SOLIB_NAME) $(LIB_OBJS) $(SOLIB_LDLIBS)
#	$(SOLIB_LDLIBS) -L/usr/lib/gcc/x86_64-redhat-linux/3.4.3/ -lg2c -lgcc -lm

$(ALIB_NAME): $(LIB_OBJS)
	ar r $(ALIB_NAME) $(LIB_OBJS)

.c.o:
	$(CC)  -c -fPIC $(LIB_CFLAGS) $(CINCL) $(CFLAGS) -o $@ $<

clean:
	-rm -f $(LIB_OBJS)
