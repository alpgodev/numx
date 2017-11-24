include ../../../../make.inc
include ../make.inc
#/*================================================================
#   NumX - Linux/MacOsx build             
#
#   Module - Simulation
#
#===================================================================*/

# library name
LIB_NAME = libnumx-simulation-java-1.0

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
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_simul -lnumx_stat -lnumx_optim -lnumx_rnd -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(LINUX_MATH_LIBS) $(LINUX_GCC_LIBS)
    LIB_SRCS := bmnjni.c    \
                pebmjni.c   pegarchjni.c    pegbmjni.c  pegbmxjni.c \
                pemrbmjni.c pesbmjni.c  testsimuljni.c  \
                $(PROTECT_PATH)/sysinfo-macosx.c    $(UTILS_PATH)/utils.c
endif
ifeq ($(PLAT),linux)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_simul -lnumx_stat -lnumx_optim -lnumx_rnd -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(LINUX_MATH_LIBS) $(LINUX_GCC_LIBS)
    LIB_SRCS := bmnjni.c    \
                pebmjni.c   pegarchjni.c    pegbmjni.c  pegbmxjni.c \
                pemrbmjni.c pesbmjni.c  testsimuljni.c  \
                $(UTILS_PATH)/utils.c
#                $(PROTECT_PATH)/decrypt-linux.c	$(PROTECT_PATH)/sysinfo-linux.c	\
               
endif

LIB_OBJS  = $(LIB_SRCS:.c=.o)

# Making Libraries
$(SOLIB_NAME):  $(LIB_OBJS)
	$(SOLIBTOOL) $(SOLIB_LDFLAGS) -o $(SOLIB_NAME) $(LIB_OBJS) $(SOLIB_LDLIBS)

$(ALIB_NAME): $(LIB_OBJS)
	ar r $(ALIB_NAME) $(LIB_OBJS)

.c.o:
	$(CC)  -c -fPIC $(LIB_CFLAGS) $(CINCL) $(CFLAGS) -o $@ $<

clean:
	-rm -f $(LIB_OBJS)
