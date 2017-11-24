include ../../../../make.inc
include ../make.inc
#/*================================================================
#   NumX - Linux/MacOsx build             
#
#   Module - Modeling
#
#===================================================================*/

# library name
LIB_NAME = libnumx-modeling-java-1.0

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
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_stat -lnumx_rnd -lnumx_optim -lnumx_data -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(MATH_LIBS) $(GCC_LIBS)
    LIB_SRCS := aptcstjni.c bayessteinjni.c calaptjni.c cokurtjni.c corexpjni.c \
                corlackjni.c    cormjni.c   coskewjni.c covexpjni.c covlackjni.c    \
                covljni.c   covmjni.c   hlretmjni.c implretjni.c    jamessteinjni.c \
                logrlackjni.c  mhoretjni.c  navemjni.c  retlackjni.c    \
                schurpcajni.c   shrinkagecovjni.c   testmodeljni.c  \
                $(PROTECT_PATH)/sysinfo-macosx.c    $(UTILS_PATH)/utils.c
endif
ifeq ($(PLAT),linux)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_stat -lnumx_rnd -lnumx_optim -lnumx_data -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(MATH_LIBS) $(GCC_LIBS)
    LIB_SRCS := aptcstjni.c bayessteinjni.c calaptjni.c cokurtjni.c corexpjni.c \
                corlackjni.c    cormjni.c   coskewjni.c covexpjni.c covlackjni.c    \
                covljni.c   covmjni.c   hlretmjni.c implretjni.c    jamessteinjni.c \
                logrlackjni.c  mhoretjni.c  navemjni.c  retlackjni.c    \
                schurpcajni.c   shrinkagecovjni.c   testmodeljni.c  \
                $(UTILS_PATH)/utils.c
#                $(PROTECT_PATH)/decrypt-linux.c	$(PROTECT_PATH)/sysinfo-linux.c	\

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
