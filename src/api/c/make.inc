####################################################################
#  make include file.                                              #
####################################################################

# library path (output)
LIB_PATH = ../../../../lib/api/$(PLAT)

# utils path
UTILS_PATH = ../../utils

# protect path
PROTECT_PATH = ../../protect

# include path
INCLUDE_PATH = ../../include

# -DTRANSPOSE
# PROTECT:
# add -DRELEASE_C for C and COM this option should be removed for VB
# add -DNOCHECKMAC to PROTECT to disable Mac Address checking
PROTECT = -DTRANSPOSE 

DEBUG = -O3 -funroll-all-loops
CFLAGS = $(PROTECT) $(DEBUG)

ifeq ($(PLAT),macosx)
    CINCL = -I. -I$(JDK_ROOT)/include/linux -I$(JDK_ROOT)/include	\
            -I$(INCLUDE_PATH)/protect	-I$(INCLUDE_PATH)/jni/macosx    \
            -I$(INCLUDE_PATH)/c -I$(INCLUDE_PATH)/utils -I$(INCLUDE_PATH)/
    SOLIB_LDFLAGS = -fPIC -dynamiclib -o $(SOLIB_NAME) --verbose
endif
ifeq ($(PLAT),linux)
    CINCL = -I. -I$(JDK_ROOT)/include/linux -I$(JDK_ROOT)/include	\
            -I$(INCLUDE_PATH)/protect   -I$(INCLUDE_PATH)/c     \
            -I$(INCLUDE_PATH)/utils
    SOLIB_LDFLAGS = -fPIC -shared -Wl,-soname,$(LIB_NAME).so --verbose
endif

AS = as
SOLIBTOOL = gcc

LDFLAGS =
