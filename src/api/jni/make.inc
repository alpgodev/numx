####################################################################
#  make include file.                                              #
####################################################################

# library path
LIB_PATH = ../../../../lib/api/$(PLAT)

# utils path
UTILS_PATH = ../../utils

# protect path
PROTECT_PATH = ../../protect

# include path
INCLUDE_PATH = ../../include

# PROTECT:
# add -DRELEASE_C
# add -DNOCHECKMAC to PROTECT to disable Mac Address checking
PROTECT = #-DRELEASE_C

DEBUG = -O3 -funroll-all-loops
CFLAGS = $(PROTECT) $(DEBUG)
LDFLAGS =

ifeq ($(PLAT),macosx)
    CINCL = -I. -I$(JDK_ROOT)/include	\
            -I$(INCLUDE_PATH)/protect	-I$(INCLUDE_PATH)/jni/macosx    \
            -I$(INCLUDE_PATH)/c    -I$(INCLUDE_PATH)/jni   \
            -I$(INCLUDE_PATH)/utils -I$(INCLUDE_PATH)/
    SOLIB_LDFLAGS = -fPIC -dynamiclib -o $(SOLIB_NAME) --verbose
endif
ifeq ($(PLAT),linux)
    CINCL = -I. -I$(JDK_ROOT)/include/linux -I$(JDK_ROOT)/include	\
            -I$(INCLUDE_PATH)/protect -I$(INCLUDE_PATH)/jni   -I$(INCLUDE_PATH)/c \
            -I$(INCLUDE_PATH)/utils
    SOLIB_LDFLAGS = -fPIC -shared -Wl,-soname,$(LIB_NAME).so --verbose
endif

AS = as
SOLIBTOOL = gcc


