include ../../../../make.inc
include ../make.inc
#/*================================================================
#   NumX - Linux/MacOsx build             
#
#   Module - Reporting
#
#===================================================================*/

# library name
LIB_NAME   = libnumx-reporting-java-1.0

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
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_stat -lnumx_rnd -lnumx_optim -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(GCC_LIBS) $(MATH_LIBS)
    LIB_SRCS := diversificationjni.c  downsidejni.c exaccjni.c    exacvarjni.c  exaipcjni.c   \
                examshjni.c   examvajni.c   exanshortfalljni.c    exanvajni.c   exaretjni.c   \
                exarirjni.c   exarkujni.c   exarskjni.c   exarvojni.c   exasrajni.c   exaterjni.c   \
                exavrijni.c   expcoefcorrjni.c  expcvarjni.c  expirajni.c   expkurjni.c   \
                expmaxlossjni.c   expmshjni.c   expmvajni.c   expnpcvarjni.c    expnpvarjni.c \
                expnshortfalljni.c    expnvajni.c   expterjni.c   expskejni.c   expsorjni.c   \
                expsrajni.c   expretjni.c   expvoljni.c   expvolexpjni.c    expvrijni.c   \
                futvaljni.c   kerneldjni.c  pastvaljni.c  retcontribjni.c   rolpiesjni.c  \
                semvoljni.c   tecontribjni.c    testanajni.c  varcontribjni.c   volcontribjni.c   \
                $(PROTECT_PATH)/sysinfo-macosx.c    $(UTILS_PATH)/utils.c
endif
ifeq ($(PLAT),linux)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_stat -lnumx_rnd -lnumx_optim -lnumx_alglin -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(GCC_LIBS) $(MATH_LIBS)
    LIB_SRCS := diversificationjni.c  downsidejni.c exaccjni.c    exacvarjni.c  exaipcjni.c   \
                examshjni.c   examvajni.c   exanshortfalljni.c    exanvajni.c   exaretjni.c   \
                exarirjni.c   exarkujni.c   exarskjni.c   exarvojni.c   exasrajni.c   exaterjni.c   \
                exavrijni.c   expcoefcorrjni.c  expcvarjni.c  expirajni.c   expkurjni.c   \
                expmaxlossjni.c   expmshjni.c   expmvajni.c   expnpcvarjni.c    expnpvarjni.c \
                expnshortfalljni.c    expnvajni.c   expterjni.c   expskejni.c   expsorjni.c   \
                expsrajni.c   expretjni.c   expvoljni.c   expvolexpjni.c    expvrijni.c   \
                futvaljni.c   kerneldjni.c  pastvaljni.c  retcontribjni.c   rolpiesjni.c  \
                semvoljni.c   tecontribjni.c    testanajni.c  varcontribjni.c   volcontribjni.c   \
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
