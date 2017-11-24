include ../../../../make.inc
include ../make.inc
#/*================================================================
#   NumX - C/C++ Linux/MacOsx build             
#
#   Module - Reporting (Performance & Risk Analysis)
#
#===================================================================*/

# Library name
LIB_NAME   = libnumx-reporting-c-1.0

# output libraries directory
SOLIB_NAME = $(LIB_PATH)/$(LIB_NAME).so
ALIB_NAME  = $(LIB_PATH)/$(LIB_NAME).a

all: $(SOLIB_NAME) $(ALIB_NAME)

# any extra libraries that your librarie may depend on.
ifeq ($(PLAT),macosx)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_mcvar -lnumx_stat -lnumx_optim -lnumx_alglin -lnumx_rnd -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(GCC_LIBS) $(MATH_LIBS)
    LIB_SRCS := diversificationc.c  downsidec.c exaccc.c    exacvarc.c  exaipcc.c   \
                examshc.c   examvac.c   exanshortfallc.c    exanvac.c   exaretc.c   \
                exarirc.c   exarkuc.c   exarskc.c   exarvoc.c   exasrac.c   exaterc.c   \
                exavric.c   expcoefcorrc.c  expcvarc.c  expirac.c   expkurc.c   \
                expmaxlossc.c   expmshc.c   expmvac.c   expnpcvarc.c    expnpvarc.c \
                expnshortfallc.c    expnvac.c   expterc.c   expskec.c   expsorc.c   \
                expsrac.c   expretc.c   expvolc.c   expvolexpc.c    expvric.c   \
                futvalc.c   kerneldc.c  pastvalc.c  retcontribc.c   rolpiesc.c  \
                semvolc.c   tecontribc.c    testanac.c  varcontribc.c   volcontribc.c   \
                $(PROTECT_PATH)/sysinfo-macosx.c    $(UTILS_PATH)/utils.c
endif
ifeq ($(PLAT),linux)
    SOLIB_LDLIBS =  -L../../../../lib/math/$(PLAT) -lnumx_mcvar -lnumx_stat -lnumx_optim -lnumx_alglin -lnumx_rnd -lnumx_calelm $(LAPACK_LIB) $(BLAS_LIB) $(GFORTRAN_LIB) $(GCC_LIBS) $(MATH_LIBS)
    LIB_SRCS := diversificationc.c  downsidec.c exaccc.c    exacvarc.c  exaipcc.c   \
                examshc.c   examvac.c   exanshortfallc.c    exanvac.c   exaretc.c   \
                exarirc.c   exarkuc.c   exarskc.c   exarvoc.c   exasrac.c   exaterc.c   \
                exavric.c   expcoefcorrc.c  expcvarc.c  expirac.c   expkurc.c   \
                expmaxlossc.c   expmshc.c   expmvac.c   expnpcvarc.c    expnpvarc.c \
                expnshortfallc.c    expnvac.c   expterc.c   expskec.c   expsorc.c   \
                expsrac.c   expretc.c   expvolc.c   expvolexpc.c    expvric.c   \
                futvalc.c   kerneldc.c  pastvalc.c  retcontribc.c   rolpiesc.c  \
                semvolc.c   tecontribc.c    testanac.c  varcontribc.c   volcontribc.c   \
#                $(PROTECT_PATH)/decrypt-linux.c	$(PROTECT_PATH)/sysinfo-linux.c	\
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
