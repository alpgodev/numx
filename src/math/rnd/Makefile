include ../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create a library for the RANDOM.
#
#  Alternatively, the commands
#       make clean all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in RANDOMLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(RANDOMLIB)
 
#---------------------------------------------------------
#  process/
#---------------------------------------------------------

PROCESS = process/pearch.o process/pebm.o process/pegarch.o	\
        process/pegbm.o process/pegbmx.o process/pemrbm.o	\
        process/pemrjbm.o process/pesbm.o
$(PROCESS): $(FRC)

#---------------------------------------------------------
#  all/
#---------------------------------------------------------

ALL =   algdiv.o alngam.o alnrel.o apser.o basym.o bcorr.o \
        betaln.o bfrac.o bgrat.o bpser.o bratio.o brcmp1.o	\
        brcomp.o bup.o cdfbet.o cdfbin.o cdfchi.o cdfchn.o	\
        cdff.o cdffnc.o cdfgam.o cdfnbn.o cdfnor.o cdfpoi.o	\
        cdft.o cumbet.o cumbin.o cumchi.o cumchn.o cumf.o	\
        cumfnc.o cumgam.o cumnbn.o cumnor.o cumpoi.o cumt.o	\
        dbetrm.o devlpl.o dexpm1.o dinvnr.o dinvr.o dlanor.o	\
        dln1mx.o dln1px.o dlnbet.o dlngam.o dstrem.o dt1.o	\
        dzror.o erf.o erfc1.o esum.o exparg.o fpser.o gam1.o	\
        gaminv.o gamln.o gamln1.o gamma.o genbet.o genchi.o	\
        genexp.o genf.o gengam.o genmn.o genmul.o gennch.o	\
        gennf.o gennor.o genprm.o genunf.o  getMNormalPdf.o	\
        getNormalPdf.o grat1.o gratio.o gsumln.o ignbin.o	\
        ignnbn.o ignpoi.o ipmpar.o lennob.o phrtsd.o psi.o	\
        ranfl.o rcomp.o rexp.o rlog.o rlog1.o sdot.o setgmn.o	\
        sexpo.o  sgamma.o snorm.o spmpar.o spofa.o stvaln.o	\
        init_seed_linux.o
$(ALL): $(FRC)

ALLOBJ=$(PROCESS) $(ALL)

$(RANDOMLIB): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o process/*.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@

.F.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@
