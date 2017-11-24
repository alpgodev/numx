include ../../../make.inc
include make.inc
#######################################################################
#  This is the makefile to create the library STATLIB.
#
#  Alternatively, the commands
#       make clean all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in STATLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(STATLIB)
 
#---------------------------------------------------------
#  apt/
#---------------------------------------------------------

APT =	apt/aptcst.o apt/aptvol.o apt/calapt.o apt/gwapt.o apt/utapt.o
$(APT): $(FRC)

#---------------------------------------------------------
#  cov/
#---------------------------------------------------------

COV =   cov/calepsmean.o    cov/calepsmedian.o  cov/chtcovm.o   \
        cov/cokurt.o    cov/corexp.o    cov/corlack.o   cov/corm.o  \
        cov/cormean.o   cov/cormvm.o    cov/corvolm.o   cov/coskew.o    \
        cov/covexp.o    cov/covfiltering.o  cov/covfiltering3.o \
        cov/covl.o  cov/covlack.o   cov/covm.o  cov/covmvm.o    \
        cov/empcov.o    cov/empcovexp.o cov/etahat.o    cov/gwcov.o \
        cov/gwshrinkagecov.o    cov/nuhat.o cov/pihat.o \
        cov/shrinkagecov.o cov/shrinkageintensity.o cov/shrinkagetarget.o
$(COV): $(FRC)

#---------------------------------------------------------
#  kernel/
#---------------------------------------------------------

KERNEL = kernel/glkern.o kernel/gwnpreg.o kernel/kernela.o  \
        kernel/kerneld.o    kernel/lokern.o kernel/npreg.o  \
        kernel/plugin.o kernel/subs.o
$(KERNEL): $(FRC)

#---------------------------------------------------------
#  perf/
#---------------------------------------------------------

PERF =  perf/bayesstein.o   perf/gwperf.o   perf/implret.o	\
        perf/jamesstein.o
$(PERF): $(FRC)

#---------------------------------------------------------
#  return/
#---------------------------------------------------------

RETURN = return/hlprim.o return/hlretm.o return/logrlack.o return/mhoret.o	\
        return/retlack.o return/vhoret.o
$(RETURN): $(FRC)

#---------------------------------------------------------
#  risk/
#---------------------------------------------------------

RISK = risk/ARETURN.o risk/DOWNSIDE.o risk/EXACVAR.o risk/EXAMSH.o risk/EXAMVA.o	\
        risk/EXANSHORTFALL.o risk/EXANVA.o risk/EXARET.o	\
        risk/EXARIR.o risk/EXARKU.o risk/EXARSK.o risk/EXARVA.o risk/EXARVO.o	\
        risk/EXASRA.o risk/EXASTARR.o	\
        risk/EXATER.o risk/EXAVRI.o risk/EXPCORRCOEF.o risk/EXPCVAR.o	\
        risk/EXPDNB.o risk/EXPDOWNC.o risk/EXPDPERCENT.o risk/EXPGMEAN.o	\
        risk/EXPGPERCENT.o risk/EXPIRA.o risk/EXPKUR.o risk/EXPLMEAN.o	\
        risk/EXPMAXLOSS.o risk/EXPMSH.o risk/EXPMVA.o risk/EXPNPVAR.o	\
        risk/EXPNSHORTFALL.o risk/EXPNVA.o risk/EXPOMEGA.o risk/EXPRET.o	\
        risk/EXPSKE.o risk/EXPSOR.o risk/EXPSRA.o risk/EXPSTARR.o	\
        risk/EXPTER.o risk/EXPUPC.o	\
        risk/EXPUPNB.o risk/EXPUPPERCENT.o risk/EXPVAR.o risk/EXPVOL.o	\
        risk/EXPVOLEXP.o risk/EXPVRI.o risk/FUTVAL.o risk/GWRISK.o risk/KURTOS.o	\
        risk/MEDIAN.o risk/MODSHA.o risk/MODVAR.o risk/MOMENT.o risk/NORMINV.o	\
        risk/PASTVAL.o risk/RETCONTRIB.o risk/ROLLINGPIES.o risk/ROLPIES.o	\
        risk/SDVHST.o risk/SDVRISK.o risk/SKEWN.o risk/VARCONTRIB.o	\
        risk/TECONTRIB.o risk/TENSORK.o risk/TENSORKK.o risk/MOMENT3.o	\
        risk/EXACC.o risk/EXAIPC.o risk/DIVERSIFICATION.o	\
        risk/VOLCONTRIB.o risk/EXPNPCVAR.o

$(RISK): $(FRC)

ALLOBJ= $(RISK) $(APT) $(COV) $(KERNEL) $(PERF) $(RETURN)

$(STATLIB): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

FRC:
	@FRC=$(FRC)

clean:
	rm -f apt/*.o cov/*.o kernel/*.o perf/*.o return/*.o risk/*.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@

.F.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@
