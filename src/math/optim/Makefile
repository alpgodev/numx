include ../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create a library for the OPTIM.
#
#  Alternatively, the commands
#       make clean all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in OPTIMLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(OPTIMLIB)
 
#---------------------------------------------------------
#  bfgs/
#---------------------------------------------------------

BFGS =  bfgs/bfgsbox.o bfgs/bfgsboxs.o  bfgs/gwbfgs.o   bfgs/qnbfgs.o
$(BFGS): $(FRC)

#---------------------------------------------------------
#  nls/
#---------------------------------------------------------

NLS =   nls/gwnls.o nls/nls.o nls/nlsb.o nls/nlsbs.o nls/nlsbx.o	\
        nls/nlss.o nls/nlssx.o nls/nlsx.o nls/rnls.o nls/rnlsb.o	\
        nls/rnlsbs.o nls/rnlsbx.o nls/rnlss.o nls/rnlssx.o	\
        nls/rnlsx.o nls/utnls.o
$(NLS): $(FRC)

#---------------------------------------------------------
#  qp/
#---------------------------------------------------------

QP = qp/buildfeas.o qp/gwqp.o qp/initfeas.o qp/qp.o
$(QP): $(FRC)

#---------------------------------------------------------
#  sdls/
#---------------------------------------------------------

SDLS =  sdls/bfgssdls.o sdls/checksdls.o    sdls/checksdlseq.o  \
        sdls/checksdlsineq.o    sdls/gwsdls.o   sdls/sdls.o \
        sdls/sdlsce.o   sdls/sdlsci.o   sdls/sdlscor.o  \
        sdls/sdlsg.o    sdls/sdlsgen.o  sdls/sdlss.o    \
        sdls/sdlssg.o   sdls/sdlstrace.o    sdls/simscor.o  \
        sdls/simsdls.o  sdls/simsdlsg.o sdls/simstrace.o    \
        sdls/solscor.o  sdls/solsdls.o  sdls/solsdlsg.o \
        sdls/solstrace.o    sdls/utsdls.o   sdls/vcsdls.o
$(SDLS): $(FRC)

#---------------------------------------------------------
#  utils/
#---------------------------------------------------------

UTILS = utils/ajour.o utils/anfm01.o utils/anfm02.o utils/anfm03.o	\
        utils/anfm04.o utils/anfm05.o utils/anfm06.o utils/anrs01.o	utils/anrs02.o	\
        utils/aux003.o utils/auxo01.o utils/bfgsd.o utils/calbx.o utils/calmaj.o	\
        utils/ctcab.o utils/ctonb.o utils/dcube.o utils/ddd2.o utils/desr03.o	\
        utils/dimp03.o utils/dipvtf.o utils/dnrm0.o utils/dogleg.o utils/dpmpar.o	\
        utils/enorm.o utils/f1qhes.o utils/f2qhes.o utils/fdjac1.o utils/fdjac2.o	\
        utils/ffinf1.o utils/fmulb1.o utils/fmuls1.o utils/fpq2.o utils/fprf2.o	\
        utils/frdf1.o utils/fremf1.o utils/fremf2.o utils/gcbd.o utils/gcp.o	\
        utils/hybrd.o utils/hybrd1.o utils/hybrj.o utils/hybrj1.o utils/icscof.o	\
        utils/icse.o utils/icse0.o utils/icse1.o utils/icse2.o utils/icsec2.o	\
        utils/icsei.o utils/lmder.o utils/lmdif.o utils/lmpar.o utils/majour.o	\
        utils/majysa.o utils/majz.o utils/n1fc1.o utils/n1fc1a.o utils/n1fc1o.o	\
        utils/n1gc2.o utils/n1gc2a.o utils/n1gc2b.o utils/n1qn1.o utils/n1qn1a.o	\
        utils/n1qn2.o utils/n1qn2a.o utils/n1qn3.o utils/n1qn3a.o utils/nlis2.o	\
        utils/nqhess.o utils/nvkt03.o utils/optr01.o utils/optr03.o utils/opvf03.o	\
        utils/pasr03.o utils/plcbas.o utils/proj.o utils/qform.o utils/qnbd.o	\
        utils/qrfac.o utils/qrsolv.o utils/r1mpyq.o utils/r1updt.o utils/rednor.o	\
        utils/relvar.o utils/rlbd.o utils/satur.o utils/shanph.o utils/simext.o	\
        utils/strang.o utils/tol03.o utils/zgcbd.o utils/zqnbd.o utils/zthz.o	\
        utils/n2qn1a.o utils/fmc11z.o utils/fajc1.o utils/fretc1.o utils/fmani1.o	\
        utils/fcomp1.o utils/fmlag1.o utils/nlis0.o utils/fcube.o utils/fuclid.o	\
        utils/fmc11a.o utils/fmc11b.o utils/fmc11e.o
$(UTILS): $(FRC)

ALLOBJ=$(BFGS) $(NLS) $(QP) $(SDLS) $(UTILS)

$(OPTIMLIB): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

FRC:
	@FRC=$(FRC)

clean:
	rm -f bfgs/*.o nls/*.o qp/*.o sdls/*.o utils/*.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@

.F.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@
