include ../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create a library for the ALGLIN.
#
#  Alternatively, the commands
#       make clean all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in ALGLINLIB,
#  which is set in the make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(ALGLINLIB)
 
#---------------------------------------------------------
#  utils
#---------------------------------------------------------

UALGLIN = utils/addmcz.o utils/addmlz.o utils/blas1.o utils/blas2.o utils/blas3.o	\
        utils/eigenvd.o utils/gweigenvd.o   utils/gwjms.o   utils/svd.o \
        utils/testsdp.o utils/utdiv.o   utils/utio.o    utils/utioi.o   \
        utils/utmat.o   utils/utmati.o  utils/utrand.o  utils/utrandil.o    \
        utils/utrandl.o utils/utstat.o

$(UALGLIN): $(FRC)

#---------------------------------------------------------
#  cmc
#---------------------------------------------------------

CMCALGLIN = cmc/calepsvar.o cmc/chol.o  cmc/chtmcm.o    cmc/greps.o cmc/gwcmc.o     \
            cmc/kat.o   cmc/katave.o    cmc/katavep.o   cmc/kateps.o    \
            cmc/katoad.o    cmc/katoblock.o cmc/katoblockx.o    cmc/katp.o  \
            cmc/katvar.o    cmc/pca.o   cmc/project.o cmc/projecta.o    \
            cmc/projects.o  cmc/rcho.o  cmc/schur.o cmc/schuro.o    cmc/schurpca.o  \
            cmc/schurs.o    cmc/schurso.o

$(CMCALGLIN): $(FRC)

ALLOBJ=$(UALGLIN) $(CMCALGLIN)

$(ALGLINLIB): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

FRC:
	@FRC=$(FRC)

clean:
	rm -f utils/*.o cmc/*.o

.f.o:
	$(FORTRAN) $(OPTS) -c $< -o $@

.F.o:
	$(FORTRAN) $(OPTS) -c $< -o $@

