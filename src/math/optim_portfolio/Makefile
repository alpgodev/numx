include ../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create a library for Portfolio Optimization
#
#  Alternatively, the commands
#       make clean all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in ALLOCLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(ALLOCLIB)
 
#---------------------------------------------------------
#  utils
#---------------------------------------------------------

UALLOC = allocit.o  allocitrfr.o    allocmdd.o  allocmv.o   \
        allocmvrfr.o    allocmvtc.o allocrb.o   allocrbrfr.o    \
        allocsr.o   allocsrrfr.o    alloctest.o alloctestrfr.o  \
        allocvar.o  bfgsmvol.o  buildrbit.o checkfeasit.o       \
        checklinbox.o   cstsec.o    gwalloc.o   gwmvol.o        \
        multivolrfr.o   mve.o   omegab.o    opit.o  opitrfr.o   \
        opmv.o  opmvrfr.o   opvar.o riskbudgetit.o      \
        simmvolrfr.o    simrb.o     testlinearcst.o utalloc.o

$(UALLOC): $(FRC)

ALLOBJ=$(UALLOC)

$(ALLOCLIB): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@

.F.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@
