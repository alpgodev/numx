include ../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create a library for the MONTECARLO.
#
#  Alternatively, the commands
#       make clean all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in MCLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(MCLIB)
 
MC = GWMC.o MCBASE.o MCBASEM.o MCBCOR.o MCBCORM.o MCBHIS.o	\
	MCBHISM.o MCFACT.o MCFACTM.o MCGEN.o RETMU.o UTMC.o 
$(MC): $(FRC)

ALLOBJ=$(MC)

$(MCLIB): $(ALLOBJ)
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
