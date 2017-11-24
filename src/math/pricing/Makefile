include ../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create a library.
#
#  Alternatively, the commands
#       make clean all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in PRICINGLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(PRICINGLIB)

PRICING = PRbinomial.o PRvanilla.o PRmonteCarlo.o PRutils.o 	
$(PRICING): $(FRC)

ALLOBJ=$(PRICING)

$(PRICINGLIB): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@
