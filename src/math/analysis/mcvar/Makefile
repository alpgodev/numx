include ../../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create a library for the Value-at-Risk module
#
#  Alternatively, the commands
#       make clean all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in MCVARLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#---------------------------------------------------------------------
#  Yann Vernaz
#######################################################################

all: clean $(MCVARLIB)
 
#---------------------------------------------------------
#  Monte-Carlo Value-at-Risk sources
#---------------------------------------------------------

VAR = mcPortfolioVaR.o mcFactor.o mcFactorBlock.o mcVaR.o mcBlockVaR.o mcBlockXVaR.o
$(VAR): $(FRC)

ALLOBJ=$(VAR) 

$(MCVARLIB): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

.c.o: 
	$(CC) $(OPTS) -I$(INCLUDE) -I$(INCLUDE_UTILS) -c $< -o $@
