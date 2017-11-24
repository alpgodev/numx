include ../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create the library DATA.
#
#  Alternatively, the commands
#       make all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in DATALIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(DATALIB)
 
#---------------------------------------------------------
#  
#---------------------------------------------------------

UDATA = gwnavem.o navem.o retem.o
$(UDATA): $(FRC)

ALLOBJ=$(UDATA)

$(DATALIB): $(ALLOBJ)
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

