include ../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create a library for the CLUSTER.
#
#  Alternatively, the commands
#       make clean all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in CLUSTERLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(CLUSTERLIB)
 
CLUSTER = 	cluster.o em.o estep.o gwcluster.o mstep.o
$(CLUSTER): $(FRC)

ALLOBJ=$(CLUSTER)

$(CLUSTERLIB): $(ALLOBJ)
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
