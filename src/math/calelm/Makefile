include ../../../make.inc
include make.inc

#######################################################################
#  This is the makefile to create a library for the CALELM.
#
#  Alternatively, the commands
#       make all
#  or
#       make
#  without any arguments creates a library.
#  The name of the library is held in CALELMLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make FRC=FRC
#
#######################################################################

all: clean $(CALELMLIB)
 
#---------------------------------------------------------
#  utils
#---------------------------------------------------------

UCALELM = arcosh.o calerf.o ccopy.o coshin.o cshep2d.o cupro.o	\
        cuproi.o cusum.o dad.o dadd.o dcabs1.o ddif.o ddpow.o	\
        ddpow1.o ddpowe.o ddrdiv.o dipow.o dipowe.o dlblks.o	\
        dlgama.o dmcopy.o dmmul.o dmmul1.o dmprod.o dmsum.o	\
        drdiv.o dsearch.o dset.o dsort.o dspfit.o dsum.o	\
        dtensbs.o dtild.o dvmul.o dwdiv.o dwpow.o dwpow1.o	\
        dwpowe.o dwrdiv.o dxlegf.o entier.o franck.o gamma.o	\
        gdcp2i.o hilber.o icopy.o idmax.o idmin.o imcopy.o	\
        infinity.o int2db.o intp.o iset.o isort.o isova0.o	\
        isoval.o ivimp.o iwamax.o kronc.o kronr.o lnblnk.o	\
        magic.o mtran.o nearfloat.o pchim.o psi.o pythag.o	\
        rat.o	rcopy.o rcsort.o rea2db.o ribesl.o rjbesl.o	\
        rkbesl.o round.o rybesl.o simple.o somespline.o	\
        unsfdcopy.o urand.o vpythag.o wacos.o wasin.o wasum.o	\
        watan.o waxpy.o wcopy.o wddiv.o wdiv.o wdotci.o	\
        wdotcr.o wdpow.o wdpow1.o wdpowe.o wdrdiv.o wipow.o	\
        wipowe.o wlog.o wmmul.o wmprod.o wmsum.o wmul.o	\
        wnrm2.o wrscal.o wscal.o wsign.o wsort.o wsqrt.o	\
        wswap.o wtan.o wvmul.o wwdiv.o wwpow.o wwpow1.o	\
        wwpowe.o wwrdiv.o
$(UCALELM): $(FRC)

ALLOBJ=$(UCALELM)

$(CALELMLIB): $(ALLOBJ)
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
