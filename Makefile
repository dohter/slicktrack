#
# Makefile for SLICKTRACK
#
#To get optimization so that Not-a-Number etc doesn't cause a crash.
#FFLAGS= -xtarget=native -O5 -libmil -fsimple=2 -dalign -xlibmopt -depend -fns -ftrap=common -pad=local -xvector=yes -xprefetch=yes -xprefetch_level=2
#FFLAGS= -O3 -libmil -fsimple=1 -dalign -xlibmopt -depend -ftrap=common -xvector=no -xprefetch=no
#No optimization so that Not-a-Number etc does cause a crash.
#FFLAGS= -g -C
#FFLAGS= -g -C +T
#-L for listing
#
FFLAGS= -O3
#FFLAGS= -O3  -m32
#FFLAGS= -O3 -m32
# -ftrace=full
#
#
# libraries needed# need option lF77 to link the F77 naglib with an F90 prgm#
#  SNAG = -lnag  -lF77 #
#  SNAG = -lnag  -L. -lF77
#  SNAG = -L/opt/products/naglib/20 -lnag -L/opt/products/gcc/3.3.3/lib -lg2c
#   PNAG = -L/afs/desy.de/i586_rhel40/products/naglib/20/ -lnag
#   RHLX = -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/32 -lg2c -m32
#   PNAG = /opt/NAG/fll3a22dfl/lib/libnag_nag.a
#    PNAG = /Users/fanglei/NAG/flmi622dcl/lib/libnag_nag.a
    PNAG = /Users/fanglei/NAG/flmi6261fl/lib/libnag_nag.a
#
# sources of SLICK_F90.MAY97
#
SRCslick=slick.f aver.f averv5sq.f avorb.f lbeambeam.f lbeambeam9.f cavity.f cfdiph.f cfdipv.f chrotn.f damint.f damper.f date.f  dcgmpr.f dertop.f diaver.f diph.f dipv.f dx66.f dx88.f  emitnc.f fixorb1.f fixorb2.f fixorb3.f fixorb4.f fixorb5.f hedge.f iparam.f jam333.f jam444.f jam666.f jam777.f jam888.f jam881.f jam999.f lattin.f linopt.f linopenopt.f mult.f mx66.f mx66damp.f mx88.f mx88damp.f mx99damp.f mxdamp.f nlbeambeam.f nlbeambeam9.f norm.f norm2.f orbit.f ortchk.f pertop.f psequil.f quartm.f renorm.f reson.f roquad.f rotatr.f rspin.f scrurita1.f scrurita2.f scrurita3.f scrurita4.f setbb.f sex88.f sexeq.f sig.f simq.f skewquad.f snake.f sol6.f sol66.f sol77.f sol8an.f sol9an.f solxyp.f spin.f ssect.f symp6.f symp8.f thiquad.f tspin.f unit.f unitm.f \
f03aaf.f f02ebf.f \
g05skf.f \
g05sqf.f \
g05kff.f \
g05rzf.f \
s15ddf.f \
rnglib.f ranlib.f \
string_trim.f90

# object modules of SLICK_F90.MAY97
#
OBJslick=slick.o aver.o averv5sq.o avorb.o lbeambeam.o lbeambeam9.o cavity.o cfdiph.o cfdipv.o chrotn.o damint.o damper.o date.o  dcgmpr.o dertop.o diaver.o diph.o dipv.o dx66.o dx88.o  emitnc.o fixorb1.o fixorb2.o fixorb3.o fixorb4.o fixorb5.o hedge.o iparam.o jam333.o jam444.o jam666.o jam777.o jam888.o jam881.o jam999.o lattin.o linopt.o linopenopt.o mult.o mx66.o mx66damp.o mx88.o mx88damp.o mx99damp.o mxdamp.o nlbeambeam.o nlbeambeam9.o norm.o norm2.o orbit.o ortchk.o pertop.o psequil.o quartm.o renorm.o reson.o roquad.o rotatr.o rspin.o scrurita1.o scrurita2.o scrurita3.o scrurita4.o setbb.o sex88.o sexeq.o sig.o simq.o skewquad.o snake.o sol6.o sol66.o sol77.o sol8an.o sol9an.o solxyp.o spin.o ssect.o symp6.o symp8.o thiquad.o tspin.o unit.o unitm.o \
f02ebf.o \
g05kff.o \
g05skf.o \
g05sqf.o \
s15ddf.o \
f03aaf.o \
g05rzf.o \
rnglib.o ranlib.o \
string_trim.o

%.o: %.f
#       gfortran also works here
#	gfortran -c $(FFLAGS) $<
#	f95 -c $(FFLAGS) $<
	gfortran -c $(FFLAGS) $<

%.o: %.f90
	gfortran -c $(FFLAGS) $<

slick:  $(OBJslick)
#	f95  -o  slick90run $(FFLAGS) $(OBJslick) $(PNAG) $(RHLX)
#       gfortran also works here
#	gfortran  -o  slick90run $(FFLAGS) $(OBJslick) $(PNAG)
#       gfortran  -o  slick90run $(FFLAGS) $(OBJslick) $(PNAG) -framework IOKit -framework CoreFoundation -lc++
#       gfortran  -o  slick90run $(FFLAGS) $(OBJslick) $(PNAG) -framework IOKit -framework CoreFoundation
#	f95  -o  slick90run $(FFLAGS) $(OBJslick) $(PNAG)
	gfortran  -o  slicktrack $(FFLAGS) $(OBJslick) -llapack

clean: ;rm -f $(OBJslick); rm -f slicktrack
