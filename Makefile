#
# Makefile for SLICKTRACK
#

FF=gfortran
FFLAGS= -O3
LFLAGS= -O3 -llapack
TARGET=slicktrack
OBJDIR=obj
SRCDIR=src

SOURCES=$(SRCDIR)/slick.f $(SRCDIR)/aver.f $(SRCDIR)/averv5sq.f $(SRCDIR)/avorb.f $(SRCDIR)/lbeambeam.f $(SRCDIR)/lbeambeam9.f $(SRCDIR)/cavity.f $(SRCDIR)/cfdiph.f $(SRCDIR)/cfdipv.f $(SRCDIR)/chrotn.f $(SRCDIR)/damint.f $(SRCDIR)/damper.f $(SRCDIR)/date.f $(SRCDIR)/dcgmpr.f $(SRCDIR)/dertop.f $(SRCDIR)/diaver.f $(SRCDIR)/diph.f $(SRCDIR)/dipv.f $(SRCDIR)/dx66.f $(SRCDIR)/dx88.f $(SRCDIR)/emitnc.f $(SRCDIR)/fixorb1.f $(SRCDIR)/fixorb2.f $(SRCDIR)/fixorb3.f $(SRCDIR)/fixorb4.f $(SRCDIR)/fixorb5.f $(SRCDIR)/hedge.f $(SRCDIR)/iparam.f $(SRCDIR)/jam333.f $(SRCDIR)/jam444.f $(SRCDIR)/jam666.f $(SRCDIR)/jam777.f $(SRCDIR)/jam888.f $(SRCDIR)/jam881.f $(SRCDIR)/jam999.f $(SRCDIR)/lattin.f $(SRCDIR)/linopt.f $(SRCDIR)/linopenopt.f $(SRCDIR)/mult.f $(SRCDIR)/mx66.f $(SRCDIR)/mx66damp.f $(SRCDIR)/mx88.f $(SRCDIR)/mx88damp.f $(SRCDIR)/mx99damp.f $(SRCDIR)/mxdamp.f $(SRCDIR)/nlbeambeam.f $(SRCDIR)/nlbeambeam9.f $(SRCDIR)/norm.f $(SRCDIR)/norm2.f $(SRCDIR)/orbit.f $(SRCDIR)/ortchk.f $(SRCDIR)/pertop.f $(SRCDIR)/psequil.f $(SRCDIR)/quartm.f $(SRCDIR)/renorm.f $(SRCDIR)/reson.f $(SRCDIR)/roquad.f $(SRCDIR)/rotatr.f $(SRCDIR)/rspin.f $(SRCDIR)/scrurita1.f $(SRCDIR)/scrurita2.f $(SRCDIR)/scrurita3.f $(SRCDIR)/scrurita4.f $(SRCDIR)/setbb.f $(SRCDIR)/sex88.f $(SRCDIR)/sexeq.f $(SRCDIR)/sig.f $(SRCDIR)/simq.f $(SRCDIR)/skewquad.f $(SRCDIR)/snake.f $(SRCDIR)/sol6.f $(SRCDIR)/sol66.f $(SRCDIR)/sol77.f $(SRCDIR)/sol8an.f $(SRCDIR)/sol9an.f $(SRCDIR)/solxyp.f $(SRCDIR)/spin.f $(SRCDIR)/ssect.f $(SRCDIR)/symp6.f $(SRCDIR)/symp8.f $(SRCDIR)/thiquad.f $(SRCDIR)/tspin.f $(SRCDIR)/unit.f $(SRCDIR)/unitm.f $(SRCDIR)/f03aaf.f $(SRCDIR)/f02ebf.f $(SRCDIR)/g05skf.f $(SRCDIR)/g05sqf.f $(SRCDIR)/g05kff.f $(SRCDIR)/g05rzf.f $(SRCDIR)/s15ddf.f $(SRCDIR)/rnglib.f $(SRCDIR)/ranlib.f $(SRCDIR)/string_trim.f90
OBJECTS=$(OBJDIR)/slick.o $(OBJDIR)/aver.o $(OBJDIR)/averv5sq.o $(OBJDIR)/avorb.o $(OBJDIR)/lbeambeam.o $(OBJDIR)/lbeambeam9.o $(OBJDIR)/cavity.o $(OBJDIR)/cfdiph.o $(OBJDIR)/cfdipv.o $(OBJDIR)/chrotn.o $(OBJDIR)/damint.o $(OBJDIR)/damper.o $(OBJDIR)/date.o $(OBJDIR)/dcgmpr.o $(OBJDIR)/dertop.o $(OBJDIR)/diaver.o $(OBJDIR)/diph.o $(OBJDIR)/dipv.o $(OBJDIR)/dx66.o $(OBJDIR)/dx88.o $(OBJDIR)/emitnc.o $(OBJDIR)/fixorb1.o $(OBJDIR)/fixorb2.o $(OBJDIR)/fixorb3.o $(OBJDIR)/fixorb4.o $(OBJDIR)/fixorb5.o $(OBJDIR)/hedge.o $(OBJDIR)/iparam.o $(OBJDIR)/jam333.o $(OBJDIR)/jam444.o $(OBJDIR)/jam666.o $(OBJDIR)/jam777.o $(OBJDIR)/jam888.o $(OBJDIR)/jam881.o $(OBJDIR)/jam999.o $(OBJDIR)/lattin.o $(OBJDIR)/linopt.o $(OBJDIR)/linopenopt.o $(OBJDIR)/mult.o $(OBJDIR)/mx66.o $(OBJDIR)/mx66damp.o $(OBJDIR)/mx88.o $(OBJDIR)/mx88damp.o $(OBJDIR)/mx99damp.o $(OBJDIR)/mxdamp.o $(OBJDIR)/nlbeambeam.o $(OBJDIR)/nlbeambeam9.o $(OBJDIR)/norm.o $(OBJDIR)/norm2.o $(OBJDIR)/orbit.o $(OBJDIR)/ortchk.o $(OBJDIR)/pertop.o $(OBJDIR)/psequil.o $(OBJDIR)/quartm.o $(OBJDIR)/renorm.o $(OBJDIR)/reson.o $(OBJDIR)/roquad.o $(OBJDIR)/rotatr.o $(OBJDIR)/rspin.o $(OBJDIR)/scrurita1.o $(OBJDIR)/scrurita2.o $(OBJDIR)/scrurita3.o $(OBJDIR)/scrurita4.o $(OBJDIR)/setbb.o $(OBJDIR)/sex88.o $(OBJDIR)/sexeq.o $(OBJDIR)/sig.o $(OBJDIR)/simq.o $(OBJDIR)/skewquad.o $(OBJDIR)/snake.o $(OBJDIR)/sol6.o $(OBJDIR)/sol66.o $(OBJDIR)/sol77.o $(OBJDIR)/sol8an.o $(OBJDIR)/sol9an.o $(OBJDIR)/solxyp.o $(OBJDIR)/spin.o $(OBJDIR)/ssect.o $(OBJDIR)/symp6.o $(OBJDIR)/symp8.o $(OBJDIR)/thiquad.o $(OBJDIR)/tspin.o $(OBJDIR)/unit.o $(OBJDIR)/unitm.o $(OBJDIR)/f03aaf.o $(OBJDIR)/f02ebf.o $(OBJDIR)/g05skf.o $(OBJDIR)/g05sqf.o $(OBJDIR)/g05kff.o $(OBJDIR)/g05rzf.o $(OBJDIR)/s15ddf.o $(OBJDIR)/rnglib.o $(OBJDIR)/ranlib.o $(OBJDIR)/string_trim.o

$(TARGET): $(OBJECTS)
	$(FF) -o $@ $(OBJECTS) $(LFLAGS)
	@echo "Linking complete!"

$(OBJDIR)/%.o : $(SRCDIR)/%.f
	$(FF) $(FFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FF) $(FFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONY: clean
clean: ;rm -f $(OBJECTS) $(OBJECTS90); rm -f $(TARGET)
