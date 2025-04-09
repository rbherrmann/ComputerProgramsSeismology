FCMP=gfortran

#CPS=/Users/rbh/PROGRAMS.310t/PROGRAMS.330
CALPLTF=$(CPS)/lib/libcalpltf.a

MCHDEP=mchdep.o
MCHCMD=mnmarg.o mgtarg.o


.f.o:
	$(FCMP) -c $<

mnmarg.o:
	cp $(CPS)/SUBS/mnmarg.f .
	$(FCMP) -c mnmarg.f

mgtarg.o:
	cp $(CPS)/SUBS/mgtarg.f .
	$(FCMP) -c mgtarg.f

mchdep.o:
	cp $(CPS)/SUBS/mchdep.f .
	$(FCMP) -c mchdep.f

tfmplot:	tfmplot.o mrs.o $(MCHDEP) $(MCHCMD)
	$(FCMP) tfmplot.o mrs.o $(MCHDEP) $(MCHCMD) $(CALPLTF) -o tfmplot

clean:
	rm -f *.o
