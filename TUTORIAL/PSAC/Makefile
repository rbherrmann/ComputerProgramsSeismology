FCMP=gfortran
CALPLTF=${CPS}/lib/libcalpltf.a

.f.o:
	$(FCMP) -c $<

grphsubf.f:
	cp ${CPS}/SUBS/grphsubf.f .
grphsubf.o: grphsubf.f
sacsubf.f:
	cp ${CPS}/SUBS/sacsubf.f .
sacsubf.o: sacsubf.f
lgstr.f:
	cp ${CPS}/SUBS/lgstr.f .
lgstr.o: lgstr.f
mnmarg.f:
	cp ${CPS}/SUBS/mnmarg.f .
mnmarg.o: mnmarg.f
mgtarg.f:
	cp ${CPS}/SUBS/mgtarg.f .
mgtarg.o: mgtarg.f

psac: psac.o lgstr.o sacsubf.o grphsubf.o mnmarg.o mgtarg.o
	$(FCMP) -o psac  psac.o lgstr.o sacsubf.o grphsubf.o mnmarg.o mgtarg.o  ${CALPLTF}

clean:
	rm -f *.o
