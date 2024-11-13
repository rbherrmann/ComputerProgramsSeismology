#####
#	LINUX
#####
CCMP=gcc -O

##### DO NOT CHANGE BELOW HERE #####

MCHDEP= mchdep.o
MCHCMD= mnmarg.o mgtarg.o
RM = rm -f
DEST=..

.f.o:
	$(FCMP) -c $<
.c.o:	
	$(CCMP) -c $<

all: sacnoise clean

sacnoise: sacnoise.o sacsubc.o
	$(CCMP) sacnoise.o sacsubc.o -lm -o sacnoise

clean:
	${RM} *.o

