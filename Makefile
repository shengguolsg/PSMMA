include SLmake.inc.intel
#include SLmake.inc.th2

F90 = $(F77) 

PMLIBS = $(TOOLSLIB) $(BLACSLIB) $(BLASLIB)

OBJS= Cauchylowrank.o Auxil_2p5D.o pscauchy_compute.o psmma_cauchy.o 

TESTDRV2 = psmma_testcauchy.o

all: test2

test2: $(TESTDRV2) $(OBJS) 
	$(F77LOADER) $(F77LOADFLAGS) -o $@ $(TESTDRV2) $(OBJS) $(PMLIBS)

%.o:%.f90
	$(F90) -c $(F77FLAGS) $<

%.o:%.f
	$(F90) -c $(F77FLAGS) $<

clean :
	rm -f *.o *.mod
	rm test2
