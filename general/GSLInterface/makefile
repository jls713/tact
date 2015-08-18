default: clean all

all:  test
forceall:	clean all

HOST=$(shell hostname)
ifeq ($(HOST), calx019.ast.cam.ac.uk)
	CPP = /opt/ioa/software/gcc/4.7.2/bin/g++
	CFLAGS = -std=c++0x -Wall -O3 -ffast-math -fno-strict-aliasing
	GSLPATH = /opt/ioa/software/gsl/1.16/
else
	CPP =  clang++
	CFLAGS = -std=c++11 -stdlib=libc++ -Wall -O3 -ffast-math -fno-strict-aliasing -Wno-deprecated-register
	EXTRACLIB =  -lc++
endif

test.o: test.cpp
	$(CPP) $(CFLAGS) -c test.cpp -I$(GSLPATH)include/ -o test.o

test.exe: test.o
	$(CPP) $(CFLAGS) -o test.exe test.o -L$(GSLPATH)lib/ -lgsl -lgslcblas -lm $(EXTRACLIB)

clean:
	rm -rf *.o; rm -rf *.so; rm -rf *.exe libGSLInterface.a
