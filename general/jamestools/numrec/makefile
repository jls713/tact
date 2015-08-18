include ../../../Makefile.inc

default: clean all

all:  libpress_cp.a
forceall:	clean all

pressstuff.o: pressstuff.cpp
	$(CCOMPILER) -c $(OPTS) pressstuff.cpp -o pressstuff.o

libpress_cp.a: pressstuff.o
	ar rcs libpress_cp.a pressstuff.o

clean:
	rm -rf *~; rm -rf *.o; rm -rf *.so; rm -rf *.exe libpress_cp.a
