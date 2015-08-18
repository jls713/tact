include ../../../Makefile.inc

default: clean all

all:  liboctint.a
forceall:	clean all

oct_int.o: oct_int.cpp
	$(CCOMPILER) $(OPTS) -c oct_int.cpp -o oct_int.o

liboctint.a: oct_int.o
	ar rcs liboctint.a oct_int.o

clean:
	rm -rf *~; rm -rf *.o; rm -rf *.so; rm -rf *.exe liboctint.a
