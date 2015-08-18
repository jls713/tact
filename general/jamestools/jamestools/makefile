include ../../../Makefile.inc

default: clean all

all:  libjamestools.a
forceall:	clean all

jamestools.o: jamestools.cpp
	$(CCOMPILER) $(OPTS) -c jamestools.cpp -o jamestools.o

libjamestools.a: jamestools.o
	ar rcs libjamestools.a jamestools.o

clean:
	rm -rf *~; rm -rf *.o; rm -rf *.so; rm -rf *.exe libjamestools.a
