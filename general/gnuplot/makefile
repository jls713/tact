include ../../Makefile.inc

.cpp.o:
	$(CCOMPILER) -c $(OPTS) $<

default: clean all

all::	libplot.a

gnuplot_i.o:	gnuplot_i.cpp

libplot.a: gnuplot_i.o
	ar rcs libplot.a gnuplot_i.o

clean:
	rm -f libplot.a gnuplot_i.o
