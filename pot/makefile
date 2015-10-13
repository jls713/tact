include ../Makefile.inc

default: clean all

all:  libpot_js.so
forceall:	clean all

#------------------------------------------------------
SRC = src/
POT = pot/
OBJ = obj/

OBJECTS = $(OBJ)potential.o $(OBJ)Multipole.o $(OBJ)orbit.o

HEADERS = $(wildcard inc/*.h)

#------------------------------------------------------

GENERALPATH = ../general/
CTPATH = ../general/coordtransforms/
GNUPLOTLIBPATH = $(GENERALPATH)gnuplot

ifeq ($(TORUS), 1)
	TORUSINCPATH = -I$(TORUSPATH)src/ -I$(TORUSPATH)src/utils/ -I$(TORUSPATH)src/pot
	TORUSLIBPATH = -L$(TORUSPATH)obj/ -L$(TORUSPATH)WDlib/obj/
else
	TORUSINCPATH =
	TORUSLIBPATH =
endif

COORDINCPATH = $(CTPATH)inc/
COORDLIBPATH = $(CTPATH)

INC = -I../../ -Iinc/ $(TORUSINCPATH) -I$(GENERALPATH) -I$(GSLPATH)include/ -I$(COORDINCPATH)

LIB = -L$(GENERALPATH) -L$(GNUPLOTLIBPATH) $(TORUSLIBPATH) -L$(GSLPATH)lib/ \
	  -L$(COORDLIBPATH) -lcoords -lPot -lWD -lplot -lgsl -lgslcblas \
	  -lm $(EXTRACLIB)

BOOSTLIB = -lboost_python -lpython2.7
LIBPY = #$(PYTHONLIBPATH) $(BOOSTLIBPATH) $(BOOSTLIB)
#------------------------------------------------------

$(OBJ)%.o : $(SRC)%.cpp ${HEADERS}
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC)

multipole_test.exe: $(OBJECTS)
	$(CCOMPILER) $(OPTS) -o multipole_test.exe $(OBJECTS) $(LIB) $(LIBPY)

libpot_js.so: $(OBJECTS)
	$(CCOMPILER) -shared -Wl,-R$(CTPATH) $(OPTS) -o libpot_js.so -Wl,-soname,libpot_js.so $(OPTS) $(OBJECTS) $(INC)

clean:
	rm -rf obj/*.o; rm -rf *.so; rm -rf *.exe
