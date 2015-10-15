
# Makefile for aa

#------------------------------------------------------
include ../Makefile.inc
#------------------------------------------------------

SRC = src/
POT = pot/
OBJ = obj/
FALCONOBJ=$(OBJ)falcONobj/
#------------------------------------------------------
# Set paths to other code

GENERALPATH = ../general/
CTPATH = $(GENERALPATH)coordtransforms/
NEMOPATH = /data/jls/nemo/
FALCONPATH = $(NEMOPATH)usr/dehnen/falcON.P/
POTPATH = ../pot/

GNUPLOTLIBPATH = $(GENERALPATH)gnuplot
NUMRECLIBPATH = $(GENERALPATH)jamestools/numrec/
JAMESLIBPATH = $(GENERALPATH)jamestools/jamestools/
POTLIBPATH = ../pot/

CUBALIBPATH = -L$(GENERALPATH)cuba
CUBA = $(CUBALIBPATH) -lcuba

ifeq ($(TORUS), 1)
	TORUSINCPATH = -I$(TORUSPATH)src/ -I$(TORUSPATH)src/utils/ -I$(TORUSPATH)src/pot
	TORUSLIBPATH = -L$(TORUSPATH)obj/ -L$(TORUSPATH)WDlib/obj/ -L$(EBFPATH)lib/
	TORUSLIBS = -lTorus -lOther -lPot -lWD -lebf_cpp
else
	TORUSINCPATH =
	TORUSLIBPATH =
	TORUSLIBS =
endif

ifeq ($(LAPACK), 1)
	LAPACKLIBPATH=$(LAPACKPATH)
	LAPACKLIBS=-llapack -lblas
else
	LAPACKLIBPATH=.
	LAPACKLIBS=
endif

#------------------------------------------------------
# Set include and library paths

INC = -I../../ -I../ -Iinc/ $(TORUSINCPATH) -I$(PYTHONINCPATH) -I$(BOOSTINCPATH) -I$(GENERALPATH) -I$(GSLPATH)include/ -I$(POTPATH)inc/ -I$(CTPATH)inc/

LIBPATHS =
 # -Wl,-R$(GENERALPATH) -Wl,-R$(NUMRECLIBPATH) -Wl,-R$(JAMESLIBPATH)\
		   # -Wl,-R$(CUBALIBPATH) -Wl,-R$(GNUPLOTLIBPATH) -Wl,-R$(TORUSLIBPATH)\
		   # -Wl,-R$(GSLPATH)lib/ -Wl,-R$(POTLIBPATH) -Wl,-R$(CTPATH)

LIB = $(LIBPATHS) -L../aa/lib/ -Wl,-rpath,$(TORUSPATH)obj/ -Wl,-rpath,$(TORUSPATH)WDlib/obj/ -L$(GENERALPATH) -L$(NUMRECLIBPATH) -L$(JAMESLIBPATH) \
	$(CUBALIBPATH)\
	  -L$(GNUPLOTLIBPATH) $(TORUSLIBPATH) -L$(GSLPATH)lib/ -L$(POTLIBPATH)\
	  -L$(CTPATH) -laa -lpot_js -lcoords \
	  -lpress_cp $(TORUSLIBS) -lplot -lgsl -lgslcblas -lm $(EXTRACLIB) -L$(LAPACKLIBPATH) $(LAPACKLIBS) -ljamestools
BOOSTLIB = -lboost_python -lpython2.7
LIBPY = $(PYTHONLIBPATH) $(BOOSTLIBPATH) $(BOOSTLIB)

#------------------------------------------------------
# Additional stuff for falcON programs

FALCONINC = -I$(FALCONPATH)inc/ -I$(FALCONPATH)inc/public/ -I$(FALCONPATH)utils/inc/ -I$(NEMOPATH)inc/ -I$(NEMOPATH)inc/max/ -I$(NEMOPATH)lib/ -DfalcON_SINGLE -DfalcON_NEMO -DfalcON_PROPER

FALCONLIB = -L$(FALCONPATH)lib/ -lfalcON.P -L$(FALCONPATH)utils/lib/ -lWDutils -L$(NEMOPATH)lib/ -lnemo++ -lnemo -DfalcON_SINGLE -DfalcON_NEMO -DfalcON_PROPER

#------------------------------------------------------
# Set files

LIBRARY = libaa.so

OBJECTS = $(OBJ)debug.o $(OBJ)spherical_aa.o $(OBJ)stackel_aa.o  $(OBJ)get_closed_Rz.o $(OBJ)uv_orb.o $(OBJ)lmn_orb_ELLz.o $(OBJ)lmn_orb.o $(OBJ)it_torus.o $(OBJ)stackel_fit.o $(OBJ)analytic_aa.o $(OBJ)genfunc_aa.o $(OBJ)adiabatic_aa.o $(OBJ)tables_aa.o

HEADERS = $(wildcard inc/*.h)

GTESTINC=$(GTESTPATH)/include/
GTESTLIB=$(GTESTPATH)/lib/.libs/

#------------------------------------------------------

default: all

all: tests.exe

forceall: clean all

#------------------------------------------------------

$(OBJ)%.o : %.cpp ${HEADERS}
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) -I$(GTESTINC) $(TORUSFLAG) $(LAPACKFLAG)

#------------------------------------------------------

## TESTS

tests.exe: $(OBJ)test.o
	$(CCOMPILER) $(OPTS) -o tests.exe $(OBJ)test.o $(LIB) $(CUBA) -L$(GTESTLIB) -lgtest $(TORUSFLAG) $(LAPACKFLAG)

.PHONY: clean
clean:
	rm -rf obj/*.o; rm -rf *.so; rm -rf *.exe
