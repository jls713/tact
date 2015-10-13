
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

TORUSINCPATH = -I$(TORUSPATH)src/ -I$(TORUSPATH)src/utils/ -I$(TORUSPATH)src/pot
TORUSLIBPATH = -L$(TORUSPATH)obj/ -L$(TORUSPATH)WDlib/obj/

CUBALIBPATH = -L$(GENERALPATH)cuba
CUBA = $(CUBALIBPATH) -lcuba

#------------------------------------------------------
# Set include and library paths

INC = -I../../ -I../ -Iinc/ $(TORUSINCPATH) -I$(PYTHONINCPATH) -I$(BOOSTINCPATH) -I$(GENERALPATH) -I$(GSLPATH)include/ -I$(POTPATH)inc/ -I$(CTPATH)inc/

LIBPATHS = -Wl,-R$(GENERALPATH) -Wl,-R$(NUMRECLIBPATH) -Wl,-R$(JAMESLIBPATH)\
		   -Wl,-R$(CUBALIBPATH) -Wl,-R$(GNUPLOTLIBPATH) -Wl,-R$(TORUSLIBPATH)\
		   -Wl,-R$(GSLPATH)lib/ -Wl,-R$(POTLIBPATH) -Wl,-R$(CTPATH)

LIB = $(LIBPATHS) -L$(GENERALPATH) -L$(NUMRECLIBPATH) -L$(JAMESLIBPATH) \
	$(CUBALIBPATH)\
	  -L$(GNUPLOTLIBPATH) $(TORUSLIBPATH) -L$(GSLPATH)lib/ -L$(POTLIBPATH)\
	  -L$(CTPATH) -lpot_js -lcoords \
	  -lpress_cp -lTorus -lOther -lPot -lWD -L$(EBFPATH) -lebf_cpp -lplot -lgsl -lgslcblas -lm $(EXTRACLIB) -L$(LAPACKPATH) -llapack -lblas -ljamestools
BOOSTLIB = -lboost_python -lpython2.7
LIBPY = $(PYTHONLIBPATH) $(BOOSTLIBPATH) $(BOOSTLIB)

#------------------------------------------------------
# Additional stuff for falcON programs

FALCONINC = -I$(FALCONPATH)inc/ -I$(FALCONPATH)inc/public/ -I$(FALCONPATH)utils/inc/ -I$(NEMOPATH)inc/ -I$(NEMOPATH)inc/max/ -I$(NEMOPATH)lib/ -DfalcON_SINGLE -DfalcON_NEMO -DfalcON_PROPER

FALCONLIB = -L$(FALCONPATH)lib/ -lfalcON.P -L$(FALCONPATH)utils/lib/ -lWDutils -L$(NEMOPATH)lib/ -lnemo++ -lnemo -DfalcON_SINGLE -DfalcON_NEMO -DfalcON_PROPER

#------------------------------------------------------
# Set files

LIBRARY = libaa.so

OBJECTS = $(OBJ)debug.o $(OBJ)spherical_aa.o $(OBJ)stackel_aa.o  $(OBJ)get_closed_Rz.o $(OBJ)uv_orb.o $(OBJ)lmn_orb.o $(OBJ)stackel_fit.o $(OBJ)analytic_aa.o $(OBJ)genfunc_aa.o $(OBJ)adiabatic_aa.o $(OBJ)tables_aa.o
ifeq ($(TORUS),1)
	$(OBJECTS) += $(OBJ)it_torus.o
endif
HEADERS = $(wildcard inc/*.h)

#------------------------------------------------------

default: all

all: lib/libaa.so mains/test_actions.exe mains/test_angles.exe mains/many_tori.exe

forceall: clean all

#------------------------------------------------------

$(OBJ)%.o : $(SRC)%.cpp ${HEADERS}
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC)

$(OBJ)%.o : mains/%.cpp ${HEADERS}
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC)

$(FALCONOBJ)%.o : $(SRC)%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) $(FALCONINC)

#------------------------------------------------------

## TESTS
analytic_test.exe: $(OBJECTS)
	$(CCOMPILER) $(OPTS) -o analytic_test.exe $(OBJECTS) $(LIB)

genfunc_test.exe: $(OBJECTS)
	$(CCOMPILER) $(OPTS) -o genfunc_test.exe $(OBJECTS) $(LIB) $(CUBA)

test.exe: $(OBJECTS) $(OBJ)test.o
	$(CCOMPILER) $(OPTS) -o test.exe $(OBJECTS) $(OBJ)test.o $(LIB) $(CUBA)

mains/test_actions.exe: $(OBJECTS) $(OBJ)test_actions.o
	$(CCOMPILER) $(OPTS) -o mains/test_actions.exe $(OBJECTS) $(OBJ)test_actions.o $(LIB) $(CUBA) -static-libstdc++

mains/test_angles.exe: $(OBJECTS) $(OBJ)test_angles.o
	$(CCOMPILER) $(OPTS) -o mains/test_angles.exe $(OBJECTS) $(OBJ)test_angles.o $(LIB) $(CUBA) -static-libstdc++

mains/many_tori.exe: $(OBJECTS) $(OBJ)many_tori.o
	$(CCOMPILER) $(OPTS) -o mains/many_tori.exe $(OBJECTS) $(OBJ)many_tori.o $(LIB) $(CUBA)

spherical_test.exe: $(OBJECTS)
	$(CCOMPILER) $(OPTS) -o spherical_test.exe $(OBJECTS) $(LIB)

stackel_fit_test.exe: $(OBJECTS)
	$(CCOMPILER) $(OPTS) -o stackel_fit_test.exe $(OBJECTS) $(LIB) $(CUBA)

lmnnew_test.exe: $(OBJECTS)
	$(CCOMPILER) $(OPTS) -o lmnnew_test.exe $(OBJECTS) $(LIB) $(CUBA)
## LIBRARY

lib/libaa.so: $(OBJECTS)
	$(CCOMPILER) -fPIC -shared $(CFL) $(OBJECTS) $(INC) -o lib/libaa.so

triax_py.so: $(OBJECTS) $(OBJ)triax_py.o
	$(CCOMPILER) -fPIC -shared $(CFL) $(OBJECTS) $(OBJ)triax_py.o $(INC) $(LIB)  -o triax_py.so

lib/aa_py.so: $(OBJECTS) $(OBJ)aa_py.o
	$(CCOMPILER) -fPIC -shared $(OPTS) $(OBJECTS) $(OBJ)aa_py.o $(INC) $(LIB) $(LIBPY) $(CUBA) -o lib/aa_py.so

## OTHER

falcON_aa.exe: $(OBJECTS) $(FALCONOBJ)falcON_aa.o
	$(CCOMPILER) $(OPTS) -o falcON_aa.exe $(OBJECTS) $(FALCONOBJ)falcON_aa.o $(FALCONLIB) $(LIB) $(LIBPY) $(CUBA)

.PHONY: clean
clean:
	rm -rf obj/*.o; rm -rf *.so; rm -rf *.exe
