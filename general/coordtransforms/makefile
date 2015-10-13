include ../../Makefile.inc

default: clean all

all:  libcoords.so
forceall:	clean all

#------------------------------------------------------
# General compiler & options
SRC = src/
POT = pot/
OBJ = obj/
OBJECTS = $(OBJ)coordtransforms.o $(OBJ)coordsys.o $(OBJ)utils.o

#------------------------------------------------------

INC = -Iinc/ -I../ -I. -I$(GSLPATH)include/

#------------------------------------------------------

$(OBJ)utils.o : ../utils.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC)

$(OBJ)%.o : $(SRC)%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC)

libcoords.so: $(OBJECTS)
	$(CCOMPILER) $(LIBFLAG) $(OPTS) -o libcoords.so -Wl,-$(LIBNAME),$(CURDIR)/libcoords.so $(OBJECTS) $(INC) -lgsl -lgslcblas

clean:
	rm -rf obj/*.o; rm -rf *.so; rm -rf *.exe
