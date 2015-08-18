include ../../Makefile.inc

default: clean all

all:  libcoords.so
forceall:	clean all

#------------------------------------------------------
# General compiler & options
SRC = src/
POT = pot/
OBJ = obj/
OBJECTS = $(OBJ)coordtransforms.o $(OBJ)coordsys.o

#------------------------------------------------------

INC = -Iinc/ -I../ -I. -I$(GSLPATH)include/

#------------------------------------------------------

$(OBJ)%.o : $(SRC)%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC)

libcoords.so: $(OBJECTS)
	$(CCOMPILER) -shared $(OPTS) -o libcoords.so -Wl,-soname,libcoords.so $(OPTS) $(OBJECTS) $(INC)

clean:
	rm -rf obj/*.o; rm -rf *.so; rm -rf *.exe
