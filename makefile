include Makefile.inc

COORD_DIR = general/coordtransforms/
CUBA_DIR = general/cuba/
POT_DIR = pot
AA_DIR = aa
NR_DIR = general/jamestools/numrec
JT_DIR = general/jamestools/jamestools
OCTINT_DIR = general/jamestools/octint
GNUPLOT_DIR = general/gnuplot

ifeq ($(TORUS),1)
	TORUSFLAG=-DTORUS
else
	TORUSFLAG=
endif

export TORUSFLAG
export TORUS

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LIBNAME =soname
    LIBFLAG=-shared
endif
ifeq ($(UNAME_S),Darwin)
    LIBNAME=install_name
    LIBFLAG=-dynamic_lib -undefined suppress -flat_namespace
endif


.PHONY: aa_code ct_code pot_code cuba_code octint_code jt_code nr_code gnuplot_code

all: cuba_code ct_code jt_code nr_code octint_code gnuplot_code pot_code aa_code

gnuplot_code:
	$(MAKE) -C $(GNUPLOT_DIR)

jt_code:
	$(MAKE) -C $(JT_DIR)

nr_code:
	$(MAKE) -C $(NR_DIR)

octint_code:
	$(MAKE) -C $(OCTINT_DIR)

cuba_code:
	$(MAKE) -C $(CUBA_DIR)

ct_code:
	$(MAKE) -C $(COORD_DIR)

pot_code:
	$(MAKE) -C $(POT_DIR)

aa_code:
	$(MAKE) -C $(AA_DIR)

docs:
	doxygen doc/doxygen.config

