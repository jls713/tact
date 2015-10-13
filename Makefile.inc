####################################################################
## Makefile
## 1. Need to install Torus (https://github.com/PaulMcMillan-Astro/Torus)
## 2. Need to install LAPACK (http://www.netlib.org/lapack/)
## 2. Need to install gsl (http://www.gnu.org/software/gsl/)
####################################################################

## You will need to change the compiler and the path to gsl

CCOMPILER=$(CXX)
CCOMPILER=g++-4.9
GSLPATH = /opt/ioa/software/gsl/1.16/
LAPACKPATH = general/lapack/
TORUSPATH = /home/jls/work/code/Torus/
EBFPATH = /home/jls/work/code/Torus/libebf_c_cpp-0.0.3/

####################################################################
DEBUG = -ggdb3 #-ffast-math -Ofast
OPTS = -std=c++0x -Wall -O2 -fopenmp -fPIC

## These paths are only important if you want to make python modules

BOOSTINCPATH =/opt/ioa/software/boost/1.55.0/include
BOOSTLIBPATH =-L/opt/ioa/software/boost/1.55.0/lib
PYTHONINCPATH=/opt/ioa/software/python/2.7.8/include/python2.7/
PYTHONLIBPATH=-L/opt/ioa/software/python/2.7.8/lib/

####################################################################

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LIBNAME =soname
    LIBFLAG=-shared
endif
ifeq ($(UNAME_S),Darwin)
    LIBNAME=install_name
    LIBFLAG=-dynamic_lib -undefined suppress -flat_namespace
endif
