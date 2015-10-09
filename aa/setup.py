#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
import numpy as np
import os

GSLPath=""
LapackPath=""
TorusPath=""
EBFPath=""
BoostPath=""

s = open("../Makefile.inc","r")
for l in s:
  if "CCOMPILER" in l:
    os.environ["CXX"]=l[10:]
    os.environ["CC"]=l[10:-3]+"gcc"
  if "GSLPATH" in l:
    GSLPath=l[8:]
  if "LAPACKPATH" in l:
    LapackPath=l[11:]
  if "TORUSPATH" in l:
    TorusPath=l[10:]
  if "EBFPATH" in l:
    EBFPath=l[8:]
  if "BOOSTINCPATH" in l:
    BoostPath=l[13:-7]

if GSLPath=="":
    print "GSLPath not found"
if TorusPath=="":
    print "TorusPath not found"
if LapackPath=="":
    print "LapackPath not found"
if EBFPath=="":
    print "EBFPath not found"
if BoostPath=="":
    print "BoostPath not found"

setup(name="aa_py",
      ext_modules=[
          Extension("aa_py", ["src/aa_py.cpp"],
                    include_dirs=['inc',
                                  '../pot/inc',
                                  '../general/GSLInterface',
				                          '../general/',
                                  '../general/coordtransforms/inc',
                                  '../general/jamestools/jamestools',
                                  '../general/jamestools/numrec',
                                  '../general/jamestools/octint',
				                          GSLPath+'/include',
				                          TorusPath + 'src',
                                  TorusPath + 'src/pot',
                                  TorusPath + 'src/utils',
                                  BoostPath+'include',
                                  np.get_include()],
                    library_dirs=['lib',
                                  TorusPath+'obj/',
                                  TorusPath+'WDlib/obj',
				                          GSLPath+'/lib',
                                  BoostPath+'lib',
                                  "../general/gnuplot/",
                                  "../general/coordtransforms/",
                                  "../general/jamestools/jamestools",
                                  "../general/jamestools/numrec",
                                  "../general/jamestools/octint",
                                  "../pot/",
                                  EBFPath+"lib"],
                    libraries=[
                        'aa','pot_js', 'plot','coords',
                        'jamestools','lapack','cuba',
                        'press_cp','octint','ebf_cpp',
                        'Torus', 'Other', 'Pot', 'WD',
                        'boost_python', 'gsl','gslcblas','m'],
                    extra_compile_args=['-fPIC', '-shared', '-std=c++0x',
                                        '-Wall', '-O3',
                                        '-ffast-math', '-fPIC', '-fopenmp'])
      ])
