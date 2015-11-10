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
  # if "CCOMPILER" in l:
  #   os.environ["CXX"]=l[10:]
  #   os.environ["CC"]=os.environ["CXX"]
  if "GSLPATH" in l and GSLPath=="":
    GSLPath=l[9:].replace('\n', '')
  if "LAPACKPATH" in l:
    LapackPath=l[11:].replace('\n', '')
  if "TORUSPATH" in l:
    TorusPath=l[10:].replace('\n', '')
  if "EBFPATH" in l:
    EBFPath=l[8:].replace('\n', '')
  if "BOOSTINCPATH" in l:
    BoostPath=l[14:-8].replace('\n', '')

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

TorusFlag=0
if 'TORUS' in os.environ:
  TorusFlag=1
LapackFlag=0
if 'LAPACK' in os.environ:
  LapackFlag=1

setup(name="aa_py",
      ext_modules=[
          Extension("aa_py", ["src/aa_py.cpp"],
                    include_dirs=['inc',
                                  '../pot/inc',
                                  '../general/GSLInterface',
				                          '../general/',
                                  '../general/cuba/',
                                  '../general/coordtransforms/inc',
                                  '../general/jamestools/jamestools',
                                  '../general/jamestools/numrec',
                                  '../general/jamestools/octint',
				                          GSLPath+'/include',
				                          TorusPath+'/src',
                                  TorusPath+'/src/pot',
                                  TorusPath+'/src/utils',
                                  BoostPath+'include',
                                  np.get_include()],
                    library_dirs=['lib',
                                  TorusPath+'/obj/',
                                  TorusPath+'/WDlib/obj',
				                          GSLPath+'/lib',
                                  BoostPath+'lib',
                                  "../general/gnuplot/",
                                  "../general/cuba/",
                                  "../general/coordtransforms/",
                                  "../general/jamestools/jamestools",
                                  "../general/jamestools/numrec",
                                  "../general/jamestools/octint",
                                  "../pot/",
                                  EBFPath+"/lib"],
                    libraries=[
                        'aa','pot_js', 'plot','coords',
                        'jamestools','lapack','cuba',
                        'press_cp','octint','ebf_cpp',
                        'Torus', 'Other', 'Pot', 'WD',
                        'boost_python', 'gsl','gslcblas','m'],
                    extra_compile_args=['-fPIC', '-shared', '-std=c++0x',
                                        '-Wall', '-O3',
                                        '-ffast-math', '-fPIC', '-fopenmp',
                                        'TORUS='+str(TorusFlag),
                                        'LAPACK='+str(LapackFlag)])
      ])
