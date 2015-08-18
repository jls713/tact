#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
import numpy as np
import os

os.environ["CC"] = "/opt/ioa/software/gcc/4.7.2/bin/gcc"
os.environ["CXX"] = "/opt/ioa/software/gcc/4.7.2/bin/g++"

TorusPath='../Torus/'
GSLPath='/opt/ioa/software/gsl/1.16'
BoostPath='/opt/ioa/software/boost/1.55.0/'
setup(name="aa_py",
      ext_modules=[
          Extension("aa_py", ["src/aa_py.cpp"],
                    include_dirs=['inc','../pot/inc','../general/GSLInterface',
				  '../general/','../general/coordtransforms/inc',
				  GSLPath+'/include',
				  TorusPath + 'src',
                                  TorusPath + 'WDlib/inc',
                                  BoostPath+'include',np.get_include()],
                    library_dirs=['lib',TorusPath+'obj/',TorusPath+'WDlib/obj',
				  GSLPath+'/lib',BoostPath+'lib'],
                    libraries=[
                        'aa','Torus', 'Other', 'Pot', 'WD',
                        'boost_python', 'gsl','gslcblas','m'],
                    extra_compile_args=['-fPIC', '-shared', '-std=c++0x',
                                        '-Wall', '-O3',
                                        '-ffast-math', '-fPIC', '-fopenmp'])
      ])
