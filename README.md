# tact

[![Build Status](https://travis-ci.org/jls713/tact.svg?branch=master)](https://travis-ci.org/jls713/tact)

Code for calculating actions, angles and frequencies in various ways

## Author

Jason Sanders -- jls at ast dot cam dot ac dot uk

Please cite the accompanying paper Sanders & Binney (2016) if you find the code useful.

## Requirements

1. [gsl](http://www.gnu.org/software/gsl/)

## Installation

* Make sure environment variable $(CXX) gives c++ compiler or specify compiler path in Makefile.inc (need a C++-11 compatible compiler, it compiles with clang 3.4 and g++ 4.8 or 4.9)
* Specify path to gsl in Makefile.inc
* Run make

This will install the basic package. To access all the features one should also install Torus and LAPACK:

* Some code uses [Torus](https://github.com/PaulMcMillan-Astro/Torus). To use this install Torus (use 'make CPP="$CXX -fPIC"' to ensure libraries are compiled with fPIC flag) , add path to Torus to Makefile.inc and run 'make TORUS=1'. Currently there are problems compiling with Torus using clang. This appears to be due to different compiler flags in Torus and tact. If users wish to use clang then they should make sure the compiler flags are the same for both.
* Some code uses [LAPACK](http://www.netlib.org/lapack/). To use this install LAPACK, add path to LAPACK to Makefile.inc and run 'make LAPACK=1'
* To do both run 'make TORUS=1 LAPACK=1'

One can also compile the code into a python module. This requires the boost library and the paths in Makefile.inc to python and boost to be correctly set. With these set either run 'make python' or use the setup.py in /aa like 'TORUS=1 LAPACK=1 python setup.py install'.

There is also test code that runs using googletest. However, there is a quick test to run detailed below.

## Methods

1. Analytic potentials (Isochrone and Harmonic oscillator)
2. General spherical potentials
3. [Cylindrical Adiabatic Approximation (CAA)](http://arxiv.org/abs/1109.4417), Schoenrich & Binney (2012)
4. Spheroidal Adiabatic Approximation (SAA) (unpublished, in my thesis)
5. [Stackel fitting](http://arxiv.org/abs/1208.2813), Sanders (2012)
6. [Axisymmetric Stackel fudge](http://arxiv.org/abs/1207.4910), Binney (2012)
7. [Interpolation using Axisymmetric Stackel fudge](http://arxiv.org/abs/1207.4910), Binney (2012)
8. [Triaxial Stackel fudge](http://arxiv.org/abs/1412.2093), Sanders & Binney (2015)
9. [Generating function from orbit (axisymmetric and triaxial, O2GF)](http://arxiv.org/abs/1401.3600), Sanders & Binney (2014)
10. [Average generating function from orbit (AvGF)](http://arxiv.org/abs/1401.2985), Bovy (2014), [Fox (2014)](http://arxiv.org/abs/1407.1688)
11. [Iterative Torus Machine (ItTC)](http://arxiv.org/abs/1412.2093), Sanders & Binney (2015)

## Docs

There is some documentation that can be produced by doxygen by running
```
make docs
```
The python library also contains doc strings.

## Test case

After successful compilation the command
```
cd aa; ./mains/test_actions.exe 8. 1. 0.2 40. 200. 50. acts.dat
```
should integrate an orbit with initial conditions X=(8. 1. 0.2) kpc and V = (40. 200. 50.)km/s in the potential

Phi(x)=Vc^2/2 log(R^2+(z/q)^2)

with Vc=220km/s and q=0.9 (or in the Piffl 2014 potential if Torus is installed) and compute the actions for each point using a variety of methods. The results are output in acts.dat with two columns per method (JR and Jz).

## Paper 

The accompanying paper is Sanders & Binney (2016). [action_comp_plots.py](aa/action_comp_plots.py) reproduces all but one of the plots in the paper. To produce the data for these plots the following scripts are available. 

1. Fig. 2 data is produced by the command 
```
./orbits.sh
```
2. Fig. 3, 4, 5, 6 data are produced by the command
```
./mains/./many_tori.exe many_tori_output.dat
```
3. Fig. 7 data are produced by
```
./orbits_converg.sh
```

