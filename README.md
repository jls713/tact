# tact

[![Build Status](https://travis-ci.org/jls713/tact.svg?branch=master)](https://travis-ci.org/jls713/tact)

Code for calculating actions, angles and frequencies in various ways

## Requirements

1. [gsl](http://www.gnu.org/software/gsl/)

## Installation

* Make sure environment variable $(CXX) gives c++ compiler or specify compiler path in Makefile.inc (need a C++-11 compatible compiler, it compiles with clang 3.4 and g++ 4.8 or 4.9)
* Specify path to gsl in Makefile.inc
* Run make

This will install the basic package. To access all the features one should also install Torus and LAPACK:

* Some code uses [Torus](https://github.com/PaulMcMillan-Astro/Torus). To use this install Torus (use 'make CFLAGS:=$(CFLAGS) -fPIC' to ensure libraries are compiled with fPIC flag) , add path to Torus to Makefile.inc and run 'make TORUS=1'
* Some code uses [LAPACK](http://www.netlib.org/lapack/). To use this install LAPACK, add path to LAPACK to Makefile.inc and run 'make LAPACK=1'
* To do both run 'make TORUS=1 LAPACK=1'

## Methods

1. Analytic potentials (Isochrone and Harmonic oscillator)
2. General spherical potentials
3. [Polar Adiabatic Approximation](http://arxiv.org/abs/1109.4417), Schoenrich & Binney (2012)
4. Spheroidal Adiabatic Approximation (unpublished, in my thesis)
5. [Stackel fitting](http://arxiv.org/abs/1208.2813), Sanders (2012)
6. [Axisymmetric Stackel fudge](http://arxiv.org/abs/1207.4910), Binney (2012)
7. [Interpolation using Axisymmetric Stackel fudge](http://arxiv.org/abs/1207.4910), Binney (2012)
8. [Triaxial Stackel fudge](http://arxiv.org/abs/1412.2093), Sanders & Binney (2014)
9. [Generating function from orbit (axisymmetric and triaxial)](http://arxiv.org/abs/1401.3600), Sanders & Binney (2014)
10. [Average generating function from orbit](http://arxiv.org/abs/1401.2985), Bovy (2014), [Fox (2014)](http://arxiv.org/abs/1407.1688)
11. [Iterative Torus Machine](http://arxiv.org/abs/1412.2093), Sanders & Binney (2014)


## Test case

After successful compilation the command
```
./mains/test_actions.exe 8. 1. 0.2 40. 200. 50. acts.dat
```
should integrate an orbit with initial conditions X=(8. 1. 0.2) kpc and V = (40. 200. 50.)km/s in the potential

Phi(x)=Vc^2/2 log(R^2+(z/q)^2)

with Vc=220km/s and q=0.9 and compute the actions for each point using a variety of methods. The results are output in tmp with two columns per method (JR and Jz).
