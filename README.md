# tact

[![Build Status](https://travis-ci.org/jls713/tact.svg?branch=master)](https://travis-ci.org/jls713/tact)

Code for calculating actions, angles and frequencies in various ways

## Requirements

1. [Torus](https://github.com/PaulMcMillan-Astro/Torus)
2. [LAPACK](http://www.netlib.org/lapack/)
3. [gsl](http://www.gnu.org/software/gsl/)

## Installation

* Specify compiler path in Makefile.inc (need a C++-11 compatible compiler, currently compiles with g++ 4.9)
* Specify path to Torus, LAPACK and gsl in Makefile.inc
* Run make

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
