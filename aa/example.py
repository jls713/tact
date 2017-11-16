### example.py
## Gives an example of how to use the python interface to the tact code
##
## Estimates actions using the variety of methods in Sanders & Binney (2016)
## in multi-component potential from Piffl et al. (2014)

import numpy as np
import sys
sys.path.append('lib/')
import aa_py

X = np.array([8.,0.1,0.3,40.,220.,30.])
Pot = aa_py.GalPot('../Torus/pot/Piffl14.Tpot')

print 'Computing actions for the phase-space point (units kpc, km/s): '
print X
print '============================================='

## Fudge
print ('Using Stackel Fudge:')
AA = aa_py.Actions_AxisymmetricStackel_Fudge(Pot,-20.)
print '\t Actions:',AA.actions(X)
print '\t Angles and Freqs:',AA.angles(X)

## Iterative Torus
print 'Using Iterative Torus:'
Tor = aa_py.IterativeTorusMachine(AA,Pot,1e-8,5,1e-5)
#print '\t Actions:',Tor.actions(X)[:3]

## Generating Function
print 'Using Generating Function -- a fuller interface is available through full_actions and full_angles -- the default parameters may need tweaking to avoid spurious results (e.g. negative vertical action for highly eccentric orbits'
AG = aa_py.Actions_Genfunc(Pot,"axisymmetric")
print '\t Actions:',AG.actions(X)
print '\t Angles and Freqs:',AG.angles(X)

## Average generating Function
print 'Averaging method'
AGav = aa_py.Actions_Genfunc_Average(Pot,"axisymmetric")
print '\t Actions:',AGav.actions(X)
print '\t Angles and Freqs:',AGav.angles(X)

## uvorb
print 'Using Stackel fudge with Delta from closed orbits'
UV = aa_py.Actions_AxisymmetricStackel_Fudge_DeltaGuess(Pot,1.,20.,10,10,"example.delta_uv")
print '\t Actions:',UV.actions(X)
print '\t Angles and Freqs:',UV.angles(X)

## Cylindrical Adiabatic
print 'Using cylindrical adiabatic approximation'
PAA = aa_py.Actions_CylindricalAdiabaticApproximation(Pot,"example.paa",True,False,1.,20.,10.,60)
print '\t Actions:',PAA.actions(X)
print '\t Angles and Freqs:',PAA.angles(X)

## Spheroidal Adiabatic
print 'Using spheroidal adiabatic approximation'
SAA = aa_py.Actions_SpheroidalAdiabaticApproximation(Pot,"example.saa",True,False,-40.,1.,20.,10.,60,10)
print '\t Actions:',SAA.actions(X)
print '\t Angles and Freqs:',SAA.angles(X)

## Stackel Fitting
print 'Using Stackel fitting'
SF = aa_py.Actions_StackelFit(Pot,1e-8)
print '\t Actions:',SF.actions(X)
print '\t Angles and Freqs:',SF.angles(X)
