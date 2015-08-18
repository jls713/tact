import numpy as np
import sys
sys.path.append('lib/')
import aa_py

X = np.array([8.,0.1,0.3,40.,220.,30.])
Pot = aa_py.GalPot('../Torus/pot/PJM11_best.Tpot')

## Fudge
AA = aa_py.Actions_AxisymmetricStackel_Fudge(Pot,-30.)
print AA.actions(X)

## Iterative Torus
Tor = aa_py.IterativeTorusMachine(AA,Pot,1e-8,5,1e-5)
print Tor.actions(X)[:3]

## Generating Function
AG = aa_py.Actions_Genfunc(Pot,"axisymmetric")
print AG.actions(X)

## Average generating Function
AGav = aa_py.Actions_Genfunc_Average(Pot,"axisymmetric")
print AGav.actions(X)

## uvorb
UV = aa_py.Actions_AxisymmetricStackel_Fudge_DeltaGuess(Pot,1.,20.,10,10,"example.delta_uv")
print UV.actions(X)

## Polar Adiabatic
PAA = aa_py.Actions_PolarAdiabaticApproximation(Pot,"example.paa",True,True)
print PAA.actions(X)

## Spheroidal Adiabatic
SAA = aa_py.Actions_SpheroidalAdiabaticApproximation(Pot,"example.saa",True,True,-30.)
print SAA.actions(X)

## Spheroidal Adiabatic
SF = aa_py.Actions_StackelFit(Pot)
print SF.actions(X)
