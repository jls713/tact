import numpy as np
import sys
sys.path.append('lib/')
import aa_py

X = np.array([8.,0.1,0.3,40.,220.,30.])
Pot = aa_py.GalPot('../Torus/pot/PJM11_best.Tpot')

## Fudge
AA = aa_py.Actions_AxisymmetricStackel_Fudge(Pot,-30.)
print AA.actions(X)

