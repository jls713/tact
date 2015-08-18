#ifndef ADIABATIC_H
#define ADIABATIC_H

#include "potential.h"
#include "aa.h"

// ==========================================================================================
// Ellipsoidal Adiabatic Approximation
// ==========================================================================================

class Actions_EllipsoidalAdiabaticApproximation : public Action_Finder{
	private:
        bool alpha_guess;
	public:
        Potential_JS *Pot;
		OblateSpheroidCoordSys *CS;
		Actions_EllipsoidalAdiabaticApproximation(Potential_JS *pot): Pot(pot){
			CS = new OblateSpheroidCoordSys(-2.);
            alpha_guess = 1;
			// alpha set for each PSP
		}
        void set_alpha(double a){CS->newalpha(a);}
		VecDoub actions(VecDoub x);
};

struct action_struct_aa{
    Actions_EllipsoidalAdiabaticApproximation *AA;
    double ENugl, lamgl, Deltagl, taubargl, Lz;
    action_struct_aa(Actions_EllipsoidalAdiabaticApproximation *AA,double ENugl, double lamgl, double Lz, double Deltagl, double taubargl)
        :AA(AA),ENugl(ENugl),lamgl(lamgl),Lz(Lz)
        ,Deltagl(Deltagl),taubargl(taubargl){}
};

#endif
