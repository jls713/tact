#ifndef IT_TORUS_H
#define IT_TORUS_H

#include "potential.h"
#include "stackel_aa.h"
#include "Torus.h"

class IterativeTorusMachine : public Action_Finder{
	private:
		Actions_AxisymmetricStackel_Fudge *AAA;
		GalPot *Phi;
		double eta;
		int MaxIterations;
		double dJ;
		VecDoub RefineGuess(VecDoub PSP, VecDoub aa, Potential *Phi, double *min, Angles &theta_new, double dJ);
		VecDoub NewActionPoint(VecDoub startpt, VecDoub oldpt, VecDoub newpt);
		VecDoub find_aa(VecDoub x);
		public:
			IterativeTorusMachine(Actions_AxisymmetricStackel_Fudge *AAA,GalPot *Phi, double eta = 1e-8, int MaxIterations=5, double dJ = 1e-5):AAA(AAA),Phi(Phi),eta(eta),MaxIterations(MaxIterations),dJ(dJ){};
			inline void set_eta(double s){eta = s;}
			inline void set_maxit(double s){MaxIterations = s;}
			inline void set_dJ(double s){dJ = s;}
			inline void set_options(double a, double b, double c)
				{eta=a;MaxIterations=b;dJ=c;}
			VecDoub actions(const VecDoub& x, void*params=nullptr);
};

VecDoub findXV(Angles theta, Torus *T, double sign);

struct Min_ItTorusMac_struct{
	Torus *T;VecDoub XV; double freqWeight;double sign;
	Min_ItTorusMac_struct(Torus *T, VecDoub XV,double freqWeight, double sign):T(T), XV(XV),freqWeight(freqWeight),sign(sign){}
};

#endif
