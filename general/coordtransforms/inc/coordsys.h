#ifndef COORDSYS_H
#define COORDSYS_H

#include "utils.h"

class UVOblateSpheroidCoordSys{
	// (u,v) coords
	private:
		double Delta, Delta2;
	public:
		UVOblateSpheroidCoordSys(double Delta): Delta(Delta),Delta2(Delta*Delta){}
		inline double delta(void){ return Delta;}
		VecDoub xv2uv(const VecDoub& x);
		VecDoub uv2Rz(const VecDoub& x);
};

class OblateSpheroidCoordSys{
	// (lam,nu) coords
	private:
		double Alpha;
		const double Gamma;
	public:
		OblateSpheroidCoordSys(double alpha): Alpha(alpha), Gamma(-1.){}
		inline double alpha(void){ return Alpha;}
		inline double gamma(void){ return Gamma;}
		inline void newalpha(double a){ Alpha = a;}
		VecDoub x2tau(const VecDoub& x);
		VecDoub tau2x(const VecDoub& x);
		VecDoub tau2polar(const VecDoub& x);
		VecDoub xv2tau(const VecDoub& x);
		VecDoub derivs(const VecDoub& x);
		VecDoub tau2p(const VecDoub& tau);
};


class ConfocalEllipsoidalCoordSys{
	// (lam,mu,nu) coords
	private:
		double Alpha, Beta, Gamma;
	public:
		ConfocalEllipsoidalCoordSys(double a, double b): Alpha(a), Beta(b), Gamma(-1.){}
		inline double alpha(void){ return Alpha;}
		inline double beta(void){ return Beta;}
		inline double gamma(void){ return Gamma;}
		inline void newalpha(double a){ Alpha = a;}
		inline void newbeta(double b){ Beta = b;}
		VecDoub x2tau(const VecDoub& x);
		VecDoub tau2x(const VecDoub& tau);
		VecDoub xv2tau(const VecDoub& x);
		VecDoub tau2p(const VecDoub& tau);
		VecDoub derivs(const VecDoub& x);
};

#endif
