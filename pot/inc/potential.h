// ============================================================================
/// \file inc/potential.h
// ============================================================================
/// \author Jason Sanders
/// \date 2014-2015
/// Institute of Astronomy, University of Cambridge (and University of Oxford)
// ============================================================================

// ============================================================================
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// ============================================================================
/*! \brief Potential classes
 *
 *  We give a list of the available potentials. Most potentials are
 *  self-documented, in that they contain a string that outputs their info.

 *  1. Potential_JS -- Base class
 *  2. SphericalPotential -- Base spherical class
 *  3. IsochronePotential -- Isochrone spherical potential
 *  4. StackelOblate_PerfectEllipsoid -- axisymmetric oblate Stackel potential perfect ellipsoid
 *  5. StackelTriaxial -- triaxial perfect ellipsoid
 *  6. Logarithmic -- triaxial logarithmic potential
 *  7. PowerLaw -- triaxial power-law potential
 *  8. Isochrone -- triaxial isochrone potential
 *  9. HarmonicOscillator -- triaxial harmonic oscillator
 *  10. Dehnen -- spherical Dehnen potential
 *  11. MiyamotoNagai_JS -- axisymmetric Miyamoto-Nagai potential
 *  11. NFW -- triaxial NFW potential
 *  12. Hernquist -- triaxial Hernquist potential
 *  13. Bulge -- triaxial Jaffe bulge potential
 *  14. GalPot -- wrapper for Torus GalPot potential
 *  15. MultiComponentPotential -- sum of potentials
 *  16. MultiComponentSphericalPotential -- sum of spherical potentials
 *  17. NFWSpherical -- spherical NFW potential
 *  18. HernquistSpherical -- spherical Hernquist potential
 *  19. PowerLawSpherical-- spherical power-law potential
 *  20. BowdenNFW -- flattened NFW from Bowden et al.(2014)
 *  21. WrapperTorusPotential -- wraps potentials from Torus code
 */
//============================================================================

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <cmath>
#include "coordsys.h"
#include "coordtransforms.h"
#include "GSLInterface/GSLInterface.h"
#ifdef TORUS
#include "falPot.h"
#include "Torus.h"
#endif

//============================================================================
/// General base class for Potential_JS
//============================================================================
class Potential_JS{
	public:
		inline virtual std::string name(void) const{
			return "You haven't written a description of this potential"; }
		inline virtual std::string params(void) const{ return "No param list"; }

		inline virtual double Phi(const VecDoub& x){
			/* potential at Cartesian x 				 */
			/* -- to be overridden by derived classes -- */
			return 0.;
		}
		inline virtual VecDoub Forces(const VecDoub& x){
			/* forces at Cartesian x 					 */
			/* -- to be overridden by derived classes -- */
			return VecDoub(0.);
		}
		inline virtual double density(const VecDoub& x){
			/* density at Cartesian x 		   			 */
			/* -- to be overridden by derived classes -- */
			return 0.;
		}
		inline virtual double Vc(double R){
			/* circular velocity  at polar R 			*/
			/* -- to be overridden by derived classes -- */
			return 0.;
		}
		inline double Mass(double r){
			/* returns mass in spherical r*/
			return r*r*(-Forces({r,0.,0.})[0])/conv::G;
		}
		inline double H(const VecDoub& xv){
			/* returns Hamiltonian at Cartesian (x,v) */
			assert(xv.size()>=6);
			VecDoub x = {xv[0],xv[1],xv[2]};
			return 0.5*(xv[3]*xv[3]+xv[4]*xv[4]+xv[5]*xv[5])+Phi(x);
		}
		inline double L(const VecDoub& xv){
			/* returns ang mom at Cartesian (x,v) */
			assert(xv.size()>=6);
			double Lx = xv[1]*xv[5]-xv[2]*xv[4];
			double Ly = xv[2]*xv[3]-xv[0]*xv[5];
			double Lz = xv[0]*xv[4]-xv[1]*xv[3];
			return sqrt(Lx*Lx+Ly*Ly+Lz*Lz);
		}
		inline VecDoub Lvec(const VecDoub& xv){
			/* returns ang mom at Cartesian (x,v) */
			assert(xv.size()>=6);
			return {xv[1]*xv[5]-xv[2]*xv[4],
					xv[2]*xv[3]-xv[0]*xv[5],
					xv[0]*xv[4]-xv[1]*xv[3]};
		}
		inline double Lz(const VecDoub& xv){
			assert(xv.size()>=6);
			/* returns z-component of angular momentum at Cartesian (x,v) */
			return xv[0]*xv[4]-xv[1]*xv[3];
		}
		inline virtual VecDoub freqs(double R){
			/* epicyclic frequencies at cylindrical R 	 */
			/* -- to be overridden by derived classes -- */
			std::cerr<<"Epicyclic frequency routine not written\n";
			return {0.,0.,0.};
		}
		VecDoub dPhidRdz(const VecDoub& Rz);
			// returns second derivative of potential wrt R and z
			// assumes axisymmetry -- takes VecDoub Rz = {R,z}
			// returns (dP/dRdR, dPdRdz, dPdzdz)
		double DeltaGuess(const VecDoub& x);

        const double Phi_max(double a = 1e5){
        	return Phi({a,a,a});
        }
        inline double E_circ(double R){
            return -Forces({R,0.,0.})[0]*.5*R+Phi({R,0.,0.});
        }
        inline double L_circ(double R){
			return sqrt(R*-Forces({R,0.,0.})[0])*R;
		}
        double R_E(double E,double r = -1.);
        double L_E(double E);
        double R_E(const VecDoub &x);
        double R_L(double E,double r = -1.);
        double R_L(const VecDoub &x);
        double L_E(const VecDoub &x);
        double torb(const VecDoub &x);
        double find_potential_intercept(double Phi0, int direction,double xmin,double xmax);
        VecDoub potential_flattening(double x){
        	return {
        		find_potential_intercept(Phi({x,1e-5,1e-5}),1,x/10.,x*2.)/x,
        		find_potential_intercept(Phi({x,1e-5,1e-5}),2,x/10.,x*2.)/x};
        }
        double Lzmax(double E,double RC);
};

class SphericalPotential: virtual public Potential_JS{
	private:
	public:
		inline virtual double Phi_r(double r){
			return 0.;
		}
		inline virtual double dPhi_r(double r){
			return 0.;
		}
		double Phi(const VecDoub& x){
			double r = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
			return Phi_r(sqrt(r));
		}
		VecDoub Forces(const VecDoub& x){
			double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
			double dP = dPhi_r(r)/r;
			return {-x[0]*dP,-x[1]*dP,-x[2]*dP};
		}
};

class IsochronePotential: public SphericalPotential{
	private:
		const std::string desc =
		"Isochrone potential (spherical):\n\tPhi(r) = -frac{GM}{b+sqrt(b^2+r^2}	\n\tTakes two parameters:\n\t\tG*mass: GM and the scale radius: b ";
		double GM, b;
	public:
		IsochronePotential(double GM, double b): GM(GM),b(b){};
		inline std::string name(void) const {return desc;}
		inline std::string params(void) const {
			return "GM = "+std::to_string(GM)+", b = "+std::to_string(b);
		}
		double Phi_r(double r){
			return -GM/(b+sqrt(b*b+r*r));
		}
		double dPhi_r(double r){
			double t = (b+sqrt(b*b+r*r));
			return GM*r/t/t/sqrt(b*b+r*r);
		}
		double JR(const VecDoub& x){
        	double E=H(x);
        	if(E>0.)
        	    return std::numeric_limits<double>::infinity();
        	double LL=L(x);
        	return GM/sqrt(-2*E)-0.5*(LL+sqrt(LL*LL+4.*GM*b));
		}
		VecDoub Actions(const VecDoub& x){
			return {JR(x),Lz(x),L(x)-fabs(Lz(x))};
		}
		VecDoub Omega(const VecDoub& x){
			double LL = L(x);
			double Omegar=pow(-2*H(x),1.5)/GM;
        	double Omegap=0.5*Omegar*(1.+LL/sqrt(LL*LL+4*GM*b));
        	return {Omegar,Omegap};
		}
		VecDoub Hessian(const VecDoub& x){
			VecDoub Om = Omega(x);
			double D_rr,D_rp,D_pp;
			double LL = L(x);
       		D_rr=-Om[0]*3.*sqrt(-2*H(x))/GM;
       		D_rp=D_rr*Om[1]/Om[0];
       		D_pp=Om[0]*2*GM*b*pow(LL*LL+4*GM*b,-1.5)
       				+0.5*(1+LL/sqrt(LL*LL+4*GM*b))*D_rp;
       		return {D_rr,D_rp,D_pp};
		}
};

//============================================================================
/// Axisymmetric Stackel potential with perfect ellipsoidal density
///
/// Requires an prolate spheroidal coordinate system instance
//============================================================================
class StackelOblate_PerfectEllipsoid: public Potential_JS{
	private:
		const std::string desc =
		"Axisymmetric oblate perfect ellipsoid.\nThis potential is of St\"ackel form and is most simply defined by its density:\n\trho(R,z) = -frac{rho0}{1+m^2}^2\n\tm^2 = (R/a)^2+(z/b)^2\n\tTakes two parameters:\t	the central density: rho0 and the parameter alpha such that ... hmm   ";
		std::unique_ptr<ProlateSpheroidCoordSys> CS;
		double Const;
	public:
		StackelOblate_PerfectEllipsoid(double Rho0, double alpha)
		 :CS(new ProlateSpheroidCoordSys(alpha)),
		  Const(2.*PI*conv::G*Rho0*(-CS->alpha())){}
		StackelOblate_PerfectEllipsoid(const StackelOblate_PerfectEllipsoid& s):
			CS(new ProlateSpheroidCoordSys(*s.CS)),Const(s.Const){}
		inline std::string name(void) const {return desc;}
		inline double alpha(){return CS->alpha();}
		inline void newalpha(double alp){CS->newalpha(alp);}
		inline double gamma(){return CS->gamma();}
		virtual double G(double tau);
		virtual double GPrime(double tau);
		virtual double BigFPrime(double t){std::cerr<<"NOT WRITTEN"<<std::endl; return 0.;}
		VecDoub Vderivs(const VecDoub& tau);
		inline VecDoub x2tau(VecDoub x){return CS->x2tau(x);}
		inline VecDoub xv2tau(VecDoub x){return CS->xv2tau(x);}
		inline VecDoub tau2x(VecDoub x){return CS->tau2x(x);}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
		double Phi_tau(const VecDoub& tau);
		VecDoub x2ints(const VecDoub& x, VecDoub *tau = NULL);
};

//============================================================================
///
///	Triaxial Stackel potential with perfect ellipsoidal density
///
///	Requires a confocal ellipsoidal coordinate system instance
///
//============================================================================
class StackelTriaxial: public Potential_JS{
	private:
		std::unique_ptr<ConfocalEllipsoidalCoordSys> CS;
		double a, b, c, Const;
		//, l, Flm, Elm, sinm;
		VecDoub Vderivs(const VecDoub& tau);
	public:
		StackelTriaxial(double Rho0, double alpha, double beta)
			:CS(new ConfocalEllipsoidalCoordSys(alpha,beta)),
			 a(sqrt(-CS->alpha())),
			 b(sqrt(-CS->beta())),
			 c(sqrt(-CS->gamma())),
			// l = acos(c/a);double sinL = (1-(c/a)*(c/a));
			// sinm = sqrt((1.-(b/a)*(b/a))/sinL);
			// Flm = ellint_first(l,sinm);
			// Elm = ellint_second(l,sinm);
			 Const(2.*PI*conv::G*Rho0*a*b*c*c){}
			// std::cout<<(-Const/a/a/sinL/sinm/sinm/3.*((1-sinm*sinm)*Flm+(2.*sinm*sinm-1.)*Elm-b*sinm*sinm/a*sqrt(sinL)*cos(l)))
			// 	<<" "<<-Const/a/a/sinL/sinm/sinm/(1-sinm*sinm)*(Elm-(1-sinm*sinm)*Flm-c*sinm*sinm*sin(l)/b)
			// 	<<" "<<-Const*0.5/a/a/3./sinL/(1-sinm*sinm)*(-(1-sinm*sinm)*Flm+(2.*(1-sinm*sinm)-1)*Elm+b*sin(l)*(b*b/c/c-(1-sinm*sinm))/c)<<std::endl;
		// virtual ~StackelTriaxial(){delete CS;}
		StackelTriaxial(const StackelTriaxial& s):
			CS(new ConfocalEllipsoidalCoordSys(*s.CS)),a(s.a),b(s.b),c(s.c),Const(s.Const){}
		inline double alpha(){return CS->alpha();}
		inline double beta(){return CS->beta();}
		inline double gamma(){return CS->gamma();}
		double G(double tau);
		double GPrime(double tau);
		inline VecDoub x2tau(VecDoub x){return CS->x2tau(x);}
		inline VecDoub xv2tau(VecDoub x){return CS->xv2tau(x);}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
		double Phi_tau(const VecDoub& tau);
		VecDoub tau2ints(const VecDoub& tau);
};

//============================================================================
///
///	Triaxial logarithmic potential
///
///	Phi = Vc^2/2 log(x^2+y^2/q1^2+z^2/q2^2)
///
//============================================================================
class Logarithmic: public Potential_JS{
	private:
		double Vc2, q1, q2, Phi0;
		const std::string desc =
		"Logarithmic potential:Phi ="
			"(1/2)V_c^2\\log(x^2+(y/qy)^2+(z/qz)^2)"
		"\n\tTakes three parameters:\n\t\t"
		"circular velocity: VC, and the axis ratios q_y, q_z";
	public:
		Logarithmic(double VC, double P, double Q): Vc2(VC*VC), q1(P*P),q2(Q*Q){Phi0=0.;Phi0=Phi({1e10,1e10,1e10});}
		inline std::string params(void) const {
			return "V_c^2 = "+std::to_string(Vc2)+
				   ", q_y^2 = "+std::to_string(q1)+
				   ", q_z^2 = "+std::to_string(q2);
		}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
};

//============================================================================
///
///	Triaxial power-law potential
///
///	Phi = -G M/(x^2+y^2/q1^2+z^2/q2^2)^{k/2}
///
//============================================================================
class PowerLaw: public Potential_JS{
	protected:
		double GM, k, q1, q2;
		const std::string desc =
		"Power-law potential:Phi ="
			"-GM/(x^2+y^2/q_y^2+z^2/q_z^2)^k"
		"\n\tTakes four parameters:\n\t\t"
		"G*mass: GM, the power: k and the axis ratios q_y, q_z";
	public:
		PowerLaw(double GM, double k, double P=1., double Q=1.): GM(GM), k(k),q1(P*P),q2(Q*Q){}
		inline std::string name(void) const {return desc;}
		inline std::string params(void) const {
			return "GM = "+std::to_string(GM)+", k = "+std::to_string(k)+
				   ", q_y^2 = "+std::to_string(q1)+", q_z^2 = "+std::to_string(q2);
		}
		inline void set(VecDoub pars){
			GM=pars[0];
			if(pars.size()>1) k=pars[1];
			if(pars.size()>2){ q1=pars[2];q2=pars[3];}
		}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
		double density(double r);
};
//============================================================================
///
///	Triaxial isochrone potential
///
///	Phi = -G M/(b+sqrt(b^2+x^2+y^2/q1^2+z^2/q2^2))
///
//============================================================================
class Isochrone: public Potential_JS{
	protected:
		double GM, b, q1, q2;
		const std::string desc =
		"Isochrone potential:Phi ="
			"-GM/(b+\\sqrt{b^2+r^2}) where r^2 = (x^2+y^2/q_y^2+z^2/q_z^2)"
		"\n\tTakes four parameters:\n\t\t"
		"G*mass: GM, the scale: b and the axis ratios q_y, q_z";
	public:
		Isochrone(double GM, double b, double P=1., double Q=1.): GM(GM), b(b),q1(P*P),q2(Q*Q){}
		inline std::string name(void) const {return desc;}
		inline std::string params(void) const {
			return "GM = "+std::to_string(GM)+", b = "+std::to_string(b)+
				   ", q_y^2 = "+std::to_string(q1)+", q_z^2 = "+std::to_string(q2);
		}
		inline void set(VecDoub pars){
			GM=pars[0];
			if(pars.size()>1) b=pars[1];
			if(pars.size()>2){ q1=pars[2];q2=pars[3];}
		}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
		double density(double r);
};

//============================================================================
///
///	Triaxial harmonic oscillator potential
///
///	Phi = .5*Om_i x_i x_i
///
//============================================================================
class HarmonicOscillator: public Potential_JS{
	protected:
		VecDoub Om;
		const std::string desc =
		"Triaxial harmonic oscillator potential:Phi ="
			"\\sum_i .5*\\omega_i x_i x_i"
		"\n\tTakes one parameter: a vector of the frequencies Om";
	public:
		HarmonicOscillator(VecDoub Om): Om(Om){}
		inline std::string name(void) const {return desc;}
		inline std::string params(void) const {
			std::string om = "\\omega_i = (";
			for(auto i:Om)om+=std::to_string(i)+","; om+=")";
			return om;
		}
		HarmonicOscillator(double Omx,double Omy,double Omz)
			:Om(VecDoub({Omx,Omy,Omz})){}
		double Phi(const VecDoub &x);
		VecDoub Forces(const VecDoub &x);
};

//============================================================================
///
///	Spherical Dehnen potential
///
///	rho = rhoS*pow(r,-gamma)*pow(1+pow(r,1./alpha),alpha*(beta-gamma))
///
//============================================================================
class Dehnen: public Potential_JS{
	private:
		double rhoS, alpha, beta, gamma, rs;
	public:
		Dehnen(double rhoS, double alpha, double beta, double gamma, double rs)
			: rhoS(rhoS), alpha(alpha),beta(beta), gamma(gamma), rs(rs){}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
		double Density(const VecDoub& x);
		double TotalMass(void);
};

//============================================================================
///
///	Axisymmetric Miyamoto-Nagai potential
///
///	Phi = -GM/sqrt(R^2+(A+sqrt(z^2+b^2))^2)
///
//============================================================================
class MiyamotoNagai_JS: public Potential_JS{
	private:
		double GM, A, Bq;
	public:
		MiyamotoNagai_JS(double gm, double a, double b)
			: GM(gm), A(a), Bq(b*b){}
		double Phi(const VecDoub& x);
		double Vc(double R);
		VecDoub Forces(const VecDoub& x);
};

//============================================================================
///
///	## Triaxial NFW potential
///		rs is the scale radius
///		GM is G times the mass M
///		qy and qz are the y and z flattenings respectively
///
///	 \f[\Phi = \frac{-GM}{\sqrt{x^2+y^2/q_y^2+z^2/q_z^2}}
///		\log\Big(1+\frac{\sqrt{x^2+y^2/q_y^2+z^2/q_z^2}}{R_s}\Big)\f]
///
//============================================================================
class NFW: public Potential_JS{
	private:
		const std::string desc =
		"NFW potential:Phi ="
			"-GM/sqrt{x^2+y^2/q_y^2+z^2/q_z^2}"
			"log(1+sqrt{x^2+y^2/q_y^2+z^2/q_z^2}/R_s)"
		"\n\tTakes four parameters:\n\t\tG*mass: GM, the scale radius: rs and the axis ratios q_y, q_z";
		double GM, rs, q1, q2;
	public:
		NFW(double gm, double RS, double qy, double qz): GM(gm), rs(RS), q1(qy*qy), q2(qz*qz){}
		inline std::string name(void) const {return desc;}
		inline std::string params(void) const {
			return "GM = "+std::to_string(GM)+", R_s = "+std::to_string(rs)+
				   ", q_y = "+std::to_string(q1)+", q_z = "+std::to_string(q2);
		}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
		double density(const VecDoub& x);
};

//============================================================================
///
///	## Triaxial Hernquist potential
///		rs is the scale radius
///		GM is G times the mass M
///		qy and qz are the y and z flattenings respectively
///
///	 \f[\Phi = \frac{-GM}{\sqrt{x^2+y^2/q_y^2+z^2/q_z^2}+R_s}\f]
///
//============================================================================
class Hernquist: public Potential_JS{
	private:
		const std::string desc =
		"Hernquist potential:Phi ="
			"-GM/(sqrt{x^2+y^2/q_y^2+z^2/q_z^2}+R_s)"
		"\n\tTakes four parameters:\n\t\tG*mass: GM, the scale radius: rs and the axis ratios q_y, q_z";
		double GM, rs, q1, q2;
	public:
		Hernquist(double gm, double RS, double qy, double qz): GM(gm), rs(RS), q1(qy*qy), q2(qz*qz){}
		inline std::string name(void) const {return desc;}
		inline std::string params(void) const {
			return "GM = "+std::to_string(GM)+", R_s = "+std::to_string(rs)+
				   ", q_y = "+std::to_string(q1)+", q_z = "+std::to_string(q2);
		}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
};

//============================================================================
/// Triaxial Jaffe bulge potential
//============================================================================
class Bulge: public Potential_JS{
	private:
		const std::string desc =
		"Bulge potential:Phi ="
			"GM/b_bulge*log(r/(r+b_bulge))"
		"\n\tTakes four parameters:\n\t\tG*mass: GM, the scale radius: b_bulge and the axis ratios q_y, q_z";
		double b_bulge, GM, q1, q2;
	public:
		Bulge(double b_bulge, double gm, double qy, double qz): b_bulge(b_bulge), GM(gm), q1(qy*qy), q2(qz*qz){}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
};
#ifdef TORUS
//============================================================================
/// Wrapper for potentials produced by the GalPot code
//============================================================================
class GalPot: virtual public Potential_JS{
	private:
		GalaxyPotential *PhiWD;
	public:
		GalPot(std::string TpotFile);
		inline std::string name(void) const{
			return "Galaxy Potential";
		}
		inline std::string params(void) const {
			std::ostringstream stream;
			PhiWD->DescribePot(stream);
			return stream.str()+" and mass at 100 kpc "
					+std::to_string(PhiWD->Mass(100.));
		}
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
		double Vc(double R);
		VecDoub freqs(double R);
		GalaxyPotential *PWD(){return PhiWD;}
};
#endif
//============================================================================

//============================================================================
/// Multicomponent potential
// //============================================================================
template<class P>
class MultiComponentPotential: public Potential_JS{
	private:
		std::vector<P*> multicomponents;
	public:
		MultiComponentPotential(std::vector<P*> multicomponents):multicomponents(multicomponents){};
		virtual ~MultiComponentPotential(void){};
		inline int ncompts(){return multicomponents.size();}
		double Phi(const VecDoub& x){
			double p = 0.;
			for(auto i: multicomponents)
				p+=i->Phi(x);
			return p;
		}
		VecDoub Forces(const VecDoub& x){
			VecDoub p(3,0);
			for(auto i: multicomponents) p=p+i->Forces(x);
			return p;
		}
};

template<class P>
class MultiComponentSphericalPotential: public SphericalPotential{
	private:
		std::vector<P*> multicomponents;
	public:
		MultiComponentSphericalPotential(std::vector<P*> multicomponents):multicomponents(multicomponents){};
		virtual ~MultiComponentSphericalPotential(void){};
		inline int ncompts(){return multicomponents.size();}
		double Phi_r(double r){
			double p = 0.;
			for(auto i: multicomponents)
				p+=i->Phi_r(r);
			return p;
		}
		double dPhi_r(double r){
			double p=0.;
			for(auto i: multicomponents)
				p+=i->dPhi_r(r);
			return p;
		}
};


//============================================================================

class AndreasPotential: public Potential_JS{
private:
	Bulge *bulge;
	MiyamotoNagai_JS *disc;
	NFW *halo;
	MultiComponentPotential<Potential_JS> *MC;
	const double Grav = 4.300918e-6;
public:
	AndreasPotential(NFW *H, Bulge *B = 0, MiyamotoNagai_JS *D = 0){
		std::cerr<<"NEED TO SET THE UNITS TO KPC_KMS_MSUN\n";
		if(B) bulge = B;
		else bulge = new Bulge(0.7,3.4e10*Grav,1.,1.);
		if(D) disc = D;
		else disc = new MiyamotoNagai_JS(1e11*Grav,6.5,0.26);
		if(H) halo = H;
		else halo = new NFW(Grav*1.81194e12,32.26,1.,0.8140);
		MC = new MultiComponentPotential<Potential_JS>({bulge,disc,halo});
	}
	double Phi(const VecDoub &x){return MC->Phi(x);}
	VecDoub Forces(const VecDoub &x){ return MC->Forces(x);}
	VecDoub NFWForces(const VecDoub &x){ return bulge->Forces(x);}
	VecDoub DiscForces(const VecDoub &x){ return halo->Forces(x);}
	VecDoub BulgeForces(const VecDoub &x){ return disc->Forces(x);}
};

class NFWSpherical: public SphericalPotential{
	private:
		const std::string desc =
		"Spherical NFW potential:Phi ="
			"-GM/r log(1+r/r_s)"
		"\n\tTakes two parameters:\n\t\tG*mass: GM and the scale radius: r_s";

		double GM, rs;
	public:
		NFWSpherical(double GM, double rs): GM(GM),rs(rs){}
		inline std::string name(void) const {return desc;}
		inline std::string params(void) const {
			return "GM = "+std::to_string(GM)+", r_s = "+std::to_string(rs);
		}
		inline double Phi_r(double r) {return -GM*log(1.+r/rs)/r;}
		inline double dPhi_r(double r){return GM*(log(1.+r/rs)/r-1./rs/(1.+r/rs))/r;}
};


class HernquistSpherical: public SphericalPotential{
	private:
		const std::string desc =
		"Spherical Hernquist potential:Phi ="
			"-GM/(r_s+r)"
		"\n\tTakes two parameters:\n\t\tG*mass: GM and the scale radius: r_s";
		double GM, rs;
	public:
		HernquistSpherical(double GM, double rs): GM(GM),rs(rs){}
		inline std::string name(void) const {return desc;}
		inline std::string params(void) const {
			return "GM = "+std::to_string(GM)+", r_s = "+std::to_string(rs);
		}
		inline double Phi_r(double r) {return -GM/(rs+r);}
		inline double dPhi_r(double r){return GM/(rs+r)/(rs+r);}
};


class PowerLawSpherical: public SphericalPotential{
	private:
		double GM, k;
		const std::string desc =
		"Power-law spherical potential:Phi ="
			"-GM/r^k"
		"\n\tTakes two parameters:\n\t\t"
		"G*mass: GM, the power: k";
	public:
		PowerLawSpherical(double GM, double k): GM(GM),k(k){}
		inline std::string name(void) const {return desc;}
		inline std::string params(void) const {
			return "GM = "+std::to_string(GM)+", k = "+std::to_string(k);
		}
		inline double Phi_r(double r) {return -GM*pow(r,-k);}
		inline double dPhi_r(double r){return k*GM*pow(r,-k-1);}
};

class BowdenNFW: public Potential_JS{
private:
	double rho0, rs, q0, qinf, p0, pinf;
	double rho1, rho2, r1,r2;
	const std::string desc =
	"Triaxial NFW from Bowden, Evans & Belokurov (2014)\n"
	"\tTakes 6 parameters: the usual central density (rho0) and scale (rs),"
	"and the z and y flattenings at zero(q0,p0) and inf(qinf,pinf)";
public:
	BowdenNFW(double rho00,double rs,double q0=1.,double qinf=1.,double p0=1.,double pinf=1.):rho0(rho00),rs(rs),q0(q0),qinf(qinf),p0(p0),pinf(pinf){
		double pinf3=pinf*pinf*pinf, qinf3=qinf*qinf*qinf;
		double pinf32qinf3=1+pinf3-2.*qinf3;
		double pinf3qinf3=1+pinf3+qinf3;
		rho0*=rs*rs*rs;
		r1 = rs*sqrt(2./3.*(1+p0+q0)*pinf32qinf3/(1+p0-2*q0)/pinf3qinf3);
		r2 = rs*sqrt(2./3.*(1+p0+q0)*(1+pinf3)/(1-p0)/pinf3qinf3);
		rho1 = rho0/12.*pinf32qinf3/pinf3qinf3;
		rho2 = rho0/4.*(1-pinf3)/pinf3qinf3;
	}
	inline std::string name(void) const {return desc;}
	inline std::string params(void) const {
		return "rho0 = "+std::to_string(rho0/pow(rs,3))+", rs = "+std::to_string(rs)+
		", q0 = "+std::to_string(q0)+", qinf = "+std::to_string(qinf)+
		", p0 = "+std::to_string(rho0)+", pinf = "+std::to_string(pinf);
	}
	double Phi(const VecDoub& x);
	VecDoub Forces(const VecDoub& x);
};

#ifdef TORUS
template<class c>
void vec2torus(VecDoub a, c &C){
	C[0]=a[0];C[2]=a[1];C[1]=a[2];
}

VecDoub torusPSPT2cartvec(PSPT FF);

class WrapperTorusPotential: public Potential{
	private:
	Potential_JS *Pot;
	public:
	WrapperTorusPotential(Potential_JS *pot):Pot(pot){}
	double operator()(const double R, const double z) const;
	double operator()(const double R, const double z, double& dPdR, double& dPdz) const;
    Frequencies KapNuOm(            // returns kappa,nu,Om
                const double R) const;  // given R at z=0
    double RfromLc(const double L_in, double* dR=0) const;
    double LfromRc(const double R, double* dR=0) const;

};

#endif

#ifdef GALPY
#include <galpy_potentials.h>

/*! Galpy wrapper -- adapted from Bovy */
class galpyPotential_JS : public Potential_JS {
	int nargs;
	struct potentialArg * potentialArgs;
	void  error(const char*) const;
	public:
		galpyPotential(int,struct potentialArg *);
		double Phi(const VecDoub& x);
		VecDoub Forces(const VecDoub& x);
};

inline galpyPotential::galpyPotential(int na,
				      struct potentialArg * inPotentialArgs) :
		      nargs(na), potentialArgs(inPotentialArgs)
{
parse_actionAngleArgs(nargs,potentialArgs,pot_type,pot_args,false);
}
#endif
// potential.h
#endif
