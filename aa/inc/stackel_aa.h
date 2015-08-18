#ifndef STACKEL_AA_H
#define STACKEL_AA_H

#include "potential.h"
#include "aa.h"

// ============================================================================
// Axisymmetric Stackel Perfect Ellipsoid Potential
// ============================================================================

struct root_struct_axi{
	StackelProlate_PerfectEllipsoid *P;
	VecDoub Ints;
	root_struct_axi(StackelProlate_PerfectEllipsoid *PP, VecDoub ints)
		:P(PP),Ints(ints){}
};

struct action_struct_axi{
	root_struct_axi RS;
	double taubargl, Deltagl;
	action_struct_axi(StackelProlate_PerfectEllipsoid *PP, VecDoub ints, double tb, double Dl)
		:RS(PP,ints),taubargl(tb),Deltagl(Dl){}
};

class Actions_AxisymmetricStackel : public Action_Finder{
	private:
		StackelProlate_PerfectEllipsoid *Pot;
		VecDoub find_limits(const VecDoub& x, const VecDoub& ints);
		std::vector<VecDoub> dtau01dint;
	public:
		Actions_AxisymmetricStackel(StackelProlate_PerfectEllipsoid *pot): Pot(pot){}
		VecDoub actions(const VecDoub& x, void *params=nullptr);
		VecDoub angles(const VecDoub& x, bool with_hess=false);
		void dtaudint(const VecDoub& limits, const VecDoub& ints);
		double dtaudint(const VecDoub& limits, int i, int j, double theta);
		double dDeltaGLdint(const VecDoub& limits, int i, int j);
		double dp2dtau(double tau, const VecDoub& ints);
		double BigFPrime(double t){return Pot->BigFPrime(t);}
};

struct hess_struct_axi{
	Actions_AxisymmetricStackel *ASS;
	VecDoub ints,limits;
	double taubargl, Deltagl;
	root_struct_axi RS;
	hess_struct_axi(Actions_AxisymmetricStackel *ASS, StackelProlate_PerfectEllipsoid *pot, VecDoub ints, VecDoub limits,double tb, double Dl)
		:ASS(ASS),ints(ints),limits(limits),taubargl(tb),Deltagl(Dl),RS(pot,ints){}
};

// ============================================================================
// Axisymmetric Stackel Fudge
// ============================================================================


class Actions_AxisymmetricStackel_Fudge : public Action_Finder{
	private:
		const double tiny_number = 1e-6;
		Potential_JS *Pot;
		double E, I2;
		VecDoub Kt;
		VecDoub find_limits(const VecDoub& x);
		void integrals(const VecDoub& tau);
	public:
		std::unique_ptr<OblateSpheroidCoordSys> CS;
		Actions_AxisymmetricStackel_Fudge(Potential_JS *pot,double a): Pot(pot){
			CS = std::unique_ptr<OblateSpheroidCoordSys>(new OblateSpheroidCoordSys(a));
			Kt.resize(2,0);
		}
		Actions_AxisymmetricStackel_Fudge(const Actions_AxisymmetricStackel_Fudge& a):Pot(a.Pot),CS(new OblateSpheroidCoordSys(*a.CS)){
			Kt.resize(2,0);
		}
		inline double Phi_tau(const VecDoub& tau){
			return Pot->Phi(CS->tau2x(tau));
		}
		inline double Phi_tau(double l, double n){
			return Pot->Phi(CS->tau2x({l,0.,n}));
		}
		inline double chi_lam(const VecDoub& tau){return -(tau[0]-tau[2])*Phi_tau(tau);}
		inline double chi_nu(const VecDoub& tau){ return -(tau[2]-tau[0])*Phi_tau(tau);}
		inline void reset(Potential_JS *pot){Pot = pot;}
		VecDoub actions(const VecDoub& x, void *params=nullptr);
		VecDoub angles(const VecDoub& x, void *params=nullptr);
};

struct root_struct_axi_fudge{
	Actions_AxisymmetricStackel_Fudge *ASF;
	VecDoub Ints;
	VecDoub tau_i;
	int swit;
	root_struct_axi_fudge(Actions_AxisymmetricStackel_Fudge *ASF, VecDoub ints, VecDoub tau_i, int swit)
		:ASF(ASF),Ints(ints),tau_i(tau_i),swit(swit){}
};

struct action_struct_axi_fudge{
	Actions_AxisymmetricStackel_Fudge *ASF;
	VecDoub Ints;
	VecDoub tau_i;
	double taubargl, Deltagl;
	int swit;
	action_struct_axi_fudge(Actions_AxisymmetricStackel_Fudge *ASF, VecDoub ints, VecDoub tau_i,double tb, double Dl,int swit)
		:ASF(ASF),Ints(ints),tau_i(tau_i),taubargl(tb),Deltagl(Dl),swit(swit){}
};

// ============================================================================
// Triaxial Stackel Perfect Ellipsoid Potential
// ============================================================================

struct root_struct_triax{
	StackelTriaxial *P;
	VecDoub Ints;
	root_struct_triax(StackelTriaxial *PP, VecDoub ints)
		:P(PP),Ints(ints){}
};

struct action_struct_triax{
	StackelTriaxial *P;
	VecDoub Ints;
	double taubargl, Deltagl;
	action_struct_triax(StackelTriaxial *PP, VecDoub ints, double tb, double Dl)
		:P(PP),Ints(ints),taubargl(tb),Deltagl(Dl){}
};

class Actions_TriaxialStackel : public Action_Finder{
	private:
		StackelTriaxial *Pot;
		VecDoub find_limits(const VecDoub& x,const VecDoub& ints);
	public:
		Actions_TriaxialStackel(StackelTriaxial *pot): Pot(pot){}
		VecDoub actions(const VecDoub& x0, void *params=nullptr);
};

// ============================================================================
// Triaxial Stackel Fudge
// ============================================================================

class Actions_TriaxialStackel_Fudge : public Action_Finder{
	private:
		const double tiny_number = 1e-6;
		Potential_JS *Pot;
		double E;
		VecDoub Jt, Kt;
		VecDoub find_limits(const VecDoub& x);
		void integrals(const VecDoub& tau);
		bool freq_yes;
	public:
		std::unique_ptr<ConfocalEllipsoidalCoordSys> CS;
		Actions_TriaxialStackel_Fudge(Potential_JS *pot,double a,double b): Pot(pot){
			CS = std::unique_ptr<ConfocalEllipsoidalCoordSys>
				(new ConfocalEllipsoidalCoordSys(a,b));
			Jt.resize(3,0);Kt.resize(3,0);
			freq_yes=false;
		}
		Actions_TriaxialStackel_Fudge(const Actions_TriaxialStackel_Fudge& s):Pot(s.Pot),
			CS(new ConfocalEllipsoidalCoordSys(*s.CS)){
			Jt.resize(3,0);Kt.resize(3,0);
			freq_yes=false;
		}
		inline double Phi_tau(const VecDoub& tau){
			return Pot->Phi(CS->tau2x(tau));
		}
		inline double Phi_tau(double l, double m, double n){
			return Pot->Phi(CS->tau2x({l,m,n}));
		}
		inline double chi_lam(const VecDoub& tau){
			return (tau[0]-tau[1])*(tau[2]-tau[0])*Phi_tau(tau);
		}
		inline double chi_mu(const VecDoub& tau){
			return (tau[1]-tau[2])*(tau[0]-tau[1])*Phi_tau(tau);
		}
		inline double chi_nu(const VecDoub& tau){
			return (tau[2]-tau[0])*(tau[1]-tau[2])*Phi_tau(tau);
		}
		double ptau_tau(double tau, const VecDoub& ints, const VecDoub& tau_i, int swit);

		inline void reset(Potential_JS *pot){Pot = pot;}
		VecDoub actions(const VecDoub& x, void *params=nullptr);
		VecDoub angles(const VecDoub& x, void *params=nullptr);
		double sos(const VecDoub& x, int comp,const std::string& outfile);
		inline void set_freq(bool s){freq_yes=s;}
};


struct root_struct_triax_fudge{
	Actions_TriaxialStackel_Fudge *ATSF;
	VecDoub Ints;
	VecDoub tau_i;
	int swit;
	root_struct_triax_fudge(Actions_TriaxialStackel_Fudge *ATSF, VecDoub ints, VecDoub tau_i, int swit)
		:ATSF(ATSF),Ints(ints),tau_i(tau_i),swit(swit){}
};

struct action_struct_triax_fudge{
	Actions_TriaxialStackel_Fudge *ATSF;
	VecDoub Ints;
	VecDoub tau_i;
	double taubargl, Deltagl;
	int swit;
	action_struct_triax_fudge(Actions_TriaxialStackel_Fudge *ATSF, VecDoub ints, VecDoub tau_i,double tb, double Dl,int swit)
		:ATSF(ATSF),Ints(ints),tau_i(tau_i),taubargl(tb),Deltagl(Dl),swit(swit){}
};
#endif
// ============================================================================
