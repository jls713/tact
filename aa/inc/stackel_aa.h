// ============================================================================
/// \file inc/stackel_aa.h
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
/// \brief Action finding in Staeckel potentials and Staeckel fudges
///
/// Four classes are implemented:
/// 1. Actions_AxisymmetricStackel: Action finding in axisymmetric Staeckel
///    potentials (currently accepts perfect ellipsoid potential)
/// 2. Actions_TriaxialStackel: Action finding in triaxial Staeckel potential
/// 3. Actions_AxisymmetricStackel_Fudge: Action estimation in axisymmetric
///    potential using Staeckel fudge (as in Binney (2012))
/// 4. Actions_TriaxialStackel_Fudge :Action estimation in triaxial potential
///    using Staeckel fudge (as in Sanders & Binney (2014))
///
//============================================================================

#ifndef STACKEL_AA_H
#define STACKEL_AA_H

#include "potential.h"
#include "aa.h"

//============================================================================
/*! Action finding in axisymmetric Staeckel potential */
class Actions_AxisymmetricStackel : public Action_Finder{
	private:
		StackelOblate_PerfectEllipsoid *Pot; /*< Staeckel potential */
		std::vector<VecDoub> dtau01dint;

		VecDoub find_limits(const VecDoub& x, const VecDoub& ints);
	public:
		//! Actions_AxisymmetricStackel constructor.
	    /*!
	      \param pot StackelOblate_PerfectEllipsoid potential in which to
	      compute the actions
	    */
		Actions_AxisymmetricStackel(StackelOblate_PerfectEllipsoid *pot): Pot(pot){}
		//! Finds actions
	    /*!
	      \param x phase-space point (x,v)
	      \param params -- doesn't do anything

		  \return actions -- 3D vector J=(J_R,J_phi,J_z)
	    */
		VecDoub actions(const VecDoub& x, void *params=nullptr);
		//! Finds angles
	    /*!
	      \param x phase-space point (x,v)
	      \param with_hess -- option to calculate hessian

		  \return angles and frequencies --
		      	6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
	    */
		VecDoub angles(const VecDoub& x, void *params=nullptr);

		// Unimportant functions that need to be accessed by integration routines
		void dtaudint(const VecDoub& limits, const VecDoub& ints);
		double dtaudint(const VecDoub& limits, int i, int j, double theta);
		double dDeltaGLdint(const VecDoub& limits, int i, int j);
		double dp2dtau(double tau, const VecDoub& ints);
		double BigFPrime(double t){return Pot->BigFPrime(t);}
};

/*! Helper structure for finding limits of action integrals for axisymmetric
	Stackel */
struct root_struct_axi{
	StackelOblate_PerfectEllipsoid *P;
	VecDoub Ints;
	root_struct_axi(StackelOblate_PerfectEllipsoid *PP, VecDoub ints)
		:P(PP),Ints(ints){}
};

/*! Helper structure for integration of actions for axisymmetric Stackel */
struct action_struct_axi{
	root_struct_axi RS;
	double taubargl, Deltagl, tiny_number;
	action_struct_axi(StackelOblate_PerfectEllipsoid *PP, VecDoub ints, double tb, double Dl, double tn)
		:RS(PP,ints),taubargl(tb),Deltagl(Dl), tiny_number(tn){}
};

/*! Helper structure for finding Hessian for axisymmetric Stackel */
struct hess_struct_axi{
	Actions_AxisymmetricStackel *ASS;
	VecDoub ints,limits;
	double taubargl, Deltagl;
	root_struct_axi RS;
	hess_struct_axi(Actions_AxisymmetricStackel *ASS, StackelOblate_PerfectEllipsoid *pot, VecDoub ints, VecDoub limits,double tb, double Dl)
		:ASS(ASS),ints(ints),limits(limits),taubargl(tb),Deltagl(Dl),RS(pot,ints){}
};

//============================================================================
/*! Action estimation in general axisymmetric potential using Staeckel fudge */
class Actions_AxisymmetricStackel_Fudge : public Action_Finder{
	private:
		Potential_JS *Pot;

		const double tiny_number = 1e-10;/*!< Tolerance for \int 1/p_tau    */
		double E, I2; 					/*!< Energy, I_2 = 0.5 L_z^2       */
		VecDoub Kt;	  					/*!< Third integrals 			   */
		VecDoub find_limits(const VecDoub& x);/*!<Find tau limits 		   */
		void integrals(const VecDoub& tau);	/*!< Find E, I2 and Kt 		   */
	public:
		std::unique_ptr<ProlateSpheroidCoordSys> CS;/*!< Coordinate system   */
		//! Actions_AxisymmetricStackel_Fudge constructor.
	    /*!
	      \param pot Potential (axisymmetric) in which to compute the actions
	      \param a   alpha value to use for coordinate system (will be
	      overridden by actions and angles if required)
	    */
		Actions_AxisymmetricStackel_Fudge(Potential_JS *pot,double a): Pot(pot){
			CS = std::unique_ptr<ProlateSpheroidCoordSys>(new ProlateSpheroidCoordSys(a));
			Kt.resize(2,0);
		}
		//! Actions_AxisymmetricStackel_Fudge copy constructor.
		Actions_AxisymmetricStackel_Fudge(const Actions_AxisymmetricStackel_Fudge& a):Pot(a.Pot),CS(new ProlateSpheroidCoordSys(*a.CS)){
			Kt.resize(2,0);
		}
		inline void reset(Potential_JS *pot){Pot = pot;}
		//! Compute potential at tau coordinate
		inline double Phi_tau(const VecDoub& tau){
			return Pot->Phi(CS->tau2x(tau));
		}
		//! Compute potential at tau =(l,0,n)
		inline double Phi_tau(double l, double n){
			return Pot->Phi(CS->tau2x({l,0.,n}));
		}
		inline double chi_lam(const VecDoub& tau){return -(tau[0]-tau[2])*Phi_tau(tau);}
		inline double chi_nu(const VecDoub& tau){ return -(tau[2]-tau[0])*Phi_tau(tau);}
		//! Finds actions
	    /*!
	      \param x phase-space point (x,v)
	      \param params -- if null, use alpha specified in constructor
	      				-- if >0, estimate using derivatives of potential eq(8) Sanders (2012)
	      				-- if <0, use value passed as new alpha

		  \return actions -- 3D vector J=(J_R,J_phi,J_z)
	    */
		VecDoub actions(const VecDoub& x, void *params=nullptr);
		//! Finds angles
	    /*!
	      \param x phase-space point (x,v)
	      \param params -- if null, use alpha specified in constructor
	      				-- if >0, estimate using derivatives of potential eq(8) Sanders (2012)
	      				-- if <0, use value passed as new alpha

		  \return angles and frequencies --
		      	6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
	    */
		VecDoub angles(const VecDoub& x, void *params=nullptr);
};

/*! Helper structure for root-finding for axisymmetric Stackel fudge */
struct root_struct_axi_fudge{
	Actions_AxisymmetricStackel_Fudge *ASF;
	VecDoub Ints;
	VecDoub tau_i;
	int swit;
	root_struct_axi_fudge(Actions_AxisymmetricStackel_Fudge *ASF, VecDoub ints, VecDoub tau_i, int swit)
		:ASF(ASF),Ints(ints),tau_i(tau_i),swit(swit){}
};

/*! Helper structure for action integrals for axisymmetric Stackel fudge */
struct action_struct_axi_fudge{
	Actions_AxisymmetricStackel_Fudge *ASF;
	VecDoub Ints;
	VecDoub tau_i;
	double taubargl, Deltagl;
	int swit;
	double tiny_number;
	action_struct_axi_fudge(Actions_AxisymmetricStackel_Fudge *ASF, VecDoub ints, VecDoub tau_i,double tb, double Dl,int swit, double tn)
		:ASF(ASF),Ints(ints),tau_i(tau_i),taubargl(tb),Deltagl(Dl),swit(swit), tiny_number(tn){}
};

//============================================================================
/*! Action finding in triaxial Staeckel potential */
class Actions_TriaxialStackel : public Action_Finder{
	private:
		StackelTriaxial *Pot; /*< Staeckel potential */
		VecDoub find_limits(const VecDoub& x,const VecDoub& ints);
	public:
		//! Actions_TriaxialStackel constructor.
	    /*!
	      \param pot Triaxial Stackel Potential (axisymmetric)
	    */
		Actions_TriaxialStackel(StackelTriaxial *pot): Pot(pot){}
		//! Finds actions
	    /*!
	      \param x phase-space point (x,v)
	      \param params -- if not null, returns frequencies

		  \return actions -- 3D vector J=(J_R,J_phi,J_z)
	    */
		VecDoub actions(const VecDoub& x0, void *params=nullptr);
};

/*! Helper structure for root-finding for triaxial Staeckel actions */
struct root_struct_triax{
	StackelTriaxial *P;
	VecDoub Ints;
	root_struct_triax(StackelTriaxial *PP, VecDoub ints)
		:P(PP),Ints(ints){}
};

/*! Helper structure for action integrals for triaxial Staeckel actions */
struct action_struct_triax{
	StackelTriaxial *P;
	VecDoub Ints;
	double taubargl, Deltagl;
	action_struct_triax(StackelTriaxial *PP, VecDoub ints, double tb, double Dl)
		:P(PP),Ints(ints),taubargl(tb),Deltagl(Dl){}
};

//============================================================================
/*! Action estimation in general triaxial potentials using Staeckel fudge */
class Actions_TriaxialStackel_Fudge : public Action_Finder{
	private:
		Potential_JS *Pot;
		const double tiny_number = 1e-6;	/*!< Tolerance for \int 1/p_tau */
		double E; 							/*!< Energy 				    */
		VecDoub Jt, Kt;						/*!< Second and third integrals */
		VecDoub find_limits(const VecDoub& x);/*!<Find tau limits 		    */
		void integrals(const VecDoub& tau);	/*!< Find E, Jt and Kt 		    */
		bool freq_yes;						/*!< Returns freq from actions  */
	public:
		std::unique_ptr<ConfocalEllipsoidalCoordSys> CS;/*!<coordinate system*/
		//! Actions_TriaxialStackel_Fudge constructor.
	    /*!
	      \param pot Potential_JS (triaxial)
	      \param a   alpha for coordinate system
	      \param b   beta for coordinate system
	    */
		Actions_TriaxialStackel_Fudge(Potential_JS *pot,double a,double b): Pot(pot){
			CS = std::unique_ptr<ConfocalEllipsoidalCoordSys>
				(new ConfocalEllipsoidalCoordSys(a,b));
			Jt.resize(3,0);Kt.resize(3,0);
			freq_yes=false;
		}
		//! Actions_TriaxialStackel_Fudge copy constructor.
		Actions_TriaxialStackel_Fudge(const Actions_TriaxialStackel_Fudge& s):Pot(s.Pot),
			CS(new ConfocalEllipsoidalCoordSys(*s.CS)){
			Jt.resize(3,0);Kt.resize(3,0);
			freq_yes=false;
		}
		//! reset -- change potential
		inline void reset(Potential_JS *pot){Pot = pot;}
		//! return freq from actions or not
		inline void set_freq(bool s){freq_yes=s;}
		//! computes potential at ellipsoidal coordinates tau
		inline double Phi_tau(const VecDoub& tau){
			return Pot->Phi(CS->tau2x(tau));
		}
		//! computes potential at ellipsoidal coordinates tau=(l,m,n)
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
		//! Finds actions
	    /*!
	      \param x phase-space point (x,v)
	      \param params -- does nothing

		  \return actions -- 3D vector J=(J_R,J_phi,J_z)
	    */
		VecDoub actions(const VecDoub& x, void *params=nullptr);
		//! Finds angles
	    /*!
	      \param x phase-space point (x,v)
	      \param params -- does nothing

		  \return angles and frequencies --
		      	6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
	    */
		VecDoub angles(const VecDoub& x, void *params=nullptr);
		//! computes surface of section
	    /*!
	      \param x phase-space point (x,v)
	      \param comp -- if 0 use x, if 1 use y, if 2 use z
	      \param outfile -- output file

		  \return nothing
	    */
		double sos(const VecDoub& x, int comp,const std::string& outfile);
};

/*! Helper structure for root-finding for triaxial Staeckel fudge */
struct root_struct_triax_fudge{
	Actions_TriaxialStackel_Fudge *ATSF;
	VecDoub Ints;
	VecDoub tau_i;
	int swit;
	root_struct_triax_fudge(Actions_TriaxialStackel_Fudge *ATSF, VecDoub ints, VecDoub tau_i, int swit)
		:ATSF(ATSF),Ints(ints),tau_i(tau_i),swit(swit){}
};
/*! Helper structure for action integrals for triaxial Staeckel fudge */
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
