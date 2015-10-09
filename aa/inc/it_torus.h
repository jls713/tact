// ============================================================================
/// \file inc/it_torus.h
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
/// \brief Action finding by iteratively constructing tori
///
///	We estimate the actions by iteratively estimating them using an approximate
/// method (here the Actions_AxisymmetricStackel_Fudge) and then constructing
/// a torus with this action and finding the error made in the action estimate.
//============================================================================

#ifndef IT_TORUS_H
#define IT_TORUS_H

#include "potential.h"
#include "stackel_aa.h"
#include "Torus.h"

//============================================================================
/*! Action finding by iteratively constructing tori*/
class IterativeTorusMachine : public Action_Finder{
	private:
		Actions_AxisymmetricStackel_Fudge *AAA;/*!< Approximate action finder*/
		Potential_JS *Phi;					   /*!< Target potential 		 */
		WrapperTorusPotential TPhi;		  /*!< Wrapper of potential for Torus*/
		Torus    *T;
		double eta;								   /*!< tolerance in ang fit */
		int MaxIterations;				 /*!< max no. of torus constructions */
		double dJ;		   /*!< relative action error for torus construction */


		VecDoub RefineGuess(VecDoub PSP, VecDoub aa, double *min, Angles &theta_new, double dJ);
		VecDoub NewActionPoint(VecDoub startpt, VecDoub oldpt, VecDoub newpt);
		VecDoub find_aa(VecDoub x);
		public:
			//! IterativeTorusMachine constructor.
		    /*!
		      \param AAA Axisymmetric fudge action finder
		      \param Phi Potential_JS (axisymmetric)
		      \param eta  tolerance in angle fit
		      \param MaxIterations  maximum number of torus fits
		      \param dJ relative error of torus fits
		    */
			IterativeTorusMachine(Actions_AxisymmetricStackel_Fudge *AAA,Potential_JS *Phi, double eta = 1e-8, int MaxIterations=5, double dJ = 1e-5):AAA(AAA),Phi(Phi),TPhi(Phi),eta(eta),MaxIterations(MaxIterations),dJ(dJ),T(new Torus){};
			//! IterativeTorusMachine destructor.
			~IterativeTorusMachine(){delete T;}
			//!< set tolerance in ang fit
			inline void set_eta(double s){eta = s;}
			//!< set max no. of torus constructions
			inline void set_maxit(double s){MaxIterations = s;}
			//!< set relative action error for torus construction
			inline void set_dJ(double s){dJ = s;}
			//!< set tolerance, max no. of tori and relative error
			/*!
		      \param a tolerance in ang. fit
		      \param b max no. of torus constructions
		      \param c relative action error for torus construction
		    */
			inline void set_options(double a, double b, double c)
				{eta=a;MaxIterations=b;dJ=c;}
			//! Finds actions
		    /*!
		      \param x phase-space point (x,v)
		      \param params -- does nothing

		      \return actions -- 3D vector J=(J_R,J_phi,J_z)
		    */
			VecDoub actions(const VecDoub& x, void*params=nullptr);
			//! Finds angles
		    /*!
		      \param x phase-space point (x,v)
		      \param params -- does nothing

		      \return angles and frequencies --
		      	6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
		    */
			VecDoub angles(const VecDoub& x, void*params=nullptr);
};

VecDoub findXV(Angles theta, Torus *T, double sign);

/*! Helper structure for passing parameters to minimization routine when finding nearest angle point in torus */
struct Min_ItTorusMac_struct{
	Torus *T;VecDoub XV; double freqWeight;double sign;
	Min_ItTorusMac_struct(Torus *T, VecDoub XV,double freqWeight, double sign):T(T), XV(XV),freqWeight(freqWeight),sign(sign){}
};
//============================================================================

#endif
