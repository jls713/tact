// ============================================================================
/// \file inc/falcON_interface.h
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
/// \brief Wrapper of Walter Dehnen's falPot potentials
///
/// The nemo accelerations available through Walter Dehnen's falPot code
/// are wrapped in a Potential_JS object such that they can be passed to
/// Action_Finder classes.
///
// ============================================================================

#ifndef FALCON_INTERFACE_H
#define FALCON_INTERFACE_H
#include "potential.h"
#include "public/basic.h"
#include "public/types.h"
#include <body.h>
#include <acceleration.h>
#include <externacc.h>

// ============================================================================
/*! Wrapper for falcON potential/acceleration classes such that they can be used by action-finding classes */
class falcON_Potential :
	public Potential_JS{
	private:
		falcON::nemo_acc *falPot;/*! Class for computing the falPot potential*/
	public:
		//! falcON_Potential constructor.
	    /*!
	      \param accname name of falPot acceleration class
	      \param accpars params for falPot acc. class
	      \param accfile file with params for falPot acc. class
	    */
		falcON_Potential(const char*accname,
	                 	 const char*accpars,
	                 	 const char*accfile){
			falPot = new falcON::nemo_acc(accname,accpars,accfile);
		}
		//! falcON_Potential destructor.
		~falcON_Potential(){delete falPot;}
		//! computes potential
	    /*!
	      \param x 3D vector position x = (x,y,z)

	      \return potential Phi
	    */
		double Phi(const VecDoub &x){
			double p=0.,m=0.;
			falcON::falcONVec<3,double> xx, vv, aa;
			for(int j=0;j<3;j++){
				xx[j]=x[j];
				vv[j]=0.;
				aa[j]=0.;
			}
			falPot->set(0.,1,&m,&xx,&vv,nullptr,&p,&aa,0);
			return p;
		}
		//! computes forces
	    /*!
	      \param x 3D vector position x = (x,y,z)

	      \return 3D Cartesian forces = -(dPhi/dx,dPhi/dy,dPhi/dz)
	    */
		VecDoub Forces(const VecDoub &x){
			VecDoub a(3,0); double p=0.,m=0.;
			falcON::falcONVec<3,double> xx, vv, aa;
			for(int j=0;j<3;j++){xx[j]=x[j];vv[j]=0.;aa[j]=0.;}
			falPot->set(0.,1,&m,&xx,&vv,nullptr,&p,&aa,0);
			for(int j=0;j<3;j++) a[j]=aa[j];
			return a;
		}
};

// ============================================================================

/*! Wrapper for falcON potential/acceleration classes such that they can be used by spherical action-finding classes */
class falcON_SphericalPotential :
	public SphericalPotential{
	private:
		falcON::nemo_acc *falPot;/*! Class for computing the falPot potential*/
	public:
		//! falcON_SphericalPotential constructor.
	    /*!
	      \param accname name of falPot acceleration class
	      \param accpars params for falPot acc. class
	      \param accfile file with params for falPot acc. class
	    */
		falcON_SphericalPotential(const char*accname,
	                 	 const char*accpars,
	                 	 const char*accfile){
			falPot = new falcON::nemo_acc(accname,accpars,accfile);
		}
		//! falcON_SphericalPotential destructor.
		~falcON_SphericalPotential(){delete falPot;}
		//! computes potential
	    /*!
	      \param r spherical radius r = sqrt(x^2+y^2+z^2)

	      \return potential Phi
	    */
		double Phi_r(double r){
			double p=0.,m=0.;
			falcON::falcONVec<3,double> xx, vv, aa;
			xx[0]=r;xx[1]=0.;xx[2]=0.;
			for(int j=0;j<3;j++){
				vv[j]=0.;
				aa[j]=0.;
			}
			falPot->set(0.,1,&m,&xx,&vv,nullptr,&p,&aa,0);
			return p;
		}
		//! computes potential deriv
	    /*!
	      \param r spherical radius r = sqrt(x^2+y^2+z^2)

	      \return potential deriv dPhi/dr
	    */
		double dPhi_r(double r){
			VecDoub a(3,0); double p=0.,m=0.;
			falcON::falcONVec<3,double> xx, vv, aa;
			xx[0]=r;xx[1]=0.;xx[2]=0.;
			for(int j=0;j<3;j++){
				vv[j]=0.;
				aa[j]=0.;
			}
			falPot->set(0.,1,&m,&xx,&vv,nullptr,&p,&aa,0);
			for(int j=0;j<3;j++) a[j]=aa[j];
			return -a[0];
		}
};
#endif
// FALCON_INTERFACE_H
// ============================================================================
