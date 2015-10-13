// ============================================================================
/// \file coordtransform/inc/coordsys.h
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
/*! \brief Coordinate systems
 *
 *  1. Implements prolate spheroidal coordinate systems (UVProlateSpheroidCoordSys and ProlateSpheroidCoordSys)
 *  2. Implements full ellipsoidal coordinate system
 */
//=============================================================================

#ifndef COORDSYS_H
#define COORDSYS_H

//=============================================================================

#include "utils.h"

//=============================================================================
/*! Prolate spheroidal coordinate system using u, v coordinates */
class UVProlateSpheroidCoordSys{
	// (u,v) coords
	private:
		double Delta, Delta2; /*! focal length and its square */
	public:
		//! UVProlateSpheroidCoordSys constructor.
		/*!
			\param Delta -- focal length
		*/
		UVProlateSpheroidCoordSys(double Delta): Delta(Delta),Delta2(Delta*Delta){}
		/*! output focal length */
		inline double delta(void){ return Delta;}
		//! convert x to prolate spheroidal coordinates
		/*!
			\param x -- input 3D Cartesian coordinate x=(x,y,z)
			\return uv = (u,phi,v)
		*/
		VecDoub xv2uv(const VecDoub& x);
		//! convert prolate spheroidal coordinates to Cartesian
		/*!
			\param uv = (u,phi,v)
			\return 3D Cartesian coordinate x=(x,y,z)
		*/
		VecDoub uv2Rz(const VecDoub& x);
};
/*! Prolate spheroidal coordinate system using lambda, nu coordinates
	Delta = sqrt(Gamma-Alpha)
*/
class ProlateSpheroidCoordSys{
	// (lam,nu) coords
	private:
		double Alpha; 		/*! coordinate system parameter Alpha 		*/
		const double Gamma; /*! coordinate system parameter Gamma = -1. */
	public:
		//! ProlateSpheroidCoordSys constructor.
		/*!
			\param Alpha -- coord sys param. focal length, Delta = sqrt(Gamma-Alpha)
		*/
		ProlateSpheroidCoordSys(double alpha): Alpha(alpha), Gamma(-1.){}
		inline double alpha(void){ return Alpha;}
		inline double gamma(void){ return Gamma;}
		/*! reset Alpha */
		inline void newalpha(double a){ Alpha = a;}
		VecDoub x2tau(const VecDoub& x);
		//! convert prolate spheroidal coordinates to Cartesian
		/*!
			\param x = (lambda,phi,nu)
			\return 3D Cartesian coordinate x=(x,y,z)
		*/
		VecDoub tau2x(const VecDoub& x);
		//! convert prolate spheroidal coordinates to cylindrical polar
		/*!
			\param x = (lambda,phi,nu)
			\return 3D cylindrical polar coordinate x=(R,phi,z)
		*/
		VecDoub tau2polar(const VecDoub& x);
		//! convert Cartesian to tau full 6D coords
		/*!
			\param x = (x,y,z,vx,vy,vz)
			\return 6D prolate spheroidal coordinates (lambda,phi,nu) and time derivs
		*/
		VecDoub xv2tau(const VecDoub& x);
		//! Calculates tau & derivatives of tau wrt to R and z
		/*!
			\param 3D Cartesian position x = (x,y,z)
			\return {lamda,nu,dlambdadR, dlambdadz, dnudR, dnudz};
		*/
		VecDoub derivs(const VecDoub& x);
		//! calculates tau momenta given tau and taudot
		/*!
			\param tau = (lambda,phi,nu,lambdadot,phidot,nudot)
			\return (p_lambda,p_nu)
		*/
		VecDoub tau2p(const VecDoub& tau);
};

//=============================================================================
/*! Ellipsoidal coordinate system using lambda, mu, nu coordinates
	Delta_1 = sqrt(Beta-Alpha), Delta_2 = sqrt(Gamma-Beta)
*/
class ConfocalEllipsoidalCoordSys{
	// (lam,mu,nu) coords
	private:
		double Alpha, Beta, Gamma;	/*!< Coordinate system parameters */
	public:
		//! ConfocalEllipsoidalCoordSys constructor.
		/*!
			\param Alpha -- coord sys param
			\param Beta -- coord sys param
			focal lengths, Delta_1 = sqrt(Beta-Alpha), Delta_2 = sqrt(Gamma-Beta)
		*/
		ConfocalEllipsoidalCoordSys(double a, double b): Alpha(a), Beta(b), Gamma(-1.){}
		inline double alpha(void){ return Alpha;}
		inline double beta(void){ return Beta;}
		inline double gamma(void){ return Gamma;}
		/*! reset Alpha */
		inline void newalpha(double a){ Alpha = a;}
		/*! reset Beta */
		inline void newbeta(double b){ Beta = b;}
		//! convert Cartesian to ellipsoidal coordinates
		/*!
			\param x = (x,y,z)
			\return (lambda,mu,nu)
		*/
		VecDoub x2tau(const VecDoub& x);
		//! convert ellipsoidal coordinates to Cartesian
		/*!
			\param x = (lambda,mu,nu)
			\return (x,y,z)
		*/
		VecDoub tau2x(const VecDoub& tau);
		//! convert Cartesian to tau full 6D coords
		/*!
			\param x = (x,y,z,vx,vy,vz)
			\return 6D prolate spheroidal coordinates (lambda,mu,nu) and time derivs
		*/
		VecDoub xv2tau(const VecDoub& x);
		//! calculates tau momenta given tau and taudot
		/*!
			\param tau = (lambda,mu,nu,lambdadot,mudot,nudot)
			\return (p_lambda,p_mu,p_nu)
		*/
		VecDoub tau2p(const VecDoub& tau);
		//! Calculates tau & derivatives of tau wrt to R and z
		/*!
			\param 3D Cartesian position x = (x,y,z)
			\return (dl/dx,dm/dx,dn/dx,dl/dy...
		*/
		VecDoub derivs(const VecDoub& x);
};

//=============================================================================

#endif

//=============================================================================
