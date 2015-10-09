// ============================================================================
/// \file inc/uv_orb.h
// ============================================================================
/// \author Jason Sanders (and James Binney 2012)
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
/// \brief Axisymmetric Staeckel fudge using Delta estimation from shells
///
/// uv_orb: Wraps axisymmetric Staeckel fudge using the Delta estimation
/// routine from Binney (2014)
///
//============================================================================

#ifndef uv_H
#define uv_H
#include "aa.h"

// ============================================================================
/*! uv_orb is an interface for using axisymmetric St\"ackel fudge apparatus.
	It handles the choice of coordinate system using the Delta estimation from
	Binney (2014). It is thread-safe as each action call creates a new instance
	of Actions_AxisymmetricStackel_Fudge */
class uv_orb: public Action_Finder{

	private:
		Potential_JS *Pot;  /*!< Potential (axisymmetric)					*/
		VecDoub E_delta;	/*!< Energy grid	`							*/
		std::vector<VecDoub> L_delta;/*!< Ang mom. grid	`					*/
		std::vector<VecDoub> Delta;/*< Delta grid	`						*/
		double Rmin; 		/*!< Minimum R intercept 						*/
		double Rmax; 		/*!< Maximum R intercept 						*/
		double Emax = 0.; 	/*!< Maximum energy								*/
		double E0 = 0.;		/*!< Core energy     							*/
		int NE;				/*!< Number of energy grid points				*/
		int NL;				/*!< Number of radial grid points 				*/
		std::string name;   /*!< Filename for storing interp grids of Delta */
	public:
		//! uv_orb constructor.
	    /*!
	      \param pot Potential_JS (axisymmetric)
	      \param Rm  minimum radial grid point
	      \param Rn  maximum radial grid point
	      \param NE  number of energy grid points
	      \param NL  number of ang. mom. grid points
	      \param name name of class
	    */
		uv_orb(Potential_JS *pot, double Rm=0.05, double Rn=80.,int NE=20, int NL = 10,std::string name="")
			:Pot(pot),Rmin(Rm),Rmax(Rn),NE(NE),NL(NL),name(name){
			E0 = Pot->E_circ(Rmin);
			Emax = Pot->E_circ(Rmax);
			fillDeltagrids(name!=""?name+".Delta_uv":"Delta_uv.tmp");
		}
		//! reset uv_orb.
	    /*!
	      Rebuilds delta grid for new potential
	      \param pot Potential_JS (axisymmetric)
	    */
		void reset(Potential_JS *pot){
			Pot = pot;
			Emax = Pot->E_circ(Rmax);
			E0 = Pot->E_circ(Rmin);
			fillDeltagrids(name!=""?name+".Delta_uv":"Delta_uv.tmp");
		}
		//! Fill delta grid.
	    /*!
	      Builds delta grid and outputs to file
	      \param f output file name string
	    */
		void fillDeltagrids(const std::string& f);
		//! Read in delta grid.
	    /*!
	      Reads delta grid from file
	      \param f input file name string
	    */
		void readDeltagrids(const std::string& f);
		//! Finds Delta
	    /*!
	      Interpolates Delta grids.
	      \param E energy
	      \param L angular momentum

	      \return Delta
	    */
		double findDelta_interp(double E,double L);
		//! Finds actions
	    /*!
	      \param x phase-space point (x,v)
	      \param params -- can pass alpha value to override Delta estimation

		  \return actions -- 3D vector J=(J_R,J_phi,J_z)
	    */
		VecDoub actions(const VecDoub& x,void*params=nullptr);
		//! Finds angles
	    /*!
	      \param x phase-space point (x,v)
	      \param params -- can pass alpha value to override Delta estimation

		  \return angles and frequencies --
		      	6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
	    */
		VecDoub angles(const VecDoub& x,void*params=nullptr);
};

#endif
// ============================================================================
