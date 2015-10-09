// ============================================================================
/// \file inc/lmn_orb.h
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
/// \brief Axisymmetric Staeckel fudge using Delta estimation from shells
///
/// lmn_orb: Wraps triaxial Staeckel fudge using the alpha, beta (or Delta_1
/// and Delta_2) estimation routine from Sanders & Binney (2014).
/// We find the closed loop orbits in the (x,y) and (y,z) planes and fit
/// ellipses to these orbits to find alpha and beta respectively.
///
//============================================================================

#ifndef LMN_H
#define LMN_H
#include "aa.h"
// ============================================================================
// Lambda_Mu_Nu orbit action finding
// ============================================================================

/*! Interface for using triaxial St\"ackel fudge apparatus
	Handles the choice of coordinate system
*/
class lmn_orb: public Action_Finder{
	private:
		VecDoub E;  						/*!< Energy grid				 */
		VecDoub Beta;						/*!< Beta grid				     */
		VecDoub Alpha;						/*!< Alpha grid				     */
		Potential_JS *Pot;  				/*!< Potential (triaxial)	     */

		double ymin; 		/*!< Maximum y intercept*/
		double ymax; 		/*!< Maximum y intercept*/
		int NE;				/*!< Size of grid 		*/
		double Emax = 0.; 	/*!< Maximum energy 	*/
		double E0 = 0.;		/*!< Core energy 		*/
		std::string name;	/*!< output file name 	*/
		bool use_log_grid;	/*!< use log-spaced grid*/
		bool use_acts;		/*!< estimate alpha, beta by minimising actions */
	public:
        //! lmn_orb constructor.
        /*!
          \param pot Potential (axisymmetric) in which to compute the actions
          \param ym -- minimum intermediate (y) intercept value
          \param yn -- maximum intermediate (y) intercept value
          \param NEE -- number of energy grid points
          \param use_log_grid -- if true, use log-spaced grid, else linear
          \param use_acts-if true, estimate alpha & beta by action minimization
          \param s -- output file
        */
		lmn_orb(Potential_JS *pot, double ym=0.05, double yn=80.,int NEE=32,bool use_log_grid=false,bool use_acts=false,std::string s="")
			:Pot(pot),ymin(ym),ymax(yn),NE(NEE),
			 name(s),use_log_grid(use_log_grid),use_acts(use_acts){
			Emax = Pot->Phi({0.,ymax,0.});
			E0 = Pot->Phi({0.,ymin,0.});
			fillDeltagrids(name!=""?name+".Delta_lmn":"Delta_lmn.tmp");
		}
		void reset(Potential_JS *pot){
			Pot = pot;
			Emax = Pot->Phi({0.,ymax,0.});
			E0 = Pot->Phi({0.,ymin,0.});
			fillDeltagrids(name!=""?name+".Delta_lmn":"Delta_lmn.tmp");
		}
		// inline void newbeta(double beta){ATSF->CS->newbeta(beta);}
		// inline void newalpha(double alpha){ATSF->CS->newalpha(alpha);}

        //! Finds actions
        /*!
          \param x phase-space point (x,v)
          \param params -- can pass pointer to VecDoub ab=(alpha,beta)
          \return actions -- 3D vector J=(J_R,J_phi,J_z)
        */
		VecDoub actions(const VecDoub& x,void*params=nullptr);
        //! Finds angles
        /*!
          \param x phase-space point (x,v)
          \param params -- does nothing

          \return angles and frequencies --
                6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
        */
		VecDoub angles(const VecDoub& x,void*params=nullptr);
		/*! get max intermediate axis (y) intercept */
		inline double get_ymax(void){return ymax;}
		/*! get min intermediate axis (y) intercept */
		inline double get_ymin(void){return ymin;}

		// Other development functions
		std::vector<int> angular_momentum(const VecDoub &x);
		int check_ang_mom(double y, double E, int swit);
		double find_closed(double E, int swit);
		VecDoub check_orbit(double y, double E, int swit, int plot=0);
		void plot_Delta2(double E);
		void plot_Delta1(double E);
		double find_beta(double E,double a, double b);
		double find_alpha(double E,double a, double b);
		VecDoub find_best_alphabeta(const VecDoub&);
		void alphabeta_grid(const VecDoub& X);
		VecDoub find_ab_from_box(double E);
		VecDoub actionSD(const std::vector<VecDoub>&,bool);
		void fillDeltagrids(const std::string& f);
		void readDeltagrids(const std::string& f);
		VecDoub findDelta_interp(double E);
		double sos(int comp,const VecDoub& x,const std::string& outfile);
};

/*! Various structures for passing to minimisation routines etc. */
struct root_struct_mindist{
	Potential_JS *Pot;
	double E;
	int swit;
	root_struct_mindist(Potential_JS *Pot, double E, int swit)
		:Pot(Pot),E(E),swit(swit){}
};

struct root_struct_action{
	Actions_TriaxialStackel_Fudge *ATSF;
	Orbit *orbit;
	int swit;
	root_struct_action(Actions_TriaxialStackel_Fudge *ATSF, Orbit *orbit, int swit)
		:ATSF(ATSF),orbit(orbit),swit(swit){}
};

struct root_struct_actions{
	Potential_JS *Pot;
	Orbit *orbit;
	root_struct_actions(Potential_JS *Pot, Orbit *orbit)
		:Pot(Pot),orbit(orbit){}
};

#endif
// ============================================================================
