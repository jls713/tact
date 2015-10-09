// ============================================================================
/// \file inc/genfunc_aa.h
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
/// \brief Action finding routines using Generating function
///
/// Computes actions via orbit integration. There are two classes given:
/// 1. Actions_Genfunc computes Fourier coefficients of the generating function
///    from a toy action system (isochrone or harmonic oscillator). This
///    procedure is detailed in Sanders \& Binney (2014).
/// 2. Actions_Genfunc_Average computes the actions by averaging the toy
///    actions over the toy angles (as in Fox (2014) and Bovy (2014)).
///
// ============================================================================
// Notes:
//      There are currently problems if z=0, vz=0
//============================================================================

#ifndef GENFUNC_AA_H
#define GENFUNC_AA_H

//============================================================================

#include "potential.h"
#include "aa.h"
#include "stackel_aa.h"
#include "orbit.h"
#include "lmn_orb.h"

//============================================================================
/*! Stores parameters for Actions_Genfunc and Actions_Genfunc_Average
    such that the signature of VecDoub actions is uniform across all
    action approximation routines */
struct Actions_Genfunc_data_structure{
	double total_T;                        /*! Total integration time */
    double stepsize;                       /*! Integration step size  */
    unsigned int N_T;                      /*! Total time samples     */
    int N_matrix;                          /*! No. of S_n coefficients*/
    double orbit_eps;             /*! Tolerance for orbit integration */
    double NTMAX;                      /*! Max number of time samples */
    double maxtimescale;                     /*! Max integration time */
    bool axisymmetric;           /*! axisymmetric (1) or triaxial (0) */
    /*! Actions_Genfunc_data_structure constructor.
        \param tt total integration time
        \param nt total time samples
        \param N_mat No. of S_n coefficients
        \param oe tolerance for orbit integration
        \param NTMAX max number of time samples
        \param maxtimescale max integration time
        \param axisymmetric -- axisymmetric (1) or triaxial (0)
    */
	Actions_Genfunc_data_structure(double tt,unsigned int nt, unsigned int N_mat,double oe, int NTMAX, double maxtimescale, bool axisymmetric=false):total_T(tt), N_T(nt),N_matrix(N_mat),orbit_eps(oe), NTMAX(NTMAX),maxtimescale(maxtimescale),axisymmetric(axisymmetric){
		stepsize=total_T/(double)N_T;
        if(N_T<(axisymmetric==true ? 1.:N_mat)*N_mat*N_mat/2.)
            std::cerr<<"Not enough time samples ("<<N_mat
                     <<") to constrain N_max="<<N_mat<<"modes"<<std::endl;
	};
};

//============================================================================
/*! Action finder using Generating function
    Finds actions by calculating the generating function components from
    a orbit integration via toy actions in either the isochrone or harmonic
    oscillator potential
    The actions routine accepts a set of parameters = (Total time, N_T, N_max)
    If no parameters are passed the routine will estimate the orbital time
    scale from T = Potential_JS::torb and use Total_T = 8*T, N_T=200 and
    N_max = 6
    Currently geared up to return J_lambda, J_mu, J_nu -- this means that
    for the inner long axis loops (identified via the Stackel fudge) have
    J_r <==> L such that J_mu is the radial action
*/
class Actions_Genfunc : public Action_Finder{
    protected:
        std::string symmetry;/**< Encodes the degree of symmetry: "axisymmetric" or "triaxial". Doubles radial action if triaxial symmetry for defining continuous action surface */
        lmn_orb *AF;
        Action_Finder *ToyAct;
        Potential_JS *TargetPot;
        bool la_switch;

        std::vector<int> angular_momentum(const VecDoub &x);
		std::vector<int> loop(const std::vector<VecDoub> &orbit_samples);
        VecDoub find_isochrone_params(const std::vector<VecDoub> &orbit_samples);
        VecDoub find_box_params(const std::vector<VecDoub> &orbit_samples);
        VecDoub find_box_params_minvar(const std::vector<VecDoub> &orbit_samples);
        VecDoub find_isochrone_params_minvar(const std::vector<VecDoub> &orbit_samples);
    public:
        /*! Actions_Genfunc constructor.
            \param Pot potential (axisymmetric or triaxial)
            \param symmetry ('axisymmetric' or 'triaxial')
            \param AF pointer to lmn_orb action finder for determining orbit class (not necessary)
            \param la_switch -- boolean for switching actions for long-axis loops
        */
        Actions_Genfunc(Potential_JS *Pot, std::string symmetry = "triaxial",lmn_orb *AF=nullptr,bool la_switch=false)
            :symmetry(symmetry),AF(AF),TargetPot(Pot),la_switch(la_switch){
                if(symmetry!="triaxial" and symmetry!="axisymmetric"){
                    std::cerr<<"Symmetry must be triaxial or axisymmetric\n";
                }
                if(symmetry=="spherical")
                    std::cerr<<"Use Actions_Spherical for spherical symmetry!\n";
            };
        //! reset potential
        inline void reset(Potential_JS *pot){TargetPot = pot;
            if(AF) AF->reset(pot);}
        //! Finds actions
        /*!
          \param x phase-space point (x,v)
          \param params -- can pass pointer to Actions_Genfunc_data_structure that contains parameters for orbit integration

          \return actions -- 3D vector J=(J_R,J_phi,J_z)
        */
        VecDoub actions(const VecDoub& x, void *params=nullptr);
        //! Finds angles
        /*!
          \param x phase-space point (x,v)
          \param params -- can pass pointer to Actions_Genfunc_data_structure that contains parameters for orbit integration

          \return angles and frequencies --
                6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
        */
        VecDoub angles(const VecDoub& x, void *params=nullptr);
        //! Finds actions - same as above but pass individual params not struct
        /*!
          \param x phase-space point (x,v)
          \param NT number of orbital times to integrate for
          \param Nsamp number of time samples
          \param Nmax number of Fourier coefficients Sn

          \return actions -- 3D vector J=(J_R,J_phi,J_z)
        */
        VecDoub full_actions(const VecDoub& x, int NT,int Nsamp,int Nmax);
        //! Finds angles
        /*!
          \param x phase-space point (x,v)
          \param NT number of orbital times to integrate for
          \param Nsamp number of time samples
          \param Nmax number of Fourier coefficients Sn

          \return angles and frequencies --
                6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
        */
        VecDoub full_angles(const VecDoub& x, int NT,int Nsamp,int Nmax);
};

//============================================================================
/*! Action finder averaging the toy actions over toy angles
    Finds actions by averaging using an orbit integration via toy actions
    in either the isochrone or harmonic oscillator potential
*/
class Actions_Genfunc_Average : public Actions_Genfunc{
    private:
    public:
        Actions_Genfunc_Average(Potential_JS *Pot,std::string symmetry,lmn_orb *AF=nullptr,bool la_switch=false):Actions_Genfunc(Pot,symmetry,AF,la_switch){};
        //! Finds actions
        /*!
          \param x phase-space point (x,v)
          \param params -- can pass pointer to Actions_Genfunc_data_structure that contains parameters for orbit integration

          \return actions -- 3D vector J=(J_R,J_phi,J_z)
        */
        VecDoub actions(const VecDoub& x, void *params=nullptr);
};

#endif
//============================================================================
