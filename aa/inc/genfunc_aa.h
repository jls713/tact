/*=============================================*/
/*        Actions from generating function     */
/*=============================================*/
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

struct Actions_Genfunc_data_structure{
    // Stores parameters for Actions_Genfunc and Actions_Genfunc_Average
    // such that the signature of VecDoub actions is uniform across all
    // action approximation routines
	double total_T, stepsize; unsigned int N_T;
    int N_matrix;
    double orbit_eps, NTMAX, maxtimescale;
    bool axisymmetric;
	Actions_Genfunc_data_structure(double tt,unsigned int nt, unsigned int N_mat,double oe, int NTMAX, double maxtimescale, bool axisymmetric=false):total_T(tt), N_T(nt),N_matrix(N_mat),orbit_eps(oe), NTMAX(NTMAX),maxtimescale(maxtimescale),axisymmetric(axisymmetric){
		stepsize=total_T/(double)N_T;
        if(N_T<(axisymmetric==true ? 1.:N_mat)*N_mat*N_mat/2.)
            std::cerr<<"Not enough time samples ("<<N_mat
                     <<") to constrain N_max="<<N_mat<<"modes"<<std::endl;
	};
};

//============================================================================

class Actions_Genfunc : public Action_Finder{
    // Finds actions by calculating the generating function components from
    // a orbit integration via toy actions in either the isochrone or harmonic
    // oscillator potential
    // The actions routine accepts a set of parameters = (Total time, N_T, N_max)
    // If no parameters are passed the routine will estimate the orbital time
    // scale from T = Potential_JS::torb and use Total_T = 8*T, N_T=200 and
    // N_max = 6
    // Currently geared up to return J_lambda, J_mu, J_nu -- this means that
    // for the inner long axis loops (identified via the Stackel fudge) have
    // J_r <==> L such that J_mu is the radial action
    protected:
        std::string symmetry;
        // Doubles radial action if triaxial symmetry for defining continuous action surface
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
        Actions_Genfunc(Potential_JS *Pot, std::string symmetry = "triaxial",lmn_orb *AF=nullptr,bool la_switch=false)
            :symmetry(symmetry),AF(AF),TargetPot(Pot),la_switch(la_switch){
                if(symmetry!="triaxial" and symmetry!="axisymmetric"){
                    std::cerr<<"Symmetry must be triaxial or axisymmetric\n";
                }
                if(symmetry=="spherical")
                    std::cerr<<"Use Actions_Spherical for spherical symmetry!\n";
            };
        //Actions_Genfunc(Potential_JS *Pot):scaling(true),AF(nullptr),TargetPot(Pot){};
        VecDoub actions(const VecDoub& x, void *params=nullptr);
        inline void reset(Potential_JS *pot){TargetPot = pot;
            if(AF) AF->reset(pot);}
        VecDoub angles(const VecDoub& x, void *params=nullptr);
        VecDoub full_actions(const VecDoub& x, int,int,int);
        VecDoub full_angles(const VecDoub& x, int,int,int);
};

//============================================================================

class Actions_Genfunc_Average : public Actions_Genfunc{
    // Finds actions by averaging using an orbit integration via toy actions
    // in either the isochrone or harmonic oscillator potential
    private:
    public:
        Actions_Genfunc_Average(Potential_JS *Pot,std::string symmetry,lmn_orb *AF=nullptr,bool la_switch=false):Actions_Genfunc(Pot,symmetry,AF,la_switch){};
        VecDoub actions(const VecDoub& x, void *params=nullptr);
};

#endif
//============================================================================
