#ifndef LMN_H
#define LMN_H
#include "aa.h"
// ============================================================================
// Lambda_Mu_Nu orbit action finding
// ============================================================================

class lmn_orb: public Action_Finder{
	// Interface for using triaxial St\"ackel fudge apparatus
	// Handles the choice of coordinate system
	private:
		VecDoub E;
		VecDoub Beta;
		VecDoub Alpha;
		Potential_JS *Pot;
		// Actions_TriaxialStackel_Fudge *ATSF;
		double ymin; 		// Maximum y intercept
		double ymax; 		// Maximum y intercept
		int NE;				// Size of grid
		double Emax = 0.; 	// Maximum energy
		double E0 = 0.;		// Core energy
		std::string name;
		bool use_log_grid, use_acts;
	public:

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
		inline double get_ymax(void){return ymax;}
		inline double get_ymin(void){return ymin;}
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
		VecDoub actions(const VecDoub& x,void*params=nullptr);
		VecDoub angles(const VecDoub& x,void*params=nullptr);
		double sos(int comp,const VecDoub& x,const std::string& outfile);
};

// Various structures for passing to minimisation routines etc.

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
