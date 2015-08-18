#ifndef LMNELLZ_H
#define LMNELLZ_H
#include "aa.h"
// ============================================================================
// Lambda_Mu_Nu orbit action finding
// ============================================================================

class lmn_orb_ELLz: public Action_Finder{
	// Interface for using triaxial St\"ackel fudge apparatus
	// Handles the choice of coordinate system
	private:
		VecDoub E;
		VecDoub Lmax;
		VecDoub IL,ILz;
		MatMatDoub Beta;
		MatMatDoub Alpha;
		Potential_JS *Pot;
		// Actions_TriaxialStackel_Fudge *ATSF;
		double ymin; 		// Maximum y intercept
		double ymax; 		// Maximum y intercept
		int NE,NL,NLz;		// Size of grid
		double Emax = 0.; 	// Maximum energy
		double E0 = 0.;		// Core energy
		double minIL,maxIL;
		std::string name;
		bool use_log_grid, use_acts;
	public:

		lmn_orb_ELLz(Potential_JS *pot, double ym=0.05, double yn=80.,int NEE=32,int NLL=32,int NLz=32,bool use_log_grid=false,double minIL=0.02,double maxIL=0.98,std::string s="")
			:Pot(pot),ymin(ym),ymax(yn),NE(NEE),NL(NLL),NLz(NLz)
			,minIL(minIL),maxIL(maxIL)
			,name(s),use_log_grid(use_log_grid){
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
		inline double get_ymax(void){return ymax;}
		inline double get_ymin(void){return ymin;}
		std::vector<int> angular_momentum(const VecDoub &x);
		int check_ang_mom(double y, double E, int swit);
		double find_closed(double E, int swit);
		VecDoub check_orbit(double y, double E, int swit, int plot=0);
		VecDoub find_best_alphabeta(const VecDoub&);
		void alphabeta_grid(const VecDoub& X);
		VecDoub find_ab_from_box(double E);
		VecDoub actionSD(const std::vector<VecDoub>&,bool);
		void fillDeltagrids(const std::string& f);
		void readDeltagrids(const std::string& f);
		VecDoub findDelta_interp(double E, double L, double Lz);
		VecDoub actions(const VecDoub& x,void*params=nullptr);
		VecDoub angles(const VecDoub& x,void*params=nullptr);
};

// Various structures for passing to minimisation routines etc.

struct root_struct_mindistELLZ{
	Potential_JS *Pot;
	double E;
	int swit;
	root_struct_mindistELLZ(Potential_JS *Pot, double E, int swit)
		:Pot(Pot),E(E),swit(swit){}
};

struct root_struct_actionELLZ{
	Actions_TriaxialStackel_Fudge *ATSF;
	Orbit *orbit;
	int swit;
	root_struct_actionELLZ(Actions_TriaxialStackel_Fudge *ATSF, Orbit *orbit, int swit)
		:ATSF(ATSF),orbit(orbit),swit(swit){}
};

struct root_struct_actionsELLZ{
	Potential_JS *Pot;
	Orbit *orbit;
	double rmax;
	root_struct_actionsELLZ(Potential_JS *Pot, Orbit *orbit, double rmax)
		:Pot(Pot),orbit(orbit),rmax(rmax){}
};

#endif
// ============================================================================
