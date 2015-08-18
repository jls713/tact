#ifndef uv_H
#define uv_H
#include "aa.h"
// ============================================================================
// uv orbit action finding
// ============================================================================

class uv_orb: public Action_Finder{
	// Interface for using axisymmetric St\"ackel fudge apparatus
	// Handles the choice of coordinate system
	private:
		VecDoub E_delta;
		std::vector<VecDoub> L_delta, Delta;
		Potential_JS *Pot;
		double Rmin; 		// Maximum R intercept
		double Rmax; 		// Maximum R intercept
		double Emax = 0.; 	// Maximum energy
		double E0 = 0.;		// Core energy
		int NE, NL;				// Size of grid
		std::string name;
	public:

		uv_orb(Potential_JS *pot, double Rm=0.05, double Rn=80.,int NE=20, int NL = 10,std::string name="")
			:Pot(pot),Rmin(Rm),Rmax(Rn),NE(NE),NL(NL),name(name){
			E0 = Pot->E_circ(Rmin);
			Emax = Pot->E_circ(Rmax);
			// Emax+=Pot->Phi({Rn,0.,Rn})-Emax;
			fillDeltagrids(name!=""?name+".Delta_uv":"Delta_uv.tmp");
		}
		void reset(Potential_JS *pot){
			Pot = pot;
			Emax = Pot->E_circ(Rmax);
			E0 = Pot->E_circ(Rmin);
			fillDeltagrids(name!=""?name+".Delta_uv":"Delta_uv.tmp");
		}
		void fillDeltagrids(const std::string& f);
		void readDeltagrids(const std::string& f);
		double findDelta_interp(double E,double L);
		VecDoub actions(const VecDoub& x,void*params=nullptr);
		VecDoub angles(const VecDoub& x,void*params=nullptr);
};

#endif
// ============================================================================
