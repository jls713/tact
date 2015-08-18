#ifndef ADIABATIC_H
#define ADIABATIC_H

#include "potential.h"
#include "aa.h"
// ============================================================================
// Polar Adiabatic Approximation
// ============================================================================

class Actions_PolarAdiabaticApproximation : public Action_Finder{
    private:
        Potential_JS *Pot;
        bool no_energy_correction;

        const int NGRID = 60;
        const double Rmin = 1., Rmax = 40.;
        VecDoub Rgrid, Ezmaxgrid;
        std::vector<VecDoub> Ezgrid, Jzgrid;
    public:
        Actions_PolarAdiabaticApproximation(Potential_JS *pot,std::string filename="",bool write=false,bool no_energy_corr=false)
            : Pot(pot), no_energy_correction(no_energy_corr){
            if(filename!="" and !write)
                load_grids(filename);
            else
                make_grids(write?filename:"");
        };
        double actions_Jz(double R, double Ez, double z, double *zlim);
        VecDoub actions(const VecDoub& x, void*params=nullptr);
        inline double Phi_z(VecDoub Polar){
            return Pot->Phi({Polar[0],0.,Polar[2]})-Pot->Phi({Polar[0],0.,0.});
        }
        inline double PhiR_eff(double R, double Lz2, double Jz){
            return Pot->Phi({R,0.,0.})+.5*Lz2/(R*R)+Ez_from_grid(R,Jz);
        }
        double find_zlimit(double Ez, VecDoub Polar);
        VecDoub find_Rlimits(double R, double Etot, double Lz2,double Jz);
        void load_grids(std::string filename);
        void make_grids(std::string filename);
        double Ez_from_grid(double R, double Jz);
};

struct PolarAA_zactions_struct{
    Actions_PolarAdiabaticApproximation *AA;
    double Ez, R, zlim;
    PolarAA_zactions_struct(Actions_PolarAdiabaticApproximation *AA,double Ez, double R,double zlim)
        :AA(AA),Ez(Ez), R(R), zlim(zlim){}
};

struct PolarAA_Ractions_struct{
    Actions_PolarAdiabaticApproximation *AA;
    double Etot, Lz2, Jz, Delta,taubar;
    PolarAA_Ractions_struct(Actions_PolarAdiabaticApproximation *AA,double Etot, double Lz2, double Jz, double Delta,double taubar)
        :AA(AA),Etot(Etot), Lz2(Lz2), Jz(Jz), Delta(Delta), taubar(taubar){}
};

// ============================================================================
// Spheroidal Adiabatic Approximation
// ============================================================================

class Actions_SpheroidalAdiabaticApproximation : public Action_Finder{
    private:
        bool alpha_guess;
    public:
        Potential_JS *Pot;
        std::unique_ptr<OblateSpheroidCoordSys> CS;
        bool no_energy_correction;

        const int NGRID = 60, NL = 5;
        const double Rmin = 1., Rmax = 40.;
        double Lmin, Lmax;
        VecDoub Rgrid, Lgrid;
        std::vector<VecDoub> Ezmaxgrid;
        std::vector<std::vector<VecDoub>> Ezgrid, Jzgrid;
    public:
        Actions_SpheroidalAdiabaticApproximation(Potential_JS *pot,std::string filename="",bool write=false,bool no_energy_corr=false,double alpha = 100.)
            : Pot(pot), no_energy_correction(no_energy_corr){

            Lmin=Pot->L_circ(Rmin);Lmax=Pot->L_circ(Rmax);

            if(alpha>0.) alpha_guess=false;
            else {alpha_guess = true; alpha=-2;}
            CS =std::unique_ptr<OblateSpheroidCoordSys>(new OblateSpheroidCoordSys(alpha));

            if(filename!="" and !write)
                load_grids(filename);
            else
                make_grids(write?filename:"");
        };
	Actions_SpheroidalAdiabaticApproximation(const Actions_SpheroidalAdiabaticApproximation& s):
		Pot(s.Pot),CS(new OblateSpheroidCoordSys(*s.CS)),no_energy_correction(s.no_energy_correction),Rgrid(s.Rgrid),Lgrid(s.Lgrid),Ezmaxgrid(s.Ezmaxgrid),Ezgrid(s.Ezgrid),Jzgrid(s.Jzgrid){
            Lmin=Pot->L_circ(Rmin);Lmax=Pot->L_circ(Rmax);
        }
        double actions_Jz(double ENu, double Lz2, VecDoub tau, double *nulim);
        VecDoub actions(const VecDoub& x, void*params=nullptr);
        inline double Phi_nu(VecDoub tau, double Lz2){
            VecDoub x = CS->tau2x(tau);
            return Pot->Phi(x)-Pot->Phi({sqrt(tau[0]+CS->alpha()),0.,0.})
                    +.5*Lz2*(1./normsq<double>({x[0],x[1]})-1./(tau[0]+CS->alpha()));
        }
        inline double Philam_eff(VecDoub tau, double Lz2, double Jz){
            double R = sqrt(tau[0]+CS->alpha());
            return Pot->Phi({R,0.,0.})+.5*Lz2/(R*R)+Ez_from_grid(Lz2,R,Jz);
        }
        double find_nulimit(double ENu, double Lz2, VecDoub tau);
        VecDoub find_lamlimits(double Etot, double Lz2,double Jz,VecDoub tau);
        void load_grids(std::string filename);
        void make_grids(std::string filename);
        double Ez_from_grid(double Lz2, double R, double Jz);
};

struct SpheroidalAA_zactions_struct{
    Actions_SpheroidalAdiabaticApproximation *AA;
    double ENu, Lz2, lam, Delta, taubar;
    SpheroidalAA_zactions_struct(Actions_SpheroidalAdiabaticApproximation *AA,double ENu, double Lz2, double lam, double Delta, double taubar)
        :AA(AA),ENu(ENu), Lz2(Lz2), lam(lam), Delta(Delta), taubar(taubar){}
};

struct SpheroidalAA_Ractions_struct{
    Actions_SpheroidalAdiabaticApproximation *AA;
    double Elam, Lz2, Jz, nu, Delta,taubar;
    SpheroidalAA_Ractions_struct(Actions_SpheroidalAdiabaticApproximation *AA,double Elam, double Lz2, double Jz, double nu, double Delta,double taubar)
        :AA(AA),Elam(Elam), Lz2(Lz2), Jz(Jz), nu(nu), Delta(Delta), taubar(taubar){}
};

#endif
