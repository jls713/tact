// ============================================================================
/// \file inc/adiabatic_aa.h
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
/// \brief Action finders for adiabatic approximations
///
/// There are two adiabatic approximations -- polar and spheroidal
/// The two classes are very similar in that the actions are estimated by
/// assuming the potential is separable in either a polar or spheroidal
/// coordinate system.
///
/// The polar method is taken from Schoenrich \& Binney (2012).
/// The spheroidal method is from Sanders \& Binney (2015).
///
/// Both classes first construct grids of the vertical energy as a function
/// of vertical action, angular momentum and radius. The properties of the
/// grids are governed by the parameters Rmin (minimum radius), Rmax (maximum
/// radius), ZMAX (maximum z height), NGRID (number of gridpoints in R and Jz),
/// NL (number of grid points in Lz)
// ============================================================================

#ifndef ADIABATIC_H
#define ADIABATIC_H

#include "potential.h"
#include "aa.h"
// ============================================================================
/*! Polar Adiabatic Approximation */
// ============================================================================
class Actions_CylindricalAdiabaticApproximation : public Action_Finder{
    private:
        Potential_JS *Pot;          /*! Potential (axisymmetric)             */
        bool no_energy_correction;  /*! if true, no radial-vertical coupling */
        double Rmin, Rmax, ZMAX;/*!min and max radius, max z height for grids*/
        int NGRID;            /*! number of grid points for Ez         */
        VecDoub Rgrid, Ezmaxgrid;   /*! radius and energy grids              */
        std::vector<VecDoub> Ezgrid;/*! vertical energy grid                 */
        std::vector<VecDoub> Jzgrid;/*! vertical action grid                 */
        std::vector<VecDoub> dJzgrid;/*! vertical freq grid                 */

        double actions_Jz(double R, double Ez, double z, double *zlim);
        double find_zlimit(double Ez, VecDoub Polar);
        VecDoub find_Rlimits(double R, double Etot, double Lz2,double Jz);
        void load_grids(std::string filename);
        void make_grids(std::string filename);
        double Ez_from_grid(double R, double Jz);
        double Jz_from_grid(double R, double Jz);
        double estimate_tiny(double Ez, double R);
        /*! computes derivative of Jz wrt R at fixed Ez */
        double dJzdR(double R, double Jz);
        double actions_dJzdEz(double R, double Ez,double z, double *zlim);
    public:
        //! Actions_CylindricalAdiabaticApproximation constructor.
        /*!
          \param pot Potential (axisymmetric) in which to compute the actions
          \param filename -- output file for Ez grid
          \param write -- output file if true
          \param no_energy_corr -- if true, no coupling between radial and vertical motion
          \param Rm -- minimum radius for grid
          \param Rn -- maximum radius for grid
          \param zmax -- maximum z height for grid
          \param NGRID -- number of grid points for Ez
        */
        Actions_CylindricalAdiabaticApproximation(Potential_JS *pot,std::string filename="",bool write=false,bool no_energy_corr=false, double Rm=1., double Rn=40.,double zmax=40.,int NGRID=60)
            : Pot(pot), no_energy_correction(no_energy_corr), Rmin(Rm), Rmax(Rn), ZMAX(zmax), NGRID(NGRID){
            if(filename!="" and !write)
                load_grids(filename);
            else
                make_grids(write?filename:"");
        };
        /*! Effective vertical potential */
        inline double Phi_z(VecDoub Polar){
            return Pot->Phi({Polar[0],0.,Polar[2]})-Pot->Phi({Polar[0],0.,0.});
        }
        /*! Effective radial potential */
        inline double PhiR_eff(double R, double Lz2, double Jz){
            return Pot->Phi({R,0.,0.})+.5*Lz2/(R*R)+Ez_from_grid(R,Jz);
        }
        /*! computes derivative of Jz wrt Ez at fixed radius */
        double dJzdEz(double R, double Jz);
        //! Finds actions
        /*!
          \param x phase-space point (x,v)
          \param params does nothing
          \return actions -- 3D vector J=(J_R,J_phi,J_z)
        */
        VecDoub actions(const VecDoub& x, void*params=nullptr);
        //! Finds angles
        /*!
          \param x phase-space point (x,v)
          \param params -- does nothing

          \return angles and frequencies --
                6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
        */
        VecDoub angles(const VecDoub& x, void *params=nullptr);
};
/*! Helper structure for Jz polar adiabatic approx integration */
struct PolarAA_zactions_struct{
    Actions_CylindricalAdiabaticApproximation *AA;
    double Ez, R, zlim, tiny_number;
    PolarAA_zactions_struct(Actions_CylindricalAdiabaticApproximation *AA,double Ez, double R,double zlim, double tn)
        :AA(AA),Ez(Ez), R(R), zlim(zlim), tiny_number(tn){}
};

/*! Helper structure for JR polar adiabatic approx integration */
struct PolarAA_Ractions_struct{
    Actions_CylindricalAdiabaticApproximation *AA;
    double Etot, Lz2, Jz, Delta,taubar, tiny_number;
    PolarAA_Ractions_struct(Actions_CylindricalAdiabaticApproximation *AA,double Etot, double Lz2, double Jz, double Delta,double taubar, double tn)
        :AA(AA),Etot(Etot), Lz2(Lz2), Jz(Jz), Delta(Delta), taubar(taubar), tiny_number(tn){}
};

// ============================================================================
/*! Spheroidal Adiabatic Approximation */
// ============================================================================
class Actions_SpheroidalAdiabaticApproximation : public Action_Finder{
    private:
        bool alpha_guess;
        Potential_JS *Pot;          /*! Potential (axisymmetric)             */
        bool no_energy_correction;  /*! if true, no radial-vertical coupling */

        double Rmin, Rmax, ZMAX;/*!min and max radius, max z height for grids*/
        int NGRID;            /*! number of grid points for Enu        */
        int NL;               /*! number of grid points for Lz         */
        double Lmin, Lmax;          /*! min and max ang mom for grids        */
        VecDoub Rgrid, Lgrid;       /*! radius and ang mom grids             */
        std::vector<VecDoub> Ezmaxgrid;/*! max vertical energy grid          */
        std::vector<std::vector<VecDoub>> Ezgrid;/*! vertical energy grid    */
        std::vector<std::vector<VecDoub>> Jzgrid;/*! vertical action grid    */

        double find_nulimit(double ENu, double Lz2, VecDoub tau);
        VecDoub find_lamlimits(double Etot, double Lz2,double Jz,VecDoub tau);
        void load_grids(std::string filename);
        void make_grids(std::string filename);
        double Ez_from_grid(double Lz2, double R, double Jz);
        double Jz_from_grid(double Lz2, double R, double Ez);
        double dJzdR(double Lz2,double R,double Ez);
        double estimate_tiny(double Ez, double l, double n, double Lz2);
        double actions_Jz(double ENu, double Lz2, VecDoub tau, double *nulim);
    public:
        //! Actions_SpheroidalAdiabaticApproximation constructor.
        /*!
          \param pot Potential (axisymmetric) in which to compute the actions
          \param filename -- output file for Ez grid
          \param write -- output file if true
          \param no_energy_corr -- if true, no coupling between radial and vertical motion
          \param alpha - alpha estimate (can be overridden in actions & angles)
          \param Rm -- minimum radius for grid
          \param Rn -- maximum radius for grid
          \param zmax -- maximum z height for grid
          \param NGRID -- number of grid points for Enu
          \param NL -- number of grid points for Lz
        */
        Actions_SpheroidalAdiabaticApproximation(Potential_JS *pot,std::string filename="",bool write=false,bool no_energy_corr=false,double alpha = 100., double Rm=1., double Rn=40.,double zmax=40., int NGRID=60, int NL=10)
            : Pot(pot), no_energy_correction(no_energy_corr), Rmin(Rm), Rmax(Rn), ZMAX(zmax), NGRID(NGRID), NL(NL){

            Lmin=Pot->L_circ(Rmin);Lmax=Pot->L_circ(Rmax);

            if(alpha>0.) alpha_guess=false;
            else {alpha_guess = true; alpha=-2;}
            CS =std::unique_ptr<ProlateSpheroidCoordSys>(new ProlateSpheroidCoordSys(alpha));

            if(filename!="" and !write)
                load_grids(filename);
            else
                make_grids(write?filename:"");
        };
        //! Actions_SpheroidalAdiabaticApproximation copy constructor.
    	Actions_SpheroidalAdiabaticApproximation(const Actions_SpheroidalAdiabaticApproximation& s):
		Pot(s.Pot),no_energy_correction(s.no_energy_correction), Rmin(s.Rmin), Rmax(s.Rmax), ZMAX(s.ZMAX),Rgrid(s.Rgrid),Lgrid(s.Lgrid),Ezmaxgrid(s.Ezmaxgrid),Ezgrid(s.Ezgrid),Jzgrid(s.Jzgrid),CS(new ProlateSpheroidCoordSys(*s.CS)){
            Lmin=Pot->L_circ(Rmin);Lmax=Pot->L_circ(Rmax);
        }
        std::unique_ptr<ProlateSpheroidCoordSys> CS;
        /*! Effective vertical potential */
        inline double Phi_nu(VecDoub tau, double Lz2){
            VecDoub x = CS->tau2x(tau);
            return Pot->Phi(x)-Pot->Phi({sqrt(tau[0]+CS->alpha()),0.,0.})
                    +.5*Lz2*(1./normsq<double>({x[0],x[1]})-1./(tau[0]+CS->alpha()));
        }
        /*! Effective radial potential */
        inline double Philam_eff(VecDoub tau, double Lz2, double Jz){
            double R = sqrt(tau[0]+CS->alpha());
            return Pot->Phi({R,0.,0.})+.5*Lz2/(R*R)+Ez_from_grid(Lz2,R,Jz);
        }
        //! Finds actions
        /*!
          \param x phase-space point (x,v)
          \param params -- if null, use alpha specified in constructor
                        -- if >0, estimate using derivatives of potential eq(8) Sanders (2012)
                        -- if <0, use value passed as new alpha

          \return actions -- 3D vector J=(J_R,J_phi,J_z)
        */
        VecDoub actions(const VecDoub& x, void*params=nullptr);
        //! Finds angles
        /*!
          \param x phase-space point (x,v)
          \param params -- if null, use alpha specified in constructor
                        -- if >0, estimate using derivatives of potential eq(8) Sanders (2012)
                        -- if <0, use value passed as new alpha

          \return angles and frequencies --
                6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
        */
        VecDoub angles(const VecDoub& x, void *params=nullptr);
        //! Derivative of Jz wrt Ez at fixed Lz2 and R
        double dJzdEz(double Lz2,double R,double Jz);
};

/*! Helper structure for Jz spheroidal adiabatic approx integration */
struct SpheroidalAA_zactions_struct{
    Actions_SpheroidalAdiabaticApproximation *AA;
    double ENu, Lz2, lam, Delta, taubar, tiny_number;
    SpheroidalAA_zactions_struct(Actions_SpheroidalAdiabaticApproximation *AA,double ENu, double Lz2, double lam, double Delta, double taubar, double tn)
        :AA(AA),ENu(ENu), Lz2(Lz2), lam(lam), Delta(Delta), taubar(taubar), tiny_number(tn){}
};

/*! Helper structure for JR spheroidal adiabatic approx integration */
struct SpheroidalAA_Ractions_struct{
    Actions_SpheroidalAdiabaticApproximation *AA;
    double Elam, Lz2, Jz, nu, Delta,taubar, tiny_number;
    SpheroidalAA_Ractions_struct(Actions_SpheroidalAdiabaticApproximation *AA,double Elam, double Lz2, double Jz, double nu, double Delta,double taubar, double tn)
        :AA(AA),Elam(Elam), Lz2(Lz2), Jz(Jz), nu(nu), Delta(Delta), taubar(taubar), tiny_number(tn){}
};

#endif
// ============================================================================
