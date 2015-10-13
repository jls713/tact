// ============================================================================
/// \file inc/stackel_aa.h
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
/// \brief Action estimation by fitting Staeckel potentials
///
/// 1. Stackel_Fitted_Potential
/// Implements a class that finds the best-fitting Staeckel potential for some
/// general axisymmetric potential using the method in Dejonghe & de Zeeuw
/// (1988).
/// 2. Actions_StackelFit
/// Fits Staeckel potential to region a given orbit probes and finds actions
/// in this best-fitting potential
//============================================================================

#ifndef STACK_FIT_H
#define STACK_FIT_H

#include "utils.h"
#include "potential.h"
#include "aa.h"

// ============================================================================
/*! Potential class that finds the best-fitting Staeckel potential to a
    axisymmetric potential */
class Stackel_Fitted_Potential: public StackelOblate_PerfectEllipsoid{
    private:
        const int DATAPOINTS = 40; /*< Number of gridpoints             */
        const int STEPMAX = 300;   /*< Max number of integrations       */
        const int NMAX = 10;       /*< Number of bisections to perform  */
        const double cl = 4.;      /*< Power of lambda weight function  */
        const double cv = 0.5;     /*< Power of nu weight function      */
        VecDoub limits;            /*< Limits of fitting regions in tau */
        Potential_JS *TruePot;     /*< Target potential                 */
        std::unique_ptr<Orbit> Orb;/*< Orbit for finding limits         */
        double *lambdagrid,*nugrid;/*< Stores tau grid for best-fit     */
        double *flamgrid,*fnugrid; /*< Stores f(tau) for best-fit       */
        double *y2Lamgrid,*y2Nugrid;/*<Stores second derivs for best-fit*/
        double L, N, chibar;

        double Vlamnu(double lambda, double nu);
        double flamINT(double lam);
        double fnuINT(double nu);
        double flamINT_DERIV(double lam);
        double fnuINT_DERIV(double nu);
        double G(double tau);
        double BigF(double tau);
        double BigFPrime(double tau);
        double GPrime(double tau);
        double fnu(double);
        double flam(double);
        double find_turning(VecDoub x, VecDoub y, double step, int lamnu);
        void find_limits(VecDoub x);
    public:
        //! Stackel_Fitted_Potential constructor.
        /*!
          \param TruePot Potential_JS (axisymmetric)

        */
        Stackel_Fitted_Potential(Potential_JS *TruePot)
            :StackelOblate_PerfectEllipsoid(10.,-30.),TruePot(TruePot){
                Orb = std::unique_ptr<Orbit>(new Orbit(TruePot,1e-8));
            };
        //! Stackel_Fitted_Potential destructor.
        ~Stackel_Fitted_Potential(){
            if(lambdagrid){
                delete lambdagrid;delete nugrid;delete flamgrid;
                delete fnugrid; delete y2Lamgrid; delete y2Nugrid;
            }
        }
        //! Fits Staeckel potential to target potential
        /*!
          \param x phase-space coordinate (x,v)
          This is the main routine!
        */
        void fit_potential(VecDoub x);
        //! chi = (lambda-nu)Phi(lambda, nu)
        double chi(double lambda, double nu);
        //! BigLambda -- nu weight function
        double BigNu(double nu);
        //! BigLambda -- lambda weight function
        double BigLambda(double lambda);
        //! chibarint -- chi*BigLambda*BigNu
        double chibarint(double lambda,double nu);
};

// ============================================================================
/*!  Action finding by first fitting a Staeckel potential to the region an
     orbit probes */
class Actions_StackelFit : public Action_Finder{
    private:
        Stackel_Fitted_Potential *SFP;
    public:
        //! Actions_StackelFit constructor.
        /*!
          \param pot Potential_JS (axisymmetric)

        */
        Actions_StackelFit(Potential_JS *Pot);
        //! Finds actions
        /*!
          \param x phase-space point (x,v)
          \param params -- can pass alpha value to override Delta estimation
        */
        VecDoub actions(const VecDoub& x, void*params=nullptr);
        //! Finds angles
        /*!
          \param x phase-space point (x,v)
          \param params -- can pass alpha value to override Delta estimation

          \return actions -- 3D vector J=(J_R,J_phi,J_z)
        */
        VecDoub angles(const VecDoub& x, void*params=nullptr);
        //! Finds angles and hessian
        /*!
          \param x phase-space point (x,v)

          \return angles and frequencies --
                6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
        */
        VecDoub angles_with_hessdet(VecDoub x);
};

#endif
// ============================================================================

