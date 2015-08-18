#ifndef STACK_FIT_H
#define STACK_FIT_H

#include "utils.h"
#include "potential.h"

// ============================================================================

class Stackel_Fitted_Potential: public StackelProlate_PerfectEllipsoid{
    private:
        const int DATAPOINTS = 40; // NUmber of gridpoints
        const int STEPMAX = 300; // Max number of integrations
        const int NMAX = 10; // Number of bisections to perform
        const double cl = 4.;
        const double cv = 0.5;
        VecDoub limits;
        Potential_JS *TruePot;
        std::unique_ptr<Orbit> Orb;
        double *lambdagrid,*nugrid,*flamgrid,*fnugrid,*y2Lamgrid,*y2Nugrid;
        double L, N, chibar;
    public:
        Stackel_Fitted_Potential(Potential_JS *TruePot)
            :StackelProlate_PerfectEllipsoid(10.,-30.),TruePot(TruePot){
                Orb = std::unique_ptr<Orbit>(new Orbit(TruePot,1e-6));
            };
        ~Stackel_Fitted_Potential(){
            if(lambdagrid){
                delete lambdagrid;delete nugrid;delete flamgrid;
                delete fnugrid; delete y2Lamgrid; delete y2Nugrid;
            }
        }

        double Vlamnu(double lambda, double nu);
        double chi(double lambda, double nu);
        double BigNu(double nu);
        double BigLambda(double lambda);
        double flamINT(double lam);
        double fnuINT(double nu);
        double chibarint(double lambda,double nu);
        double flamINT_DERIV(double lam);
        double fnuINT_DERIV(double nu);
        double G(double tau);
        double BigF(double tau);
        double BigFPrime(double tau);
        double GPrime(double tau);
        // double Phi(VecDoub x);
        double fnu(double);
        double flam(double);
        // VecDoub Forces(VecDoub x);
        double find_turning(VecDoub x, VecDoub y, double step, int lamnu);
        void find_limits(VecDoub x);
        void fit_potential(VecDoub x);
};

// ============================================================================

class Actions_StackelFit : public Action_Finder{

    private:
        Stackel_Fitted_Potential *SFP;
    public:
        Actions_StackelFit(Potential_JS *Pot);
        VecDoub actions(const VecDoub& x, void*params=nullptr);
        VecDoub angles(VecDoub x, bool with_hess=false);
        VecDoub angles_with_hessdet(VecDoub x);
};

#endif
// ============================================================================

