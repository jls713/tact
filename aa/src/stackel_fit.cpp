// ============================================================================
/// \file src/stackel_fit.cpp
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

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "jamestools/numrec/press.h"
#include "GSLInterface/GSLInterface.h"
#include "utils.h"
#include "coordtransforms.h"
#include "potential.h"
#include "spherical_aa.h"
#include "debug.h"
#include "orbit.h"
#include "stackel_fit.h"
#include "cuba/cuba.h"
#include "stackel_aa.h"
#include "aa.h"

struct chi_integrals{
    Stackel_Fitted_Potential *SFP;
    double tau;
    VecDoub x2min,x2max;
    chi_integrals(Stackel_Fitted_Potential *SFP, double tau,VecDoub x2min={0.}, VecDoub x2max={0.}): SFP(SFP), tau(tau),x2min(x2min),x2max(x2max){};
};

static double *dmatrix(int n){// creates array of length n and fills with zeroes
    double *m1 = new double[n];
    for(int i=0;i<n;i++) m1[i]=0;
    return m1;
}

static int ChiBarIntCubature(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata) {
    chi_integrals * SFP = (chi_integrals *) fdata;
    double y2[2];
    for(int i=0;i<2;i++)
        y2[i]=(SFP->x2max[i]-SFP->x2min[i])*y[i]+SFP->x2min[i];
     // interface to integrand for chibar for use with cuba code
    fval[0] = SFP->SFP->chibarint(y2[0],y2[1]);
    return 0;
}

static double chiNu(double nu, void *params){
    chi_integrals * SFP = (chi_integrals *) params;
    return SFP->SFP->chi(SFP->tau,nu)*SFP->SFP->BigNu(nu);
}

static double chiLam(double lambda, void *params){
    chi_integrals * SFP = (chi_integrals *) params;
    return SFP->SFP->chi(lambda,SFP->tau)*SFP->SFP->BigLambda(lambda);
}

double Stackel_Fitted_Potential::Vlamnu(double lambda, double nu){
// fitpot to be fit to as a function of tau coordinates
    VecDoub x = tau2x({lambda,0.,nu});
    return TruePot->Phi(x);
}

double Stackel_Fitted_Potential::chibarint(double lambda,double nu){
    return chi(lambda,nu)*BigLambda(lambda)*BigNu(nu);
    }


double Stackel_Fitted_Potential::chi(double lambda, double nu){
    double V = Vlamnu(lambda,nu);
    return -(lambda-nu)*V;
}

double Stackel_Fitted_Potential::BigNu(double nu){
    //Nu weight - see Paper 1 with d=0
    return pow((nu+gamma())/fabs(gamma()-alpha()),cv);
}

double Stackel_Fitted_Potential::BigLambda(double lambda){
//Lambda weight - see Paper 1 with d=0
    return pow((1.0/(lambda+alpha()-gamma())),cl);
}

double Stackel_Fitted_Potential::flam(double lambda){
    chi_integrals SFP(this,lambda);
    return ((GaussLegendreQuad(&chiNu,limits[2],limits[3],&SFP,50)/N)-0.5*chibar);
}

double Stackel_Fitted_Potential::fnu(double nu){
    if(nu+gamma()<TINY) nu+=TINY;
    chi_integrals SFP(this,nu);
    return (-(GaussLegendreQuad(&chiLam,limits[0],limits[1],&SFP,50)/L)+0.5*chibar);
}

double Stackel_Fitted_Potential::flamINT(double lam){
    if(lam<lambdagrid[0] || lam>lambdagrid[DATAPOINTS-1])
        return flam(lam);
    else return splint(lambdagrid,flamgrid,y2Lamgrid,DATAPOINTS,lam);
}

double Stackel_Fitted_Potential::fnuINT(double nu){
    if(nu<nugrid[0] || nu>nugrid[DATAPOINTS-1])
        return fnu(nu);
    else return splint(nugrid,fnugrid,y2Nugrid,DATAPOINTS,nu);
}

double Stackel_Fitted_Potential::flamINT_DERIV(double lam){
    if(lam<lambdagrid[0] || lam>lambdagrid[DATAPOINTS-1]){
        return (flamINT(lam+0.01)-flamINT(lam-0.01))/0.02;
        }
    else{
    double flamInter = splintp(lambdagrid,flamgrid,y2Lamgrid,DATAPOINTS,lam);
    return flamInter;
    }
}

double Stackel_Fitted_Potential::fnuINT_DERIV(double nu){
    if(nu<nugrid[0] || nu>nugrid[DATAPOINTS-1]){
        return (fnuINT(nu+0.001)-fnuINT(nu-0.001))/0.002;
        }
    else{
    double fnuInter = splintp(nugrid,fnugrid,y2Nugrid,DATAPOINTS,nu);
    return fnuInter;
    }
}

double Stackel_Fitted_Potential::G(double tau){
  // de Zeeuw G(tau)
        if(tau==-gamma()){tau=tau+TINY;}
        if((tau>-alpha() and tau<-gamma() and -gamma()>-alpha()) or (tau>-alpha() and -   alpha()>-gamma()))
            {   return flamINT(tau)/(tau+gamma());}
        else{   return fnuINT(tau)/(tau+gamma());}
}


  double Stackel_Fitted_Potential::BigF(double tau){
  // de Zeeuw F(tau)
        if(tau==-gamma()){tau=tau+TINY;}
        if((tau>-alpha() and tau<-gamma() and -gamma()>-alpha()) or (tau>-alpha() and -   alpha()>-gamma()))
                {   return flamINT(tau)*(tau+alpha());}
            else{return fnuINT(tau)*(tau+alpha());}
        }

double Stackel_Fitted_Potential::BigFPrime(double tau){
// de Zeeuw F'(tau)
  if((tau>-alpha() and tau<-gamma() and -gamma()>-alpha()) or (tau>-alpha() and -alpha()> -gamma()))
      {return flamINT_DERIV(tau)*(tau+alpha())+flamINT(tau);}
  else{return fnuINT_DERIV(tau)*(tau+alpha())+fnuINT(tau);}
  }

double Stackel_Fitted_Potential::GPrime(double tau){
    //    computes G'(tau) - necessary for forces
    return (BigFPrime(tau)-G(tau)*(2.0*tau+alpha()+gamma()))/((tau+gamma())*( tau+alpha())) ;
}

// double Stackel_Fitted_Potential::Phi(VecDoub x){
//     return 0.;
// }

// VecDoub Stackel_Fitted_Potential::Stackel_Fitted_Potential::Forces(VecDoub x){
//     return {0.,0.,0.};
// };

double Stackel_Fitted_Potential::find_turning(VecDoub x, VecDoub y, double step, int lamnu){
    int index = 3;
    if(lamnu) index = 5;

    double STEP = step/2.0; int N = 0;
    VecDoub ylower=x, taulower = xv2tau(ylower);
    VecDoub yupper=y, tauupper = xv2tau(ylower);
    // Carry out N bisections
    while(N<NMAX){
        // Carry out integration over half the interval
        yupper = Orb->integrate(ylower,STEP,STEP);
        tauupper = xv2tau(yupper);
        // Check which half root occurs in
        if(tauupper[index]*taulower[index]>0.0){
            ylower=yupper;
            taulower=tauupper;
        }
        STEP = STEP/2.0;
        N++;
    }
    return tauupper[index-3];
}

void Stackel_Fitted_Potential::find_limits(VecDoub x){
    // This function begins by integrating the initial phase space point with an adaptive timestep
    // to an accuracy of 0.01% and finding an estimate of alpha using the first few points.
    // We then determine the three points which define the fitting region - lambdaminus, lambdaplus
    // and nuplus which are defined by taudot=0.
    // We perform the potential fit over the region defined by these three points.
    // At the end, we fill the arrays to be used in the interpolation of the f functions.
    limits = {0.,0.,0.,0.};
    VecDoub y = x, ynew, tau, taunew;
    double step = 0.01*TruePot->torb(x);
    int steps=0; // counts no. of int. steps
    // flags for checking when three edge points are found
    int alldone1=0,alldone2=0,alldone3=0,alldone4=0;
    newalpha(-TruePot->DeltaGuess(y)+gamma()); // Initial guess of alpha
    if(alpha()>gamma() or alpha()!=alpha() )newalpha(gamma()-0.1);
    tau = xv2tau(y);
    int n = 2; // counts no. of alpha guesses+1
    while (alldone1*alldone2*alldone3*alldone4==0){

        ynew = Orb->integrate(y,step,step);
        // Refine our guess of alpha
        if(steps%5==0 && steps<=STEPMAX){
            double new_alpha = -TruePot->DeltaGuess(ynew)+gamma();
            if(new_alpha>gamma()){
                new_alpha=gamma()-0.1;
                if(debug_NegativeDelta)
                    std::cerr<<"Negative Delta at R="<<sqrt(y[0]*y[0]+y[1]*y[1])<<", z="<<y[2]<<std::endl;
            }
            newalpha(((n-1)*alpha()+new_alpha)/(double)n);
            n++;
        }
        taunew = xv2tau(ynew);
        double nudotdot = (taunew[5]-tau[5]);
        double lamdotdot = (taunew[3]-tau[3]);
        // Edge is defined by points when lambdadot or nudot change sign

        // Inner Lambda Edge
        if(taunew[3]*tau[3]<0.0
           && lamdotdot>0.0
           && steps>STEPMAX
           && alldone1==0){
            limits[0] = find_turning(y,ynew,step,0);alldone1=1;
        }

        // Outer Lambda Edge - same algorithm as above
        if(taunew[3]*tau[3]<0.0
           && lamdotdot<0.0
           && steps>STEPMAX
           && alldone2==0){
            limits[1] = find_turning(y,ynew,step,0);alldone2=1;
        }

        if(taunew[5]*tau[5]<0.0
           && nudotdot<0.0
           && steps>STEPMAX
           && alldone3==0){
            limits[3] = find_turning(y,ynew,step,1); alldone3=1;
        }

        if(taunew[5]*tau[5]<0.0
           && nudotdot>0.0
           && steps>STEPMAX
           && alldone4==0){
            limits[2] = find_turning(y,ynew,step,1);alldone4=1;
        }
        y = ynew; tau = taunew;
        steps++;
    }
    // if(y[0]>300.) limits = {-1.,-1.,-1.,-1.};
    return;
}

void Stackel_Fitted_Potential::fit_potential(VecDoub x){

    find_limits(x);
    // limits[3]*=1.02; limits[0]*=0.98; limits[1]*=1.02;
    if(limits[0]>limits[1]){
        double tautmp=limits[0];
        limits[0]=limits[1];limits[1]=tautmp;
    }
    if(fabs(limits[2]-limits[3])<SMALL) limits[3]+=SMALL;
    // std::cout<<alpha()<<std::endl;
    // Fit //
    // Perform integrals which are indep. of R,z
    // can do integrals analytically rather than numerically
    double B = alpha()-gamma();
    L = (1./(cl-1.))*(pow(limits[0]+B,-(cl-1))-pow(limits[1]+B,-(cl-1)));
    N = (pow((limits[3]+gamma()),cv+1))/((cv+1.)*pow(fabs(gamma()-alpha()),cv));
    VecDoub xmin = {limits[0],limits[2]}, xmax = {limits[1],limits[3]};
    double integral[1],error[1],prob[1];

    chi_integrals SFP(this,0.,xmin,xmax);
    double prod = 1.;int neval,fail,nregions;
    for(int i=0;i<2;i++) prod*=(SFP.x2max[i]-SFP.x2min[i]);

    Cuhre(2,1,ChiBarIntCubature,&SFP,1,1e-6,0,0,
    MINEVAL, MAXEVAL, 0, STATEFILE,SPIN,
    &nregions, &neval, &fail, integral, error, prob);

    integral[0]*=prod;
    chibar = integral[0]/(L*N);

    // Grids for interpolation
    lambdagrid=dmatrix(DATAPOINTS);
    nugrid=dmatrix(DATAPOINTS);
    flamgrid=dmatrix(DATAPOINTS);
    fnugrid=dmatrix(DATAPOINTS);

    // Fill grids
    for(int i=0;i<DATAPOINTS;i++){
        lambdagrid[i] = limits[0]+(double)i*(limits[1]-limits[0])/(double)(DATAPOINTS-1);
        nugrid[i] = limits[2]+(double)i*(limits[3]-limits[2])/(double)(DATAPOINTS-1);
        flamgrid[i] = flam(lambdagrid[i]);
        fnugrid[i] = fnu(nugrid[i]);
    }

    // Create splines for interpolation
    y2Lamgrid = dmatrix(DATAPOINTS);
    y2Nugrid = dmatrix(DATAPOINTS);
    spline(lambdagrid, flamgrid, DATAPOINTS, 1e30, 1e30, y2Lamgrid);
    spline(nugrid, fnugrid, DATAPOINTS, 1e30, 1e30, y2Nugrid);
}

Actions_StackelFit::Actions_StackelFit(Potential_JS *Pot, double eps){
    SFP = new Stackel_Fitted_Potential(Pot,eps);
};

VecDoub Actions_StackelFit::actions(const VecDoub& x, void*params){
    VecDoub acts(3,0.);
    if(action_check(x,acts,SFP->potential())) return acts;
    SFP->fit_potential(x);
    Actions_AxisymmetricStackel AS(SFP);
    return AS.actions(x);
}
VecDoub Actions_StackelFit::angles(const VecDoub& x, void*params){
    VecDoub angs(6,0.);
    if(angle_check(x,angs,SFP->potential())) return angs;
    SFP->fit_potential(x);
    Actions_AxisymmetricStackel AS(SFP);
    return AS.angles(x,params);
}

VecDoub Actions_StackelFit::angles_with_hessdet(VecDoub x){
    SFP->fit_potential(x);
    Actions_AxisymmetricStackel AS(SFP);
    int bh=1;
    VecDoub aa = AS.angles(x,&bh);
    VecDoub bb;
    for(int i=0;i<6;i++)
        bb.push_back(aa[i]);
    bb.push_back(aa[6]*det2(aa[10],aa[11],aa[13],aa[14])
                - aa[7]*det2(aa[9],aa[11],aa[12],aa[14])
                + aa[8]*det2(aa[9],aa[10],aa[12],aa[13]));
    return bb;
}

// int main(){

//     GalPot L("/home/jls/work/code/Torus/pot/PJM11_best.Tpot");
//     // StackelProlate_PerfectEllipsoid L(.01,-20.);
//     VecDoub x = {3.756020000e+00,1.489170000e+01,-5.674940000e-01,-2.520260000e-01,4.498560000e-03,1.339060000e-01};
//     for(unsigned i=3;i<6;++i)x[i]*=977.775;
//     // Actions_AxisymmetricStackel AASS(&L);
//     // printVector(AASS.actions(x));
//     // Actions_AxisymmetricStackel_Fudge AAS(&L,-20.);
//     Actions_StackelFit AS(&L);
//     Orbit orb(&L);
//     double torb = L.torb(x);
//     VecDoub y = orb.integrate(x,10.*torb,0.1*torb);
//     for(auto i:orb.results()){
//         // printVector(concatVectors(AAS.actions(x),AS.actions(x)));
//         // printVector(AS.angles(x,true));
//         // AS.angles(x,false);
//         printVector(AS.actions(i)*(1./977.775));
//     }
//     return 0;
// }
