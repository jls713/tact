// ============================================================================
// Action calculation for a few analytic cases
// ============================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface/GSLInterface.h"
#include "utils.h"
#include "potential.h"
#include "analytic_aa.h"
#include "debug.h"
#include "orbit.h"

// ============================================================================
// Triaxial harmonic oscillator
// ============================================================================

VecDoub Actions_HarmonicOscillator::actions(const VecDoub& x, void *params){
    VecDoub J(3,0);
    for(int i=0;i<3;i++)
        J[i] = (x[i+3]*x[i+3]+Om[i]*Om[i]*x[i]*x[i])/(2.*Om[i]);
    return J;
}

VecDoub Actions_HarmonicOscillator::angles(const VecDoub& x, void *params){
    VecDoub Theta(3,0);
    for(int i=0;i<3;i++){
        Theta[i] = atan(-x[3+i]/Om[i]/x[i])+(x[i]>0.?PI:0.);
        if(Theta[i]<0.)Theta[i]+=2.*PI;
        if(Theta[i]>2.*PI)Theta[i]-=2.*PI;
    }
    return Theta;
}

// ============================================================================
// Isochrone
// ============================================================================

VecDoub Actions_Isochrone::actions(const VecDoub& x, void *params){
    double E=H(x), L_ = L(x), Lz_ = Lz(x);
    if(E>0.)
        return {std::numeric_limits<double>::infinity(),0,0};
    return {GM/sqrt(-2*E)-0.5*(L_+sqrt(L_*L_+4.*GM*b)),Lz_,L_-fabs(Lz_)};

}

VecDoub Actions_Isochrone::freq(const VecDoub& x){
    double LL = L(x);
    double Omegar=pow(-2*H(x),1.5)/GM;
    double Omegap=0.5*Omegar*(1.+LL/sqrt(LL*LL+4*GM*b));
    return {Omegar,Omegap};
}

VecDoub Actions_Isochrone::Hessian(const VecDoub& x){
    VecDoub Om = freq(x);
    double D_rr,D_rp,D_pp;
    double LL = L(x);
    D_rr=-Om[0]*3.*sqrt(-2*H(x))/GM;
    D_rp=D_rr*Om[1]/Om[0];
    D_pp=Om[0]*2*GM*b*pow(LL*LL+4*GM*b,-1.5)
            +0.5*(1+LL/sqrt(LL*LL+4*GM*b))*D_rp;
    return {D_rr,D_rp,D_pp};
}

static double F(double x, double y){
    if (y>.5*PI) return .5*PI-atan(tan(.5*(PI-y))/x);
    else if (y<-.5*PI) return -.5*PI+atan(tan(.5*(PI+y))/x);
    else return atan(x*tan(0.5*y));
}

VecDoub Actions_Isochrone::angles(const VecDoub& x, void *params){
    VecDoub Sph = conv::CartesianToSphericalPolar(x);
    VecDoub Theta(3,0);
    double E=H(x), L_ = L(x), Lz_ = Lz(x);
    double c=GM/(-2*E)-b;
    double e=sqrt(1-L_*L_*(1+b/c)/GM/c);
    double r = Sph[0];
    double eta=atan2(r*Sph[3]/sqrt(-2.*E),b+c-sqrt(b*b+r*r));
    double OmR=pow(-2*E,1.5)/GM;
    double Omp=0.5*OmR*(1+L_/sqrt(L_*L_+4*GM*b));
    Theta[0]=eta-e*c*sin(eta)/(c+b);

    double psi = PI/2.;
    if(fabs(Sph[5])>1e-10) psi=atan2(cos(Sph[2]),-sin(Sph[2])*r*Sph[5]/L_);
    double a=sqrt((1+e)/(1-e));
    double ap=sqrt((1+e+2*b/c)/(1-e+2*b/c));

    Theta[2]=psi+Omp*Theta[0]/OmR-F(a,eta)-F(ap,eta)/sqrt(1+4*GM*b/L_/L_);

    double LR=Lz_/L_;
    double sinu = LR/sqrt(1.-LR*LR)/tan(Sph[2]);
    double u = 0;
    if(sinu>1.)
        u=.5*PI;
    else if(sinu<-1.)
        u = -.5*PI;
    else
        u = asin(sinu);
    if(Sph[5]>0.)
        u=PI-u;
    Theta[1]=Sph[1]-u+SIGN(Lz_)*Theta[2];
    for(int i=0;i<3;i++){
        while(Theta[i]<0.) Theta[i]+=2.*PI;
        while(Theta[i]>2.*PI) Theta[i]-=2.*PI;
    }
    return Theta;
}

// ============================================================================
// ============================================================================


/*

def deltaH_ho(omega,xsamples):
    if(np.any(omega<1e-5)):
        return np.nan
    H = 0.5*np.sum(xsamples.T[3:]**2,axis=0)+0.5*np.sum((omega[:3]*xsamples.T[:3].T)**2,axis=1)
    return H-np.mean(H)

def Jac_deltaH_ho(omega,xsamples):
    dHdparams = omega[:3]*xsamples.T[:3].T**2
    return dHdparams-np.mean(dHdparams,axis=0)

def findbestparams_ho(xsamples):
    """ Minimize sum of square differences of H_sho-<H_sho> for timesamples """
    return np.abs(leastsq(deltaH_ho,np.array([10.,10.,10.]), Dfun = Jac_deltaH_ho, args=(xsamples,))[0])[:3]
// */

// int main(){
//     VecDoub X = {8.,0.,0.,0.,200.,40.};
//     Actions_Isochrone ActI(2e6,6.);
//     printVector(ActI.actions(X));
//     Isochrone Iso(2e6,6.);
//     Orbit O(&Iso);
//     VecDoub Q = O.integrate(X,1.,1e-6);
//     printVector(ActI.actions(X));
//     return 0;
// }
