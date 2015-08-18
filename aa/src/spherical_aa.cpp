/*======================================*/
/*Actions, Angles, Frequencies & Hessian*/
/*======================================*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface/GSLInterface.h"
#include "utils.h"
#include "coordtransforms.h"
#include "potential.h"
#include "spherical_aa.h"
#include "debug.h"
#include "orbit.h"

static double p_r(double r, void *params){
    // need to set taubargl and Deltagl
    Actions_Spherical_limits_struct * AS = (Actions_Spherical_limits_struct *) params;
    return 2.*AS->E-2.*AS->Pot->Phi_r(r)-AS->L*AS->L/r/r;
}

static double J_r(double theta,void *params){
    // need to set taubargl and Deltagl
    Actions_Spherical_data_struct * AS = (Actions_Spherical_data_struct *) params;
    double r = sin(theta)*AS->Delta+AS->taubar;
    return cos(theta)*sqrt(MAX(0.,2.*AS->E-2.*AS->Pot->Phi_r(r)-AS->L*AS->L/r/r));
}

static double dJrdH(double theta,void *params){
    // need to set taubargl and Deltagl
    Actions_Spherical_data_struct * AS = (Actions_Spherical_data_struct *) params;
    double r = sin(theta)*AS->Delta+AS->taubar;
    return cos(theta)/sqrt(MAX(TINY,2.*AS->E-2.*AS->Pot->Phi_r(r)-AS->L*AS->L/r/r));
}

static double dJrdL(double theta,void *params){
    // need to set taubargl and Deltagl
    Actions_Spherical_data_struct * AS = (Actions_Spherical_data_struct *) params;
    double r = sin(theta)*AS->Delta+AS->taubar;
    return -cos(theta)*AS->L/r/r/sqrt(MAX(TINY,2.*AS->E-2.*AS->Pot->Phi_r(r)-AS->L*AS->L/r/r));
}

static double dLdL(double theta,void *params){
    // need to set taubargl and Deltagl
    Actions_Spherical_data_struct * AS = (Actions_Spherical_data_struct *) params;
    double st = sin(theta);
    return AS->L/sqrt(MAX(TINY,AS->L*AS->L-AS->Lz*AS->Lz/st/st));
}

static double d2JdH2(double theta,void *params){
    // need to set taubargl and Deltagl
    Actions_Spherical_data_struct * AS = (Actions_Spherical_data_struct *) params;
    double r = sin(theta)*AS->Delta+AS->taubar;
    double L2r2 = AS->L*AS->L/r/r;
    double P = sqrt(MAX(TINY,2.*AS->E-2.*AS->Pot->Phi_r(r)-L2r2));
    // double drdH = AS->dDelta*sin(theta)+AS->dtaubar;
    return cos(theta)*(-1./*+(AS->Pot->dPhi_r(r)-L2r2/r)*drdH*/)/pow(P,3);
}

static double d2JdHdL(double theta,void *params){
    // need to set taubargl and Deltagl
    Actions_Spherical_data_struct * AS = (Actions_Spherical_data_struct *) params;
    double r = sin(theta)*AS->Delta+AS->taubar;
    double L2r2 = AS->L*AS->L/r/r;
    double P = sqrt(MAX(TINY,2.*AS->E-2.*AS->Pot->Phi_r(r)-L2r2));
    double drdH = AS->dDelta*sin(theta)+AS->dtaubar;
    return cos(theta)*(AS->L/r/r+(AS->Pot->dPhi_r(r)-L2r2/r)*drdH)/pow(P,3);
}


VecDoub Actions_Spherical::find_limits(double r, double E, double L){
    VecDoub limits;
    Actions_Spherical_limits_struct Act(Pot,E,L);
    double r_in=r, r_out=r;
    root_find RF(SMALL,100);
    if(p_r(0.,&Act)>0) r_in=0.;
    else while(p_r(r_in,&Act)>=0.0) r_in*=0.9;
    while(p_r(r_out,&Act)>=0.0) r_out*=1.1;
    limits.push_back(RF.findroot(&p_r,r_in,r,&Act));
    limits.push_back(RF.findroot(&p_r,r,r_out,&Act));
    return limits;
}

VecDoub Actions_Spherical::actions(const VecDoub &x, void *params){
    double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    double E = Pot->H(x), L = Pot->L(x), Lz = Pot->Lz(x);
    if(E>0.){
        return {std::numeric_limits<double>::infinity(),Lz,L-fabs(Lz)};
    }
    VecDoub limits = find_limits(r,E,L);
    double taubar = 0.5*(limits[0]+limits[1]);
    double Delta = 0.5*(limits[1]-limits[0]);
    Actions_Spherical_data_struct Act(Pot,E,L,Lz,Delta,taubar);
    double JR = Delta*GaussLegendreQuad(&J_r,-PI/2.,PI/2.,&Act,8)/PI;
    return {JR,Lz,L-fabs(Lz)};
}

VecDoub Actions_Spherical::angles_and_freqs(const VecDoub &x){
    // call actions before
    double E = Pot->H(x), L = Pot->L(x), Lz = Pot->Lz(x);
    if(E>0.){
        return {std::numeric_limits<double>::infinity(),L};
    }
    double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    VecDoub limits = find_limits(r,E,L);
    double taubar = 0.5*(limits[0]+limits[1]);
    double Delta = 0.5*(limits[1]-limits[0]);
    Actions_Spherical_data_struct Act(Pot,E,L,Lz,Delta,taubar);

    double OmegaR = Delta*GaussLegendreQuad(&dJrdH,-PI/2.,PI/2.,&Act,8)/PI;
    OmegaR = 1./OmegaR;
    double Omegap = Delta*GaussLegendreQuad(&dJrdL,-PI/2.,PI/2.,&Act,8)/PI;
    Omegap*=-OmegaR;

    VecDoub SPol = conv::CartesianToSphericalPolar(x);

    double thetaLim = asin(MAX(-1.,MIN(1.,(SPol[0]-taubar)/Delta)));
    double dSRdH=sign(SPol[3])*Delta*GaussLegendreQuad(&dJrdH,-PI/2.,thetaLim,&Act,8);
    double dSRdL=sign(SPol[3])*Delta*GaussLegendreQuad(&dJrdL,-PI/2.,thetaLim,&Act);

    double ThetaR = dSRdH*OmegaR;

    double dStdL=sign(SPol[5])*GaussLegendreQuad(&dLdL,PI/2.,SPol[2],&Act,8);
    // printVector(x);
    double Thetap=dSRdL+dStdL+dSRdH*Omegap;
    double LR=fabs(Lz)/L;
    double sinu = LR/sqrt(1.-LR*LR)/tan(SPol[2]);
    double u = 0.;
    if(sinu>1.)
        u=PI/2.;
    else if(sinu<-1.)
        u = -PI/2.;
    else
        u = asin(sinu);
    if(SPol[5]>0.)
        u=PI-u;
    double Thetat=SPol[1]-u+sign(Lz)*Thetap;

    if(SPol[5]>0.)
        Thetat-=PI;

    if(SPol[5]<0.) Thetap+=PI;
    if(ThetaR>2.*PI) ThetaR-=2.*PI;
    if(ThetaR<0.) ThetaR+=2.*PI;
    if(Thetap>2.*PI) Thetap-=2.*PI;
    if(Thetap<0.) Thetap+=2.*PI;
    if(Thetat>2.*PI) Thetat-=2.*PI;
    if(Thetat<0.) Thetat+=2.*PI;
    if(Thetat<0.) Thetat+=2.*PI;

    return {ThetaR,Thetat,Thetap,OmegaR,sign(Lz)*Omegap,Omegap};
}

double Actions_Spherical::dr0dH(double r, double L){
    return 1./(Pot->dPhi_r(r)-L*L/r/r/r);
}
double Actions_Spherical::dr0dL(double r, double L){
    return -L/r/r/(Pot->dPhi_r(r)-L*L/r/r/r);
}

VecDoub Actions_Spherical::Hessian(const VecDoub &x){
    std::cerr<<"Currently not working!!"<<std::endl;
    VecDoub Freq(2,0), AF = angles_and_freqs(x);
    Freq[0]=AF[3];Freq[1]=AF[4];
    double E = Pot->H(x), L = Pot->L(x), Lz = Pot->Lz(x);
    double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    VecDoub limits = find_limits(r,E,L);
    double taubar = 0.5*(limits[0]+limits[1]);
    double Delta = 0.5*(limits[1]-limits[0]);
    double drpdH = dr0dH(limits[0],L);
    double dradH = dr0dH(limits[1],L);
    double dDelta=0.5*(dradH-drpdH);
    double dtaubar=0.5*(dradH+drpdH);
    Actions_Spherical_data_struct Act(Pot,E,L,Lz,Delta,taubar,dDelta,dtaubar);
    double d2JdH = Delta*GaussLegendreQuad(&d2JdH2,-PI/2.,PI/2.,&Act)/PI;
    double drpdL = dr0dL(limits[0],L);
    double dradL = dr0dL(limits[1],L);
    dDelta=0.5*(dradL-drpdL);
    dtaubar=0.5*(dradL+drpdL);
    // double d2JdHdL2 = Delta*GaussLegendreQuad(&d2JdHdL,-PI/2.,PI/2.,&Act)/PI;
    // d2JdH+=Delta*(dr0dH(limits[1])*dJrdH(PI/2.,&Act)-dr0dH(limits[0])*dJrdH(-PI/2.,&Act))/PI;
    double dOmegadJR = -Freq[0]*Freq[0]*Freq[0]*d2JdH;
    return {dOmegadJR};
}

void simple_isochrone_orbit_for_movie(void){
    IsochronePotential Iso(1000000.,3.64);
    Actions_Spherical AS(&Iso);
    VecDoub X ={2.67519,1.04903,-3.08583,175.546,-216.936,63.3228};
    Orbit orbit(&Iso);
    VecDoub QQ=X;double PPtmp = 0., fac = 0.;
    for(int i=0;i<250;i++){
        QQ=orbit.integrate(QQ, 0.004,0.004);
        VecDoub PP2 = AS.actions(QQ);
        VecDoub PP = AS.angles_and_freqs(QQ);
        if(PP[0]<PPtmp){PPtmp = PP[0];fac+=1.;}
        else PPtmp = PP[0];
        std::cout<<i*0.004<<" "<<QQ[0]<<" "<<QQ[1]<<" "<<PP2[0]<<" "<<PP2[1]<<" "<<PP[0]+fac*2*PI<<" "<<PP[1]<<std::endl;
    }
}

// int main(){
//     IsochronePotential Iso(1000000.,3.64);
//     // NFWSpherical Iso(1000000.,3.64);
//     Actions_Spherical AS(&Iso);
//     VecDoub X ={2.67519,1.04903,-3.08583,10.,50.,10.};
//     // VecDoub X = {-0.779383,0.586109,0.00117222,-0.263395,-0.443456 };
//     // printVector(AS.actions(X));
//     // printVector(AS.angles_and_freqs(X));
//     // printVector(Iso.Omega(X));
//     // printVector(Iso.Hessian(X));
//     // printVector(AS.Hessian(Iso.Omega(X)));
//     // return 0;
//     // printVector(AS.angles_and_freqs(X));
//     Orbit orbit(&Iso);
//     VecDoub QQ=X;double PPtmp = 0., fac = 0.;
//     for(int i=0;i<250;i++){
//         QQ=orbit.integrate(QQ, 0.004,0.004);
//         // printVector(AS.actions(QQ));
//         VecDoub PP2 = AS.actions(QQ);
//         VecDoub PP = AS.angles_and_freqs(QQ);
//         if(PP[0]<PPtmp){PPtmp = PP[0];fac+=1.;}
//         else PPtmp = PP[0];
//         std::cout<<i*0.004<<" "<<QQ[0]<<" "<<QQ[1]<<" "<<PP2[0]<<" "<<PP2[1]<<" "<<PP[3]<<" "<<PP[4]<<std::endl;
//         // printVector(AS.Hessian({QQ[3],QQ[4]}));
//     }
// }
