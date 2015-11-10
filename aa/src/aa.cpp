#include "aa.h"
#include "spherical_aa.h"
#include <limits>

VecDoub planar_sphr_actions(const VecDoub &x,Potential_JS* Pot){
    PlanarAxisymPotential PP(Pot);
    Actions_Spherical AS(&PP);
    return AS.actions(x);
}
VecDoub planar_sphr_angles(const VecDoub &x,Potential_JS* Pot){
    PlanarAxisymPotential PP(Pot);
    Actions_Spherical AS(&PP);
    return AS.angles(x);
}

/*!
    Checks input for actions
*/
int action_check(const VecDoub &x, VecDoub &acts, Potential_JS *Pot){
    // Check length of x
    if(x.size()<6){
        std::cerr<<"Must pass 6D vector to action routine\n";
        for(unsigned i=0;i<3;++i) acts[i]=std::numeric_limits<double>::infinity();
        return 1;
    }
    // Check if unbound
    if(Pot->H(x)>Pot->Phi({1.e7*x[0],1.e7*x[1],1.e7*x[2]})){
        std::cerr<<"Orbit passed to action routine is unbound\n";
        for(unsigned i=0;i<3;++i)
            acts[i]=std::numeric_limits<double>::infinity();
        return 1;
    }

    VecDoub Polar = conv::CartesianToPolar(x);
    if(fabs(Polar[0])<SMALL){
        acts[0]=0.;acts[1]=0.;acts[2]=0.;
        return 1;
    }
    acts[1]=Polar[0]*Polar[4];
    if(fabs(x[2])<SMALL and fabs(x[5])<SMALL){
        acts[2]=0.;
        if(fabs(Polar[3])<SMALL){
            acts[0]=0.;
            return 1;
        }
        else{
            acts=planar_sphr_actions(x,Pot);
            acts[2]=0.;
            return 1;
        }
    }
    return 0;
}
/*!
    Checks input for angles
*/
int angle_check(const VecDoub &x, VecDoub &angs, Potential_JS *Pot){
    // Check length of x
    if(x.size()<6){
        std::cerr<<"Must pass 6D vector to angle routine\n";
        for(unsigned i=0;i<6;++i) angs[i]=std::numeric_limits<double>::infinity();
        return 1;
    }
    // Check if unbound
    if(Pot->H(x)>Pot->Phi({1.e7*x[0],1.e7*x[1],1.e7*x[2]})){
        std::cerr<<"Orbit passed to angle routine is unbound\n";
        for(unsigned i=0;i<6;++i)
            angs[i]=std::numeric_limits<double>::infinity();
        return 1;
    }
    VecDoub Polar = conv::CartesianToPolar(x);
    if(fabs(Polar[0])<SMALL){
        angs[0]=0.;angs[1]=0.;angs[2]=0.;
        angs[3]=0.;angs[4]=0.;angs[5]=0.;
        return 1;
    }
    if(fabs(x[2])<SMALL and fabs(x[5])<SMALL){
        angs[2]=0.;
        angs[5]=0.;
        if(fabs(Polar[3])<SMALL){
            angs[0]=0.;
            angs[3]=0.;
            angs[1]=Polar[1];
            angs[4]=Polar[4]/Polar[0];
            return 1;
        }
        else{
            angs=planar_sphr_angles(x,Pot);
            angs[2]=0.;
            angs[5]=0.;
            return 1;
        }
    }
    return 0;
    }
