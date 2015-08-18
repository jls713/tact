#ifndef AA_H
#define AA_H

#include "utils.h"
#include "potential.h"

class Action_Finder{
    /* General base class for Action_Finder */
    public:
        inline virtual void reset(Potential_JS *Pot){std::cerr<<"You shouldn't be in Action_Finder reset aa.h\n";};
        inline virtual void partial_reset(Potential_JS *Pot){std::cerr<<"You shouldn't be in Action_Finder partial reset aa.h\n";};
        inline virtual void reset_sph(SphericalPotential *Pot){std::cerr<<"You shouldn't be in Action_Finder reset aa.h\n";};
        inline virtual VecDoub actions(const VecDoub& x, void *params=nullptr){
            /* actions for Cartesian position x */
            std::cerr<<"Action routine not overridden in aa.h"<<std::endl;
            return VecDoub(3,0);
        }
	inline VecDoub actions_nopars(const VecDoub& x){return actions(x);}
    inline VecDoub angles_nopars(const VecDoub& x){return angles(x);}
        inline virtual VecDoub angles(const VecDoub& x, void *params=nullptr){
            /* angles for Cartesian position x */
            std::cerr<<"Action routine not overridden in aa.h"<<std::endl;
            return VecDoub(6,0);
        }
        inline virtual double L_circ(double R){
            std::cout<<"You shouldn't be in L_circ"<<std::endl;return 0.;}
};
#endif
