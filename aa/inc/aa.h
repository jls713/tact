// ============================================================================
/// \file inc/aa.h
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
/// \brief Defines base class for action finders
///
/// Action_Finder takes a potential class and a minimum Action_Finder class
/// has routines
/// VecDoub actions(const VecDoub& x, void *params=nullptr)
/// VecDoub angles(const VecDoub& x, void *params=nullptr)
/// where params gives the possibility of passing options to the routines.
/// These routines accept a 6D phase-space point (x) and
/// return J=(JR,J_\phi,Jz) from actions and
/// return A=(\theta_R,\theta_\phi,\theta_z,
///           \Omega_R\Omega_\phi,\Omega_z) from angles
///
// ============================================================================

#ifndef AA_H
#define AA_H

#include "utils.h"
#include "potential.h"

// ============================================================================

/*! General base class for Action_Finder */
class Action_Finder{
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

// ============================================================================
