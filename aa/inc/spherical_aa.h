// ============================================================================
/// \file inc/spherical_aa.h
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
/// \brief Action finding in spherical potentials
///
/// In all spherical potentials the actions are J_R, L_z, L-|L_z| where
/// J_R is a 1D quadrature. Actions_Spherical provides routines for performing
/// the 1D quadratures.
///
//============================================================================

#ifndef SPHERICAL_AA
#define SPHERICAL_AA

#include "potential.h"
#include "aa.h"

/*! Action finding in spherical potential */
class Actions_Spherical : public Action_Finder{
    private:
        VecDoub find_limits(double r, double E, double L);
        SphericalPotential *Pot;
        double dr0dH(double r, double L);
        double dr0dL(double r, double L);
    public:
        //! Actions_Spherical constructor
        /*!
          \param Pot Spherical potential (SphericalPotential)
        */
        Actions_Spherical(SphericalPotential *Pot): Pot(Pot){};
        //! reset potential
        inline void reset_sph(SphericalPotential *pot){Pot=pot;}
        //! Finds actions
        /*!
          \param x phase-space point (x,v)
          \param params -- does nothing

          \return actions -- 3D vector J=(J_R,J_phi=L_z,J_z=L-|L_z|)
        */
        VecDoub actions(const VecDoub &x, void *params=nullptr);
        //! Finds angles
        /*!
          \param x phase-space point (x,v)
          \param params -- does nothing

          \return angles --
                3D vector (theta_R,theta_phi,theta_z)
        */
        VecDoub angles(const VecDoub &x, void *params=nullptr){
            VecDoub X = angles_and_freqs(x);
            X.resize(3);
            return X;
        }
        //! Finds angles and frequencies
        /*!
          \param x phase-space point (x,v)

          \return angles and frequencies --
                6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)
        */
        VecDoub angles_and_freqs(const VecDoub &x);
        //! Finds Hessian
        /*!
          \param x phase-space point (x,v)
            // currently not working
        */
        VecDoub Hessian(const VecDoub &x);
};

/*! Helper structure for finding limits of orbit in spherical potential */
struct Actions_Spherical_limits_struct{
    SphericalPotential *Pot;
    double E, L;
    Actions_Spherical_limits_struct(SphericalPotential *Pot, double E, double L)
    : Pot(Pot), E(E), L(L){}
};

/*! Helper structure for integrating action integrand in spherical potential */
struct Actions_Spherical_data_struct{
    SphericalPotential *Pot;
    double E, L, Lz, Delta, taubar, dDelta, dtaubar;
    Actions_Spherical_data_struct(SphericalPotential *Pot, double E, double L, double Lz, double Delta, double taubar, double dDelta=0., double dtaubar=0.)
    : Pot(Pot), E(E), L(L), Lz(Lz), Delta(Delta), taubar(taubar), dDelta(dDelta), dtaubar(dtaubar){}
};


#endif

//============================================================================
