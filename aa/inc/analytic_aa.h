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
/// \brief Action finders for potentials in which actions are analytic
///
/// Defines action-finding routines for the two potentials in which the actions
/// , angles and frequencies are analytic. These are
/// 1. Isochrone potential given by Phi(r) = -GM/(b+sqrt(b**2+r**2))
/// 2. 3D harmonic oscillator given by Phi(\bs{x})=0.5 \omega_i^2 x^2
// ============================================================================

#ifndef ANALYTIC_AA_H
#define ANALYTIC_AA_H

#include "potential.h"
#include "aa.h"

// ============================================================================
/*! Action finder for 3D harmonic oscillator
    Finds actions in potential Phi(x,y,z)=0.5*((Om_x*x)^2+(Om_y*y)^2+(Om_z*z)^2
*/
class Actions_HarmonicOscillator : public Action_Finder, HarmonicOscillator{
    private:
    public:
        //! Actions_HarmonicOscillator constructor.
        /*!
          \param Om 3D vector of frequencies of oscillator Om=(Om_x,Om_y,Om_z)
        */
        Actions_HarmonicOscillator(VecDoub Om): HarmonicOscillator(Om){};
        //! Finds actions
        /*!
          \param x phase-space point (x,v)
          \param params -- does nothing
        */
        VecDoub actions(const VecDoub &x, void *params=nullptr);
        //! Finds angles
        /*!
          \param x phase-space point (x,v)
          \param params -- does nothing
        */
        VecDoub angles(const VecDoub &x, void *params=nullptr);
};

/*! Action finder for isochrone potential
Finds actions in spherical potential Phi(r)=-GM/(b+sqrt(b^2+r^2))
*/
class Actions_Isochrone : public Action_Finder, public Isochrone{
    private:
    public:
        //! Actions_Isochrone constructor.
        /*!
          \param GM 'mass' parameter of potential
          \param b scale parameter of potential
        */
        Actions_Isochrone(double GM, double b): Isochrone(GM,b){};
        //! Finds actions
        /*!
          \param x phase-space point (x,v)
          \param params -- does nothing
        */
        VecDoub actions(const VecDoub &x, void *params=nullptr);
        //! Finds angles
        /*!
          \param x phase-space point (x,v)
          \param params -- does nothing
        */
        VecDoub angles(const VecDoub &x, void *params=nullptr);
        //! Finds frequencies
        /*!
          \param x phase-space point (x,v)
        */
        VecDoub freq(const VecDoub &x);
        //! Finds Hessian
        /*!
          \param x phase-space point (x,v)
        */
        VecDoub Hessian(const VecDoub &x);
};

#endif
// ============================================================================
