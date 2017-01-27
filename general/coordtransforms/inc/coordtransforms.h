// ============================================================================
/// \file coordtransform/inc/coordtransforms.h
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
/*! \brief Coordinate Transformations
 *
 *  1. Contains definitions of constants and conventions
 *  2. Contains Transformations between different coordinate systems
 */
//=============================================================================

#ifndef COORDTRANSFORMS_H
#define COORDTRANSFORMS_H

// ============================================================================
#include "utils.h"
// ============================================================================

namespace conv{

#define NATURAL_UNITS 1

// #ifndef KPC_KMS_MSUN_UNITS
// #define KPC_KMS_MSUN_UNITS 0
// #endif

// #ifndef KPC_KMS_11MSUN_UNITS
// #define KPC_KMS_11MSUN_UNITS 0
// #endif

// ============================================================================
// Conventions
// ============================================================================
const VecDoub SunPeculiar {11.1,12.24,7.25}; // peculiar motion of Sun
const VecDoub StandardSolar {8.0,0.014,11.1,232.24,7.25}; // in km/s with v_c = 220 km/s
const VecDoub StandardSolar2 {8.29,0.014,0.01135,0.257055,0.0074148}; // in kpc/Myr with v_c=239.1km/s
const VecDoub StandardSolarPAUL {8.29,0.014,11.1,251.34,7.25}; // in km/s with v_c = 239.1 km/s
const double G =
#ifdef KPC_KMS_11MSUN_UNITS
	430091.5694; /* in units kpc, km/s and 10^11 M_solar */
#elif KPC_KMS_MSUN_UNITS
	4.300918e-6; // In Units kpc (km/s) / M_sun
#elif NATURAL_UNITS
	1.;
#endif
const double FPG = 4.*PI*G; // In Units kpc (km/s) / M_sun
const double masyr2radMyr = 4.8481368111e-3;
const double deg2rad = 0.017453292;
const double kpcMyr2kms = 977.775;
const double kpcMyr2kmsSq = kpcMyr2kms*kpcMyr2kms;
const double kms2kpcMyr = 1./kpcMyr2kms;
const double kms2kpcGyr = 1000.*kms2kpcMyr;
const double RA_GP=3.36603292,decGP=0.473477282,lCP=2.145566725;
const double PM_Const = 4.74057170372;

const double sdGP = sin(decGP), cdGP = cos(decGP);

// ============================================================================
// Cartesian = (x,y,z,vx,vy,vz):
// 		 x away from centre of Galaxy towards Sun
//		 y opposite to Galactic rotation
//		 z up towards North Galactic Pole
// Polar = (R,\phi,z,vR,v\phi,vz)
// Galactic = (l,b,s,vlos,\mu_l\cos(b),\mu_b)
// Equatorial =(\alpha,\delta,s,vlos,\mu_\alpha\cos(\delta),mu_\delta)
// Units: velocities in km/s, positions kpc and proper motions in mas/yr
// ============================================================================
/*! Cartesian ==> SphericalPolar */
VecDoub CartesianToSphericalPolar(const VecDoub& Cartesian);
/*! Cartesian <== SphericalPolar */
VecDoub SphericalPolarToCartesian(const VecDoub& Polar);
// ============================================================================
/*! Cartesian ==> Polar */
VecDoub CartesianToPolar(const VecDoub& Cartesian);
/*! Cartesian <== Polar */
VecDoub PolarToCartesian(const VecDoub& Polar);

// ============================================================================
/*! Galactic ==> Cartesian */
VecDoub GalacticToCartesian(const VecDoub &Galactic,
								      const VecDoub& SolarPosition=StandardSolar);
/*! Galactic <== Cartesian */
VecDoub CartesianToGalactic(const VecDoub &Cartesian,
									const VecDoub& SolarPosition=StandardSolar);

// ============================================================================
/*! Galactic ==> Polar */
VecDoub PolarToGalactic(const VecDoub &Polar,
									const VecDoub& SolarPosition=StandardSolar);
/*! Galactic <== Polar */
VecDoub GalacticToPolar(const VecDoub &Galactic,
									const VecDoub& SolarPosition=StandardSolar);
// ============================================================================
/*! Equatorial ==> Galactic */
VecDoub EquatorialToGalactic(const VecDoub &Eq);
/*! Equatorial ==> Galactic including errors */
std::vector<VecDoub> EquatorialToGalacticwithErrors(const VecDoub &E,const VecDoub &EqE);
/*! Equatorial <== Galactic */
VecDoub GalacticToEquatorial(const VecDoub &Galactic);
}
#endif
// ============================================================================
