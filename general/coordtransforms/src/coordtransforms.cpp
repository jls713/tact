// ============================================================================
/// \file coordtransform/src/coordtransforms.cpp
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
/// \brief Coordinate Transformations
///
/// 1. Contains definitions of constants and conventions
/// 2. Contains Transformations between different coordinate systems
///
//=============================================================================

/*==============================================*/
/* 			Coordinate Transformations 			*/
/*==============================================*/

// ============================================================================
#include <vector>
#include <cmath>
#include <iostream>
#include "utils.h"
#include "coordtransforms.h"
// ============================================================================

namespace conv{
VecDoub CartesianToSphericalPolar(const VecDoub& Cartesian){
	double r = sqrt(Cartesian[0]*Cartesian[0]+Cartesian[1]*Cartesian[1]+Cartesian[2]*Cartesian[2]);
    VecDoub SPolar = {r,atan2(Cartesian[1],Cartesian[0]),acos(Cartesian[2]/r)};
    if(Cartesian.size()==3)	return SPolar;
    SPolar.push_back((Cartesian[3]*cos(SPolar[1])+Cartesian[4]*sin(SPolar[1]))*sin(SPolar[2])+cos(SPolar[2])*Cartesian[5]);
    SPolar.push_back(-Cartesian[3]*sin(SPolar[1])+Cartesian[4]*cos(SPolar[1]));
    SPolar.push_back((Cartesian[3]*cos(SPolar[1])+Cartesian[4]*sin(SPolar[1]))*cos(SPolar[2])-sin(SPolar[2])*Cartesian[5]);
    return SPolar;
}

VecDoub SphericalPolarToCartesian(const VecDoub& Spherical){
    VecDoub SPolar = {Spherical[0]*sin(Spherical[2])*cos(Spherical[1]),
    					Spherical[0]*sin(Spherical[2])*sin(Spherical[1]),
    					Spherical[0]*cos(Spherical[2])};
    if(Spherical.size()==3)	return SPolar;
    return SPolar;
}

// ======================================================================================
// Cartesian <==> Polar
VecDoub CartesianToPolar(const VecDoub& Cartesian){
	// X,Y,Z -> R,phi,z
	VecDoub Polar {	sqrt(Cartesian[0]*Cartesian[0]+Cartesian[1]*Cartesian[1]),
					atan2(Cartesian[1],Cartesian[0]),
					Cartesian[2]};
	if(Cartesian.size()==3)	return Polar;
	// vx,vy,vz -> vR,vphi,vz
	else{
		double cp = cos(Polar[1]), sp = sin(Polar[1]);
		VecDoub PolarVel {	Cartesian[3]*cp+Cartesian[4]*sp,Cartesian[4]*cp-Cartesian[3]*sp,
					        Cartesian[5]};
		for (	VecDoub::iterator it = PolarVel.begin();
				it != PolarVel.end(); ++it) Polar.push_back(*it);
		return Polar;
		}
}

VecDoub PolarToCartesian(const VecDoub& Polar){
	// R,phi,z -> X,Y,Z
	double cp = cos(Polar[1]), sp = sin(Polar[1]);
	VecDoub Cartesian {	Polar[0]*cp,
						Polar[0]*sp,
						Polar[2]};
	if(Polar.size()==3) return Cartesian;
	// vR,vphi,vz -> vx,vy,vz
	else{
		VecDoub CartVel {Polar[3]*cp-Polar[4]*sp,Polar[4]*cp+Polar[3]*sp,Polar[5]};
		for (	VecDoub::iterator it = CartVel.begin();
				it != CartVel.end(); ++it) Cartesian.push_back(*it);
		return Cartesian;
		}
}

// ======================================================================================
// Galactic <==> Cartesian
VecDoub GalacticToCartesian(const VecDoub &Galactic,
								      const VecDoub& SolarPosition){
	// l,b,s->X,Y,Z
	double 	cl = cos(Galactic[0]), sl = sin(Galactic[0]),
			cb = cos(Galactic[1]), sb = sin(Galactic[1]);

	double x = Galactic[2]*cb*cl;
	double z = Galactic[2]*sb;
	// Need to rotate to account for the height of the Sun above the plane
	double h = sqrt(SolarPosition[0]*SolarPosition[0]
	                +SolarPosition[1]*SolarPosition[1]);
	double ct = SolarPosition[0]/h, st = SolarPosition[1]/h;

	VecDoub Cartesian {	SolarPosition[0]-ct*x-st*z,
						-Galactic[2]*cb*sl,
						-st*x+ct*z+SolarPosition[1]};
	if(Galactic.size()==3)return Cartesian;
	// vlos,mu_lcos(b),mu_b -> vx,vy,vz
	// in units km/s, mas/yr -> km/s
	else{
		double vl = PM_Const*Galactic[2]*Galactic[4];
		double vb = PM_Const*Galactic[2]*Galactic[5];
		double tmp = cb*Galactic[3]-sb*vb;

		double vx = cl*tmp-sl*vl+SolarPosition[2];
		double vy = sl*tmp+cl*vl+SolarPosition[3];
		double vz = sb*Galactic[3]+cb*vb+SolarPosition[4];
		VecDoub CartVel{-(vx*ct+vz*st),-vy,-vx*st+vz*ct};
	  	for (	VecDoub::iterator it = CartVel.begin();
				it != CartVel.end(); ++it) Cartesian.push_back(*it);
			return Cartesian;
	}
}

VecDoub GalacticToCartesian(const VecDoub &Galactic){
  return GalacticToCartesian(Galactic,StandardSolar);
}

VecDoub CartesianToGalactic(const VecDoub &Cartesian,
									const VecDoub& SolarPosition){
	// X,Y,Z->l,b,s
	double tmp1 = SolarPosition[0]-Cartesian[0];
	double tmp2 = -Cartesian[1];
	double tmp3 = Cartesian[2]-SolarPosition[1];
	// Need to rotate to account for the height of the Sun above the plane
	double h = sqrt(SolarPosition[0]*SolarPosition[0]
	                +SolarPosition[1]*SolarPosition[1]);
	double ct = SolarPosition[0]/h, st = SolarPosition[1]/h;

	double x = tmp1*ct-tmp3*st, z = tmp1*st+tmp3*ct;

	double Distance = norm<double>({x,tmp2,z});

	VecDoub Galactic {	atan2(tmp2,x),
						asin(z/Distance),
						Distance};
	if(Cartesian.size()==3)return Galactic;
	// vx,vy,vz -> vlos,mu_lcos(b),mu_b
	// in units km/s -> km/s mas/yr
	else{ 	double vx=-Cartesian[3]*ct-Cartesian[5]*st-SolarPosition[2];
			double vy = -Cartesian[4]-SolarPosition[3];
			double vz = Cartesian[5]*ct+Cartesian[3]*st-SolarPosition[4];
			double 	cl = cos(Galactic[0]), sl = sin(Galactic[0]),
			cb = cos(Galactic[1]), sb = sin(Galactic[1]);
			VecDoub GalVel {vx*cl*cb+vy*sl*cb+vz*sb,(-vx*sl+vy*cl)/(PM_Const*Distance),
				        	(-vx*cl*sb-vy*sl*sb+vz*cb)/(PM_Const*Distance)};
			for (	VecDoub::iterator it = GalVel.begin();
					it != GalVel.end(); ++it) Galactic.push_back(*it);
			return Galactic;
		}
}

VecDoub CartesianToGalactic(const VecDoub &Cartesian){
  return CartesianToGalactic(Cartesian,StandardSolar);
}
// ======================================================================================
// Galactic <==> Polar
VecDoub PolarToGalactic(const VecDoub &Polar,
									const VecDoub& SolarPosition){
	// R,phi,z -> l,b,s; SolarPosition=(R0,z0)
	VecDoub Cart = PolarToCartesian(Polar);
	VecDoub Galactic = CartesianToGalactic(Cart,SolarPosition);
	return Galactic;
	}

VecDoub PolarToGalactic(const VecDoub &Polar){
  return PolarToGalactic(Polar,StandardSolar);
}

VecDoub GalacticToPolar(const VecDoub &Galactic,
									const VecDoub& SolarPosition){
	// l,b,s->R,phi,z; SolarPosition=(R0,z0)
	VecDoub Cart = GalacticToCartesian(Galactic,SolarPosition);
	VecDoub Polar = CartesianToPolar(Cart);
	return Polar;
	}

VecDoub GalacticToPolar(const VecDoub &Galactic){
  return GalacticToPolar(Galactic,StandardSolar);
}
// ======================================================================================
// Equatorial <==> Galactic

VecDoub EquatorialToGalactic(const VecDoub &Equatorial){
	//alpha, dec, s => l,b,s
	double alpha = Equatorial[0], delta = Equatorial[1];
	double cd = cos(delta), sd = sin(delta);
	double dalpha = alpha-RA_GP;
	double b=asin(sdGP*sd+cdGP*cd*cos(dalpha));
	double l=lCP-atan2(cd*sin(dalpha),cdGP*sd-sdGP*cd*cos(dalpha));
	if(l<0.)l+=2.*PI;
	VecDoub Galactic {l,b,Equatorial[2]};
	if(Equatorial.size()==3)return Galactic;
	else{
		//vlos, ma_cos(d), md => vlos, ml_cos(b), mb
		double cb = cos(b), sb = sin(b);
		double dl = lCP-l;
		double A11=(sdGP*cd-cdGP*sd*cos(dalpha))/cb;
		double A12=-cdGP*sin(dalpha)/cb;
		double A21,A22;
		if(fabs(cos(dl))>fabs(sin(dl))){
			A21= (sd*sin(dalpha)-sb*sin(dl)*A11)/cos(dl);
			A22=-(   cos(dalpha)+sb*sin(dl)*A12)/cos(dl);
		}else{
			A21=(cdGP*cd+sdGP*sd*cos(dalpha)+sb*cos(dl)*A11)/sin(dl);
			A22=(sdGP*sin(dalpha)+sb*cos(dl)*A12)/sin(dl);
		}

		VecDoub GalVel {Equatorial[3],A21*Equatorial[5]+A22*Equatorial[4],
						A11*Equatorial[5]+A12*Equatorial[4]};
		for (	VecDoub::iterator it = GalVel.begin();
				it != GalVel.end(); ++it) Galactic.push_back(*it);
		return Galactic;
		}
}
std::vector<VecDoub> EquatorialToGalacticwithErrors(const VecDoub &Equatorial,
	const VecDoub &Errors){
	//alpha, dec, s => l,b,s
	double alpha = Equatorial[0], delta = Equatorial[1];
	double cd = cos(delta), sd = sin(delta);
	double dalpha = alpha-RA_GP;
	double b=asin(sdGP*sd+cdGP*cd*cos(dalpha));
	double l=lCP-atan2(cd*sin(alpha-RA_GP),cdGP*sd-sdGP*cd*cos(dalpha));
	if(l<0.)l+=2.*PI;
	if(Equatorial.size()==3){	std::vector<VecDoub> Galactic {{l,b,Equatorial[2]},Errors};
								return Galactic;}
	else{
		//vlos, ma_cos(d), md => vlos, ml_cos(b), mb
		double cb = cos(b), sb = sin(b);
		double A11=(sdGP*cd-cdGP*sd*cos(dalpha))/cb;
		double A12=-cdGP*sin(dalpha)/cb;
		double A21,A22;
		double dl = lCP-l;
		if(fabs(cos(dl))>fabs(sin(dl))){
			A21=(sd*sin(dalpha)-sb*sin(dl)*A11)/cos(dl);
			A22=-(cos(dalpha)+sb*sin(dl)*A12)/cos(dl);
		}else{
			A21=(cdGP*cd+sdGP*sd*cos(dalpha)+sb*cos(dl)*A11)/sin(dl);
			A22=(sdGP*sin(dalpha)+sb*cos(dl)*A12)/sin(dl);
		}

		std::vector<VecDoub> Galactic {
		{l,b,Equatorial[2],Equatorial[3],
		A21*Equatorial[5]+A22*Equatorial[4],A11*Equatorial[5]+A12*Equatorial[4]}
		,{Errors[0],Errors[1],Errors[2],Errors[3],
		sqrt(A21*A21*Errors[5]*Errors[5]+A22*A22*Errors[4]*Errors[4]),
		sqrt(A11*A11*Errors[5]*Errors[5]+A12*A12*Errors[4]*Errors[4])}};
		return Galactic;
		}
}
VecDoub GalacticToEquatorial(const VecDoub &Galactic){
	//l,b,s => alpha, dec, s
	double l = Galactic[0], b = Galactic[1];
	double cb = cos(b),sb = sin(b);
	double dl = lCP-l;
	double delta=asin(cdGP*cb*cos(-dl)+sb*sdGP);
	double alpha=RA_GP+atan2(cb*sin(dl),sb*cdGP-cb*sdGP*cos(-dl));
	if(alpha>2.*PI)alpha-=2.*PI;
	VecDoub Equatorial {alpha,delta,Galactic[2]};
	if(Galactic.size()==3)return Equatorial;
	else{
		double dalpha = alpha-RA_GP;
		//vlos, ml_cos(b), mb => vlos, ma_cos(d), md
		double cd = cos(delta), sd = sin(delta);
		double A11=(sdGP*cd-cdGP*sd*cos(dalpha))/cb;
		double A12=-cdGP*sin(dalpha)/cb;
		double A21,A22;
		if(fabs(cos(dl))>fabs(sin(dl))){
			A21=(sd*sin(dalpha)-sb*sin(dl)*A11)/cos(dl);
			A22=-(cos(dalpha)+sb*sin(dl)*A12)/cos(dl);
		}else{
			A21=(cdGP*cd+sdGP*sd*cos(dalpha)+sb*cos(dl)*A11)/sin(dl);
			A22=(sdGP*sin(dalpha)+sb*cos(dl)*A12)/sin(dl);
		}
		double Prod = A11*A22-A12*A21;
		VecDoub EqVel {Galactic[3],(A11*Galactic[4]-A21*Galactic[5])/Prod,
						(A22*Galactic[5]-A12*Galactic[4])/Prod};
		for (	VecDoub::iterator it = EqVel.begin();
				it != EqVel.end(); ++it) Equatorial.push_back(*it);
		return Equatorial;
		}
}
}
// ============================================================================
