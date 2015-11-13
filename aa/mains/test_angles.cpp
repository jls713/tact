// ============================================================================
/// \file src/test_angles.cpp
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
/// \brief Computes the angles and frequencies using most of the available methods for axisymmetric potentials
///
/// Must pass phase-space point x y z vx vy vz, output file and no of orbital
/// times to integrate for.
/// e.g. mains/./test_angles.exe 8.29 0.1 0.1 30.22 211.1 19.22 thin 10.
// ============================================================================
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface/GSLInterface.h"
#include "gnuplot/gnuplot_i.h"
#include <gsl/gsl_poly.h>
#include "utils.h"
#include "coordsys.h"
#include "coordtransforms.h"
#include "potential.h"
#include "orbit.h"
#include "aa.h"
#include "stackel_aa.h"
#include "spherical_aa.h"
#include "genfunc_aa.h"
#include "adiabatic_aa.h"
#include "uv_orb.h"
#include "lmn_orb.h"
#include "stackel_fit.h"
#ifdef TORUS
#include "falPot.h"
#include "it_torus.h"
#include "PJM_cline.h"
#endif
int main(int argc, char*argv[]){

	#ifdef TORUS
	GalPot Pot("pot/Piffl14.Tpot");
	WrapperTorusPotential TPot(&Pot);
	std::cerr<<"Using Piffl14 GalPot"<<std::endl;
	#else
	Logarithmic Pot(220.,1.,0.9);
	std::cerr<<"Using logarithmic potential -- to use GalPot need TORUS flag to make"<<std::endl;
	#endif

	if(argc<8)
		std::cerr<<"Need to pass phase-space point and filename\n";
	VecDoub X(6,0.);
	for(unsigned i=0;i<6;++i)
		X[i]=atof(argv[i+1]);
	Orbit O(&Pot);
	// Fudge
	Actions_AxisymmetricStackel_Fudge AA(&Pot,100.);

	// Iterative Torus
	#ifdef TORUS
	IterativeTorusMachine Tor(&AA,&Pot,1e-8,5,1e-3);
	#endif
	// Generating Function
	Actions_Genfunc AG(&Pot,"axisymmetric");
	VecDoub acts = AG.actions(X);
	double Rmin=1.,Rmax=40.,zmax=30.;
	#ifdef TORUS
	Actions J;J[0]=acts[0]/conv::kpcMyr2kms;
	J[2]=acts[1]/conv::kpcMyr2kms;J[1]=acts[2]/conv::kpcMyr2kms;
	Torus T; T.AutoFit(J,&TPot,1e-5);
	Rmin = .9*T.minR();
	Rmax = 1.1*T.maxR();
	zmax = 1.1*T.maxz();
	#endif

	// Average generating Function
	Actions_Genfunc_Average AGav(&Pot,"axisymmetric");

	// uvorb
	uv_orb UV(&Pot,3.,30.,10,10,"example.delta_uv");

	// Cylindrical Adiabatic
	Actions_CylindricalAdiabaticApproximation PAA(&Pot,"example.paa",true,false,4.,30.,15.,60);

	std::cout<<"Action methods loaded"<<std::endl;
	// Spheroidal Adiabatic
	Actions_SpheroidalAdiabaticApproximation SAA(&Pot,"example.saa",true,false,100.,4.,30.,15.);

	// Spheroidal Adiabatic
	Actions_StackelFit SF(&Pot);


	double tt = 10.;
	double tstep = 0.01*Pot.torb(X);
	if(argc>8) tt=atof(argv[8]);
	O.integrate(X,tt*Pot.torb(X),tstep);
	//O.plot(0,1);

	std::ofstream outfile;
	outfile.open(argv[7]);
	outfile<<"# Fudgev1 ";
	#ifdef TORUS
	outfile<<"ItTC ";
	#endif
	outfile<<"O2GF AvGF Fudgev2 CAA SAA Fit\n";

	int guess_alpha=1; int N=0;
	VecDoub Fudge, ITorus, Genfunc, GenfuncAv, uvAct, paaAct, saaAct, fitAct;
	for(auto i:O.results()){
		Fudge = AA.angles(i,&guess_alpha);
		#ifdef TORUS
		ITorus = Tor.angles(i);
		#endif
		Genfunc = AG.angles(i);
		GenfuncAv = AGav.angles(i);
		uvAct = UV.angles(i);
		paaAct = PAA.angles(i);
		saaAct = SAA.angles(i,&guess_alpha);
		fitAct = SF.angles(i);
		outfile << tstep*N<<" ";
		for(auto k: Fudge) outfile<<k<<" ";
		#ifdef TORUS
		for(auto k: ITorus) outfile<<k<<" ";
		#endif
		for(auto k: Genfunc) outfile<<k<<" ";
		for(auto k: GenfuncAv) outfile<<k<<" ";
		for(auto k: uvAct) outfile<<k<<" ";
		for(auto k: paaAct) outfile<<k<<" ";
		for(auto k: saaAct) outfile<<k<<" ";
		for(auto k: fitAct) outfile<<k<<" ";
		for(auto j:i) outfile<<j<<" ";
		outfile <<std::endl;
	++N;
	}
	outfile.close();
}

// ============================================================================
