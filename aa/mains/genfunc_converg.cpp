// ============================================================================
/// \file src/genfunc_converg.cpp
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
/// \brief Computes the errors in the actions as a function of the parameters in the O2GF method
///
/// Must pass phase-space point x y z vx vy vz, output file and no of orbital
/// times to integrate for.
/// e.g. mains/./genfunc_converg.exe 8.29 0.1 0.1 30.22 211.1 19.22 thin 10.
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
#include "genfunc_aa.h"
#include <ctime>
#include <ratio>
#include <chrono>

#ifdef TORUS
#include "falPot.h"
#include "it_torus.h"
#include "PJM_cline.h"
#endif

using namespace std::chrono;
int main(int argc, char*argv[]){

	#ifdef TORUS
	GalPot Pot("pot/Piffl14.Tpot");
	WrapperTorusPotential TPot(&Pot);
	std::cout<<TPot.KapNuOm(8.29)*conv::kpcMyr2kms<<std::endl;
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

	// Generating Function
	Actions_Genfunc AG(&Pot,"axisymmetric");

	double tt = 10.;
	double tstep = 0.01*Pot.torb(X);
	if(argc>7) tt=atof(argv[8]);
	O.integrate(X,tt*Pot.torb(X),tstep);

	VecDoub range={2,4,6,8,10,12,14,16,18};

	std::ofstream outfile;
	outfile.open(argv[7]);
	outfile<<"# DeltaJr DeltaJz Time\n";

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	auto timeT=duration_cast<nanoseconds>(t1-t1);
	VecDoub Genfunc;int N=O.results().size();
	std::vector<VecDoub> R = std::vector<VecDoub>(N,VecDoub(2,0.));
	for(auto j: range){
		int n=0;
		t1 = high_resolution_clock::now();
		for(auto i:O.results()){
			Genfunc = AG.full_actions(i,24,2400,j);
			R[n][0]=Genfunc[0];
			R[n][1]=Genfunc[2];
			++n;
		}
		timeT=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
		VecDoub CSD=columnSD(R);
		outfile<<j<<" "<<CSD[0]<<" "<<CSD[1]<<" "<<timeT.count()/N<<std::endl;
	}
	outfile.close();
}

// ============================================================================
