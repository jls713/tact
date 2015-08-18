#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface/GSLInterface.h"
#include "gnuplot/gnuplot_i.h"
#include <gsl/gsl_poly.h>
#include "falPot.h"
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
#include "it_torus.h"
#include "stackel_fit.h"

int main(int argc, char*argv[]){
	// GalPot Pot("../../Torus/pot/PJM11.Tpot");
	Logarithmic Pot(220.,1.,0.9);

	if(argc<8){
		std::cerr<<"Need to pass phase-space point and filename\n";
		return 0;
	}
	VecDoub X(6,0.);
	for(unsigned i=0;i<6;++i)
		X[i]=atof(argv[i+1]);
	Orbit O(&Pot);
	// Fudge
	Actions_AxisymmetricStackel_Fudge AA(&Pot,-30.);

	// Iterative Torus
	// IterativeTorusMachine Tor(&AA,&Pot,1e-8,5,1e-3);

	// Generating Function
	Actions_Genfunc AG(&Pot,"axisymmetric");

	// Average generating Function
	Actions_Genfunc_Average AGav(&Pot,"axisymmetric");

	// uvorb
	uv_orb UV(&Pot,1.,20.,10,10,"example.delta_uv");

	// Polar Adiabatic
	Actions_PolarAdiabaticApproximation PAA(&Pot,"example.paa",true,false);

	// Spheroidal Adiabatic
	Actions_SpheroidalAdiabaticApproximation SAA(&Pot,"example.saa",true,false,-30.);

	// Spheroidal Adiabatic
	Actions_StackelFit SF(&Pot);

	double tt = 10.;
	if(argc>8) tt=atof(argv[8]);
	O.integrate(X,tt*Pot.torb(X),0.01*Pot.torb(X));
	// O.plot(0,2);

	std::ofstream outfile;
	outfile.open(argv[7]);
	outfile<<"# Fudge Genfunc GenfuncAv uvOrb PAA SAA FIT\n";

	int guess_alpha=1;
	VecDoub Fudge, ITorus, Genfunc, GenfuncAv, uvAct, paaAct, saaAct, fitAct;
	for(auto i:O.results()){
		Fudge = AA.actions(i,&guess_alpha);
		// ITorus = Tor.actions(i);
		Genfunc = AG.actions(i);
		GenfuncAv = AGav.actions(i);
		uvAct = UV.actions(i);
		paaAct = PAA.actions(i);
		saaAct = SAA.actions(i,&guess_alpha);
		fitAct = SF.actions(i);
		outfile
				<<Fudge[0]<<" "<<Fudge[2]<<" "
				// <<ITorus[0]<<" "<<ITorus[2]<<" "
				<<Genfunc[0]<<" "<<Genfunc[2]<<" "
				<<GenfuncAv[0]<<" "<<GenfuncAv[2]<<" "
				<<uvAct[0]<<" "<<uvAct[2]<<" "
				<<paaAct[0]<<" "<<paaAct[2]<<" "
				<<saaAct[0]<<" "<<saaAct[2]<<" "
				<<fitAct[0]<<" "<<fitAct[2]<<" ";
		for(auto j:i) outfile<<j<<" ";
		outfile <<std::endl;
	}
	outfile.close();
}
