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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <ctime>
#include <ratio>
#include <chrono>

MatDoub dvdJ(VecDoub X, double dv, Actions_Genfunc *AG){
	VecDoub Y = X, uJ,dJ;MatDoub mat(3,VecDoub(3,0.));
	for(int j=0;j<3;++j){
		Y[j+3]+=   dv; uJ = AG->actions(Y);
		Y[j+3]-=2.*dv; dJ = AG->actions(Y);
		for(unsigned i=0;i<3;++i) mat[i][j]=(uJ[i]-dJ[i])/(2.*dv);
		Y[j+3]+=dv;
	}
	MatDoub inv = inverse3D(mat);
	return inv;
}

using namespace std::chrono;

int main(int argc, char*argv[]){

	GalPot Pot("../Torus/pot/PJM11_best.Tpot");
	// Logarithmic Pot(220.,1.,0.9);
	VecDoub X(6,1e-4);

	X[0]=conv::StandardSolarPAUL[0];X[2]=conv::StandardSolarPAUL[1];
	X[4]=sqrt(X[0]*-Pot.Forces(X)[0]);
	printVector(X);
	Orbit O(&Pot);
	// Fudge
	Actions_AxisymmetricStackel_Fudge AA(&Pot,-30.);

	// Iterative Torus
	IterativeTorusMachine Tor(&AA,&Pot,1e-8,5,1e-3);

	// Generating Function
	Actions_Genfunc AG(&Pot,"axisymmetric");

	// Average generating Function
	Actions_Genfunc_Average AGav(&Pot,"axisymmetric");

	// uvorb
	uv_orb UV(&Pot,1.,40.,10,10,"example.delta_uv");

	// Polar Adiabatic
	Actions_PolarAdiabaticApproximation PAA(&Pot,"example.paa",true,false);

	// Spheroidal Adiabatic
	Actions_SpheroidalAdiabaticApproximation SAA(&Pot,"example.saa",true,false,-30.);

	// Spheroidal Adiabatic
	Actions_StackelFit SF(&Pot);

	std::ofstream outfile;
	outfile.open(argv[1]);
	outfile<<"# JRJzLz Fudge ItTorus Genfunc GenfuncAv uvOrb PAA SAA FIT\n";
	double VMax = sqrt(2.*(Pot.Phi({40.,0.,40.})-Pot.Phi(X)));

	VecDoub range = create_range(0.01*VMax,0.3*VMax,30);

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	for(auto j: range){
		X[3]=sqrt(2.)*j;
		X[5]=sqrt(2.)*j;
		O.integrate(X,10.*Pot.torb(X),0.2*Pot.torb(X));

		int guess_alpha=1;
		MatDoub FResults,ITResults,GResults,GAvResults,UVResults,PAAResults,SAAResults,FITResults;
		VecDoub Fudge, ITorus, Genfunc, GenfuncAv, uvAct, paaAct, saaAct, fitAct;
		double FudgeV=0., ITorusV=0., GenfuncV=0., GenfuncAvV=0., uvActV=0., paaActV=0., saaActV=0., fitActV=0.;
		MatDoub dvdJ_e;
		t1 = high_resolution_clock::now();
		std::vector<nanoseconds> times(8,duration_cast<nanoseconds>(t1-t1));
		for(auto i:O.results()){
			t1 = high_resolution_clock::now();
			Genfunc = AG.actions(i);
			times[2]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);GenfuncAv.resize(3);
			Genfunc.resize(3);
			GResults.push_back({Genfunc[0],Genfunc[1],Genfunc[2]});
		}
		VecDoub acts = columnMean(GResults);
		VecDoub GF_SD = columnSD(GResults);
		outfile<<acts[0]<<" "<<acts[2]<<" "<<(acts[0]+acts[2])/fabs(acts[1])<<" ";
		int N=0;
		for(auto i:O.results()){
			// t1 = high_resolution_clock::now();
			// Fudge = AA.actions(i,&guess_alpha);
			// times[0]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			// Fudge.resize(3);Fudge=Fudge-acts;
			// t1 = high_resolution_clock::now();
			// ITorus = Tor.actions(i);ITorus.resize(3);ITorus=ITorus-acts;
			// times[1]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			// t1 = high_resolution_clock::now();
			// GenfuncAv = AGav.actions(i);
			// times[3]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			// GenfuncAv.resize(3);GenfuncAv=GenfuncAv-acts;
			t1 = high_resolution_clock::now();
			uvAct = UV.actions(i);
			times[4]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			uvAct.resize(3);uvAct=uvAct; //-acts;
			// t1 = high_resolution_clock::now();
			// paaAct = PAA.actions(i);
			// times[5]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			// paaAct.resize(3);paaAct=paaAct-acts;
			// t1 = high_resolution_clock::now();
			// saaAct = SAA.actions(i,&guess_alpha);
			// times[6]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			// saaAct.resize(3);saaAct=saaAct-acts;
			// t1 = high_resolution_clock::now();
			// fitAct = SF.actions(i);
			// times[7]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			// fitAct.resize(3);fitAct=fitAct-acts;

			// dvdJ_e=dvdJ(i,.5,&AG);
			// FudgeV += normsq(dvdJ_e*Fudge);
			// ITorusV += normsq(dvdJ_e*ITorus);
			// GenfuncV += normsq(dvdJ_e*GResults[N]);
			// GenfuncAvV += normsq(dvdJ_e*GenfuncAv);
			// uvActV += normsq(dvdJ_e*uvAct);
			// paaActV += normsq(dvdJ_e*paaAct);
			// saaActV += normsq(dvdJ_e*saaAct);
			// fitActV += normsq(dvdJ_e*fitAct);

			// FResults.push_back({Fudge[0],Fudge[2]});
			// ITResults.push_back({ITorus[0],ITorus[2]});
			// GAvResults.push_back({GenfuncAv[0],GenfuncAv[2]});
			UVResults.push_back({uvAct[0],uvAct[2]});
			// PAAResults.push_back({paaAct[0],paaAct[2]});
			// SAAResults.push_back({saaAct[0],saaAct[2]});
			// FITResults.push_back({fitAct[0],fitAct[2]});
			++N;
		}
		// for(auto k:columnSD(FResults)) outfile<<k<<" "; outfile<<sqrt(FudgeV/(double)N)<<" "<<times[0].count()<<" ";
		// for(auto k:columnSD(ITResults)) outfile<<k<<" ";outfile<<sqrt(ITorusV/(double)N)<<" ";
		// for(auto k:columnSD(GResults)) outfile<<k<<" ";outfile<<sqrt(GenfuncV/(double)N)<<" "<<times[2].count()<<" ";
		// for(auto k:columnSD(GAvResults)) outfile<<k<<" ";outfile<<sqrt(GenfuncAvV/(double)N)<<" "<<times[3].count()<<" ";
		for(auto k:columnSD(UVResults)) outfile<<k<<" ";outfile<<sqrt(uvActV/(double)N)<<" "<<times[4].count()<<" ";
		// for(auto k:columnSD(PAAResults)) outfile<<k<<" ";outfile<<sqrt(paaActV/(double)N)<<" "<<times[5].count()<<" ";
		// for(auto k:columnSD(SAAResults)) outfile<<k<<" ";outfile<<sqrt(saaActV/(double)N)<<" "<<times[6].count()<<" ";
		// for(auto k:columnSD(FITResults)) outfile<<k<<" ";outfile<<sqrt(fitActV/(double)N)<<" "<<times[7].count()<<" ";
		outfile<<std::endl;
	}
	outfile.close();
}
