// ============================================================================
/// \file src/many_tori.cpp
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
/// \brief Computes the variance of the actions, angles and frequencies using various methods for axisymmetric potential and computes total time taken.
///
/// Must pass output file
/// e.g. mains/./many_tori.exe output
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <ctime>
#include <ratio>
#include <chrono>

#ifdef TORUS
#include "it_torus.h"
#include "falPot.h"
#endif

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

	#ifdef TORUS
	GalPot Pot("pot/Piffl14.Tpot");
	WrapperTorusPotential TPot(&Pot);
	// GalPot Pot("../Torus/pot/PJM11.Tpot");
	std::cout<<TPot.KapNuOm(8.29)*conv::kpcMyr2kms<<std::endl;
	#else
	Logarithmic Pot(220.,1.,0.9);
	#endif

	VecDoub X(6,1e-4);

	X[0]=conv::StandardSolarPAUL[0];X[2]=conv::StandardSolarPAUL[1];
	X[4]=sqrt(X[0]*-Pot.Forces(X)[0]);
	printVector(X);
	Orbit O(&Pot,1e-8);
	// Fudge
	Actions_AxisymmetricStackel_Fudge AA(&Pot,-30.);

	// Iterative Torus
	#ifdef TORUS
	IterativeTorusMachine Tor(&AA,&Pot,1e-4,5,5e-3);
	#endif

	// Generating Function
	Actions_Genfunc AG(&Pot,"axisymmetric");

	// Average generating Function
	Actions_Genfunc_Average AGav(&Pot,"axisymmetric");

	// uvorb
	uv_orb UV(&Pot,4.,30.,50,50,"example.delta_uv");

	// Polar Adiabatic
	Actions_PolarAdiabaticApproximation PAA(&Pot,"example.paa",true,false,4.,30.,15.);

	// Spheroidal Adiabatic
	Actions_SpheroidalAdiabaticApproximation SAA(&Pot,"example.saa",true,false,100.,4.,30.,15.);

	// Spheroidal Adiabatic
	Actions_StackelFit SF(&Pot,1e-5);

	std::ofstream outfile;
	outfile.open(argv[1]);
	outfile<<"# JR Lz Jz JRJzLz ";
	#ifdef TORUS
	outfile<<"Rperi Rapo Zmax ";
	#endif
	outfile<<"OmR Omp Omz Fudge ItTorus Genfunc GenfuncAv uvOrb PAA SAA FIT\n";
	double VMax = sqrt((Pot.Phi({50.,0.,50.})-Pot.Phi(X))-.5*X[4]*X[4]);
	VecDoub range = create_log_range(0.03*VMax,0.8*VMax,500);
	int count=0;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(auto j: range){
		count+=1;
		if(count<480)
			continue;
		X[3]=j;
		X[5]=j*.8;
		O.integrate(X,10.*Pot.torb(X),0.204*Pot.torb(X));
		int guess_alpha=1;
		MatDoub FResults,ITResults,GResults,GAvResults,UVResults,PAAResults,SAAResults,FITResults;
		VecDoub Fudge, ITorus, Genfunc, GenfuncAv, uvAct, paaAct, saaAct, fitAct,Energy;
		MatDoub dvdJ_e;
		t1 = high_resolution_clock::now();
		std::vector<nanoseconds> times(8,duration_cast<nanoseconds>(t1-t1));
		for(auto i:O.results()){
			t1 = high_resolution_clock::now();
			Genfunc = AG.actions(i);
			times[2]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);GenfuncAv.resize(3);
			VecDoub aa = AG.angles(i);
			GResults.push_back({Genfunc[0],Genfunc[2],0.,0.,0.,aa[3],aa[4],aa[5]});
			Energy.push_back(Pot.H(i));
		}
		VecDoub acts = {columnMean(GResults)[0],Pot.Lz(X),columnMean(GResults)[1],columnMean(GResults)[5],columnMean(GResults)[6],columnMean(GResults)[7]};

		VecDoub GF_SD = columnSD(GResults);
		outfile<<acts[0]<<" "<<acts[1]<<" "<<acts[2]<<" "<<(acts[0]+acts[2])/fabs(acts[1])<<" ";
		#ifdef TORUS
		Actions J;J[0]=acts[0]/conv::kpcMyr2kms;
		J[2]=acts[1]/conv::kpcMyr2kms;J[1]=acts[2]/conv::kpcMyr2kms;
		Torus T; T.AutoFit(J,&TPot,1e-5);
		outfile<<T.minR()<<" "<<T.maxR()<<" "<<" "<<T.maxz()<<" ";
		#endif
		outfile<<acts[3]<<" "<<acts[4]<<" "<<acts[5]<<" "<<SD(Energy)/Mean(Energy)<<" ";

		int N=0;
		for(auto i:O.results()){
			VecDoub ang = AG.angles(i);
			t1 = high_resolution_clock::now();
			Fudge = AA.actions(i,&guess_alpha);
			times[0]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			VecDoub ang2 = AA.angles(i,&guess_alpha);
			FResults.push_back({Fudge[0]-acts[0],Fudge[2]-acts[2],ang2[0]-ang[0],ang2[1]-ang[1],ang2[2]-ang[2],ang2[3]-ang[3],ang2[4]-ang[4],ang2[5]-ang[5]});
			for(unsigned k=2;k<5;++k){
				if(FResults[N][k]>PI) FResults[N][k] = 2.*PI-FResults[N][k];
				if(FResults[N][k]<-PI) FResults[N][k] = 2.*PI+FResults[N][k];
			}
			t1 = high_resolution_clock::now();
			#ifdef TORUS
			ITorus = Tor.actions(i);
			times[1]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ITResults.push_back({ITorus[0]-acts[0],ITorus[2]-acts[2],ITorus[6]-ang[0],ITorus[7]-ang[1],ITorus[8]-ang[2],ITorus[3]-ang[3],ITorus[4]-ang[4],ITorus[5]-ang[5]});
			for(unsigned k=2;k<5;++k){
				if(ITResults[N][k]>PI) ITResults[N][k]=2.*PI-ITResults[N][k];
				if(ITResults[N][k]<-PI) ITResults[N][k]=2.*PI+ITResults[N][k];
			}
			#endif
			t1 = high_resolution_clock::now();
			GenfuncAv = AGav.actions(i);
			times[3]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			GAvResults.push_back({GenfuncAv[0]-acts[0],GenfuncAv[2]-acts[2],0.,0.,0.,0.,0.,0.});
			t1 = high_resolution_clock::now();
			uvAct = UV.actions(i);
			times[4]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ang2 = UV.angles(i);
			UVResults.push_back({uvAct[0]-acts[0],uvAct[2]-acts[2],ang2[0]-ang[0],ang2[1]-ang[1],ang2[2]-ang[2],ang2[3]-ang[3],ang2[4]-ang[4],ang2[5]-ang[5]});
			for(unsigned k=2;k<5;++k){
				if(UVResults[N][k]>PI) UVResults[N][k]=2.*PI-UVResults[N][k];
				if(UVResults[N][k]<-PI) UVResults[N][k]=2.*PI+UVResults[N][k];
			}
			t1 = high_resolution_clock::now();
			paaAct = PAA.actions(i);
			times[5]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ang2 = PAA.angles(i);
			PAAResults.push_back({paaAct[0]-acts[0],paaAct[2]-acts[2],ang2[0]-ang[0],ang2[1]-ang[1],ang2[2]-ang[2],ang2[3]-ang[3],ang2[4]-ang[4],ang2[5]-ang[5]});
			for(unsigned k=2;k<5;++k){
				if(PAAResults[N][k]>PI)PAAResults[N][k]=2.*PI-PAAResults[N][k];
				if(PAAResults[N][k]<-PI)PAAResults[N][k]=2.*PI+PAAResults[N][k];
			}
			t1 = high_resolution_clock::now();
			saaAct = SAA.actions(i,&guess_alpha);
			times[6]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ang2 = SAA.angles(i,&guess_alpha);
			SAAResults.push_back({saaAct[0]-acts[0],saaAct[2]-acts[2],ang2[0]-ang[0],ang2[1]-ang[1],ang2[2]-ang[2],ang2[3]-ang[3],ang2[4]-ang[4],ang2[5]-ang[5]});
			for(unsigned k=2;k<5;++k){
				if(SAAResults[N][k]>PI)SAAResults[N][k]=2.*PI-SAAResults[N][k];
				if(SAAResults[N][k]<-PI)SAAResults[N][k]=2.*PI+SAAResults[N][k];
			}
			t1 = high_resolution_clock::now();
			fitAct = SF.actions(i);
			times[7]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ang2 = SF.angles(i);
			FITResults.push_back({fitAct[0]-acts[0],fitAct[2]-acts[2],ang2[0]-ang[0],ang2[1]-ang[1],ang2[2]-ang[2],ang2[3]-ang[3],ang2[4]-ang[4],ang2[5]-ang[5]});
			for(unsigned k=2;k<5;++k){
				if(FITResults[N][k]>PI)FITResults[N][k]=2.*PI-FITResults[N][k];
				if(FITResults[N][k]<-PI)FITResults[N][k]=2.*PI+FITResults[N][k];
			}
			++N;
		}
		for(auto k:columnRMS(FResults)) outfile<<k<<" ";
		#ifdef TORUS
		for(auto k:columnRMS(ITResults)) outfile<<k<<" ";
		#endif
		for(auto k:columnSD(GResults)) outfile<<k<<" ";
		for(auto k:columnRMS(GAvResults)) outfile<<k<<" ";
		for(auto k:columnRMS(UVResults)) outfile<<k<<" ";
		for(auto k:columnRMS(PAAResults)) outfile<<k<<" ";
		for(auto k:columnRMS(SAAResults)) outfile<<k<<" ";
		for(auto k:columnRMS(FITResults)) outfile<<k<<" ";
		for(unsigned N=0;N<8;++N) outfile<<times[N].count()/range.size()<<" ";
		outfile<<std::endl;
	}
	outfile.close();
}

// ============================================================================
