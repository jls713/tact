// ============================================================================
/// \file src/many_tori_james.cpp
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
/// Differs from many_torus.cpp as computes standard deviations in actions and frequencies and variance
/// of angles wrt those calculated from theta0+Omega t for the given method. In many_torus.cpp the variance
/// is computed wrt to O2GF method
///
/// Must pass output file
/// e.g. mains/./many_tori_james.exe many_torus_output_james.dat
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

MatDoub dOmdJ(Actions J, Actions dJ, WrapperTorusPotential *Pot){
	Torus T; Actions Jp = J; Frequencies Omu, Omd;
	MatDoub mat(3,VecDoub(3,0.));
	Jp[0]+=dJ[0];
	T.AutoFit(Jp,Pot,1e-5); Omu = T.omega();
	Jp=J; Jp[0]-=dJ[0];
	T.AutoFit(Jp,Pot,1e-5); Omd = T.omega();
	mat[0][0]=(Omu[0]-Omd[0])/2./dJ[0];
	mat[1][0]=(Omu[2]-Omd[2])/2./dJ[0];
	mat[2][0]=(Omu[1]-Omd[1])/2./dJ[0];
	Jp=J; Jp[1]+=dJ[1];
	T.AutoFit(Jp,Pot,1e-5); Omu = T.omega();
	Jp=J; Jp[2]-=dJ[2];
	T.AutoFit(Jp,Pot,1e-5); Omd = T.omega();
	mat[1][1]=(Omu[2]-Omd[2])/2./dJ[2];
	mat[2][1]=(Omu[1]-Omd[1])/2./dJ[2];
	Jp=J; Jp[1]+=dJ[1];
	T.AutoFit(Jp,Pot,1e-5); Omu = T.omega();
	Jp=J; Jp[1]-=dJ[1];
	T.AutoFit(Jp,Pot,1e-5); Omd = T.omega();
	mat[2][2]=(Omu[1]-Omd[1])/2./dJ[1];
	mat[0][1]=mat[1][0];
	mat[0][2]=mat[2][0];
	mat[1][2]=mat[2][1];
	return mat;
}
# endif

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

	if(argc>2)
		X[0]=atof(argv[2]);
	else
		X[0]=conv::StandardSolarPAUL[0];
	X[2]=0.001;
	X[4]=sqrt(X[0]*-Pot.Forces(X)[0]);
	printVector(X);
	Orbit O(&Pot,1e-8);
	// Fudge
	Actions_AxisymmetricStackel_Fudge AA(&Pot,-30.);

	// Iterative Torus
	#ifdef TORUS
	IterativeTorusMachine Tor(&AA,&Pot,1e-8,5,1e-3);
	#endif

	// Generating Function
	Actions_Genfunc AG(&Pot,"axisymmetric");

	// Average generating Function
	Actions_Genfunc_Average AGav(&Pot,"axisymmetric");

	// uvorb
	uv_orb UV(&Pot,4.,30.,50,50,"example.delta_uv");

	// Cylindrical Adiabatic
	Actions_CylindricalAdiabaticApproximation PAA(&Pot,"example.paa",true,false,4.,30.,15.,100);

	// Spheroidal Adiabatic
	Actions_SpheroidalAdiabaticApproximation SAA(&Pot,"example.saa",true,false,100.,4.,30.,15.,100);

	// Spheroidal Adiabatic
	Actions_StackelFit SF(&Pot,1e-8);

	std::ofstream outfile;
	outfile.open(argv[1]);
	outfile<<"# JR Lz Jz JRJzLz ";
	#ifdef TORUS
	outfile<<"Rperi Rapo Zmax ";
	#endif
	outfile<<"OmR Omp Omz Fudgev1 ItTC O2GF AvGF Fudgev2 CAA SAA Fit\n";
	double VMax = sqrt((Pot.Phi({50.,0.,50.})-Pot.Phi(X))-.5*X[4]*X[4]);
	int number = 500;
	if(argc>3)
		number=atoi(argv[3]);
	VecDoub range = create_log_range(0.03*VMax,0.8*VMax,number);
	int count=0;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(auto j: range){
		count+=1;
		// if(count<480)
		// 	continue;
		X[3]=j;
		X[5]=j*.8;
		printVector(X);
		double Torb = Pot.torb(X), tstep=0.204*Torb, tmax=10.*Torb;
		O.integrate(X,tmax,tstep);
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
			GResults.push_back({Genfunc[0],Genfunc[2],aa[0],aa[1],aa[2],aa[3],aa[4],aa[5]});
			Energy.push_back(Pot.H(i));
		}
		VecDoub acts = {columnMean(GResults)[0],Pot.Lz(X),columnMean(GResults)[1],columnMean(GResults)[5],columnMean(GResults)[6],columnMean(GResults)[7]};

		VecDoub GF_SD = columncarefulSD(GResults);
		outfile<<acts[0]<<" "<<acts[1]<<" "<<acts[2]<<" "<<(acts[0]+acts[2])/fabs(acts[1])<<" ";
		#ifdef TORUS
		Actions J;J[0]=acts[0]/conv::kpcMyr2kms;
		J[2]=acts[1]/conv::kpcMyr2kms;J[1]=acts[2]/conv::kpcMyr2kms;
		Torus T; T.AutoFit(J,&TPot,1e-5);
		outfile<<T.minR()<<" "<<T.maxR()<<" "<<" "<<T.maxz()<<" ";
		MatDoub Hess = dOmdJ(J,.1*J,&TPot);
		#endif
		outfile<<acts[3]<<" "<<acts[4]<<" "<<acts[5]<<" "<<carefulSD(Energy)/Mean(Energy)<<" ";

		int N=0;
		for(auto i:O.results()){
			VecDoub ang = AG.angles(i);
			t1 = high_resolution_clock::now();
			Fudge = AA.actions(i,&guess_alpha);
			times[0]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			VecDoub ang2 = AA.angles(i,&guess_alpha);
			FResults.push_back({Fudge[0],Fudge[2],ang2[0],ang2[1],ang2[2],ang2[3],ang2[4],ang2[5]});
			t1 = high_resolution_clock::now();
			#ifdef TORUS
			ITorus = Tor.actions(i);
			times[1]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ITResults.push_back({ITorus[0],ITorus[2],ITorus[6],ITorus[7],ITorus[8],ITorus[3],ITorus[4],ITorus[5]});
			#endif
			t1 = high_resolution_clock::now();
			GenfuncAv = AGav.actions(i);
			times[3]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			GAvResults.push_back({GenfuncAv[0],GenfuncAv[2],0.,0.,0.,0.,0.,0.});
			t1 = high_resolution_clock::now();
			uvAct = UV.actions(i);
			times[4]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ang2 = UV.angles(i);
			UVResults.push_back({uvAct[0],uvAct[2],ang2[0],ang2[1],ang2[2],ang2[3],ang2[4],ang2[5]});
			t1 = high_resolution_clock::now();
			paaAct = PAA.actions(i);
			times[5]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ang2 = PAA.angles(i);
			PAAResults.push_back({paaAct[0],paaAct[2],ang2[0],ang2[1],ang2[2],ang2[3],ang2[4],ang2[5]});
			t1 = high_resolution_clock::now();
			saaAct = SAA.actions(i,&guess_alpha);
			times[6]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ang2 = SAA.angles(i,&guess_alpha);
			SAAResults.push_back({saaAct[0],saaAct[2],ang2[0],ang2[1],ang2[2],ang2[3],ang2[4],ang2[5]});
			t1 = high_resolution_clock::now();
			fitAct = SF.actions(i);
			times[7]+=duration_cast<nanoseconds>(high_resolution_clock::now()-t1);
			ang2 = SF.angles(i);
			FITResults.push_back({fitAct[0],fitAct[2],ang2[0],ang2[1],ang2[2],ang2[3],ang2[4],ang2[5]});
			++N;
		}

		double timeT=tstep;VecDoub freqs;
		for(int i=1;i<N;++i){
			freqs=columnMean(FResults)*timeT;
			FResults[i][2]-=FResults[0][2]+freqs[5];
			FResults[i][3]-=FResults[0][3]+freqs[6];
			FResults[i][4]-=FResults[0][4]+freqs[7];
			freqs=columnMean(ITResults)*timeT;
			ITResults[i][2]-=ITResults[0][2]+freqs[5];
			ITResults[i][3]-=ITResults[0][3]+freqs[6];
			ITResults[i][4]-=ITResults[0][4]+freqs[7];
			freqs=columnMean(GResults)*timeT;
			GResults[i][2]-=GResults[0][2]+freqs[5];
			GResults[i][3]-=GResults[0][3]+freqs[6];
			GResults[i][4]-=GResults[0][4]+freqs[7];
			freqs=columnMean(UVResults)*timeT;
			UVResults[i][2]-=UVResults[0][2]+freqs[5];
			UVResults[i][3]-=UVResults[0][3]+freqs[6];
			UVResults[i][4]-=UVResults[0][4]+freqs[7];
			freqs=columnMean(PAAResults)*timeT;
			PAAResults[i][2]-=PAAResults[0][2]+freqs[5];
			PAAResults[i][3]-=PAAResults[0][3]+freqs[6];
			PAAResults[i][4]-=PAAResults[0][4]+freqs[7];
			freqs=columnMean(SAAResults)*timeT;
			SAAResults[i][2]-=SAAResults[0][2]+freqs[5];
			SAAResults[i][3]-=SAAResults[0][3]+freqs[6];
			SAAResults[i][4]-=SAAResults[0][4]+freqs[7];
			freqs=columnMean(FITResults)*timeT;
			FITResults[i][2]-=FITResults[0][2]+freqs[5];
			FITResults[i][3]-=FITResults[0][3]+freqs[6];
			FITResults[i][4]-=FITResults[0][4]+freqs[7];
			timeT+=tstep;
			for(unsigned k=2;k<5;++k){
				while(FITResults[i][k]<-PI)FITResults[i][k]+=2.*PI;
				while(GResults[i][k]<-PI)GResults[i][k]+=2.*PI;
				while(FResults[i][k]<-PI)FResults[i][k]+=2.*PI;
				while(ITResults[i][k]<-PI)ITResults[i][k]+=2.*PI;
				while(UVResults[i][k]<-PI)UVResults[i][k]+=2.*PI;
				while(SAAResults[i][k]<-PI)SAAResults[i][k]+=2.*PI;
				while(PAAResults[i][k]<-PI)PAAResults[i][k]+=2.*PI;
			}
		}
		for(unsigned k=2;k<5;++k) FResults[0][k]=0.;
		for(unsigned k=2;k<5;++k) UVResults[0][k]=0.;
		for(unsigned k=2;k<5;++k) GResults[0][k]=0.;
		for(unsigned k=2;k<5;++k) ITResults[0][k]=0.;
		for(unsigned k=2;k<5;++k) FITResults[0][k]=0.;
		for(unsigned k=2;k<5;++k) PAAResults[0][k]=0.;
		for(unsigned k=2;k<5;++k) SAAResults[0][k]=0.;

		for(auto k:columncarefulSD(FResults)) outfile<<k<<" ";
		#ifdef TORUS
		for(auto k:columncarefulSD(ITResults)) outfile<<k<<" ";
		#endif
		for(auto k:columncarefulSD(GResults)) outfile<<k<<" ";
		for(auto k:columncarefulSD(GAvResults)) outfile<<k<<" ";
		for(auto k:columncarefulSD(UVResults)) outfile<<k<<" ";
		for(auto k:columncarefulSD(PAAResults)) outfile<<k<<" ";
		for(auto k:columncarefulSD(SAAResults)) outfile<<k<<" ";
		for(auto k:columncarefulSD(FITResults)) outfile<<k<<" ";
		for(unsigned N=0;N<8;++N) outfile<<times[N].count()/range.size()<<" ";
		#ifdef TORUS
		for(unsigned kk=0;kk<3;++kk)
			for(unsigned pp=0;pp<3;++pp)
				outfile<<Hess[kk][pp]<<" ";
		#endif
		outfile<<std::endl;
	}
	outfile.close();
}

// ============================================================================
