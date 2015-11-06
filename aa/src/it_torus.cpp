// ============================================================================
/// \file src/it_torus.cpp
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
/// \brief Action finding by iteratively constructing tori
///
///	We estimate the actions by iteratively estimating them using an approximate
/// method (here the Actions_AxisymmetricStackel_Fudge) and then constructing
/// a torus with this action and finding the error made in the action estimate.
//============================================================================

/*======================================*/
/*	  Iterative Torus action finding	*/
/*======================================*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include "GSLInterface/GSLInterface.h"
#include "Fit.h"
#include "Torus.h"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "coordtransforms.h"
#include "it_torus.h"


VecDoub IterativeTorusMachine::find_aa(VecDoub x){
	for(int i=3;i<6;i++)x[i]*=conv::kpcMyr2kms;
	VecDoub V = AAA->actions(x);
	VecDoub P = AAA->angles(x);
	for(int i=0;i<3;i++)V[i]/=conv::kpcMyr2kms;
	for(int i=3;i<6;i++)P[i]/=conv::kpcMyr2kms;
	V.insert(V.end(),P.begin(),P.end());
	VecDoub vtmp = {V[3],V[4],V[5]};
	V[3]=V[6];V[4]=V[7];V[5]=V[8];
	V[6]=vtmp[0];V[7]=vtmp[1];V[8]=vtmp[2];
	return V;
}

VecDoub findXV(Angles theta, Torus *T, double sign){
	PSPT Results;
	Results = T->Map3D(theta);
	VecDoub X = {Results[0]*cos(Results[2])
				,Results[0]*sin(Results[2])
				,Results[1]
				,Results[3]*cos(Results[2])-sign*Results[5]*sin(Results[2])
				,Results[3]*sin(Results[2])+sign*Results[5]*cos(Results[2])
				,Results[4]};
	return X;
}

double PhaseOffset(const gsl_vector *v, void *p){
	Min_ItTorusMac_struct *P = (Min_ItTorusMac_struct *) p;
	Angles theta;
	for(int i=0;i<3;i++) theta[i]=gsl_vector_get(v,i);
	VecDoub XV1 = findXV(theta,P->T,P->sign);
	VecDoub Diff = P->XV-XV1;
	double Off = P->freqWeight*
				 (Diff[0]*Diff[0]+Diff[1]*Diff[1]+Diff[2]*Diff[2])
				 +(Diff[3]*Diff[3]+Diff[4]*Diff[4]+Diff[5]*Diff[5]);
	return Off;
}

VecDoub IterativeTorusMachine::RefineGuess(VecDoub PSP, VecDoub aa, double *min, Angles &theta_new, double dJ){
	// Construct torus
	int flag,err=0;
	Actions  J,oJ;	Angles theta;
	J[0] = aa[0];	J[1] = aa[2];	J[2] = aa[1];
	flag = T->AutoFit(J,&TPhi,dJ,700,300,15,5,24,200,24,err);
	if(flag!=0) std::cerr<<"Non-zero flag in IterativeTorusMachine::RefineGuess\n";
	// Use theta as initial guess
	theta[0]=aa[6];	theta[1]=aa[8];	theta[2]=aa[7];
	double sign = SIGN(PSP[0]*PSP[4]-PSP[1]*PSP[3]);
	VecDoub xout = findXV(theta,T,sign);
	// Minimise Omega^2(x_i x_i)+v_i v_i w.r.t theta
	double freqWeight =  T->omega(0)*T->omega(0)
						+T->omega(1)*T->omega(1)
						+T->omega(2)*T->omega(2);
	Min_ItTorusMac_struct Min(T,PSP,freqWeight,sign);
	VecDoub sizes = {0.05,0.05,0.05};
	VecDoub thetav = {theta[0],theta[1],theta[2]};
	minimiser M(&PhaseOffset,thetav,sizes,1e-3,&Min);
	VecDoub results;
	(*min) = M.minimise(&results,100,0);
	for(int i=0;i<3;i++) theta_new[i]=fmod(results[i],2.*PI);
	return findXV(theta_new,T,sign);
}

VecDoub IterativeTorusMachine::NewActionPoint(VecDoub startpt, VecDoub oldpt, VecDoub newpt){
	for(int i=0;i<3;i++)newpt[i]=startpt[i]+(oldpt[i]-newpt[i]);
	return newpt;
}

VecDoub IterativeTorusMachine::actions(const VecDoub &X, void*params){
	// ,int MaxIt, double et, double dJ,bool freq_yes){
	VecDoub x=X;
	for(int i=3;i<6;i++)x[i]/=conv::kpcMyr2kms;
	double min=1.; Angles theta;
	// Use Stackel to initial guess
	VecDoub aa_init = find_aa(x);
	VecDoub aaGuess=aa_init, aaNEW, aaFirst;
	int iter=0;
	while(min>eta and iter<MaxIterations){
		// Find closest PSP to x on this torus
		VecDoub xNEW = RefineGuess(x, aaGuess, &min, theta,dJ);
		// Use Stackel to find actions of this guess
		aaNEW = find_aa(xNEW);
		// Guess new based on offset from Stackel
		aaGuess = NewActionPoint(aa_init, aaGuess, aaNEW);
		iter++;
		if(iter==1)aaFirst = aaGuess;
	}

	// find frequencies of best torus
	int flag,err=0; Actions J;
	J[0] = aaGuess[0];	J[1] = aaGuess[2];	J[2] = aaGuess[1];
	flag = T->AutoFit(J,&TPhi,dJ,700,300,15,5,24,200,24,err);
	if(flag!=0) std::cerr<<"Non-zero flag in IterativeTorusMachine::actions:"<<flag<<"\n";
	for(int i=0;i<3;i++) aaGuess[i+3]=T->omega(i);

	// Return actions in right units
	for(int i=0;i<6;i++) aaGuess[i]*=conv::kpcMyr2kms;

	return {aaGuess[0],aaGuess[1],aaGuess[2], // Action results
			aaGuess[3],aaGuess[5],aaGuess[4], // Freq. results
			theta[0],theta[2],theta[1],		  // Angle results
			aa_init[0],aa_init[2],            // After first iter. guess
			aa_init[3],aa_init[5],aa_init[6],aa_init[7],aa_init[8],
			aaFirst[0],aaFirst[2],            // First action guess
			aaFirst[3],aaFirst[5],aaFirst[6],aaFirst[7],aaFirst[8],
			(double)iter
			};
}

VecDoub IterativeTorusMachine::angles(const VecDoub &X, void*params){
	VecDoub aa = actions(X,params);
	return {aa[6],aa[7],aa[8],aa[3],aa[4],aa[5]};
}
