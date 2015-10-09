// ============================================================================
/// \file src/lmn_orb.cpp
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
/// \brief Axisymmetric Staeckel fudge using Delta estimation from shells
///
/// lmn_orb: Wraps triaxial Staeckel fudge using the alpha, beta (or Delta_1
/// and Delta_2) estimation routine from Sanders & Binney (2014).
/// We find the closed loop orbits in the (x,y) and (y,z) planes and fit
/// ellipses to these orbits to find alpha and beta respectively.
///
//============================================================================

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
#include "stackel_aa.h"
#include "lmn_orb.h"
#include "debug.h"
// What is gsl_set_error_handler_off() doing? l125
// ============================================================================

static double min_distance(double y, void *params){

	root_struct_mindist *RS = (root_struct_mindist *) params;
	Orbit orbit(RS->Pot);
	VecDoub X = {0.,y,0.,0.,0.,0.};
	double p = sqrt(2.*(RS->E-RS->Pot->Phi(X)));
	if(RS->swit==0) // Short axis
		X[3]=p;
	else if(RS->swit==1) // Long axis
		X[5]=p;
	VecDoub QQ = X; VecDoub QQ2;
	// Choose appropriate time-step
	double torb = RS->Pot->torb(X);
	double step = 1e-3*torb;
	int i=0, maxsteps = 10000;
	while(i<maxsteps){
		QQ2=orbit.integrate(QQ, step,step);
		if(RS->swit==0 and QQ[0]*QQ2[0]<0. and QQ2[1]<0. and QQ2[3]<0.) break;
		if(RS->swit==1 and QQ[2]*QQ2[2]<0. and QQ2[1]<0. and QQ2[5]<0.) break;
		QQ = QQ2;
		i++;
	}
	double min = 4.*PI*PI/torb/torb*(-QQ[1]-X[1])*(-QQ2[1]-X[1]);
	if(RS->swit==0)min+=(-QQ[3]-X[3])*(-QQ[3]-X[3]);
	else if(RS->swit==1)min+=(-QQ[5]-X[5])*(-QQ[5]-X[5]);
	return min;
}

std::vector<int> lmn_orb::angular_momentum(const VecDoub &x){
	VecDoub xx = {x[0],x[1],x[2]}, vv = {x[3],x[4],x[5]};
	VecDoub ll = cross_product<double>(xx,vv);
	return {sign(ll[0]),sign(ll[1]),sign(ll[2])};
}

int lmn_orb::check_ang_mom(double y, double E, int swit){
	double P = Pot->Phi({0.,y,0.});
	double p = sqrt(2.*(E-P));
	Orbit orbit(Pot);
	VecDoub X = {0.,y,0.,0.,0.,0.};
	int index = 0; if(swit==0) index = 2;
	if(swit==0) X[3]=p;
	else if(swit==1) X[5]=p;
	double step = 1e-2*Pot->torb(X);
	VecDoub QQ=orbit.integrate(X, 10000*step, step);
	orbit.plot(0,1);
	int l=1, l0=angular_momentum(orbit.results()[0])[index];
    int result = 1;
	for(auto it = std::begin(orbit.results())+1;
	    	 it!= std::end(orbit.results());
	    	 ++it)
	{
		l = angular_momentum(*it)[index];
		if(l!=l0){ result = 0; break;}
	}
	return result;
}

VecDoub lmn_orb::check_orbit(double y, double E, int swit, int plot){
	// Plots closed orbit to check if it is closed

	double p = sqrt(2.*(E-Pot->Phi({0.,y,0.})));
	Orbit orbit(Pot);
	VecDoub X = {0.,y,0.,0.,0.,0.};
	if(swit==0) X[3]=p;
	else if(swit==1) X[5]=p;
	double step = 1e-4*Pot->torb(X);
	VecDoub QQ=orbit.integrate(X, 100000*step, step);
	double ymax = Max<double>(transpose(orbit.results())[swit]);
	double zmax = Max<double>(transpose(orbit.results())[swit+1]);
	if(debug_find_Delta){
		orbit.output2file("orbit.tmp");
		std::cin.ignore();
		// if(swit==0 and plot==1)orbit.plot(0,1);
		// if(swit==1) orbit.plot(1,2);
	}
	return {ymax,zmax};
}

static double EminusPot(double x, void *p){
	root_struct_mindist *RS = (root_struct_mindist *) p;
	return RS->E-RS->Pot->Phi({0.,x,0.});
}


double lmn_orb::find_closed(double E, int swit){
	// At a fixed energy finds the short (swit=0) or long (swit=1)
	// axis closed loop orbit
	root_struct_mindist RS(Pot,E,swit);
	double yMax=ymax,yMin=1e-3, mid=0.;
	root_find RF(1e-4,100);
	yMax = RF.findroot(&EminusPot,1e-5,ymax,&RS)*.99;
	mid = yMax*0.7;
	double Dmax=min_distance(yMax, &RS),Dmin=min_distance(yMin, &RS),Dmid=min_distance(mid, &RS);

	// // Check if only boxes
	// if(fabs((Dmax-Dmin)/(yMax-yMin))-fabs((Dmid-Dmin)/(mid-yMin))<.5)
	//    return -1.;

	bracketer bra;
	int brr = bra.bracket(&yMax,&mid,&yMin,&Dmax,&Dmid,&Dmin,&min_distance,&RS);
	if(yMin<0. or yMax<0. or mid<0. or brr==0)
		return -1.;
	if(debug_find_Delta)
		std::cout<<yMin<<" "<<mid<<" "<<yMax<<" "<<Dmin<<" "<<Dmid<<" "<<Dmax<<std::endl;
	int status=0;
	minimiser1D min(&min_distance,mid,yMin,yMax,1e-4,0.,&status,&RS);
	if(status!=0) return -1.;
	double y = min.minimise(100);
	// if(swit==1) check_orbit(y,E,swit);
	return y;
}

static double min_action(double b, void *params){
	root_struct_action *RS = (root_struct_action *) params;
	if(RS->swit==0)RS->ATSF->CS->newalpha(b);
	else if(RS->swit==1)RS->ATSF->CS->newbeta(b);
	double L = (double)RS->orbit->results().size();
	double means=0.,vars=0.;int index = RS->swit+1;
	for(auto Y: RS->orbit->results()){
		VecDoub i = RS->ATSF->actions(Y);
		if(i[index]==i[index]){means+=i[index];	vars+=i[index]*i[index];}
	}
	double f = sqrt(MAX(0.,vars/L-means*means/L/L));
	return f;
}

static double min_lambda(double b, void *params){
	root_struct_action *RS = (root_struct_action *) params;
	if(RS->swit==0)RS->ATSF->CS->newalpha(b);
	else if(RS->swit==1)RS->ATSF->CS->newbeta(b);
	double L = (double)RS->orbit->results().size();
	double means=0.,vars=0.;int index = 0;
	for(auto Y: RS->orbit->results()){
		VecDoub i = RS->ATSF->CS->xv2tau(Y);
		i[index]+=RS->ATSF->CS->alpha();
		if(i[index]==i[index]){means+=i[index];	vars+=i[index]*i[index];}
	}
	double f = sqrt(MAX(0.,vars/L-means*means/L/L));
	return f;
}

void lmn_orb::plot_Delta2(double E){
	double y = find_closed(E,1);
	double P = Pot->Phi({0.,y,0.});
	double p = sqrt(2.*(E-P));
	Orbit orbit(Pot);
	VecDoub X = {0.0,y,0.0,0.0,0.0,p};
	VecDoub QQ=orbit.integrate(X, 2.5, 0.005);
	orbit.output2file("orbit_Delta2_loweralpha.dat");
	Actions_TriaxialStackel_Fudge ATSF(Pot,-80.,-20.);
	root_struct_action RS(&ATSF,&orbit,1);
	std::ofstream outFile; outFile.open("Delta2_loweralpha.dat");
	for(double b = ATSF.CS->alpha()+0.01;b<ATSF.CS->gamma()-0.01;b+=1.)
		outFile<<sqrt(ATSF.CS->gamma()-b)<<" "<<min_action(b,&RS)<<std::endl;
	outFile.close();
	return;
}

void lmn_orb::plot_Delta1(double E){
	double y = find_closed(E,0);
	double P = Pot->Phi({0.,y,0.});
	double p = sqrt(2.*(E-P));
	Orbit orbit(Pot);
	VecDoub X = {0.0,y,0.0,p,0.0,0.0};
	VecDoub QQ=orbit.integrate(X, 2.5, 0.005);
	orbit.output2file("orbit_Delta1.dat");
	Actions_TriaxialStackel_Fudge ATSF(Pot,-80.,-20.);
	root_struct_action RS(&ATSF,&orbit,0);
	std::ofstream outFile; outFile.open("Delta1.dat");
	for(double a = ATSF.CS->alpha()+1.;a<ATSF.CS->beta()-1.;a+=0.5)
		outFile<<sqrt(ATSF.CS->beta()-a)<<" "<<min_action(a,&RS)<<std::endl;
	outFile.close();
	return;
}

double lmn_orb::find_beta(double E, double a_init, double b_init){
	if(E<E0){
		std::cerr<<"find_beta: Energy too low"<<std::endl;
		return 0.;
	}
	double y = find_closed(E,1);
	if(y<0.) return 10.; // We have failed to find loop orbit at this energy

	// Use orbit shape -- assume elliptical -- default
	if(!use_acts){
		VecDoub orbit_shape = check_orbit(y,E,1,1);
		if(orbit_shape[1]<orbit_shape[0])
			std::cout<<"Orbit wrong shape!!!: "
					 <<orbit_shape[0]<<" "<<orbit_shape[1]<<" "
					 <<Pot->Phi({0.,orbit_shape[0],0.})<<" "
					 <<Pot->Phi({0.,0.,orbit_shape[1]})<<std::endl;
		return -1-(orbit_shape[1]*orbit_shape[1]-orbit_shape[0]*orbit_shape[0]);
	} // or minimise actions
	else{
		double P = Pot->Phi({0.,y,0.});
		double p = sqrt(2.*(E-P));
		Orbit orbit(Pot);
		VecDoub X = {1e-8,y,1e-8,0.,0.,p};
		double torb = Pot->torb(X);
		double step = 1e-4*torb;
		int maxsteps = 100000;
		VecDoub QQ=orbit.integrate(X, maxsteps*step, step);
		Actions_TriaxialStackel_Fudge ATSF(Pot, a_init, b_init);
		root_struct_action RS(&ATSF,&orbit,1);
		double a = -1.-y*y;//
		// a = ATSF.CS->alpha()+0.0001;
		double g = ATSF.CS->gamma()-1e-6;

		auto min_func = /*&min_action;//*/&min_lambda;

		double down = min_func(a,&RS), up = min_func(g,&RS);
		double av = (a+g)/2.;
		double avM = min_func(av,&RS);

		if(debug_find_Delta){
			std::cerr<<"Beta finder at energy "<<E<<": ";
			std::cerr<<"beta = "<<a<<", spread = "<<down;
			std::cerr<<", beta = "<<av<<", spread = "<<avM;
			std::cerr<<", beta = "<<g<<", spread = "<<up<<std::endl;
		}
		while(avM>down){
			av = (av+a)/2.;avM = min_func(av,&RS);
		}
		while(avM>up){
			av = (av+g)/2.;avM = min_func(av,&RS);
		}
		if(avM==up) return -1.001;
		else{
			int status = 0;
			minimiser1D min(min_func,av,a,g,1e-3,0.,&status,&RS);
			if(status!=0) std::cerr<<"Minimum not encompassed in find_beta\n";
			return min.minimise(100);
		}
	}
}

double lmn_orb::find_alpha(double E, double alpha_i, double beta){

	if(E<E0){
		std::cerr<<"find_alpha: Energy too low"<<std::endl;
		return 0.;
	}

	double y = find_closed(E,0);
	if(y<0.) return 10.; // We have failed to find loop orbit at this energy


	// We either use the shape of the orbit -- assume elliptical --default
	if(!use_acts){
		VecDoub orbit_shape = check_orbit(y,E,0,1);
		if(orbit_shape[1]<orbit_shape[0])
			std::cout<<"Orbit wrong shape!!!: "
					 <<orbit_shape[0]<<" "<<orbit_shape[1]<<" "
					 <<Pot->Phi({orbit_shape[0],0.,0.})<<" "
					 <<Pot->Phi({0.,orbit_shape[1],0.})<<std::endl;
		return beta-(orbit_shape[1]*orbit_shape[1]-orbit_shape[0]*orbit_shape[0]);
	}// or minimising the actions
	else{
		double P = Pot->Phi({0.,y,0.});
		double p = sqrt(2.*(E-P));
		Orbit orbit(Pot);
		VecDoub X = {1e-8,y,1e-8,p,0.,0.};
		double torb = Pot->torb(X);
		double step = 1e-4*torb;
		int maxsteps = 100000;
		VecDoub QQ=orbit.integrate(X, maxsteps*step, step);
		Actions_TriaxialStackel_Fudge ATSF(Pot, alpha_i, beta);
		root_struct_action RS(&ATSF,&orbit,0);
		double b = ATSF.CS->beta()-1e-5;
		double a1 = 10.*b;
		double av = b-1.;

		auto min_func = /*&min_action;//*/&min_lambda;

		double down = min_func(a1,&RS), up = min_func(b,&RS), avM = min_func(av,&RS);
		while(avM>down){
			av = (av+a1)/2.;avM = min_func(av,&RS);
		}
		while(avM>up){
			av = (av+b)/2.;avM = min_func(av,&RS);
		}
		if(avM==down or avM==up) return -1.002;
		else{
			int status = 0;
			minimiser1D min(min_func,av,a1,b,5e-2,0.,&status,&RS);
			if(status!=0) std::cerr<<"Minimum not encompassed in find_alpha\n";
			return min.minimise(100);
		}
	}
}

VecDoub lmn_orb::find_ab_from_box(double E){
	root_struct_mindist RS(Pot,E,0);
	root_find RF(1e-3,100);
	double y = RF.findroot(&EminusPot,1e-5,ymax,&RS);
	double p = sqrt(2.*(E-Pot->Phi({0.,y*.7,0.})));
	return find_best_alphabeta({0.,y*.7,0.,p/sqrt(2.),0.,p/sqrt(2.)});
}

static VecDoub actionSD(Actions_TriaxialStackel_Fudge *ATSF, const std::vector<VecDoub> & results, bool with_freq=false){

	if(with_freq) ATSF->set_freq(true);
	std::vector<VecDoub> Results;
	for(auto Y: results) Results.push_back(ATSF->actions(Y));

	return concatVectors<double>(columnMean<double>(Results),
	                             columnSD<double>(Results),
	                             columnMedian<double>(Results));
}

static double min_actions(const gsl_vector *v, void *params){
	root_struct_actions *RS = (root_struct_actions *) params;
	if(gsl_vector_get(v,1)<gsl_vector_get(v,0) or gsl_vector_get(v,1)>-1. or gsl_vector_get(v,0)>-1.)return 1e10;
	Actions_TriaxialStackel_Fudge ATSF(RS->Pot,gsl_vector_get(v,0),gsl_vector_get(v,1));
	VecDoub f = actionSD(&ATSF,RS->orbit->results());
	return sqrt(f[3]*f[3]+f[4]*f[4]+f[5]*f[5]);
}

void lmn_orb::alphabeta_grid(const VecDoub& X){
	std::ofstream outfile; outfile.open("grid.grid");
	double En = Pot->H(X);
	findDelta_interp(En);
	Orbit orbit(Pot);
	VecDoub QQ=orbit.integrate(X, 5., 0.1);
	root_struct_actions RS(Pot,&orbit);
	gsl_vector *v; v = gsl_vector_alloc(2);
	for(double alpha = -20.;alpha<-1.1;alpha+=0.2){
		gsl_vector_set(v,0,alpha);
		for(double beta = alpha+0.001;beta<-1.01;beta+=0.2){
			gsl_vector_set(v,1,beta);
			outfile<<sqrt(beta-alpha)<<" "<<sqrt(-1.-beta)<<" "<<min_actions(v,&RS)<<std::endl;
		}
	}
	outfile.close();
	gsl_vector_free(v);
	return;
}

VecDoub lmn_orb::find_best_alphabeta(const VecDoub& X){
	double torb = Pot->torb(X);
	double step = 5e-1*torb;
	Orbit orbit(Pot);
	VecDoub QQ=orbit.integrate(X, 100*step, step);
	orbit.plot(1,2);
	root_struct_actions RS(Pot,&orbit);
	VecDoub a1 = {-1.5,-1.1};
	VecDoub sizes ={.5,.1};
	minimiser min(&min_actions,a1,sizes,1e-5,&RS);
	VecDoub results;
	min.minimise(&results,100,0);
	return results;
}



void lmn_orb::readDeltagrids(const std::string& file){
	std::ifstream infile; infile.open(file);
	if(!infile.is_open())std::cerr<<"Problem: "<<file<<" doesn't exist."<<std::endl;
	double tmp,tmp2,tmp3;
	while(infile>>tmp>>tmp2>>tmp3){
		E.push_back(tmp);
		Beta.push_back(tmp2);
		Alpha.push_back(tmp3);
	}
	infile.close();
}

void lmn_orb::fillDeltagrids(const std::string& file){
	E=VecDoub(NE,0.);Beta=E;Alpha=E;

	if(use_log_grid)
		E = create_log_range(E0,Emax,NE);
	else
		E = create_range(E0,Emax,NE);

	#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<NE;i++){
		double alpha=-1.2, beta=-1.05;
		beta = find_beta(E[i],alpha,beta);
		if(beta>1.){
			// These are energies where there are no loop orbits
			beta = 1.;alpha = 1.;
		}
		else alpha = find_alpha(E[i],beta-0.001, beta);
		// just check they are not the same
		if(beta>-1.){
			// double xin = find_closed(E[i],1);
			beta=-1.05;
			// while(sqrt(-1.-beta)<.5*xin) beta=-1.-(-1.-beta)*2.;
		}
		if(alpha>=beta){
			// double xin = find_closed(E[i],0);
			alpha=beta-0.1;//-0.05;
			// while(sqrt(beta-alpha)<.5*xin) alpha=beta-(beta-alpha)*2.;
		}
		Alpha[i] = alpha;
		Beta[i] = beta;
		std::cerr<<E[i]<<" "<<beta<<" "<<alpha<<std::endl;
	}
	// Extrapolate back for those with no loop orbits
	int n_a=-1, n_b=-1;
	for(unsigned i=0;i<E.size();i++)
		if(Beta[i]<0.){ n_b = i; break; }
	for(unsigned i=0;i<E.size();i++)
		if(Alpha[i]<0.){ n_a = i; break; }
	if(n_a<0 and n_b<0)
		std::cerr<<"ALL ENERGIES DISALLOW LOOP ORBITS!!!!"<<std::endl;

	for(int i=0;i<n_a;i++)
		Alpha[i] = (Alpha[n_a+1]-Alpha[n_a])/(E[n_a+1]-E[n_a])*(E[i]-E[n_a])+Alpha[n_a];
	for(int i=0;i<n_b;i++)
		Beta[i] = (Beta[n_b+1]-Beta[n_b])/(E[n_b+1]-E[n_b])*(E[i]-E[n_b])+Beta[n_b];

	std::ofstream outfile; outfile.open(file);
	for(unsigned i=0;i<E.size();i++)
		outfile<<E[i]<<" "<<Beta[i]<<" "<<Alpha[i]<<std::endl;
	outfile.close();

}

VecDoub lmn_orb::findDelta_interp(double En){
	double a, b;
	// if(En>Emax){
	// 	double interp = (log(-Alpha[NE-1])-log(-Alpha[NE-2]))/(log(-E[NE-1])-log(-E[NE-2]));
	// 	a = pow(En/E[NE-1],interp)*Alpha[NE-1];
	// 	interp = (log(-Beta[NE-1])-log(-Beta[NE-2]))/(log(-E[NE-1])-log(-E[NE-2]));
	// 	b = pow(En/E[NE-1],interp)*Beta[NE-1];
	// }
	if(En>Emax){
		a = Alpha[NE-1];
		b = Beta[NE-1];
	}
	else if(En<E0){
		a = Alpha[0];
		b = Beta[0];
	}
	else{
		int i=0;
		while(E[i]<En and i<NE-1)	i++;
		a = (Alpha[i]-Alpha[i-1])/(E[i]-E[i-1])*(En-E[i-1])+Alpha[i-1];
		b = (Beta[i]-Beta[i-1])/(E[i]-E[i-1])*(En-E[i-1])+Beta[i-1];
	}
	return {a,b};
}

VecDoub lmn_orb::actions(const VecDoub& x,void*params){

	double En = Pot->H(x);
	if(params!=nullptr){
		double *alphabeta = (double*)params;
		Actions_TriaxialStackel_Fudge ATSF(Pot,alphabeta[0],alphabeta[1]);
		return ATSF.actions(x);
	}
	else{
		VecDoub ab = findDelta_interp(En);


		// // if alpha>beta or beta>gamma then the axes are not in the right order
		// // so we must flip the axes and change alpha and beta
		// if(ab[0]>ab[1]){
		// 	if(ab[0]<-1.){
		// 		// This is the case when z is minor, x inter, y major
		// 		VecDoub X = x;
		// 		double tmp=ab[1];ab[1]=ab[0];ab[0]=tmp;
		// 		tmp=X[1];X[1]=X[0];X[0]=tmp;
		// 		Actions_TriaxialStackel_Fudge ATSF(Pot,ab[0],ab[1]);
		// 		return ATSF.actions(x);
		// 	}
		// 	else if(ab[1]>-1.)
		// 		// This is the case when x is minor, y inter, z major
		// 		VecDoub X = x;
		// 		double tmp=ab[2];ab[2]=ab[0];ab[0]=tmp;
		// 		tmp=X[2];X[2]=X[0];X[0]=tmp;
		// 	}
		// 	else{
		// 		// This is the case when x is minor, z inter, y major
		// 		VecDoub X = x;
		// 		double tmp=ab[1];ab[1]=ab[0];ab[0]=tmp;
		// 		ab[1]=-2.-ab[1];
		// 		tmp=X[1];X[1]=X[2];X[2]=tmp;
		// 		tmp=X[1];X[1]=X[0];X[0]=tmp;
		// 		Actions_TriaxialStackel_Fudge ATSF(Pot,ab[0],ab[1]);
		// 		return ATSF.actions(x);
		// 	}
		// }
		// else if(ab[1]>-1.){
		// 	if(ab[0]<-1.){
		// 		// This is the case when y is minor, z inter, x major
		// 		VecDoub X = x;
		// 		ab[1]=-2.-ab[1];
		// 		double tmp=X[2];X[2]=X[1];X[1]=tmp;
		// 		Actions_TriaxialStackel_Fudge ATSF(Pot,ab[0],ab[1]);
		// 		return ATSF.actions(x);
		// 	}
		// 	else{
		// 		// This is the case when y is minor, x inter, z major
		// 		VecDoub X = x;
		// 	}
		// }
		Actions_TriaxialStackel_Fudge ATSF(Pot,ab[0],ab[1]);
		return ATSF.actions(x);
	}
}

VecDoub lmn_orb::angles(const VecDoub& x,void *params){
	double En = Pot->H(x);
	if(params!=nullptr){
		double *alphabeta = (double*)params;
		Actions_TriaxialStackel_Fudge ATSF(Pot,alphabeta[0],alphabeta[1]);
		return ATSF.angles(x);
	}
	else{
		VecDoub ab = findDelta_interp(En);
		Actions_TriaxialStackel_Fudge ATSF(Pot,ab[0],ab[1]);
		return ATSF.angles(x);
	}
}

double lmn_orb::sos(int comp, const VecDoub& x, const std::string& outfile){
	double En = Pot->H(x);
	VecDoub alphabeta = findDelta_interp(En);
	Actions_TriaxialStackel_Fudge ATSF(Pot,alphabeta[0],alphabeta[1]);
	return ATSF.sos(x,comp,outfile);
}
// ============================================================================
