// ============================================================================
/// \file src/orbit.cpp
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
/// \brief Orbit integration classes
///
/// Orbit_base implements a base class for common functionality amongst orbit
/// integration class and there are two inherited classes
/// 1. Orbit -- integrates using GSL ode routines (default Dortmund-Prince 8)
/// 2. SymplecticOrbit -- integrates using 2nd order kick-drift-kick
///
// ============================================================================

/*======================================*/
/* 				  Orbits	 			*/
/*======================================*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface/GSLInterface.h"
#include "utils.h"
#include "gnuplot/gnuplot_i.h"
#include "potential.h"
#include "orbit.h"

// ============================================================================
// Orbit base class -- implements 3D orbit integration in a Potential_JS
// ============================================================================

double Orbit_base::dynamic_time(const VecDoub& x){
	VecDoub F = Pot->Forces(x);
	double t = 0., n=0.;
	for(int i=0;i<3;i++){t+=-F[i]*x[i]; n+=x[i]*x[i];}
	return sqrt(n/t);
}

void Orbit_base::plot(int i, int j, std::string name){
	std::vector<std::string> labels = {"x","y","z","vx","vy","vz"};
	// return;
// }
	/* Plot the results of the orbit integration  */
	/* i,j give the indices of the arrays to plot */
	VecDoub x,y; std::string x_s=labels[i], y_s=labels[j];
	for(unsigned int k=0;k<results_.size();k++){
		x.push_back(results_[k][i]);y.push_back(results_[k][j]);
	}
	Gnuplot G("lines ls 1");
	G.set_xrange(1.1*Min(x),1.1*Max(x)).set_yrange(1.1*Min(y),1.1*Max(y));
	G.set_xlabel(x_s).set_ylabel(y_s);
	G.savetotex(name).plot_xy(x,y);
	G.outputpdf(name);
}

void Orbit_base::output2file(std::string file){
	std::ofstream outFile; outFile.open(file);
	for(auto i: results_){
		for(auto j:i)outFile<<j<<" ";
		outFile<<std::endl;
	}
}

// ============================================================================
// Orbit integration using GSL ode integrators
// ============================================================================

int derivs(double t, const double y[], double f[],void *params){
	/* derivatives for orbit integration */
	Potential_JS *P = (Potential_JS*) params;
	VecDoub x = {y[0],y[1],y[2]};
	VecDoub F = P->Forces(x);
	for(int i=0;i<3;i++){ f[i]=y[i+3]; f[i+3]=F[i]; }
	return GSL_SUCCESS;
}

Orbit::Orbit(Potential_JS *P, double eps,std::string type):Orbit_base(P){
	/* Create instance of orbit with associated Potential_JS P */
	O = std::unique_ptr<ode>(new ode(derivs,6,eps,type,P));
}

VecDoub Orbit::integrate(const VecDoub& x, double t_interval, double step_size, bool adaptive){
	/* Integrate from time t=0 to t=t_interval at regular steps step_size */
	/* and store results in the array results							  */
	double y[6]; VecDoub xf(6); for(int i=0;i<6;i++)y[i]=x[i];
	double start=0.;
	results_.clear();
	results_.push_back(x);
	double tnow = 1., step = step_size, tnorm=1.;
	if(adaptive)
		tnorm=Pot->torb(x)/10.;
	while((t_interval-start)>0.){
		if(adaptive)
			tnow = dynamic_time(results_.back());
		step=step_size*tnow/tnorm;
		O->step(start,start+step,y,step);
		for(int i=0;i<6;i++) xf[i]=y[i];
		results_.push_back(xf);
		start+=step;
	}
	return xf;
}

void Orbit::SoS(int comp, VecDoub Init, std::string file){

	// comp = 0 for x, 1 for y, 2 for z
	int ind1, ind2;
	if(comp==0){ind1=1;ind2=2;}
	if(comp==1){ind1=0;ind2=2;}
	if(comp==2){ind1=0;ind2=1;}

	results_.clear(); VecDoub xf(6,0);
	double y[6], ydown[6], yup[6], ytmp[6];
	for(int i=0;i<6;i++)y[i]=Init[i];
	double start=0., step_size_i=0.005, step_size=step_size_i;
	int direction = -1, N=0, NMAX = 10, total=0;

	while(total<200){

		for(int i=0; i<6; i++) ydown[i]=y[i];
		O->step(start,start+step_size,y,step_size);
		for(int i=0; i<6; i++){ yup[i]=y[i]; ytmp[i]=y[i];}

		if(yup[ind1]*ydown[ind1]<0.){
			while(yup[ind2]*ydown[ind2]<0. and N<NMAX){
				step_size_i/=2.;
				O->step(start,start+direction*step_size_i,y,direction*step_size_i);
				// std::cout<<ydown[ind1]<<" "<<yup[ind1]<<" "<<y[ind1]<<std::endl;
				if(y[ind1]*yup[ind1]<0.){ for(int i=0; i<6; i++) ydown[i]=y[i]; direction=1;}
				else if(y[ind1]*ydown[ind1]<0.){ for(int i=0; i<6; i++) yup[i]=y[i]; direction = -1;}
				N++;
			}
			if(N==NMAX){
				VecDoub xf = {fabs(y[comp]),fabs(y[comp+3]),y[ind1],y[ind2]};
				results_.push_back(xf);
				total++;
			}
			N=0;
			for(int i=0; i<6; i++) y[i]=ytmp[i];
		}

		start+=step_size; step_size_i=step_size; direction = -1;

	}

	std::sort(results_.begin(), results_.end(), cmp_by_first());
	std::ofstream outfile; outfile.open(file);
	if(comp==0)outfile<<"# x px\n";
	if(comp==1)outfile<<"# y py\n";
	if(comp==2)outfile<<"# z pz\n";
	for(auto i: results_) outfile<<i[0]<<" "<<i[1]<<" "<<i[2]<<" "<<i[3]<<std::endl;
	outfile.close();

}

// ============================================================================
// Symplectic orbit integration
// ============================================================================

void SymplecticOrbit::kdk2(VecDoub&x,VecDoub &a, double dt){
	for(int i=0;i<3;++i) x[i+3]+=a[i]*.5*dt;
	for(int i=0;i<3;++i) x[i]+=x[i+3]*dt;
	a = Pot->Forces(x);
	for(int i=0;i<3;++i) x[i+3]+=a[i]*.5*dt;
}

VecDoub SymplecticOrbit::integrate(const VecDoub& x, double t_interval, double step_size, bool adaptive){
	/* Integrate from time t=0 to t=t_interval at regular steps step_size */
	/* and store results in the array results							  */
	VecDoub xf=x,acc=Pot->Forces(xf);
	double start=0.,t=0.;
	results_.clear();
	results_.push_back(x);
	double dt = step_size/fac;
	double adotx=0.,xdotx=0.;
	while(start<t_interval){
		while(t<=step_size and dt!=0.){
			if(adaptive){
				adotx=0.,xdotx=0.;
				for(int i=0;i<3;i++){adotx+=-acc[i]*xf[i]; xdotx+=xf[i]*xf[i];}
				dt = sqrt(xdotx/adotx)/fac;
			}
			if(step_size-t<dt) dt=step_size-t;
			kdk2(xf,acc,dt);
			t+=dt;
		}
		if(adaptive){
			adotx=0.,xdotx=0.;
			for(int i=0;i<3;i++){adotx+=-acc[i]*xf[i]; xdotx+=xf[i]*xf[i];}
			dt = sqrt(xdotx/adotx)/fac;
		}
		else dt=step_size/fac;
		t=0.;
		start+=step_size;
		results_.push_back(xf);
	}
	return xf;
}

// ============================================================================
