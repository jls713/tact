// ============================================================================
/// \file src/adiabatic_aa.cpp
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
/// \brief Action finders for adiabatic approximations
///
/// There are two adiabatic approximations -- cylindrical and spheroidal
/// The two classes are very similar in that the actions are estimated by
/// assuming the potential is separable in either a cylindrical or spheroidal
/// coordinate system.
///
/// The cylindrical method is taken from Schoenrich \& Binney (2012).
/// The spheroidal method is from Sanders \& Binney (2015).
///
/// Both classes first construct grids of the vertical energy as a function
/// of vertical action, angular momentum and radius. The properties of the
/// grids are governed by the parameters Rmin (minimum radius), Rmax (maximum
/// radius), ZMAX (maximum z height), NGRID (number of gridpoints in R and Jz),
/// NL (number of grid points in Lz)
// ============================================================================

/*==============================================*/
/* 			Adiabatic Approximations 			*/
/*==============================================*/

// ============================================================================
#include "coordtransforms.h"
#include "adiabatic_aa.h"
#include "debug.h"
// ============================================================================
/*======================================================*/
/* 			Spheroidal Adiabatic Approximation 			*/
/*======================================================*/

static double pnu_squared(double nu, void *params){
	SpheroidalAA_zactions_struct *RS = (SpheroidalAA_zactions_struct *) params;
	return RS->ENu-RS->AA->Phi_nu({RS->lam,0.,nu},RS->Lz2);
}

static double Jnuint(double theta, void *params){
	SpheroidalAA_zactions_struct *RS = (SpheroidalAA_zactions_struct *) params;
	double nu=RS->taubar+RS->Delta*sin(theta);
	return sqrt(MAX(0.,0.5*(nu-RS->lam)/((nu+RS->AA->CS->alpha())*(MAX(TINY,nu+RS->AA->CS->gamma())))*pnu_squared(nu,RS)))*cos(theta);
}

static double plam_squared(double lam, void *params){
	SpheroidalAA_Ractions_struct *RS = (SpheroidalAA_Ractions_struct *) params;
	return RS->Elam-RS->AA->Philam_eff({lam,0.,RS->nu},RS->Lz2,RS->Jz);
}

static double Jlamint(double theta, void *params){
	SpheroidalAA_Ractions_struct *RS = (SpheroidalAA_Ractions_struct *) params;
	double lam=RS->taubar+RS->Delta*sin(theta);
	return sqrt(MAX(0.,0.5*(lam-RS->nu)/((lam+RS->AA->CS->alpha())*(lam+RS->AA->CS->gamma()))*plam_squared(lam, RS)))*cos(theta);
}

static double dJnudLzint(double theta, void *params){
	SpheroidalAA_zactions_struct *RS = (SpheroidalAA_zactions_struct *) params;
	double nu=RS->taubar+RS->Delta*sin(theta), ct=cos(theta); ct*=ct;
	double R = RS->AA->CS->tau2x({RS->lam,0.,nu})[0];
	return -.5*sqrt(RS->Lz2)*(1./(R*R)-1./(RS->lam+RS->AA->CS->alpha()))*sqrt(MAX(0.,0.5*(nu-RS->lam)*ct/((nu+RS->AA->CS->alpha())*(MAX(TINY,nu+RS->AA->CS->gamma())))/MAX(RS->tiny_number,pnu_squared(nu,RS))));
}

static double dJnudEnuint(double theta, void *params){
	SpheroidalAA_zactions_struct *RS = (SpheroidalAA_zactions_struct *) params;
	double nu=RS->taubar+RS->Delta*sin(theta), ct=cos(theta); ct*=ct;
	return .5*sqrt(MAX(0.,0.5*(nu-RS->lam)*ct/((nu+RS->AA->CS->alpha())*(MAX(TINY,nu+RS->AA->CS->gamma())))/MAX(RS->tiny_number,pnu_squared(nu,RS))));
}

static double dJlamdEint(double theta, void *params){
	SpheroidalAA_Ractions_struct *RS = (SpheroidalAA_Ractions_struct *) params;
	double lam=RS->taubar+RS->Delta*sin(theta), ct=cos(theta); ct*=ct;
	return .5*sqrt(MAX(0.,0.5*(lam-RS->nu)*ct/((lam+RS->AA->CS->alpha())*(lam+RS->AA->CS->gamma()))/MAX(RS->tiny_number,plam_squared(lam, RS))));
}

static double dJlamdLzint(double theta, void *params){
	SpheroidalAA_Ractions_struct *RS = (SpheroidalAA_Ractions_struct *) params;
	double lam=RS->taubar+RS->Delta*sin(theta), ct=cos(theta); ct*=ct;
	double R = sqrt(lam+RS->AA->CS->alpha());
	//RS->AA->CS->tau2x({lam,0.,RS->nu})[0];
	return -.5*sqrt(RS->Lz2)/(R*R)*
			sqrt(MAX(0.,0.5*(lam-RS->nu)*ct/((lam+RS->AA->CS->alpha())*(lam+RS->AA->CS->gamma()))/MAX(RS->tiny_number,plam_squared(lam, RS))));
}

static double dJlamdEnuint(double theta, void *params){
	SpheroidalAA_Ractions_struct *RS = (SpheroidalAA_Ractions_struct *) params;
	double lam=RS->taubar+RS->Delta*sin(theta), ct=cos(theta); ct*=ct;
	return -.5*sqrt(MAX(0.,0.5*(lam-RS->nu)*ct/((lam+RS->AA->CS->alpha())*(lam+RS->AA->CS->gamma()))/MAX(RS->tiny_number,plam_squared(lam, RS))));
}
static double dJlamdJzint(double theta, void *params){
	SpheroidalAA_Ractions_struct *RS = (SpheroidalAA_Ractions_struct *) params;
	double lam=RS->taubar+RS->Delta*sin(theta), ct=cos(theta); ct*=ct;
	return -(1./RS->AA->dJzdEz(RS->Lz2,sqrt(lam+RS->AA->CS->alpha()),RS->Jz))*.5*sqrt(MAX(0.,0.5*(lam-RS->nu)*ct/((lam+RS->AA->CS->alpha())*(lam+RS->AA->CS->gamma()))/MAX(RS->tiny_number,plam_squared(lam, RS))));
}
void Actions_SpheroidalAdiabaticApproximation::load_grids(std::string filename){
	std::ifstream in; in.open(filename);
	int NG,NLL;	in>>NG>>NLL; assert(NG==NGRID); assert(NLL==NL);

	Rgrid = VecDoub(NGRID,0.);
	Ezmaxgrid = std::vector<VecDoub>(NL,VecDoub(NGRID,0.));
	Ezgrid = std::vector<std::vector<VecDoub>>(NL,std::vector<VecDoub>(NGRID,VecDoub(NGRID,0.)));
	Jzgrid = std::vector<std::vector<VecDoub>>(NL,std::vector<VecDoub>(NGRID,VecDoub(NGRID,0.)));

	for(int i=0;i<NGRID;i++) in>>Rgrid[i];
	for(int n=0;n<NL;n++)
		for(int i=0;i<NGRID;i++) in>>Ezmaxgrid[n][i];
	for(int n=0;n<NL;n++)
		for(int i=0;i<NGRID;i++) for(int j=0;j<NGRID;j++) in>>Ezgrid[n][i][j];
	for(int n=0;n<NL;n++)
		for(int i=0;i<NGRID;i++) for(int j=0;j<NGRID;j++) in>>Jzgrid[n][i][j];
	in.close();
}

void Actions_SpheroidalAdiabaticApproximation::make_grids(std::string filename){
	for(int n=0; n<NGRID; n++)
		Rgrid.push_back(.75*Rmin+n*(1.5*Rmax-.75*Rmin)/((double)(NGRID-1)));

	for(int n=0; n<NL; n++)
		Lgrid.push_back(Lmin+(double)n*(Lmax-Lmin)/(double)(NL-1));

	Ezmaxgrid = std::vector<VecDoub>(NL,VecDoub(NGRID,0.));
	Ezgrid = std::vector<std::vector<VecDoub>>(NL,std::vector<VecDoub>(NGRID,VecDoub(NGRID,0.)));
	Jzgrid = std::vector<std::vector<VecDoub>>(NL,std::vector<VecDoub>(NGRID,VecDoub(NGRID,0.)));

	for(int nL=0; nL<NL; nL++){
		double Lz2 = Lgrid[nL];Lz2*=Lz2;
		for(int nR=0; nR<NGRID; nR++){
			VecDoub xx = {Rgrid[nR],0.,.1};
			double a = -Pot->DeltaGuess(xx)+CS->gamma();
			CS->newalpha(a>CS->gamma()?CS->gamma()-0.1:a);
			double lam = Rgrid[nR]*Rgrid[nR]-CS->alpha();
			double zmax = MIN(ZMAX,sqrt(lam+CS->gamma())*0.9);
			double Rmax = sqrt(MAX((1.0-zmax*zmax/(Rgrid[nR]*Rgrid[nR]-CS->alpha()+CS->gamma()))*Rgrid[nR]*Rgrid[nR],0.0));
			if(Rmax==0.0){Rmax = 0.5; zmax = sqrt((1.0-Rmax*Rmax/(Rgrid[nR]*Rgrid[nR]))*(Rgrid[nR]*Rgrid[nR]-CS->alpha()+CS->gamma()));}
			Ezmaxgrid[nL][nR] = Pot->Phi({Rmax,0.,zmax})-Pot->Phi({Rgrid[nR],0.,0.})+.5*Lz2*(1./(Rmax*Rmax)-1./(Rgrid[nR]*Rgrid[nR]));
			Ezgrid[nL][nR][0]=Jzgrid[nL][nR][0]=0;
			double nulast=-CS->gamma()+SMALL, nulim = 0.;

			for(int i=1; i<NGRID; i++){
				Ezgrid[nL][nR][i]=((double)i)*sqrt(Ezmaxgrid[nL][nR])/(double)(NGRID-1);
				Ezgrid[nL][nR][i]*=Ezgrid[nL][nR][i];
				Jzgrid[nL][nR][i]=actions_Jz(Ezgrid[nL][nR][i],Lz2,{Rgrid[nR]*Rgrid[nR]-CS->alpha(),0.,nulast},&nulim);
				nulast = nulim;
				// std::cout<<CS->alpha()<<" "<<Rgrid[nR]<<" "<<nulim<<" "<<Lz2<<" "<<Ezgrid[nL][nR][i]<<" "<<Jzgrid[nL][nR][i]<<std::endl;
			}
		}
	}
	if(filename!=""){
		std::ofstream out; out.open(filename);
		out<<NGRID<<" "<<NL<<std::endl;
		for(auto i:Rgrid)out<<i<<" ";out<<std::endl;
		for(auto i:Ezmaxgrid){for(auto j:i)out<<j<<" ";out<<std::endl;}
		for(auto i:Ezgrid) for(auto j:i){ for(auto k:j) out<<k<<" "; out<<std::endl;}
		for(auto i:Jzgrid) for(auto j:i){ for(auto k:j) out<<k<<" "; out<<std::endl;}
		out.close();
	}
}

double Actions_SpheroidalAdiabaticApproximation::Ez_from_grid(double Lz2,double R,double Jz){
	//Ez(R,Jz) by interpolation
	if(no_energy_correction)
		return 0;
	else{
		int botR,topR;
		if(R<Rgrid[0]) R = Rgrid[0];
		if(R>Rgrid[NGRID-1]) R = Rgrid[NGRID-1];
		topbottom<double>(Rgrid,R,&botR,&topR,"SAA Ez_from_grid R");

		double Lz = sqrt(Lz2);
		int botL,topL;
		if(Lz<Lgrid[0]) Lz = Lgrid[0];
		if(Lz>Lgrid[NL-1]) Lz = Lgrid[NL-1];
		topbottom<double>(Lgrid,Lz,&botL,&topL,"SAA Ez_from_grid Lz");

		double fac=(R-Rgrid[botR])/(Rgrid[topR]-Rgrid[botR]);
		double fac2=(Lz-Lgrid[botL])/(Lgrid[topL]-Lgrid[botL]);

		// std::cout<<Jz<<std::endl;
		double Ezb=linterp<double>(Jzgrid[botL][botR],Ezgrid[botL][botR],MIN(Jz,Jzgrid[botL][botR][NGRID-1]));
		double Ezt=linterp(Jzgrid[botL][topR],Ezgrid[botL][topR],MIN(Jz,Jzgrid[botL][topR][NGRID-1]));

		double Ezb2=linterp<double>(Jzgrid[topL][botR],Ezgrid[topL][botR],MIN(Jz,Jzgrid[topL][botR][NGRID-1]));
		double Ezt2=linterp(Jzgrid[topL][topR],Ezgrid[topL][topR],MIN(Jz,Jzgrid[topL][topR][NGRID-1]));

		return (1-fac2)*((1-fac)*Ezb+fac*Ezt)+fac2*((1-fac)*Ezb2+fac*Ezt2);
	}
}


double Actions_SpheroidalAdiabaticApproximation::Jz_from_grid(double Lz2,double R,double Ez){
	//Jz(R,Ez) by interpolation
	if(no_energy_correction)
		return 0;
	else{
		int botR,topR;
		if(R<Rgrid[0]) R = Rgrid[0];
		if(R>Rgrid[NGRID-1]) R = Rgrid[NGRID-1];
		topbottom<double>(Rgrid,R,&botR,&topR,"SAA Jz_from_grid R");

		double Lz = sqrt(Lz2);
		int botL,topL;
		if(Lz<Lgrid[0]) Lz = Lgrid[0];
		if(Lz>Lgrid[NL-1]) Lz = Lgrid[NL-1];
		topbottom<double>(Lgrid,Lz,&botL,&topL,"SAA Jz_from_grid Lz");

		double fac=(R-Rgrid[botR])/(Rgrid[topR]-Rgrid[botR]);
		double fac2=(Lz-Lgrid[botL])/(Lgrid[topL]-Lgrid[botL]);

		double Jzb=linterp<double>(Ezgrid[botL][botR],Jzgrid[botL][botR],MIN(Ez,Ezgrid[botL][botR][NGRID-1]));
		double Jzt=linterp(Ezgrid[botL][topR],Jzgrid[botL][topR],MIN(Ez,Ezgrid[botL][topR][NGRID-1]));

		double Jzb2=linterp<double>(Ezgrid[topL][botR],Jzgrid[topL][botR],MIN(Ez,Ezgrid[topL][botR][NGRID-1]));
		double Jzt2=linterp(Ezgrid[topL][topR],Jzgrid[topL][topR],MIN(Ez,Ezgrid[topL][topR][NGRID-1]));

		return (1-fac2)*((1-fac)*Jzb+fac*Jzt)+fac2*((1-fac)*Jzb2+fac*Jzt2);
	}
}

double Actions_SpheroidalAdiabaticApproximation::dJzdEz(double Lz2,double R,double Jz){
	// dJzdEz at fixed lam and Lz
	double dJ = 0.1*Jz;
	double Eu = Ez_from_grid(Lz2,R,Jz+dJ);
	double Ed = Ez_from_grid(Lz2,R,Jz-dJ);
	return 2.*dJ/(Eu-Ed);
}

double Actions_SpheroidalAdiabaticApproximation::dJzdR(double Lz2,double R,double Ez){
	// dJzdR at fixed Ez and Lz
	double dR = 0.1*R;
	double Ju = Jz_from_grid(Lz2,R+dR,Ez);
	double Jd = Jz_from_grid(Lz2,R-dR,Ez);
	return (Ju-Jd)/2./dR;
}

double Actions_SpheroidalAdiabaticApproximation::estimate_tiny(double ENu, double lam, double nu, double Lz2){
	return 2.*(ENu-Phi_nu({lam,0.,nu},Lz2))/1.e8;
}

double Actions_SpheroidalAdiabaticApproximation::find_nulimit(double ENu, double Lz2, VecDoub tau){
	SpheroidalAA_zactions_struct RS(this,ENu,Lz2,tau[0],0.,0.,0.);
	double nu=tau[2],nutry = tau[2];double tiny_number=TINY;
	root_find RF(TINY,100);
	if(pnu_squared(nu, &RS)>=0.0){
		while(pnu_squared(nutry, &RS)>=0.
		      and -(nutry+CS->alpha())>tiny_number)
			nutry+=0.1*(-CS->alpha()-nutry);
		if(-(nutry+CS->alpha())>tiny_number)
			return RF.findroot(&pnu_squared,nu,nutry,&RS);
		else return -CS->alpha();
	}
	else return tau[2]+tiny_number;
}

VecDoub Actions_SpheroidalAdiabaticApproximation::find_lamlimits(double Elam, double Lz2,double Jz,VecDoub tau){

	SpheroidalAA_Ractions_struct RS(this,Elam,Lz2,Jz,tau[2],0.,0.,0.);
	// find roots of p^2(lam)
	VecDoub limits;
	double lambda=tau[0], tiny_number=TINY;
	root_find RF(TINY,100);
	double laminner=lambda, lamouter=lambda;
	if(plam_squared(lambda, &RS)>=0.0){
		while(plam_squared(laminner, &RS)>=0.0
		      and (laminner+CS->alpha())>tiny_number)
			laminner-=.1*(laminner+CS->alpha());
		while(plam_squared(lamouter, &RS)>=0.)	lamouter*=1.1;

		if((laminner+CS->alpha())>tiny_number)
			limits.push_back(
			    RF.findroot(&plam_squared,laminner,lambda,&RS));
		else limits.push_back(-CS->alpha());
		limits.push_back(RF.findroot(&plam_squared,lambda,lamouter,&RS));
		if(fabs(limits[0]-limits[1])<TINY){
		    if(plam_squared(lambda-5*TINY, &RS)>0.0)
			limits[0]=RF.findroot(&plam_squared,laminner,lambda-5*TINY,&RS);
	            else
			limits[1]=RF.findroot(&plam_squared,lambda+5*TINY,lamouter,&RS);
		}
	}
	else{
		limits.push_back(lambda-tiny_number);
		limits.push_back(lambda+tiny_number);
	}
	return limits;
}

double Actions_SpheroidalAdiabaticApproximation::actions_Jz(double ENu, double Lz2, VecDoub tau, double *nulim){
	// Calculates vertical action -- only require R and Ez but for speed reasons it is efficient to pass the input z value and output the found z limit
	*nulim = find_nulimit(ENu,Lz2,tau);
	double Delta = .5*(*nulim+CS->gamma()), taubar = .5*(*nulim-CS->gamma());
	SpheroidalAA_zactions_struct PA(this,ENu,Lz2,tau[0],Delta,taubar,0.);
	return 2.*Delta*GaussLegendreQuad(&Jnuint,-.5*PI,.5*PI,&PA)/PI;

}

VecDoub Actions_SpheroidalAdiabaticApproximation::actions(const VecDoub& x, void*params){

	VecDoub acts(3,0.);
    if(action_check(x,acts,Pot)) return acts;
	if(fabs(x[2])>ZMAX)
		std::cerr<<"z outside grid range: Pass larger ZMAX to constructor"<<std::endl;
	if(params){
		double *deltaguess = (double *) params;
		if(*deltaguess>0.) CS->newalpha(-Pot->DeltaGuess(x)+CS->gamma());
		else CS->newalpha(*deltaguess);
	}
	else{CS->newalpha(-Pot->DeltaGuess(x)+CS->gamma());}

	if(CS->alpha()>CS->gamma()){
	    if(debug_NegativeDelta)
	    	std::cerr<<"Negative Delta at R="<<sqrt(x[0]*x[0]+x[1]*x[1])<<", z="<<x[2]<<std::endl;
        CS->newalpha(CS->gamma()-0.1);
	}
	VecDoub tau = CS->xv2tau(x);
	// If exactly on the plane we integrate a bit
	if(x[2]==0)
		tau = CS->xv2tau(integrate_a_bit(x,Pot));
	double Lz = Pot->Lz(x), Lz2 = Lz*Lz;
	// if vz is small use Enu \approx Ez = 1/2 vz^2 + Pot(R,z)-Pot(R,0)
	double ENu = (fabs(tau[5])>TINY?(tau[2]-tau[0])/(8.*(tau[2]+CS->alpha())
	               *(tau[2]+CS->gamma()))*tau[5]*tau[5]:.5*x[5]*x[5]+Pot->Phi(x)-Pot->Phi({x[0],0.,0.}))
				 +Phi_nu(tau,Lz2);

	double nulim = tau[2];
	double Jz = actions_Jz(ENu, Lz2, tau, &nulim);
	double Elam = (tau[0]-tau[2])*tau[3]*tau[3]/(8.0*(tau[0]+CS->alpha())*(tau[0]+CS->gamma()))+Philam_eff(tau,Lz2,Jz);
	VecDoub lamlims = find_lamlimits(Elam,Lz2,Jz,tau);
	double Delta = .5*(lamlims[1]-lamlims[0]), taubar = .5*(lamlims[1]+lamlims[0]);
	SpheroidalAA_Ractions_struct RA(this,Elam,Lz2,Jz,-CS->gamma(),Delta,taubar, 0.);
	double JR = Delta*GaussLegendreQuad(&Jlamint,-.5*PI,.5*PI,&RA)/PI;

	return {JR,Lz,Jz};
}

VecDoub Actions_SpheroidalAdiabaticApproximation::angles(const VecDoub& x, void *params){

	VecDoub angs(6,0.);
    if(angle_check(x,angs,Pot)) return angs;

	if(fabs(x[2])>ZMAX) std::cerr<<"z outside grid range"<<std::endl;

	if(params){
		double *deltaguess = (double *) params;
		if(*deltaguess>0.) CS->newalpha(-Pot->DeltaGuess(x)+CS->gamma());
		else CS->newalpha(*deltaguess);
	}
	else{CS->newalpha(-Pot->DeltaGuess(x)+CS->gamma());}

	if(CS->alpha()>CS->gamma()){
	    if(debug_NegativeDelta)
	    	std::cerr<<"Negative Delta at R="<<sqrt(x[0]*x[0]+x[1]*x[1])<<", z="<<x[2]<<std::endl;
        CS->newalpha(CS->gamma()-0.1);
	}

	VecDoub angles(6,0.);

	VecDoub tau = CS->xv2tau(x);
	// If exactly on the plane we integrate a bit
	if(x[2]==0.)
		tau = CS->xv2tau(integrate_a_bit(x,Pot));
	double Lz = Pot->Lz(x), Lz2 = Lz*Lz;
	// if vz is small use Enu \approx Ez = 1/2 vz^2 + Pot(R,z)-Pot(R,0)
	double ENu = (fabs(tau[5])>TINY?(tau[2]-tau[0])/(8.*(tau[2]+CS->alpha())
	               *(tau[2]+CS->gamma()))*tau[5]*tau[5]:.5*x[5]*x[5]+Pot->Phi(x)-Pot->Phi({x[0],0.,0.}))
				 +Phi_nu(tau,Lz2);

	double tn = estimate_tiny(ENu, tau[0], tau[2], Lz2);
	std::cout<<tn<<" "<<ENu<<std::endl;
	if(tn==0.)tn=1e-15;
	double nulim = tau[2];
	double Jz = actions_Jz(ENu, Lz2, tau, &nulim);
	double Delta = .5*(nulim+CS->gamma()), taubar = .5*(nulim-CS->gamma());
	SpheroidalAA_zactions_struct PA(this,ENu,Lz2,tau[0],Delta,taubar, tn);
	VecDoub dJdIn(3,0.), dJdIl(3,0.), dSdI(3,0.);
	dSdI[1]+=SIGN(tau[4])*tau[1];

	double anglim = asin((tau[2]-taubar)/Delta);
	double lsign = SIGN(tau[5]);
	dJdIn[0] = 0.;//2.*Delta*GaussLegendreQuad(&dJnudEint,-.5*PI,.5*PI,&PA)/PI;
	dJdIn[1] = 2.*Delta*GaussLegendreQuad(&dJnudLzint,-.5*PI,.5*PI,&PA)/PI;
	dJdIn[2] = 1.;//2.*Delta*GaussLegendreQuad(&dJnudEnuint,-.5*PI,.5*PI,&PA)/PI;
	dSdI[0] += 0.;//lsign*Delta*GaussLegendreQuad(&dJnudEint,-.5*PI,anglim,&PA)/PI;
	dSdI[1] += lsign*Delta*GaussLegendreQuad(&dJnudLzint,-.5*PI,anglim,&PA);
	dSdI[2] += lsign*Delta*GaussLegendreQuad(&dJnudEnuint,-.5*PI,anglim,&PA)/(2.*Delta*GaussLegendreQuad(&dJnudEnuint,-.5*PI,.5*PI,&PA)/PI);

	// /dJzdEz(Lz2,sqrt(tau[0]+CS->alpha()),Jz);

	double Elam = (tau[0]-tau[2])*tau[3]*tau[3]/(8.0*(tau[0]+CS->alpha())*(tau[0]+CS->gamma()))+Philam_eff(tau,Lz2,Jz);
	VecDoub lamlims = find_lamlimits(Elam,Lz2,Jz,tau);

	Delta = .5*(lamlims[1]-lamlims[0]);
	taubar = .5*(lamlims[1]+lamlims[0]);
	SpheroidalAA_Ractions_struct RA(this,Elam,Lz2,Jz,-CS->gamma(),Delta,taubar, tn);

	anglim = asin((tau[0]-taubar)/Delta);
	lsign = SIGN(tau[3]);
	dJdIl[0] = Delta*GaussLegendreQuad(&dJlamdEint,-.5*PI,.5*PI,&RA)/PI;
	dJdIl[1] = Delta*GaussLegendreQuad(&dJlamdLzint,-.5*PI,.5*PI,&RA)/PI;
	dJdIl[2] = Delta*GaussLegendreQuad(&dJlamdJzint,-.5*PI,.5*PI,&RA)/PI;
	dSdI[0] += lsign*Delta*GaussLegendreQuad(&dJlamdEint,-.5*PI,anglim,&RA);
	dSdI[1] += lsign*Delta*GaussLegendreQuad(&dJlamdLzint,-.5*PI,anglim,&RA);
	dSdI[2] += lsign*Delta*GaussLegendreQuad(&dJlamdJzint,-.5*PI,anglim,&RA);

	VecDoub dJdIp = {0.,1.,0.};

	double Determinant = dJdIp[1]*(dJdIl[0]*dJdIn[2]-dJdIl[2]*dJdIn[0]);

	double dIdJ[3][3];
	dIdJ[0][0] = det2(dJdIp[1],dJdIp[2],dJdIn[1],dJdIn[2])/Determinant;
	dIdJ[0][1] = -det2(dJdIl[1],dJdIl[2],dJdIn[1],dJdIn[2])/Determinant;
	dIdJ[0][2] = det2(dJdIl[1],dJdIl[2],dJdIp[1],dJdIp[2])/Determinant;
	dIdJ[1][0] = -det2(dJdIp[0],dJdIp[2],dJdIn[0],dJdIn[2])/Determinant;
	dIdJ[1][1] = det2(dJdIl[0],dJdIl[2],dJdIn[0],dJdIn[2])/Determinant;
	dIdJ[1][2] = -det2(dJdIl[0],dJdIl[2],dJdIp[0],dJdIp[2])/Determinant;
	dIdJ[2][0] = det2(dJdIp[0],dJdIp[1],dJdIn[0],dJdIn[1])/Determinant;
	dIdJ[2][1] = -det2(dJdIl[0],dJdIl[1],dJdIn[0],dJdIn[1])/Determinant;
	dIdJ[2][2] = det2(dJdIl[0],dJdIl[1],dJdIp[0],dJdIp[1])/Determinant;

	angles[0]=dSdI[0]*dIdJ[0][0]+dSdI[1]*dIdJ[1][0]+dSdI[2]*dIdJ[2][0];
	angles[1]=dSdI[0]*dIdJ[0][1]+dSdI[1]*dIdJ[1][1]+dSdI[2]*dIdJ[2][1];
	angles[2]=dSdI[0]*dIdJ[0][2]+dSdI[1]*dIdJ[1][2]+dSdI[2]*dIdJ[2][2];
	angles[3]=det2(dJdIp[1],dJdIp[2],dJdIn[1],dJdIn[2])/Determinant;
	angles[4]=det2(dJdIn[1],dJdIn[2],dJdIl[1],dJdIl[2])/Determinant;
	angles[5]=det2(dJdIl[1],dJdIl[2],dJdIp[1],dJdIp[2])/Determinant;

	angles[4]*=SIGN(Lz);

	if(tau[5]<0.0){angles[2]+=PI;}
	if(x[2]<0.0){angles[2]+=PI;}
	if(angles[2]>2.*PI)	angles[2]-=2.0*PI;
	if(angles[0]<0.0)	angles[0]+=2.0*PI;
	if(angles[1]<0.0)	angles[1]+=2.0*PI;
	if(angles[2]<0.0)	angles[2]+=2.0*PI;

	return angles;

}
// double J_integrand_AxiSym(double theta, void *params){
// 	// need to set taubargl and Deltagl
// 	action_struct_axi * AS = (action_struct_axi *) params;
// 	double tau=AS->taubargl+AS->Deltagl*sin(theta);
// 	return sqrt(MAX(0.,ptau2_AxiSym(tau,&(AS->RS))))*cos(theta);
// }

// double PhiNu(double nu, action_struct_aa *RS){
// 	double Alpha = RS->AA->CS->alpha(), Gamma = RS->AA->CS->gamma();
// 	double Lz2 = RS->Lz*RS->Lz;
// 	double R = sqrt((RS->lamgl+Alpha)*(nu+Alpha)/(Alpha-Gamma));
// 	double z = sqrt((RS->lamgl+Gamma)*(nu+Gamma)/(Gamma-Alpha));
// 	return RS->AA->Pot->Phi({R,0.,z})+Lz2*(1./(2.0*R*R)
// 	        -1./(2.0*(RS->lamgl+RS->AA->CS->alpha())));
// }

// double nuturnfn(double nu, void *params){
// 	action_struct_aa * AS = (action_struct_aa *) params;
// 	return AS->ENugl-PhiNu(nu,AS);
// }

// double JNuint(double theta, void *params){
// 	action_struct_aa * AS = (action_struct_aa *) params;
// 	double Alpha = AS->AA->CS->alpha(), Gamma = AS->AA->CS->gamma();
// 	double nu=AS->taubargl+AS->Deltagl*sin(theta);
// 	return sqrt(MAX(0.,0.5*(nu-AS->lamgl)/((nu+Alpha)*(nu+Gamma))
// 			*(AS->ENugl-PhiNu(nu,AS))))*cos(theta);
// }

// VecDoub Actions_EllipsoidalAdiabaticApproximation::find_limits(VecDoub tau, action_struct_aa *RS){

// 	double lambda = tau[0], nu = tau[2];
// 	VecDoub limits;
// 	limits.push_back(-CS->gamma()+TINY);
// 	// find root of p^2(nu)
// 	double nuouter=nu;
// 	while(nuturnfn(nuouter, &RS)<0.)	nuouter+=0.1*(-CS->alpha()-nuouter);
// 	root_find RF(SMALL,100);
// 	limits.push_back(RF.findroot(&nuturnfn,nu,nuouter,&RS));
// 	return limits;
// }

// VecDoub Actions_EllipsoidalAdiabaticApproximation::actions(VecDoub x){
// 	CS->newalpha(Pot->DeltaGuess(x)+CS->gamma());
// 	double Lz = fabs(Pot->Lz(x));
// 	VecDoub tau = CS->xv2tau(x);
// 	double ENugl = (tau[2]-tau[0])/(4.*(tau[2]+CS->alpha())
// 	               *(tau[2]+CS->gamma()))*tau[5]*tau[5];
// 	action_struct_aa RS(this,ENugl,tau[0],Lz,0.,0.);
// 	ENugl+=PhiNu(tau[2],&RS);
// 	// Jnu
// }

// ============================================================================
/*======================================================*/
/* 		 Cylindrical Adiabatic Approximation 			*/
/*======================================================*/
static double pz_squared(double z, void *params){
	PolarAA_zactions_struct *RS = (PolarAA_zactions_struct *) params;
	return RS->Ez-RS->AA->Phi_z({RS->R,0.,z});
}

static double Jzint(double theta, void *params){
	PolarAA_zactions_struct *RS = (PolarAA_zactions_struct *) params;
	double z=RS->zlim*sin(theta);
	return sqrt(MAX(0.,2.*pz_squared(z,RS)))*cos(theta);
}

static double pR_squared(double R, void *params){
	PolarAA_Ractions_struct *RS = (PolarAA_Ractions_struct *) params;
	return RS->Etot-RS->AA->PhiR_eff(R,RS->Lz2,RS->Jz);
}

static double JRint(double theta, void *params){
	PolarAA_Ractions_struct *RS = (PolarAA_Ractions_struct *) params;
	double R=RS->taubar+RS->Delta*sin(theta);
	return sqrt(MAX(0.,2.*pR_squared(R, RS)))*cos(theta);
}

static double dJzdEzint(double theta, void *params){
	PolarAA_zactions_struct *RS = (PolarAA_zactions_struct *) params;
	double z=RS->zlim*sin(theta), ct = cos(theta); ct*=ct;
	return sqrt(MAX(0.,ct/(2.*pz_squared(z,RS))));
}

static double dJRdEint(double theta, void *params){
	PolarAA_Ractions_struct *RS = (PolarAA_Ractions_struct *) params;
	double R=RS->taubar+RS->Delta*sin(theta), ct = cos(theta); ct*=ct;
	return sqrt(MAX(0.,ct/(2.*pR_squared(R,RS))));
}

static double dJRdLzint(double theta, void *params){
	PolarAA_Ractions_struct *RS = (PolarAA_Ractions_struct *) params;
	double R=RS->taubar+RS->Delta*sin(theta), ct = cos(theta); ct*=ct;
	return -sqrt(MAX(0.,RS->Lz2*ct/(2.*pR_squared(R,RS))))/R/R;
}

static double dJRdJzint(double theta, void *params){
	PolarAA_Ractions_struct *RS = (PolarAA_Ractions_struct *) params;
	double R=RS->taubar+RS->Delta*sin(theta), ct = cos(theta); ct*=ct;
	return -(1./RS->AA->dJzdEz(R,RS->Jz))*sqrt(MAX(0.,ct/(2.*pR_squared(R,RS))));
}

static double dJRdEzint(double theta, void *params){
	PolarAA_Ractions_struct *RS = (PolarAA_Ractions_struct *) params;
	double R=RS->taubar+RS->Delta*sin(theta), ct = cos(theta); ct*=ct;
	return -sqrt(ct/(2.*pR_squared(R,RS)));
}
void Actions_CylindricalAdiabaticApproximation::load_grids(std::string filename){

	std::ifstream in; in.open(filename);
	int NG;	in>>NG; assert(NG==NGRID);

	Rgrid = VecDoub(NGRID,0.); Ezmaxgrid = VecDoub(NGRID,0.);
	Ezgrid = std::vector<VecDoub>(NGRID,VecDoub(NGRID,0.));
	Jzgrid= std::vector<VecDoub>(NGRID,VecDoub(NGRID,0.));

	for(int i=0;i<NGRID;i++) in>>Rgrid[i];
	for(int i=0;i<NGRID;i++) in>>Ezmaxgrid[i];
	for(int i=0;i<NGRID;i++) for(int j=0;j<NGRID;j++) in>>Ezgrid[i][j];
	for(int i=0;i<NGRID;i++) for(int j=0;j<NGRID;j++) in>>Jzgrid[i][j];
	in.close();
}

void Actions_CylindricalAdiabaticApproximation::make_grids(std::string filename){

	for(int n=0; n<NGRID; n++){
		Rgrid.push_back(.75*Rmin+n*(1.5*Rmax-.75*Rmin)/((double)(NGRID-1)));
		Ezmaxgrid.push_back(Pot->Phi({Rgrid[n],0.,ZMAX})-Pot->Phi({Rgrid[n],0.,0.}));
		Ezgrid.push_back(VecDoub(NGRID,0.));
		Jzgrid.push_back(VecDoub(NGRID,0.));
		dJzgrid.push_back(VecDoub(NGRID,0.));
		double zlast=.005, zlim = 0.;

		for(int i=1; i<NGRID; i++){
			Ezgrid[n][i]=((double)i)*sqrt(Ezmaxgrid[n])/(double)(NGRID-1);
			Ezgrid[n][i]*=Ezgrid[n][i];
			Jzgrid[n][i]=actions_Jz(Rgrid[n],Ezgrid[n][i],zlast,&zlim);
			dJzgrid[n][i]=actions_dJzdEz(Rgrid[n],Ezgrid[n][i],zlast,&zlim);
			zlast = zlim;
		}
	}
	if(filename!=""){
		std::ofstream out; out.open(filename);
		out<<NGRID<<std::endl;
		for(auto i:Rgrid)out<<i<<" ";out<<std::endl;
		for(auto i:Ezmaxgrid)out<<i<<" ";out<<std::endl;
		for(auto i:Ezgrid){ for(auto j:i) out<<j<<" "; }out<<std::endl;
		for(auto i:Jzgrid){ for(auto j:i) out<<j<<" "; }out<<std::endl;
		out.close();
	}
}

double Actions_CylindricalAdiabaticApproximation::Ez_from_grid(double R,double Jz){
	//Ez(R,Jz) by interpolation
	if(no_energy_correction)
		return 0;
	else{
		int botR,topR;
		if(R<Rgrid[0]){ R = Rgrid[0]; std::cerr<<"Outside R grid in CAA:"<<R<<". Pass smaller Rmin to constructor.\n";}
		if(R>Rgrid[NGRID-1]){ R = Rgrid[NGRID-1]; std::cerr<<"Outside R grid in CAA:"<<R<<". Pass larger Rmax to constructor.\n";}
		topbottom<double>(Rgrid,R,&botR,&topR);
		double fac=(R-Rgrid[botR])/(Rgrid[topR]-Rgrid[botR]);
		double Ezb=linterp<double>(Jzgrid[botR],Ezgrid[botR],MIN(Jz,Jzgrid[botR][NGRID-1]));
		double Ezt=linterp(Jzgrid[topR],Ezgrid[topR],MIN(Jz,Jzgrid[topR][NGRID-1]));
		return (1-fac)*Ezb+fac*Ezt;
	}
}


double Actions_CylindricalAdiabaticApproximation::Jz_from_grid(double R,double Ez){
	//Jz(R,Ez) by interpolation
	if(no_energy_correction)
		return 0;
	else{
		int botR,topR;
		if(R<Rgrid[0]){ R = Rgrid[0]; std::cerr<<"Outside R grid in CAA:"<<R<<". Pass smaller Rmin to constructor.\n";}
		if(R>Rgrid[NGRID-1]){ R = Rgrid[NGRID-1]; std::cerr<<"Outside R grid in CAA:"<<R<<". Pass larger Rmax to constructor.\n";}
		topbottom<double>(Rgrid,R,&botR,&topR);
		double fac=(R-Rgrid[botR])/(Rgrid[topR]-Rgrid[botR]);
		double Jzb=linterp<double>(Ezgrid[botR],Jzgrid[botR],MIN(Ez,Ezgrid[botR][NGRID-1]));
		double Jzt=linterp(Ezgrid[topR],Jzgrid[topR],MIN(Ez,Ezgrid[topR][NGRID-1]));
		return (1-fac)*Jzb+fac*Jzt;
	}
}

double Actions_CylindricalAdiabaticApproximation::dJzdEz(double R, double Jz){
	// dJzdEz at fixed R0
	double dJ = 0.1*Jz;
	double Eu = Ez_from_grid(R,Jz+dJ);
	double Ed = Ez_from_grid(R,Jz-dJ);
	return 2.*dJ/(Eu-Ed);

}

double Actions_CylindricalAdiabaticApproximation::dJzdR(double R, double Ez){
	// dJzdR at fixed Ez
	double dR = 0.1*R;
	double Ju = Jz_from_grid(R+dR,Ez);
	double Jd = Jz_from_grid(R-dR,Ez);
	return (Ju-Jd)/2./dR;
}

double Actions_CylindricalAdiabaticApproximation::find_zlimit(double Ez, VecDoub Polar){
	PolarAA_zactions_struct RS(this,Ez,Polar[0],0.,0.);
	double ztry = fabs(Polar[2]);
	if(ztry==0.)ztry+=SMALL;
	if(pz_squared(ztry, &RS)<0.)return ztry;
	while(pz_squared(ztry, &RS)>0.) ztry*=1.1;
	root_find RF(TINY,100);
	return RF.findroot(&pz_squared,fabs(Polar[2]),ztry,&RS);
}

VecDoub Actions_CylindricalAdiabaticApproximation::find_Rlimits(double R, double Etot, double Lz2,double Jz){

	PolarAA_Ractions_struct RS(this,Etot,Lz2,Jz,0.,0.,0.);
	double Rin = R, Rout = R;
	while(pR_squared(Rout, &RS)>=0.) Rout*=1.01;
	while(pR_squared(Rin, &RS)>=0.)  Rin*=.99;
	root_find RF(TINY,100);
	return {Lz2>0.?RF.findroot(&pR_squared,Rin,Rin/.99,&RS):0.,
				RF.findroot(&pR_squared,Rout/1.01,Rout,&RS)};
}

double Actions_CylindricalAdiabaticApproximation::estimate_tiny(double Ez, double R){
	return 2.*(Ez-Phi_z({R,0.,0.}))/1.e5;
}

double Actions_CylindricalAdiabaticApproximation::actions_Jz(double R, double Ez,double z, double *zlim){
	// Calculates vertical action -- only require R and Ez but for speed reasons it is efficient to pass the input z value and output the found z limit
	*zlim = find_zlimit(Ez, {R,0.,z});
	PolarAA_zactions_struct PA(this,Ez,R,*zlim,0.);
	return 2.*(*zlim)*GaussLegendreQuad(&Jzint,0.,.5*PI,&PA)/PI;

}

double Actions_CylindricalAdiabaticApproximation::actions_dJzdEz(double R, double Ez,double z, double *zlim){
	PolarAA_zactions_struct PA(this,Ez,R,*zlim,0.);
	return 2.*(*zlim)*GaussLegendreQuad(&dJzdEzint,0.,.5*PI,&PA)/PI;

}

VecDoub Actions_CylindricalAdiabaticApproximation::actions(const VecDoub& x, void*params){
	VecDoub acts(3,0.);
    if(action_check(x,acts,Pot)) return acts;
	VecDoub Polar = conv::CartesianToPolar(x);

	double vz = Polar[5], R = Polar[0], z = Polar[2];
	double Ez = .5*vz*vz+Phi_z(Polar);
	double zlim = z;
	double Jz = actions_Jz(R, Ez, Polar[2], &zlim);
	double Lz = Pot->Lz(x), Lz2 = Lz*Lz;
	double vR = Polar[3];
	double Etot = .5*vR*vR+PhiR_eff(R, Lz2, Jz);
	VecDoub Rlims = find_Rlimits(R, Etot, Lz2,Jz);
	double Delta = .5*(Rlims[1]-Rlims[0]), taubar = .5*(Rlims[1]+Rlims[0]);
	PolarAA_Ractions_struct RA(this,Etot,Lz2,Jz,Delta,taubar,0.);
	double JR = Delta*GaussLegendreQuad(&JRint,-.5*PI,.5*PI,&RA)/PI;

	return {JR,Lz,Jz};
}


VecDoub Actions_CylindricalAdiabaticApproximation::angles(const VecDoub& x, void *params){

    VecDoub angs(6,0.);
    if(angle_check(x,angs,Pot)) return angs;

	VecDoub Polar = conv::CartesianToPolar(x);
	VecDoub angles(6,0.);

	double vz = Polar[5], R = Polar[0], z = Polar[2];
	double Ez = .5*vz*vz+Phi_z(Polar);
	double zlim;
	double Jz = actions_Jz(R, Ez, z, &zlim);
	double Delta = .5*zlim, taubar = .5*zlim;
	double tn = estimate_tiny(Ez,R);
	PolarAA_zactions_struct PA(this,Ez,R,zlim,tn);
	VecDoub dJdIn(3,0.), dJdIl(3,0.), dSdI(3,0.);
	dSdI[1]+=SIGN(Polar[4])*Polar[1];

	double anglim = asin((fabs(z)-taubar)/Delta);
	double lsign = SIGN(x[5]*x[2]);
	dJdIn[0] = 0.;
	dJdIn[1] = 0.;
	dJdIn[2] = 1.;//2.*Delta*GaussLegendreQuad(&dJzdEzint,-.5*PI,.5*PI,&PA)/PI;
	dSdI[0] += 0.;
	dSdI[1] += 0.;
	dSdI[2] += lsign*Delta*GaussLegendreQuad(&dJzdEzint,-.5*PI,anglim,&PA)/(2.*Delta*GaussLegendreQuad(&dJzdEzint,-.5*PI,.5*PI,&PA)/PI);

	double Lz = Pot->Lz(x), Lz2 = Lz*Lz;
	double vR = Polar[3];
	double Etot = .5*vR*vR+PhiR_eff(R, Lz2, Jz);
	VecDoub Rlims = find_Rlimits(R, Etot, Lz2,Jz);
	Delta = .5*(Rlims[1]-Rlims[0]);
	taubar = .5*(Rlims[1]+Rlims[0]);
	PolarAA_Ractions_struct RA(this,Etot,Lz2,Jz,Delta,taubar,tn);

	anglim = asin((R-taubar)/Delta);
	lsign = SIGN(vR);
	dJdIl[0] = Delta*GaussLegendreQuad(&dJRdEint,-.5*PI,.5*PI,&RA,8)/PI;
	dJdIl[1] = Delta*GaussLegendreQuad(&dJRdLzint,-.5*PI,.5*PI,&RA,8)/PI;
	dJdIl[2] = Delta*GaussLegendreQuad(&dJRdJzint,-.5*PI,.5*PI,&RA,8)/PI;
	dSdI[0] += lsign*Delta*GaussLegendreQuad(&dJRdEint,-.5*PI,anglim,&RA,8);
	dSdI[1] += lsign*Delta*GaussLegendreQuad(&dJRdLzint,-.5*PI,anglim,&RA,8);
	dSdI[2] += lsign*Delta*GaussLegendreQuad(&dJRdJzint,-.5*PI,anglim,&RA,8);

	VecDoub dJdIp = {0.,1.,0.};

	double Determinant = dJdIp[1]*(dJdIl[0]*dJdIn[2]-dJdIl[2]*dJdIn[0]);

	double dIdJ[3][3];
	dIdJ[0][0] = det2(dJdIp[1],dJdIp[2],dJdIn[1],dJdIn[2])/Determinant;
	dIdJ[0][1] = -det2(dJdIl[1],dJdIl[2],dJdIn[1],dJdIn[2])/Determinant;
	dIdJ[0][2] = det2(dJdIl[1],dJdIl[2],dJdIp[1],dJdIp[2])/Determinant;
	dIdJ[1][0] = -det2(dJdIp[0],dJdIp[2],dJdIn[0],dJdIn[2])/Determinant;
	dIdJ[1][1] = det2(dJdIl[0],dJdIl[2],dJdIn[0],dJdIn[2])/Determinant;
	dIdJ[1][2] = -det2(dJdIl[0],dJdIl[2],dJdIp[0],dJdIp[2])/Determinant;
	dIdJ[2][0] = det2(dJdIp[0],dJdIp[1],dJdIn[0],dJdIn[1])/Determinant;
	dIdJ[2][1] = -det2(dJdIl[0],dJdIl[1],dJdIn[0],dJdIn[1])/Determinant;
	dIdJ[2][2] = det2(dJdIl[0],dJdIl[1],dJdIp[0],dJdIp[1])/Determinant;

	angles[0]=dSdI[0]*dIdJ[0][0]+dSdI[1]*dIdJ[1][0]+dSdI[2]*dIdJ[2][0];
	angles[1]=dSdI[0]*dIdJ[0][1]+dSdI[1]*dIdJ[1][1]+dSdI[2]*dIdJ[2][1];
	angles[2]=dSdI[0]*dIdJ[0][2]+dSdI[1]*dIdJ[1][2]+dSdI[2]*dIdJ[2][2];
	angles[3]=det2(dJdIp[1],dJdIp[2],dJdIn[1],dJdIn[2])/Determinant;
	angles[4]=det2(dJdIn[1],dJdIn[2],dJdIl[1],dJdIl[2])/Determinant;
	angles[5]=det2(dJdIl[1],dJdIl[2],dJdIp[1],dJdIp[2])/Determinant;

	angles[4]*=SIGN(Lz);

	if(x[5]<0.0){angles[2]+=PI;}
	if(x[2]<0.0){angles[2]+=2.*PI;}
	if(angles[2]>2.*PI)	angles[2]-=2.0*PI;
	if(angles[0]<0.0)	angles[0]+=2.0*PI;
	if(angles[1]<0.0)	angles[1]+=2.0*PI;
	if(angles[2]<0.0)	angles[2]+=2.0*PI;

	return angles;
}
// ============================================================================
