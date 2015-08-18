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
#include "lmn_orb_ELLz.h"
#include "lmn_orb.h"
#include "debug.h"

// ============================================================================

static double min_distance(double y, void *params){

	root_struct_mindistELLZ *RS = (root_struct_mindistELLZ *) params;
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

std::vector<int> lmn_orb_ELLz::angular_momentum(const VecDoub &x){
	VecDoub xx = {x[0],x[1],x[2]}, vv = {x[3],x[4],x[5]};
	VecDoub ll = cross_product<double>(xx,vv);
	return {sign(ll[0]),sign(ll[1]),sign(ll[2])};
}

int lmn_orb_ELLz::check_ang_mom(double y, double E, int swit){
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

VecDoub lmn_orb_ELLz::check_orbit(double y, double E, int swit, int plot){
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
	root_struct_mindistELLZ *RS = (root_struct_mindistELLZ *) p;
	return RS->E-RS->Pot->Phi({0.,x,0.});
}


double lmn_orb_ELLz::find_closed(double E, int swit){
	// At a fixed energy finds the short (swit=0) or long (swit=1)
	// axis closed loop orbit
	root_struct_mindistELLZ RS(Pot,E,swit);
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
    gsl_set_error_handler_off();
	minimiser1D min(&min_distance,mid,yMin,yMax,1e-4,0.,&status,&RS);
	if(status!=0) return -1.;
	double y = min.minimise(100);
	// if(swit==1) check_orbit(y,E,swit);
	return y;
}

static double min_action(double b, void *params){
	root_struct_actionELLZ *RS = (root_struct_actionELLZ *) params;
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
	root_struct_actionELLZ *RS = (root_struct_actionELLZ *) params;
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


static VecDoub actionSD(Actions_TriaxialStackel_Fudge *ATSF, const std::vector<VecDoub> & results, bool with_freq=false){

	if(with_freq) ATSF->set_freq(true);
	std::vector<VecDoub> Results;
	for(auto Y: results) Results.push_back(ATSF->actions(Y));

	return concatVectors<double>(columnMean<double>(Results),
	                             columnSD<double>(Results),
	                             columnMedian<double>(Results));
}

static double min_actions(const gsl_vector *v, void *params){
	root_struct_actionsELLZ *RS = (root_struct_actionsELLZ *) params;
	double Delta1=gsl_vector_get(v,0);
	double Delta2=gsl_vector_get(v,1);
	// if(Delta2<Delta1) return 1e10;
	if(Delta1>2.*RS->rmax or Delta2>2.*RS->rmax or fabs(Delta2-Delta1)>5.*Delta1 or fabs(Delta2-Delta1)>5.*Delta2) return 1e10;
	double beta = -1.-Delta2*Delta2;
	double alpha = beta-Delta1*Delta1;
	Actions_TriaxialStackel_Fudge ATSF(RS->Pot,alpha,beta);
	VecDoub f = actionSD(&ATSF,RS->orbit->results());
	return sqrt(f[3]*f[3]+f[4]*f[4]+f[5]*f[5]);
}


VecDoub lmn_orb_ELLz::find_best_alphabeta(const VecDoub& X){
	double torb = Pot->torb(X);
	double step = 1e-1*torb;
	Orbit orbit(Pot);
	VecDoub QQ=orbit.integrate(X, 100*step, step);
	//orbit.plot(1,2);
	double norm = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
	root_struct_actionsELLZ RS(Pot,&orbit,ymax);
	VecDoub a1 = {norm/10.,norm/5.};
	VecDoub sizes ={norm/10.,norm/10.};
	minimiser min(&min_actions,a1,sizes,1e-5,&RS);
	VecDoub results;
	results.push_back(min.minimise(&results,100,0));
	return results;
}



void lmn_orb_ELLz::readDeltagrids(const std::string& file){
	std::ifstream infile; infile.open(file);
	if(!infile.is_open())std::cerr<<"Problem: "<<file<<" doesn't exist."<<std::endl;
	double tmp,tmp2,tmp3,tmp4,tmp5,tmp6;
	E=VecDoub(NE,0.);
	Lmax=VecDoub(NE,0.);
	IL=VecDoub(NL,0.);
	ILz=VecDoub(NLz,0.);
	for(int i=0;i<NL;++i)
		IL[i] = (maxIL-minIL)*i/(NL-1)+minIL;
	for(int i=0;i<NLz;++i)
		ILz[i] = (maxIL-minIL)*i/(NLz-1)+minIL;
	Beta=MatMatDoub(NE,MatDoub(NL,VecDoub(NLz,0.)));
	Alpha=MatMatDoub(NE,MatDoub(NL,VecDoub(NLz,0.)));
	int n=0, nE=-1, nL=0, nLz=0;
	while(infile>>tmp>>tmp2>>tmp3>>tmp4>>tmp5>>tmp6){
		if(n%(NL*NLz)==0){nE++;E[nE]=tmp;Lmax[nE]=tmp2;nL=-1;}
		if(n%NL==0){ nL++; nLz=0;}
		Beta[nE][nL][nLz]=tmp6;
		Alpha[nE][nL][nLz]=tmp5;
		++n;++nLz;
	}
	infile.close();
}

void lmn_orb_ELLz::fillDeltagrids(const std::string& file){
	E=VecDoub(NE,0.);
	Lmax=VecDoub(NE,0.);
	IL=VecDoub(NL,0.);
	ILz=VecDoub(NLz,0.);
	for(int i=0;i<NL;++i)
		IL[i] = (maxIL-minIL)*i/(NL-1)+minIL;
	for(int i=0;i<NLz;++i)
		ILz[i] = (maxIL-minIL)*i/(NLz-1)+minIL;
	Beta=MatMatDoub(NE,MatDoub(NL,VecDoub(NLz,0.)));
	Alpha=MatMatDoub(NE,MatDoub(NL,VecDoub(NLz,0.)));

	if(use_log_grid)
		E = create_log_range(E0,Emax,NE);
	else
		E = create_range(E0,Emax,NE);

	#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<NE;i++){
		double y = find_closed(E[i],0);
		double p = sqrt(2.*(E[i]-Pot->Phi({0.,y,0.})));
		Lmax[i]=y*p;
		double Lz,L;
		for(int j=0;j<NL;++j){
			L=Lmax[i]*IL[j];
			double pxpz = L/y, px, pz;
			double py = sqrt(p*p-pxpz*pxpz);
			for(int k=0;k<NLz;++k){
				Lz = L*ILz[k];
				px = Lz/y;
				pz = sqrt(pxpz*pxpz-px*px);
				VecDoub ab = find_best_alphabeta({1e-4,y,1e-4,px,py,pz});
				double beta = -1.-ab[1]*ab[1];
				double alpha = beta-ab[0]*ab[0];
				Alpha[i][j][k]=alpha;
				Beta[i][j][k]=beta;
				std::cout<<E[i]<<" "<<L<<" "<<Lz<<" "<<alpha<<" "<<beta<<" "<<ab[2]<<" "<<y<<" "<<px<<" "<<py<<" "<<pz<<std::endl;
			}
		}
	}
	std::ofstream outfile; outfile.open(file);
	for(unsigned i=0;i<E.size();i++)
		for(unsigned j=0;j<IL.size();j++)
			for(unsigned k=0;k<ILz.size();k++)
				outfile<<E[i]<<" "<<Lmax[i]<<" "<<IL[j]<<" "<<ILz[k]<<" "<<Beta[i][j][k]<<" "<<Alpha[i][j][k]<<std::endl;
	outfile.close();
}

VecDoub lmn_orb_ELLz::findDelta_interp(double En,double L, double Lz){

	int E_bot=0, E_top=0, L_bot1=0, L_bot2=0, L_top1=0, L_top2=0;
	int Lz_bot=0, Lz_top=0;
	double Lrat = Lz/L;

	if(Lrat<minIL)Lrat=minIL;
	if(Lrat>maxIL)Lrat=maxIL;
    topbottom<double>(ILz,Lrat,&Lz_bot,&Lz_top,"Find Delta E grid");

    if(En<E[0]) En=E[0];
    else if(En>E[NE-1]) En = E[NE-1];
    topbottom<double>(E,En,&E_bot,&E_top,"Find Delta E grid");

    double Lintb = Lmax[E_bot]/L;
	if(Lintb<minIL)Lintb=minIL;
	if(Lintb>maxIL)Lintb=maxIL;
    topbottom<double>(IL,Lintb,&L_bot1,&L_top1,"Find Delta E grid");
    double Lintt = Lmax[E_top]/L;
    if(Lintt<minIL)Lintt=minIL;
	if(Lintt>maxIL)Lintt=maxIL;
    topbottom<double>(IL,Lintt,&L_bot2,&L_top2,"Find Delta E grid");

    double fac = (Lrat-ILz[Lz_bot])/(ILz[Lz_top]-ILz[Lz_bot]);

    double alpha_bb = Alpha[E_bot][L_bot1][Lz_bot]+fac*(Alpha[E_bot][L_bot1][Lz_top]-Alpha[E_bot][L_bot1][Lz_bot]);
    double beta_bb = Beta[E_bot][L_bot1][Lz_bot]+fac*(Beta[E_bot][L_bot1][Lz_top]-Beta[E_bot][L_bot1][Lz_bot]);

    double alpha_bt = Alpha[E_bot][L_top1][Lz_bot]+fac*(Alpha[E_bot][L_top1][Lz_top]-Alpha[E_bot][L_top1][Lz_bot])/(ILz[Lz_top]-ILz[Lz_bot]);
    double beta_bt = Beta[E_bot][L_top1][Lz_bot]+fac*(Beta[E_bot][L_top1][Lz_top]-Beta[E_bot][L_top1][Lz_bot]);

    double alpha_tb = Alpha[E_top][L_bot2][Lz_bot]+fac*(Alpha[E_top][L_bot2][Lz_top]-Alpha[E_top][L_bot2][Lz_bot]);
    double beta_tb = Beta[E_top][L_bot2][Lz_bot]+fac*(Beta[E_top][L_bot2][Lz_top]-Beta[E_top][L_bot2][Lz_bot]);

    double alpha_tt = Alpha[E_top][L_top2][Lz_bot]+fac*(Alpha[E_top][L_top2][Lz_top]-Alpha[E_top][L_top2][Lz_bot]);
    double beta_tt = Beta[E_top][L_top2][Lz_bot]+fac*(Beta[E_top][L_top2][Lz_top]-Beta[E_top][L_top2][Lz_bot]);

    fac = (Lintb-IL[L_bot1])/(IL[L_top1]-IL[L_bot1]);
    double alpha_b = alpha_bb+(alpha_bt-alpha_bb)*fac;
    double beta_b = beta_bb+(beta_bt-beta_bb)*fac;

    fac = (Lintt-IL[L_bot2])/(IL[L_top2]-IL[L_bot2]);
    double alpha_t = alpha_tb+(alpha_tt-alpha_tb)*fac;
    double beta_t = beta_tb+(beta_tt-beta_tb)*fac;

    double alpha = alpha_b+(alpha_t-alpha_b)*(En-E[E_bot])/(E[E_top]-E[E_bot]);
    double beta = beta_b+(beta_t-beta_b)*(En-E[E_bot])/(E[E_top]-E[E_bot]);
    return {alpha,beta};
}

VecDoub lmn_orb_ELLz::actions(const VecDoub& x,void*params){
	double En = Pot->H(x);
	double L = Pot->L(x);
	double Lz = fabs(Pot->Lz(x));

	if(params!=nullptr){
		double *alphabeta = (double*)params;
		Actions_TriaxialStackel_Fudge ATSF(Pot,alphabeta[0],alphabeta[1]);
		return ATSF.actions(x);
	}
	else{
		VecDoub ab = findDelta_interp(En,L,Lz);
		Actions_TriaxialStackel_Fudge ATSF(Pot,ab[0],ab[1]);
		return ATSF.actions(x);
	}
}

VecDoub lmn_orb_ELLz::angles(const VecDoub& x,void *params){
	double En = Pot->H(x);
	double L = Pot->L(x);
	double Lz = fabs(Pot->Lz(x));
	if(params!=nullptr){
		double *alphabeta = (double*)params;
		Actions_TriaxialStackel_Fudge ATSF(Pot,alphabeta[0],alphabeta[1]);
		return ATSF.angles(x);
	}
	else{
		VecDoub ab = findDelta_interp(En,L,Lz);
		Actions_TriaxialStackel_Fudge ATSF(Pot,ab[0],ab[1]);
		return ATSF.angles(x);
	}
}

// #include "Multipole.h"

// int main(){
// 	NFW Pot(1.,1.,0.9,0.8);
// 	TestDensity_Hernquist rho(1.,1.,{1.,0.9,0.8});
// 	// MultipoleExpansion ME("/data/jls/self_con/May15/hernquist_we/isotropic_centre_lmax4_NA4_actmin_multiply2_D1_0.05.me");
// 	lmn_orb_ELLz LMN(&Pot,0.05,60.,12,6,6,true,0.02,0.98,"ELLz.tmp");
// 	// lmn_orb LMN2(&Pot,0.05,60,24,true,true);
// 	// printVector(LMN.actions({10.,2.,2.,0.2,0.3,0.2}));
// }

// ============================================================================
