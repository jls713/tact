/*==============================================*/
/* 			Adiabatic Approximations 			*/
/*==============================================*/

// ============================================================================
#include "coordtransforms.h"
#include "adiabatic_aa.h"
// ============================================================================
// Ugh, I don't think the spheroidal Staeckel bits are thread-safe as we are
// altering the coordinate system state every time.

static double pnu_squared(double nu, void *params){
	SpheroidalAA_zactions_struct *RS = (SpheroidalAA_zactions_struct *) params;
	return RS->ENu-RS->AA->Phi_nu({RS->lam,0.,nu},RS->Lz2);
}

static double Jnuint(double theta, void *params){
	SpheroidalAA_zactions_struct *RS = (SpheroidalAA_zactions_struct *) params;
	double nu=RS->taubar+RS->Delta*sin(theta);
	return sqrt(MAX(0.,0.5*(nu-RS->lam)/((nu+RS->AA->CS->alpha())*(nu+RS->AA->CS->gamma()))*pnu_squared(nu,RS)))*cos(theta);
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
			VecDoub xx = {Rgrid[nR],0.,0.1};
			double a = Pot->DeltaGuess(xx)+CS->gamma();
			CS->newalpha(a>CS->gamma()?CS->gamma()-0.1:a);
			double zmax = 1.5+0.25*Rgrid[nR];
			double Rmax = sqrt(MAX((1.0-zmax*zmax/(Rgrid[nR]*Rgrid[nR]-CS->alpha()+CS->gamma()))*Rgrid[nR]*Rgrid[nR],0.0));
			if(Rmax==0.0){Rmax = 0.5; zmax = sqrt((1.0-Rmax*Rmax/(Rgrid[nR]*Rgrid[nR]))*(Rgrid[nR]*Rgrid[nR]-CS->alpha()+CS->gamma()));}
			Ezmaxgrid[nL][nR] = Pot->Phi({Rmax,0.,zmax})-Pot->Phi({Rgrid[nR],0.,0.})+.5*Lz2*(1./(Rmax*Rmax)-1./(Rgrid[nR]*Rgrid[nR]));
			Ezgrid[nL][nR][0]=Jzgrid[nL][nR][0]=0;
			double nulast=-CS->gamma()+SMALL, nulim = 0.;

			for(int i=1; i<NGRID; i++){
				Ezgrid[nL][nR][i]=((double)i)*Ezmaxgrid[nL][nR]/(double)(NGRID-1);
				Jzgrid[nL][nR][i]=actions_Jz(Ezgrid[nL][nR][i],Lz2,{Rgrid[nR]*Rgrid[nR]-CS->alpha(),0.,nulast},&nulim);
				nulast = nulim;
				std::cout<<Rgrid[nR]<<" "<<nulim<<" "<<Lz2<<" "<<Ezgrid[nL][nR][i]<<" "<<Jzgrid[nL][nR][i]<<std::endl;
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
		topbottom<double>(Rgrid,R,&botR,&topR);

		double Lz = sqrt(Lz2);
		int botL,topL;
		if(Lz<Lgrid[0]) Lz = Lgrid[0];
		if(Lz>Lgrid[NL-1]) Lz = Lgrid[NL-1];
		topbottom<double>(Lgrid,Lz,&botL,&topL);

		double fac=(R-Rgrid[botR])/(Rgrid[topR]-Rgrid[botR]);
		double fac2=(Lz-Lgrid[botL])/(Lgrid[topL]-Lgrid[botL]);

		double Ezb=linterp<double>(Jzgrid[botL][botR],Ezgrid[botL][botR],MIN(Jz,Jzgrid[botL][botR][NGRID-1]));
		double Ezt=linterp(Jzgrid[botL][topR],Ezgrid[botL][topR],MIN(Jz,Jzgrid[botL][topR][NGRID-1]));

		double Ezb2=linterp<double>(Jzgrid[topL][botR],Ezgrid[topL][botR],MIN(Jz,Jzgrid[topL][botR][NGRID-1]));
		double Ezt2=linterp(Jzgrid[topL][topR],Ezgrid[topL][topR],MIN(Jz,Jzgrid[topL][topR][NGRID-1]));

		return (1-fac2)*((1-fac)*Ezb+fac*Ezt)+fac2*((1-fac)*Ezb2+fac*Ezt2);
	}
}

double Actions_SpheroidalAdiabaticApproximation::find_nulimit(double ENu, double Lz2, VecDoub tau){
	SpheroidalAA_zactions_struct RS(this,ENu,Lz2,tau[0],0.,0.);
	double nutry = tau[2];
	while(pnu_squared(nutry, &RS)>0.) nutry+=0.1*(-CS->alpha()-nutry);
	root_find RF(TINY,20);
	return RF.findroot(&pnu_squared,tau[2],nutry,&RS);
}

VecDoub Actions_SpheroidalAdiabaticApproximation::find_lamlimits(double Elam, double Lz2,double Jz,VecDoub tau){
	SpheroidalAA_Ractions_struct RS(this,Elam,Lz2,Jz,tau[2],0.,0.);
	double lamin = tau[0], lamout = tau[0];
	while(plam_squared(lamout, &RS)>0.) lamout*=1.1;;
	while(plam_squared(lamin, &RS)>0.)  lamin-=.1*(lamin+CS->alpha());
	root_find RF(TINY,20);
	return {RF.findroot(&plam_squared,lamin,tau[0],&RS),
			RF.findroot(&plam_squared,tau[0],lamout,&RS)};
}

double Actions_SpheroidalAdiabaticApproximation::actions_Jz(double ENu, double Lz2, VecDoub tau, double *nulim){
	// Calculates vertical action -- only require R and Ez but for speed reasons it is efficient to pass the input z value and output the found z limit
	*nulim = find_nulimit(ENu,Lz2,tau);
	double Delta = .5*(*nulim+CS->gamma()), taubar = .5*(*nulim-CS->gamma());
	SpheroidalAA_zactions_struct PA(this,ENu,Lz2,tau[0],Delta,taubar);
	return 2.*Delta*GaussLegendreQuad(&Jnuint,-.5*PI,.5*PI,&PA)/PI;

}

VecDoub Actions_SpheroidalAdiabaticApproximation::actions(const VecDoub& x, void*params){

	if(params){
		double *deltaguess = (double *) params;
		if(*deltaguess>0.) CS->newalpha(-Pot->DeltaGuess(x)-1.);
		else CS->newalpha(*deltaguess);
	}

	if(CS->alpha()>CS->gamma())CS->newalpha(CS->gamma()-0.1);

	VecDoub tau = CS->xv2tau(x);
	double Lz = Pot->Lz(x), Lz2 = Lz*Lz;
	double ENu = (tau[2]-tau[0])/(8.*(tau[2]+CS->alpha())
	               *(tau[2]+CS->gamma()))*tau[5]*tau[5]
				 +Phi_nu(tau,Lz2);

	double nulim = tau[2];
	double Jz = actions_Jz(ENu, Lz2, tau, &nulim);

	double Elam = (tau[0]-tau[2])*tau[3]*tau[3]/(8.0*(tau[0]+CS->alpha())*(tau[0]+CS->gamma()))+Philam_eff(tau,Lz2,Jz);
	VecDoub lamlims = find_lamlimits(Elam,Lz2,Jz,tau);
	double Delta = .5*(lamlims[1]-lamlims[0]), taubar = .5*(lamlims[1]+lamlims[0]);
	SpheroidalAA_Ractions_struct RA(this,Elam,Lz2,Jz,tau[2],Delta,taubar);
	double JR = Delta*GaussLegendreQuad(&Jlamint,-.5*PI,.5*PI,&RA)/PI;

	return {JR,Lz,Jz};
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

void Actions_PolarAdiabaticApproximation::load_grids(std::string filename){

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

void Actions_PolarAdiabaticApproximation::make_grids(std::string filename){

	for(int n=0; n<NGRID; n++){
		Rgrid.push_back(.75*Rmin+n*(1.5*Rmax-.75*Rmin)/((double)(NGRID-1)));
		Ezmaxgrid.push_back(Pot->Phi({Rgrid[n],0.,30.})-Pot->Phi({Rgrid[n],0.,0.}));
		Ezgrid.push_back(VecDoub(NGRID,0.));
		Jzgrid.push_back(VecDoub(NGRID,0.));
		double zlast=.005, zlim = 0.;

		for(int i=1; i<NGRID; i++){
			Ezgrid[n][i]=((double)i)*Ezmaxgrid[n]/(double)(NGRID-1);
			Jzgrid[n][i]=actions_Jz(Rgrid[n],Ezgrid[n][i],zlast,&zlim);
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

double Actions_PolarAdiabaticApproximation::Ez_from_grid(double R,double Jz){
	//Ez(R,Jz) by interpolation
	if(no_energy_correction)
		return 0;
	else{
		int botR,topR;
		if(R<Rgrid[0]){ R = Rgrid[0]; std::cerr<<"Outside R grid in PAA:"<<R<<"\n";}
		if(R>Rgrid[NGRID-1]){ R = Rgrid[NGRID-1]; std::cerr<<"Outside R grid in PAA:"<<R<<"\n";}
		topbottom<double>(Rgrid,R,&botR,&topR);
		double fac=(R-Rgrid[botR])/(Rgrid[topR]-Rgrid[botR]);
		double Ezb=linterp<double>(Jzgrid[botR],Ezgrid[botR],MIN(Jz,Jzgrid[botR][NGRID-1]));
		double Ezt=linterp(Jzgrid[topR],Ezgrid[topR],MIN(Jz,Jzgrid[topR][NGRID-1]));
		return (1-fac)*Ezb+fac*Ezt;
	}
}

double Actions_PolarAdiabaticApproximation::find_zlimit(double Ez, VecDoub Polar){
	PolarAA_zactions_struct RS(this,Ez,Polar[0],0.);
	double ztry = fabs(Polar[2]);
	while(pz_squared(ztry, &RS)>0.) ztry*=1.1;
	root_find RF(TINY,20);
	return RF.findroot(&pz_squared,fabs(Polar[2]),ztry,&RS);
}

VecDoub Actions_PolarAdiabaticApproximation::find_Rlimits(double R, double Etot, double Lz2,double Jz){
	PolarAA_Ractions_struct RS(this,Etot,Lz2,Jz,0.,0.);
	double Rin = R, Rout = R;
	while(pR_squared(Rout, &RS)>0.) Rout*=1.01;
	while(pR_squared(Rin, &RS)>0.)  Rin*=.99;
	root_find RF(TINY,20);
	return {RF.findroot(&pR_squared,Rin,R,&RS),
			RF.findroot(&pR_squared,R,Rout,&RS)};
}

double Actions_PolarAdiabaticApproximation::actions_Jz(double R, double Ez,double z, double *zlim){
	// Calculates vertical action -- only require R and Ez but for speed reasons it is efficient to pass the input z value and output the found z limit
	*zlim = find_zlimit(Ez, {R,0.,z});
	PolarAA_zactions_struct PA(this,Ez,R,*zlim);
	return 2.*(*zlim)*GaussLegendreQuad(&Jzint,0.,.5*PI,&PA)/PI;

}

VecDoub Actions_PolarAdiabaticApproximation::actions(const VecDoub& x, void*params){

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
	PolarAA_Ractions_struct RA(this,Etot,Lz2,Jz,Delta,taubar);
	double JR = Delta*GaussLegendreQuad(&JRint,-.5*PI,.5*PI,&RA)/PI;

	return {JR,Lz,Jz};
}
