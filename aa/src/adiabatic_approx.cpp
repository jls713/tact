#include "adiabatic_approx.h"

double J_integrand_AxiSym(double theta, void *params){
	// need to set taubargl and Deltagl
	action_struct_axi * AS = (action_struct_axi *) params;
	double tau=AS->taubargl+AS->Deltagl*sin(theta);
	return sqrt(MAX(0.,ptau2_AxiSym(tau,&(AS->RS))))*cos(theta);
}

double PhiNu(double nu, action_struct_aa *RS){
	double Alpha = RS->AA->CS->alpha(), Gamma = RS->AA->CS->gamma();
	double Lz2 = RS->Lz*RS->Lz;
	double R = sqrt((RS->lamgl+Alpha)*(nu+Alpha)/(Alpha-Gamma));
	double z = sqrt((RS->lamgl+Gamma)*(nu+Gamma)/(Gamma-Alpha));
	return RS->AA->Pot->Phi({R,0.,z})+Lz2*(1./(2.0*R*R)
	        -1./(2.0*(RS->lamgl+RS->AA->CS->alpha())));
}

double nuturnfn(double nu, void *params){
	action_struct_aa * AS = (action_struct_aa *) params;
	return AS->ENugl-PhiNu(nu,AS);
}

double JNuint(double theta, void *params){
	action_struct_aa * AS = (action_struct_aa *) params;
	double Alpha = AS->AA->CS->alpha(), Gamma = AS->AA->CS->gamma();
	double nu=AS->taubargl+AS->Deltagl*sin(theta);
	return sqrt(MAX(0.,0.5*(nu-AS->lamgl)/((nu+Alpha)*(nu+Gamma))
			*(AS->ENugl-PhiNu(nu,AS))))*cos(theta);
}

VecDoub Actions_EllipsoidalAdiabaticApproximation::find_limits(VecDoub tau, action_struct_aa *RS){

	double lambda = tau[0], nu = tau[2];
	VecDoub limits;
	limits.push_back(-CS->gamma()+TINY);
	// find root of p^2(nu)
	double nuouter=nu;
	while(nuturnfn(nuouter, &RS)<0.)	nuouter+=0.1*(-CS->alpha()-nuouter);
	root_find RF(SMALL,100);
	limits.push_back(RF.findroot(&nuturnfn,nu,nuouter,&RS));
	return limits;
}

VecDoub Actions_EllipsoidalAdiabaticApproximation::actions(VecDoub x){
	CS->newalpha(Pot->DeltaGuess(x)+CS->gamma());
	double Lz = fabs(Pot->Lz(x));
	VecDoub tau = CS->xv2tau(x);
	double ENugl = (tau[2]-tau[0])/(4.*(tau[2]+CS->alpha())
	               *(tau[2]+CS->gamma()))*tau[5]*tau[5];
	action_struct_aa RS(this,ENugl,tau[0],Lz,0.,0.);
	ENugl+=PhiNu(tau[2],&RS);
	// Jnu
}
