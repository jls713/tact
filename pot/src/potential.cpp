/*======================================*/
/* 			    Potential_JS 			*/
/*======================================*/
#include <Python.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "coordsys.h"
#include "coordtransforms.h"
#include "potential.h"

// ============================================================================
struct RE_st{
	double E; Potential_JS *Pot;
	RE_st(double E, Potential_JS *Pot):E(E), Pot(Pot){}
};

double RE_fn(double r, void*params){
	RE_st *RE = (RE_st*)params;
	return -0.5*r*RE->Pot->Forces({r,0.,0.})[0]+RE->Pot->Phi({r,0.,0.})-RE->E;
}
struct RL_st{
	double Lz2; Potential_JS *Pot;
	RL_st(double Lz, Potential_JS *Pot):Lz2(Lz*Lz), Pot(Pot){}
};

double RL_fn(double R, void*params){
	RL_st *RL = (RL_st*)params;
	return RL->Pot->Forces({R,0.,0.})[0]+RL->Lz2/pow(R,3);
}

double Potential_JS::R_E(double E,double r){
	// Can pass initial guess for speed
	RE_st RES(E,this);
	root_find RF(1e-6,200);
	double up = 1e5,down=1e-5;
	if(r>0.){up=10.*r;down=0.1*r;}
    return RF.findroot(&RE_fn,down,up,&RES);
}

double Potential_JS::R_L(double Lz,double r){
	// Can pass initial guess for speed
	RL_st RLS(Lz,this);
	root_find RF(1e-6,200);
	double up = 1e5,down=1e-5;
	if(r>0.){up=10.*r;down=0.1*r;}
    return RF.findroot(&RL_fn,down,up,&RLS);
}

double Lzmax_fn(double r, void*params){
	RE_st *RE = (RE_st*)params;
	return -2.*r*r*(RE->E-RE->Pot->Phi({r,0.,0.}));
}

double Potential_JS::Lzmax(double E,double RC){
	RE_st RES(E,this);int status;
	minimiser1D min(&Lzmax_fn,RC/2.,SMALL,RC,1e-4,0.,&status,&RES);
    double r = min.minimise(100);
    return -Lzmax_fn(r,&RES);
}

double Potential_JS::R_L(const VecDoub &x){
	assert(x.size()==6);
	return R_L(Lz(x),norm<double>({x[0],x[1],x[2]}));
}

double Potential_JS::R_E(const VecDoub &x){
	assert(x.size()==6);
	return R_E(H(x),norm<double>({x[0],x[1],x[2]}));
}

double Potential_JS::L_E(double E){
	double R_e = R_E(E);
	return sqrt(pow(R_e,3.)*-Forces({R_e,0.,0.})[0]);
}

double Potential_JS::L_E(const VecDoub &x){
	double R_e = R_E(x);
	return sqrt(pow(R_e,3.)*-Forces({R_e,0.,0.})[0]);
}

double Potential_JS::torb(const VecDoub &x){
	double R = R_E(x);
	return 2.*PI*sqrt(R/-Forces({R,0.,0.})[0]);
}

VecDoub Potential_JS::dPhidRdz(const VecDoub& Rz){
	double Delta = 0.001;

	VecDoub Rplus =  Forces({Rz[0]+Delta,0.,Rz[1]});
	VecDoub Rminus = Forces({Rz[0]-Delta,0.,Rz[1]});
	VecDoub zplus =  Forces({Rz[0],0.,Rz[1]+Delta});
	VecDoub zminus = Forces({Rz[0],0.,Rz[1]-Delta});
	return {0.5*(Rminus[0]-Rplus[0])/Delta,
			0.5*(Rminus[2]-Rplus[2])/Delta,
			0.5*(zminus[2]-zplus[2])/Delta};
}

double Potential_JS::DeltaGuess(const VecDoub& x){
	// Returns a guess of Gamma-Alpha=Delta^2 assuming dldv((l-v)V)=0
	double R = norm<double>({x[0],x[1]}), z = x[2];
	VecDoub F = Forces({R,0.,z});
	VecDoub d2P = dPhidRdz({R,z});
	return z*z-R*R+(-3.0*z*F[0]+3.0*R*F[2]+R*z*(d2P[0]-d2P[2]))/d2P[1];
}

struct pot_int_struc{
	Potential_JS *Pot;
	double intercept;
	int coeff;
	pot_int_struc(Potential_JS *Pot,double intercept,int coeff):Pot(Pot), intercept(intercept),coeff(coeff){}
};

static double find_potential_int(double x, void *params){
	pot_int_struc *RS = (pot_int_struc *) params;
	VecDoub Phi(3,1e-5);
	Phi[RS->coeff]=x;
	return RS->Pot->Phi(Phi)-RS->intercept;
}

double Potential_JS::find_potential_intercept(double Phi0, int direction,double xmin,double xmax){
	root_find RF(1e-4,100);
	pot_int_struc PP(this,Phi0,direction);
	return RF.findroot(&find_potential_int,xmin,xmax,&PP);
}

// ============================================================================
// Prolate Stackel Perfect Ellipsoid Potential
// ============================================================================

double StackelProlate_PerfectEllipsoid::G(double tau){
	/* de Zeeuw's G(tau) function */
	double Gamma = CS->gamma(), sqG =sqrt(-Gamma/(tau+Gamma));
	return Const*sqG*atan(1./sqG);
}
double StackelProlate_PerfectEllipsoid::GPrime(double tau){
	/* derivative of G wrt tau */
	double Gamma = CS->gamma(), sqG =sqrt(-Gamma/(tau+Gamma));
	return 0.5*Const*sqG*sqG*(sqG*atan(1./sqG)/Gamma+1./tau);
}

VecDoub StackelProlate_PerfectEllipsoid::Vderivs(const VecDoub& tau){
	/* Calculates the derivatives of the Potential_JS wrt tau given tau=(l,v) */
	double Gl = G(tau[0]), Gv = G(tau[1]), Gamma = CS->gamma();
	double dVdl = 	(-GPrime(tau[0])*(tau[0]+Gamma)-Gl
			+(Gl*(tau[0]+Gamma)-Gv*(tau[1]+Gamma))/(tau[0]-tau[1]))
			/(tau[0]-tau[1]);
	double dVdv = 	(GPrime(tau[1])*(tau[1]+Gamma)+Gv
			-(Gl*(tau[0]+Gamma)-Gv*(tau[1]+Gamma))/(tau[0]-tau[1]))
			/(tau[0]-tau[1]);
	VecDoub derivs = {dVdl, dVdv};
	return derivs;
}


VecDoub StackelProlate_PerfectEllipsoid::Forces(const VecDoub& x){
	/* forces at Cartesian x */
	VecDoub derivs = CS->derivs(x);
	VecDoub tau = {derivs[0],derivs[1]};
	VecDoub Vderiv = Vderivs(tau);

	double dvdR = -Vderiv[0]*derivs[2]-Vderiv[1]*derivs[4];
	double R = norm<double>({x[0],x[1]});
	VecDoub result ={ 	x[0]*dvdR/R, x[1]*dvdR/R,
			  			-Vderiv[0]*derivs[3]-Vderiv[1]*derivs[5]};
	return result;
}

double StackelProlate_PerfectEllipsoid::Phi_tau(const VecDoub& tau){
	/* Potential at tau */
	double Gamma = CS->gamma();
	return -((tau[0]+Gamma)*G(tau[0])-(tau[2]+Gamma)*G(tau[2]))/(tau[0]-tau[2]);
}

double StackelProlate_PerfectEllipsoid::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	VecDoub tau = CS->x2tau(x);
	return Phi_tau(tau);
}

VecDoub StackelProlate_PerfectEllipsoid::x2ints(const VecDoub& x, VecDoub *tau){
	VecDoub Ints = {H(x), 0.5*pow(Lz(x),2.)};
	if(!tau) (*tau) = CS->xv2tau(x);
	Ints.push_back(
	 ((*tau)[0]+CS->gamma())*
	 	(Ints[0]-(Ints[1]/((*tau)[0]+CS->alpha()))+G((*tau)[0]))
	 -(pow(((*tau)[3]*((*tau)[0]-(*tau)[2])),2.0))
	 	/(8.0*((*tau)[0]+CS->alpha())*((*tau)[0]+CS->gamma())));
	Ints.push_back(
	 ((*tau)[2]+CS->gamma())*
	 	(Ints[0]-(Ints[1]/((*tau)[2]+CS->alpha()))+G((*tau)[2]))
	 -(pow(((*tau)[5]*((*tau)[0]-(*tau)[2])),2.0))
	 	/(8.0*((*tau)[2]+CS->alpha())*((*tau)[2]+CS->gamma())));
	// Ints[3]=Ints[2];
	return Ints;
}

// ============================================================================
// Triaxial Stackel Perfect Ellipsoid Potential
// ============================================================================

struct TriaxialStackel_GIntegrand_struct{
	double taugl, acgl, bcgl, c2gl;
	TriaxialStackel_GIntegrand_struct(double taugl, double acgl, double bcgl, double c2gl):taugl(taugl),acgl(acgl),bcgl(bcgl),c2gl(c2gl){};
};

static double G_integrand(double s,void* params){
	TriaxialStackel_GIntegrand_struct *TG =
		(TriaxialStackel_GIntegrand_struct* )params;
	return sqrt(TG->c2gl+TG->bcgl*s*s)/sqrt(TG->c2gl+TG->acgl*s*s)/(TG->c2gl+(TG->taugl-TG->c2gl)*s*s);
}

double StackelTriaxial::G(double tau){
	/* de Zeeuw's perfect ellipsoid G function 			*/
	/* calculated using GL integration of eq B9 of dZ85 */
	TriaxialStackel_GIntegrand_struct P(tau,a*a-c*c,b*b-c*c,c*c);
	return Const*GaussLegendreQuad(&G_integrand,0.,1.,&P);
}

static double GP_integrand(double s,void* params){
	TriaxialStackel_GIntegrand_struct *TG =
		(TriaxialStackel_GIntegrand_struct* )params;
	double p = TG->c2gl+(TG->taugl-TG->c2gl)*s*s;
	return -sqrt(TG->c2gl+TG->bcgl*s*s)/sqrt(TG->c2gl+TG->acgl*s*s)*s*s/p/p;
}
double StackelTriaxial::GPrime(double tau){
	/* de Zeeuw's perfect ellipsoid GPrime function 	*/
	/* calculated using GL integration of eq B9 of dZ85 */
	TriaxialStackel_GIntegrand_struct P(tau,a*a-c*c,b*b-c*c,c*c);
	return Const*GaussLegendreQuad(&GP_integrand,0.,1.,&P);
}

/*
// Attempts at analytic expressions for above -- doesn't work

double StackelTriaxial::G(double tau){
	// std::cout<<l<<" "<<sinm<<" "<<-(tau+CS->alpha())/(CS->alpha()-CS->gamma())<<std::endl;
	// std::cout<<ellint_third(l,sinm,-(tau+CS->alpha())/(CS->alpha()-CS->gamma()))<<std::endl;
	std::cout<<ellint_third(l,sinm,0.)<<" "<<Flm<<" "<<Elm<<std::endl;
	return Const/(tau+CS->alpha())*
		   ((tau+CS->beta())*ellint_third(l,sinm,-(tau+CS->alpha())/(CS->alpha()-CS->gamma()))
		   +(CS->alpha()-CS->beta())*Flm);
}

double StackelTriaxial::GPrime(double tau){
	return 0.5*Const/(tau+CS->alpha())/(tau+CS->gamma())*
		   (((tau+CS->alpha())*(tau+CS->gamma())-(tau+CS->alpha())*(tau+CS->beta())-(tau+CS->gamma())*(tau+CS->beta()))
		   	*ellint_third(l,sinm,-(tau+CS->alpha())/(CS->alpha()-CS->gamma()))/(tau+CS->alpha())
		   +((tau+CS->alpha())*(CS->beta()-CS->gamma())+(tau+CS->gamma())*(CS->beta()-CS->alpha()))
		    *Flm/(tau+CS->alpha())
		   -(tau+CS->alpha())*b*c*sin(l)/CS->alpha()+(CS->gamma()-CS->alpha())*Elm);
}
*/

VecDoub StackelTriaxial::Vderivs(const VecDoub& tau){
	/* Calculates the derivatives of the Potential_JS wrt tau given tau=(l,m,v) */
	double Gl = G(tau[0]), Gm = G(tau[1]), Gn = G(tau[2]), Alpha = CS->alpha(), Gamma = CS->gamma();
	double  F_l = (tau[0]+Alpha)*(tau[0]+Gamma)*Gl,
			F_m = (tau[1]+Alpha)*(tau[1]+Gamma)*Gm,
			F_n = (tau[2]+Alpha)*(tau[2]+Gamma)*Gn;
	double dlm = tau[0]-tau[1], dln = tau[0]-tau[2];
	double dml = tau[1]-tau[0], dmn = tau[1]-tau[2];
	double dnl = tau[2]-tau[0], dnm = tau[2]-tau[1];

	double dVdl = -F_l*GPrime(tau[0])/Gl/dlm/dln-(2.*tau[0]+Alpha+Gamma)*Gl/dlm/dln+F_l/dlm/dln*(1./dln+1./dlm)
					-F_m/dmn/dml/dml-F_n/dnl/dnl/dnm;
	double dVdm = -F_m*GPrime(tau[1])/Gm/dml/dmn-(2.*tau[1]+Alpha+Gamma)*Gm/dmn/dml+F_m/dml/dmn*(1./dmn+1./dml)
					-F_l/dlm/dlm/dln-F_n/dnl/dnm/dnm;
	double dVdn = -F_n*GPrime(tau[2])/Gn/dnl/dnm-(2.*tau[2]+Alpha+Gamma)*Gn/dnm/dnl+F_n/dnl/dnm*(1./dnm+1./dnl)
					-F_l/dlm/dln/dln-F_m/dmn/dmn/dml;
	VecDoub derivs = {dVdl, dVdm, dVdn};
	return derivs;
}

VecDoub StackelTriaxial::Forces(const VecDoub& x){
	/* forces at Cartesian x */
	VecDoub derivs = CS->derivs(x);
	VecDoub tau = {derivs[0],derivs[1],derivs[2]};
	VecDoub Vderiv = Vderivs(tau);
	VecDoub result ={ 	-Vderiv[0]*derivs[3]-Vderiv[1]*derivs[4]-Vderiv[2]*derivs[5],
						-Vderiv[0]*derivs[6]-Vderiv[1]*derivs[7]-Vderiv[2]*derivs[8],
			  			-Vderiv[0]*derivs[9]-Vderiv[1]*derivs[10]-Vderiv[2]*derivs[11]};
	return result;
}

double StackelTriaxial::Phi_tau(const VecDoub& tau){
	/* Potential at tau */
	double  F_l = (tau[0]+CS->alpha())*(tau[0]+CS->gamma())*G(tau[0]),
			F_m = (tau[1]+CS->alpha())*(tau[1]+CS->gamma())*G(tau[1]),
			F_n = (tau[2]+CS->alpha())*(tau[2]+CS->gamma())*G(tau[2]);
	return -F_l/(tau[0]-tau[1])/(tau[0]-tau[2])-F_m/(tau[1]-tau[2])/(tau[1]-tau[0])-F_n/(tau[2]-tau[0])/(tau[2]-tau[1]);
}

double StackelTriaxial::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	VecDoub tau = CS->x2tau(x);
	return Phi_tau(tau);
}

VecDoub StackelTriaxial::tau2ints(const VecDoub& tau){
	VecDoub pp = CS->tau2p(tau);
	double X = 0.5*pp[0]-(tau[0]+CS->alpha())*(tau[0]+CS->gamma())*G(tau[0])/(tau[0]-tau[1])/(tau[0]-tau[2]);
	double Y = 0.5*pp[1]-(tau[1]+CS->alpha())*(tau[1]+CS->gamma())*G(tau[1])/(tau[1]-tau[0])/(tau[1]-tau[2]);
	double Z = 0.5*pp[2]-(tau[2]+CS->alpha())*(tau[2]+CS->gamma())*G(tau[2])/(tau[2]-tau[1])/(tau[2]-tau[0]);
	VecDoub Ints = {X+Y+Z};
	double J =(tau[1]+tau[2])*X+(tau[2]+tau[0])*Y+(tau[0]+tau[1])*Z;
	double K = tau[1]*tau[2]*X+tau[2]*tau[0]*Y+tau[0]*tau[1]*Z;
	Ints.push_back((CS->alpha()*(CS->alpha()*Ints[0]+J)+K)/(CS->alpha()-CS->gamma()));
	Ints.push_back((CS->gamma()*(CS->gamma()*Ints[0]+J)+K)/(CS->gamma()-CS->alpha()));
	return Ints;
}
// ============================================================================
// PowerLaw Potential
// ============================================================================
double PowerLaw::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	double r = x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2;
	return -GM*pow(r,-.5*k);
}
VecDoub PowerLaw::Forces(const VecDoub& x){
	/* Forces at Cartesian x */
	double r = x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2;
	r = -k*GM*pow(r,-.5*k-1);
	return {x[0]*r,x[1]*r/q1,x[2]*r/q2};
}
double PowerLaw::density(double r){
	return GM/conv::FPG*k*(1-k)*pow(r,-k-2.);
}
// ============================================================================
// Isochrone Potential
// ============================================================================
double Isochrone::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	double r = x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2;
	return -GM/(b+sqrt(r+b*b));
}
VecDoub Isochrone::Forces(const VecDoub& x){
	/* Forces at Cartesian x */
	double r = x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2;
	r = GM/(b+sqrt(r+b*b))/(b+sqrt(r+b*b))/sqrt(r+b*b);
	return {-x[0]*r,-x[1]*r/q1,-x[2]*r/q2};
}
double Isochrone::density(double r){
	double a = sqrt(b*b+r*r);
	return GM/conv::FPG*(3.*(b+a)*a*a-r*r*(b+3.*a))/pow(a*(b+a),3.);
}

// ============================================================================
// Logarithmic Potential
// ============================================================================
double Logarithmic::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	// Phi0 is to make potential negative so E<0 for orbits
	return Vc2/2.*log(x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2)-Phi0;
}
VecDoub Logarithmic::Forces(const VecDoub& x){
	/* Forces at Cartesian x */
	double r = Vc2/(x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2);
	VecDoub Forces = {x[0],x[1]/q1,x[2]/q2};
	Forces = Forces*-r;
	return Forces;
}


// ============================================================================
// HarmonicOscillator Potential
// ============================================================================
double HarmonicOscillator::Phi(const VecDoub &x){
	/* potential at Cartesian x */
	double P = 0.;
	for(int i=0;i<3;i++) P+=Om[i]*Om[i]*x[i]*x[i];
	return P*.5;
}
VecDoub HarmonicOscillator::Forces(const VecDoub &x){
	/* Forces at Cartesian x */
	VecDoub F(3,0);
	for(int i=0;i<3;i++) F[i]=-Om[i]*Om[i]*x[i];
	return F;
}

// ============================================================================
// Dehnen Potential
// ============================================================================
double Dehnen::Phi(const VecDoub &x){
	/* potential at Cartesian x */
	assert(x.size()==3);
	double r = norm<double>(x);
	double chi = pow(r/rs,1./alpha);
	chi=chi/(1+chi);
	return -conv::FPG*rhoS*rs*rs*alpha*
	(rs/r*incomplete_beta(alpha*(3-gamma),alpha*(beta-3),chi)
	+incomplete_beta(alpha*(beta-2),alpha*(2-gamma),1-chi));
}
VecDoub Dehnen::Forces(const VecDoub &x){
 	/* Forces at Cartesian x */
 	assert(x.size()==3);
 	double r = norm<double>(x);
	double chi = pow(r/rs,1./alpha);
	double dchi = chi/r/alpha/(1+chi)*(1.-chi/(1+chi));
	chi = chi/(1+chi);
	r = -conv::FPG*rhoS*rs*rs*alpha*
	(-rs/r/r*incomplete_beta(alpha*(3-gamma),alpha*(beta-3),chi)
	+rs/r*pow(chi,alpha*(3-gamma)-1)*pow(1-chi,alpha*(beta-3)-1)*dchi
	-pow(1-chi,alpha*(beta-2)-1)*pow(chi,alpha*(2-gamma)-1)*dchi);
 	VecDoub F = x*-r;
 	return F;
}
double Dehnen::Density(const VecDoub& x){
	assert(x.size()==3);
	double r = norm<double>(x)/rs;
	return rhoS*pow(r,-gamma)*pow(1.+pow(r,1./alpha),(gamma-beta)*alpha);

}
double Dehnen::TotalMass(){
	return 4.*PI*rhoS*rs*rs*rs*alpha
			*complete_beta(alpha*(3-gamma),alpha*(beta-3));
}

// ============================================================================
// Miyamoto-Nagai Potential
// ============================================================================
double MiyamotoNagai_JS::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	double AZB=A+sqrt(x[2]*x[2]+Bq);
	return -GM/sqrt(x[0]*x[0]+x[1]*x[1]+AZB*AZB);
}

double MiyamotoNagai_JS::Vc(double R){
	double t = A+sqrt(Bq);
	double F = 1./(R*R+t*t);
	return sqrt(GM*R*R*F*sqrt(F));
}

VecDoub MiyamotoNagai_JS::Forces(const VecDoub& x){
	/* Forces at Cartesian x */
	double  ZB=sqrt(x[2]*x[2]+Bq), AZB = A + ZB;
	double f = 1./(x[0]*x[0]+x[1]*x[1]+AZB*AZB),rtF=sqrt(f);
	VecDoub Force = {-GM*f*rtF*x[0], -GM*x[1]*f*rtF, ZB? -GM*x[2]*f*rtF*AZB/ZB : 0.};
	return Force;
}
// ============================================================================
// Jaffe Bulge Potential
// ============================================================================

double Bulge::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	double r = sqrt(x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2);
	return GM/b_bulge*log(r/(r+b_bulge));
}
VecDoub Bulge::Forces(const VecDoub& x){
	// Forces at Cartesian x
	double r = sqrt(x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2);
	double dpdr = GM/b_bulge*(1./r-1./(r+b_bulge))/r;
	VecDoub Force = {x[0], x[1]/q1, x[2]/q2};
	Force = Force*-dpdr;
	return Force;
}

// ============================================================================
// NFW Potential
// ============================================================================

double NFW::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	double r = sqrt(x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2);
	return -GM*log(1.+r/rs)/r;
}
VecDoub NFW::Forces(const VecDoub& x){
	/* Forces at Cartesian x */
	double r = sqrt(x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2);
	double dpdr = GM*(log(1.+r/rs)/r-1./rs/(1.+r/rs))/r;
	VecDoub Force = {x[0], x[1]/q1, x[2]/q2};
	Force = Force*(-dpdr/r);
	return Force;
}
double NFW::density(const VecDoub& x){
	double Delta = 0.005;
	VecDoub xtmp=x;
	xtmp[0]+=Delta;		VecDoub plus = Forces(xtmp);
	xtmp[0]-=2.*Delta; 	VecDoub minus = Forces(xtmp);
	double d2p = (minus[0]-plus[0])/2./Delta;
	xtmp[0]=x[0];
	xtmp[1]+=Delta;		plus = Forces(xtmp);
	xtmp[1]-=2.*Delta; 	minus = Forces(xtmp);
	d2p += (minus[1]-plus[1])/2./Delta;
	xtmp[1]=x[1];
	xtmp[2]+=Delta;		plus = Forces(xtmp);
	xtmp[2]-=2.*Delta; 	minus = Forces(xtmp);
	d2p += (minus[2]-plus[2])/2./Delta;
	xtmp[2]=x[2];
	return d2p/4./PI;
}
// ============================================================================
// Hernquist Potential
// ============================================================================

double Hernquist::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	double r = sqrt(x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2);
	return -GM/(rs+r);
}
VecDoub Hernquist::Forces(const VecDoub& x){
	/* Forces at Cartesian x */
	double r = sqrt(x[0]*x[0]+x[1]*x[1]/q1+x[2]*x[2]/q2);
	double dpdr = GM/(rs+r)/(rs+r);
	VecDoub Force = {x[0], x[1]/q1, x[2]/q2};
	Force = Force*(-dpdr/r);
	return Force;
}
// ============================================================================
// GalPot Potential for interface with Walter Dehnen's code
// ============================================================================
GalPot::GalPot(std::string TpotFile){
	if(TpotFile.substr(TpotFile.find_last_of(".") + 1)!="Tpot")
		std::cerr<<"No functionality implemented for passing anything but Tpot file to Galpot."<<std::endl;
	std::ifstream file;
	file.open(TpotFile);
	if(!file.is_open())
		std::cerr<<TpotFile<<" cannot be opened."<<std::endl;
	PhiWD=new GalaxyPotential(file);
	file.close();
}

double GalPot::Phi(const VecDoub& x){
	/* potential at Cartesian x */
	double R = norm<double>({x[0],x[1]});
	return conv::kpcMyr2kmsSq*(*PhiWD)(R,x[2]);
}
VecDoub GalPot::Forces(const VecDoub& x){
	/* Forces at Cartesian x */
	double R = norm<double>({x[0],x[1]});
	double dR, dz;
	(*PhiWD)(R,x[2],dR,dz);
	VecDoub f = {x[0]*dR/R,x[1]*dR/R, dz};
	f = f*-conv::kpcMyr2kmsSq;
	return f;
}

VecDoub GalPot::freqs(double R){
	/* frequencies: kappa, omega_c and nu at polar R */
	Frequencies freqs=(*PhiWD).KapNuOm(R);
	VecDoub freqs_v = {freqs(0),freqs(1),freqs(2)};
	freqs_v = freqs_v*conv::kpcMyr2kms;
	return freqs_v;
}

double GalPot::Vc(double R){
	return sqrt(R*-Forces({R,0.,0.})[0]);
}

// ============================================================================
// Bowden NFW Potential
// ============================================================================

double BowdenNFW::Phi(const VecDoub &x){
	double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	double ct = x[2]/r;ct*=ct;
	double pot = -rho0*log(1+r/rs)/r;
	pot+=rho1*r/pow(r+r1,2.)*(3*ct-1.);
	double R2 = (x[0]*x[0]+x[1]*x[1]);
	double iR2 = 1./R2;
	double c2p = (R2>0.?(x[0]*x[0]-x[1]*x[1])*iR2:0.); // cos(2\phi)
	ct=1-ct;// sin^2(theta)
	pot-=rho2*r/pow(r+r2,2.)*ct*c2p;
	return pot*conv::FPG;
}


VecDoub BowdenNFW::Forces(const VecDoub &x){
	double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	double R2 = (x[0]*x[0]+x[1]*x[1]);
	double iR2 = 1./R2;
	double ct = x[2]/r, st = sqrt(1-ct*ct);
	double dPdr = rho0/r*(-log(1+r/rs)/r+1./(rs+r));
	dPdr+=rho1*(2.*r/(r+r1)-1)/pow(r+r1,2.)*(3*ct*ct-1.);
	double dPdt=rho1*r/pow(r+r1,2)*6.*st*ct;
	double c2p = (R2>0.?(x[0]*x[0]-x[1]*x[1])*iR2:0.); // cos(2\phi)
	double s2p = (R2>0.?2.*x[0]*x[1]*iR2:0.); // sin(2\phi)
	dPdr-=rho2*(2.*r/(r+r2)-1)/pow(r+r2,2.)*st*st*c2p;
	dPdt+=rho2*r/pow(r+r2,2.)*2.*st*ct*c2p;
	double dPdp=-rho2*r/pow(r+r2,2)*st*st*2.*s2p;
	VecDoub f = {0.,0.,0.};
    for(int i=0;i<3;i++) f[i] = dPdr*x[i]/r;
    double ist=1./st;
    f[0]+= (R2>0.?-x[1]*iR2*dPdp:0.)+(st!=0.?x[0]*x[2]/pow(r,3.)*dPdt*ist:0.);
    f[1]+= (R2>0.?x[0]*iR2*dPdp:0.) +(st!=0.?x[1]*x[2]/pow(r,3.)*dPdt*ist:0.);
    f[2]+=-(st!=0.?(1.-ct*ct)/r*dPdt*ist:0.);
	return f*(conv::FPG);
}
// ============================================================================

VecDoub torusPSPT2cartvec(PSPT FF){
	VecDoub X = {FF[0]*cos(FF[2]),FF[0]*sin(FF[2]),FF[1],FF[3]*cos(FF[2])-FF[5]*sin(FF[2]),FF[3]*sin(FF[2])+FF[5]*cos(FF[2]),FF[4]};
	for(int i=3;i<6;++i)X[i]*=conv::kpcMyr2kms;
	return X;
}

double WrapperTorusPotential::operator()(const double R, const double z) const{
	return Pot->Phi({R,0.,z})/conv::kpcMyr2kmsSq;
}

double WrapperTorusPotential::operator()(const double R, const double z, double& dPdR, double& dPdz) const{
	VecDoub F = Pot->Forces({R,0.,z});
	dPdR = -F[0]/conv::kpcMyr2kmsSq;
	dPdz = -F[2]/conv::kpcMyr2kmsSq;
	return Pot->Phi({R,0.,z})/conv::kpcMyr2kmsSq;
}
double WrapperTorusPotential::RfromLc(const double L_in, double* dR) const
{
  bool more=false;
  double R,lR=0.,dlR=0.001,dPR,dPz,LcR,oldL,L=fabs(L_in);
  R=exp(lR);
  (*this)(R,0.,dPR,dPz);
  LcR=sqrt(R*R*R*dPR);
  if(LcR == L) return R;
  if(L>LcR) more=true;
  oldL=LcR;

  for( ; ; ) {
    lR += (more)? dlR : -dlR;
    R=exp(lR);
    (*this)(R,0.,dPR,dPz);
    LcR=sqrt(R*R*R*dPR);
    if(LcR == L) return R;
    if((L< LcR && L>oldL) ||(L>LcR && L<oldL)){
  R=(more)? exp(lR-0.5*dlR) : exp(lR+0.5*dlR);
  return R;}
    oldL=LcR;
  }
}
double WrapperTorusPotential::LfromRc(const double R, double* dR) const
{
  double dPR,dPz;
  (*this)(R,0.,dPR,dPz);
  return sqrt(R*R*R*dPR);
}
Frequencies WrapperTorusPotential::KapNuOm(            // returns kappa,nu,Om
            const double R) const  // given R at z=0
{
  Frequencies epi;
  double dPR,dPz, tmp, dPR2, dPz2;
  (*this)(R,0.,dPR,dPz);
  tmp = dPR/R;
  double delz = 2e-3, delR = 0.01*R;
  (*this)(R+delR,0.,dPR,dPz2);
  (*this)(R-delR,0.,dPR2,dPz2);
  (*this)(R,delz,dPR,dPz2);
  epi[2] = sqrt(tmp);
  epi[1] = sqrt((dPz2-dPz)/delz);
  epi[0] = sqrt(.5*(dPR-dPR2)/delR+3*tmp);
  return epi;
}

// ============================================================================
// potential.cpp
