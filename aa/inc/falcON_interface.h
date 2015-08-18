// ============================================================================
// Simple interface to Walter's falcON potential codes
// ============================================================================
#ifndef FALCON_INTERFACE_H
#define FALCON_INTERFACE_H
#include "potential.h"
#include "public/basic.h"
#include "public/types.h"
#include <body.h>
#include <acceleration.h>
#include <externacc.h>

class falcON_Potential :
	public Potential_JS{
	private:
		falcON::nemo_acc *falPot;
	public:

		falcON_Potential(const char*accname,
	                 	 const char*accpars,
	                 	 const char*accfile){
			falPot = new falcON::nemo_acc(accname,accpars,accfile);
		}

		~falcON_Potential(){delete falPot;}

		double Phi(const VecDoub &x){
			double p=0.,m=0.;
			falcON::falcONVec<3,double> xx, vv, aa;
			for(int j=0;j<3;j++){
				xx[j]=x[j];
				vv[j]=0.;
				aa[j]=0.;
			}
			falPot->set(0.,1,&m,&xx,&vv,nullptr,&p,&aa,0);
			return p;
		}
		VecDoub Forces(const VecDoub &x){
			VecDoub a(3,0); double p=0.,m=0.;
			falcON::falcONVec<3,double> xx, vv, aa;
			for(int j=0;j<3;j++){xx[j]=x[j];vv[j]=0.;aa[j]=0.;}
			falPot->set(0.,1,&m,&xx,&vv,nullptr,&p,&aa,0);
			for(int j=0;j<3;j++) a[j]=aa[j];
			return a;
		}
};
class falcON_SphericalPotential :
	public SphericalPotential{
	private:
		falcON::nemo_acc *falPot;
	public:

		falcON_SphericalPotential(const char*accname,
	                 	 const char*accpars,
	                 	 const char*accfile){
			falPot = new falcON::nemo_acc(accname,accpars,accfile);
		}

		~falcON_SphericalPotential(){delete falPot;}

		double Phi_r(double r){
			double p=0.,m=0.;
			falcON::falcONVec<3,double> xx, vv, aa;
			xx[0]=r;xx[1]=0.;xx[2]=0.;
			for(int j=0;j<3;j++){
				vv[j]=0.;
				aa[j]=0.;
			}
			falPot->set(0.,1,&m,&xx,&vv,nullptr,&p,&aa,0);
			return p;
		}
		double dPhi_r(double r){
			VecDoub a(3,0); double p=0.,m=0.;
			falcON::falcONVec<3,double> xx, vv, aa;
			xx[0]=r;xx[1]=0.;xx[2]=0.;
			for(int j=0;j<3;j++){
				vv[j]=0.;
				aa[j]=0.;
			}
			falPot->set(0.,1,&m,&xx,&vv,nullptr,&p,&aa,0);
			for(int j=0;j<3;j++) a[j]=aa[j];
			return -a[0];
		}
};
#endif
// FALCON_INTERFACE_H
// ============================================================================
