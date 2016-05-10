// ============================================================================
/// \file GSLInterface/test.cpp
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
/// \brief test file for GSLInterface -- Performs tests on algorithms
///
//============================================================================

#include "GSLInterface.h"
#include <iostream>

//=============================================================================
// Test Functions
double Offset=25.0;

double quad(double x,void*){
	return x*x-Offset;
	}

double cube(double x,void*p){
	double *P = (double *) p;
	return x*x*x+2.0*x*x-5.*x+(*P);
}

int shm(double t, const double y[], double f[],void *params){
  f[0] = y[1];
  f[1] = -y[0];
  return GSL_SUCCESS;
}

double rosenbrock(const gsl_vector *v, void *params){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	return (1.-x)*(1.-x)+100.*(y-x*x)*(y-x*x);
}

double exact = 1.3932039296856768591842462603255;

double g (double *k, size_t dim, void *params){
  	double A = 1.0 / (M_PI * M_PI * M_PI);
 	return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
	// A/=8.0;
	// return 3.7796447*sqrt(A)*exp(-k[0]*k[0]/0.25/2.-k[1]*k[1]/0.7/2.-k[2]*k[2]/0.4/2.);
}

//=============================================================================

int main(){

	// 1. Random Numbers

		// 1.1 Uniform
		unsigned long int Seed = 329567593472;
		rand_uniform rn(Seed);
		std::cout << "Uniform random number: "<< rn.nextnumber() << ", ";

		// 1.2 Gaussian
		double sigma = 5.0;
		rand_gaussian rn2(sigma,Seed+1000);
		std::cout << "Gaussian random number with sigma=5: "<< rn2.nextnumber() << ", ";

		// 1.2.1 Pass new sigma
		sigma = 0.2;
		rn2.newsigma(sigma);
		std::cout << "Gaussian random number with sigma=0.2: "<< rn2.nextnumber() << std::endl
		<< std::endl;

	// 2. Root Finding

		// 2.1 Find root of quad in interval [1,50]
		double xlow = 1.0; double xhigh = 50.; double eps = 1e-8; int maxiter = 100;
		root_find RF(eps,maxiter);
		std::cout << "Root of quad in interval [1,50] = "<< RF.findroot(&quad,xlow,xhigh) << ", ";

		// 2.2 Find root of cube in interval [1,50]
		double a = 1.;
		std::cout << "Root of cube in interval [1,50] = "<< RF.findroot(&cube,xlow,xhigh,&a) << std::endl
		<< std::endl;

	// 3. Integration

		// 3.1 Integrate quad from 0 to 10
		integrator Int(1000);
		xlow = 0.0; xhigh = 10.0;
		std::cout << "Integral of quad from 0 to 10 = " << Int.integrate(&quad,xlow,xhigh,eps) << ", ";

		// 3.2 Integrate cube from 0 to 10
		std::cout << "Integral of cube from 0 to 10 = " << Int.integrate(&cube,xlow,xhigh,eps) <<
		" with error =" << Int.error() << std::endl << std::endl;

		GaussLegendreIntegrator GL(50);
		xlow = 0.0; xhigh = 10.0;
		std::cout << "Gauss Legendre Integral of quad from 0 to 10 = " << GL.integrate(&quad,xlow,xhigh) << std::endl << std::endl;

		// 3.3 Monte Carlo integration of g from [0,0,0] to [pi,pi,pi]
		MCintegrator Int2(3);
		double xl[3]={0.0,0.0,0.0}; double xh[3]={1.,1.,1.}; double err;
		std::cout << "Integral of g from [0,0,0] to [pi,pi,pi] = "<<Int2.integrate(g,xl,xh,100000,&err)<<
		". True result = "<<exact<<"."<<std::endl<<std::endl;

	// 4. Interpolation

		// 4.1 Create spline
		int N = 10;
		double x[N],y[N];
		for(int i=0;i<N;i++){x[i]=(double)i;y[i]=cube(x[i],&a);}
		interpolator interp(x,y,N);

		// 4.2 Interpolate at xtest
		double xtest = .5*(x[3]+x[4]);
		std::cout<<"True value = "<< cube(xtest,&a)<<", Interpolated value = "
		<<interp.interpolate(xtest)<<std::endl << std::endl;

	// 5. Sorting

		// 5.1 Sort list of random numbers
		std::cout<<"List: ";
		for(int i=0;i<N;i++){y[i]=rn.nextnumber();std::cout<<y[i]<<" ";}
		std::cout<<std::endl<<"Sorted list: ";
		sorter S;S.sort(y,N);
		for(int i=0;i<N;i++){std::cout<<y[i]<<" ";x[i]=rn.nextnumber();}
		std::cout<<std::endl<<"Sorted by x: ";

		// 5.2 Sort list of random numbers (x) and arrange second list (y) in same order
		S.sort2(x,N,y);
		for(int i=0;i<N;i++){std::cout<<y[i]<<" ";}
		std::cout<<std::endl << std::endl;

	// 6. ODE

		// 6.1 Setup ode integrator for simple harmonic motion (shm)
		int n = 2;double stepsize = 1e-8;
		ode ODE(shm,n,eps);
		double Initial[2] = {1.,0.};
		double tstart = 0., tfinish = 1.;

		// 6.2 Integrate by 1 time unit
		ODE.step(tstart,tfinish,Initial,stepsize);
		std::cout<<"End values (x,v) = " << Initial[0]<<" "<<Initial[1]<<", True values = "
		<<cos(1.)<<" "<<-sin(1.)<<std::endl;

		// 6.3 Integrate back to starting point
		ODE.step(tfinish,tstart,Initial,-stepsize);
		std::cout<<"Recovered initial values (x,v) = "<<Initial[0]<<" "<<Initial[1]<<std::endl<< std::endl;

	// 7. MultiDim Minimisation

		// 7.1 Setup minimiser for rosenbrock
		int NParam = 2;
		double xy[2] = {2.234,9.123};
		double sizes[2]={1.,1.};
		minimiser M(&rosenbrock,xy,NParam,sizes,eps,NULL);

		// 7.2 Minimise with a max number of 1000 steps and put results in array
		double results[2];
		M.minimise(results,1000,0);
		std::cout<<"Minimum found at ("<<results[0]<<", "<<results[1]<<"), true minimum at (1,1)."
		<< std::endl <<std::endl;

	// 7.1 1D Minimisation

		// 7.1 Setup minimiser for quad
		double bracket[2] = {-10.,10.};
		double initial = 2.;
		int status;
		minimiser1D M1D(&quad,initial,bracket[0],bracket[1],0.,eps,&status,NULL);

		// 7.2 Minimise with a max number of 1000 steps
		double result = M1D.minimise(1000);
		std::cout<<"Minimum of quad found at "<<result<<", true minimum at x=0."
		<< std::endl <<std::endl;

	// 8. Special Functions

		// 8.1 Error Function
		std::cout<<"Erf(0.5) = "<<erf(0.5);
		// 8.2 Bessel Function
		std::cout<<", J_0(0.2) = "<<besselJ(0.2,0);
		// 8.3 Modified Bessel Function
		std::cout<<", I_1(0.3) = "<<besselI(0.3,1);
		// 8.4 Gamma Function
		std::cout<<", Gamma(0.5) = "<<gamma_fn(0.5)<<std::endl;
		// 8.5 Elliptic Integrals
		// 8.5.1 of the First Kind
		std::cout<<"F(0.5,0.1) = "<<ellint_first(0.5,0.1);
		// 8.5.2 of the Second Kind
		std::cout<<", E(0.5,0.1) = "<<ellint_second(0.5,0.1);
		// 8.5.3 of the Third Kind
		std::cout<<", Pi(0.5,0.1,2) = "<<ellint_third(0.5,0.1,2.)<<std::endl;

}
//=============================================================================
