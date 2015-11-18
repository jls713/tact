// ============================================================================
/// \file utils.h
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
/// \brief Utility functions
//============================================================================

/*======================================*/
/*		  A few utility functions 	  	*/
/*======================================*/
#ifndef UTILS_H
#define UTILS_H

// ============================================================================

#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include "GSLInterface/GSLInterface.h"

// ============================================================================

typedef std::vector<int> VecInt;
typedef std::vector<double> VecDoub;
typedef std::vector<std::vector<double> > MatDoub;
typedef std::vector<std::vector<std::vector<double> > >MatMatDoub;

const double PI = 3.14159265359;
const double TPI = 6.283185307;
const double PIH = 1.57079632679, HPI=PIH;
const double SMALL = 1e-5;
const double TINY = 1e-10;

#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))
#define SIGN(A) ((A)>0.?(1.0):(-1.0))

extern rand_uniform *rn;
extern rand_gaussian *rnGauss;
extern rand_exponential *rnExp;

// ============================================================================
#define OUTPUT(x) #x<<" = "<<(x)<<", "
#define OUTPUTE(x) #x<<" = "<<(x)<<std::endl

inline double sign(double a){
    if(a>0.)return 1.;
    else return -1.;
}

inline int sign(int a){
    if(a>0.)return 1.;
    else return -1.;
}

template<class c>
c Max(const std::vector<c> &a){
	c R = a[0];
	for_each(begin(a),end(a),[&R](c p){if(p>R)R=p;});
	return R;
}

template<class c>
c Min(const std::vector<c> &a){
	c R = a[0];
	for_each(begin(a),end(a),[&R](c p){if(p<R)R=p;});
	return R;
}

template<class c>
c normsq(const std::vector<c> &a){
    c R = 0.;
    for(auto i:a)R+=i*i;
    return R;
}

template<class c>
c norm(const std::vector<c> &a){
    return sqrt(normsq<c>(a));
}

template<class c>
void printVector(const std::vector<c> &a){
	for ( unsigned int i=0;i<a.size();i++)std::cout<<a[i]<<" ";
	std::cout<<std::endl;
}

template<class c>
void printVector_tofile(const std::vector<c> &a, std::ofstream &outputfile){
    if(!outputfile.is_open())
        std::cerr<<"Output file in printVector not open\n";
    for (unsigned int i=0;i<a.size();i++)outputfile<<a[i]<<" ";
}
template<class c>
std::vector<c> string2vec(std::string str) {
  std::vector<c> internal;
  std::istringstream ss(str); // Turn the string into a stream.
  c tmp;
  while(ss>>tmp)
    internal.push_back(tmp);

  return internal;
}

template<class c>
void printMatrix(const std::vector<std::vector<c>> &a){
    for(unsigned i=0;i<a.size();++i){
        for(unsigned j=0;j<a[0].size();++j)
            std::cout<<a[i][j]<<" ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}


template<class c>
void printVector_c(const std::vector<c> &a){
    for ( unsigned int i=0;i<a.size();i++)std::cout<<a[i]<<",";
    std::cout<<std::endl;
}

template<class c>
void print(const c &a){
	std::cout<<a<<std::endl;
}

template<class c>
std::vector<c> concatVectors(const std::vector<c> &a, const std::vector<c> &b){
	unsigned int p = a.size();
	std::vector<c> Result (p+b.size());
	for(unsigned int i=0;i<p;i++)Result[i]=a[i];
	for(unsigned int i=0;i<b.size();i++)Result[i+p]=b[i];
	return Result;
}

template<class c>
std::vector<c> concatVectors(const std::vector<c> &a, const std::vector<c> &b, const std::vector<c> &d){
	return concatVectors(concatVectors(a,b),d);
}

template<class c>
double Mean(const std::vector<c> &a){
	double R=0.;
	for_each(begin(a),end(a),[&R](c p){R+=p;});
	return ((double)R/(double)a.size());
}

template<class c>
double SD(const std::vector<c> &a){
	double R1=0., R2=0.;
	for_each(begin(a),end(a),[&R1,&R2](c p){R1+=p;R2+=p*p;});
	int N=a.size();R1/=(double)N;
	return sqrt(R2/N-R1*R1);
}
template<class c>
double carefulSD(const std::vector<c> &a){
    double R1=0., R2=0., delta; unsigned n;
    for(n=1;n<a.size()+1;++n){
        delta=a[n-1]-R1;
        R1+=delta/n;
        R2+=delta*(a[n-1]-R1);
    }

    return sqrt(R2/(n-1));
}

template<class c>
double RMS(const std::vector<c> &a){
    double R2=0.;
    for_each(begin(a),end(a),[&R2](c p){R2+=p*p;});
    int N=a.size();
    return sqrt(R2/N);
}

template<class c>
double Median(const std::vector<c> &a){
	std::vector<c> a2 = a;
	size_t midIndex = a2.size()/2;
	std::nth_element(a2.begin(), a2.begin() + midIndex, a2.end());
	return a2[midIndex];
}

template<class c>
std::vector<c> rowMean(const std::vector<std::vector<c>> &a){
	std::vector<c> rM(a.size());
	for(unsigned int i=0;i<a.size();i++) rM[i] = Mean<c>(a[i]);
	return rM;
}

template<class c>
std::vector<c> rowSD(const std::vector<std::vector<c>> &a){
	std::vector<c> rM(a.size());
	for(unsigned int i=0;i<a.size();i++) rM[i] = SD<c>(a[i]);
	return rM;
}

template<class c>
std::vector<c> rowcarefulSD(const std::vector<std::vector<c>> &a){
    std::vector<c> rM(a.size());
    for(unsigned int i=0;i<a.size();i++) rM[i] = carefulSD<c>(a[i]);
    return rM;
}

template<class c>
std::vector<c> rowRMS(const std::vector<std::vector<c>> &a){
    std::vector<c> rM(a.size());
    for(unsigned int i=0;i<a.size();i++) rM[i] = RMS<c>(a[i]);
    return rM;
}
template<class c>
std::vector<c> rowMedian(const std::vector<std::vector<c>> &a){
	std::vector<c> rM(a.size());
	for(unsigned int i=0;i<a.size();i++) rM[i] = Median<c>(a[i]);
	return rM;
}

template<class c>
std::vector<std::vector<c> > transpose(const std::vector<std::vector<c>> &a) {
    std::vector<std::vector<c>> result(a[0].size(),std::vector<c>(a.size()));
    for (unsigned int i = 0; i < a[0].size(); i++)
        for (unsigned int j = 0; j < a.size(); j++) {
            result[i][j] = a[j][i];
        }
    return result;
}

template<class c>
std::vector<c> columnMean(const std::vector<std::vector<c>> &a){
	return rowMean<c>(transpose(a));
}

template<class c>
std::vector<c> columnSD(const std::vector<std::vector<c>> &a){
	return rowSD<c>(transpose(a));
}

template<class c>
std::vector<c> columncarefulSD(const std::vector<std::vector<c>> &a){
    return rowcarefulSD<c>(transpose(a));
}

template<class c>
std::vector<c> columnRMS(const std::vector<std::vector<c>> &a){
    return rowRMS<c>(transpose(a));
}

template<class c>
std::vector<c> columnMedian(const std::vector<std::vector<c>> &a){
	return rowMedian<c>(transpose(a));
}


struct cmp_by_first {
  template<typename T>
  bool operator()(const T& x, const T& y) const { return x[0] < y[0]; }
};


inline double det2(double a, double b, double c, double d){
	// |a b|
	// |c d|
	return a*d-b*c;
	}

inline VecDoub operator*(const VecDoub& v, double a){
	VecDoub v2(v);
	for (VecDoub::iterator i = v2.begin(); i != v2.end(); ++i)
		(*i)*=a;
	return v2;
}

inline VecDoub operator+(const VecDoub& v1, const VecDoub& v2){
	if(v1.size()!=v2.size())
		std::cerr<<"Vector lengths incompatible in VecDoub operator +\n";
	VecDoub v(v1);
    for(unsigned int i=0;i<v2.size();i++)
    	v[i]+=v2[i];
   	return v;
}


inline VecDoub operator-(const VecDoub& v1, const VecDoub& v2){
	if(v1.size()!=v2.size())
		std::cerr<<"Vector lengths incompatible in VecDoub operator -\n";
    VecDoub v(v1);
    for(unsigned int i=0;i<v2.size();i++)
    	v[i]-=v2[i];
   	return v;
}

inline double operator*(const VecDoub& v1, const VecDoub& v2){
    if(v1.size()!=v2.size())
        std::cerr<<"Vector lengths incompatible in VecDoub operator *\n";
    double v=0.;
    for(unsigned int i=0;i<v2.size();i++)
        v+=v1[i]*v2[i];
    return v;
}

inline VecDoub operator*(const MatDoub& v1, const VecDoub& v2){
    if(v1[0].size()!=v2.size())
        std::cerr<<"Vec & Matrix sizes incompatible in Matrix operator *\n";
    VecDoub R(v1.size(),0.);
    for(unsigned i=0;i<R.size();++i)
        for(unsigned k=0;k<v1[0].size();++k)
            R[i]+=v1[i][k]*v2[k];
    return R;
}

inline MatDoub operator*(const MatDoub& v1, double a){
    MatDoub R = v1;
    for(unsigned i=0;i<R.size();++i)
        for(unsigned j=0;j<R[0].size();++j)
            R[i][j]*=a;
    return R;
}

inline MatDoub operator+(const MatDoub& v1, const MatDoub& v2){
    assert(v1.size()==v2.size());
    assert(v1[0].size()==v2[0].size());
    MatDoub R = v1;
    for(unsigned i=0;i<R.size();++i)
        for(unsigned j=0;j<R[0].size();++j)
            R[i][j]+=v2[i][j];
    return R;
}

inline MatDoub operator-(const MatDoub& v1, const MatDoub& v2){
    assert(v1.size()==v2.size());
    assert(v1[0].size()==v2[0].size());
    MatDoub R = v1;
    for(unsigned i=0;i<R.size();++i)
        for(unsigned j=0;j<R[0].size();++j)
            R[i][j]-=v2[i][j];
    return R;
}

inline MatDoub operator*(const MatDoub& v1, const MatDoub& v2){
    if(v1[0].size()!=v2.size())
        std::cerr<<"Matrix sizes incompatible in Matrix operator *\n";
    MatDoub R(v1.size(),VecDoub(v2[0].size(),0.));
    for(unsigned i=0;i<R.size();++i)
        for(unsigned j=0;j<R[0].size();++j)
            for(unsigned k=0;k<v1[0].size();++k)
                R[i][j]+=v1[i][k]*v2[k][j];
    return R;
}

inline MatDoub inverse3D(const MatDoub& A){
    assert(A[0].size()==3);
    assert(A.size()==3);
    MatDoub R = A;
    double det = A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2])
                        -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                        +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
    double invdet = 1/det;
    R[0][0] =  (A[1][1]*A[2][2]-A[2][1]*A[1][2])*invdet;
    R[0][1] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invdet;
    R[0][2] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;
    R[1][0] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invdet;
    R[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invdet;
    R[1][2] = -(A[0][0]*A[1][2]-A[1][0]*A[0][2])*invdet;
    R[2][0] =  (A[1][0]*A[2][1]-A[2][0]*A[1][1])*invdet;
    R[2][1] = -(A[0][0]*A[2][1]-A[2][0]*A[0][1])*invdet;
    R[2][2] =  (A[0][0]*A[1][1]-A[1][0]*A[0][1])*invdet;
    return R;
}


template<class c>
inline void topbottom(const std::vector<c> &xm,c x,int *botx,int *topx, std::string name=""){
    int n=xm.size();
    int m, top=n-1, bot=0;
    if((xm[top]-x)*(x-xm[bot])<0.){
        std::cout<<"improper x input to topbottom: "<<xm[bot]<<" "<<xm[top]<<" "<<x<<" "<<name<<std::endl;
        exit(0);
    }
    while(top-bot>1){
        m=(top+bot)/2;
        if((xm[top]-x)*(x-xm[m])>=0) bot=m;
        else top=m;
    }
    *topx=top; *botx=bot;
}

template<class c>
c linterp(const std::vector<c> &x, const std::vector<c> &y, c x0){ // linear interpolation
    int bot, top;
    topbottom<c>(x,x0,&bot,&top);
    return y[bot]+(x0-x[bot])/(x[top]-x[bot])*(y[top]-y[bot]);
}

template<class c>
c quad_extrapolate(c x0, std::vector<c> x, std::vector<c> y){
    c dy21 = y[2]-y[1], dy10 = y[1]-y[0];
    c dx21 = x[2]-x[1], dx10 = x[1]-x[0];
    c dx102= x[1]*x[1]-x[0]*x[0], dx212= x[2]*x[2]-x[1]*x[1];
    c b = (dy21*dx102-dy10*dx212)/(dx21*dx102-dx10*dx212);
    c a = (dy21-b*dx21)/dx212;
    return a*(x0*x0-x[0]*x[0])+b*(x0-x[0])+y[0];
}


template<class c>
void spline(const std::vector<c> & x, const std::vector<c> & y, c yp1, c ypn, std::vector<c> & y2)
{
    int i,k;
    c p,qn,sig,un,*u;

    int n = x.size();

    u=(c*)calloc(n,sizeof(c));
    if (yp1 > 0.99e30)
        y2[0]=u[0]=0.0;
    else {
        y2[0] = -0.5;
        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    for (i=1;i<n-1;i++) {
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)
        qn=un=0.0;
    else {
        qn=0.5;
        un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (k=n-2;k>=0;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];
    free(u);
}

template<class c>
void output_to_tmpfile(std::vector<std::vector<c>> A,std::string f){
	std::ofstream outfile; outfile.open(f);
	for(std::vector<c> i:A){
		for(c j:i)outfile<<j<<" ";outfile<<std::endl;
	}
	outfile.close();
}

template<class c>
std::vector<c> create_range(c min, c max, int length){
    std::vector<c> P(length);
    for(int i=0;i<length;i++) P[i]=min+(max-min)*i/(c)(length-1);
    return P;
}

template<class c>
std::vector<c> create_log_range(c min, c max, int length){
    std::vector<c> P(length);
    c lmm = log(max/min);
    for(int i=0;i<length;i++) P[i]=min*exp(lmm*i/(c)(length-1));
    return P;
}

template<class c>
std::vector<c> cross_product(std::vector<c> &a,std::vector<c> &b){
    assert(a.size()==b.size());
    assert(a.size()==3);
    return {a[1]*b[2]-a[2]*b[1],-a[0]*b[2]+a[2]*b[0],a[0]*b[1]-a[1]*b[0]};
}

template<class c>
c dot_product(std::vector<c> &a,std::vector<c> &b){
    assert(a.size()==b.size());
    c sum = 0.;
    for(unsigned int i=0;i<a.size();i++) sum+=a[i]*b[i];
    return sum;
}

template<class c>
double dot_product_int(std::vector<int> &a,std::vector<c> &b){
    assert(a.size()==b.size());
    c sum = 0.;
    for(unsigned int i=0;i<a.size();i++) sum+=a[i]*b[i];
    return sum;
}


extern const double weights4[4];
extern const double abscissa4[4];
extern const double weights5[5];
extern const double abscissa5[5];
extern const double weights8[8];
extern const double abscissa8[8];
extern const double abscissa10[10];
extern const double weights10[10];
extern const double abscissa16[16];
extern const double weights16[16];
extern const double abscissa32[32];
extern const double weights32[32];
extern const double abscissa50[50];
extern const double weights50[50];

inline double GaussLegendreQuad(double (*func)(double, void*),double a,double b, void *p = NULL, int order = 10){
    // Adapted from numerical recipes routine qgaus
    auto abscissa = abscissa50, weights = weights50;
    if(order==8){abscissa=abscissa8;weights=weights8;}
    if(order==4){abscissa=abscissa4;weights=weights4;}
    if(order==16){abscissa=abscissa16;weights=weights16;}
    if(order==10){abscissa=abscissa10;weights=weights10;}
    if(order!=50  and order!=8 and order!=4 and order!=10 and order!=16)
        std::cerr<<"GaussLegendre Order not available"<<std::endl;
    double xm=0.5*(b+a);
    double xr=0.5*(b-a);
    double ss=0.0;
    double dx=0.;
    for (int j=0; j<order; j++){
        dx=xr*abscissa[j];
        ss=ss+weights[j]*(func(xm+dx,p)+func(xm-dx,p));
    }
    return ss*xr;
}

inline double GaussLegendreQuad_2D(double (*func)(double,double, void*),VecDoub a, VecDoub b, void *p=NULL){
    double rr=0.0,ss=0.,dx,dy;
    double xm=0.5*(b[0]+a[0]);
    double xr=0.5*(b[0]-a[0]);
    double ym=0.5*(b[1]+a[1]);
    double yr=0.5*(b[1]-a[1]);
    for (int j=0; j<8; j++){
        dx=xr*abscissa8[j];
        ss = 0.;
        for (int k=0; k<8; k++){
            dy=yr*abscissa8[k];
            ss=ss+weights8[k]*(func(xm-dx,ym+dy,p)+func(xm-dx,ym-dy,p)+func(xm+dx,ym+dy,p)+func(xm+dx,ym-dy,p));
        }
        rr=rr+weights8[j]*ss;
    }
    return rr*xr*yr;
}

inline VecDoub GaussianQuad_Weights_8(int i){
    return {weights8[i],abscissa8[i]};
}

inline double GFunction(double Foff,double sigF){return exp(-Foff*Foff/(2.0*sigF*sigF))/(sqrt(TPI)*sigF);}


#endif

//============================================================================
