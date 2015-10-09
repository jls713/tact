// ============================================================================
/// \file inc/potential.h
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


static double weights4[4] = {0.3626837833783620,0.3137066458778873,0.2223810344533745,0.1012285362903763};
static double abscissa4[4] = {0.1834346424956498,0.5255324099163290,0.7966664774136267,0.9602898564975363};
static double weights5[5] = {0.2955242247147529,0.2692667193099963,0.2190863625159820,0.1494513491505806,0.0666713443086881};
static double abscissa5[5] = {0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717};
static double weights8[8]={0.1894506104550685,0.1826034150449236,0.1691565193950025,0.1495959888165767,0.1246289712555339,0.0951585116824928,0.0622535239386479,0.0271524594117541};
static double abscissa8[8]={0.0950125098376374,0.2816035507792589,0.4580167776572274,0.6178762444026438,0.7554044083550030,0.8656312023878318,0.9445750230732326,0.9894009349916499};
static double abscissa10[10] = {0.0765265211334973337546404,0.2277858511416450780804962,0.3737060887154195606725482,0.5108670019508270980043641,0.6360536807265150254528367,0.7463319064601507926143051,0.8391169718222188233945291,0.9122344282513259058677524,0.9639719272779137912676661,0.9931285991850949247861224};
static double weights10[10] = {0.1527533871307258506980843,0.1491729864726037467878287,0.1420961093183820513292983,0.1316886384491766268984945,0.1181945319615184173123774,0.1019301198172404350367501,0.0832767415767047487247581,0.0626720483341090635695065,0.0406014298003869413310400,0.0176140071391521183118620};
static double abscissa16[16] = {0.0483076656877383,0.1444719615827965,0.2392873622521371,0.3318686022821277,0.4213512761306353,0.5068999089322294,0.5877157572407623,0.6630442669302152,0.7321821187402897,0.7944837959679424,0.8493676137325700,0.8963211557660521,0.9349060759377397,0.9647622555875064,0.9856115115452684,0.9972638618494816};
static double weights16[16]={0.0965400885147278,0.0956387200792749,0.0938443990808046,.0911738786957639,0.0876520930044038,0.0833119242269467,0.0781938957870703,0.0723457941088485,0.0658222227763618,0.0586840934785355,0.0509980592623762,0.0428358980222267,0.0342738629130214,0.0253920653092621,0.0162743947309057,0.0070186100094701};
static double abscissa32[32]={0.0243502926634244325089558,0.0729931217877990394495429,0.1214628192961205544703765,0.1696444204239928180373136,0.2174236437400070841496487,0.2646871622087674163739642,0.3113228719902109561575127,0.3572201583376681159504426,0.4022701579639916036957668,0.4463660172534640879849477,0.4894031457070529574785263,0.5312794640198945456580139,0.5718956462026340342838781,0.6111553551723932502488530,0.6489654712546573398577612,0.6852363130542332425635584,0.7198818501716108268489402,0.7528199072605318966118638,0.7839723589433414076102205,0.8132653151227975597419233,0.8406292962525803627516915,0.8659993981540928197607834,0.8893154459951141058534040,0.9105221370785028057563807,0.9295691721319395758214902,0.9464113748584028160624815,0.9610087996520537189186141,0.9733268277899109637418535,0.9833362538846259569312993,0.9910133714767443207393824,0.9963401167719552793469245,0.9993050417357721394569056};
static double weights32[32]={0.0486909570091397203833654,0.0485754674415034269347991,0.0483447622348029571697695,0.0479993885964583077281262,0.0475401657148303086622822,0.0469681828162100173253263,0.0462847965813144172959532,0.0454916279274181444797710,0.0445905581637565630601347,0.0435837245293234533768279,0.0424735151236535890073398,0.0412625632426235286101563,0.0399537411327203413866569,0.0385501531786156291289625,0.0370551285402400460404151,0.0354722132568823838106931,0.0338051618371416093915655,0.0320579283548515535854675,0.0302346570724024788679741,0.0283396726142594832275113,0.0263774697150546586716918,0.0243527025687108733381776,0.0222701738083832541592983,0.0201348231535302093723403,0.0179517157756973430850453,0.0157260304760247193219660,0.0134630478967186425980608,0.0111681394601311288185905,0.0088467598263639477230309,0.0065044579689783628561174,0.0041470332605624676352875,0.0017832807216964329472961};

static double abscissa50[50]={0.0156289844215430828722167,0.0468716824215916316149239,0.0780685828134366366948174,0.1091892035800611150034260,0.1402031372361139732075146,0.1710800805386032748875324,0.2017898640957359972360489,0.2323024818449739696495100,0.2625881203715034791689293,0.2926171880384719647375559,0.3223603439005291517224766,0.3517885263724217209723438,0.3808729816246299567633625,0.4095852916783015425288684,0.4378974021720315131089780,0.4657816497733580422492166,0.4932107892081909335693088,0.5201580198817630566468157,0.5465970120650941674679943,0.5725019326213811913168704,0.5978474702471787212648065,0.6226088602037077716041908,0.6467619085141292798326303,0.6702830156031410158025870,0.6931491993558019659486479,0.7153381175730564464599671,0.7368280898020207055124277,0.7575981185197071760356680,0.7776279096494954756275514,0.7968978923903144763895729,0.8153892383391762543939888,0.8330838798884008235429158,0.8499645278795912842933626,0.8660146884971646234107400,0.8812186793850184155733168,0.8955616449707269866985210,0.9090295709825296904671263,0.9216092981453339526669513,0.9332885350430795459243337,0.9440558701362559779627747,0.9539007829254917428493369,0.9628136542558155272936593,0.9707857757637063319308979,0.9778093584869182885537811,0.9838775407060570154961002,0.9889843952429917480044187,0.9931249370374434596520099,0.9962951347331251491861317,0.9984919506395958184001634,0.9997137267734412336782285};
static double weights50[50]={0.0312554234538633569476425,0.0312248842548493577323765,0.0311638356962099067838183,0.0310723374275665165878102,0.0309504788504909882340635,0.0307983790311525904277139,0.0306161865839804484964594,0.0304040795264548200165079,0.0301622651051691449190687,0.0298909795933328309168368,0.0295904880599126425117545,0.0292610841106382766201190,0.0289030896011252031348762,0.0285168543223950979909368,0.0281027556591011733176483,0.0276611982207923882942042,0.0271926134465768801364916,0.0266974591835709626603847,0.0261762192395456763423087,0.0256294029102081160756420,0.0250575444815795897037642,0.0244612027079570527199750,0.0238409602659682059625604,0.0231974231852541216224889,0.0225312202563362727017970,0.0218430024162473863139537,0.0211334421125276415426723,0.0204032326462094327668389,0.0196530874944353058653815,0.0188837396133749045529412,0.0180959407221281166643908,0.0172904605683235824393442,0.0164680861761452126431050,0.0156296210775460027239369,0.0147758845274413017688800,0.0139077107037187726879541,0.0130259478929715422855586,0.0121314576629794974077448,0.0112251140231859771172216,0.0103078025748689695857821,0.0093804196536944579514182,0.0084438714696689714026208,0.0074990732554647115788287,0.0065469484508453227641521,0.0055884280038655151572119,0.0046244500634221193510958,0.0036559612013263751823425,0.0026839253715534824194396,0.0017093926535181052395294,0.0007346344905056717304063};

inline double GaussLegendreQuad(double (*func)(double, void*),double a,double b, void *p = NULL, int order = 50){
    // Adapted from numerical recipes routine qgaus
    auto abscissa = abscissa50, weights = weights50;
    if(order==8){abscissa=abscissa8;weights=weights8;}
    if(order==4){abscissa=abscissa4;weights=weights4;}
    if(order!=50  and order!=8 and order!=4)
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
