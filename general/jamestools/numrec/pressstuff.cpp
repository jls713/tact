#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "press.h"

using namespace std;

#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define SIGN(A,B) ((B)>=0 ?fabs(A):-fabs(A))

float linterp(float *xm,float *ym,int np,float x){ // linear interpolation
	int n, top=np-1, bot=0;
	if((xm[top]-x)*(x-xm[bot])<0.){
		printf("improper input to linterp: %f %f %f\n",xm[bot],xm[top],x);
		exit(0);}
	while(top-bot>1){
		n=(top+bot)/2;
		if((xm[top]-x)*(x-xm[n])>=0) bot=n;
		else top=n;
	}
	return ym[bot]+(x-xm[bot])/(xm[top]-xm[bot])*(ym[top]-ym[bot]);
}

void topbottom(double *xm,int n,double x,int *botx,int *topx){
	int m, top=n-1, bot=0;
	if((xm[top]-x)*(x-xm[bot])<0.){
		printf("improper x input to topbottom: %f %f %f\n",xm[bot],xm[top],x);
		exit(0);
	}
	while(top-bot>1){
		m=(top+bot)/2;
		if((xm[top]-x)*(x-xm[m])>=0) bot=m;
		else top=m;
	}
	*topx=top; *botx=bot;
}
double linterp(double *xm,double *ym,int np,double x){ // linear interpolation
	int top,bot;
	topbottom(xm,np,x,&bot,&top);
	return ym[bot]+(x-xm[bot])/(xm[top]-xm[bot])*(ym[top]-ym[bot]);
}
double dlinterp(double *xm,double *ym,int np,double x){ // linear interp for derivative
	int top,bot;
	topbottom(xm,np,x,&bot,&top);
	return (ym[top]-ym[bot])/(xm[top]-xm[bot]);
}
void linterp2(double *xm,double **ym,int np,double x,double *y){
// linear interp of 2 dependent variables
	int top,bot;
	topbottom(xm,np,x,&bot,&top);
	y[0]=ym[bot][0]+(x-xm[bot])/(xm[top]-xm[bot])*(ym[top][0]-ym[bot][0]);
	y[1]=ym[bot][1]+(x-xm[bot])/(xm[top]-xm[bot])*(ym[top][1]-ym[bot][1]);
}
double linterp2d(double **f,double *xm,int nx,double *ym,int ny,
		 double x,double y){ // 2d linear interpolation
	int topx,botx,topy,boty;
	topbottom(xm,nx,x,&botx,&topx);
	topbottom(ym,ny,y,&boty,&topy);
	double dx=(x-xm[botx])/(xm[topx]-xm[botx]);
	double dy=(y-ym[boty])/(ym[topy]-ym[boty]);
	double fa=f[botx][boty]*(1-dx)+f[topx][boty]*dx;
	double fb=f[botx][topy]*(1-dx)+f[topx][topy]*dx;
	return fa*(1-dy)+fb*dy;
}
double linterp3d(double ***f,double *xm,int nx,double *ym,int ny,double *zm,
		 int nz,double x,double y,double z){//3d linear interpolation
	int topx,botx,topy,boty,topz,botz;
	topbottom(xm,nx,x,&botx,&topx);
	topbottom(ym,ny,y,&boty,&topy);
	topbottom(zm,nz,z,&botz,&topz);
	double dx=(x-xm[botx])/(xm[topx]-xm[botx]), dxb=1-dx;
	double dy=(y-ym[boty])/(ym[topy]-ym[boty]), dyb=1-dy;
	double dz=(z-zm[botz])/(zm[topz]-zm[botz]), dzb=1-dz;
	return f[botx][boty][botz]*dxb*dyb*dzb + f[botx][boty][topz]*dxb*dyb*dz
		+f[botx][topy][botz]*dxb*dy*dzb + f[botx][topy][topz]*dxb*dy*dz
		+f[topx][boty][botz]*dx*dyb*dzb + f[topx][boty][topz]*dx*dyb*dz
		+f[topx][topy][botz]*dx*dy*dzb + f[topx][topy][topz]*dx*dy*dz;
}
double trapzd(double (*func)(double), double a, double b,double *s,int *it, int n){
	double x,tnm,sum,del;
	int j;
	if (n == 0) {
		*s=0.5*(b-a)*((*func)(a)+(*func)(b));
		(*it)=1;
	} else {
		del=(b-a)/(*it);
		x=a+0.5*del; sum=0;
		for (j=0;j<(*it);j++,x+=del) sum += (*func)(x);
		*s=0.5*(*s+(b-a)*sum/(*it));
		(*it)*=2;
	}
	return *s;
}

#define EPS 1.0e-5
#define JMAX 20

double qsimp(double (*func)(double), double a, double b){
	int j,it;
	double s,st,ost,os;
	ost = os = -1.0e30;
	for (j=0;j<JMAX;j++) {
		s=(4.0*trapzd(func,a,b,&st,&it,j)-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		if (s == 0.0 && os == 0.0 && j > 6) return s;
		os=s;
		ost=st;
	}
	printf("Too many steps in routine qsimp\n");
	return s;
}
#undef EPS
#undef JMAX

double trapzdi(double (*func)(double), double a, double b,double *s,int *it, int n){
	double x,tnm,sum,del;
	int j;
	if (n == 0) {
		*s=0.5*(b-a)*((*func)(a)+(*func)(b));
		(*it)=1;
	} else {
		del=(b-a)/(*it);
		x=a+0.5*del; sum=0;
		for (j=0;j<(*it);j++,x+=del) sum += (*func)(x);
		*s=0.5*(*s+(b-a)*sum/(*it));
		(*it)*=2;
	}
	return *s;
}

#define EPS 1.0e-6
#define JMAX 20

double qsimpi(double (*func)(double), double a, double b){
	int j,it;
	double s,st,ost,os;
	ost = os = -1.0e30;
	for (j=0;j<JMAX;j++) {
		s=(4.0*trapzdi(func,a,b,&st,&it,j)-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		if (s == 0.0 && os == 0.0 && j > 6) return s;
		os=s;
		ost=st;
	}
	printf("Too many steps in routine qsimp\n");
	return s;
}
#undef EPS
#undef JMAX


void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
	int i,k;
	double p,qn,sig,un,*u;

	u=(double*)calloc(n,sizeof(double));
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

double splint(double *xa, double *ya, double *y2a, int np, double x)
{
	int n, khi=np-1, klo=0;
	if((xa[khi]-x)*(x-xa[klo])<0.){
		printf("improper input to splint: %f %f %f\n",xa[klo],xa[khi],x);
		exit(0);}
	while(khi-klo>1){
		n=(khi+klo)/2;
		if((xa[khi]-x)*(x-xa[n])>=0) klo=n;
		else khi=n;
	}
	double h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad xa input to routine splint\n");
	double a=(xa[khi]-x)/h;
	double b=(x-xa[klo])/h;
	return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

double splintp(double *xa, double *ya, double *y2a, int np, double x)
// returns first deriv
{
	int n, khi=np-1, klo=0;
	if((xa[khi]-x)*(x-xa[klo])<0.){
		printf("improper input to splintp: %f %f %f\n",xa[klo],xa[khi],x);
		exit(0);}
	while(khi-klo>1){
		n=(khi+klo)/2;
		if((xa[khi]-x)*(x-xa[n])>=0) klo=n;
		else khi=n;
	}
	double h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad xa input to routine splintp\n");
	double a=(xa[khi]-x)/h;
	double b=(x-xa[klo])/h;
	return (ya[khi]-ya[klo])/h+((3.*b*b-1.)*y2a[khi]-(3.*a*a-1.)*y2a[klo])*h/6.0;
}

double splintp2(double *xa, double *ya, double *y2a, int np, double x)
// returns second deriv
{
	int n, khi=np-1, klo=0;
	if((xa[khi]-x)*(x-xa[klo])<0.){
		printf("improper input to splintp2: %f %f %f\n",xa[klo],xa[khi],x);
		exit(0);}
	while(khi-klo>1){
		n=(khi+klo)/2;
		if((xa[khi]-x)*(x-xa[n])>=0) klo=n;
		else khi=n;
	}
	double h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad xa input to routine splintp\n");
	double a=(xa[khi]-x)/h;
	double b=(x-xa[klo])/h;
	return y2a[klo]*a + y2a[khi]*b;
}

void spline2(double *x1a, double *x2a, double **ya, int m, int n, double **y2a)
{
	for(int j=0;j<m;j++)spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
}

double splint2(double *x1a, double *x2a, double **ya, double **y2a, int m, int n, double x1, double x2)
{
	double *ytmp=new double[m];for(int i=0;i<m;i++)ytmp[i]=0;
	double *yytmp=new double[m];for(int i=0;i<m;i++)yytmp[i]=0;
	for (int j=0;j<m;j++)yytmp[j]=splint(x2a,ya[j],y2a[j],n,x2);
	spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
	double y = splint(x1a,yytmp,ytmp,m,x1);
	delete ytmp; delete yytmp;
	return y;
}

double **matrix(long nrow, long ncol)
/* allocate a matrix*/
{
	long i;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow)*sizeof(double*)));
	if (!m) printf("allocation failure 1 in strvector()");

	/* allocate rows and set pointers to them */
	m[0]=(double *) malloc((size_t)((nrow*ncol)*sizeof(double)));
	if (!m[0]) printf("allocation failure 2 in strvector()");

	for(i=1;i<nrow;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **fmatrix(long nrow, long ncol)
/* allocate a float vector that contains nrow strings of length ncol*/
{
	long i;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow)*sizeof(float*)));
	if (!m) printf("allocation failure 1 in strvector()");

	/* allocate rows and set pointers to them */
	m[0]=(float *) malloc((size_t)((nrow*ncol)*sizeof(float)));
	if (!m[0]) printf("allocation failure 2 in fmatrix()");

	for(i=1;i<nrow;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

char **smatrix(long nrow, long ncol)
/* allocate a float vector that contains nrow strings of length ncol*/
{
	long i;
	char **m;

	/* allocate pointers to rows */
	m=(char **) malloc((size_t)((nrow)*sizeof(char*)));
	if (!m) printf("allocation failure 1 in smatrix()");

	/* allocate rows and set pointers to them */
	m[0]=(char *) malloc((size_t)((nrow*ncol)*sizeof(char)));
	if (!m[0]) printf("allocation failure 2 in smatrix()");

	for(i=1;i<nrow;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void rk4(double y[], double dydx[],int n, double x, double h,
		 void (*derivs)(double, double [], double [],void*),void*p)
{
	int i;
	double xh,hh,h6,dym[6],dyt[6],yt[6];
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt,p);
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym,p);
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt,p);
	for (i=0;i<n;i++)
		y[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

void rkck(double y[], double dydx[],int N, double x, double h, double yout[],
		  double yerr[], void (*derivs)(double, double [], double [],void*),void *params = NULL)
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
	b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
	b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
	b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
	b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
	c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
	dc5 = -277.00/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
	dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double ak2[6],ak3[6],ak4[6],ak5[6],ak6[6],ytemp[6];

	for (i=0;i<N;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,ak2,params);
	for (i=0;i<N;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,ak3,params);
	for (i=0;i<N;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,ak4,params);
	for (i=0;i<N;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,ak5,params);
	for (i=0;i<N;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)(x+a6*h,ytemp,ak6,params);
	for (i=0;i<N;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=0;i<N;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqs(double *y,double *dydx,int N,double *x,double htry,double eps,
		  double *yscal, double *hdid, double *hnext,
		  void (*derivs)(double, double*, double*,void*), void *params)
{
	int i;
	double errmax,h,htemp,xnew,yerr[6],ytemp[6];
	h=htry;
	for (;;) {
		rkck(y,dydx,N,*x,h,ytemp,yerr,derivs,params);
		errmax=0.0;
		for (i=0;i<N;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x)
		{
			cerr << "At x=" << xnew;
			cerr << "stepsize underflow in rkqs";
		}
	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=0;i<N;i++) y[i]=ytemp[i];
}

#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

#define EPS 3.0e-8

double zbrent(double (*func)(double,void*), double x1, double x2,double fa,double fb,double tol,int itmax,void* params){
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
		printf("Root must be bracketed in zbrent:\nx1,x2,f1,f2: %g %g %g %g",x1,x2,fa,fb);
		exit(0);}
	fc=fb;
	for (iter=1;iter<=itmax;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a; fc=fa; e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b; b=c; c=a; fa=fb; fb=fc; fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s; q=1.0-s;
			} else {
				q=fa/fc; r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d; d=p/q;
			} else {
				d=xm; e=d;
			}
		} else {
			d=xm; e=d;
		}
		a=b; fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(b,params);
	}
//	printf("Too many iters in zbrent: %g %g ",b,fb);
//	exit(0);
	return b;
}

#undef EPS

#define FACTOR 1.2
#define NTRY 25

int zbrac(double (*func)(double), double *x1, double *x2,double *f1,double *f2,double xmax)
{//modified so MAX(*x1,*x2)<xmax
	int j=2;
	double fall[NTRY],xall[NTRY],top=1.e20,bot=-1.e20;

	if (*x1 == *x2 || MAX(*x1,*x2)>xmax){
		printf("Bad initial range in zbrac"); exit(0);}
	*f1=(*func)(*x1); fall[0]=*f1; xall[0]=*x1;
	*f2=(*func)(*x2); fall[1]=*f2; xall[1]=*x2;
	while((*f1)*(*f2) > 0 && j<NTRY) {
		if (fabs(*f1) < fabs(*f2)){
			*x1 += MIN(FACTOR*(*x1-*x2),.95*(xmax-*x1));
			*f1=(*func)(*x1); fall[j]=*f1; xall[j]=*x1;
		}else{
			*x2 += MIN(FACTOR*(*x2-*x1),.95*(xmax-*x2));
			*f2=(*func)(*x2); fall[j]=*f2; xall[j]=*x2;
		}
		j++;
	}
	for(int i=0;i<j;i++){
		if(fall[i]>=0 && fall[i]<top){top=fall[i]; *f1=fall[i], *x1=xall[i];}
		if(fall[i] <0 && fall[i]>bot){bot=fall[i]; *f2=fall[i], *x2=xall[i];}
	}
	if(top*bot<=0) return 1;
	return 0;
}
#undef FACTOR
#undef NTRY

double gammln(double xx)//returns ln(Gamma(xx)) for xx>0
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	24.01409824083091,-1.231739572450155,
	0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#define ITMAX 100
#define EPS 3.0e-7
void gser(double *gamser, double a, double x, double *gln){
	//needed for incomplete gamma fn
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) printf("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		printf("a too large, ITMAX too small in routine gser");
		return;
	}
}

#define FPMIN 1.0e-30
void gcf(double *gammcf, double a, double x, double *gln){
	//needed for incomplete gamma fn
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) printf("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN

double gammp(double a, double x){
	//Returns incomplete Gamma fn P(a,x)=1/Gamma(a)\int_0^x dt e^{-t}t^{a-1}
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) printf("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

double gammq(double a, double x){
	//returns complement of incomplete gamma: Q(a,x)=1-P(a,x)
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) printf("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

double erff(double x){
	//returns error fn erf(x)
	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
double gasdev(long *idum)
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp
#define M 7
#define NSTACK 50

void sort(unsigned long n, double *arr)
{
	unsigned long i,ir=n,j,k,l=1;
	int jstack=0,*istack;
	double a,temp;
	arr-=1;

	istack=(int*)calloc(NSTACK,sizeof(int));
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1]);
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir]);
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir]);
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l]);
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) printf("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free(istack);
}
void sort2(unsigned long n, double *arr, double *brr)
{
	unsigned long i,ir=n,j,k,l=1;
	int jstack=0,*istack;
	double a,b,temp;
	arr-=1;brr-=1;

	istack=(int*)calloc(NSTACK,sizeof(int));
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];b=brr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];brr[i+1]=brr[i];
				}
				arr[i+1]=a;brr[i+1]=b;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1]);SWAP(brr[k],brr[l+1]);
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir]);SWAP(brr[l+1],brr[ir]);
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir]);SWAP(brr[l],brr[ir]);
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l]);SWAP(brr[l+1],brr[l]);
			}
			i=l+1;
			j=ir;
			a=arr[l];b=brr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);SWAP(brr[i],brr[j]);
			}
			arr[l]=arr[j];brr[l]=brr[j];
			arr[j]=a;brr[j]=b;
			jstack += 2;
			if (jstack > NSTACK) printf("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free(istack);
}

#undef M
#undef NSTACK
#undef SWAP

#define EPS1 0.001
#define EPS2 1.0e-8
double probks(double alam){
	int j;
	double a2,fac=2.0,sum=0.0,term,termbf=0.0;
	a2 = -2.0*alam*alam;
	for (j=1;j<=100;j++) {
		term=fac*exp(a2*j*j);
		sum += term;
		if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
		fac = -fac;
		termbf=fabs(term);
	}
	return 1.0;
}
#undef EPS1
#undef EPS2

double ksone(double *data, unsigned long n, double (*func)(double), double *d){
	double probks(double alam);
	void sort(unsigned long n, double *arr);
	unsigned long j;
	double dt,en,ff,fn,fo=0.0;

	sort(n,data);
	en=n;
	*d=0.0;
	for (j=0;j<n;j++) {
		fn=j/en;
		ff=(*func)(data[j]);
		dt=MAX(fabs(fo-ff),fabs(fn-ff));
		if (dt > *d) *d=dt;
		fo=fn;
	}
	en=sqrt(en);
	return probks((en+0.12+0.11/en)*(*d));
}

double kstwo(double *data1, unsigned long n1, double *data2, unsigned long n2, double *d){
	double probks(double alam);
	void sort(unsigned long n, double *arr);
	unsigned long j1=0,j2=0;
	double dt,en1,en2,d1,d2,fn1,fn2;

	sort(n1,data1);sort(n2,data2);
	en1=n1;en2=n2;
	*d=0.0;
	while(j1<n1&&j2<n2){
		if((d1=data1[j1])<=(d2=data2[j2]))fn1=j1++/en1;
		if(d2<=d1)fn2=j2++/en2;
		if((dt=fabs(fn2-fn1))>*d)*d=dt;
	}
	double en=sqrt(en1*en2/(en1+en2));
	return probks((en+0.12+0.11/en)*(*d));
}

double bessj0(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
						+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
					   +y*(-0.6911147651e-5+y*(0.7621095161e-6
			-y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}
double bessj1(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
			+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
			+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
					   +y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
				      +y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}
double bessi0(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}
double bessi1(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	} else {
		y=3.75/ax;
		ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2));
		ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		ans *= (exp(ax)/sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}
double bessk0(double x)
{
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
			+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
			+y*(0.10750e-3+y*0.74e-5))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
			+y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return ans;
}
double bessk1(double x)
{
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
	return ans;
}


double qgaus(double (*func)(double),double a,double b){
	double x[5] = {.1488743389,.4333953941,.6794095682,.8650633666,.9739065285};
	double w[5] = {.2955242247,.2692667193,.2190863625,.1494513491,.0666713443};
	double xm=0.5*(b+a);
	double xr=0.5*(b-a);
	double ss=0.0;
	for (int j=0; j<5; j++){
		double dx=xr*x[j];
		ss=ss+w[j]*(func(xm+dx)+func(xm-dx));
	}
	return ss*xr;
}

