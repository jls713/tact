#ifndef NUMREC_H
#define NUMREC_H

#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))

double bessj0(double);
double bessj1(double);
double bessi0(double);
double bessi1(double);
double bessk0(double);
double bessk1(double);

double qsimp(double (*)(double),double,double);
double qsimpi(double (*)(double),double,double);
double qgaus(double (*)(double),double,double);
void topbottom(double*,int,double,int*,int*);
float linterp(float *,float *,int,float);
double linterp(double *,double *,int,double);
double dlinterp(double *,double *,int,double);
double linterp2d(double **,double *,int,double *,int,double,double);
double linterp3d(double***,double*,int,double*,int,double*,int,double,double,double);
void linterp2(double *,double **,int,double,double *);
void spline(double *,double *,int,double,double,double *);
double splint(double *, double *, double *,int,double);
double splintp(double *,double *,double *,int,double);
void spline2(double *, double *, double**, int, int, double**);
double splint2(double *, double *, double**, double**, int, int, double ,double);
double **matrix(long,long);
float **fmatrix(long,long);
char **smatrix(long,long);
void rk4(double *,double *,int,double,double, void(*)(double,double *,double *,void*),void*p=nullptr);
void rkqs(double *,double *,int,double *,double,double,
		  double *,double *, double *,void (*)(double, double*, double*,void*),void *params = nullptr);
double zbrent(double (*)(double,void*),double,double,double,double,double,int,void* params = nullptr);
int zbrac(double (*)(double),double *,double *,double *,double *,double);
double gammln(double);
double erff(double);
double ran1(long *);
double gasdev(long *);
void sort(unsigned long, double*);
void sort2(unsigned long, double*,double*);
double probks_NR(double);
double ksone(double*,unsigned long,double (*)(double),double*);

void mrqmin(double*,double*,double*,int,double*,int*,int,double**,
	    double**,double*,void(*)(double,double*,double*,double*,int),double*);

void gaussj(double**, int, double**, int);

void chebft(double,double,double*,int,double (*)(double));
double chebev(double,double,double*,int,double);
double polint(double*,double*,int,double,double*);
void polyfit(double,double,double*,int,double*,double*,int);

void four1(float*,unsigned long,int);
void four1(double*,unsigned long,int);
void fourn(float*,unsigned long*,int,int);
void fourn(double*,unsigned long*,int,int);
#endif
