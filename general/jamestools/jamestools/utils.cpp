#include <stdlib.h>
#include "GSLInterface/GSLInterface.h"
rand_uniform *rn;
rand_gaussian *rnGauss;
rand_exponential *rnExp;

double *dmatrix(int n){
	double *m1 = new double[n];
	for(int i=0;i<n;i++) m1[i]=0;
	return m1;
}
int *imatrix(int n){
	int *m1 = new int[n];
	for(int i=0;i<n;i++) m1[i]=0;
	return m1;
}
double **dmatrix(int n,int m){
	double **m1 = new double*[n];
	for(int i=0; i<n; i++) m1[i] = dmatrix(m);
	return m1;
}
double ***dmatrix(int n,int m,int l){
	double ***m1 = new double**[n];
	for(int i=0; i<n; i++) m1[i] = dmatrix(m,l);
	return m1;
}
double ****dmatrix(int n,int m,int l,int k){
	double ****m1 = new double***[n];
	for(int i=0; i<n; i++) m1[i] = dmatrix(m,l,k);
	return m1;
}
void delmatrix(double **m1,int n){
	for(int i=0;i<n;i++) delete[] m1[i];
	delete [] m1;
}
void delmatrix(double ***m1,int n,int m){
	for(int i=0;i<n;i++) delmatrix(m1[i],m);
	delete [] m1;
}
void delmatrix(double ****m1,int n,int m,int l){
	for(int i=0;i<n;i++) delmatrix(m1[i],m,l);
	delete [] m1;
}
