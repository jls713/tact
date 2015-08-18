// jamestools.cpp
// Contains all of JJB's functions which do stuff to files: compress etc.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

namespace iocompress{

#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define ABS(A) ((A)<0?(-A):(A))
#define SIGN(A,B) ((B)>0?(A):(-A))
#define MOD(A,B) ((A)-(B)*((int)(A)/(B)))

float *scratch;
double *dscratch;

void compress(FILE *tmpf,float *x,int npt){

	int i,j,k,lnx,rem,i1 = 90,i2 = 8100,i3 = 729000;
	float at2t23 = 8388608., big=1.e28,small=1.e-28,aln2;
	char m[6];
	m[5]='\0';

	aln2 = log(2.);
	for(k=0 ; k<npt ; k++){
		if (x[k]==0.){
			m[4] = 80;
			m[3] = 0;
			m[2] = 0;
			m[1] = 0;
			m[0] = 0;
		}else{
			if (ABS(x[k])>big) x[k]= SIGN(big, x[k]);
			if (ABS(x[k])<small) x[k] = SIGN(small, x[k]);
			if (ABS(x[k])>=1.)
				lnx = (int)(log(ABS(x[k]))/aln2) + 129;
			else
				lnx = (int)(log(ABS(x[k]))/aln2) + 128;
			j = lnx/i1;
			m[3] = lnx - j*i1;
			lnx = lnx - 129;
			if (lnx>0)
				rem = (int)((ABS(x[k])/pow(2.,(double)(lnx)) - 1.)*at2t23);
			else
				rem = (int)((ABS(x[k])*pow(2.,(double)(ABS(lnx))) - 1.)*at2t23);
			m[4] = rem/i3;
			rem = rem - m[4]*i3;
			m[2] = rem/i2;
			rem = rem - m[2]*i2;
			m[1] = rem/i1;
			m[0] = rem - m[1]*i1;
			m[4] = m[4] + j*12;
			if (x[k]<0.) m[4] = m[4] + 40;
		}
		for (i=0;i<5;i++){
			m[i] = m[i] + 33;
			if (m[i]>=94) m[i] = m[i] + 1;
			if (m[i]>=123) m[i] = m[i] + 1;
		}
		fprintf(tmpf,"%s",m);
		if(15*((k+1)/15)==k+1) fprintf(tmpf,"\n");
	}
	if(15*((npt)/15)!=npt) fprintf(tmpf,"\n");
}

void compres2(FILE *tmpf,float *array,int nx,int ny)
{
	int i,j;
	scratch=(float*)calloc(nx,sizeof(float));
	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++) scratch[j]=array[j+nx*i];
		compress(tmpf,scratch,nx);}
	free(scratch);
	scratch=NULL;
}
void compress(FILE *tmpf,double *x,int npt){

	int i,j,k,lnx,rem,i1 = 90,i2 = 8100,i3 = 729000;
	float at2t23 = 8388608., big=1.e28,small=1.e-28,aln2;
	char m[6];
	m[5]='\0';

	aln2 = log(2.);
	for(k=0 ; k<npt ; k++){
		if (x[k]==0.){
			m[4] = 80;
			m[3] = 0;
			m[2] = 0;
			m[1] = 0;
			m[0] = 0;
		}else{
			if (ABS(x[k])>big) x[k]= SIGN(big, x[k]);
			if (ABS(x[k])<small) x[k] = SIGN(small, x[k]);
			if (ABS(x[k])>=1.)
				lnx = (int)(log(ABS(x[k]))/aln2) + 129;
			else
				lnx = (int)(log(ABS(x[k]))/aln2) + 128;
			j = lnx/i1;
			m[3] = lnx - j*i1;
			lnx = lnx - 129;
			if (lnx>0)
				rem = (int)((ABS(x[k])/pow(2.,(double)(lnx)) - 1.)*at2t23);
			else
				rem = (int)((ABS(x[k])*pow(2.,(double)(ABS(lnx))) - 1.)*at2t23);
			m[4] = rem/i3;
			rem = rem - m[4]*i3;
			m[2] = rem/i2;
			rem = rem - m[2]*i2;
			m[1] = rem/i1;
			m[0] = rem - m[1]*i1;
			m[4] = m[4] + j*12;
			if (x[k]<0.) m[4] = m[4] + 40;
		}
		for (i=0;i<5;i++){
			m[i] = m[i] + 33;
			if (m[i]>=94) m[i] = m[i] + 1;
			if (m[i]>=123) m[i] = m[i] + 1;
		}
		fprintf(tmpf,"%s",m);
		if(15*((k+1)/15)==k+1) fprintf(tmpf,"\n");
	}
	if(15*((npt)/15)!=npt) fprintf(tmpf,"\n");
}

void get(FILE *tmpf,float *xg,int npt)
{
	float expon,at2t23= 8388608.;
	long int i1=90,i,j,k,line,nline;
	unsigned char lin[80],m[5];
	unsigned int np,kp;

	nline=npt/15;
	if(15*nline!=npt) nline=nline+1;
	for(line=0;line<nline;line++){
		fscanf(tmpf,"%s",lin);
		np=strlen((char*)lin)/5;
		if((5*np)!=strlen((char*)lin)) np=np+1;
		for(kp=0;kp<np;kp++){
			k=line*15+kp;
			for(i = 0;i<5;i++){
				m[i] = lin[5*kp+i];
				if (m[i]>=124) m[i] = m[i] - 1;
				if (m[i]>=95) m[i] = m[i] - 1;
				m[i] = m[i] - 33;}
			if (m[4]==80)
				xg[k] = 0.;
			else{
				if (m[4]>=40){
					m[4] = m[4] - 40;
					xg[k] = -1.;}
				else xg[k] = 1.;
				j = m[4]/12;
				m[4] = m[4] - j*12;
				expon=j*i1+m[3]-129;
				xg[k]=xg[k]*((((m[4]*i1+m[2])*i1+m[1])*i1
					      +m[0])/at2t23+ 1.);
				if (expon>0.)
					xg[k]=xg[k]*pow(2.,expon);
				else
					xg[k]=xg[k]/pow(2.,-expon);
			}
		}
	}
}

void get2(FILE *tmpf,float *array,int nx,int ny)
{
	int i,j;
	scratch=(float*)calloc(nx,sizeof(float));
	for(i=0;i<ny;i++){
		get(tmpf,scratch,nx);
		for(j=0;j<nx;j++) array[j+nx*i]=scratch[j];}
	free(scratch);
	scratch=NULL;
}
void compres2(FILE *tmpf,double *array,int nx,int ny)
{
	int i,j;
	dscratch=(double*)calloc(nx,sizeof(double));
	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++) dscratch[j]=array[j+nx*i];
		compress(tmpf,dscratch,nx);}
	free(dscratch);
	dscratch=NULL;
}
void compres2(FILE *tmpf,double **array,int nx,int ny){
	for(int i=0;i<nx;i++) compress(tmpf,array[i],ny);
}
void compress(FILE *tmpf,double **array,int nx,int ny){//same as compres2
	for(int i=0; i<nx;i++) compress(tmpf,array[i],ny);
}
void compress(FILE *tmpf,double ***array,int nx,int ny,int nz){//might be compres3
	for(int i=0; i<nx; i++) compress(tmpf,array[i],ny,nz);
}
void get(FILE *tmpf,double *xg,int npt)
{
	float expon,at2t23= 8388608.;
	long int i1=90,i,j,k,line,nline;
	unsigned char lin[80],m[5];
	unsigned int np,kp;

	nline=npt/15;
	if(15*nline!=npt) nline=nline+1;
	for(line=0;line<nline;line++){
		fscanf(tmpf,"%s",lin);
		np=strlen((char*)lin)/5;
		if((5*np)!=strlen((char*)lin)) np=np+1;
		for(kp=0;kp<np;kp++){
			k=line*15+kp;
			for(i = 0;i<5;i++){
				m[i] = lin[5*kp+i];
				if (m[i]>=124) m[i] = m[i] - 1;
				if (m[i]>=95) m[i] = m[i] - 1;
				m[i] = m[i] - 33;}
			if (m[4]==80)
				xg[k] = 0.;
			else{
				if (m[4]>=40){
					m[4] = m[4] - 40;
					xg[k] = -1.;}
				else xg[k] = 1.;
				j = m[4]/12;
				m[4] = m[4] - j*12;
				expon=j*i1+m[3]-129;
				xg[k]=xg[k]*((((m[4]*i1+m[2])*i1+m[1])*i1
					      +m[0])/at2t23+ 1.);
				if (expon>0.)
					xg[k]=xg[k]*pow(2.,expon);
				else
					xg[k]=xg[k]/pow(2.,-expon);
			}
		}
	}
}

void get2(FILE *tmpf,double *array,int nx,int ny)
{
	int i,j;
	dscratch=(double*)calloc(nx,sizeof(double));
	for(i=0;i<ny;i++){
		get(tmpf,dscratch,nx);
		for(j=0;j<nx;j++) array[j+nx*i]=dscratch[j];}
	free(dscratch);
	dscratch=NULL;
}
void get2(FILE *tmpf,double **array,int nx,int ny){
	for(int i=0; i<nx; i++) get(tmpf,array[i],ny);
}
void get(FILE *tmpf,double **array,int nx,int ny){//same as get2
	for(int i=0; i<nx; i++) get(tmpf,array[i],ny);
}
void get(FILE *tmpf,double ***array,int nx,int ny,int nz){//might be get3
	for(int i=0; i<nx; i++) get(tmpf,array[i],ny,nz);
}

void compress(FILE *tmpf,std::vector<double> x){

	unsigned npt = x.size();

	int i,j,lnx,rem,i1 = 90,i2 = 8100,i3 = 729000;
	float at2t23 = 8388608., big=1.e28,small=1.e-28,aln2;
	char m[6];
	m[5]='\0';

	aln2 = log(2.);
	for(unsigned k=0 ; k<npt ; ++k){
		if (x[k]==0.){
			m[4] = 80;
			m[3] = 0;
			m[2] = 0;
			m[1] = 0;
			m[0] = 0;
		}else{
			if (ABS(x[k])>big) x[k]= SIGN(big, x[k]);
			if (ABS(x[k])<small) x[k] = SIGN(small, x[k]);
			if (ABS(x[k])>=1.)
				lnx = (int)(log(ABS(x[k]))/aln2) + 129;
			else
				lnx = (int)(log(ABS(x[k]))/aln2) + 128;
			j = lnx/i1;
			m[3] = lnx - j*i1;
			lnx = lnx - 129;
			if (lnx>0)
				rem = (int)((ABS(x[k])/pow(2.,(double)(lnx)) - 1.)*at2t23);
			else
				rem = (int)((ABS(x[k])*pow(2.,(double)(ABS(lnx))) - 1.)*at2t23);
			m[4] = rem/i3;
			rem = rem - m[4]*i3;
			m[2] = rem/i2;
			rem = rem - m[2]*i2;
			m[1] = rem/i1;
			m[0] = rem - m[1]*i1;
			m[4] = m[4] + j*12;
			if (x[k]<0.) m[4] = m[4] + 40;
		}
		for (i=0;i<5;i++){
			m[i] = m[i] + 33;
			if (m[i]>=94) m[i] = m[i] + 1;
			if (m[i]>=123) m[i] = m[i] + 1;
		}
		fprintf(tmpf,"%s",m);
		if(15*((k+1)/15)==k+1) fprintf(tmpf,"\n");
	}
	if(15*((npt)/15)!=npt) fprintf(tmpf,"\n");
}

void compress(FILE *tmpf,std::vector<std::vector<double>> xg){
	for(unsigned i=0; i<xg.size(); ++i) compress(tmpf,xg[i]);
}
void compress(FILE *tmpf,std::vector<std::vector<std::vector<double>>> xg){
	for(unsigned i=0; i<xg.size(); ++i) compress(tmpf,xg[i]);
}
void get(FILE *file,std::vector<double>& xg){

	unsigned npt = xg.size();
	float expon,at2t23= 8388608.;
	long int i1=90,i,j,k,line,nline;
	unsigned char lin[80],m[5];
	unsigned int np,kp;
	std::string s;

	nline=npt/15;
	if(15*nline!=npt) nline=nline+1;
	for(line=0;line<nline;line++){
		fscanf(file,"%s",lin);
		np=strlen((char*)lin)/5;
		if((5*np)!=strlen((char*)lin)) np=np+1;
		for(kp=0;kp<np;kp++){
			k=line*15+kp;
			for(i = 0;i<5;i++){
				m[i] = lin[5*kp+i];
				if (m[i]>=124) m[i] = m[i] - 1;
				if (m[i]>=95) m[i] = m[i] - 1;
				m[i] = m[i] - 33;}
			if (m[4]==80)
				xg[k] = 0.;
			else{
				if (m[4]>=40){
					m[4] = m[4] - 40;
					xg[k] = -1.;}
				else xg[k] = 1.;
				j = m[4]/12;
				m[4] = m[4] - j*12;
				expon=j*i1+m[3]-129;
				xg[k]=xg[k]*((((m[4]*i1+m[2])*i1+m[1])*i1
					      +m[0])/at2t23+ 1.);
				if (expon>0.)
					xg[k]=xg[k]*pow(2.,expon);
				else
					xg[k]=xg[k]/pow(2.,-expon);
			}
		}
	}
}

void get(FILE *tmpf,std::vector<std::vector<double>>& xg){
	for(unsigned i=0; i<xg.size(); ++i) get(tmpf,xg[i]);
}

void get(FILE *tmpf,std::vector<std::vector<std::vector<double>>>& xg){
	for(unsigned i=0; i<xg.size(); ++i) get(tmpf,xg[i]);
}

}
