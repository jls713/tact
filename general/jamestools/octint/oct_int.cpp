//for f(x,y,z) oct_int evaluates rho, or rho & rho*ybar, or rho, rho*ybar, rho*y^2bar or rho, rho*ybar, and all three rho*x_i^2bar
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//using namespace std;

#define MAX(A,B) ((A)>(B)?(A):(B))
#define SMALL 1.e-13

#include "oct_int.h"

void eval_node1(double (*fn)(double,void*),node1 parent,double DX,double *I,int lmax ,void *params=NULL){
	double I0,I1,m[2];
	double dx,c0[3];
	double BottomLeft[2],zero=parent.GetLoc();
	int l1=parent.GetLevel()+1;
	dx=.5*parent.Getdx();// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<2;i++){
		BottomLeft[i]=zero+i*dx;
		m[i]=(*fn)(BottomLeft[i]+.5*dx,params);
	}
	node1 kiddie0(BottomLeft[0],l1,DX,c0[0],c0[2],m[0]);
	node1 kiddie1(BottomLeft[1],l1,DX,c0[2],c0[1],m[1]);
	I1=kiddie0.GetI()+kiddie1.GetI();
	if(fabs(2*I0-I1)<SMALL || l1>lmax) *I+=I1*dx;
	else{
		eval_node1(fn,kiddie0,DX,I,lmax,params);
		eval_node1(fn,kiddie1,DX,I,lmax,params);
	}
}
double oct_int(double (*fn)(double,void*),double x0,double x1,int lmax, void *params){
	double I=0, DX=x1-x0;
	node1 domain(x0,0,DX);
	domain.set_values(fn,params);
	eval_node1(fn,domain,DX,&I,lmax, params);
	return I;
}
void eval_node2(double (*fn)(double,double),node2 parent,location2 DIFF,double *I,int lmax){
	double f0,f1,f2,f3,I0,I1,dx,dy,c0[5],m[4];
	location2 BottomLeft[4],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<4;i++){//bottom-left corners & central values of children
		BottomLeft[i].x=zero.x+fix2[i].x*dx;// clockwise from origin
		BottomLeft[i].y=zero.y+fix2[i].y*dy;
		m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy);
	}
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y);   //fn values clockwise from
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y);
	node2 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,c0[4],f3,m[0]);
	node2 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,c0[4],m[1]);
	node2 kiddie2(BottomLeft[2],l1,DIFF,c0[4],f1,c0[2],f2,m[2]);
	node2 kiddie3(BottomLeft[3],l1,DIFF,f3,c0[4],f2,c0[3],m[3]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI();
//	if(l1>3 && fabs(4*I1-I0)<SMALL) I+=I1*dx*dy;
//	if(l1>7) printf("%d %g (%g,%g) ",l1,4*I1-I0,zero.x,zero.y);
	if(fabs(4*I0-I1)<SMALL || l1>lmax) *I+=I1*dx*dy;
	else{
		eval_node2(fn,kiddie0,DIFF,I,lmax); eval_node2(fn,kiddie1,DIFF,I,lmax);
		eval_node2(fn,kiddie2,DIFF,I,lmax); eval_node2(fn,kiddie3,DIFF,I,lmax);
	}
}
double oct_int(double (*fn)(double,double),double x0,double x1,double y0,double y1,int lmax){
	double I=0;
	location2 DIFF={x1-x0,y1-y0};
	node2 domain(x0,y0,0,DIFF);
	domain.set_values(fn);
	eval_node2(fn,domain,DIFF,&I,lmax);
	return I;
}
void eval_node3(double (*fn)(double,double,double),node3 parent,location3 DIFF,double *I,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
		m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	}
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if(fabs(8*I0-I1)<SMALL || l1>lmax) *I+=I1*dx*dy*dz;
	else{
		eval_node3(fn,kiddie0,DIFF,I,lmax); eval_node3(fn,kiddie1,DIFF,I,lmax);
		eval_node3(fn,kiddie2,DIFF,I,lmax); eval_node3(fn,kiddie3,DIFF,I,lmax);
		eval_node3(fn,kiddie4,DIFF,I,lmax); eval_node3(fn,kiddie5,DIFF,I,lmax);
		eval_node3(fn,kiddie6,DIFF,I,lmax); eval_node3(fn,kiddie7,DIFF,I,lmax);
	}
}
void eval_node3(double (*fn)(double,double,double),node3 parent,location3 DIFF,double *I,double *J,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
		m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	}
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if(fabs(8*I0-I1)<SMALL || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
	}else{
		eval_node3(fn,kiddie0,DIFF,I,J,lmax); eval_node3(fn,kiddie1,DIFF,I,J,lmax);
		eval_node3(fn,kiddie2,DIFF,I,J,lmax); eval_node3(fn,kiddie3,DIFF,I,J,lmax);
		eval_node3(fn,kiddie4,DIFF,I,J,lmax); eval_node3(fn,kiddie5,DIFF,I,J,lmax);
		eval_node3(fn,kiddie6,DIFF,I,J,lmax); eval_node3(fn,kiddie7,DIFF,I,J,lmax);
	}
}
void eval_node3(double (*fn)(double,double,double),node3 parent,location3 DIFF,double *I,double *J,double* L,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,L1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
		m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	}
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if(fabs(8*I0-I1)<SMALL || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
		L1=kiddie0.GetIysq()+kiddie1.GetIysq()+kiddie2.GetIysq()+kiddie3.GetIysq()+kiddie4.GetIysq()+
		   kiddie5.GetIysq()+kiddie6.GetIysq()+kiddie7.GetIysq();
		*L+=L1*dx*dy*dz;
	}else{
		eval_node3(fn,kiddie0,DIFF,I,J,L,lmax); eval_node3(fn,kiddie1,DIFF,I,J,L,lmax);
		eval_node3(fn,kiddie2,DIFF,I,J,L,lmax); eval_node3(fn,kiddie3,DIFF,I,J,L,lmax);
		eval_node3(fn,kiddie4,DIFF,I,J,L,lmax); eval_node3(fn,kiddie5,DIFF,I,J,L,lmax);
		eval_node3(fn,kiddie6,DIFF,I,J,L,lmax); eval_node3(fn,kiddie7,DIFF,I,J,L,lmax);
	}
}
void eval_node3(double (*fn)(double,double,double),node3 parent,location3 DIFF,double *I,double *J,double *K,double *L,double *M,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,K1,L1,M1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
		m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	}
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if(fabs(8*I0-I1)<SMALL || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
		K1=kiddie0.GetIxsq()+kiddie1.GetIxsq()+kiddie2.GetIxsq()+kiddie3.GetIxsq()+kiddie4.GetIxsq()+
		   kiddie5.GetIxsq()+kiddie6.GetIxsq()+kiddie7.GetIxsq();
		*K+=K1*dx*dy*dz;
		L1=kiddie0.GetIysq()+kiddie1.GetIysq()+kiddie2.GetIysq()+kiddie3.GetIysq()+kiddie4.GetIysq()+
		   kiddie5.GetIysq()+kiddie6.GetIysq()+kiddie7.GetIysq();
		*L+=L1*dx*dy*dz;
		M1=kiddie0.GetIzsq()+kiddie1.GetIzsq()+kiddie2.GetIzsq()+kiddie3.GetIzsq()+kiddie4.GetIzsq()+
		   kiddie5.GetIzsq()+kiddie6.GetIzsq()+kiddie7.GetIzsq();
		*M+=M1*dx*dy*dz;
	}else{
		eval_node3(fn,kiddie0,DIFF,I,J,K,L,M,lmax); eval_node3(fn,kiddie1,DIFF,I,J,K,L,M,lmax);
		eval_node3(fn,kiddie2,DIFF,I,J,K,L,M,lmax); eval_node3(fn,kiddie3,DIFF,I,J,K,L,M,lmax);
		eval_node3(fn,kiddie4,DIFF,I,J,K,L,M,lmax); eval_node3(fn,kiddie5,DIFF,I,J,K,L,M,lmax);
		eval_node3(fn,kiddie6,DIFF,I,J,K,L,M,lmax); eval_node3(fn,kiddie7,DIFF,I,J,K,L,M,lmax);
	}
}
double oct_int(double (*fn)(double,double,double),double x0,double x1,
		 double y0,double y1,double z0,double z1,int lmax){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(fn);
	eval_node3(fn,domain,DIFF,&I,lmax);
	return I;
}
double oct_int(double (*fn)(double,double,double),double x0,double x1,
		double y0,double y1,double z0,double z1,int lmax,double *J){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(fn);
	eval_node3(fn,domain,DIFF,&I,J,lmax);
	return I;
}
double oct_int(double (*fn)(double,double,double),double x0,double x1,
		double y0,double y1,double z0,double z1,int lmax,double *J,double *L){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0; *L=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(fn);
	eval_node3(fn,domain,DIFF,&I,J,L,lmax);
	return I;
}
double oct_int(double (*fn)(double,double,double),double x0,double x1,
		double y0,double y1,double z0,double z1,int lmax,double *J,double *K,double *L,double *M){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0; *K=0; *L=0; *M=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(fn);
	eval_node3(fn,domain,DIFF,&I,J,K,L,M,lmax);
	return I;
}
