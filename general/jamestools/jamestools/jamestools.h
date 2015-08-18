#ifndef JAMESTOOLS_H
#define JAMESTOOLS_H

#include <stdio.h>
#include <fstream>
#include <vector>
//============================================================================
/// James' IO tools
//============================================================================

namespace iocompress{
	void compress(FILE *,float *,int);
	void compress(FILE *,double *,int);
	void compress(FILE *,double **,int,int);
	void compress(FILE *,double ***,int,int,int);
	void compres2(FILE *,float *,int,int);
	void compres2(FILE *,double *,int,int);
	void compres2(FILE *,double **,int,int);
	void get(FILE *,float *,int);
	void get(FILE *,double *,int);
	void get(FILE *,double **,int,int);
	void get(FILE *,double ***,int,int,int);
	void get2(FILE *,float *,int,int);
	void get2(FILE *,double *,int,int);
	void get2(FILE *,double **,int,int);
	// new c++ vector
	void compress(FILE *tmpf,std::vector<double> xg);
	void compress(FILE *tmpf,std::vector<std::vector<double> >  xg);
	void compress(FILE *tmpf,std::vector<std::vector<std::vector<double>>> xg);
	void get(FILE *tmpf,std::vector<double>& xg);
	void get(FILE *tmpf,std::vector<std::vector<double>>& xg);
	void get(FILE *tmpf,std::vector<std::vector<std::vector<double>>>& xg);
}
//============================================================================

#endif
