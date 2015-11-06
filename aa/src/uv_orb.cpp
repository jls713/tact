// ============================================================================
/// \file src/uv_orb.cpp
// ============================================================================
/// \author Jason Sanders (and James Binney 2012)
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
/// \brief Axisymmetric Staeckel fudge using Delta estimation from shells
///
/// uv_orb: Wraps axisymmetric Staeckel fudge using the Delta estimation
/// routine from Binney (2014)
///
//============================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface/GSLInterface.h"
#include "gnuplot/gnuplot_i.h"
#include <gsl/gsl_poly.h>
#include "utils.h"
#include "coordsys.h"
#include "coordtransforms.h"
#include "potential.h"
#include "orbit.h"
#include "stackel_aa.h"
#include "get_closed_Rz.h"
#include "uv_orb.h"
#include "debug.h"

// ============================================================================

void uv_orb::readDeltagrids(const std::string& file){
	std::ifstream infile; infile.open(file);
	if(!infile.is_open())std::cerr<<"Problem: "<<file<<" doesn't exist."<<std::endl;
	E_delta = VecDoub(NE,0.);
	L_delta = std::vector<VecDoub>(NE,VecDoub(NL,0.));
	Delta = L_delta;
	for(int i=0;i<NE;++i) infile>>E_delta[i];
	for(int i=0;i<NE;++i)for(int j=0;j<NL;++j) infile>>L_delta[i][j];
	for(int i=0;i<NE;++i)for(int j=0;j<NL;++j) infile>>Delta[i][j];
	infile.close();
}

void uv_orb::fillDeltagrids(const std::string& file){

	E_delta.clear();L_delta.clear();Delta.clear();
    E_delta = create_range(E0,Emax,NE);
    L_delta = std::vector<VecDoub>(NE,VecDoub(NL,0.));
    Delta = std::vector<VecDoub>(NE,VecDoub(NL,0.));
    #pragma omp parallel for schedule (dynamic)
    for(int i=0; i<NE;i++){
        double R = Pot->R_E(E_delta[i]);
        for(double j=0;j<NL;j++){
            L_delta[i][j] = Pot->L_circ(R)*(j*0.8/(double)(NL-1)+0.001);
            find_best_delta DD(Pot, E_delta[i], L_delta[i][j]);
            Delta[i][j]=DD.delta(R*.9);
            if(Delta[i][j]<0. or Delta[i][j]!=Delta[i][j])
            	Delta[i][j]=R/5.;
            if(debug_deltaGrid)
		std::cerr<<"DeltaGrid: "
            	<<OUTPUT(E_delta[i])
            	<<OUTPUT(L_delta[i][j])
            	<<OUTPUT(R)
            	<<OUTPUTE(Delta[i][j]);
        }
    }
    if(debug_deltaGrid) std::cerr<<"DeltaGrid calculated: NE = "<<NE<<", NL = "<<NL<<std::endl;

    std::ofstream outfile; outfile.open(file);
	for(int i=0;i<NE;++i)
		for(int j=0;j<NL;++j)
		outfile<<E_delta[i]<<" "<<L_delta[i][j]<<" "<<Delta[i][j]<<std::endl;
	outfile.close();

}

double uv_orb::findDelta_interp(double E, double L){
    int E_bot=0, E_top=0, L_bot1=0, L_bot2=0, L_top1=0, L_top2=0;
    VecDoub t = L_delta[0];
    double delta1, delta2;
    if(E<E_delta[0]){ E_bot=0; E_top=0;}
    else if(E>E_delta[NE-1]){ E_bot = NE-1; E_top = NE-1;}
    else topbottom<double>(E_delta,E,&E_bot,&E_top,"Find Delta E grid");

    if(L<L_delta[E_bot][0]) delta1 = Delta[E_bot][0];
    else if(L>L_delta[E_bot][NL-1]) delta1 = Delta[E_bot][NL-1];
    else{
        topbottom<double>(L_delta[E_bot],L,&L_bot1,&L_top1);
        delta1 = Delta[E_bot][L_bot1]+(L-L_delta[E_bot][L_bot1])*(Delta[E_bot][L_top1]-Delta[E_bot][L_bot1])/(L_delta[E_bot][L_top1]-L_delta[E_bot][L_bot1]);
    }

    if(L<L_delta[E_top][0]) delta2 = Delta[E_top][0];
    else if(L>L_delta[E_top][NL-1]) delta2 = Delta[E_top][NL-1];
    else{
        topbottom<double>(L_delta[E_top],L,&L_bot2,&L_top2);
        delta2 = Delta[E_top][L_bot2]+(L-L_delta[E_top][L_bot2])*(Delta[E_top][L_top2]-Delta[E_top][L_bot2])/(L_delta[E_top][L_top2]-L_delta[E_top][L_bot2]);
    }

    if(E_bot!=E_top) delta1+=(delta2-delta1)*(E-E_delta[E_bot])/(E_delta[E_top]-E_delta[E_bot]);
    return delta1;
}

VecDoub uv_orb::actions(const VecDoub& x,void*params){

    VecDoub acts(3,0.);
    if(action_check(x,acts,Pot)) return acts;

	double En = Pot->H(x), Lz = Pot->Lz(x);
	if(params!=nullptr){
		double *alphabeta = (double*)params;
		Actions_AxisymmetricStackel_Fudge ATSF(Pot,alphabeta[0]);
		return ATSF.actions(x);
	}
	else{
		double dd = findDelta_interp(En,fabs(Lz));
		Actions_AxisymmetricStackel_Fudge ATSF(Pot,-1.-dd*dd);
		return ATSF.actions(x);
	}
}

VecDoub uv_orb::angles(const VecDoub& x,void *params){

    VecDoub angs(6,0.);
    if(angle_check(x,angs,Pot)) return angs;

	double En = Pot->H(x), Lz = Pot->Lz(x);
	if(params!=nullptr){
		double *alphabeta = (double*)params;
		Actions_AxisymmetricStackel_Fudge ATSF(Pot,alphabeta[0]);
		return ATSF.angles(x);
	}
	else{
		double dd = findDelta_interp(En,fabs(Lz));
		Actions_AxisymmetricStackel_Fudge ATSF(Pot,-1.-dd*dd);
		return ATSF.angles(x);
	}
}

// ============================================================================
