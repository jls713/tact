// ============================================================================
/// \file src/tables_aa.cpp
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
/// \brief Interpolation grids for axisymmetric Staeckel fudge
///
/// Actions_AxisymmetricFudge_InterpTables: Constructs grids for interpolating
/// the results of the axisymmetric Staeckel fudge
///
//============================================================================

#include "utils.h"
#include "potential.h"
#include "uv_orb.h"
#include "jamestools/jamestools/jamestools.h"
#include "coordsys.h"
#include "GSLInterface/GSLInterface.h"
#include "tables_aa.h"

//============================================================================
/// Interpolation tables for Axisymmetric fudge
//============================================================================

double Actions_AxisymmetricFudge_InterpTables::Phiuv(const VecDoub& uv,double Delta){
    UVProlateSpheroidCoordSys UVV(Delta);
    VecDoub Rz = UVV.uv2Rz(uv);
    return Pot->Phi({Rz[0],0.,Rz[1]});
}
double Actions_AxisymmetricFudge_InterpTables::dU(double u, double v, double sv2, double Phi1, double sh1sq, double Delta){
    //returns U(u)-U(u0)
    return (pow(sinh(u),2)+sv2)*Phiuv({u,v},Delta)-(sh1sq+sv2)*Phi1;
}

double Actions_AxisymmetricFudge_InterpTables::dV(double v,double u,double shu2,double Delta){//returns V(v)-V(PIH)
    return (1+shu2)*Phiuv({u,PIH},Delta)-(shu2+pow(sin(v),2))*Phiuv({u,v},Delta);
}

double Actions_AxisymmetricFudge_InterpTables::dPhiu(const VecDoub& uv,double *f,double Delta){
    //returns Phiu(u) and dPhi/dR etc
    UVProlateSpheroidCoordSys UVV(Delta);
    VecDoub Rz = UVV.uv2Rz(uv);
    VecDoub F = Pot->Forces({Rz[0],0.,Rz[1]});
    f[0]=-F[0];f[1]=-F[2];
    return Pot->Phi({Rz[0],0.,Rz[1]});
}

static double duturnfn(double u,void*p){//derivative of uturnfn
    uturn_st *P = (uturn_st *) p;
    double shu=sinh(u),chu=cosh(u),shu2=shu*shu,f[2];
    double Phi0=P->AAFT->dPhiu({u,P->v},f,P->Delta);
    double dPhidu=P->Delta*(f[0]*chu*sqrt(P->sv2)+f[1]*shu*sqrt(1-P->sv2));
    return 2*(P->E-Phi0)*shu*chu-(shu2+P->sv2)*dPhidu+P->Lzsq*chu/(P->Delta2*shu2*shu);
}

static double negduturnfn(double u,void*p){return -duturnfn(u,p);}

double Actions_AxisymmetricFudge_InterpTables::Getu0(double utry,double E,double Delta,double Iu,double Lzsq,double sh1sq, double v, double sv2){
    uturn_st UVV(this,E,Delta,Iu,Lzsq,sh1sq,v,sv2);
    double ui,uo=0.,dut=duturnfn(utry,&UVV);
    if(dut>0){//inside
        while(dut>0){
            ui=utry; utry*=1.2; dut=duturnfn(utry,&UVV);
        } uo=utry;
    }else{//outside
        while(dut<0){
            uo=utry; utry*=.8; dut=duturnfn(utry,&UVV);
        } ui=utry;
    }
    if(ui<TINY_UV and uo<TINY_UV){
        int status;
        minimiser1D min(&negduturnfn,1e-1,1e-5,100.,10.,1e-3,&status,&UVV);
        return min.minimise(100);
    }
    if(dut!=dut or std::isinf(dut) or std::isnan(dut)){
        std::cerr<<"duturnfn is nan in Getu0\n";
        return us;
    }
    root_find RF(TINY_UV,100);
    return RF.findroot(&duturnfn,ui,uo,&UVV);
}

double Actions_AxisymmetricFudge_InterpTables::find_umid(VecDoub X,double Delta,double E,double Lzsq){
    CartesianToUVCoords PT(X,Delta);
    double sh1sq=pow(sinh(us),2), Phiu1=Phiuv({us,PT.v},Delta);
    double IuS=E*(PT.shu2-sh1sq)-.5*pow(PT.pu,2)/PT.Delta2-.5*Lzsq/PT.Delta2*(1/PT.shu2-1/sh1sq)-dU(PT.u,PT.v,PT.sv2,Phiu1,sh1sq,Delta);
    return Getu0(us,E,Delta,IuS,Lzsq,sh1sq,PT.v,PT.sv2);
}

VecDoub Actions_AxisymmetricFudge_InterpTables::IuIv(VecDoub X,double Delta,double E,double Lzsq){
    CartesianToUVCoords PT(X,Delta);
    double sh1sq=pow(sinh(us),2), Phiu1=Phiuv({us,PT.v},Delta);
    double IuS=E*(PT.shu2-sh1sq)-.5*pow(PT.pu,2)/PT.Delta2-.5*Lzsq/PT.Delta2*(1/PT.shu2-1/sh1sq)-dU(PT.u,PT.v,PT.sv2,Phiu1,sh1sq,Delta);
    double umid = Getu0(us,E,Delta,IuS,Lzsq,sh1sq,PT.v,PT.sv2);
    sh1sq=pow(sinh(umid),2); Phiu1=Phiuv({umid,PT.v},Delta);
    IuS=E*(PT.shu2-sh1sq)-.5*pow(PT.pu,2)/PT.Delta2-.5*Lzsq/PT.Delta2*(1/PT.shu2-1/sh1sq)-dU(PT.u,PT.v,PT.sv2,Phiu1,sh1sq,Delta);
    IuS=IuS/pow(cosh(umid),2);
    double IvS=.5*pow(PT.pv,2)/PT.Delta2-E*PT.sv2+.5*Lzsq/(PT.Delta2*PT.sv2)-dV(PT.v,PT.u,PT.shu2,Delta);
    IvS=(IvS+E-.5*Lzsq/PT.Delta2)/pow(cosh(umid),2);
    return {IuS,IvS};
}

Actions_AxisymmetricFudge_InterpTables::Actions_AxisymmetricFudge_InterpTables(Potential_JS *Pot, std::string name, bool tab, double Rm, double Rn,int NR, int NGRID,int NED,int NEL)
    :Pot(Pot),name(name)
    ,Rmin(Rm),Rmax(Rn),NR(NR),NGRID(NGRID),UV(new uv_orb(Pot,Rm,Rn,NED,NEL)){

    no_table=false;

    Rgrid = VecDoub(NR); Lzgrid = VecDoub(NR); Vmax = VecDoub(NR);
    Lcgrid = VecDoub(NR);
    E0grid = VecDoub(NR); Unigrid = VecDoub(NGRID+1);
    Iu0 = VecDoub(NR); Iv0 = VecDoub(NR);
    kappagrid = VecDoub(NR); nugrid = VecDoub(NR); Omegacgrid = VecDoub(NR);
    Jr_LzEIu=std::vector<std::vector<VecDoub>>(NR,std::vector<VecDoub>(NGRID,VecDoub(NGRID)));
    Iugrid=Jr_LzEIu; Jz_LzEIv=Jr_LzEIu; Ivgrid=Jr_LzEIu;
    Iumin=std::vector<VecDoub>(NR,VecDoub(NGRID)); Iumax=Iumin;
    Ivmin=std::vector<VecDoub>(NR,VecDoub(NGRID)); Ivmax=Ivmin;

    for(int i=0; i<=NGRID; i++)
        Unigrid[i]=((double)i/(double)(NGRID));//   standard grid
    for(int i=0; i<NR; i++)
        Rgrid[i]=Rmin+i*(Rmax-Rmin)/((double)(NR-1));

    if(tab) tab_acts();
    else get_acts();
}

void Actions_AxisymmetricFudge_InterpTables::fillLzEVmaxgrids(bool refill_freq_grid
    // option to not refill freq grid for adiabatic contraction
    ){
    double R = 0.;
    for(int i=0; i<NR; i++){
        R = Rgrid[i];
        Lzgrid[i]=Pot->L_circ(R);
        E0grid[i]=Pot->E_circ(R);
        Vmax[i]=sqrt(2*(Pot->Phi({4.*Rmax,0.,4.*Rmax})-E0grid[i]));//Range of Egrid at R
        if(refill_freq_grid){
            VecDoub f = Pot->freqs(R);
            kappagrid[i]=f[0];Omegacgrid[i]=f[2];nugrid[i]=f[1];
            Lcgrid[i]=Lzgrid[i];
        }
    }
}

void Actions_AxisymmetricFudge_InterpTables::tab_acts(bool refill_freq_grid
    // option to not refill freq grid for adiabatic contraction
    ){
    FILE *acts_file;
    acts_file = fopen(name.c_str(),"w");
    if(acts_file==nullptr){
        std::cout<<"I cannot open "<<name<<"\n";
        exit(0);
    }

    fillLzEVmaxgrids(refill_freq_grid);

    double E,speed,vc2,R,Lz,Phi0,Delta,alpha;
    VecDoub X(6,0.),Iuv,Acts;X[3]=SMALL;
    for(int nL=0; nL<NR; nL++){
        Lz=Lzgrid[nL]; R=Rgrid[nL];
        for(int nE=-1; nE<NGRID; nE++){
            E=E0grid[nL]+.5*pow(Vmax[nL]*(nE+1)/(double)(NGRID),2);
            Phi0=Pot->Phi({R,0.,0.}); vc2=pow(Lz/R,2);
            speed=sqrt(MAX(0,2*(E-Phi0)-vc2));
            Delta = UV->findDelta_interp(E,Lz);
            alpha = -1-Delta*Delta;
            // find minimum of u potential
            double utry=asinh(R/Delta);
            X[0]=R; X[3]=0; X[4] = Lz/R; X[5]=speed;
            double umid = find_umid(X,Delta,E,Lz*Lz);
            while(fabs(utry-umid)>SMALL){
                X[0]=R=Delta*sinh(umid);
                Phi0=Pot->Phi({R,0.,0.});
                vc2=pow(Lz/R,2);
                X[4]=Lz/R;
                X[5]=speed=sqrt(MAX(0,2*(E-Phi0)-vc2));
                utry=umid;
            }
            if(nE<0){//At this E only circ orbit
                Iuv = IuIv(X,Delta,E,Lz*Lz);
                Iu0[nL]=Iuv[0]; Iv0[nL]=Iuv[1];
            }else{//Now fire off at each angle to plane
                for(int nz=0; nz<NGRID; nz++){
                    double theta=1e-2+(double)nz*(PIH-2e-2)/(double)(NGRID-1);
                    X[3]=speed*cos(theta); X[5]=speed*sin(theta);
                    Iuv = IuIv(X,Delta,E,Lz*Lz);
                    Acts = UV->actions(X,&alpha);
                    Jr_LzEIu[nL][nE][nz]=Acts[0];
                    Iugrid[nL][nE][nz]=Iuv[0];
                    Jz_LzEIv[nL][nE][nz]=Acts[2];
                    Ivgrid[nL][nE][nz]=Iuv[1];
                }
                Iumin[nL][nE]=Iugrid[nL][nE][0];
                Iumax[nL][nE]=Iugrid[nL][nE][NGRID-1];
                Ivmin[nL][nE]=Ivgrid[nL][nE][0];
                Ivmax[nL][nE]=Ivgrid[nL][nE][NGRID-1];
                for(int nz=0; nz<NGRID; nz++){//renormalise
                    Iugrid[nL][nE][nz]-=Iumin[nL][nE];
                    Iugrid[nL][nE][nz]/=(Iumax[nL][nE]-Iumin[nL][nE]);
                    Ivgrid[nL][nE][nz]-=Ivmin[nL][nE];
                    Ivgrid[nL][nE][nz]/=(Ivmax[nL][nE]-Ivmin[nL][nE]);
                }
            }
        }
    }
    printf("\n");
    iocompress::compress(acts_file,Iu0);
    iocompress::compress(acts_file,Iugrid);
    iocompress::compress(acts_file,Iumin);
    iocompress::compress(acts_file,Iumax);
    iocompress::compress(acts_file,Jr_LzEIu);
    iocompress::compress(acts_file,Iv0);
    iocompress::compress(acts_file,Ivgrid);
    iocompress::compress(acts_file,Ivmin);
    iocompress::compress(acts_file,Ivmax);
    iocompress::compress(acts_file,Jz_LzEIv);
    fclose(acts_file);
    std::cerr<<"Interpolation grids calculated and stored for NR = "<<NR<<
        ", NGRID = "<<NGRID<<std::endl;
}

void Actions_AxisymmetricFudge_InterpTables::get_acts(void){//read in previously tabulated values

    fillLzEVmaxgrids();

    FILE *acts_file;
    acts_file = fopen(name.c_str(),"r");
    if(acts_file==nullptr){
        std::cout<<"I cannot open "<<name<<"\n";
        exit(0);
    }
    iocompress::get(acts_file,Iu0);
    iocompress::get(acts_file,Iugrid);
    iocompress::get(acts_file,Iumin);
    iocompress::get(acts_file,Iumax);
    iocompress::get(acts_file,Jr_LzEIu);
    iocompress::get(acts_file,Iv0);
    iocompress::get(acts_file,Ivgrid);
    iocompress::get(acts_file,Ivmin);
    iocompress::get(acts_file,Ivmax);
    iocompress::get(acts_file,Jz_LzEIv);
    fclose(acts_file);
    std::cerr<<"Interpolation grids read in for NR = "<<NR<<
        ", NGRID = "<<NGRID<<std::endl;

}

int Actions_AxisymmetricFudge_InterpTables::int_acts(double Lz,double E,double Iu,double Iv,double *Jr,double *Jz){
    int bx,tx,by,ty;
    if(Lz<Lzgrid[0] || Lz>Lzgrid[NR-1]) return 0;
    topbottom<double>(Lzgrid,Lz,&bx,&tx,"Lz Interp Grids");
    double dx=(Lz-Lzgrid[bx])/(Lzgrid[tx]-Lzgrid[bx]), dxb=1-dx;
    double Ebot=dxb*E0grid[bx]+dx*E0grid[tx];
    double e=sqrt(2*(E-Ebot))/(dxb*Vmax[bx]+dx*Vmax[tx]);
    if(e<0 || e>1) return 0;
    double iu,jubb,jubt,jutb,jutt,iv,jvbb,jvbt,jvtb,jvtt;
    topbottom<double>(Unigrid,e,&by,&ty,"E Interp Grids");
    double dy=(e-Unigrid[by])/(Unigrid[ty]-Unigrid[by]), dyb=1-dy;
    by-=1; ty-=1;
    if(by<0){
        double Iubot=dxb*(dyb*Iu0[bx]+dy*Iumin[bx][ty])
                 +dx*(dyb*Iu0[tx]+dy*Iumin[tx][ty]);
        double Iutop=dxb*(dyb*Iu0[bx]+dy*Iumax[bx][ty])
                 +dx*(dyb*Iu0[tx]+dy*Iumax[tx][ty]);
        iu=(Iu-Iubot)/(Iutop-Iubot);//renormalised
        if(iu<0 || iu>1) return 0;
        jubt=linterp<double>(Iugrid[bx][0],Jr_LzEIu[bx][0],iu);
        jutt=linterp<double>(Iugrid[tx][0],Jr_LzEIu[tx][0],iu);
        *Jr=dxb*dy*jubt+dx*dy*jutt;
        double Ivbot=dxb*(dyb*Iv0[bx]+dy*Ivmin[bx][ty])
                 +dx*(dyb*Iv0[tx]+dy*Ivmin[tx][ty]);
        double Ivtop=dxb*(dyb*Iv0[bx]+dy*Ivmax[bx][ty])
                 +dx*(dyb*Iv0[tx]+dy*Ivmax[tx][ty]);
        iv=(Iv-Ivbot)/(Ivtop-Ivbot);//renormalised
        if(iv<0 || iv>1) return 0;
        jvbt=linterp<double>(Ivgrid[bx][0],Jz_LzEIv[bx][0],iv);
        jvtt=linterp<double>(Ivgrid[tx][0],Jz_LzEIv[tx][0],iv);
        *Jz=dxb*dy*jvbt+dx*dy*jvtt;
    }else{
        double Iubot=dxb*(dyb*Iumin[bx][by]+dy*Iumin[bx][ty])
                 +dx*(dyb*Iumin[tx][by]+dy*Iumin[tx][ty]);
        double Iutop=dxb*(dyb*Iumax[bx][by]+dy*Iumax[bx][ty])
                 +dx*(dyb*Iumax[tx][by]+dy*Iumax[tx][ty]);
        iu=(Iu-Iubot)/(Iutop-Iubot);//renormalised
        if(iu<0 || iu>1) return 0;
        jubb=linterp<double>(Iugrid[bx][by],Jr_LzEIu[bx][by],iu);
        jubt=linterp<double>(Iugrid[bx][ty],Jr_LzEIu[bx][ty],iu);
        jutb=linterp<double>(Iugrid[tx][by],Jr_LzEIu[tx][by],iu);
        jutt=linterp<double>(Iugrid[tx][ty],Jr_LzEIu[tx][ty],iu);
        *Jr=dxb*(dyb*jubb+dy*jubt)+dx*(dyb*jutb+dy*jutt);//combine estimates    of 4 cols
        double Ivbot=dxb*(dyb*Ivmin[bx][by]+dy*Ivmin[bx][ty])
                 +dx*(dyb*Ivmin[tx][by]+dy*Ivmin[tx][ty]);
        double Ivtop=dxb*(dyb*Ivmax[bx][by]+dy*Ivmax[bx][ty])
                 +dx*(dyb*Ivmax[tx][by]+dy*Ivmax[tx][ty]);
        iv=(Iv-Ivbot)/(Ivtop-Ivbot);//renormalised
        if(iv<0 || iv>1) return 0;
        jvbb=linterp<double>(Ivgrid[bx][by],Jz_LzEIv[bx][by],iv);
        jvbt=linterp<double>(Ivgrid[bx][ty],Jz_LzEIv[bx][ty],iv);
        jvtb=linterp<double>(Ivgrid[tx][by],Jz_LzEIv[tx][by],iv);
        jvtt=linterp<double>(Ivgrid[tx][ty],Jz_LzEIv[tx][ty],iv);
        *Jz=dxb*(dyb*jvbb+dy*jvbt)+dx*(dyb*jvtb+dy*jvtt);//combine estimates    of 4 cols
    }
    return 1;
}

VecDoub Actions_AxisymmetricFudge_InterpTables::actions(const VecDoub &XX, void* with_f){
    // THIS IS A HORRIBLE FUDGE TO STOP J_R=0
    VecDoub X = XX;
    if(fabs(X[3])<0.001)X[3]=0.001;

    bool wf;
    if(with_f)
        wf=(bool *)with_f;

    double E = Pot->H(X),Lz = Pot->Lz(X),JR,Jz;
    double Delta = UV->findDelta_interp(E,fabs(Lz));
    VecDoub Iuv = IuIv(X,Delta,E,Lz*Lz);
    if(int_acts(Lz,E,Iuv[0],Iuv[1],&JR,&Jz)==0 or no_table){
        double alpha = -1.-Delta*Delta;
        VecDoub JJ = UV->actions(X,&alpha);
        JR=JJ[0]; Jz=JJ[2];
    }
    VecDoub JJ = {JR,Lz,Jz};
    if(wf){
        if(Lz<Lcgrid[0]) JJ.push_back(Pot->R_L(Lz,Rgrid[0]));
        else if(Lz>Lcgrid[NR-1]) JJ.push_back(Pot->R_L(Lz,Rgrid[NR-1]));
        else{
            int bot,top;
            topbottom(Lcgrid, Lz, &bot, &top,"actions");
            JJ.push_back(Rgrid[bot]+(Lz-Lcgrid[bot])/(Lcgrid[top]-Lcgrid[bot])*(Rgrid[top]-Rgrid[bot]));
        }
        VecDoub freqs = PotentialFreq(JJ[3]);
        for(auto i:freqs)JJ.push_back(i);
    }
    return JJ;
}

VecDoub Actions_AxisymmetricFudge_InterpTables::PotentialFreq(double R){
    if(R<Rmin)
        return {kappagrid[0],nugrid[0],Omegacgrid[0]};
    else if(R>Rmax)
        return {kappagrid[NR-1],nugrid[NR-1],Omegacgrid[NR-1]};
    else{
        int bot,top;
        topbottom(Rgrid, R, &bot, &top,"PotentialFreqs");
        VecDoub ff(3,0);
        ff[0]=kappagrid[bot]+(R-Rgrid[bot])/(Rgrid[top]-Rgrid[bot])*(kappagrid[top]-kappagrid[bot]);
        ff[1]=nugrid[bot]+(R-Rgrid[bot])/(Rgrid[top]-Rgrid[bot])*(nugrid[top]-nugrid[bot]);
        ff[2]=Omegacgrid[bot]+(R-Rgrid[bot])/(Rgrid[top]-Rgrid[bot])*(Omegacgrid[top]-Omegacgrid[bot]);
        return ff;
    }
}

double Actions_AxisymmetricFudge_InterpTables::R_circ(double L){
    if(L<Lcgrid[0])      return Pot->R_L(L,Rgrid[0]);
    else if(L>Lcgrid[NR-1]) return Pot->R_L(L,Rgrid[NR-1]);
    else{
        int bot,top;
        topbottom(Lcgrid, L, &bot, &top,"R_circ");
        return Rgrid[bot]+(L-Lcgrid[bot])/(Lcgrid[top]-Lcgrid[bot])*(Rgrid[top]-Rgrid[bot]);
    }
}


double Actions_AxisymmetricFudge_InterpTables::L_circ(double R){
    if(R<Rmin) return Lcgrid[0];
    else if(R>Rmax) return Lcgrid[NR-1];
    else{
        int bot,top;
        topbottom(Rgrid, R, &bot, &top,"L_circ");
        return Lcgrid[bot]+(R-Rgrid[bot])/(Rgrid[top]-Rgrid[bot])*(Lcgrid[top]-Lcgrid[bot]);
    }
}


double Actions_AxisymmetricFudge_InterpTables::L_circ_current(double R){
    if(R<Rmin) return Pot->L_circ(R);
    else if(R>Rmax) return Pot->L_circ(R);
    else{
        int bot,top;
        topbottom(Rgrid, R, &bot, &top,"L_circ");
        return Lzgrid[bot]+(R-Rgrid[bot])/(Rgrid[top]-Rgrid[bot])*(Lzgrid[top]-Lzgrid[bot]);
    }
}

void Actions_AxisymmetricFudge_InterpTables::testActions(void){
    VecDoub X ={8.,0.,0.,0.2,2.1,0.3};
    bool withfreq = true;
    for(auto i:actions(X,&withfreq))std::cout<<i<<" ";
    for(auto i:UV->actions(X))std::cout<<i<<" ";
    std::cout<<std::endl;
}
//============================================================================
