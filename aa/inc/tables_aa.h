#ifndef TABLES_H
#define TABLES_H
#include "uv_orb.h"
#include "potential.h"
#include "coordsys.h"
#include "utils.h"

//============================================================================
/// Interpolation tables for Axisymmetric fudge
//============================================================================

class Actions_AxisymmetricFudge_InterpTables: public Action_Finder{
private:
    const double us=1.58, TINY_UV = SMALL;
    Potential_JS *Pot;
    std::string name;
    double Rmin,Rmax;
    int NR,NGRID;
    uv_orb *UV;
    bool no_table;
    // Grids
    VecDoub Rgrid,Lzgrid,E0grid,Vmax,Unigrid,Iu0,Iv0;
    VecDoub kappagrid,nugrid,Omegacgrid,Lcgrid;
    std::vector<VecDoub> Egrid,Iumin,Iumax,Ivmin,Ivmax;
    std::vector<std::vector<VecDoub>> Iugrid,Ivgrid,Jr_LzEIu,Jz_LzEIv;

    double Phiuv(const VecDoub& uv,double Delta);
    double dU(double u, double v, double sv2, double Phi1, double sh1sq,double Delta);
    double dV(double v,double u,double shu,double Delta2);
    double Getu0(double utry,double E,double Delta,double Iu,double Lzsq,double sh1sq, double v, double sv2);
    VecDoub IuIv(VecDoub X,double Delta,double E,double Lzsq);
    double find_umid(VecDoub X,double Delta,double E,double Lzsq);
public:

    Actions_AxisymmetricFudge_InterpTables(Potential_JS *Pot, std::string name="", bool tab=false,double Rm=1., double Rn=20.,int NR=60, int NGRID = 120,int NED=20,int NEL=5);
    ~Actions_AxisymmetricFudge_InterpTables(){}
    void use_table(){no_table=false;}
    void dont_use_table(){no_table=true;}
    double dPhiu(const VecDoub& uv,double *f,double Delta);
    void tab_acts(bool reset_freq=true);
    void get_acts(void);
    int int_acts(double Lz,double E,double Iu,double Iv,double *Jr,double *Jz);
    VecDoub actions(const VecDoub& X,void *params=nullptr);
    void testActions(void);
    VecDoub PotentialFreq(double);
    void fillLzEVmaxgrids(bool reset_freq=true);
    inline void reset(Potential_JS *pot){Pot = pot;tab_acts();}
    inline void partial_reset(Potential_JS *pot){Pot = pot;tab_acts(false);UV->reset(pot);}
    double L_circ(double R);
    double R_circ(double L);
    double L_circ_current(double R);
};

struct uturn_st{
    Actions_AxisymmetricFudge_InterpTables *AAFT;
    double E,Delta,Delta2,Iu,Lzsq,sh1sq,v,sv2;
    uturn_st(Actions_AxisymmetricFudge_InterpTables *AAFT,double E,double Delta,double Iu,double Lzsq,double sh1sq,double v,double sv2)
        :AAFT(AAFT),E(E),Delta(Delta),Delta2(Delta*Delta)
        ,Iu(Iu),Lzsq(Lzsq),sh1sq(sh1sq),v(v),sv2(sv2){}
};

// James' struct for coordinate transformation -- could be moved
class CartesianToUVCoords{
public:
    double u,v,Delta2;
    double shu2,sv2,pu,pv;
    CartesianToUVCoords(VecDoub X,double Delta){
        Delta2 = Delta*Delta;
        UVOblateSpheroidCoordSys UVV(Delta);
        VecDoub cc = UVV.xv2uv(X);
        u=cc[0];v=cc[2];
        pu=cc[3];pv=cc[5];
        shu2=cc[6];sv2=cc[7];
    }
};
//============================================================================

#endif
