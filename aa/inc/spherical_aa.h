#ifndef SPHERICAL_AA
#define SPHERICAL_AA

#include "potential.h"
#include "aa.h"

class Actions_Spherical : public Action_Finder{
    private:
        VecDoub find_limits(double r, double E, double L);
        SphericalPotential *Pot;
        double dr0dH(double r, double L);
        double dr0dL(double r, double L);
    public:
        Actions_Spherical(SphericalPotential *Pot): Pot(Pot){};
        inline void reset_sph(SphericalPotential *pot){Pot=pot;}
        VecDoub actions(const VecDoub &x, void *params=nullptr);
        VecDoub angles_and_freqs(const VecDoub &x);
        VecDoub Hessian(const VecDoub &x);
};

struct Actions_Spherical_limits_struct{
    SphericalPotential *Pot;
    double E, L;
    Actions_Spherical_limits_struct(SphericalPotential *Pot, double E, double L)
    : Pot(Pot), E(E), L(L){}
};

struct Actions_Spherical_data_struct{
    SphericalPotential *Pot;
    double E, L, Lz, Delta, taubar, dDelta, dtaubar;
    Actions_Spherical_data_struct(SphericalPotential *Pot, double E, double L, double Lz, double Delta, double taubar, double dDelta=0., double dtaubar=0.)
    : Pot(Pot), E(E), L(L), Lz(Lz), Delta(Delta), taubar(taubar), dDelta(dDelta), dtaubar(dtaubar){}
};


#endif
