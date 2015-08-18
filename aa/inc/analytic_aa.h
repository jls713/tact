#ifndef ANALYTIC_AA_H
#define ANALYTIC_AA_H

#include "potential.h"
#include "aa.h"

class Actions_HarmonicOscillator : public Action_Finder, HarmonicOscillator{
    private:
    public:
        Actions_HarmonicOscillator(VecDoub Om): HarmonicOscillator(Om){};
        VecDoub actions(const VecDoub &x, void *params=nullptr);
        VecDoub angles(const VecDoub &x, void *params=nullptr);
};


class Actions_Isochrone : public Action_Finder, Isochrone{
    private:
    public:
        Actions_Isochrone(double GM, double b): Isochrone(GM,b){};
        VecDoub actions(const VecDoub &x, void *params=nullptr);
        VecDoub angles(const VecDoub &x, void *params=nullptr);
        VecDoub freq(const VecDoub &x);
        VecDoub Hessian(const VecDoub &x);
};

#endif
