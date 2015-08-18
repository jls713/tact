#ifndef GET_CLOSED_H
#define GET_CLOSED_H

class delta_st{
    // Class to find delta from a (R,z) series of a shell orbit
    public:
        int N;
        double *Rp, *zp;
        double R0,D2,sqDR,sv,cv,v0;
        delta_st(double R0,double *Ri, double *zi,int N);
        ~delta_st(){delete[]Rp;delete[]zp;}
        double s2(double v, double R, double z);
        void get_v0(double Ri,double zi);
        double dv0dD2(double R,double z);
        double ds2dD2(double R,double z);
        double get_delta2(void);
};


struct closed_approach_st{
    // Structure to store parameters for root-finding routines
    delta_st *DS;
    double R, z;
    closed_approach_st(delta_st *DS,double R, double z):DS(DS),R(R),z(z){};
};

class find_best_delta{
    // Given an energy and angular momentum
    // will find the shell orbit at this E,Lz and fit an ellipse to it
    // producing an estimate of Delta(E,Lz)
    private:
        const int nmaxRi = 100;
    public:
        Potential_JS *Pot;
        double E, Lzsq;
        find_best_delta(Potential_JS *Pot, double E, double Lz)
            :Pot(Pot),E(E),Lzsq(Lz*Lz){};
        double delta(double x0);
        double dPhi_eff(double x);
        double vsqx(double x);
        double Ez(double *x2);
        int go_up(double x,double *Ri,double *zi,int nmaxR);
        double sorted(double x0);
};

#endif
