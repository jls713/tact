// ============================================================================
/// \file inc/get_closed_Rz.h
// ============================================================================
/// \author Jason Sanders, based on James Binney (2014)
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
/// \brief Finding best coordinate system from shell orbits
///
/// Routines for determining Delta = gamma-alpha = a^2-c^2 by fitting ellipses
/// in the meridional plane (R,z) to shell orbits

//============================================================================

#ifndef GET_CLOSED_H
#define GET_CLOSED_H

/*! Class to find delta from a (R,z) series of a shell orbit */
class delta_st{
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
        double find_focus(void);
};


/*! Structure to store parameters for root-finding routines */
struct closed_approach_st{
    delta_st *DS;
    double R, z;
    closed_approach_st(delta_st *DS,double R, double z):DS(DS),R(R),z(z){};
};

/*! Given an energy and angular momentum will find the shell orbit
at this E,Lz and fit an ellipse to it producing an estimate of Delta(E,Lz)
*/
class find_best_delta{
    private:
        const int nmaxRi = 100;      /*! Number of stored integration points */
        Potential_JS *Pot;                      /*! Potential (axisymmetric) */
        double E, Lzsq;        /*! Energy and z-component of ang mom squared */
        double dPhi_eff(double x);
        double Ez(double *x2);
        int go_up(double x,double *Ri,double *zi,int nmaxR);
        double sorted(double x0);
    public:
        /*! find_best_delta constructor.
            \param Pot Potential (axisymmetric)
            \param E energy
            \param Lz z-component of ang. mom.
        */
        find_best_delta(Potential_JS *Pot, double E, double Lz)
            :Pot(Pot),E(E),Lzsq(Lz*Lz){};
        /*! find delta
            \param x0 -- guess at radius of shell orbit

            \return Delta
        */
        VecDoub delta(double x0, int iter=0);

        // Functions that must be accessed by root-finding but not public
        double vsqx(double x);
        inline Potential_JS* pot(void){return Pot;}
        inline double Lz2(void){return Lzsq;}
};

#endif
//============================================================================
