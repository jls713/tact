// ============================================================================
/// \file inc/tables_aa.h
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
/*! \brief Interpolation grids for axisymmetric Staeckel fudge
 *
 *  Actions_AxisymmetricFudge_InterpTables: Constructs grids for interpolating
 *  the results of the axisymmetric Staeckel fudge
 */
//============================================================================

#ifndef TABLES_H
#define TABLES_H
#include "uv_orb.h"
#include "potential.h"
#include "coordsys.h"
#include "utils.h"

//============================================================================
/*! Interpolation tables for Axisymmetric fudge
    This class is thread-safe as each call will call the uv_orb class which
    in turn creates new Actions_AxisymmetricStaeckel_Fudge instances if
    necessary.
*/
class Actions_AxisymmetricFudge_InterpTables: public Action_Finder{
private:
    const double us=1.58, TINY_UV = SMALL; /*!< Parameters for root-finding  */
    Potential_JS *Pot;                     /*!< Target potential             */
    std::string name;                      /*!< File for storing interp grids*/
    double Rmin,Rmax;                      /*!< Min and max grid radii       */
    int NR;                                /*!< Number of radial grid points */
    int NGRID;                             /*!< Number of energy grid points */
    uv_orb *UV;                            /*!< uv_orb finds the actions     */
    bool no_table;                         /*!< Option to not interpolate    */
    // Grids
    VecDoub Rgrid,Lzgrid;    /*!< Radius and ang. mom. of circular orbit grid*/
    VecDoub E0grid;                             /*!< Energy of circular orbit*/
    VecDoub Vmax;                                 /*!< Max velocity grid     */
    VecDoub Unigrid;                         /*!< Grid from 0 to 1 for interp*/
    VecDoub Iu0,Iv0;                          /*!< Iu,Iv circular orbit grid */
    VecDoub kappagrid,nugrid,Omegacgrid,Lcgrid;    /*!< Epicyclic frequencies*/
    /// There is the option to leave the epicyclic frequency grids and the
    /// corresponding angular momentum grid (Lcgrid) fixed when changing the
    /// potential -- this is useful for adiabatic relaxation of DFs
    std::vector<VecDoub> Egrid;                     /*!< Energy grid         */
    std::vector<VecDoub> Iumin,Iumax,Ivmin,Ivmax;   /*!< Min, max I3 grids   */
    std::vector<std::vector<VecDoub>> Iugrid,Ivgrid;/*!< I3 grids            */
    std::vector<std::vector<VecDoub>>Jr_LzEIu,Jz_LzEIv;/*!< action grids     */

    double Phiuv(const VecDoub& uv,double Delta);
    double dU(double u, double v, double sv2, double Phi1, double sh1sq,double Delta);
    double dV(double v,double u,double shu,double Delta2);
    double Getu0(double utry,double E,double Delta,double Iu,double Lzsq,double sh1sq, double v, double sv2);
    VecDoub IuIv(VecDoub X,double Delta,double E,double Lzsq);
    double find_umid(VecDoub X,double Delta,double E,double Lzsq);
    void tab_acts(bool reset_freq=true);
    void get_acts(void);
    int int_acts(double Lz,double E,double Iu,double Iv,double *Jr,double *Jz);
    void fillLzEVmaxgrids(bool reset_freq=true);
public:
    //! Actions_AxisymmetricFudge_InterpTables constructor
    /*!
      \param pot Potential_JS (axisymmetric)
      \param name filename for storing or reading interpolation grids
      \param tab  if true, tabulate the actions. if false, read from file
      \param Rm  minimum radial grid point
      \param Rn  maximum radial grid point
      \param NR  number of ang. mom./radial grid points
      \param NGRID  number of energy grid points
      \param NED  number of energy grid points for Delta estimation
      \param NEL  number of ang. mom. grid points for Delta estimation
    */
    Actions_AxisymmetricFudge_InterpTables(Potential_JS *Pot, std::string name="", bool tab=false,double Rm=1., double Rn=20.,int NR=60, int NGRID = 120,int NED=20,int NEL=5);
    //! Actions_AxisymmetricFudge_InterpTables destructor
    ~Actions_AxisymmetricFudge_InterpTables(){}
     //! turn on use of interpolation table
    void use_table(){no_table=false;}
     //! turn off use of interpolation table
    void dont_use_table(){no_table=true;}

    VecDoub PotentialFreq(double);
    //! dPhiu
    /*!
      Computes gradient of effective u potential
    */
    double dPhiu(const VecDoub& uv,double *f,double Delta);
    //! reset
    /*!
      Retabulates actions as a function of E, Lz and I3
    */
    inline void reset(Potential_JS *pot){Pot = pot;tab_acts();}
    //! partial_reset
    /*!
      Retabulates actions as a function of E, Lz and I3 but doesn't update
      the epicyclic frequencies grid.
    */
    inline void partial_reset(Potential_JS *pot){Pot = pot;tab_acts(false);UV->reset(pot);}
    //! L_circ
    /*!
      Ang. mom. of circular orbit from epicyclic frequency grids (so could be using old potential)
      \param R radius of orbit

      \return ang. mom. of circular orbit
    */
    double L_circ(double R);
    //! R_circ
    /*!
      Radius of circular orbit from epicyclic frequency grids (so could be using old potential)
      \param L ang. mom. of orbit

      \return radius of circular orbit
    */
    double R_circ(double L);
    //! L_circ_current
    /*!
      Ang. mom. of circular orbit in current potential
      \param R radius of orbit

      \return ang. mom. of circular orbit
    */
    double L_circ_current(double R);
    //! Finds actions
    /*!
      \param x phase-space point (x,v)
      \param params -- can pass alpha value to override Delta estimation

      \return actions -- 3D vector J=(J_R,J_phi,J_z)
    */
    VecDoub actions(const VecDoub& X,void *params=nullptr);
    // angles not implemented as they are not interpolated -- could add functionality to call uv_orb ?

    //! test Actions -- test routine for class //
    void testActions(void);
};

/*! Helper structure for finding turning points in InterpTables */
struct uturn_st{
    Actions_AxisymmetricFudge_InterpTables *AAFT;
    double E,Delta,Delta2,Iu,Lzsq,sh1sq,v,sv2;
    uturn_st(Actions_AxisymmetricFudge_InterpTables *AAFT,double E,double Delta,double Iu,double Lzsq,double sh1sq,double v,double sv2)
        :AAFT(AAFT),E(E),Delta(Delta),Delta2(Delta*Delta)
        ,Iu(Iu),Lzsq(Lzsq),sh1sq(sh1sq),v(v),sv2(sv2){}
};

/*! James' struct for coordinate transformation -- could be moved */
class CartesianToUVCoords{
public:
    double u,v,Delta2;
    double shu2,sv2,pu,pv;
    CartesianToUVCoords(VecDoub X,double Delta){
        Delta2 = Delta*Delta;
        UVProlateSpheroidCoordSys UVV(Delta);
        VecDoub cc = UVV.xv2uv(X);
        u=cc[0];v=cc[2];
        pu=cc[3];pv=cc[5];
        shu2=cc[6];sv2=cc[7];
    }
};
//============================================================================

#endif
