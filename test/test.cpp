// ============================================================================
/// \file test/test.cpp
// ============================================================================
/// \author Jason Sanders
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
/// \brief test functions
///
//============================================================================

#include "../general/utils.h"
#include "../pot/inc/potential.h"
#include "../pot/inc/Multipole.h"
#include "../pot/inc/orbit.h"
#include "../aa/inc/spherical_aa.h"
#include "../aa/inc/analytic_aa.h"
#include "../aa/inc/adiabatic_aa.h"
#include "../aa/inc/stackel_aa.h"
#include "../aa/inc/stackel_fit.h"
#include "../aa/inc/genfunc_aa.h"
#include "../aa/inc/uv_orb.h"
#ifdef TORUS
#include "../aa/inc/it_torus.h"
#endif
#include "gtest/gtest.h"

namespace {

TEST(MultipoleT,Sph){
    for(int p=0;p<3;++p){
    TestDensity_Hernquist rho(1.,1.,{1.,1.,1.});
    MultipoleExpansion_Spherical ME(&rho,100,1.,0.01,200.);
    VecDoub X = {1e-5,1e-5,1e-5};

    int NMAX = 200;
    VecDoub xx(NMAX,0), exact(NMAX,0), triaxial(NMAX,0), multipole(NMAX,0);

    // #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
        double x = 0.0001*(double)xn+.0001;
        xx[xn] = x;
        VecDoub X2 = X;
        X2[p]=x;
        exact[xn] = rho.Phi(X2);
        multipole[xn] = ME.Phi(X2);
    }
    for(int xn = 0; xn<NMAX; xn++)
      EXPECT_NEAR(exact[xn],multipole[xn],1e-3*fabs(exact[xn]));
  }
}
TEST(MultipoleT,Axisym){
    Miyamoto_NagaiDensity rho(1.,1.,0.7);
    MultipoleExpansion_Axisymmetric ME2(&rho,1000,30,20,1.,0.001,1000.);

    int NMAX = 50;
    for(int p=0;p<3;++p){
    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
        VecDoub X = {0.001,0.001,0.01};
        double centre2 = rho.Phi(X);
        double centre3 = ME2.Phi(X);
        double x = 0.1*(double)xn;
        VecDoub X2=X;
        X2[p]=x;
        double triaxial = (rho.Phi(X2)-centre2);
        double multipole = (ME2.Phi(X2)-centre3);
        EXPECT_NEAR(triaxial,multipole,5e-3*fabs(triaxial));
        triaxial = rho.Forces(X2)[p];
        multipole = ME2.Forces(X2)[p];
        if(triaxial==0.) continue;
        EXPECT_NEAR(triaxial,multipole,1e-3*fabs(triaxial));
    }
  }
}
TEST(MultipoleT,Stackel){
    TestDensity_Stackel rho(1.,-30.,-10.);
    MultipoleExpansion ME(&rho,5000,12,12,-1,1.,0.0001,10000.,false,true,true);

    for(auto qq: {"No","general"}){
    for(int p=0;p<3;++p){
    VecDoub X = {1.,1.,1.};
    double centre  = rho.Phi(X);
    double centre3 = ME.Phi(X);

    int NMAX = 100;

    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
        double exact,multipole;
        double x = (double)xn+.1;
        VecDoub X2 = X;
        X2[p]=x;
        if(qq=="general"){
            X[0]=x/2.;X[1]=x/2.;X[2]=x/2.;
        }
        exact = (rho.Phi(X2)-centre);
        multipole = (ME.Phi(X2)-centre3);
        EXPECT_NEAR(exact,multipole,5e-3*fabs(exact));
    }
  }}
}
TEST(MultipoleT,StackelForces){
    TestDensity_Stackel rho(1.,-30.,-10.);
    MultipoleExpansion ME(&rho,5000,12,12,-1,1.,0.0001,10000.,false,true,true);

    int NMAX = 50;

    for(auto qq: {"general"}){
    for(int q=0;q<1;++q){
    for(int p=0;p<1;++p){

    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
      double exact, multipole;
        double x = (double)xn+1.;
        VecDoub X(3,1);
        X[p]=x;
        if(qq=="general"){
            X[0]=x/2.;X[1]=x/3.;X[2]=x;
        }
        VecDoub ex = rho.Forces(X);
        exact = rho.Forces(X)[q];
        multipole = ME.Forces(X)[q];
        EXPECT_NEAR(exact,multipole,5e-3*norm(ex));
    }
    }}}
}
TEST(MultipoleT,StackelTP){
    TestDensity_Stackel rho(1.,-30.,-10.);
    TriaxialPotential TP(&rho,1e8);
    MultipoleExpansion ME(&rho,5000,12,12,-1,1.,0.0001,10000.,false,true,true);
    double centre  = TP.Phi({1.,1.,1.});
    double centre3 = ME.Phi({1.,1.,1.});

    int NMAX = 50;

    for(auto qq: {"general"}){
    for(int q=0;q<1;++q){
    for(int p=0;p<1;++p){

    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
      double exact, multipole;
        double x = (double)xn+1.;
        VecDoub X(3,1);
        X[p]=x;
        if(qq=="general"){
            X[0]=x/2.;X[1]=x/3.;X[2]=x;
        }
        VecDoub ex = TP.Forces(X);
        exact = TP.Forces(X)[q];
        multipole = ME.Forces(X)[q];
        EXPECT_NEAR(exact,multipole,5e-3*norm(ex));
        exact = TP.Phi(X)-centre;
        multipole = ME.Phi(X)-centre3;
        EXPECT_NEAR(exact,multipole,5e-3*fabs(exact));
    }
    }}}
}

TEST(MultipoleT,StackelTriax){
    TestDensity_Stackel rho(1.,-30.,-10.);
    TriaxialPotential TP(&rho,1e8);
    double centre=0.,centre3=0.;
    int NMAX = 50;

    for(auto qq: {"general"}){
    for(int q=1;q<3;++q){
    for(int p=0;p<1;++p){

    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
      double exact, triaxial;
        double x = (double)xn+1.;
        VecDoub X(3,1);
        X[p]=x;
        if(qq=="general"){
            X[0]=x/2.;X[1]=x/3.;X[2]=x;
        }
        VecDoub ex = TP.Forces(X);
        exact = rho.Forces(X)[q];
        triaxial = TP.Forces(X)[q];
        EXPECT_NEAR(exact,triaxial,5e-3*norm(ex));
        exact = rho.Phi(X)-centre;
        triaxial = TP.Phi(X)-centre3;
        EXPECT_NEAR(exact,triaxial,5e-3*fabs(exact));
    }
    }}}
}

TEST(MultipoleT,Hernquist){
    TestDensity_Hernquist rho(1.,10.,{1.,0.6,0.3});
    TriaxialPotential TP(&rho,1e6);
    MultipoleExpansion ME(&rho,5000,12,12,-1,1.,0.0001,10000.,false,true,true);
    double centre  = 0.;//TP.Phi({1.,1.,1.});
    double centre3 = 0.;//ME.Phi({1.,1.,1.});

    int NMAX = 50;

    for(auto qq: {"No","general"}){
    for(int q=0;q<3;++q){
    for(int p=0;p<3;++p){

    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
      double exact, multipole;
        double x = (double)xn+1.;
        VecDoub X(3,1);
        X[p]=x;
        if(qq=="general"){
            X[0]=x/2.;X[1]=x/3.;X[2]=x;
        }
        VecDoub ex = TP.Forces(X);
        exact = TP.Forces(X)[q];
        multipole = ME.Forces(X)[q];
        EXPECT_NEAR(exact,multipole,5e-3*norm(ex));
        exact = TP.Phi(X)-centre;
        multipole = ME.Phi(X)-centre3;
        EXPECT_NEAR(exact,multipole,5e-3*fabs(exact));
    }
    }}}
}
TEST(GusTest,Number1){
  VecDoub x = {0.184035,0.184035,0.476411,-97.8458,70.3121,15.3668};
  PowerLaw pot(430107.95338751527*.5,0.5,1.,0.9);
  Actions_AxisymmetricStackel_Fudge af(&pot,5.);
  printVector(af.actions(x));
}

TEST(PotentialTest,Isochrone){
  IsochronePotential ISO(1.,1.);
  EXPECT_DOUBLE_EQ(-1./3.,ISO.Phi({1.,1.,1.}));
  EXPECT_DOUBLE_EQ(-1./3.,ISO.Phi_r(sqrt(3.)));
  EXPECT_DOUBLE_EQ(-1./18.,ISO.Forces({1.,1.,1.})[0]);
  EXPECT_DOUBLE_EQ(1./18.*sqrt(3.),ISO.dPhi_r(sqrt(3.)));
  EXPECT_DOUBLE_EQ(0.,ISO.Lz({1.,1.,1.,1.,1.,1.}));
  EXPECT_DOUBLE_EQ(1.5-1./3.,ISO.H({1.,1.,1.,1.,1.,1.}));
  EXPECT_DOUBLE_EQ(0.,ISO.L({1.,1.,1.,1.,1.,1.}));
  EXPECT_DOUBLE_EQ(0.,ISO.Lvec({1.,1.,1.,1.,1.,1.})[0]);
  EXPECT_DOUBLE_EQ(0.12132034355964258,ISO.Mass(1.));
  EXPECT_DOUBLE_EQ(-5.7734693586591491e-06,ISO.Phi_max());
  EXPECT_DOUBLE_EQ(-0.35355339059327379,ISO.E_circ(1.));
  EXPECT_DOUBLE_EQ(0.34831069974900652,ISO.L_circ(1.));
}

TEST(PotentialTest,Logarithmic){
  Logarithmic Log(1.,1.,1.);
  EXPECT_DOUBLE_EQ(.5*(log(3.)-log(3e20)),Log.Phi({1.,1.,1.}));
  EXPECT_DOUBLE_EQ(-1./3.,Log.Forces({1.,1.,1.})[0]);
}

TEST(PotentialTest,IsochroneActions){
  IsochronePotential ISO(1.,1.);
  EXPECT_DOUBLE_EQ(std::numeric_limits<double>::infinity(),ISO.JR({1.,1.,1.,1e10,1e10,1e10}));
  EXPECT_DOUBLE_EQ(std::numeric_limits<double>::infinity(),ISO.Actions({1.,1.,1.,1e10,1e10,1e10})[0]);
  EXPECT_DOUBLE_EQ(0.25326797943307078,ISO.JR({1.,1.,1.,0.1,0.1,0.1}));
  EXPECT_DOUBLE_EQ(0.25326797943307078,ISO.Actions({1.,1.,1.,0.1,0.1,0.1})[0]);
  EXPECT_DOUBLE_EQ(0.,ISO.Actions({1.,1.,1.,1.,1.,1.})[1]);
  EXPECT_DOUBLE_EQ(0.,ISO.Actions({1.,1.,1.,1.,1.,1.})[2]);
  EXPECT_DOUBLE_EQ(0.50800521286330935,ISO.Omega({1.,1.,1.,0.1,0.1,0.1})[0]);
  EXPECT_DOUBLE_EQ(0.25400260643165468,ISO.Omega({1.,1.,1.,0.1,0.1,0.1})[1]);
  EXPECT_DOUBLE_EQ(-1.2160333333333331,ISO.Hessian({1.,1.,1.,0.1,0.1,0.1})[0]);
  EXPECT_DOUBLE_EQ(-0.60801666666666654,ISO.Hessian({1.,1.,1.,0.1,0.1,0.1})[1]);
  EXPECT_DOUBLE_EQ(-0.17700703011750593,ISO.Hessian({1.,1.,1.,0.1,0.1,0.1})[2]);
}
//=============================================================================
TEST(ActionTest_Spherical,Spherical){
  IsochronePotential Pot(1.,1.);
  Actions_Spherical AA(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({radius,0.,0.})[0]);
  VecDoub X = {radius,0.,0.,0.,Vc,0.};
  Actions_Isochrone Iso(1.,1.);
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub Acts = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acts[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acts[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acts[2],0.001);
}

//=============================================================================

TEST(ActionTest_TriaxialStack,StackelFudge){
  StackelTriaxial Pot(1.,-30.,-20.);
  Actions_TriaxialStackel_Fudge AA(&Pot,-30.,-20.);
  Actions_TriaxialStackel AAS(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
  VecDoub ActsTrue = AAS.actions(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
}

TEST(ActionTest_TriaxialStack,StackelFudgeWrongAlphaBeta){
  StackelTriaxial Pot(1.,30.,20.);
  ASSERT_THROW(Actions_TriaxialStackel_Fudge AA(&Pot,30.,20.),std::invalid_argument);
}
TEST(ActionTest_Iso,StackelFudge){
  Isochrone Pot(1.,1.,0.99999);
  Actions_AxisymmetricStackel_Fudge AA(&Pot,1.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
  Actions_Isochrone Iso(1.,1.);
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub FreqsTrue = Iso.freq(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
  EXPECT_NEAR(FreqsTrue[0],Angs[3],0.001);
  EXPECT_NEAR(FreqsTrue[1],Angs[4],0.001);
  EXPECT_NEAR(FreqsTrue[2],Angs[5],0.004);
}

// NANs in Stackel forces on axis

TEST(ActionTest_Stack,StackelFudge){
  StackelOblate_PerfectEllipsoid Pot(1.,-30.);
  Actions_AxisymmetricStackel_Fudge AA(&Pot,-30.);
  Actions_AxisymmetricStackel AAS(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
  VecDoub ActsTrue = AAS.actions(X);
  VecDoub FreqsTrue = AAS.angles(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
  EXPECT_NEAR(FreqsTrue[3],Angs[3],0.001);
  EXPECT_NEAR(FreqsTrue[4],Angs[4],0.001);
  EXPECT_NEAR(FreqsTrue[5],Angs[5],0.004);
}

TEST(ActionTest_Zeros,StackelFudge){
  Logarithmic Pot(1.,1.,0.9);
  Actions_AxisymmetricStackel_Fudge AA(&Pot,1.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.,0.,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(0.,Acs[0],1e-10);
  EXPECT_DOUBLE_EQ(Vc,Acs[1]);
  EXPECT_DOUBLE_EQ(0.,Acs[2]);
  VecDoub Angs = AA.angles(X);
  EXPECT_DOUBLE_EQ(0.,Angs[0]);
  EXPECT_DOUBLE_EQ(0.,Angs[1]);
  EXPECT_DOUBLE_EQ(0.,Angs[2]);
  EXPECT_DOUBLE_EQ(0.,Angs[3]);
  EXPECT_DOUBLE_EQ(Vc,Angs[4]);
  EXPECT_DOUBLE_EQ(0.,Angs[5]);

  X = {40.,0.,40.,0.,0.,0.0};
  Acs = AA.actions(X);
  EXPECT_EQ(Acs[0]>10.,1);
  printVector(Acs);
  X = {40.,0.,0.,0.,0.,0.};
  EXPECT_EQ(AA.actions(X)[0]>10.,1);
  X = {40.,0.,0.,0.,0.,.2};
  Acs = AA.actions(X);
  EXPECT_EQ(Acs[2]>0.,1);
  X = {40.,0.,0.,0.,0.,0.};
  Acs = AA.actions(X);
  EXPECT_EQ(Acs[2]==0.,1);
}


TEST(ActionTest_Planar,StackelFudge){
  IsochronePotential Pot(1.,1.);
  Actions_AxisymmetricStackel_Fudge AA(&Pot,1.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.,0.,0.1*Vc,Vc,0.};
  VecDoub Acs = AA.actions(X);
  Actions_Isochrone Iso(1.,1.);
  VecDoub ActsTrue = Iso.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],1e-10);
  EXPECT_DOUBLE_EQ(ActsTrue[1],Acs[1]);
  EXPECT_DOUBLE_EQ(ActsTrue[2],Acs[2]);
  VecDoub Angs = AA.angles(X);
  VecDoub AngsTrue = Iso.angles(X);
  EXPECT_NEAR(AngsTrue[0],Angs[0],2e-4);
  AngsTrue = Iso.freq(X);
  EXPECT_NEAR(AngsTrue[0],Angs[3],1e-4);
  EXPECT_NEAR(AngsTrue[1],Angs[4],1e-4);
}

TEST(ActionTest_Zero,StackelFudge){
  IsochronePotential Pot(1.,1.);
  Actions_AxisymmetricStackel_Fudge AA(&Pot,1.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {0.,0.,0.,0.1*Vc,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_DOUBLE_EQ(0.,Acs[0]);
  EXPECT_DOUBLE_EQ(0.,Acs[1]);
  EXPECT_DOUBLE_EQ(0.,Acs[2]);
  VecDoub Angs = AA.angles(X);
  EXPECT_DOUBLE_EQ(0.,Acs[0]);
  EXPECT_DOUBLE_EQ(0.,Angs[1]);
  EXPECT_DOUBLE_EQ(0.,Angs[2]);
  EXPECT_DOUBLE_EQ(0.,Angs[3]);
  EXPECT_DOUBLE_EQ(0.,Angs[4]);
  EXPECT_DOUBLE_EQ(0.,Angs[5]);
}

TEST(ActionTest_Small,StackelFudge){
  Logarithmic Pot(1.,1.,0.9);
  Actions_AxisymmetricStackel_Fudge AA(&Pot,1.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  std::cout<<"small z, zero vR, zero vz"<<std::endl;
  VecDoub X = {1.,2e-5,2e-5,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  VecDoub Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"zero z, small vR, small vz"<<std::endl;
  X = {1.,0.,0.,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vR, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,.3*Vc};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
}
//=============================================================================

TEST(ActionTest_Iso,StackelCAA){
  Isochrone Pot(1.,1.,0.99999);
  Actions_CylindricalAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
  Actions_Isochrone Iso(1.,1.);
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub FreqsTrue = Iso.freq(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
  EXPECT_NEAR(FreqsTrue[0],Angs[3],0.003);
  EXPECT_NEAR(FreqsTrue[1],Angs[4],0.001);
  EXPECT_NEAR(FreqsTrue[2],Angs[5],0.008);
}

TEST(ActionTest_Stack,StackelCAA){
  StackelOblate_PerfectEllipsoid Pot(1.,-30.);
  Actions_CylindricalAdiabaticApproximation AA(&Pot,"",false,false,0.2,3.,2.);
  Actions_AxisymmetricStackel AAS(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
  VecDoub ActsTrue = AAS.actions(X);
  VecDoub FreqsTrue = AAS.angles(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
  EXPECT_NEAR(FreqsTrue[3],Angs[3],0.001);
  EXPECT_NEAR(FreqsTrue[4],Angs[4],0.001);
  EXPECT_NEAR(FreqsTrue[5],Angs[5],0.004);
}

TEST(ActionTest_Zeros,StackelCAA){
  Logarithmic Pot(1.,1.,0.9);
  Actions_CylindricalAdiabaticApproximation AA(&Pot,"",false,false,0.00001,3.,4.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.,0.,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(0.,Acs[0],1e-10);
  EXPECT_DOUBLE_EQ(Vc,Acs[1]);
  EXPECT_DOUBLE_EQ(0.,Acs[2]);
  VecDoub Angs = AA.angles(X);
  EXPECT_DOUBLE_EQ(0.,Angs[0]);
  EXPECT_DOUBLE_EQ(0.,Angs[1]);
  EXPECT_DOUBLE_EQ(0.,Angs[2]);
  EXPECT_DOUBLE_EQ(0.,Angs[3]);
  EXPECT_DOUBLE_EQ(Vc,Angs[4]);
  EXPECT_DOUBLE_EQ(0.,Angs[5]);
  X = {1.,0.,1.,0.,0.,0.};
  Acs = AA.actions(X);
  EXPECT_EQ(Acs[0]>0.,1);
  X = {1.,0.,0.,0.,0.,0.2};
  Acs = AA.actions(X);
  EXPECT_EQ(Acs[2]>0.,1);
  X = {1.,0.,0.,0.,0.,0.};
  EXPECT_EQ(AA.actions(X)[0]>0.,1);
}


TEST(ActionTest_Planar,StackelCAA){
  IsochronePotential Pot(1.,1.);
  Actions_CylindricalAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.,0.,0.1*Vc,Vc,0.};
  VecDoub Acs = AA.actions(X);
  Actions_Isochrone Iso(1.,1.);
  VecDoub ActsTrue = Iso.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],1e-10);
  EXPECT_DOUBLE_EQ(ActsTrue[1],Acs[1]);
  EXPECT_DOUBLE_EQ(ActsTrue[2],Acs[2]);
  VecDoub Angs = AA.angles(X);
  VecDoub AngsTrue = Iso.angles(X);
  EXPECT_NEAR(AngsTrue[0],Angs[0],2e-4);
  AngsTrue = Iso.freq(X);
  EXPECT_NEAR(AngsTrue[0],Angs[3],1e-4);
  EXPECT_NEAR(AngsTrue[1],Angs[4],1e-4);
}

TEST(ActionTest_Zero,StackelCAA){
  IsochronePotential Pot(1.,1.);
  Actions_CylindricalAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {0.,0.,0.,0.1*Vc,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_DOUBLE_EQ(0.,Acs[0]);
  EXPECT_DOUBLE_EQ(0.,Acs[1]);
  EXPECT_DOUBLE_EQ(0.,Acs[2]);
  VecDoub Angs = AA.angles(X);
  EXPECT_DOUBLE_EQ(0.,Acs[0]);
  EXPECT_DOUBLE_EQ(0.,Angs[1]);
  EXPECT_DOUBLE_EQ(0.,Angs[2]);
  EXPECT_DOUBLE_EQ(0.,Angs[3]);
  EXPECT_DOUBLE_EQ(0.,Angs[4]);
  EXPECT_DOUBLE_EQ(0.,Angs[5]);
}

TEST(ActionTest_Small,StackelCAA){
  IsochronePotential Pot(1.,1.);
  Actions_CylindricalAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  std::cout<<"small z, zero vR, zero vz"<<std::endl;
  VecDoub X = {1.,2e-5,2e-5,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  VecDoub Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"zero z, small vR, small vz"<<std::endl;
  X = {1.,0.,0.,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vR, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,.3*Vc};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
}


TEST(ActionTest_Small,StackelSF){
  Logarithmic Pot(1.,1.,0.9);
  Actions_AxisymmetricStackel_Fudge AA(&Pot,1.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  std::cout<<"small z, zero vR, zero vz"<<std::endl;
  VecDoub X = {1.,2e-5,2e-5,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  VecDoub Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"zero z, small vR, small vz"<<std::endl;
  X = {1.,0.,0.,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vR, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,.3*Vc};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
}


#ifdef TORUS
// Memory corruption problem in Torus?
// TEST(ActionTest_Small,ItTorus){
//   Logarithmic Pot(1.,1.,0.9);
//   Actions_AxisymmetricStackel_Fudge AAF(&Pot,100.);
//   IterativeTorusMachine AA(&AAF,&Pot,1e-8,5,1e-2);
//   double radius = 1.;
//   double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
//   std::cout<<"small z, zero vR, zero vz"<<std::endl;
//   VecDoub X = {1.,2e-5,2e-5,0.,Vc,0.};
//   VecDoub Acs = AA.actions(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
//   VecDoub Angs = AA.angles(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
//   std::cout<<"zero z, small vR, small vz"<<std::endl;
//   X = {1.,0.,0.,2e-5,Vc,2e-5};
//   Acs = AA.actions(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
//   Angs = AA.angles(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
//   std::cout<<"small z, small vR, small vz"<<std::endl;
//   X = {1.,2e-5,2e-5,2e-5,Vc,2e-5};
//   Acs = AA.actions(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
//   Angs = AA.angles(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
//   std::cout<<"small z, small vz"<<std::endl;
//   X = {1.,2e-5,2e-5,.2*Vc,Vc,2e-5};
//   Acs = AA.actions(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
//   Angs = AA.angles(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
//   std::cout<<"small z"<<std::endl;
//   X = {1.,2e-5,2e-5,.2*Vc,Vc,.3*Vc};
//   Acs = AA.actions(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
//   Angs = AA.angles(X);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
//   EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
// }
#endif

TEST(ActionTest_Small,Genfunc){
  Isochrone Pot(1.,1.,1.,0.9);
  Actions_Genfunc AA(&Pot,"axisymmetric");
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  std::cout<<"small z, zero vR, zero vz"<<std::endl;
  VecDoub X = {1.,2e-5,2e-5,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  VecDoub Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"zero z, small vR, small vz"<<std::endl;
  X = {1.,0.,0.,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vR, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,.3*Vc};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);

  X = {40.,0.,40.,0.,0.,0.};
  EXPECT_EQ(AA.actions(X)[0]>0.,1);
  Actions_Isochrone Iso(1.,1.);
  radius = 1.;
  Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  X = {1.,0.,0.,0.,0.,0.};
  VecDoub ActsTrue = Iso.actions(X);
  X = {1.,0.,0.,0.,0.,0.};
  EXPECT_EQ(AA.actions(X)[0]>0.,1);
}
//=============================================================================

TEST(ActionTest,StackelSAA){
  Logarithmic Pot(1.,1.,0.8);
  Actions_SpheroidalAdiabaticApproximation AA(&Pot,"",false,false,1.,0.02,3.,1.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.4*Vc,.6*Vc,0.3*Vc};
  Actions_Genfunc Iso(&Pot,"axisymmetric");
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub AngsTrue = Iso.angles(X);
  double guess_alpha=1.;
  VecDoub Angs = AA.angles(X,&guess_alpha);
  VecDoub Acs = AA.actions(X,&guess_alpha);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.006);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.0075);
  EXPECT_NEAR(AngsTrue[3],Angs[3],0.02);
  EXPECT_NEAR(AngsTrue[4],Angs[4],0.04);
  EXPECT_NEAR(AngsTrue[5],Angs[5],0.065);
}
TEST(ActionTest,StackelCAA){
  Logarithmic Pot(1.,1.,0.8);
  Actions_CylindricalAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.4*Vc,.6*Vc,0.3*Vc};
  Actions_Genfunc Iso(&Pot,"axisymmetric");
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub AngsTrue = Iso.angles(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.02);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.0075);
  EXPECT_NEAR(AngsTrue[3],Angs[3],0.02);
  EXPECT_NEAR(AngsTrue[4],Angs[4],0.02);
  EXPECT_NEAR(AngsTrue[5],Angs[5],0.065);
}

TEST(ActionTest,StackelUV){
  Logarithmic Pot(1.,1.,0.8);
  uv_orb AA(&Pot,0.1,3.,10,10,"");
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.4*Vc,.6*Vc,0.3*Vc};
  Actions_Genfunc Iso(&Pot,"axisymmetric");
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub AngsTrue = Iso.angles(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.003);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.005);
  EXPECT_NEAR(AngsTrue[3],Angs[3],0.02);
  EXPECT_NEAR(AngsTrue[4],Angs[4],0.02);
  EXPECT_NEAR(AngsTrue[5],Angs[5],0.065);
  X = {40.,0.,40.,0.,0.,0.};
  EXPECT_EQ(AA.actions(X)[0]>10.,1);
  X = {40.,0.,0.,0.,0.,0.};
  EXPECT_EQ(AA.actions(X)[0]>10.,1);
}

TEST(DeltaTest,StackelUV){
  VecDoub X = {0.15,0.01,0.01,0.02,0.02,0.2};
  TestDensity_Plummer Plum(1.,1.,{1.,1.,0.9});
  MultipoleExpansion_Axisymmetric Pot(&Plum,150,12,-1,1.,0.001,1000.);
  uv_orb AA(&Pot,0.01,100.,10,10,"");
  printVector(AA.actions(X));
}

TEST(ActionTest,StackelF){
  Logarithmic Pot(1.,1.,0.8);
  Actions_AxisymmetricStackel_Fudge AA(&Pot,100.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.4*Vc,.6*Vc,0.3*Vc};
  Actions_Genfunc Iso(&Pot,"axisymmetric");
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub AngsTrue = Iso.angles(X);
  double guess_alpha=1.;
  VecDoub Angs = AA.angles(X,&guess_alpha);
  VecDoub Acs = AA.actions(X,&guess_alpha);
  printVector(Acs);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.003);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.005);
  EXPECT_NEAR(AngsTrue[3],Angs[3],0.02);
  EXPECT_NEAR(AngsTrue[4],Angs[4],0.02);
  EXPECT_NEAR(AngsTrue[5],Angs[5],0.131);
}
TEST(ActionTest,StackelFit){
  Logarithmic Pot(1.,1.,0.8);
  Actions_StackelFit AA(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.4*Vc,.6*Vc,0.3*Vc};
  Actions_Genfunc Iso(&Pot,"axisymmetric");
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub AngsTrue = Iso.angles(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.003);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.003);
  EXPECT_NEAR(AngsTrue[3],Angs[3],0.02);
  EXPECT_NEAR(AngsTrue[4],Angs[4],0.02);
  EXPECT_NEAR(AngsTrue[5],Angs[5],0.065);
}
#ifdef TORUS
TEST(ActionTest,ItTorus){
  Logarithmic Pot(1.,1.,0.8);
  Actions_AxisymmetricStackel_Fudge AAF(&Pot,100.);
  IterativeTorusMachine AA(&AAF,&Pot,1e-8,5,1e-3);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.4*Vc,.6*Vc,0.3*Vc};
  Actions_Genfunc Iso(&Pot,"axisymmetric");
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub AngsTrue = Iso.angles(X);
  double guess_alpha=1.;
  VecDoub Angs = AA.angles(X,&guess_alpha);
  VecDoub Acs = AA.actions(X,&guess_alpha);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.002);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.0075);
  EXPECT_NEAR(AngsTrue[3],Angs[3],0.02);
  EXPECT_NEAR(AngsTrue[4],Angs[4],0.02);
  EXPECT_NEAR(AngsTrue[5],Angs[5],0.065);
}
#endif
TEST(ActionTest,AvGenfunc){
  Logarithmic Pot(1.,1.,0.8);
  Actions_Genfunc_Average AA(&Pot,"axisymmetric");
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.4*Vc,.6*Vc,0.3*Vc};
  Actions_Genfunc Iso(&Pot,"axisymmetric");
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub AngsTrue = Iso.angles(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.002);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.0075);
  EXPECT_NEAR(AngsTrue[3],Angs[3],0.02);
  EXPECT_NEAR(AngsTrue[4],Angs[4],0.02);
  EXPECT_NEAR(AngsTrue[5],Angs[5],0.065);
}

//=============================================================================
TEST(ActionTest_Iso,StackelFit){
  Isochrone Pot(1.,1.,0.99999);
  Actions_StackelFit AA(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
  Actions_Isochrone Iso(1.,1.);
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub FreqsTrue = Iso.freq(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
  EXPECT_NEAR(FreqsTrue[0],Angs[3],0.001);
  EXPECT_NEAR(FreqsTrue[1],Angs[4],0.001);
  EXPECT_NEAR(FreqsTrue[2],Angs[5],0.004);
}

TEST(ActionTest_Stack,StackelFit){
  StackelOblate_PerfectEllipsoid Pot(100.,-30.);
  Actions_StackelFit AA(&Pot);
  Actions_AxisymmetricStackel AAS(&Pot);
  // Actions_Genfunc AAS(&Pot,"axisymmetric");
  double radius = 10.;
  double Vc = sqrt(radius*-Pot.Forces({radius,0.1,0.1})[0]);
  VecDoub X = {radius,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
  VecDoub ActsTrue = AAS.actions(X);
  VecDoub FreqsTrue = AAS.angles(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
  EXPECT_NEAR(FreqsTrue[3],Angs[3],0.004);
  EXPECT_NEAR(FreqsTrue[4],Angs[4],0.004);
  EXPECT_NEAR(FreqsTrue[5],Angs[5],0.004);
}


TEST(ActionTest_Zeros,StackelFit){
  Logarithmic Pot(1.,1.,0.9);
  Actions_StackelFit AA(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.,0.,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(0.,Acs[0],1e-10);
  EXPECT_DOUBLE_EQ(Vc,Acs[1]);
  EXPECT_DOUBLE_EQ(0.,Acs[2]);
  VecDoub Angs = AA.angles(X);
  EXPECT_DOUBLE_EQ(0.,Angs[0]);
  EXPECT_DOUBLE_EQ(0.,Angs[1]);
  EXPECT_DOUBLE_EQ(0.,Angs[2]);
  EXPECT_DOUBLE_EQ(0.,Angs[3]);
  EXPECT_DOUBLE_EQ(Vc,Angs[4]);
  EXPECT_DOUBLE_EQ(0.,Angs[5]);
  X = {40.,0.,40.,0.,0.,0.};
  Acs = AA.actions(X);
  EXPECT_EQ(Acs[0]>.1,1);
  X = {40.,0.,0.,0.,0.,0.2};
  Acs = AA.actions(X);
  EXPECT_EQ(Acs[2]>0.,1);
  X = {40.,0.,0.,0.,0.,0.};
  EXPECT_EQ(AA.actions(X)[0]>10.,1);
}


TEST(ActionTest_Planar,StackelFit){
  IsochronePotential Pot(1.,1.);
  Actions_StackelFit AA(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.,0.,0.1*Vc,Vc,0.};
  VecDoub Acs = AA.actions(X);
  Actions_Isochrone Iso(1.,1.);
  VecDoub ActsTrue = Iso.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],1e-10);
  EXPECT_DOUBLE_EQ(ActsTrue[1],Acs[1]);
  EXPECT_DOUBLE_EQ(ActsTrue[2],Acs[2]);
  VecDoub Angs = AA.angles(X);
  VecDoub AngsTrue = Iso.angles(X);
  EXPECT_NEAR(AngsTrue[0],Angs[0],2e-4);
  AngsTrue = Iso.freq(X);
  EXPECT_NEAR(AngsTrue[0],Angs[3],1e-4);
  EXPECT_NEAR(AngsTrue[1],Angs[4],1e-4);
}

TEST(ActionTest_Zero,StackelFit){
  IsochronePotential Pot(1.,1.);
  Actions_StackelFit AA(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {0.,0.,0.,0.1*Vc,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_DOUBLE_EQ(0.,Acs[0]);
  EXPECT_DOUBLE_EQ(0.,Acs[1]);
  EXPECT_DOUBLE_EQ(0.,Acs[2]);
  VecDoub Angs = AA.angles(X);
  EXPECT_DOUBLE_EQ(0.,Acs[0]);
  EXPECT_DOUBLE_EQ(0.,Angs[1]);
  EXPECT_DOUBLE_EQ(0.,Angs[2]);
  EXPECT_DOUBLE_EQ(0.,Angs[3]);
  EXPECT_DOUBLE_EQ(0.,Angs[4]);
  EXPECT_DOUBLE_EQ(0.,Angs[5]);
}

TEST(ActionTest_Small,StackelFit){
  Logarithmic Pot(1.,1.,0.99);
  Actions_StackelFit AA(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  std::cout<<"small z, zero vR, zero vz"<<std::endl;
  VecDoub X = {1.,2e-5,2e-5,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  VecDoub Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"zero z, small vR, small vz"<<std::endl;
  X = {1.,0.,0.,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vR, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  // std::cout<<"small z"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,.3*Vc};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
}
//=============================================================================

TEST(ActionTest_Iso,StackelSAA){
  Isochrone Pot(1.,1.,0.99999);
  Actions_CylindricalAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
  Actions_Isochrone Iso(1.,1.);
  VecDoub ActsTrue = Iso.actions(X);
  VecDoub FreqsTrue = Iso.freq(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
  EXPECT_NEAR(FreqsTrue[0],Angs[3],0.003);
  EXPECT_NEAR(FreqsTrue[1],Angs[4],0.001);
  EXPECT_NEAR(FreqsTrue[2],Angs[5],0.008);
}

TEST(ActionTest_Stack,StackelSAA){
  StackelOblate_PerfectEllipsoid Pot(1.,-30.);
  Actions_CylindricalAdiabaticApproximation AA(&Pot,"",false,false,0.2,3.,2.);
  Actions_AxisymmetricStackel AAS(&Pot);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
  VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
  VecDoub ActsTrue = AAS.actions(X);
  VecDoub FreqsTrue = AAS.angles(X);
  VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
  EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
  EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
  EXPECT_NEAR(FreqsTrue[3],Angs[3],0.001);
  EXPECT_NEAR(FreqsTrue[4],Angs[4],0.001);
  EXPECT_NEAR(FreqsTrue[5],Angs[5],0.004);
}

TEST(ActionTest_Zeros,StackelSAA){
  Logarithmic Pot(1.,1.,0.9);
  Actions_SpheroidalAdiabaticApproximation AA(&Pot,"",false,false,100.,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.,0.,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(0.,Acs[0],1e-10);
  EXPECT_DOUBLE_EQ(Vc,Acs[1]);
  EXPECT_DOUBLE_EQ(0.,Acs[2]);
  VecDoub Angs = AA.angles(X);
  EXPECT_DOUBLE_EQ(0.,Angs[0]);
  EXPECT_DOUBLE_EQ(0.,Angs[1]);
  EXPECT_DOUBLE_EQ(0.,Angs[2]);
  EXPECT_DOUBLE_EQ(0.,Angs[3]);
  EXPECT_DOUBLE_EQ(Vc,Angs[4]);
  EXPECT_DOUBLE_EQ(0.,Angs[5]);
  X = {1.,0.,1.,0.,0.,0.};
  Acs = AA.actions(X);
  EXPECT_EQ(Acs[0]>0.,1);
  X = {1.,0.,0.,0.,0.,0.2};
  Acs = AA.actions(X);
  EXPECT_EQ(Acs[2]>0.,1);
  X = {1.,0.,0.,0.,0.,0.};
  EXPECT_EQ(AA.actions(X)[0]>0.,1);
}

TEST(ActionTest_Zero,StackelSAA){
  Logarithmic Pot(1.,1.,0.9);
  Actions_SpheroidalAdiabaticApproximation AA(&Pot,"",false,false,100.,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {0.,0.,0.,0.1*Vc,Vc,0.};
  VecDoub Acs = AA.actions(X);
  EXPECT_DOUBLE_EQ(0.,Acs[0]);
  EXPECT_DOUBLE_EQ(0.,Acs[1]);
  EXPECT_DOUBLE_EQ(0.,Acs[2]);
  VecDoub Angs = AA.angles(X);
  EXPECT_DOUBLE_EQ(0.,Acs[0]);
  EXPECT_DOUBLE_EQ(0.,Angs[1]);
  EXPECT_DOUBLE_EQ(0.,Angs[2]);
  EXPECT_DOUBLE_EQ(0.,Angs[3]);
  EXPECT_DOUBLE_EQ(0.,Angs[4]);
  EXPECT_DOUBLE_EQ(0.,Angs[5]);
}

TEST(ActionTest_Planar,StackelSAA){
  IsochronePotential Pot(1.,1.);
  Actions_SpheroidalAdiabaticApproximation AA(&Pot,"",false,false,-1.1,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  VecDoub X = {1.,0.,0.,0.1*Vc,Vc,0.};
  VecDoub Acs = AA.actions(X);
  Actions_Isochrone Iso(1.,1.);
  VecDoub ActsTrue = Iso.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],1e-10);
  EXPECT_DOUBLE_EQ(ActsTrue[1],Acs[1]);
  EXPECT_DOUBLE_EQ(ActsTrue[2],Acs[2]);
  VecDoub Angs = AA.angles(X);
  VecDoub AngsTrue = Iso.angles(X);
  EXPECT_NEAR(AngsTrue[0],Angs[0],2e-4);
  AngsTrue = Iso.freq(X);
  EXPECT_NEAR(AngsTrue[0],Angs[3],1e-4);
  EXPECT_NEAR(AngsTrue[1],Angs[4],1e-4);
}
TEST(ActionTest_Small,StackelSAA){
  Logarithmic Pot(1.,1.,0.9);
  Actions_SpheroidalAdiabaticApproximation AA(&Pot,"",false,false,100.,0.02,3.,2.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  std::cout<<"small z, zero vR, zero vz"<<std::endl;
  VecDoub X = {1.,2e-5,2e-5,0.,Vc,0.};
  VecDoub Acs = AA.actions(X);
  printVector(Acs);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  VecDoub Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"zero z, small vR, small vz"<<std::endl;
  X = {1.,0.,0.,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  printVector(Acs);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  printVector(Angs);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z, small vR, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);

  std::cout<<"small z, small vz"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);
  std::cout<<"small z"<<std::endl;
  X = {1.,2e-5,2e-5,.2*Vc,Vc,.3*Vc};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0] or Acs[0]!=Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1] or Acs[1]!=Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2] or Acs[2]!=Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0] or Angs[0]!=Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1] or Angs[1]!=Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2] or Angs[2]!=Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3] or Angs[3]!=Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4] or Angs[4]!=Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5] or Angs[5]!=Angs[5],false);

}}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
