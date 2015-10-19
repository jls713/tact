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
#include "../pot/inc/orbit.h"
#include "../aa/inc/analytic_aa.h"
#include "../aa/inc/adiabatic_aa.h"
#include "../aa/inc/stackel_aa.h"
#include "../aa/inc/stackel_fit.h"
#include "gtest/gtest.h"

namespace {

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
  Logarithmic Pot(1.,1.,0.99);
  Actions_AxisymmetricStackel_Fudge AA(&Pot,1.);
  double radius = 1.;
  double Vc = sqrt(radius*-Pot.Forces({1.,0.,0.})[0]);
  // small z, zero vR, zero vz
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
  // zero z, small vR, small vz
  X = {1.,0.,0.,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5],false);
  // small z, small vR, small vz
  X = {1.,2e-5,2e-5,2e-5,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5],false);
  // small z, small vz
  X = {1.,2e-5,2e-5,.2*Vc,Vc,2e-5};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5],false);
  // small z
  X = {1.,2e-5,2e-5,.2*Vc,Vc,.3*Vc};
  Acs = AA.actions(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Acs[2],false);
  Angs = AA.angles(X);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[0],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[1],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[2],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[3],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[4],false);
  EXPECT_EQ(std::numeric_limits<double>::infinity()==Angs[5],false);
}
//=============================================================================

TEST(ActionTest_Iso,StackelPAA){
  Isochrone Pot(1.,1.,0.99999);
  Actions_PolarAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
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
  EXPECT_NEAR(FreqsTrue[2],Angs[5],0.004);
}

// NANs in Stackel forces on axis

// TEST(ActionTest_Stack,StackelPAA){
//   StackelOblate_PerfectEllipsoid Pot(1.,-30.);
//   Actions_PolarAdiabaticApproximation AA(&Pot,"",false,false,0.2,3.,2.);
//   Actions_AxisymmetricStackel AAS(&Pot);
//   double radius = 1.;
//   double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
//   VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
//   VecDoub ActsTrue = AAS.actions(X);
//   VecDoub FreqsTrue = AAS.angles(X);
//   VecDoub Angs = AA.angles(X);
//   VecDoub Acs = AA.actions(X);
//   EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
//   EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
//   EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
//   EXPECT_NEAR(FreqsTrue[3],Angs[3],0.001);
//   EXPECT_NEAR(FreqsTrue[4],Angs[4],0.001);
//   EXPECT_NEAR(FreqsTrue[5],Angs[5],0.004);
// }

TEST(ActionTest_Zeros,StackelPAA){
  Logarithmic Pot(1.,1.,0.9);
  Actions_PolarAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
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
}


TEST(ActionTest_Planar,StackelPAA){
  IsochronePotential Pot(1.,1.);
  Actions_PolarAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
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

TEST(ActionTest_Zero,StackelPAA){
  IsochronePotential Pot(1.,1.);
  Actions_PolarAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
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

TEST(ActionTest_Small,StackelPAA){
  IsochronePotential Pot(1.,1.);
  Actions_PolarAdiabaticApproximation AA(&Pot,"",false,false,0.02,3.,2.);
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
//=============================================================================
// TEST(ActionTest_Iso,StackelFit){
//   Isochrone Pot(1.,1.,0.99999);
//   Actions_StackelFit AA(&Pot);
//   double radius = 1.;
//   double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
//   VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
//   Actions_Isochrone Iso(1.,1.);
//   VecDoub ActsTrue = Iso.actions(X);
//   VecDoub FreqsTrue = Iso.freq(X);
//   VecDoub Angs = AA.angles(X);
//   VecDoub Acs = AA.actions(X);
//   EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
//   EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
//   EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
//   EXPECT_NEAR(FreqsTrue[0],Angs[3],0.001);
//   EXPECT_NEAR(FreqsTrue[1],Angs[4],0.001);
//   EXPECT_NEAR(FreqsTrue[2],Angs[5],0.004);
// }

// TEST(ActionTest_Stack,StackelFit){
//   StackelOblate_PerfectEllipsoid Pot(100.,-30.);
//   Actions_StackelFit AA(&Pot);
//   Actions_AxisymmetricStackel AAS(&Pot);
//   double radius = 1.;
//   double Vc = sqrt(radius*-Pot.Forces({1.,0.1,0.1})[0]);
//   VecDoub X = {1.,0.1,0.1,0.1*Vc,Vc,0.1*Vc};
//   VecDoub ActsTrue = AAS.actions(X);
//   VecDoub FreqsTrue = AAS.angles(X);
//   VecDoub Angs = AA.angles(X);
//   VecDoub Acs = AA.actions(X);
//   EXPECT_NEAR(ActsTrue[0],Acs[0],0.001);
//   EXPECT_NEAR(ActsTrue[1],Acs[1],0.001);
//   EXPECT_NEAR(ActsTrue[2],Acs[2],0.001);
//   EXPECT_NEAR(FreqsTrue[0],Angs[3],0.001);
//   EXPECT_NEAR(FreqsTrue[1],Angs[4],0.001);
//   EXPECT_NEAR(FreqsTrue[2],Angs[5],0.004);
// }


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

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
