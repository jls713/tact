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

TEST(ActionTest,StackelFit){
  // Logarithmic Pot(220.,1.,0.9);
  Isochrone Pot(1.,1.);
  Actions_StackelFit AA(&Pot);
  VecDoub X = {8.29,0.1,0.1,40.,200.,50.};
  VecDoub ActsTrue = {29.357028436047674,1653.9999999999998,44.218614307037832};
  VecDoub AngsTrue = {1.8104951133863227,6.0951452351556137,6.1521997033776064,39.443007528837896,27.741995155048556,30.810000871423856};
  // VecDoub Angs = AA.angles(X);
  VecDoub Acs = AA.actions(X);
  EXPECT_NEAR(ActsTrue[0],Acs[0],1e-10);
  EXPECT_NEAR(ActsTrue[1],Acs[1],1e-10);
  EXPECT_NEAR(ActsTrue[2],Acs[2],1e-10);
  // EXPECT_DOUBLE_EQ(AngsTrue[0],Angs[0]);
  // EXPECT_DOUBLE_EQ(AngsTrue[1],Angs[1]);
  // EXPECT_DOUBLE_EQ(AngsTrue[2],Angs[2]);
  // EXPECT_DOUBLE_EQ(AngsTrue[3],Angs[3]);
  // EXPECT_DOUBLE_EQ(AngsTrue[4],Angs[4]);
  // EXPECT_DOUBLE_EQ(AngsTrue[5],Angs[5]);
}


}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
