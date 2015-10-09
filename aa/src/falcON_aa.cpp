// ============================================================================
/// \file src/falcON_aa.cpp
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
/// \brief Computing actions from nemo snapshot
///
/// Reads in snapshot and computes actions (designed for spherical potentials
/// but simple to extend)
///
// ============================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "main.h"
#include "GSLInterface/GSLInterface.h"
#include "utils.h"
#include "coordsys.h"
#include "coordtransforms.h"
#include "potential.h"
#include "spherical_aa.h"
#include "orbit.h"
#include "falcON_interface.h"
#include "fields.h"
#include "body.h"
// #include "getparam.h"
#include "nemo++.h"
#include <cstdlib>

#define falcON_VERSION   "0.8.2"
#define falcON_VERSION_D "13-apr-2010 Walter Dehnen                          "
//////////////////////////////////////////////////////////////////////////////
const char*defv[] = {
  "in=???\n           input snapshot                                     ",
  "out=???\n          output snapshot                                    ",
  "times=all\n        output time                                        ",
  "accname=???\n      name of gravitational potential                    ",
  "accpars=\n         parameters of gravitational potential              ",
  "accfile=\n         data file for gravitational potential              ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*usage =
  "Find actions from a snapshot given a ***spherical*** nemo acceleration"
  ;
//-----------------------------------------------------------------------------
namespace falcON{
class SphericalActionCalculator
{
  public:
  falcON_SphericalPotential *pot;
  Actions_Spherical *AS;
  std::ofstream     Out, LCfile;
  fieldset  Need;

  SphericalActionCalculator():
    pot(hasvalue("accname")?
       new falcON_SphericalPotential(
       getparam_z("accname"),
       getparam_z("accpars"),
       getparam_z("accfile")) : 0),
    AS(new Actions_Spherical(pot)),
    Need(fieldset::basic){}

  void operator()(snapshot&shot)
    {
      double Hmin=1e50,Hmax=-1e50, H;
      if(!Out.is_open()) Out.open(getparam("out"));
      VecDoub actions(2,0),angles(6,0),x(6,0);
      if(shot.N_bodies()) {
        LoopAllBodies(&shot,B){
          for(int i=0;i<3;i++){x[i]=pos(B)[i];x[i+3]=vel(B)[i];}
          actions = AS->actions(x);
          angles = AS->angles_and_freqs(x);
          Out<<mass(B)<<" ";
          for(int i=0;i<3;i++)Out<<pos(B)[i]<<" ";
          for(int i=0;i<3;i++)Out<<vel(B)[i]<<" ";
          // output JR and L
          Out<<actions[0]<<" "<<actions[2]+fabs(actions[1])<<" ";
          H=pot->H(x);
          for(auto i:angles)Out<<i<<" ";
          if(H>Hmax)Hmax=H;
          if(H<Hmin)Hmin=H;
          Out<<H<<std::endl;
        }
      }
      Out.close();
      LCfile.open(std::string(getparam("out"))+".Lc");
      if(LCfile.is_open()){
        for(auto i:create_range<double>(Hmin*1.005, Hmax*0.995, 200)){
          LCfile<<i<<" "<<pot->L_E(i)<<std::endl;
        }
        LCfile.close();
      }
    }
};
}
void falcON::main() falcON_THROWING
{
  nemo_in  In(getparam("in"));
  fieldset Read;
  snapshot Shot;
  SphericalActionCalculator ACF;

  if(!In.has_snapshot())
    falcON_THROW("no snapshots found in input file\n");

  if(0==strcmp(getparam("times"),"first")) {
    // special case times=first
    Shot.read_nemo(In,Read,ACF.Need,0,0);
    ACF(Shot);
  } else if(0==strcmp(getparam("times"),"last")) {
    // special case times=last
    while(In.has_snapshot())
      Shot.read_nemo(In,Read,ACF.Need,0,0);
    Shot.del_fields(~Read);
    ACF(Shot);
  } else {
    // general case for times
    while(In.has_snapshot())
      if(Shot.read_nemo(In,Read,ACF.Need,getparam("times"),0)) {
        Shot.del_fields(~Read);
        ACF(Shot);
      }
  }
  if(!ACF.Out.is_open())
    falcON_Warning("no snapshot matching \"times=%s\" found in input\n",
       getparam("times"));
}

// ============================================================================
