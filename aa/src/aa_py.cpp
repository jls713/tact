// ============================================================================
/// \file src/aa_py.cpp
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
/// \brief Boost python wrappers for action-finding routines
///
/// Wraps the action-finding routines such that they can be used in python
//============================================================================

#include <Python.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface/GSLInterface.h"
#include "gnuplot/gnuplot_i.h"
#include <gsl/gsl_poly.h>
#include "utils.h"
#include "coordsys.h"
#include "coordtransforms.h"
#include "potential.h"
#include "Multipole.h"
#include "orbit.h"
#include "aa.h"
#include "analytic_aa.h"
#include "stackel_aa.h"
#include "spherical_aa.h"
#include "genfunc_aa.h"
#include "adiabatic_aa.h"
#include "uv_orb.h"
#include "lmn_orb.h"
#include "stackel_fit.h"
#include "tables_aa.h"
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/list.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/docstring_options.hpp>
#include <numpy/arrayobject.h>
#ifdef TORUS
#include "falPot.h"
#include "it_torus.h"
#endif
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

using namespace boost::python;

// ==========================================================================
// Define converters for std::vector <==> numpy array

template <class T>
struct vector_to_ndarray {
  static PyObject* convert(const std::vector<T>& v) {
    list l;
    typename std::vector<T>::const_iterator p;
    for(p=v.begin();p!=v.end();++p){
      l.append(object(*p));
    }
    return incref(numeric::array(l).ptr());
  }
};

template <typename T>
struct vector_from_ndarray {
  vector_from_ndarray() {
    converter::registry::push_back(&vector_from_ndarray<T>::convertible,
                                   &vector_from_ndarray<T>::construct,
                                   type_id<std::vector<T>>());
  }

  // Determine if obj_ptr can be converted in a std::vector<T>
  static void* convertible(PyObject* obj_ptr) {
    if (!PyArray_Check(obj_ptr)) {
      std::cerr<<"You have passed a non-numpy array"<<std::endl;
      return 0;
    }
    return obj_ptr;
  }

  // Convert obj_ptr into a std::vector<T>
  static void construct(PyObject* obj_ptr,
                        converter::rvalue_from_python_stage1_data* data) {

    list l(handle<>(borrowed(obj_ptr)));
    // Grab pointer to memory into which to construct the new std::vector<T>
    void* storage = ((converter::rvalue_from_python_storage<std::vector<T>>*)
                     data)->storage.bytes;
    // in-place construct the new std::vector<T> using the character data
    // extraced from the python object
    std::vector<T>& v = *(new (storage) std::vector<T>());
    // populate the vector from list contains !!!
    int le = len(l);
    v.resize(le);
    for (int i = 0; i != le; ++i) {
      v[i] = extract<T>(l[i]);
    }
    // Stash the memory chunk pointer for later use by boost.python
    data->convertible = storage;
  }
};

// ==========================================================================
// For creating potentials in python

struct Potential_JS_PythonCallback : Potential_JS {
  Potential_JS_PythonCallback(PyObject* p) : self(p) {}
  double Phi(const VecDoub& x) {
    return call_method<double>(self, "Phi", static_cast<numeric::array>(x));
  }
  VecDoub Forces(const VecDoub& x) {
    VecDoub F(3, 0.);
    auto g = call_method<numeric::array>(
        self, "Forces", numeric::array(make_tuple(x[0], x[1], x[2])));
    for (unsigned i = 0; i < 3; ++i) F[i] = extract<double>(g[i]);
    return F;
  }
  PyObject* self;
};

// ==========================================================================
// Now define members of python library

BOOST_PYTHON_MODULE_INIT(aa_py) {
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  docstring_options doc_options(true);

// ============================================================================
// Potentials
// ============================================================================
  double (Potential_JS::*R_E)(const VecDoub&) = &Potential_JS::R_E;
  double (Potential_JS::*R_L)(const VecDoub&) = &Potential_JS::R_L;

  class_<Potential_JS, boost::noncopyable,
         boost::shared_ptr<Potential_JS_PythonCallback>>("Potential_JS",
            "Base potential class\n"
            "\n"
            "All potential classes must have a Phi(x) and Forces(x) routine in C++. These are called by __call__ and Forces respectively in python",init<>())
      .def("__call__", &Potential_JS::Phi,
           "Compute potential Phi(x)\n"
            "\n"
            "Args:\n"
            "    param1: np.array of Cartesian position x.\n"
            "\n"
            "Returns:\n"
            "    Potential at x\n"
        "")
      .def("energy", &Potential_JS::H,
           "Compute energy = \frac{1}{2}(v_x^2+v_y^2+v_z^2)+Phi(x)\n"
           "\n"
            "Args:\n"
            " param1: np.array of Cartesian position X.\n"
            "\n"
            "Returns:\n"
            " Potential at x\n"
        "")
      .def("Forces", &Potential_JS::Forces,
           "Compute Forces at x\n"
           "\n"
            "Args:\n"
            " param1: np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz).\n"
            "\n"
            "Returns:\n"
            " np.array of Forces at x\n"
        "")
      .def("torb", &Potential_JS::torb,
           "Compute orbital time-scale of circular orbit with energy H(x)\n"
           "\n"
            "Args:\n"
            " param1: np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz).\n"
            "\n"
            "Returns:\n"
            "    Orbital time-scale for x\n"
        "")
      .def("DeltaGuess", &Potential_JS::DeltaGuess,
           "Compute focal length (squared) for local approximation to St\"ackel potential \n"
           "\n"
            "Args:\n"
            " param1: np.array of Cartesian position X = (x,y,z).\n"
            "\n"
            "Returns:\n"
            "    focal length Delta^2 = gamma-alpha\n"
        "")
      .def("RCircular_Energy", R_E,
           "Compute radius of circular orbit with equivalent energy \n"
           "\n"
            "Args:\n"
            " param1: np.array of Cartesian position X = (x,y,z).\n"
            "\n"
            "Returns:\n"
            "    radius of circular orbit of energy\n"
        "")
      .def("RCircular_AngMom", R_L,
           "Compute radius of circular orbit with equivalent angular momentum \n"
           "\n"
            "Args:\n"
            " param1: np.array of Cartesian position X = (x,y,z).\n"
            "\n"
            "Returns:\n"
            "    radius of circular orbit of angular momentum\n"
        "");

  // class_<Density, boost::noncopyable>("Density",init<double>());

  // class_<TriaxialDensity, bases<Density>>("TriaxialDensity",init<VecDoub>());
  // class_<Density_NFW, bases<TriaxialDensity>>("Density_NFW",init<double, double,VecDoub,double>());

  // class_<TriaxialPotential, bases<Potential_JS>>("TriaxialPotential",init<TriaxialDensity*,double>());

  class_<SphericalPotential, bases<Potential_JS>>("SphericalPotential",
                                                  "Base class for spherical potentials",
                                                  no_init);
  class_<IsochronePotential, bases<SphericalPotential>>(
    "IsochronePotential",
    "Isochrone potential (spherical):\n\tPhi(r) = -frac{GM}{b+sqrt(b^2+r^2} \n\tTakes two parameters:\n\t\tG*mass: GM and the scale radius: b ",
    init<double, double>());
  class_<NFWSpherical, bases<SphericalPotential>>(
      "NFWSpherical",
      "Spherical NFW potential:Phi ="
      "-GM/r log(1+r/r_s)"
      "\n\tTakes two parameters:\n\t\tG*mass: GM and the scale radius: r_s",
                                                  init<double, double>());
  class_<HernquistSpherical, bases<SphericalPotential>>(
    "HernquistSpherical",
    "Spherical Hernquist potential:Phi ="
      "-GM/(r_s+r)"
    "\n\tTakes two parameters:\n\t\tG*mass: GM and the scale radius: r_s",
    init<double, double>());
  class_<PowerLawSpherical, bases<SphericalPotential>>(
         "PowerLawSpherical",
         "Power-law spherical potential:Phi ="
          "-GM/r^k"
        "\n\tTakes two parameters:\n\t\t"
        "G*mass: GM, the power: k",
        init<double, double>());
  class_<PowerLawSphericalScale, bases<SphericalPotential>>(
         "PowerLawSphericalScale",
         "Power-law spherical potential with scale: Phi ="
          "-v0^2/(1+(r/rs)^2)^(alpha/2)"
        "\n\tTakes three parameters:\n\t\t"
        "the amplitude: v0, the power: alpha and the scale: rs.",
        init<double, double, double>());

  class_<StackelOblate_PerfectEllipsoid,
	bases<Potential_JS>>("StackelOblate_PerfectEllipsoid","Axisymmetric oblate perfect ellipsoid.\nThis potential is of St\"ackel form and is most simply defined by its density:\n\trho(R,z) = -frac{rho0}{1+m^2}^2\n\tm^2 = (R/a)^2+(z/b)^2\n\tTakes two parameters:\t the central density: rho0 and the parameter alpha such that ... hmm   ",init<double,double>());

  class_<StackelTriaxial,
	bases<Potential_JS>>("StackelTriaxial","Potential of perfect ellipsoid: takes three parameters: central density, alpha and beta", init<double,double,double>());

  class_<Logarithmic, bases<Potential_JS>>("Logarithmic",
                                           "Logarithmic potential:Phi ="
      "(1/2)V_c^2\\log(x^2+(y/qy)^2+(z/qz)^2)"
    "\n\tTakes three parameters:\n\t\t"
    "circular velocity: VC, and the axis ratios q_y, q_z",
                                           init<double, double, double>());

  class_<PowerLaw, bases<Potential_JS>>("PowerLaw",
                                        "Power-law potential:Phi ="
      "-GM/(x^2+y^2/q_y^2+z^2/q_z^2)^k"
    "\n\tTakes four parameters:\n\t\t"
    "G*mass: GM, the power: k and the axis ratios q_y, q_z",
                                        init<double, double, double, double>());
  class_<Isochrone, bases<Potential_JS>>(
      "Isochrone",
      "Isochrone potential:Phi ="
      "-GM/(b+\\sqrt{b^2+r^2}) where r^2 = (x^2+y^2/q_y^2+z^2/q_z^2)"
    "\n\tTakes four parameters:\n\t\t"
    "G*mass: GM, the scale: b and the axis ratios q_y, q_z"
    ,init<double, double, double, double>());

  class_<HarmonicOscillator, bases<Potential_JS>>(
      "HarmonicOscillator",
      "Triaxial harmonic oscillator potential:Phi ="
      "\\sum_i .5*\\omega_i x_i x_i"
    "\n\tTakes three parameter: the frequencies Om in each direction",
    init<double, double, double>());

#ifdef TORUS
  class_<GalPot, bases<Potential_JS>>("GalPot","Wrapper of Torus galaxy potential -- pass input Tpot file.", init<std::string>());
#endif

  class_<BowdenNFW, bases<Potential_JS>>("BowdenNFW",
  "Triaxial NFW from Bowden, Evans & Belokurov (2014)\n"
  "\tTakes 6 parameters: the usual central density (rho0) and scale (rs),"
  "and the z and y flattenings at zero(q0,p0) and inf(qinf,pinf)",init<double,double,double,double,double,double>());
  class_<MWPotential2014, bases<Potential_JS>>("MWPotential2014",
  "Bovy's MWPotential2014 (2014)\n"
  "\tPower Law bulge with exp cut off, Miyamoto-Nagai disc, NFW halo",init<>());

// ============================================================================
// Action Finders
// ============================================================================
  class_<Action_Finder>("Action_Finder",
        "Base action finding class\n"
        "\n"
        "    All action finding classes must have actions(x) and angles(x) routine that returns the actions and (angles,frequencies) respectively.\n"
        "\n"
        "",no_init)
	.def("actions",&Action_Finder::actions_nopars,
       "Compute actions(X)\n"
       "\n"
       "  Args:\n"
       "    param1: np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz).\n"
       "\n"
       "  Returns:\n"
       "    np.array of actions = (J_R,J_phi,J_z)\n"
       "")
  .def("angles",&Action_Finder::angles_nopars,
       "Compute angles(X)\n"
       "\n"
       "  Args:\n"
       "    param1: np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz).\n"
       "\n"
       "  Returns:\n"
       "    np.array of angles and frequencies = (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)\n"
       "");

  class_<Actions_Spherical, bases<Action_Finder> >(
	"Actions_Spherical","Action finding in spherical potential", init<SphericalPotential*>());

  class_<Actions_Isochrone, bases<Action_Finder, Isochrone> >(
  "Actions_Isochrone",
  "Action finder for isochrone potential\n"
    "Args:\n"
      "param1: potential param GM\n"
      "param2: potential param b\n"
    "Finds actions in spherical potential Phi(r)=-GM/(b+sqrt(b^2+r^2))"
    , init<double,double>());

  class_<Actions_CylindricalAdiabaticApproximation, bases<Action_Finder> >(
	"Actions_CylindricalAdiabaticApproximation",
  "Cylindrical Adiabatic Approximation\n"
    "Args:\n"
     " param1: Potential_JS potential\n"
      "param2: output file\n"
      "param3: write to file bool (if 1 write)\n"
      "param4: no_energy_corr (if true don't couple radial and vertical energies)\n"
      "param5 Rm: minimum radius for grid\n"
      "param6 Rn: maximum radius for grid\n"
      "param7 zmax: maximum z height for grid\n"
      "param8 NGRID: number of grid points for Ez\n"
      "", init<Potential_JS*, std::string,bool,bool,double,double,double,int>());

  class_<Actions_SpheroidalAdiabaticApproximation, bases<Action_Finder> >(
	"Actions_SpheroidalAdiabaticApproximation",
  "Spheroidal Adiabatic Approximation:"
    "can pass alpha guess to actions and angles, or if alpha>1 uses local guess\n"
    "Args:\n"
      "param1: Potential_JS potential\n"
      "param2: output file\n"
     " param3: write to file bool (if 1 write)\n"
      "param4: no_energy_corr (if true don't couple radial and vertical energies)\n"
      "param5: alpha estimate\n"
      "param6 Rm -- minimum radius for grid\n"
      "param7 Rn -- maximum radius for grid\n"
      "param8 zmax -- maximum z height for grid\n"
      "param9 NGRID -- number of grid points for Enu\n"
      "param10 NL -- number of grid points for Lz\n"
      "", init<Potential_JS*, std::string,bool,bool,double,double,double,double,int,int>());

  class_<Actions_AxisymmetricStackel, bases<Action_Finder> >(
	"Actions_AxisymmetricStackel","Action finding in axisymmetric Staeckel potential:\n"
      "params: StackelOblate_PerfectEllipsoid potential\n"
      "", init<StackelOblate_PerfectEllipsoid*>());

  class_<Actions_AxisymmetricStackel_Fudge, bases<Action_Finder> >(
      "Actions_AxisymmetricStackel_Fudge"," Action estimation in general axisymmetric potential using Staeckel fudge (can pass alpha guess to actions and angles, or if alpha>1 uses local guess:\n"
        "Args:\n"
          "param1 Potential\n"
          "param2 alpha guess"
      , init<Potential_JS*, double>());

  class_<Actions_TriaxialStackel, bases<Action_Finder> >(
	"Actions_TriaxialStackel",
  "Action finding in triaxial Staeckel potential:\n"
     " Args:\n"
      "  param1: StackelTriaxial potential\n"
      "", init<StackelTriaxial*>());

  class_<Actions_TriaxialStackel_Fudge, bases<Action_Finder> >(
	"Actions_TriaxialStackel_Fudge",
  "Action estimation in general triaxial potentials using Staeckel fudge\n"
   " Args:\n"
      "param1: Potential_JS potential\n"
      "param2: alpha\n"
     " param3: beta\n"
  "", init<Potential_JS*, double, double>());

  class_<Actions_Genfunc, bases<Action_Finder> >(
	"Actions_Genfunc",
  ""
    "Action finder using Generating function\n"
    "Args:\n"
      "param1: Potential_JS potential\n"
      "param2: symmetry string (axisymmetric or triaxial)\n"
    "Finds actions by calculating the generating function components from a orbit integration via toy actions in either the isochrone or harmonic oscillator potential\n"
    "The actions routine accepts a set of parameters = (Total time, N_T, N_max). If no parameters are passed the routine will estimate the orbital time scale from T = Potential_JS::torb and use Total_T = 8*T, N_T=200 and N_max = 6\n"
    "Currently geared up to return J_lambda, J_mu, J_nu -- this means that for the inner long axis loops (identified via the Stackel fudge) have J_r <==> L such that J_mu is the radial action\n"
  "", init<Potential_JS*,std::string>())
  .def("full_actions",&Actions_Genfunc::full_actions,
       ""
       "Finds actions\n"
        "Args:\n"
          "param1 np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz).\n"
          "param2 NT number of orbital times to integrate for\n"
          "param3 Nsamp number of time samples\n"
         " param4 Nmax number of Fourier coefficients Sn\n"
         "\n"
        "Returns: actions -- 3D vector J=(J_R,J_phi,J_z)\n"
        "\n")
  .def("full_angles",&Actions_Genfunc::full_angles,
       ""
       "Finds angles\n"
        "Args:\n"
          "param1 np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz).\n"
          "param2 NT number of orbital times to integrate for\n"
          "param3 Nsamp number of time samples\n"
         " param4 Nmax number of Fourier coefficients Sn\n"
         "\n"
        "Returns: angles and frequencies -- 6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)\n"
        "\n");

  class_<Actions_Genfunc_Average, bases<Action_Finder> >(
  "Actions_Genfunc_Average\n",
  ""
    "Action finder averaging the toy actions over toy angles\n"
    "Args:\n"
      "param1: Potential_JS potential\n"
      "param2: symmetry string (axisymmetric or triaxial)\n"
    "Finds actions by averaging using an orbit integration via toy actions in either the isochrone or harmonic oscillator potential\n"
    "\n", init<Potential_JS*,std::string>());

  class_<Actions_StackelFit, bases<Action_Finder> >(
  "Actions_StackelFit\n",
  ""
    "Action finding by first fitting a Staeckel potential to the region an orbit probes\n"
    "Args:\n"
      "param1: Potential_JS potential\n"
      "param2: tolerance of orbit integration\n"
  "", init<Potential_JS*,double>());

#ifdef TORUS
  class_<IterativeTorusMachine, bases<Action_Finder> >(
      "IterativeTorusMachine\n",
      "Action finding by iteratively constructing tori\n"
      "Args:\n"
       " param1 Actions_AxisymmetricStackel_Fudge object\n"
        "param2 GalPot potential\n"
        "param3 tolerance in angle fit\n"
       " param4 max number of torus fits\n"
        "param5 relative error of torus fits\n",
      init<Actions_AxisymmetricStackel_Fudge*, GalPot*, double, int, double>()
      ).def("set_maxit",&IterativeTorusMachine::set_maxit,
      "set max no. of torus constructions");
#endif

  class_<uv_orb, bases<Action_Finder> >(
      "Actions_AxisymmetricStackel_Fudge_DeltaGuess",
      "uv_orb is an interface for using axisymmetric St\"ackel fudge apparatus. It handles the choice of coordinate system using the Delta estimation from Binney (2014). It is thread-safe as each action call creates a new instance of Actions_AxisymmetricStackel_Fudge\n"
"  Args:\n"
    "param pot Potential_JS (axisymmetric)\n"
    "param Rm  minimum radial grid point\n"
    "param Rn  maximum radial grid point\n"
   " param NE  number of energy grid points\n"
    "param NL  number of ang. mom. grid points\n"
    "param name name of class\n"
    "\n",
      init<Potential_JS*, double, double, int, int, std::string>());

  class_<Actions_AxisymmetricFudge_InterpTables, bases<Action_Finder> >(
      "Actions_AxisymmetricFudge_InterpTables",
      "Interpolation tables for actions using scheme from Binney (2014)\n"
"  Args:\n"
    "param pot Potential_JS (axisymmetric)\n"
    "param name filename for tables\n"
    "param tab tabulate actions (False) or read from file (True)\n"
    "param Rm  minimum radial grid point\n"
    "param Rn  maximum radial grid point\n"
   " param NE  number of radial grid points\n"
    "param NGRID  number of energy grid points\n"
   " param NED  number of delta energy grid points\n"
    "param NEL  number of delta ang. mom. grid points\n"
    "\n",
      init<Potential_JS*, std::string, optional<bool, double, double, int, int, int, int>>());

  class_<lmn_orb, bases<Action_Finder> >(
      "Actions_TriaxialStackel_Fudge_DeltaGuess",
      ""
      "Interface for using triaxial St\"ackel fudge apparatus\n"
      "Handles the choice of coordinate system\n"
     " Args:\n"
          "param pot Potential (axisymmetric) in which to compute the actions\n"
          "param ym -- minimum intermediate (y) intercept value\n"
         " param yn -- maximum intermediate (y) intercept value\n"
          "param NEE -- number of energy grid points\n"
          "param use_log_grid -- if true, use log-spaced grid, else linear\n"
          "param use_acts-if true, estimate alpha & beta by action minimization\n"
          "param s -- output file\n"
      "\n",
      init<Potential_JS*, double, double, int, bool,bool, std::string>());


  class_<Orbit, boost::noncopyable>("Orbit",
                                    "Orbit integrator using GSL ode routines"
                                    "\nTakes two args: potential and accuracy", init<Potential_JS*, double>())
      .def("integrate", &Orbit::integrate,
           "Integrate phase-space point.\n"
              "Args:\n"
             " param1 np.array x -- starting phase-space point x=(x,y,z,vx,vy,vz)\n"
              "param2 t_interval -- total integration time\n"
              "param3 step_size -- integration output step size\n"
             " param4 adaptive -- adaptive step-sizes\n"
             ""
              "Returns: final phase-space point\n"
              ""
              "intermediate results stored in results\n"
            "");


  def("GalacticToPolar",conv::GalacticToPolar,
      "\n"
        "in: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)\n"
       " out: Polar = (R,phi,z,vR,vphi,vz)\n"
"\n"
        "Units: velocities in km/s, positions kpc and proper motions in mas/yr\n"
      "\n");
  def("GalacticToCartesian",conv::GalacticToCartesian,
      "\n"
        "in: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)\n"
        "out: Cartesian = (x,y,z,vx,vy,vz):\n"
           "  x away from centre of Galaxy towards Sun\n"
            " y opposite to Galactic rotation\n"
             "z up towards North Galactic Pole\n"
"\n"
        "Units: velocities in km/s, positions kpc and proper motions in mas/yr\n"
      "\n");
  def("GalacticToEquatorial",conv::GalacticToEquatorial,
      "\n"
       " in: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)\n"
       " out: Equatorial =(alpha,delta,s,vlos,mu_alphacos(delta),mu_delta)\n"
      "\n");
  def("CartesianToPolar",conv::CartesianToPolar,
      "\n"
       " in:  Cartesian = (x,y,z,vx,vy,vz)\n"
       " out: Polar = (R,phi,z,vR,vphi,vz)\n"
      "\n");
  def("PolarToCartesian",conv::PolarToCartesian,
      "\n"
       " in:  Polar = (R,phi,z,vR,vphi,vz)\n"
       " out: Cartesian = (x,y,z,vx,vy,vz)\n"
      "\n");
  def("CartesianToGalactic",conv::CartesianToGalactic,
      "\n"
       " in:  Cartesian = (x,y,z,vx,vy,vz):\n"
       "      x away from centre of Galaxy towards Sun\n"
        "     y opposite to Galactic rotation\n"
       "      z up towards North Galactic Pole\n"
      "  out: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)\n"
"\n"
      "  Units: velocities in km/s, positions kpc and proper motions in mas/yr\n"
      "\n");
  def("PolarToGalactic",conv::PolarToGalactic,
      "\n"
      "  in:  Polar = (R,phi,z,vR,vphi,vz)\n"
      "  out: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)\n"
      "\n"
      "  Units: velocities in km/s, positions kpc and proper motions in mas/yr\n"
      "\n");
  def("EquatorialToGalactic",conv::EquatorialToGalactic,
      "\n"
      "  in: Equatorial =(alpha,delta,s,vlos,mu_alphacos(delta),mu_delta)\n"
      "  out: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)\n"
      "\n");
  def("EquatorialToGalacticwithErrors",conv::EquatorialToGalacticwithErrors,
      "\n"
      "  in: Equatorial =(alpha,delta,s,vlos,mu_alphacos(delta),mu_delta)\n"
      "  in: EqErrs =(e_alpha,e_delta,e_s,e_vlos,e_mu_alphacos(delta),e_mu_delta)\n"
      "  out: Galactic = ((l,b,s,vlos,mu_lcos(b),mu_b),(e_l,e_b,e_s,e_vlos,e_mu_lcos(b),e_mu_b))\n"
      "\n");


  to_python_converter<VecDoub, vector_to_ndarray<double>>();
  vector_from_ndarray<double>();
  to_python_converter<std::vector<VecDoub>, vector_to_ndarray<VecDoub>>();
  vector_from_ndarray<VecDoub>();
  to_python_converter<std::vector<Potential_JS*>, vector_to_ndarray<Potential_JS*>>();
  vector_from_ndarray<Potential_JS*>();
  import_array();
}

// ==========================================================================

