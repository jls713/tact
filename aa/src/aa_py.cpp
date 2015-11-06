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

#include <python2.7/Python.h>
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
  class_<Potential_JS, boost::noncopyable,
         boost::shared_ptr<Potential_JS_PythonCallback>>("Potential_JS",
            "Base potential class"
            ""
            "All potential classes must have a Phi(x) and Forces(x) routine in C++. These are called by __call__ and Forces respectively in python",init<>())
      .def("__call__", &Potential_JS::Phi,
           "Compute potential Phi(x)"
            ""
            "Args:"
            "    param1: np.array of Cartesian position x."
            ""
            "Returns:"
            "    Potential at x"
        "")
      .def("energy", &Potential_JS::H,
           "Compute energy = \frac{1}{2}(v_x^2+v_y^2+v_z^2)+Phi(x)"
           ""
            "Args:"
            " param1: np.array of Cartesian position X."
            ""
            "Returns:"
            " Potential at x"
        "")
      .def("Forces", &Potential_JS::Forces,
           "Compute Forces at x"
           ""
            "Args:"
            " param1: np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz)."
            ""
            "Returns:"
            " np.array of Forces at x"
        "")
      .def("torb", &Potential_JS::torb,
           "Compute orbital time-scale of circular orbit with energy H(x)"
           ""
            "Args:"
            " param1: np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz)."
            ""
            "Returns:"
            "    Orbital time-scale for x"
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

// ============================================================================
// Action Finders
// ============================================================================
  class_<Action_Finder>("Action_Finder",
        "Base action finding class"
        ""
        "    All action finding classes must have actions(x) and angles(x) routine that returns the actions and (angles,frequencies) respectively."
        ""
        "",no_init)
	.def("actions",&Action_Finder::actions_nopars,
       "Compute actions(X)"
       ""
       "  Args:"
       "    param1: np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz)."
       ""
       "  Returns:"
       "    np.array of actions = (J_R,J_phi,J_z)"
       "")
  .def("angles",&Action_Finder::angles_nopars,
       "Compute angles(X)"
       ""
       "  Args:"
       "    param1: np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz)."
       ""
       "  Returns:"
       "    np.array of angles and frequencies = (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)"
       "");

  class_<Actions_Spherical, bases<Action_Finder> >(
	"Actions_Spherical","Action finding in spherical potential", init<SphericalPotential*>());

  class_<Actions_Isochrone, bases<Action_Finder, Isochrone> >(
  "Actions_Isochrone",
  "Action finder for isochrone potential"
    "Args:"
      "param1: potential param GM"
      "param2: potential param b"
    "Finds actions in spherical potential Phi(r)=-GM/(b+sqrt(b^2+r^2))"
    , init<double,double>());

  class_<Actions_PolarAdiabaticApproximation, bases<Action_Finder> >(
	"Actions_PolarAdiabaticApproximation",
  "Polar Adiabatic Approximation"
    "Args:"
     " param1: Potential_JS potential"
      "param2: output file"
      "param3: write to file bool (if 1 write)"
      "param4: no_energy_corr (if true don't couple radial and vertical energies)"
      "param5 Rm: minimum radius for grid"
      "param6 Rn: maximum radius for grid"
      "param7 zmax: maximum z height for grid"
      "param8 NGRID: number of grid points for Ez"
      "", init<Potential_JS*, std::string,bool,bool,double,double,double,int>());

  class_<Actions_SpheroidalAdiabaticApproximation, bases<Action_Finder> >(
	"Actions_SpheroidalAdiabaticApproximation",
  "Spheroidal Adiabatic Approximation:"
    "can pass alpha guess to actions and angles, or if alpha>1 uses local guess"
    "Args:"
      "param1: Potential_JS potential"
      "param2: output file"
     " param3: write to file bool (if 1 write)"
      "param4: no_energy_corr (if true don't couple radial and vertical energies)"
      "param5: alpha estimate"
      "param6 Rm -- minimum radius for grid"
      "param7 Rn -- maximum radius for grid"
      "param8 zmax -- maximum z height for grid"
      "param9 NGRID -- number of grid points for Enu"
      "param10 NL -- number of grid points for Lz"
      "", init<Potential_JS*, std::string,bool,bool,double,double,double,double,int,int>());

  class_<Actions_AxisymmetricStackel, bases<Action_Finder> >(
	"Actions_AxisymmetricStackel","Action finding in axisymmetric Staeckel potential:"
      "params: StackelOblate_PerfectEllipsoid potential"
      "", init<StackelOblate_PerfectEllipsoid*>());

  class_<Actions_AxisymmetricStackel_Fudge, bases<Action_Finder> >(
      "Actions_AxisymmetricStackel_Fudge"," Action estimation in general axisymmetric potential using Staeckel fudge (can pass alpha guess to actions and angles, or if alpha>1 uses local guess:"
        "Args:"
          "param1 Potential"
          "param2 alpha guess"
      , init<Potential_JS*, double>());

  class_<Actions_TriaxialStackel, bases<Action_Finder> >(
	"Actions_TriaxialStackel",
  "Action finding in triaxial Staeckel potential:"
     " Args:"
      "  param1: StackelTriaxial potential"
      "", init<StackelTriaxial*>());

  class_<Actions_TriaxialStackel_Fudge, bases<Action_Finder> >(
	"Actions_TriaxialStackel_Fudge",
  "Action estimation in general triaxial potentials using Staeckel fudge"
   " Args:"
      "param1: Potential_JS potential"
      "param2: alpha"
     " param3: beta"
  "", init<Potential_JS*, double, double>());

  class_<Actions_Genfunc, bases<Action_Finder> >(
	"Actions_Genfunc",
  ""
    "Action finder using Generating function"
    "Args:"
      "param1: Potential_JS potential"
      "param2: symmetry string (axisymmetric or triaxial)"
    "Finds actions by calculating the generating function components from a orbit integration via toy actions in either the isochrone or harmonic oscillator potential"
    "The actions routine accepts a set of parameters = (Total time, N_T, N_max). If no parameters are passed the routine will estimate the orbital time scale from T = Potential_JS::torb and use Total_T = 8*T, N_T=200 and N_max = 6"
    "Currently geared up to return J_lambda, J_mu, J_nu -- this means that for the inner long axis loops (identified via the Stackel fudge) have J_r <==> L such that J_mu is the radial action"
  "", init<Potential_JS*,std::string>())
  .def("full_actions",&Actions_Genfunc::full_actions,
       ""
       "Finds actions"
        "Args:"
          "param1 np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz)."
          "param2 NT number of orbital times to integrate for"
          "param3 Nsamp number of time samples"
         " param4 Nmax number of Fourier coefficients Sn"
         ""
        "Returns: actions -- 3D vector J=(J_R,J_phi,J_z)"
        "")
  .def("full_angles",&Actions_Genfunc::full_angles,
       ""
       "Finds angles"
        "Args:"
          "param1 np.array of Cartesian position and velocity X = (x,y,z,vx,vy,vz)."
          "param2 NT number of orbital times to integrate for"
          "param3 Nsamp number of time samples"
         " param4 Nmax number of Fourier coefficients Sn"
         ""
        "Returns: angles and frequencies -- 6D vector (theta_R,theta_phi,theta_z,Omega_R,Omega_phi,Omega_z)"
        "");

  class_<Actions_Genfunc_Average, bases<Action_Finder> >(
  "Actions_Genfunc_Average",
  ""
    "Action finder averaging the toy actions over toy angles"
    "Args:"
      "param1: Potential_JS potential"
      "param2: symmetry string (axisymmetric or triaxial)"
    "Finds actions by averaging using an orbit integration via toy actions in either the isochrone or harmonic oscillator potential"
    "", init<Potential_JS*,std::string>());

  class_<Actions_StackelFit, bases<Action_Finder> >(
  "Actions_StackelFit",
  ""
    "Action finding by first fitting a Staeckel potential to the region an orbit probes"
    "Args:"
      "param1: Potential_JS potential"
      "param2: tolerance of orbit integration"
  "", init<Potential_JS*,double>());

#ifdef TORUS
  class_<IterativeTorusMachine, bases<Action_Finder> >(
      "IterativeTorusMachine",
      "Action finding by iteratively constructing tori"
      "Args:"
       " param1 Actions_AxisymmetricStackel_Fudge object"
        "param2 GalPot potential"
        "param3 tolerance in angle fit"
       " param4 max number of torus fits"
        "param5 relative error of torus fits",
      init<Actions_AxisymmetricStackel_Fudge*, GalPot*, double, int, double>()
      ).def("set_maxit",&IterativeTorusMachine::set_maxit,
      "set max no. of torus constructions");
#endif

  class_<uv_orb, bases<Action_Finder> >(
      "Actions_AxisymmetricStackel_Fudge_DeltaGuess",
      "uv_orb is an interface for using axisymmetric St\"ackel fudge apparatus. It handles the choice of coordinate system using the Delta estimation from Binney (2014). It is thread-safe as each action call creates a new instance of Actions_AxisymmetricStackel_Fudge"
"  Args:"
    "param pot Potential_JS (axisymmetric)"
    "param Rm  minimum radial grid point"
    "param Rn  maximum radial grid point"
   " param NE  number of energy grid points"
    "param NL  number of ang. mom. grid points"
    "param name name of class"
    "",
      init<Potential_JS*, double, double, int, int, std::string>());

  class_<lmn_orb, bases<Action_Finder> >(
      "Actions_TriaxialStackel_Fudge_DeltaGuess",
      ""
      "Interface for using triaxial St\"ackel fudge apparatus"
      "Handles the choice of coordinate system"
     " Args:"
          "param pot Potential (axisymmetric) in which to compute the actions"
          "param ym -- minimum intermediate (y) intercept value"
         " param yn -- maximum intermediate (y) intercept value"
          "param NEE -- number of energy grid points"
          "param use_log_grid -- if true, use log-spaced grid, else linear"
          "param use_acts-if true, estimate alpha & beta by action minimization"
          "param s -- output file"
      "",
      init<Potential_JS*, double, double, int, bool,bool, std::string>());


  class_<Orbit, boost::noncopyable>("Orbit","Orbit integrator using GSL ode routines", init<Potential_JS*, double>())
      .def("integrate", &Orbit::integrate,
           "Integrate phase-space point."
              "Args:"
             " param1 np.array x -- starting phase-space point x=(x,y,z,vx,vy,vz)"
              "param2 t_interval -- total integration time"
              "param3 step_size -- integration output step size"
             " param4 adaptive -- adaptive step-sizes"
             ""
              "Returns: final phase-space point"
              ""
              "intermediate results stored in results"
            "");


  def("GalacticToPolar",conv::GalacticToPolar,
      ""
        "in: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)"
       " out: Polar = (R,phi,z,vR,vphi,vz)"
""
        "Units: velocities in km/s, positions kpc and proper motions in mas/yr"
      "");
  def("GalacticToCartesian",conv::GalacticToCartesian,
      ""
        "in: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)"
        "out: Cartesian = (x,y,z,vx,vy,vz):"
           "  x away from centre of Galaxy towards Sun"
            " y opposite to Galactic rotation"
             "z up towards North Galactic Pole"
""
        "Units: velocities in km/s, positions kpc and proper motions in mas/yr"
      "");
  def("GalacticToEquatorial",conv::GalacticToEquatorial,
      ""
       " in: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)"
       " out: Equatorial =(alpha,delta,s,vlos,mu_alphacos(delta),mu_delta)"
      "");
  def("CartesianToPolar",conv::CartesianToPolar,
      ""
       " in:  Cartesian = (x,y,z,vx,vy,vz)"
       " out: Polar = (R,phi,z,vR,vphi,vz)"
      "");
  def("CartesianToGalactic",conv::CartesianToGalactic,
      ""
       " in:  Cartesian = (x,y,z,vx,vy,vz):"
       "      x away from centre of Galaxy towards Sun"
        "     y opposite to Galactic rotation"
       "      z up towards North Galactic Pole"
      "  out: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)"
""
      "  Units: velocities in km/s, positions kpc and proper motions in mas/yr"
      "");
  def("PolarToGalactic",conv::PolarToGalactic,
      ""
      "  in:  Polar = (R,phi,z,vR,vphi,vz)"
      "  out: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)"
      ""
      "  Units: velocities in km/s, positions kpc and proper motions in mas/yr"
      "");
  def("EquatorialToGalactic",conv::EquatorialToGalactic,
      ""
      "  in: Equatorial =(alpha,delta,s,vlos,mu_alphacos(delta),mu_delta)"
      "  out: Galactic = (l,b,s,vlos,mu_lcos(b),mu_b)"
      "");



  to_python_converter<VecDoub, vector_to_ndarray<double>>();
  vector_from_ndarray<double>();
  import_array();
}

// ==========================================================================

