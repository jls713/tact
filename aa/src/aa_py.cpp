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
#include "falPot.h"
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
#include "it_torus.h"
#include "stackel_fit.h"
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/list.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <numpy/arrayobject.h>
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

  // Potentials

  class_<Potential_JS, boost::noncopyable,
         boost::shared_ptr<Potential_JS_PythonCallback>>("Potential_JS",
                                                         init<>())
      .def("__call__", &Potential_JS::Phi)
      .def("energy", &Potential_JS::H)
      .def("Forces", &Potential_JS::Forces)
      .def("torb", &Potential_JS::torb);

  class_<Density, boost::noncopyable>("Density",init<double>());

  class_<TriaxialDensity, bases<Density>>("TriaxialDensity",init<VecDoub>());
  class_<Density_NFW, bases<TriaxialDensity>>("Density_NFW",init<double, double,VecDoub,double>());

  class_<TriaxialPotential, bases<Potential_JS>>("TriaxialPotential",init<TriaxialDensity*,double>());

  class_<SphericalPotential, bases<Potential_JS>>("SphericalPotential",
                                                  no_init);
  class_<IsochronePotential, bases<SphericalPotential>>("IsochronePotential",
                                                        init<double, double>());
  class_<NFWSpherical, bases<SphericalPotential>>("NFWSpherical",
                                                  init<double, double>());
  class_<HernquistSpherical, bases<SphericalPotential>>("HernquistSpherical",
                                                        init<double, double>());
  class_<PowerLawSpherical, bases<SphericalPotential>>("PowerLawSpherical",
                                                       init<double, double>());
  class_<StackelProlate_PerfectEllipsoid,
	bases<Potential_JS>>("StackelProlate_PerfectEllipsoid",init<double,double>());
  class_<StackelTriaxial,
	bases<Potential_JS>>("StackelTriaxial",init<double,double,double>());
  class_<Logarithmic, bases<Potential_JS>>("Logarithmic",
                                           init<double, double, double>());
  class_<PowerLaw, bases<Potential_JS>>("PowerLaw",
                                        init<double, double, double, double>());
  class_<Isochrone, bases<Potential_JS>>(
      "Isochrone", init<double, double, double, double>());
  class_<HarmonicOscillator, bases<Potential_JS>>(
      "HarmonicOscillator", init<double, double, double>());

  class_<GalPot, bases<Potential_JS>>("GalPot", init<std::string>());

  class_<BowdenNFW, bases<Potential_JS>>("BowdenNFW",init<double,double,double,double,double,double>());

  // Action Finders

  class_<Action_Finder>("Action_Finder",no_init)
	.def("actions",&Action_Finder::actions_nopars)
  .def("angles",&Action_Finder::angles_nopars);

  class_<Actions_Spherical, bases<Action_Finder> >(
	"Actions_Spherical", init<SphericalPotential*>());

  class_<Actions_Isochrone, bases<Action_Finder, Isochrone> >(
  "Actions_Isochrone", init<double,double>());

  class_<Actions_PolarAdiabaticApproximation, bases<Action_Finder> >(
	"Actions_PolarAdiabaticApproximation", init<Potential_JS*, std::string,bool,bool>());

  class_<Actions_SpheroidalAdiabaticApproximation, bases<Action_Finder> >(
	"Actions_SpheroidalAdiabaticApproximation", init<Potential_JS*, std::string,bool,bool,double>());

  class_<Actions_AxisymmetricStackel, bases<Action_Finder> >(
	"Actions_AxisymmetricStackel", init<StackelProlate_PerfectEllipsoid*>());

  class_<Actions_AxisymmetricStackel_Fudge, bases<Action_Finder> >(
      "Actions_AxisymmetricStackel_Fudge", init<Potential_JS*, double>());

  class_<Actions_TriaxialStackel, bases<Action_Finder> >(
	"Actions_TriaxialStackel", init<StackelTriaxial*>());

  class_<Actions_TriaxialStackel_Fudge, bases<Action_Finder> >(
	"Actions_TriaxialStackel_Fudge", init<Potential_JS*, double, double>());

  class_<Actions_Genfunc, bases<Action_Finder> >(
	"Actions_Genfunc", init<Potential_JS*,std::string>())
  .def("full_actions",&Actions_Genfunc::full_actions)
  .def("full_angles",&Actions_Genfunc::full_angles);

  class_<Actions_Genfunc_Average, bases<Action_Finder> >(
  "Actions_Genfunc_Average", init<Potential_JS*,std::string>());

  class_<Actions_StackelFit, bases<Action_Finder> >(
  "Actions_StackelFit", init<Potential_JS*>());


  class_<IterativeTorusMachine, bases<Action_Finder> >(
      "IterativeTorusMachine",
      init<Actions_AxisymmetricStackel_Fudge*, GalPot*, double, int, double>()).def("set_maxit",&IterativeTorusMachine::set_maxit);

  class_<uv_orb, bases<Action_Finder> >(
      "Actions_AxisymmetricStackel_Fudge_DeltaGuess",
      init<Potential_JS*, double, double, int, int, std::string>());

  class_<lmn_orb, bases<Action_Finder> >(
      "Actions_TriaxialStackel_Fudge_DeltaGuess",
      init<Potential_JS*, double, double, int, bool,bool, std::string>());


  class_<Orbit, boost::noncopyable>("Orbit", init<Potential_JS*, double>())
      .def("integrate", &Orbit::integrate);


  def("GalacticToPolar",conv::GalacticToPolar);
  def("GalacticToCartesian",conv::GalacticToCartesian);
  def("GalacticToEquatorial",conv::GalacticToEquatorial);
  def("CartesianToPolar",conv::CartesianToPolar);
  def("CartesianToGalactic",conv::CartesianToGalactic);
  def("PolarToGalactic",conv::PolarToGalactic);
  def("EquatorialToGalactic",conv::EquatorialToGalactic);



  to_python_converter<VecDoub, vector_to_ndarray<double>>();
  vector_from_ndarray<double>();
  import_array();
}

// ==========================================================================

