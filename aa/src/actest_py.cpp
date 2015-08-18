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
#include "orbit.h"
#include "stackel_aa.h"
#include "lmn_orb.h"
#include "it_torus.h"
#include "cubature/cubature.h"
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(actest_py)
{
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    // Potentials

    class_<Potential_JS>("Potential_JS",no_init)
        .def("__call__",&Potential_JS::Phi_py)
        .def("energy",&Potential_JS::energy_py)
        .def("forces",&Potential_JS::Forces_py);

    class_<SphericalPotential, bases<Potential_JS> >("SphericalPotential",no_init);
    class_<IsochronePotential, bases<SphericalPotential> >("IsochronePotential",init<double,double>());
    class_<NFWSpherical, bases<SphericalPotential> >("NFWSpherical",init<double,double>());
    class_<HernquistSpherical, bases<SphericalPotential> >("HernquistSpherical",init<double,double>());
    class_<PowerLawSpherical, bases<SphericalPotential> >("PowerLawSpherical",init<double,double>());

    // class_<StackelProlate_PerfectEllipsoid, bases<Potential_JS>>("StackelProlate_PerfectEllipsoid",init<double,double>());
    // class_<StackelTriaxial, bases<Potential_JS>>("StackelTriaxial",init<double,double,double>());
    class_<Logarithmic, bases<Potential_JS>>("Logarithmic",init<double,double,double>());
    class_<PowerLaw, bases<Potential_JS>>("PowerLaw",init<double,double,double,double>());
    class_<Isochrone, bases<Potential_JS>>("Isochrone",init<double,double,double,double>());
    class_<HarmonicOscillator, bases<Potential_JS>>("HarmonicOscillator",init<double,double,double>());

    class_<GalPot, bases<Potential_JS> >("GalPot", init<std::string>());

    // class_<Actions_AxisymmetricStackel_Fudge>
    //     ("Actions_AxisymmetricStackel_Fudge",init<Potential_JS*,double>())
    //     .def("actions",&Actions_AxisymmetricStackel_Fudge::actions_py);

    // class_<IterativeTorusMachine>
    //     ("IterativeTorusMachine",
    //      init<Actions_AxisymmetricStackel_Fudge*,GalPot*,double,int,double>())
    //     .def("actions",&IterativeTorusMachine::actions_py)
    //     .def("actions_and_freqs",&IterativeTorusMachine::actionsandfreqs_py);

    // class_<Orbit>
    //     ("Orbit", init<Potential_JS*,double>())
    //     .def("integrate",&Orbit::integrate_py);

    // class_<test>("test",init<Potential_JS*>()).def("printer",&test::printer);
    // class_<Orbit>("Orbit", init<std::string>())
    // .def(init<double, double>())
    // .def("greet", &World::greet)
    // .def("set", &World::set)
;
}
