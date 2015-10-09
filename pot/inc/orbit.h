// ============================================================================
/// \file inc/orbit.h
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
/// \brief Orbit integration classes
///
/// Orbit_base implements a base class for common functionality amongst orbit
/// integration class and there are two inherited classes
/// 1. Orbit -- integrates using GSL ode routines (default Dortmund-Prince 8)
/// 2. SymplecticOrbit -- integrates using 2nd order kick-drift-kick
///
//============================================================================

#ifndef ORBIT_H
#define ORBIT_H
// ============================================================================

#include "GSLInterface/GSLInterface.h"

// ============================================================================
class Orbit_base{
	protected:
		Potential_JS *Pot;
		std::vector<VecDoub> results_;  /*!< Stores orbit integration results*/

	public:
		Orbit_base(Potential_JS *P):Pot(P){}

		std::vector<VecDoub> const &results() const { return results_; }
		void plot(int i, int j, std::string name = "orbit");
		void output2file(std::string file);
		virtual VecDoub integrate(const VecDoub& x, double t_interval, double step_size, bool adaptive=false) = 0;
		double dynamic_time(const VecDoub& x);
};

/*! Orbit integrator using GSL ode routines */
class Orbit: public Orbit_base{
	/* class to integrate a phase-space point in a given Potential_JS */
	private:
		std::unique_ptr<ode> O;
	public:
	    //! Orbit constructor.
	    /*!
	        \param P -- potential
	        \param eps -- tolerance for orbit integration
			\param type -- type of ode integration routine from GSL
	    */
		Orbit(Potential_JS *P, double eps = 1e-10,std::string type="rkdp8");
		//! Integrate phase-space point.
	    /*!
	        \param x -- starting phase-space point x=(x,y,z,vx,vy,vz)
	        \param t_interval -- total integration time
	        \param step_size -- integration output step size
	        \param adaptive -- adaptive step-sizes

	        \return final phase-space point

	        intermediate results stored in results
	    */
		VecDoub integrate(const VecDoub& x, double t_interval, double step_size, bool adaptive=false);
		void SoS(int comp, VecDoub Init,std::string outfile);

};

class SymplecticOrbit: public Orbit_base{
	/* class to integrate a phase-space point in a given Potential_JS */
	private:
		double fac; /*!< govern accuracy */
	public:
		//! SymplecticOrbit constructor.
	    /*!
	        \param P -- potential
	        \param fac -- governs accuracy
	    */
		SymplecticOrbit(Potential_JS *P, double fac = 100.):Orbit_base(P),fac(fac){};
		//! Integrate phase-space point.
	    /*!
	        \param x -- starting phase-space point x=(x,y,z,vx,vy,vz)
	        \param t_interval -- total integration time
	        \param step_size -- integration output step size
	        \param adaptive -- adaptive step-sizes

	        \return final phase-space point

	        intermediate results stored in results
	    */
		VecDoub integrate(const VecDoub& x, double t_interval, double step_size, bool adaptive=false);
		void kdk2(VecDoub&x,VecDoub &a, double dt);

};
#endif
// ============================================================================
