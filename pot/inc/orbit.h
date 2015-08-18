#ifndef ORBIT_H
#define ORBIT_H
// ============================================================================
#include "GSLInterface/GSLInterface.h"

// ============================================================================
class Orbit_base{
	protected:
		Potential_JS *Pot;
		std::vector<VecDoub> results_;

	public:
		Orbit_base(Potential_JS *P):Pot(P){}

		std::vector<VecDoub> const &results() const { return results_; }
		void plot(int i, int j, std::string name = "orbit");
		void output2file(std::string file);
		virtual VecDoub integrate(const VecDoub& x, double t_interval, double step_size, bool adaptive=false) = 0;
		double dynamic_time(const VecDoub& x);
};

class Orbit: public Orbit_base{
	/* class to integrate a phase-space point in a given Potential_JS */
	private:
		std::unique_ptr<ode> O;
	public:
		Orbit(Potential_JS *P, double eps = 1e-10,std::string type="rkdp8");
		VecDoub integrate(const VecDoub& x, double t_interval, double step_size, bool adaptive=false);
		void SoS(int comp, VecDoub Init,std::string outfile);

};

class SymplecticOrbit: public Orbit_base{
	/* class to integrate a phase-space point in a given Potential_JS */
	private:
		double fac;
	public:
		SymplecticOrbit(Potential_JS *P, double fac = 100.):Orbit_base(P),fac(fac){};
		VecDoub integrate(const VecDoub& x, double t_interval, double step_size, bool adaptive=false);
		void kdk2(VecDoub&x,VecDoub &a, double dt);

};
#endif
// ============================================================================
