#include <iostream>
#include <vector>
//#include <boost/python.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>

//this example of odeint should work, but apparently I need a full copy of boost 1.53.0

typedef std::vector<double> state_type;

void sys(const state_type& x, state_type& dxdt, const double t) {
	dxdt[0] = -x[0];
}

int main() {
	namespace ode = boost::numeric::odeint;
	//namespace bp = boost::python;

	ode::runge_kutta_cash_karp54<state_type> stepper;
	state_type error = {0.0};
	state_type x = {1.0};
	std::cout << "xo = " << x << std::endl; 
	stepper.do_step(sys, x, 0.0, 1.0, error); //(ode_fun, x, t, dt, xerr)
	std::cout << "x1 = " << x << "xerr = " << error << std::endl;

	return 0;
}

