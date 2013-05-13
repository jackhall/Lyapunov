/*
	Lyapunov: a library for integrating nonlinear dynamical systems
	Copyright (C) 2013  John Wendell Hall

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
	The author may be reached at jackhall@utexas.edu.
*/

#include <boost/python.hpp>
#include "cash_karp.h"

BOOST_PYTHON_MODULE(lyapunov) {
	using namespace boost::python;

	class_<Stepper>("Stepper", init<object>())
		.def("step", &Stepper::step)
		.def("revert", &Stepper::revert)
		.add_property("state", &Stepper::get_state, &Stepper::set_state)
		.add_property("time", &Stepper::get_time, &Stepper::set_time)
		.add_property("error", &Stepper::get_error)
		.add_property("step_size", &Stepper::get_step_size)
		.add_property("system", &Stepper::get_system, &Stepper::set_system);
}
