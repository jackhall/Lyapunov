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

#include <vector>
#include <boost/python.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

namespace lyapunov {
	
	//needed for element-wise list adding
	void LengthError() { PyErr_SetString(PyExc_IndexError, "List lengths don't match."); }
	void RootFindError() { PyErr_SetString(PyExc_RuntimeError, "Multiple roots detected."); }
	void RevertError() { PyErr_SetString(PyExc_RuntimeError, "No state saved for necessary revert."); }

	struct Interval {
		double lower, upper;
		double length() const { return upper - lower; }
		double midpoint() const { return (upper + lower) / 2.0; }
	};

	struct vector_to_python_tuple {
		vector_to_python_tuple() {
			using namespace boost::python;
			to_python_converter<std::vector<double>, vector_to_python_tuple>();
		}

		static PyObject* convert(const std::vector<double>& x) {
			using namespace boost::python;
			list new_tuple; //is there a way around this list middle stage?
			for(auto i : x) new_tuple.append(i);
			return incref(tuple(new_tuple).ptr()); 
		}
	};
	struct vector_from_python_tuple {
		vector_from_python_tuple() {
			using namespace boost::python;
			converter::registry::push_back(&convertible, &construct, 
										   type_id<std::vector<double>>());
		}

		static void* convertible(PyObject* obj_ptr) {
			if (!PyTuple_Check(obj_ptr)) return 0;
			return obj_ptr;
		}

		static void construct(PyObject* obj_ptr,
			    boost::python::converter::rvalue_from_python_stage1_data* data) {
			using namespace boost::python;
			assert(PyTuple_Check(obj_ptr));
			unsigned int length = PyTuple_Size(obj_ptr);
			void* storage = ((converter::rvalue_from_python_storage< std::vector<double> >*)data)->storage.bytes;
			new (storage) std::vector<double>(length);
			for(unsigned int i=0; i<length; ++i)
				static_cast< std::vector<double>* >(storage)->at(i) 
					= PyFloat_AsDouble(PyTuple_GetItem(obj_ptr, i));
			data->convertible = storage;
		}
	};

	template<typename stepper_type>
	boost::python::list find_root(stepper_type& stepper, double min_step_size) {
		//implements a simple bisection rootfinder
		//make this an independent function, pass in system, min_step_size?
		//should events be a property of system? should they be a set?
		//should system be an iterator?
		namespace bp = boost::python;

		//Initialize interval
		Interval interval = {stepper.previous_time, 
							 bp::extract<double>(stepper.system.attr("time"))};
		//Revert (need to be able to step back across the boundary)
		if(!stepper.revert()) RevertError();

		//Loop over events
			//call and store event function values to vector
		bp::ssize_t num_events = bp::len(stepper.events);
		std::vector<double> starting_values(num_events), test_values(num_events);
		for(bp::ssize_t i=0; i<num_events; ++i) 
			starting_values[i] = bp::extract<double>(stepper.events[i]());

		while(interval.length() > min_step_size) {
			//step to midpoint of interval
			stepper.step(interval.length() / 2.0);
			//check for sign changes
			for(bp::ssize_t i=0; i<num_events; ++i) {
				if(starting_values[i]*bp::extract<double>(stepper.events[i]()) < 0) {
					stepper.revert();
					interval.upper = interval.midpoint();
					break;
				}
				if(i == num_events-1) interval.lower = interval.midpoint();
			}
		}

		//find out which events changed sign and return them
		stepper.step(interval.length());
		bp::list flagged;
		for(bp::ssize_t i=0; i<num_events; ++i) {
			if(starting_values[i]*bp::extract<double>(stepper.events[i]()) < 0) 
				flagged.append(stepper.events[i]);
		}
		stepper.revert();
		return flagged; //return interval.length()?
	}

	/*template<typename stepper_type>
	struct std_step {
		namespace bp = boost::python;
		typedef typename stepper_type::state_type state_type;
		typedef typename stepper_type::value_type num_type; //should be double

		struct sys_functor {
			bp::object system;

			sys_functor() = delete;
			sys_functor(bp::object sys) : system(sys) {}

			void operator()(const state_type& x, state_type& dx, const num_type t) {
				system.attr("time") = t;
				system.attr("state") = x;
				dx = system();
			}
		};

		static void step(stepper_type& self, bp::object system, num_type step_size) {
			//get current state and time
			state_type next_state = system.attr("state"); //will be altered in place
			num_type time = system.attr("time");

			//get the next state
			sys_functor sys(system);
			self.do_step(sys_functor, next_state, system.attr("time"), step_size);
			
			//make sure system object has completed the step
			system.attr("time") = time + step_size;
			system.attr("state") = next_state;
		}
	}; */

	template<typename stepper_type>
	class explicit_stepper_wrapper {
		typedef typename stepper_type::state_type state_type;
		typedef typename stepper_type::value_type num_type; //should be double

		boost::python::object system;
		bool revert_possible;
		stepper_type stepper;
		state_type saved_state, temporary; 
		num_type saved_time;

	public:
		explicit_stepper_wrapper() = delete;
		explicit_stepper_wrapper(boost::python::object sys) : system(sys), revert_possible(false) {}
		
		boost::python::object get_system() const { return system; }
		void set_system(boost::python::object new_system) { system = new_system; }

		//system function can't be passed to do_step ... because it's a method?
		void system_function(const state_type& x, state_type& dx, const num_type t) {
			system.attr("time") = t;
			system.attr("state") = x;
			dx = system();
		}
		void step(num_type step_size) {
			namespace bp = boost::python;
			saved_time = bp::extract<double>(system.attr("time"));
			temporary = bp::extract<state_type>(system.attr("state"));
			saved_state = temporary;
			revert_possible = true;
			//can I alter system.state in place? if not, use a temporary vector?
			stepper.do_step(system_function, temporary, saved_time, step_size);
			system.attr("time") = saved_time + step_size;
			system.attr("state") = temporary;
		}
		bool revert() {
			if(!revert_possible) return false;
			system.attr("time") = saved_time;
			system.attr("state") = saved_state;
			revert_possible = false;
		}
	};
}

//#include "cash_karp.h"

BOOST_PYTHON_MODULE(solvers) {
	using namespace boost::python;
	namespace ode = boost::numeric::odeint;
	using namespace lyapunov;

	vector_to_python_tuple vec2tup;
	vector_from_python_tuple tup2vec;
	//to_python_converter<std::vector<double>, vector_to_python_tuple>();

	//def("find_root", &find_root<object>);

	typedef std::vector<double> state_type;
	typedef ode::runge_kutta4<state_type> rk4;
	class_< explicit_stepper_wrapper<rk4> >("runge_kutta4")
		.def("step", &explicit_stepper_wrapper<rk4>::step)
		.def("revert", &explicit_stepper_wrapper<rk4>::revert)
		.add_property("system", &explicit_stepper_wrapper<rk4>::get_system, 
								&explicit_stepper_wrapper<rk4>::set_system);

	/*class_<Stepper>("Stepper", init<object>())
		.def("step", &Stepper::step)
		.def("euler_step", &Stepper::euler_step)
		.def("find_root", &Stepper::find_root)
		.def("revert", &Stepper::revert)
		.def("save", &Stepper::save)
		.def("__len__", &Stepper::get_num_states)
		.add_property("error", &Stepper::get_error)
		.add_property("step_size", &Stepper::get_step_size)
		.add_property("system", &Stepper::get_system, &Stepper::set_system); */
}

