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
#include <functional>
#include <boost/python.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>

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

	class stepper_wrapper {
	protected: //makes derived code much easier to read
		typedef double num_type;
		typedef std::vector<num_type> state_type;
		boost::python::object system;
		std::function<void(const state_type&, state_type, num_type)> system_function;

	public:
		stepper_wrapper() = delete;
		//trouble creating system_function!
		explicit stepper_wrapper(boost::python::object sys) 
			: system(sys), system_function([this](const state_type& x, state_type& dx, 
												 const num_type t) { 
				namespace bp = boost::python;
				system.attr("time") = t;
				system.attr("state") = x;
				dx = bp::extract<state_type>(system()); } ) {}
		virtual void step(num_type step_size) = 0;
		virtual bool revert() = 0;
		boost::python::object get_system() const { return system; }
		void set_system(boost::python::object new_system) { system = new_system; }
	};

	boost::python::list find_root(stepper_wrapper& stepper, boost::python::list events, 
								  typename stepper_wrapper::num_type min_step_size) {
		//implements a simple bisection rootfinder
		//make this an independent function, pass in system, min_step_size?
		//should events be a property of system? should they be a set?
		//should system be an iterator?
		namespace bp = boost::python;
		typedef typename stepper_wrapper::num_type num_type;
		typedef typename stepper_wrapper::state_type state_type;

		//Revert (need to be able to step back across the boundary)
		num_type limit_time = bp::extract<num_type>(stepper.get_system().attr("time"));
		if(!stepper.revert()) RevertError();
		//Initialize interval
		Interval interval = {bp::extract<num_type>(stepper.get_system().attr("time")),
							 limit_time};

		//Loop over events
			//call and store event function values to vector
		bp::ssize_t num_events = bp::len(events);
		state_type starting_values(num_events), test_values(num_events);
		for(bp::ssize_t i=0; i<num_events; ++i) 
			starting_values[i] = bp::extract<num_type>(events[i]());

		while(interval.length() > min_step_size) {
			//step to midpoint of interval
			stepper.step(interval.length() / 2.0);
			//check for sign changes
			for(bp::ssize_t i=0; i<num_events; ++i) {
				if(starting_values[i]*bp::extract<num_type>(events[i]()) < 0) {
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
			if(starting_values[i]*bp::extract<num_type>(events[i]()) < 0) 
				flagged.append(events[i]);
		}
		stepper.revert();
		return flagged; //return interval.length()?
	}
	
	template<typename stepper_type>
	class explicit_stepper_wrapper : public stepper_wrapper {
		//typedef typename stepper_type::state_type state_type;
		//typedef typename stepper_type::value_type num_type; //should be double
	protected:
		bool revert_possible;
		stepper_type stepper;
		state_type saved_state, temporary; 
		num_type saved_time;

	public:
		explicit_stepper_wrapper() = delete;
		explicit explicit_stepper_wrapper(boost::python::object sys) 
			: stepper_wrapper(sys), revert_possible(false) {
			namespace bp = boost::python;
			auto num_states = bp::len(sys.attr("state"));
			saved_state.resize(num_states);
			temporary.resize(num_states);
		}
	
		void step(num_type step_size) {
			namespace bp = boost::python;
			saved_time = bp::extract<num_type>(system.attr("time"));
			saved_state = bp::extract<state_type>(system.attr("state"));
			revert_possible = true;
			//can I alter system.state in place? if not, use a temporary vector?
			stepper.do_step(system_function, 
							saved_state, 
							saved_time, 
							temporary, 
							step_size); //(sys, xin, tin, xout, h)
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

	template<typename stepper_type>
	class error_stepper_wrapper : public explicit_stepper_wrapper<stepper_type> {
		typedef explicit_stepper_wrapper<stepper_type> base_type;
		typedef typename base_type::state_type state_type;
		typedef typename base_type::num_type num_type;
		state_type error;

	public:
		error_stepper_wrapper() = delete;
		explicit error_stepper_wrapper(boost::python::object sys) 
			: base_type(sys), error(base_type::saved_state.size()) {}
		
		void step(num_type step_size) {
			namespace bp = boost::python;
			base_type::saved_time = bp::extract<num_type>(base_type::system.attr("time"));
			base_type::saved_state = bp::extract<state_type>(base_type::system.attr("state"));
			base_type::revert_possible = true;
			//can I alter system.state in place? if not, use a temporary vector?
			//there is trouble accessing system within this lambda...
			base_type::stepper.do_step(base_type::system_function, 
									   base_type::saved_state, 
									   base_type::saved_time, 
									   base_type::temporary, 
									   step_size, 
									   error); //(sys, xin, tin, xout, h, e)
			base_type::system.attr("time") = base_type::saved_time + step_size;
			base_type::system.attr("state") = base_type::temporary;
		}
		boost::python::tuple get_error() const { return error; } 
		//error may not have been set when get_error() is called... raise exception?
	};

	template<typename stepper_type>
	class multistepper_wrapper : public explicit_stepper_wrapper<stepper_type> {
		typedef explicit_stepper_wrapper<stepper_type> base_type;
		typedef typename base_type::state_type state_type;
		typedef typename base_type::num_type num_type;

	public:
		multistepper_wrapper() = delete;
		explicit multistepper_wrapper(boost::python::object sys)
			: base_type(sys) {}
		
		void reset() { 
			base_type::stepper.reset();
			base_type::revert_possible = false;
		}
	};
}

//macros for exposing wrapped steppers
//assumes 'using namespace boost::python'
//no semicolon afterwards
#define LYAPUNOV_EXPOSE_SIMPLE_STEPPER(STEPPER) { \
class_< explicit_stepper_wrapper< STEPPER > >(#STEPPER, init<object>()) \
	.def("step", &explicit_stepper_wrapper< STEPPER >::step) \
	.def("revert", &explicit_stepper_wrapper< STEPPER >::revert) \
	.add_property("system", &explicit_stepper_wrapper< STEPPER >::get_system, &explicit_stepper_wrapper< STEPPER >::set_system); }

#define LYAPUNOV_EXPOSE_ERROR_STEPPER(STEPPER) { \
class_< error_stepper_wrapper< STEPPER > >(#STEPPER, init<object>()) \
	.def("step", &error_stepper_wrapper< STEPPER >::step) \
	.def("revert", &error_stepper_wrapper< STEPPER >::revert) \
	.add_property("system", &error_stepper_wrapper< STEPPER >::get_system, &error_stepper_wrapper< STEPPER >::set_system) \
	.add_property("error", &error_stepper_wrapper< STEPPER >::get_error); }

#define LYAPUNOV_EXPOSE_MULTISTEPPER(STEPPER) { \
class_< multistepper_wrapper< STEPPER > >(#STEPPER, init<object>()) \
	.def("step", &multistepper_wrapper< STEPPER >::step) \
	.def("revert", &multistepper_wrapper< STEPPER >::revert) \
	.def("reset", &multistepper_wrapper< STEPPER >::reset) \
	.add_property("system", &multistepper_wrapper< STEPPER >::get_system, &multistepper_wrapper< STEPPER >::set_system); }

//macro to help with variable order solvers
//no semicolon afterwards
#define LYAPUNOV_EXPOSE_VARIABLE_ORDER_SOLVER(VARSTEPPER) { \
typedef ode::VARSTEPPER<1, state_type> VARSTEPPER1; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER1); \
typedef ode::VARSTEPPER<2, state_type> VARSTEPPER2; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER2); \
typedef ode::VARSTEPPER<3, state_type> VARSTEPPER3; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER3); \
typedef ode::VARSTEPPER<4, state_type> VARSTEPPER4; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER4); \
typedef ode::VARSTEPPER<5, state_type> VARSTEPPER5; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER5); \
typedef ode::VARSTEPPER<6, state_type> VARSTEPPER6; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER6); }

BOOST_PYTHON_MODULE(solvers) {
	using namespace boost::python;
	namespace ode = boost::numeric::odeint;
	using namespace lyapunov;

	vector_to_python_tuple vec2tup;
	vector_from_python_tuple tup2vec;
	//to_python_converter<std::vector<double>, vector_to_python_tuple>();

	def("find_root", &find_root);
	
	typedef stepper_wrapper::num_type num_type;
	typedef stepper_wrapper::state_type state_type;

	typedef ode::runge_kutta4<state_type> runge_kutta4;
	LYAPUNOV_EXPOSE_SIMPLE_STEPPER(runge_kutta4)

	typedef ode::runge_kutta_cash_karp54<state_type> cash_karp;
	LYAPUNOV_EXPOSE_ERROR_STEPPER(cash_karp)
	typedef ode::runge_kutta_fehlberg78<state_type> fehlberg78;
	//LYAPUNOV_EXPOSE_ERROR_STEPPER(fehlberg78)

	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_SOLVER(adams_bashforth)
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_SOLVER(adams_moulton)
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_SOLVER(adams_bashforth_moulton)

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

