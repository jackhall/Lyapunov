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

#include <iostream>
#include <vector>
#include <functional>
#include <boost/python.hpp>
#include <boost/numeric/odeint.hpp>
//#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
//#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
//#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>

namespace lyapunov {
	
	//needed for element-wise list adding
	void NotImplementedError() {
		PyErr_SetString(PyExc_NotImplementedError, "Feature not provided."); 
		boost::python::throw_error_already_set();
	}
	void LengthError() { 
		PyErr_SetString(PyExc_IndexError, "List lengths don't match."); 
		boost::python::throw_error_already_set();
	}
	void RuntimeError(const char* error_string) {
		PyErr_SetString(PyExc_RuntimeError, error_string);
		boost::python::throw_error_already_set();
	}
	void StopIteration() { 
		PyErr_SetString(PyExc_StopIteration, "Simulation finished."); 
		boost::python::throw_error_already_set();
	}

	struct Interval {
		double lower, upper;
		double length() const { return upper - lower; }
		double midpoint() const { return (upper + lower) / 2.0; }
	};

	struct vector_to_python_tuple {
		vector_to_python_tuple() {
			using namespace boost::python;
			to_python_converter<const std::vector<double>, vector_to_python_tuple>();
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

	boost::python::object pass_through(const boost::python::object& obj) { return obj; }
	
	//should step_across, next, and set_events be redefined for multistepper? probably
	//use bp::arg for keyword arguments
	class stepper_wrapper {
	public:
		typedef double num_type;
		typedef std::vector<num_type> state_type;

	protected: //makes derived code much easier to read
		enum stepper_state_type {NOTHING, LAST, BOUNDARY};

		boost::python::object system, steps, next_time_obj;
		state_type saved_state, event_signs;
		num_type saved_time; 
		stepper_state_type saved_information;
		std::function<void(const state_type&, state_type&, num_type)> system_function;
		bool events_occurred() const {
			namespace bp = boost::python;
			if(tracking_events) {
				bp::object event, iter = system.attr("events").attr("__iter__")();
				for(auto x : event_signs) {
					event = iter.attr("next")();
					if(x * bp::extract<num_type>( event() ) < 0) return true;
				}
			} 
			return false;
		}
		void update_signs() {
			namespace bp = boost::python;
			if(tracking_events) {
				event_signs.resize( bp::len(system.attr("events")) );
				bp::object event, iter = system.attr("events").attr("__iter__")();
				for(auto& x : event_signs) {
					//must not throw StopIteration early!
					event = iter.attr("next")(); 
					x = bp::extract<num_type>( event() );
				}
			}
		}

	public:
		num_type event_tolerance;
		bool tracking_events;

		stepper_wrapper() = delete;
		stepper_wrapper(boost::python::object sys, 
					    boost::python::object time)
			: system(sys), 
			  steps( time.attr("__iter__")() ), 
			  saved_state( boost::python::len(sys.attr("state")[1]) ),
			  event_signs(),
			  saved_time(0.0),
			  saved_information(NOTHING),
			  next_time_obj(), 
			  system_function([this](const state_type& x, state_type& dx, 
												 const num_type t) { 
				namespace bp = boost::python;
				system.attr("state") = bp::make_tuple(t, x);
				dx = bp::extract<state_type>(system()); } ),
			  event_tolerance(0.001),
			  tracking_events(PyObject_HasAttrString(sys.ptr(), "events")) {
			update_signs(); //passes through if !tracking_events
		}
		virtual void step(num_type next_time) = 0;
		virtual void step_back() = 0;
		virtual void reset() { saved_information = NOTHING; }
		void step_across() { 
			namespace bp = boost::python;
			if(saved_information != BOUNDARY) RuntimeError("No boundary to step across.");
			system.attr("state") = bp::make_tuple(saved_time, saved_state);
			saved_information = NOTHING;
		}
		boost::python::object get_system() const { return system; }
		void set_system(boost::python::object new_system) { system = new_system; }
		void use_times(boost::python::object time) {
			//stepper does not store the container, only the iterator
			steps = time.attr("__iter__")(); 
			next_time_obj = boost::python::object();
		}
		boost::python::object next() {
			namespace bp = boost::python;
			//get time for this next step
			if( next_time_obj.is_none() ) 
				next_time_obj = steps.attr("next")(); //may throw StopIteration
			num_type next_time = bp::extract<num_type>(next_time_obj);
			//take the step
			step(next_time);
			//check for event function sign changes, if any
			if( events_occurred() ) {
				auto flagged = find_root();
				//next_time_obj is not reset!
				return bp::make_tuple(system.attr("state")[0], flagged);
			} else {
				//reset next_time_obj to NoneType so the next call
				//will continue iterating through steps
				next_time_obj = bp::object();
				if(tracking_events) 
					return bp::make_tuple(system.attr("state")[0], bp::list());
				else return system.attr("state")[0];
			}
		}
		boost::python::object find_root() {
			//implements a simple bisection rootfinder
			namespace bp = boost::python;

			//initialize interval and step size goal
			Interval interval = {saved_time, 
								 bp::extract<num_type>(system.attr("state")[0]) };
			num_type tolerance = event_tolerance * interval.length();

			//Revert (need to be able to step back across the boundary)
			step_back();

			//main rootfinding loop	
			while(interval.length() > tolerance) {
				//step to midpoint of interval
				step(interval.midpoint());

				//check for sign changes
				if( events_occurred() ) {
					step_back();
					interval.upper = interval.midpoint();
				} else interval.lower = interval.midpoint();
			}

			//find out which events changed sign and return them
			step(interval.upper);
			num_type across_time = bp::extract<num_type>(system.attr("state")[0]);
			state_type across_state = bp::extract<state_type>(system.attr("state")[1]);
			bp::list flagged;
			bp::object event, iter = system.attr("events").attr("__iter__")();
			for(auto x : event_signs) {
				event = iter.attr("next")();
				if(x * bp::extract<num_type>(event()) < 0) 
					flagged.append(event);
			}
			step_back();

			//save state for step_across
			saved_time = across_time;
			saved_state = std::move(across_state);
			saved_information = BOUNDARY;
			return flagged; 
		}
		boost::python::str get_status() const {
			switch(saved_information) {
				case NOTHING:
					return "nothing";
				case LAST:
					return "last step";
				case BOUNDARY:
					return "boundary";
				default:
			}
		}
	};
	template<typename stepper_type>
	class explicit_stepper_wrapper : public stepper_wrapper {
	protected:
		stepper_type stepper;
		state_type temporary; 

	public:
		using stepper_wrapper::stepper_wrapper;
		
		void step(num_type next_time) {
			namespace bp = boost::python;
			bp::object state_tup = system.attr("state");
			saved_time = bp::extract<num_type>(state_tup[0]);
			saved_state = bp::extract<state_type>(state_tup[1]);
			saved_information = LAST;
			stepper.do_step(system_function, 
							saved_state, 
							saved_time, 
							temporary, 
							next_time - saved_time); //(sys, xin, tin, xout, h)
			system.attr("state") = bp::make_tuple(next_time, temporary);
		}
		void step_back() {
			namespace bp = boost::python;
			if(saved_information != LAST) RuntimeError("Last state not saved.");
			system.attr("state") = bp::make_tuple(saved_time, saved_state);
			saved_information = NOTHING;
		}
	};
	template<typename stepper_type, unsigned int stepper_order, unsigned int error_order>
	class error_stepper_wrapper : public explicit_stepper_wrapper<stepper_type> {
		typedef explicit_stepper_wrapper<stepper_type> base_type;
	public:
		typedef typename base_type::state_type state_type;
		typedef typename base_type::num_type num_type;
		num_type minimum_step_size, absolute_tolerance, relative_tolerance;
		int steps_taken;

	private:
		state_type error;
		num_type step_size;
		using base_type::stepper;
		using base_type::system;
		using base_type::system_function;
		using base_type::saved_state;
		using base_type::temporary;
		using base_type::saved_time;
		using base_type::saved_information;
		using base_type::LAST;

	public:
		error_stepper_wrapper() = delete;
		error_stepper_wrapper(boost::python::object sys,
							  boost::python::object time)
			: base_type(sys, time), 
			  minimum_step_size(),
			  absolute_tolerance(),
			  relative_tolerance(),
			  steps_taken(0),
			  error(saved_state.size()),
	   		  step_size() {}
		
		virtual void step(num_type next_time) {
			namespace bp = boost::python;
			bp::object state_tup = system.attr("state");
			saved_time = bp::extract<num_type>(state_tup[0]);
			saved_state = bp::extract<state_type>(state_tup[1]);
			saved_information = LAST;
			//note: a private do_step subroutine is probably a good idea
			//while current_time < next_time
				//if step_size < (next_time - current_time)
					//step with step_size
					//if error_index<0.5
						//calculate new step size
				//else
					//step with (next_time - current_time)
				//if error_index>1 
					//calculate new step size
			//set system state
			stepper.do_step(system_function, 
							saved_state, 
							saved_time, 
							temporary, 
							next_time - saved_time, 
							error); //(sys, xin, tin, xout, h, e)
			system.attr("state") = bp::make_tuple(next_time, temporary);
		}
		virtual void reset() {}

		}
	};
	template<typename stepper_wrapper_type>
	struct multistepper_wrapper : public stepper_wrapper_type {
		using stepper_wrapper_type::state_type;
		using stepper_wrapper_type::num_type;
		
		using stepper_wrapper_type::stepper_wrapper_type;
		void reset() { 
			stepper_wrapper_type::stepper.reset(); //use initialize instead?
			stepper_wrapper_type::saved_information = stepper_wrapper_type::NOTHING;
		}
	};
	template<typename stepper_type>
	class controlled_stepper_wrapper : public stepper_wrapper {
		typedef stepper_wrapper base_type;
		stepper_type stepper; //must be a Controlled Stepper
		state_type absolute_tolerance, relative_tolerance;

	public:
		void step(num_type next_time) {
			
		}
		bool step_back() {
			return true;
		}
	};
}

//macros for exposing wrapped steppers
//assumes 'using namespace boost::python'
//no semicolon afterwards
//TO DO: make steppers iterable (see std_iterator above for reference)
#define LYAPUNOV_EXPOSE_SIMPLE_STEPPER(STEPPER) { \
class_< explicit_stepper_wrapper< STEPPER > >(#STEPPER, init<object, object>()) \
	.def("step_back", &explicit_stepper_wrapper< STEPPER >::step_back) \
	.def("reset", &explicit_stepper_wrapper< STEPPER >::reset) \
	.def("__iter__", pass_through) \
	.def("next", &explicit_stepper_wrapper< STEPPER >::next) \
	.def("step_across", &explicit_stepper_wrapper< STEPPER >::step_across) \
	.def_readwrite("event_tolerance", &explicit_stepper_wrapper< STEPPER >::event_tolerance) \
	.def_readwrite("tracking_events", &explicit_stepper_wrapper< STEPPER >::tracking_events) \
	.def("use_times", &explicit_stepper_wrapper< STEPPER >::use_times) \
	.add_property("system", &explicit_stepper_wrapper< STEPPER >::get_system, &explicit_stepper_wrapper< STEPPER >::set_system); }

#define LYAPUNOV_EXPOSE_STEPPER(STEPPER, WRAPPER) { \
class_< WRAPPER< STEPPER > >(#STEPPER, init<object, object>()) \
	.def("step_back", &WRAPPER< STEPPER >::step_back) \
	.def("reset", &WRAPPER< STEPPER >::reset) \
	.def("__iter__", pass_through) \
	.def("next", &WRAPPER< STEPPER >::next) \
	.def("step_across", &WRAPPER< STEPPER >::step_across) \
	.def_readwrite("event_tolerance", &WRAPPER< STEPPER >::event_tolerance) \
	.def_readwrite("tracking_events", &WRAPPER< STEPPER >::tracking_events) \
	.def("use_times", &WRAPPER< STEPPER >::use_times) \
	.add_property("system", &WRAPPER< STEPPER >::get_system, &WRAPPER< STEPPER >::set_system) \
	.add_property("error", &WRAPPER< STEPPER >::get_error); }

//macros to help with variable order solvers
//no semicolon afterwards
#define LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(VARSTEPPER) { \
typedef ode::VARSTEPPER<1, state_type> VARSTEPPER##1; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER##1); \
typedef ode::VARSTEPPER<2, state_type> VARSTEPPER##2; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER##2); \
typedef ode::VARSTEPPER<3, state_type> VARSTEPPER##3; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER##3); \
typedef ode::VARSTEPPER<4, state_type> VARSTEPPER##4; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER##4); \
typedef ode::VARSTEPPER<5, state_type> VARSTEPPER##5; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER##5); \
typedef ode::VARSTEPPER<6, state_type> VARSTEPPER##6; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER##6); \
typedef ode::VARSTEPPER<7, state_type> VARSTEPPER##7; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER##7); \
typedef ode::VARSTEPPER<8, state_type> VARSTEPPER##8; \
LYAPUNOV_EXPOSE_MULTISTEPPER(VARSTEPPER##8); }

BOOST_PYTHON_MODULE(solvers) {
	using namespace boost::python;
	namespace ode = boost::numeric::odeint;
	using namespace lyapunov;

	vector_to_python_tuple vec2tup;
	vector_from_python_tuple tup2vec;
	//to_python_converter<std::vector<double>, vector_to_python_tuple>();

	typedef stepper_wrapper::num_type num_type;
	typedef stepper_wrapper::state_type state_type;

	//add saved_information observation (return a string)

	typedef ode::euler<state_type> euler;
	LYAPUNOV_EXPOSE_STEPPER(euler, explicit_stepper_wrapper)
	typedef ode::modified_midpoint<state_type> modified_midpoint;
	LYAPUNOV_EXPOSE_STEPPER(modified_midpoint, explicit_stepper_wrapper)
	typedef ode::runge_kutta4<state_type> runge_kutta4;
	LYAPUNOV_EXPOSE_STEPPER(runge_kutta4, explicit_stepper_wrapper)

	typedef ode::runge_kutta_cash_karp54<state_type> cash_karp;
	LYAPUNOV_EXPOSE_STEPPER(cash_karp, error_stepper_wrapper)
	typedef ode::runge_kutta_fehlberg78<state_type> fehlberg87; //order 7 error est.
	LYAPUNOV_EXPOSE_STEPPER(fehlberg87, error_stepper_wrapper)

	//turn error_stepper_wrapper into a variable_stepper_wrapper
	//remove error observation

	//typedef ode::runge_kutta_dopri5<state_type> dormand_prince;
	//typedef error_multistepper_wrapper<dormand_prince> dopri_wrap;
	//class_<dopri_wrap>("dormand_prince", 
	//				   init<object, object, optional<object, double> >())
	//	.def("step", &dopri_wrap::step) \
	//	.def("reset", &dopri_wrap::reset) \
	//	.def("step_back", &dopri_wrap::step_back) \
	//	.def("__iter__", pass_through) \
	//	.def("next", &dopri_wrap::next) \
	//	.def("step_across", &dopri_wrap::step_across) \
	//	.def("use_times", &dopri_wrap::use_times) \
	//	.def_readwrite("tolerance", &dopri_wrap::tolerance) \
	//	.add_property("events", &dopri_wrap::get_events, &dopri_wrap::set_events) \
	//	.add_property("system", &dopri_wrap::get_system, &dopri_wrap::set_system) \
	//	.add_property("error", &dopri_wrap::get_error); 

	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_bashforth)
	//ode::adams_moulton has a weird extra argument for do_step (a buffer?)
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_moulton)
	//ode::adams_bashforth_moulton lacks a reset method for some reason (a bug?)
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_bashforth_moulton)
	
	//think about Bulirsch-Stoer solver too!
}

