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

	//pass_through so that __iter__ can return self
	//template<typename T, typename I>
	//struct std_iterator {
	//	static T next(I& obj) {
	//		T* result = obj.next();
	//		if(!result) {
	//			StopIteration();
	//			boost::python::throw_error_already_set();
	//		}
	//		return *result;
	//	}
	//	
	//	static I pass_through(const I& obj) { return obj; }
	//};

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
	template<typename stepper_type>
	void std_revert(stepper_type& stepper) {
		if(!stepper.revert()) RuntimeError("No saved state to revert to.");
	}
	class stepper_wrapper;
	boost::python::object find_root(stepper_wrapper& stepper, 
									boost::python::object events, 
									double tolerance);

	//should step_through, next, and set_events be redefined for multistepper? probably
	//use bp::arg for keyword arguments
	class stepper_wrapper {
	public:
		typedef double num_type;
		typedef std::vector<num_type> state_type;
		num_type tolerance;

	protected: //makes derived code much easier to read
		boost::python::object system, steps, events, next_time_obj;
		state_type event_function_values;
		bool tracking_events;
		std::function<void(const state_type&, state_type&, num_type)> system_function;
		bool events_occurred() const {
			namespace bp = boost::python;
			if(tracking_events) {
				bp::object event, iter = events.attr("__iter__")();
				for(auto x : event_function_values) {
					event = iter.attr("next")();
					if(x * bp::extract<num_type>( event() ) < 0) return true;
				}
			} 
			return false;
		}
		void update_signs() {
			namespace bp = boost::python;
			if(tracking_events) {
				event_function_values.resize( bp::len(events) );
				bp::object event, iter = events.attr("__iter__")();
				for(auto& x : event_function_values) {
					//must not throw StopIteration early!
					event = iter.attr("next")(); 
					x = bp::extract<num_type>( event() );
				}
			}
		}

	public:
		stepper_wrapper() = delete;
		stepper_wrapper(boost::python::object sys, 
					    boost::python::object time,
						boost::python::object eventlist=boost::python::object(),
						num_type min_step_size=0.00001) //nearest 10us
			: system(sys), 
			  steps( time.attr("__iter__")() ), 
			  events(),
			  event_function_values(),
			  next_time_obj(), 
			  tracking_events(false),
			  tolerance(min_step_size), 
			  system_function([this](const state_type& x, state_type& dx, 
												 const num_type t) { 
				namespace bp = boost::python;
				system.attr("state") = bp::make_tuple(t, x);
				dx = bp::extract<state_type>(system()); } ) {
			set_events(eventlist); //sets events, event_function_values
		}
		virtual void step(num_type next_time) = 0;
		virtual bool revert() = 0;
		void step_through() { 
			//assumes events haven't been reset
			namespace bp = boost::python;
			if( tracking_events && !next_time_obj.is_none() ) {
				//use tolerance*2 to be sure of passing through 
				//check event functions instead?
				int counter = 1;
				num_type time;
				do {
					time = bp::extract<num_type>(system.attr("state")[0]);
					step(time + tolerance);
					++counter;
				} while( !events_occurred() or counter < 10 );
				if(counter >= 10) { RuntimeError("Could not find boundary."); }
			} 
		}
		boost::python::object get_system() const { return system; }
		void set_system(boost::python::object new_system) { system = new_system; }
		void use_times(boost::python::object time) {
			steps = time.attr("__iter__")(); 
			next_time_obj = boost::python::object();
			if(tracking_events) update_signs();
		}
		boost::python::object get_events() const { return events; }
		void set_events(boost::python::object new_events) {
			events = new_events;
			tracking_events = !events.is_none();
			update_signs(); 
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
				//std::cout << "Rootfinding...";
				auto flagged = find_root(*this, events, tolerance);
				//note that next_time_obj is not reset!
				//std::cout << "Root found." << std::endl;
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
	};
	template<typename stepper_type>
	class explicit_stepper_wrapper : public stepper_wrapper {
	public:
		typedef typename stepper_type::state_type state_type;
		typedef typename stepper_type::value_type num_type; //should be double
	protected:
		bool revert_possible;
		stepper_type stepper;
		state_type saved_state, temporary; 
		num_type saved_time;

	public:
		explicit_stepper_wrapper() = delete;
		explicit_stepper_wrapper(boost::python::object sys,
								 boost::python::object time,
								 boost::python::object eventlist=boost::python::list(),
								 num_type min_step_size=0.00001) 
			: stepper_wrapper(sys, time, eventlist, min_step_size), 
			  revert_possible(false) {
			namespace bp = boost::python;
			auto num_states = bp::len(sys.attr("state")[1]);
			saved_state.resize(num_states);
			temporary.resize(num_states);
		}
		
		void step(num_type next_time) {
			namespace bp = boost::python;
			bp::object state_tup = system.attr("state");
			saved_time = bp::extract<num_type>(state_tup[0]);
			saved_state = bp::extract<state_type>(state_tup[1]);
			revert_possible = true;
			//can I alter system.state in place? if not, use a temporary vector?
			stepper.do_step(system_function, 
							saved_state, 
							saved_time, 
							temporary, 
							next_time - saved_time); //(sys, xin, tin, xout, h)
			system.attr("state") = bp::make_tuple(next_time, temporary);
		}
		bool revert() {
			namespace bp = boost::python;
			if(!revert_possible) return false;
			system.attr("state") = bp::make_tuple(saved_time, saved_state);
			revert_possible = false;
			return true;
		}
	};
	template<typename stepper_type>
	class error_stepper_wrapper : public explicit_stepper_wrapper<stepper_type> {
		typedef explicit_stepper_wrapper<stepper_type> base_type;
	public:
		typedef typename base_type::state_type state_type;
		typedef typename base_type::num_type num_type;
	private:
		state_type error;
	public:
		error_stepper_wrapper() = delete;
		error_stepper_wrapper(boost::python::object sys,
							  boost::python::object time,
							  boost::python::object eventlist=boost::python::object(),
							  num_type min_step_size=0.00001) 
			: base_type(sys, time, eventlist, min_step_size), 
			  error(base_type::saved_state.size()) {}
		
		void step(num_type next_time) {
			namespace bp = boost::python;
			bp::object state_tup= base_type::system.attr("state");
			base_type::saved_time = bp::extract<num_type>(state_tup[0]);
			base_type::saved_state = bp::extract<state_type>(state_tup[1]);
			base_type::revert_possible = true;
			//can I alter system.state in place? if not, use a temporary vector?
			base_type::stepper.do_step(base_type::system_function, 
									   base_type::saved_state, 
									   base_type::saved_time, 
									   base_type::temporary, 
									   next_time - base_type::saved_time, 
									   error); //(sys, xin, tin, xout, h, e)
			base_type::system.attr("state") = bp::make_tuple(next_time, base_type::temporary);
		}
		boost::python::object get_error() { 
			if(base_type::revert_possible) 
				return boost::python::tuple(error); 
			else RuntimeError("No error estimate has been calcuated.");
		} 
	};
	template<typename stepper_type>
	class multistepper_wrapper : public explicit_stepper_wrapper<stepper_type> {
		typedef explicit_stepper_wrapper<stepper_type> base_type;
	public:
		typedef typename base_type::state_type state_type;
		typedef typename base_type::num_type num_type;
		
		using base_type::explicit_stepper_wrapper;
		void reset() { 
			base_type::stepper.reset();
			base_type::revert_possible = false;
		}
	};
	template<typename stepper_type>
	class error_multistepper_wrapper : public error_stepper_wrapper<stepper_type> {
		typedef error_stepper_wrapper<stepper_type> base_type;
	public:
		typedef typename base_type::state_type state_type;
		typedef typename base_type::num_type num_type;

		using base_type::error_stepper_wrapper;
		void reset() { 
			base_type::stepper.reset();
			base_type::revert_possible = false;
		}
	};

	boost::python::object find_root(stepper_wrapper& stepper, 
									boost::python::object events, 
								  	double tolerance=0.00001) {
		//implements a simple bisection rootfinder
		//make this an independent function, pass in system, min_step_size?
		//should events be a property of system? should they be a set?
		//should system be an iterator?
		namespace bp = boost::python;
		typedef typename stepper_wrapper::num_type num_type;
		typedef typename stepper_wrapper::state_type state_type;

		//Loop over events
			//call and store event function values to vector
		state_type starting_values( bp::len(events) );
		bp::object event, iter = events.attr("__iter__")();
		for(auto& x : starting_values) {
			event = iter.attr("next")();
			x = bp::extract<num_type>(event());
		}

		//Revert (need to be able to step back across the boundary)
		num_type limit_time = bp::extract<num_type>(stepper.get_system().attr("state")[0]);
		std_revert(stepper); //raises a python exception if revert fails
		//Initialize interval
		Interval interval = {bp::extract<num_type>(stepper.get_system().attr("state")[0]),
							 limit_time};

		//Loop over events again
			//check for changed signs - if none, raise an exception
			//call and store event function values to vector
		iter = events.attr("__iter__")();
		num_type new_value;
		bool event_occurred = false;
		for(auto& x : starting_values) {
			event = iter.attr("next")();
			new_value = bp::extract<num_type>(event());
			if(x*new_value < 0) event_occurred = true;
			x = bp::extract<num_type>(event());
		}
		if(!event_occurred) RuntimeError("No boundary was crossed.");

		bool boundary_crossed = false;
		while(interval.length() > tolerance) {
			//step to midpoint of interval
			stepper.step(interval.midpoint());
			//check for sign changes
			iter = events.attr("__iter__")();
			for(auto x : starting_values) {
				event = iter.attr("next")();
				boundary_crossed = x * bp::extract<num_type>(event()) < 0;
				if(boundary_crossed) break;
			}
			if(boundary_crossed) {
				stepper.revert();
				interval.upper = interval.midpoint();
				boundary_crossed = false;
			} else interval.lower = interval.midpoint();
		}

		//find out which events changed sign and return them
		stepper.step(interval.upper);
		bp::list flagged;
		iter = events.attr("__iter__")();
		for(auto x : starting_values) {
			event = iter.attr("next")();
			if(x * bp::extract<num_type>(event()) < 0) 
				flagged.append(event);
		}
		stepper.revert();
		return flagged; //return interval.length()?
	}
}

//macros for exposing wrapped steppers
//assumes 'using namespace boost::python'
//no semicolon afterwards
//TO DO: make steppers iterable (see std_iterator above for reference)
#define LYAPUNOV_EXPOSE_SIMPLE_STEPPER(STEPPER) { \
class_< explicit_stepper_wrapper< STEPPER > >(#STEPPER, init<object, object, optional<object, double> >()) \
	.def("step", &explicit_stepper_wrapper< STEPPER >::step) \
	.def("revert", &explicit_stepper_wrapper< STEPPER >::revert) \
	.def("__iter__", pass_through) \
	.def("next", &explicit_stepper_wrapper< STEPPER >::next) \
	.def("step_through", &explicit_stepper_wrapper< STEPPER >::step_through) \
	.def_readwrite("tolerance", &explicit_stepper_wrapper< STEPPER >::tolerance) \
	.def("use_times", &explicit_stepper_wrapper< STEPPER >::use_times) \
	.add_property("events", &explicit_stepper_wrapper< STEPPER >::get_events, &explicit_stepper_wrapper< STEPPER >::set_events) \
	.add_property("system", &explicit_stepper_wrapper< STEPPER >::get_system, &explicit_stepper_wrapper< STEPPER >::set_system); }

#define LYAPUNOV_EXPOSE_ERROR_STEPPER(STEPPER) { \
class_< error_stepper_wrapper< STEPPER > >(#STEPPER, init<object, object, optional<object, double> >()) \
	.def("step", &error_stepper_wrapper< STEPPER >::step) \
	.def("revert", &error_stepper_wrapper< STEPPER >::revert) \
	.def("__iter__", pass_through) \
	.def("next", &error_stepper_wrapper< STEPPER >::next) \
	.def("step_through", &error_stepper_wrapper< STEPPER >::step_through) \
	.def("use_times", &error_stepper_wrapper< STEPPER >::use_times) \
	.def_readwrite("tolerance", &error_stepper_wrapper< STEPPER >::tolerance) \
	.add_property("events", &error_stepper_wrapper< STEPPER >::get_events, &error_stepper_wrapper< STEPPER >::set_events) \
	.add_property("system", &error_stepper_wrapper< STEPPER >::get_system, &error_stepper_wrapper< STEPPER >::set_system) \
	.add_property("error", &error_stepper_wrapper< STEPPER >::get_error); }

#define LYAPUNOV_EXPOSE_MULTISTEPPER(STEPPER) { \
class_< multistepper_wrapper< STEPPER > >(#STEPPER, init<object, object, optional<object, double> >()) \
	.def("step", &multistepper_wrapper< STEPPER >::step) \
	.def("revert", &multistepper_wrapper< STEPPER >::revert) \
	.def("reset", &multistepper_wrapper< STEPPER >::reset) \
	.def("__iter__", pass_through) \
	.def("next", &multistepper_wrapper< STEPPER >::next) \
	.def("step_through", &multistepper_wrapper< STEPPER >::step_through) \
	.def("use_times", &multistepper_wrapper< STEPPER >::use_times) \
	.def_readwrite("tolerance", &explicit_stepper_wrapper< STEPPER >::tolerance) \
	.add_property("events", &multistepper_wrapper< STEPPER >::get_events, &multistepper_wrapper< STEPPER >::set_events) \
	.add_property("system", &multistepper_wrapper< STEPPER >::get_system, &multistepper_wrapper< STEPPER >::set_system); }

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

	def("find_root", &find_root);
	
	typedef stepper_wrapper::num_type num_type;
	typedef stepper_wrapper::state_type state_type;

	typedef ode::euler<state_type> euler;
	LYAPUNOV_EXPOSE_SIMPLE_STEPPER(euler)
	typedef ode::modified_midpoint<state_type> modified_midpoint;
	LYAPUNOV_EXPOSE_SIMPLE_STEPPER(modified_midpoint)
	typedef ode::runge_kutta4<state_type> runge_kutta4;
	LYAPUNOV_EXPOSE_SIMPLE_STEPPER(runge_kutta4)

	typedef ode::runge_kutta_cash_karp54<state_type> cash_karp;
	LYAPUNOV_EXPOSE_ERROR_STEPPER(cash_karp)
	typedef ode::runge_kutta_fehlberg78<state_type> fehlberg87; //order 7 error est.
	LYAPUNOV_EXPOSE_ERROR_STEPPER(fehlberg87)

	typedef ode::runge_kutta_dopri5<state_type> dormand_prince;
	typedef error_multistepper_wrapper<dormand_prince> dopri_wrap;
	class_<dopri_wrap>("dormand_prince", 
					   init<object, object, optional<object, double> >())
		.def("step", &dopri_wrap::step) \
		.def("reset", &dopri_wrap::reset) \
		.def("revert", &dopri_wrap::revert) \
		.def("__iter__", pass_through) \
		.def("next", &dopri_wrap::next) \
		.def("step_through", &dopri_wrap::step_through) \
		.def("use_times", &dopri_wrap::use_times) \
		.def_readwrite("tolerance", &dopri_wrap::tolerance) \
		.add_property("events", &dopri_wrap::get_events, &dopri_wrap::set_events) \
		.add_property("system", &dopri_wrap::get_system, &dopri_wrap::set_system) \
		.add_property("error", &dopri_wrap::get_error); 

	LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_bashforth)
	//ode::adams_moulton has a weird extra argument for do_step (a buffer?)
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_moulton)
	//ode::adams_bashforth_moulton lacks a reset method for some reason (a bug?)
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_bashforth_moulton)
	
	//think about Bulirsch-Stoer solver too!
}

