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

#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/numeric/odeint.hpp>
//#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
//#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
//#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>

namespace lyapunov {

    //convenience functions to call when I want to throw a particular kind of error
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
	void ValueError(const char* error_string) {
		PyErr_SetString(PyExc_ValueError, error_string);
		boost::python::throw_error_already_set();
	}
	void StopIteration() { 
		PyErr_SetString(PyExc_StopIteration, "Simulation finished."); 
		boost::python::throw_error_already_set();
	}

	struct Interval {
        /* A utility to handle the contracting interval of a bisection rootfinder.
         */
		double lower, upper;
		double length() const { return upper - lower; }
		double midpoint() const { return (upper + lower) / 2.0; }
	};

    //random utilites
	struct vector_to_python_tuple {
		vector_to_python_tuple() {
			using namespace boost::python;
			to_python_converter<const std::vector<double>, vector_to_python_tuple>();
		}

		static PyObject* convert(const std::vector<double>& x) {
            auto new_tuple = PyTuple_New(x.size());
            for(int i=x.size()-1; i>=0, --i) {
                PyTuple_SetItem(new_tuple, i, PyFloat_FromDouble(x[i]));
            }
            return boost::python::incref(new_tuple);
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
    double minimum(double x, double y) { 
        if(x>y) return y;
        else return x;
    }
    double maximum(double x, double y) {
        if(x<y) return y;
        else return x;
    }

    //should step_across, next, and set_events be defined for multistepper? probably
	//use bp::arg for keyword arguments
	class stepper_wrapper {
    /* This class represents both a wrapper for boost::numeric::odeint steppers and a 
     * state machine for a proper simulation with rootfinding. It appears to python as an
     * iterator.
     */
	public:
		typedef double num_type;
		typedef std::vector<num_type> state_type;

	protected: //makes derived code much easier to read
        // NOTHING: no saved point
        // LAST: the previous point is saved
        // BOUNDARY: the point on the other side of an event boundary is saved
		enum stepper_state_type {NOTHING, LAST, BOUNDARY};

		boost::python::object system, steps, next_time_obj;
		state_type temporary, saved_state, event_signs;
		num_type saved_time; 
		stepper_state_type saved_information; // which point is stored: previous or boundary
		std::function<void(const state_type&, state_type&, num_type)> system_function;
        bool signs_verified; // true if all event_signs are nonzero
		bool events_occurred() const {
            /* Check to see whether any event functions crossed zero since the last step.
             */
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
            /* Record the value of each event function at the current state. The only
             * important aspect of these numbers is their sign. If any of the values 
             * are exactly zero, try to record a value at the next step.
             */
			namespace bp = boost::python;
            event_signs.clear();
            sign_check = 1.0; //any nonzero value will do
            if(PyObject_HasAttrString(system.ptr(), "events")) {
                for(bp::stl_input_iterator<bp::object> iter( system.attr("events") ), end;
                      iter != end ; ++iter) {
                    event_signs.push_back(bp::extract<num_type>( (*iter)() ));
                    sign_check *= event_signs.back();
                }
            }
            if(sign_check != 0) signs_verified = true;
            if(event_signs.size() == 0) tracking_events = false;

		}
        void verify_signs() {
            /* Ensure that all event signs are nonzero;
             * that is: make sure that event functions equal to
             * zero at the initial conditions still get detected.
             */
            namespace bp = boost::python;
            zero_sign_count = 0;
            for(int i=event_signs.size()-1; i>=0; --i) {
                if(event_signs[i] == 0.0) {
                    event_signs[i] = bp:extract<num_type>(system.attr("events")[i]);
                    if(event_signs[i] == 0) ++zero_sign_count;
                }
            }
            if(zero_sign_count == 0) signs_verified = true;
        }
		void save_last() {
            /* Save the current point in anticipation of stepping to the next.
             */
			namespace bp = boost::python;
			bp::object state_tup = system.attr("state");
			saved_time = bp::extract<num_type>(state_tup[0]);
			saved_state = bp::extract<state_type>(state_tup[1]);
			saved_information = LAST;
		}
		
	public:
		num_type time_tolerance; //smallest unit of time the solver/rootfinder can see
		bool tracking_events; // true if the event list is nonempty

		stepper_wrapper() = delete;
		stepper_wrapper(boost::python::object sys, 
					    boost::python::object time)
			: system(sys), 
			  steps(), //actually set in use_times
			  next_time_obj(), //actually set in use_times
			  temporary( boost::python::len(sys.attr("state")[1]) ),
			  saved_state( temporary.size() ),
			  event_signs(), //actually set in update signs
			  saved_time(0.0),
			  saved_information(NOTHING),
			  system_function([this](const state_type& x, state_type& dx, const num_type t) { 
                  namespace bp = boost::python;
                  system.attr("state") = bp::make_tuple(t, x);
                  dx = bp::extract<state_type>(system()); } ),
			  time_tolerance(0.00001), //10 microseconds ... make this relative?
			  tracking_events(true), //actually set in update_signs
              signs_verified(false) { //actually set in update_signs (and verify_signs)
			use_times(time);
			update_signs(); 
		}
		virtual void step(num_type next_time) = 0;
		virtual void step_back() {
            /* Return the stepper to the previous point, if saved. Note that there is
             * no way to move back two steps in a row, since only one point is saved 
             * at a time.
             */
			namespace bp = boost::python;
			if(saved_information != LAST) RuntimeError("Last state not saved.");
			system.attr("state") = bp::make_tuple(saved_time, saved_state);
			reset();
		}
        void reset_with_events() {
            /* Make sure no point is saved by the solver, and recompute the sign of
             * each event function.
             */
            reset();
            update_signs();
        }
		virtual void reset() { 
            saved_information = NOTHING; 
        }
		void step_across(boost::python::object new_state = boost::python::object()) { 
            /* Move the system across the event boundary identified by the rootfinder.
             * This does not invoke the solver, since the point across the boundary is 
             * already known. After moving, recompute event function signs.
             */
			namespace bp = boost::python;
			if(saved_information != BOUNDARY) RuntimeError("No boundary to step across.");
			bp::object state_tup;
			if(new_state.is_none()) state_tup = bp::make_tuple(saved_time, saved_state);
			else state_tup = bp::make_tuple(saved_time, new_state);
			system.attr("state") = state_tup;
			reset_with_events(); 
		}
		num_type get_step_size() const { 
            /* Compute the time difference between this point and the saved previous point.
             * Doesn't work without a saved previous point.
             */
			namespace bp = boost::python;
			if(saved_information != LAST) RuntimeError("Last state not saved.");
			return bp::extract<num_type>(system.attr("state")[0]) - saved_time;
		}
		num_type get_time_tolerance() const { return time_tolerance; }
		void set_time_tolerance(num_type new_tolerance) {
			if(new_tolerance >= 0) time_tolerance = new_tolerance;
			else ValueError("Positive values only.");
		}
		boost::python::object get_system() const { return system; }
		void set_system(boost::python::object new_system) { system = new_system; }
		void use_times(boost::python::object time) {
            /* Set the iterable which generates time steps.
             */
			steps = time.attr("__iter__")(); 
			next_time_obj = boost::python::object();
		}
		boost::python::object next() {
            /* This method makes the wrapper an iterable in python. It numerically 
             * estimates the state at the next time step and checks to see if any
             * event functions changed sign. If so, it calls the rootfinder to move
             * the system right up to the event boundary without crossing. 
             */
			namespace bp = boost::python;
            //if an event occurred, step through it
            if(saved_information == BOUNDARY) step_across();
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
				return bp::make_tuple(system.attr("state")[0], bp::list());
			}
            if(not signs_verified) verify_signs();  //check that signs are nonzero
		}
		boost::python::object find_root() {
			/* Implements a simple bisection rootfinder.
             * Moves the system to within time_tolerance of the event boundary
             * and saves a point just across the boundary.
             *
             * If more than one event is detected, the rootfinder will stop at the
             * first one it finds. If the event boundaries are within time_tolerance
             * of each other, both get detected.
             */
			namespace bp = boost::python;

			//initialize interval and step size goal
			Interval interval = {saved_time, 
								 bp::extract<num_type>(system.attr("state")[0]) };

			//Revert (need to be able to step back across the boundary)
			step_back();

			//main rootfinding loop	
			while(interval.length() > time_tolerance) {
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
            /* Returns what point the stepper is currently storing.
             */
			switch(saved_information) {
				case NOTHING:
					return "nothing";
				case LAST:
					return "last step";
				case BOUNDARY:
					return "boundary";
				default:
					return "status undefined!";
			}
		}
	};
	class variable_stepper_wrapper : public stepper_wrapper {
    /* As a wrapper of variable step size solvers, this class extends stepper_wrapper
     * to (optionally) independently choose its own step size, halting at a given 
     * final time.
     */
	protected:
		typedef stepper_wrapper base_type;
		num_type absolute_tolerance, relative_tolerance;
		state_type temporary, current_state, error;
		num_type step_size, final_time;
        bool final_time_reached;

	public:
		variable_stepper_wrapper() = delete;
		variable_stepper_wrapper(boost::python::object sys,
								 boost::python::object time) 
			: base_type(sys, boost::python::list()), 
			  absolute_tolerance(0.000001),
			  relative_tolerance(0.001),
			  temporary(saved_state.size()),
			  current_state(saved_state.size()),
			  error(saved_state.size()),
	   		  step_size(-1),  // flagged to recompute
	   		  final_time(0.0),  // actually set in use_times
              final_time_reached(false) {
			use_times(time);
		}
		virtual void step_back() {
            /* See step_back docstring for stepper_wrapper.
             */
            final_time_reached = false;
            base_type::step_back();
        }
        num_type get_step_size() const { return step_size; }
		num_type get_relative_tolerance() const { return relative_tolerance; }
		void set_relative_tolerance(num_type new_tolerance) {
			if(new_tolerance >= 0) relative_tolerance = new_tolerance;
			else ValueError("Positive values only.");
		}
		num_type get_absolute_tolerance() const { return absolute_tolerance; }
		void set_absolute_tolerance(num_type new_tolerance) {
			if(new_tolerance >= 0) absolute_tolerance = new_tolerance;
			else ValueError("Positive values only.");
		}
		void use_times(boost::python::object time) {
            /*
             */
			namespace bp = boost::python;
            PyObject* time_ptr = time.ptr();
			if(!PySequence_Check(time_ptr) && PyNumber_Check(time_ptr)) {
				steps = boost::python::object();
				final_time = bp::extract<num_type>(time);
                if(step_size < 0) {
				    step_size = 0.001*(final_time - 
							bp::extract<num_type>(system.attr("state")[0]));
                }
			} else base_type::use_times(time);
            final_time_reached = false;
		}
		virtual void try_step(num_type current_time, num_type current_step_size) = 0;
        virtual void try_free_step(num_type current_step_size) = 0;
        double compute_error_index() {
            /*
             */
            double error_index = 0.0;
            for(int i=error.size()-1; i>=0; --i) {
				error_index = maximum(error_index, std::abs(error[i]) /
							  (absolute_tolerance + relative_tolerance*std::abs(temporary[i])));
			}
            return error_index;
        }
		virtual void increase_step(num_type error_index) = 0; 		
        virtual bool decrease_step(num_type error_index, num_type original_step_size) = 0;
		bool adjust_step(num_type error_index, num_type original_step_size) {
            /*
             */
			increase_step(error_index);
			return decrease_step(error_index, original_step_size);
		}
		void free_step() {
            /*
             */
			namespace bp = boost::python;
            if(final_time_reached) StopIteration();
			save_last();
			bool step_successful = false;
			num_type error_index, current_step_size;
			if(step_size < (final_time - saved_time)) {
				while(!step_successful) {
                    try_free_step(step_size);
				    current_step_size = step_size;
                   	error_index = compute_error_index(); 
					step_successful = adjust_step(error_index, step_size); 
				}
			} else {
				current_step_size = final_time - saved_time;
                final_time_reached = true;
				while(!step_successful) {
                    try_free_step(current_step_size);
				    error_index = compute_error_index(); 
					step_successful = decrease_step(error_index, current_step_size); 
                    if(!step_successful) {
                        current_step_size = step_size;
                        final_time_reached = false;
                    }
				}
			}
			system.attr("state") = bp::make_tuple(saved_time + current_step_size, 
												  temporary);
		}
		virtual void step(num_type next_time) {
            /*
             */
			namespace bp = boost::python;
			save_last();
			auto current_time = saved_time; 
            current_state = saved_state;
			num_type current_step_size, error_index=0.0;
			if(step_size < 0) step_size = next_time - current_time;
            if(std::abs(step_size) < std::numeric_limits<num_type>::epsilon()) {
                step_size = -1;
                return;
            }
			bool last_step = false;
			while(true) {
				current_step_size = next_time - current_time;
				if(step_size < current_step_size) {
					current_step_size = step_size;
					last_step = false;
				} else last_step = true; //problem here?
			    
                try_step(current_time, current_step_size);
               	error_index = compute_error_index();

				if(error_index > 1.0) {
					//calculate smaller step size
					decrease_step(error_index, current_step_size);
					continue;
				} else if(last_step) {
					break;
				} else if(error_index < 0.5) {
					//calculate larger step size
					increase_step(error_index);
				}
				current_time += current_step_size;
				std::swap(current_state, temporary);
			}
			system.attr("state") = bp::make_tuple(next_time, temporary);
		}
		boost::python::object next() {
            /*
             */
            // there is quite a bit of redundancy with base_type::next
			namespace bp = boost::python;
			if(steps.is_none()) {
                if(saved_information == BOUNDARY) step_across();
				free_step();
				if( events_occurred() ) {
					auto flagged = find_root();
					return bp::make_tuple(system.attr("state")[0], flagged);
				} else {
					return bp::make_tuple(system.attr("state")[0], bp::list());
				}
                if(not signs_verified) verify_signs();
			} else {
				return base_type::next();
			}
		}	
	};	
	template<typename stepper_type>
	class simple_stepper_instance : public stepper_wrapper {
		stepper_type stepper;

	public:
		using stepper_wrapper::stepper_wrapper;
		void step(num_type next_time) {
			save_last();
			namespace bp = boost::python;
			stepper.do_step(system_function, 
							saved_state, 
							saved_time, 
							temporary, 
							next_time - saved_time); //(sys, xin, tin, xout, h)
			system.attr("state") = bp::make_tuple(next_time, temporary);
		}
	};
	template<typename stepper_type>
	class simple_multistepper_instance : public stepper_wrapper {
		stepper_type stepper;
	
	public:
		using stepper_wrapper::stepper_wrapper;
		//'step' method is duplicated in simple_stepper_instance
		void step(num_type next_time) {
			save_last();
			namespace bp = boost::python;
			stepper.do_step(system_function, 
							saved_state, 
							saved_time, 
							temporary, 
							next_time - saved_time); //(sys, xin, tin, xout, h)
			system.attr("state") = bp::make_tuple(next_time, temporary);
		}
		virtual void reset() {
			stepper.reset();
			stepper_wrapper::reset();
		}
	};
	template<typename stepper_type, unsigned int stepper_order, unsigned int error_order>
	class variable_stepper_instance : public variable_stepper_wrapper {
		stepper_type stepper;

	public:
		using variable_stepper_wrapper::variable_stepper_wrapper;
		virtual void try_step(num_type current_time, num_type current_step_size) {
			stepper.do_step(system_function,
							current_state,
							current_time,
							temporary,
							current_step_size,
							error);
		}
		virtual void try_free_step(num_type current_step_size) {
            stepper.do_step(system_function,
                            saved_state,
                            saved_time,
                            temporary,
                            current_step_size,
                            error);
        }
        virtual void increase_step(num_type error_index) {
			if(error_index < 0.5) {
				step_size *= minimum(5, 0.9*pow(error_index, -1.0/stepper_order));
            }
		}
		virtual bool decrease_step(num_type error_index, num_type original_step_size) {
			if(error_index > 1.0) { //problems here
				//decrease step size
				step_size = original_step_size * maximum(0.2, 0.9 *
							pow(error_index, -1.0/(error_order - 1.0)));
				step_size = maximum(step_size, time_tolerance);
				return false;
			} 
			return true;
		}
	};
	template<typename stepper_type, unsigned int stepper_order, unsigned int error_order>
	class variable_multistepper_instance : public variable_stepper_wrapper {
		stepper_type stepper;

	public:
		using variable_stepper_wrapper::variable_stepper_wrapper;
		//'try_step' and 'adjust_step' methods are duplicated in variable_stepper_instance
		virtual void try_step(num_type current_time, num_type current_step_size) {
			stepper.do_step(system_function,
							current_state,
							current_time,
							temporary,
							current_step_size,
							error);
		}
		virtual void try_free_step(num_type current_step_size) {
            stepper.do_step(system_function,
                            saved_state,
                            saved_time,
                            temporary,
                            current_step_size,
                            error);
        }
        virtual void increase_step(num_type error_index) {
			if(error_index < 0.5) {
				step_size *= minimum(5, 0.9*pow(error_index, -1.0/stepper_order));
            }
		}
		virtual bool decrease_step(num_type error_index, num_type original_step_size) {
			if(error_index > 1.0) { //problems here
				//decrease step size
				step_size = original_step_size * maximum(0.2, 0.9 *
							pow(error_index, -1.0/(error_order - 1.0)));
				step_size = maximum(step_size, time_tolerance);
				return false;
			} 
			return true;
		}
		virtual void reset() {
			stepper.reset();
			variable_stepper_wrapper::reset();
		}
	};
}

//macros for exposing wrapped steppers
//assumes 'using namespace boost::python'
//no semicolon afterwards
#define LYAPUNOV_EXPOSE_STEPPER(STEPPER, WRAPPER) { \
class_< WRAPPER< STEPPER > >(#STEPPER, init<object, object>()) \
	.def("step_across", &WRAPPER< STEPPER >::step_across, (arg("new_state")=object()) ) \
	.def("reset", &WRAPPER< STEPPER >::reset_with_events) \
	.def("__iter__", pass_through) \
	.def("next", &WRAPPER< STEPPER >::next) \
	.def("use_times", &WRAPPER< STEPPER >::use_times) \
	.add_property("time_tolerance", &WRAPPER< STEPPER >::get_time_tolerance, &WRAPPER< STEPPER >::set_time_tolerance) \
	.add_property("system", &WRAPPER< STEPPER >::get_system, &WRAPPER< STEPPER >::set_system) \
	.def_readonly("tracking_events", &WRAPPER< STEPPER >::tracking_events) \
	.def_readonly("step_size", &WRAPPER< STEPPER >::get_step_size) \
	.def_readonly("status", &WRAPPER< STEPPER >::get_status); }

#define LYAPUNOV_EXPOSE_VARIABLE_STEPPER(STEPPER, WRAPPER, ORDER, ERROR) { \
class_< WRAPPER< STEPPER, ORDER, ERROR > >(#STEPPER, init<object, object>()) \
	.def("step_across", &WRAPPER< STEPPER, ORDER, ERROR >::step_across, (arg("new_state")=object()) ) \
	.def("reset", &WRAPPER< STEPPER, ORDER, ERROR >::reset_with_events) \
	.def("__iter__", pass_through) \
	.def("next", &WRAPPER< STEPPER, ORDER, ERROR >::next) \
	.def("use_times", &WRAPPER< STEPPER, ORDER, ERROR >::use_times) \
	.add_property("time_tolerance", &WRAPPER< STEPPER, ORDER, ERROR >::get_time_tolerance, &WRAPPER< STEPPER, ORDER, ERROR >::set_time_tolerance) \
	.add_property("relative_tolerance", &WRAPPER< STEPPER, ORDER, ERROR >::get_relative_tolerance, &WRAPPER< STEPPER, ORDER, ERROR >::set_relative_tolerance) \
	.add_property("absolute_tolerance", &WRAPPER< STEPPER, ORDER, ERROR >::get_absolute_tolerance, &WRAPPER< STEPPER, ORDER, ERROR >::set_absolute_tolerance) \
	.add_property("system", &WRAPPER< STEPPER, ORDER, ERROR >::get_system, &WRAPPER< STEPPER, ORDER, ERROR >::set_system) \
	.def_readonly("tracking_events", &WRAPPER< STEPPER, ORDER, ERROR >::tracking_events) \
	.def_readonly("step_size", &WRAPPER< STEPPER, ORDER, ERROR >::get_step_size) \
	.def_readonly("status", &WRAPPER< STEPPER, ORDER, ERROR >::get_status); }

//macro to help with variable order solvers
//no semicolon afterwards
#define LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(VARSTEPPER) { \
typedef ode::VARSTEPPER<1, state_type> VARSTEPPER##1; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##1, simple_multistepper_instance); \
typedef ode::VARSTEPPER<2, state_type> VARSTEPPER##2; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##2, simple_multistepper_instance); \
typedef ode::VARSTEPPER<3, state_type> VARSTEPPER##3; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##3, simple_multistepper_instance); \
typedef ode::VARSTEPPER<4, state_type> VARSTEPPER##4; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##4, simple_multistepper_instance); \
typedef ode::VARSTEPPER<5, state_type> VARSTEPPER##5; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##5, simple_multistepper_instance); \
typedef ode::VARSTEPPER<6, state_type> VARSTEPPER##6; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##6, simple_multistepper_instance); \
typedef ode::VARSTEPPER<7, state_type> VARSTEPPER##7; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##7, simple_multistepper_instance); \
typedef ode::VARSTEPPER<8, state_type> VARSTEPPER##8; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##8, simple_multistepper_instance); }

BOOST_PYTHON_MODULE(solvers) {
	using namespace boost::python;
	namespace ode = boost::numeric::odeint;
	using namespace lyapunov;

	vector_to_python_tuple vec2tup;
	vector_from_python_tuple tup2vec;
	//to_python_converter<std::vector<double>, vector_to_python_tuple>();

	typedef stepper_wrapper::state_type state_type;

	typedef ode::euler<state_type> euler;
	LYAPUNOV_EXPOSE_STEPPER(euler, simple_stepper_instance)
	typedef ode::modified_midpoint<state_type> modified_midpoint;
	LYAPUNOV_EXPOSE_STEPPER(modified_midpoint, simple_stepper_instance)
	typedef ode::runge_kutta4<state_type> runge_kutta4;
	LYAPUNOV_EXPOSE_STEPPER(runge_kutta4, simple_stepper_instance)

	typedef ode::runge_kutta_cash_karp54<state_type> cash_karp;
	LYAPUNOV_EXPOSE_VARIABLE_STEPPER(cash_karp, variable_stepper_instance, 5, 4)
	typedef ode::runge_kutta_fehlberg78<state_type> fehlberg87; //order 7 error est.
	LYAPUNOV_EXPOSE_VARIABLE_STEPPER(fehlberg87, variable_stepper_instance, 8, 7)

	typedef ode::runge_kutta_dopri5<state_type> dormand_prince;
	LYAPUNOV_EXPOSE_VARIABLE_STEPPER(dormand_prince, variable_multistepper_instance, 5, 4)

	LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_bashforth)
	//ode::adams_moulton has a weird extra argument for do_step (a buffer?)
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_moulton)
	//ode::adams_bashforth_moulton lacks a reset method for some reason (a bug?)
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_bashforth_moulton)
	
	//think about Bulirsch-Stoer solver too!
	//write a dense_stepper_wrapper?
}

