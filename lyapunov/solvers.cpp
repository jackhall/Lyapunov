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
            for(int i=x.size()-1; i>=0; --i) {
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
    class stepper_iterator {
        /* Represents the state machine aspect of a stepper. As an abstract base class, it
         * needs to be subclassed with a particular stepper so that step and reset_stepper
         * can be defined.
         */
    public:
        typedef std::vector<double> state_type;
    
    protected:
        state_type current_state, saved_state;
        double desired_time, current_time, saved_time, time_tolerance;

        state_type event_signs;
        bool signs_verified;  // true if all event_signs are nonzero

        enum stepper_state_type {CLEAN, LAST, BOUNDARY};
        stepper_state_type stepper_state;

		boost::python::object system, steps;
		std::function<void(const state_type&, state_type&, double)> system_function;

    public:
        stepper_iterator() = delete;
		stepper_iterator(boost::python::object sys, 
					     boost::python::object time) 
            : time_tolerance(0.00001) { //10 microseconds ... make this relative?
            set_system(sys);
            use_times(time);
        }

		boost::python::str get_status() const {
            /* Returns what point the stepper is currently storing.
             */
			switch(stepper_state) {
				case CLEAN:
					return "nothing";
				case LAST:
					return "last step";
				case BOUNDARY:
					return "boundary";
				default:
					return "status undefined!";
			}
		}
        bool tracking_events() const {
            return not event_signs.empty();
        }
        double get_time_tolerance() const { return time_tolerance; }
		void set_time_tolerance(double new_tolerance) {
			if(new_tolerance >= 0) time_tolerance = new_tolerance;
			else ValueError("Positive values only.");
		}
		boost::python::object get_system() const { return system; }
        void set_system(boost::python::object new_system) {
            /* Tell the stepper to integrate the given system.
             */
            namespace bp = boost::python;
            system = new_system;
            system_function = [this](const state_type& x, state_type& dx, const double t) { 
                system.attr("state") = bp::make_tuple(t, x);
                dx = bp::extract<state_type>(system()); 
            };
            saved_time = 0.0;
            saved_state.resize(bp::len(system.attr("state")[1]));
            stepper_state = CLEAN;
            // State and event information will be updated when next is next called. 
        }
        virtual void use_times(boost::python::object time) {
            /* Take user input for time data. For a fixed stepper, this means an iterable
             * of times. For a variable stepper, it could also mean a single desired end 
             * time. If a variable stepper is given
             */
			if(PyObject_HasAttrString(time.ptr(), "__iter__")) {
                steps = time.attr("__iter__")();
                desired_time = boost::python::extract<double>(steps.attr("next"));
            } else ValueError("Provide a sequence of times at which to approximate the ODE.");
        }
        double get_step_size() const { return desired_time - current_time; }
        virtual bool error_acceptable() { return true; }
        virtual bool goal_achieved() { 
            /* Fetch a new time from the array thereof.
             * Always returns true for fixed-step solvers.
             */
            desired_time = boost::python::extract<double>( steps.attr("next")() );
            return true;
        }
		bool event_occurred() const {
            /* Check to see whether any event functions crossed zero since the last step.
             */
			namespace bp = boost::python;
			if(tracking_events()) {
				bp::object event, iter = system.attr("events").attr("__iter__")();
				for(auto x : event_signs) {
					event = iter.attr("next")();
					if(x * bp::extract<double>( event() ) < 0) return true;
				}
			} 
			return false;
		}
        boost::python::list find_root() {
			/* Implements a simple bisection rootfinder.
             * Moves the system to within time_tolerance of the event boundary
             * and saves a point just across the boundary.
             *
             * If more than one event is detected, the rootfinder will stop at the
             * first one it finds. If the event boundaries are within time_tolerance
             * of each other, both get detected.
             */
			namespace bp = boost::python;

			// Initialize interval and step size goal.
			Interval interval = {current_time, saved_time};

			// main rootfinding loop	
            std::vector<double> lower_state = current_state;
			while(interval.length() > time_tolerance) {
				// Step to midpoint of interval.
				step(interval.midpoint());

				// Check for sign changes.
				if( event_occurred() ) {
                    // Narrow focus to the first half of the interval and jump back.
					interval.upper = interval.midpoint();
                    saved_time = interval.lower;
                    saved_state = lower_state;
                    swap_states();
				} else {
                    // Narrow focus to the second half of the interval.
                    interval.lower = interval.midpoint();
                    lower_state = current_state;
                }
			}

			// Find out which events changed sign and return them.
            swap_states();
			bp::list flagged;
			bp::object event, iter = system.attr("events").attr("__iter__")();
			for(auto x : event_signs) {
				event = iter.attr("next")();
				if(x * bp::extract<double>(event()) < 0) 
					flagged.append(event);
			}

			// Save state for step_across.
            swap_states();
			return flagged; 
		}
        double next_time() const { return desired_time; }
        virtual void step(double new_time) = 0;
        virtual void reset_stepper() {}
        void reset() {
            /* Erase saved event signs and recheck the number of events.
             * Also wipe any internal states in the stepper.
             */
            reset_stepper();
            signs_verified = false;
            event_signs.clear();
            if(PyObject_HasAttrString(system.ptr(), "events"))
                event_signs.resize(len(system.attr("events")), 0);
            stepper_state = CLEAN;
        }
        void save_current() {
            /* Save the current point in anticipation of stepping to the next.
             */
			namespace bp = boost::python;
			bp::object state_tup = system.attr("state");
			saved_time = bp::extract<double>(state_tup[0]);
			saved_state = bp::extract<state_type>(state_tup[1]);
            stepper_state = LAST;
		}
        void swap_states() {
            /* Move the system and the stepper to the currently saved time.
             */
            if(stepper_state == CLEAN)
                RuntimeError("Internal Error: No saved state to swap.");
            std::swap(current_state, saved_state);
            std::swap(current_time, saved_time);
            update_system();
        }
        void update_system() {
            // Move the system to the current state and time.
            system.attr("state") = make_tuple(current_time, current_state);
        }
        void step_across(boost::python::object new_state = boost::python::object()) {
            /* Jump the system across the event boundary discontinuously. Such a 
             * jump means that no state will be saved after the jump, event detection
             * needs to be reset, and any internal stepper state needs to be reset. 
             * The optional state argument is for state only - not time.
             */
            if(stepper_state != BOUNDARY)
                RuntimeError("No boundary to step across.");  
                
            // Allow the user to jump across state space, but not time.
            if(not new_state.is_none())
                saved_state = boost::python::extract<state_type>(new_state);

            // Move the system and stepper across the boundary.
            swap_states();
            reset();
        }
        boost::python::tuple next() {
            /* This method is exposed to python for iteration. It includes most
             * of the important state machine logic.
             *
             * The state machine always maintains a second system state. During normal 
             * solving, this means that the stepper can always undo a step. When an 
             * event occurs the stepper already has a bounded interval. The rootfinder 
             * then narrows this interval until it's no wider than time_tolerance. The 
             * current state is considered to be the one just before the event occurs. 
             * The saved state is just after. 
             */
            namespace bp = boost::python;
            if(stepper_state == BOUNDARY)  // if user didn't do this manually
                step_across(); 
            
            do {
                // Ensure the starting sign of each event is known.
                if(not signs_verified) 
                    verify_signs();
               
                // Take a step.
                save_current();
                step(next_time());

                // Check for events and call the rootfinder if necessary.
                if(event_occurred()) {
                    // Jump back to the safe side of the boundary.
                    swap_states();
                    stepper_state = BOUNDARY;

                    // Figure out which event(s) occurred first and where.
                    auto flagged = find_root();
				    return bp::make_tuple(system.attr("state")[0], flagged);
                }
                // These conditions are always met for simple steppers, but they
                // will be important for variable steppers. 
            } while(not error_acceptable() or not goal_achieved()); 

			return bp::make_tuple(system.attr("state")[0], boost::python::list());
        }
        void verify_signs() {
            /* Ensure that all event signs are nonzero;
             * that is: make sure that event functions equal to
             * zero at the initial conditions still get detected.
             */
            namespace bp = boost::python;
            unsigned int zero_sign_count = 0;
            for(unsigned int i=0; i<event_signs.size(); ++i) {
                if(event_signs[i] == 0) {
                    event_signs[i] = bp::extract<double>(system.attr("events")[i]);
                    if(event_signs[i] == 0) ++zero_sign_count;
                }
            }
            if(zero_sign_count == 0) 
                signs_verified = true;
        }
    };
    class variable_stepper_iterator : public stepper_iterator {
    protected:
		double absolute_tolerance, relative_tolerance;
		state_type error;
		double step_size, final_time;
        bool free_iteration;

        double get_error_index() {
            double error_index = 0.0;
            for(int i=error.size()-1; i>=0; --i) {
				error_index = maximum(error_index, std::abs(error[i]) /
							  (absolute_tolerance + 
                               relative_tolerance*std::abs(current_state[i])));
			}
            return error_index;
        }
        virtual unsigned int get_stepper_order() const = 0;
        virtual unsigned int get_error_order() const = 0;

    public:
        using stepper_iterator::stepper_iterator;

		virtual void use_times(boost::python::object time) {
            /* For variable_steppers, simulation time can be specified either as
             * a list (as with simple_steppers) or a float. List calls are simply
             * delegated to the stepper_wrapper version. Float calls give the stepper
             * a 'goal' time at which they'll stop simulating, and otherwise let it
             * choose step sizes freely. 
             */
			namespace bp = boost::python;
            PyObject* time_ptr = time.ptr();
			if(!PySequence_Check(time_ptr) && PyNumber_Check(time_ptr)) {
				steps = boost::python::object();
				final_time = bp::extract<double>(time);
                if(step_size < 0) {
				    step_size = 0.001*(final_time - 
							bp::extract<double>(system.attr("state")[0]));
                }
			} else stepper_iterator::use_times(time);
		}
        double get_step_size() const { return step_size; }
        virtual bool error_acceptable() {
            /* Computes an error index and checks that it's within bounds. If the error
             * is too large, update the step size, back up, and return false. If it's 
             * too small, update the step size and return true. 
             */
            double error_index = get_error_index();
			if(error_index < 0.5)  // possibly increase step size
				step_size *= minimum(5, 0.9*pow(error_index, -1.0/get_stepper_order()));
            else if(error_index > 1.0) {  // possibly decrease step size
				step_size = step_size * maximum(0.2, 0.9 * pow(error_index, 
                    -1.0/(get_error_order() - 1.0)));
				step_size = maximum(step_size, time_tolerance);

                // revert to saved step
                return false;
            }
			return true;
        }
        virtual bool goal_achieved() {
            /* If no time iterable is set, always return true. If one is set,
             * return true only if the stepper has reached the current next step, and
             * assign a new next step. 
             */
            if(free_iteration) {
                if(current_time >= final_time) StopIteration();
                return true;
            } else if(current_time >= desired_time) {
                return stepper_iterator::goal_achieved();
            } else return false;
        }
        double next_time() const {
            /* Returns the next time to step to, abbreviating the step if needed.
             */
            return minimum(step_size + current_time, desired_time);
        }
        double get_relative_tolerance() const { return relative_tolerance; }
        void set_relative_tolerance(double new_relative_tolerance) {
            if(new_relative_tolerance > 0.0) 
                relative_tolerance = new_relative_tolerance;
            else ValueError("Tolerances cannot be negative.");
        }
        double get_absolute_tolerance() const { return absolute_tolerance; }
        void set_absolute_tolerance(double new_absolute_tolerance) {
            if(new_absolute_tolerance > 0.0) 
                absolute_tolerance = new_absolute_tolerance;
            else ValueError("Tolerances cannot be negative.");
        }
    };
	template<typename stepper_type>
	struct simple_stepper : public stepper_iterator {
        /* For basic solvers that:
         *  - do not vary their own step size 
         *  - do not make use of prior state information
         */
		stepper_type stepper;

		using stepper_iterator::stepper_iterator;
		virtual void step(double new_time) {
            namespace bp = boost::python;
			stepper.do_step(system_function, 
							saved_state, 
							saved_time, 
							current_state, 
							new_time - saved_time); //(sys, xin, tin, xout, h)
			system.attr("state") = bp::make_tuple(new_time, current_state);
		}
	};
	template<typename stepper_type>
	struct simple_multistepper : public simple_stepper<stepper_type> {
        /* For solvers that:
         *  - do not vary their own step size 
         *  - make use of prior state information
         */
		using simple_stepper<stepper_type>::simple_stepper;
		virtual void reset_stepper() { simple_stepper<stepper_type>::stepper.reset(); }
	};
	template<typename stepper_type, unsigned int stepper_order, unsigned int error_order>
	struct variable_stepper : public variable_stepper_iterator {
        /* For solvers that:
         *  - vary their own step size 
         *  - do not make use of prior state information
         */
		stepper_type stepper;

		using variable_stepper_iterator::variable_stepper_iterator;
		virtual void step(double new_time) {
            namespace bp = boost::python;
			stepper.do_step(system_function,
							saved_state,
							saved_time,
							current_state,
							new_time - saved_time,
							error); // (sys, xin, tin, xout, h, err)
            system.attr("state") = bp::make_tuple(new_time, current_state);
		}
        virtual unsigned int get_stepper_order() const { return stepper_order; }
        virtual unsigned int get_error_order() const { return error_order; }
	};
	template<typename stepper_type, unsigned int stepper_order, unsigned int error_order>
	struct variable_multistepper : 
      public variable_stepper<stepper_type, stepper_order, error_order> {
        /* For solvers that both:
         *  - vary their own step size 
         *  - make use of prior state information
         */
        typedef variable_stepper<stepper_type, stepper_order, error_order> 
                variable_stepper_type;

		using variable_stepper_type::variable_stepper;
		virtual void reset() { variable_stepper_type::stepper.reset(); }
	};
}

//macros for exposing wrapped steppers
//assumes 'using namespace boost::python'
//no semicolon afterwards
#define LYAPUNOV_EXPOSE_STEPPER(STEPPER, WRAPPER) { \
class_< WRAPPER< STEPPER > >(#STEPPER, init<object, object>()) \
	.def("step_across", &WRAPPER< STEPPER >::step_across, (arg("new_state")=object()) ) \
	.def("reset", &WRAPPER< STEPPER >::reset) \
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
	.def("reset", &WRAPPER< STEPPER, ORDER, ERROR >::reset) \
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
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##1, simple_multistepper); \
typedef ode::VARSTEPPER<2, state_type> VARSTEPPER##2; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##2, simple_multistepper); \
typedef ode::VARSTEPPER<3, state_type> VARSTEPPER##3; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##3, simple_multistepper); \
typedef ode::VARSTEPPER<4, state_type> VARSTEPPER##4; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##4, simple_multistepper); \
typedef ode::VARSTEPPER<5, state_type> VARSTEPPER##5; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##5, simple_multistepper); \
typedef ode::VARSTEPPER<6, state_type> VARSTEPPER##6; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##6, simple_multistepper); \
typedef ode::VARSTEPPER<7, state_type> VARSTEPPER##7; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##7, simple_multistepper); \
typedef ode::VARSTEPPER<8, state_type> VARSTEPPER##8; \
LYAPUNOV_EXPOSE_STEPPER(VARSTEPPER##8, simple_multistepper); }

BOOST_PYTHON_MODULE(solvers) {
	using namespace boost::python;
	namespace ode = boost::numeric::odeint;
	using namespace lyapunov;

	vector_to_python_tuple vec2tup;
	vector_from_python_tuple tup2vec;
	//to_python_converter<std::vector<double>, vector_to_python_tuple>();

	typedef stepper_iterator::state_type state_type;

	typedef ode::euler<state_type> euler;
	LYAPUNOV_EXPOSE_STEPPER(euler, simple_stepper)
	typedef ode::modified_midpoint<state_type> modified_midpoint;
	LYAPUNOV_EXPOSE_STEPPER(modified_midpoint, simple_stepper)
	typedef ode::runge_kutta4<state_type> runge_kutta4;
	LYAPUNOV_EXPOSE_STEPPER(runge_kutta4, simple_stepper)

	typedef ode::runge_kutta_cash_karp54<state_type> cash_karp;
	LYAPUNOV_EXPOSE_VARIABLE_STEPPER(cash_karp, variable_stepper, 5, 4)
	typedef ode::runge_kutta_fehlberg78<state_type> fehlberg87; //order 7 error est.
	LYAPUNOV_EXPOSE_VARIABLE_STEPPER(fehlberg87, variable_stepper, 8, 7)

	typedef ode::runge_kutta_dopri5<state_type> dormand_prince;
	LYAPUNOV_EXPOSE_VARIABLE_STEPPER(dormand_prince, variable_multistepper, 5, 4)

	LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_bashforth)
	//ode::adams_bashforth_moulton lacks a reset method for some reason (a bug?)
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_bashforth_moulton)
	
	//think about Bulirsch-Stoer solver too!
	//write a dense_stepper_wrapper?
}

