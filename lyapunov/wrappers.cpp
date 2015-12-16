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
#include <boost/python/stl_iterator.hpp>

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

    template<typename state_type>
    class python_system {
        boost::python::object system;

    public:
        python_system(boost::python::object sys) 
            : system(sys) {}
        
        std::pair<double, state_type> get() const {
            namespace bp = boost::python;
            return std::make_pair(current_time, current_state);
        }
    };
}

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
	//ode::adams_bashforth_moulton lacks a reset method until next boost release
	//LYAPUNOV_EXPOSE_VARIABLE_ORDER_STEPPER(adams_bashforth_moulton)
	
	//think about Bulirsch-Stoer solver too!
	//write a dense_stepper_wrapper?
}

