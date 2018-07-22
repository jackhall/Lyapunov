#ifndef COMMON_CPP
#define COMMON_CPP

#include <boost/python.hpp>

namespace lyapunov {
    namespace bp = boost::python

    // convenience functions to call when I want to throw a particular kind of error
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


    // conversions between C++ vectors and python tuples
    struct vector_to_python_tuple {
        vector_to_python_tuple() {
            bp::to_python_converter<const std::vector<double>, vector_to_python_tuple>();
        }

        static PyObject* convert(const std::vector<double>& x) {
            auto new_tuple = PyTuple_New(x.size());
            for(int i=x.size()-1; i>=0; --i) {
                PyTuple_SetItem(new_tuple, i, PyFloat_FromDouble(x[i]));
            }
            return bp::incref(new_tuple);
        }
    };

    struct vector_from_python_tuple {
        vector_from_python_tuple() {
            bp::converter::registry::push_back(&convertible, &construct, bp::type_id<std::vector<double>>());
        }

        static void* convertible(PyObject* obj_ptr) {
            if (!PyTuple_Check(obj_ptr)) return 0;
            return obj_ptr;
        }

        static void construct(PyObject* obj_ptr,
            bp::converter::rvalue_from_python_stage1_data* data) {
            assert(PyTuple_Check(obj_ptr));
            unsigned int length = PyTuple_Size(obj_ptr);
            void* storage = ((bp::converter::rvalue_from_python_storage< std::vector<double> >*)data)->storage.bytes;
            new (storage) std::vector<double>(length);
            for(unsigned int i=0; i<length; ++i)
                static_cast< std::vector<double>* >(storage)->at(i)
                    = PyFloat_AsDouble(PyTuple_GetItem(obj_ptr, i));
            data->convertible = storage;
        }
    };


    // TODO: Wrap python system objects and iterables in C++ so only wrappers have to mess with boost python.
    struct ODEIterator {
    /* Abstract base class for anything wrapped as a stepper.
     */
        virtual bp::str get_status() const = 0;

        virtual double get_tolerance() const = 0;
        virtual void set_tolerance(double new_tolerance) = 0;

        virtual bp::object get_system() const = 0;
        virtual void set_system(bp::object new_system) = 0;

        virtual double get_step_size() const = 0;
        virtual void use_times(bp::object time) = 0;

        virtual void reset() = 0;
        void step_across(bp::object new_state = bp::object()) {
            _step_across(new_state);
        }
        virtual bp::tuple next() = 0;
    };
}

#endif
