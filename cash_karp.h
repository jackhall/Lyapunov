#ifndef cash_karp_h
#define cash_karp_h

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

#include <algorithm>
#include <array>
#include <type_traits>
#include <boost/python.hpp>
namespace lyapunov {
	//Butcher Tableau for Cash-Karp Runge-Kutta method
	template<unsigned int row>
	double a(unsigned int column) {
		static_assert(row < 7 && row > 1, "Row out of bounds.");
		return a<row-1>(column + row - 2);
	}

	constexpr std::array<double, 15>
	elems = {{0.2,
			  (3.0/40.0), (9.0/40.0), 
			  0.3, -0.9, 1.2,
			  -(11.0/54.0), 2.5, -(70.0/27.0), (35.0/27.0),
			  (1631.0/55296.0), (175.0/512.0), (575.0/13824.0), (44275.0/110592.0), (253.0/4096.0)}};

	template<>
	double a<2>(unsigned int column) { return elems[column-1]; }

	static constexpr std::array<double, 6>  
	c = {{0.0, 0.2, 0.3, 0.6, 1.0, (7.0/8.0)}};

	static constexpr std::array<double, 6>
	b5 = {{(37.0/378.0), 0.0, (250.0/621.0), (125.0/594.0), 0.0, (512.0/1771.0)}};

	static constexpr std::array<double, 6>
	b4 = {{(2825.0/27648.0), 0.0, (18575.0/48384.0), (13525.0/55296.0), (277.0/14336.0), 0.25}};

	//needed for element-wise list adding
	void LengthError() { PyErr_SetString(PyExc_IndexError, "List lengths don't match."); }
	void RootFindError() { PyErr_SetString(PyExc_RuntimeError, "Multiple roots detected."); }
	void RevertError() { PyErr_SetString(PyExc_RuntimeError, "No state saved for necessary revert."); }

	//rvalue converters for std::vector<double> and boost::python::tuple
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
	

	struct Interval {
		double lower, upper;
		double length() const { return upper - lower; }
		double midpoint() const { return (upper + lower)/2.0; }
	}; 

	
	class Stepper {
		boost::python::object system;
		bool saved, state_is_property;
		unsigned int num_states;
		double previous_time;
		std::vector<double> previous_state, previous_error, current_error;
		std::vector<double> k1, k2, k3, k4, k5, k6;

		void scale_and_add(boost::python::tuple& x, 
						   const std::vector<double>& y, double factor) {
			for(unsigned int i=0; i<num_states; ++i) x[i] += y[i] * factor;
		}
		void scale_and_add(std::vector<double>& x, 
						   const std::vector<double>& y, double factor) {
			for(unsigned int i=0; i<num_states; ++i) x[i] += y[i] * factor;
		}

	public:
		Stepper() = delete;
		Stepper(boost::python::object sys) 
		  : system(sys), 
			saved(false), 
			num_states( len(system.attr("state")) ), 
			previous_time(0), 
			previous_state(num_states), 
			previous_error(num_states), 
			current_error(num_states, 0.0), 
			k1(num_states), k2(num_states), k3(num_states), 
			k4(num_states), k5(num_states), k6(num_states) {} 
		Stepper(const Stepper& rhs) 
		  : system(rhs.system), 
			saved(rhs.saved), 
			num_states(rhs.num_states),
			previous_time(rhs.previous_time), 
			previous_state(rhs.previous_state), 
			previous_error(rhs.previous_error), 
			current_error(rhs.current_error), 
			k1(num_states), k2(num_states), k3(num_states), 
			k4(num_states), k5(num_states), k6(num_states) {}
		Stepper& operator=(const Stepper& rhs) {
			if(this != &rhs) {
				set_system(rhs.system);
				saved = rhs.saved;
				previous_time = rhs.previous_time;
				previous_state = rhs.previous_state;
				previous_error = rhs.previous_error;
				current_error = rhs.current_error;
				num_states = rhs.num_states;
			}
		}
		~Stepper() = default;
		
		boost::python::object get_system() const { return system; }
		void set_system(boost::python::object new_system) { system = new_system; }
		void find_root(boost::python::list events, double min_step_size) {
			//implements a simple bisection rootfinder
			//  Initialize Interval
			//  Revert
			//  Loop over events: 
			//    call and store the result to a vector
			//  Loop while Inverval.length() > min_step_size and signchanges.count(True) > 1
			//    call step to midpoint of Interval
			//    Loop over events: 
			//      check for sign changes, store as seq of bools
			//    If a sign changed: 
			//    	Revert
			//    	set interval.upper to midpoint
			//    Else:
			//      set interval.lower to midpoint
			//  If event.step_through is True, step to end of Interval
			//  return event
			//this algorithm will loop forever if two events occur simultaneously!
			//If two events still occur simultaneously, call both flag methods? 
			//If the two events differ on step_through, what to do?
			using namespace boost::python;
			//initialize interval for bisection method
			object time_obj = system.attr("time");
			double result_time = extract<double>(time_obj);
			Interval interval = {previous_time, result_time};

			//save state after crossing
			//should system.attr("state") be called to check every time?
			tuple result_state = extract<tuple>(system.attr("state"));
			
			//need to be able to step back across the boundary
			if(!revert()) RevertError();

			//get mode before crossing
			std::string old_mode = extract<std::string>(mode_obj);
			std::string current_mode;
			while(interval.length() > min_step_size) {
				step(interval.length() / 2.0);
				current_mode = extract<std::string>(mode_obj);
				if(current_mode == old_mode) 
					interval.lower = interval.midpoint();
				else if(current_mode == new_mode) {
					//We're closer to the boundary than before, on the proper side.
					interval.upper = interval.midpoint();
					result_time = extract<double>(time_obj);
					result_state = extract<tuple>(system.attr("state"));
					revert();
				} else RootFindError();
			}

			//in case the rootfind didn't end on the proper side of the boundary...
			current_mode = extract<std::string>(mode_obj);
			if(current_mode == old_mode) {
				system.attr("state") = result_state;
				time_obj = object(result_time);
			}
		}
		void step(double h) {
			//h is the step size
			using namespace boost::python;
			//extract state directly to vector?
			std::vector<double> state = 
				extract< std::vector<double> >(system.attr("state")); 
			previous_state = state;
			std::swap(previous_error, current_error);
			double t = extract<double>(system.attr("time"));

			//first slope
			//rename conversion routines!
			k1 = extract< std::vector<double> >(system());

			//second slope
			for(unsigned int i=0; i<num_states; ++i)
				state[i] += h*a<2>(1) * k1[i];
			system.attr("state") = state;
			system.attr("time") = t + c[1]*h;
			k2 = extract< std::vector<double> >(system());

			//third slope
			for(unsigned int i=0; i<num_states; ++i)
				state[i] = previous_state[i] + h*(a<3>(1)*k1[i] + a<3>(2)*k2[i]);
			system.attr("state") = state;
			system.attr("time") = t + c[2]*h;
			k3 = extract< std::vector<double> >(system());

			//fourth slope
			for(unsigned int i=0; i<num_states; ++i)
				state[i] = previous_state[i] + h*(a<4>(1)*k1[i] + a<4>(2)*k2[i] 
												+ a<4>(3)*k3[i]);
			system.attr("state") = state;
			system.attr("time") = t + c[3]*h;
			k4 = extract< std::vector<double> >(system());

			//fifth slope
			for(unsigned int i=0; i<num_states; ++i)
				state[i] = previous_state[i] + h*(a<5>(1)*k1[i] + a<5>(2)*k2[i] 
												+ a<5>(3)*k3[i] + a<5>(4)*k4[i]);
			system.attr("state") = state;
			system.attr("time") = t + c[4]*h;
			k5 = extract< std::vector<double> >(system());

			//sixth slope
			for(unsigned int i=0; i<num_states; ++i)
				state[i] = previous_state[i] + h*(a<6>(1)*k1[i] + a<6>(2)*k2[i] 
												+ a<6>(3)*k3[i] + a<6>(4)*k4[i]
												+ a<6>(5)*k5[i]);
			system.attr("state") = state;
			system.attr("time") = t + c[5]*h;
			k6 = extract< std::vector<double> >(system());

			//state update
			for(unsigned int i=0; i<num_states; ++i)
				state[i] = previous_state[i] + h*(b5[0]*k1[i] + b5[1]*k2[i] 
												+ b5[2]*k3[i] + b5[3]*k4[i]
												+ b5[4]*k5[i] + b5[5]*k6[i]);
			system.attr("state") = state;
			system.attr("time") = t + h;

			//error update
			for(unsigned int i=0; i<num_states; ++i) {
				current_error[i] = -previous_state[i] - h*(b4[0]*k1[i] + b4[1]*k2[i] 
														 + b4[2]*k3[i] + b4[3]*k4[i]
														 + b4[4]*k5[i] + b4[5]*k6[i]);
				current_error[i] += state[i];
			}

			//flag prior state/error as saved
			saved = true;
		}
		void euler_step(double h) {
			using namespace boost::python;
			//h is the step size
			std::vector<double> state = extract< std::vector<double> >(system.attr("state"));
			previous_state = state;
			std::swap(previous_error, current_error);
			double previous_time = extract<double>(system.attr("time"));

			k1 = extract< std::vector<double> >(system());

			for(unsigned int i=0; i<num_states; ++i) 
				state[i] += k1[i] * h;
			system.attr("state") = state;
			system.attr("time") = previous_time + h;
			saved = true;
		}
		bool revert() {
			using namespace boost::python;
			if(!saved) return false;
			system.attr("state") = previous_state;
			system.attr("time") = previous_time;
			std::swap(current_error, previous_error);
			saved = false;
			return true;
		}
		void save() {
			using namespace boost::python;
			std::swap(previous_error, current_error);
			previous_state = extract< std::vector<double> >(system.attr("state"));
			previous_time = extract<double>(system.attr("time"));
			saved = true;
		}
		boost::python::tuple get_error() const { return boost::python::tuple(current_error); }
		double get_step_size() const {
			if(saved) 
				return boost::python::extract<double>(system.attr("time")) - previous_time;
			else return 0.0;
		}
		unsigned int get_num_states() const { return num_states; }
	};

} //namespace lyapunov

#endif

