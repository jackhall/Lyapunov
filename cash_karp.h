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


	struct Interval {
		double lower, upper;
		double length() const { return upper - lower; }
		double midpoint() const { return (upper + lower)/2.0; }
	}; 


	class Stepper {
		boost::python::object system;
		bool saved;
		unsigned int num_states;
		double previous_time;
		std::vector<double> previous_state, previous_error, current_error;
		std::vector<double> k1, k2, k3, k4, k5, k6;

		/*void scale(std::vector<double>& x, double factor) {
			//multiply each element in x by factor
			for(auto& i : x)
				i += factor;
		}
		
		void scale_and_add(list& x, const list& y, double factor) {
			//multiply each element in y by factor and add it to the
			//corresponding element in x
			using namespace boost::python;
			if(len(x) != len(y) or len(x) != num_states) LengthError();
			for(unsigned int i=0; i<num_states; ++i)
				x[i] += extract<double>(y[i]) * factor;
		}*/

		void scale_and_add(boost::python::list& x, 
						   const std::vector<double>& y, double factor) {
			for(unsigned int i=0; i<num_states; ++i) x[i] += y[i] * factor;
		}
		void scale_and_add(std::vector<double>& x, 
						   const std::vector<double>& y, double factor) {
			for(unsigned int i=0; i<num_states; ++i) x[i] += y[i] * factor;
		}


	public:
		Stepper() = delete;
		Stepper(boost::python::object sys) : system(sys), saved(false), 
			num_states( len(system) ), previous_time(0), previous_state(num_states),
			previous_error(num_states), current_error(num_states, 0.0),
			k1(num_states), k2(num_states), k3(num_states), k4(num_states),
			k5(num_states), k6(num_states) {} 
		Stepper(const Stepper& rhs) : system(rhs.system), saved(rhs.saved), 
			previous_time(rhs.previous_time), previous_state(rhs.previous_state),
			previous_error(rhs.previous_error), current_error(rhs.current_error),
			num_states(rhs.num_states), k1(num_states), k2(num_states), 
			k3(num_states), k4(num_states), k5(num_states), k6(num_states) {}
		Stepper& operator=(const Stepper& rhs) {
			if(this != &rhs) {
				system = boost::python::object(rhs.system);
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
		void find_root(double step_size, double min_step_size) {
			//implements a simply bisection rootfinder which passes
			//the system _through_ the boundary
			using namespace boost::python;
			//initialize interval for bisection method
			double result_time = extract<double>(system.attr("time"));
			Interval interval = {result_time-step_size, result_time};

			//save state after crossing
			std::vector<double> result_state(num_states);
			list state = extract<list>(system.attr("state"));
			list_to_vector(result_state, state);
			std::string new_mode = extract<std::string>(system.attr("mode"));
			
			//need to be able to step back across the boundary
			if(!revert()) RevertError();

			//get mode before crossing
			std::string old_mode = extract<std::string>(system.attr("mode"));
			std::string current_mode;
			while(interval.length() > min_step_size) {
				step(interval.length() / 2.0);
				current_mode = extract<std::string>(system.attr("mode"));
				if(current_mode == old_mode) 
					interval.lower = interval.midpoint();
				else if(current_mode == new_mode) {
					//We're closer to the boundary than before, on the proper side.
					interval.upper = interval.midpoint();
					result_time = extract<double>(system.attr("time"));
					list_to_vector(result_state, state);
					revert();
				} else RootFindError();
			}

			//in case the rootfind didn't end on the proper side of the boundary...
			current_mode = extract<std::string>(system.attr("mode"));
			if(current_mode == old_mode) {
				vector_to_list(state, result_state);
				system.attr("time") = result_time;
			}
		}
		void step(double h) {
			//h is the step size
			using namespace boost::python;
			list state = extract<list>(system.attr("state"));
			list_to_vector(previous_state, state);
			std::swap(previous_error, current_error);
			double t = extract<double>(system.attr("time"));

			//first slope
			list_to_vector(k1, extract<list>(system()));

			//second slope
			scale_and_add(state, k1, h*a<2>(1));
			system.attr("time") = t + c[1]*h;
			list_to_vector(k2, extract<list>(system()));

			//third slope
			vector_to_list(state, previous_state);
			scale_and_add(state, k1, h*a<3>(1));
			scale_and_add(state, k2, h*a<3>(2));
			system.attr("time") = t + c[2]*h;
			list_to_vector(k3, extract<list>(system()));

			//fourth slope
			vector_to_list(state, previous_state);
			scale_and_add(state, k1, h*a<4>(1));
			scale_and_add(state, k2, h*a<4>(2));
			scale_and_add(state, k3, h*a<4>(3));
			system.attr("time") = t + c[3]*h;
			list_to_vector(k4, extract<list>(system()));

			//fifth slope
			vector_to_list(state, previous_state);
			scale_and_add(state, k1, h*a<5>(1));
			scale_and_add(state, k2, h*a<5>(2));
			scale_and_add(state, k3, h*a<5>(3));
			scale_and_add(state, k4, h*a<5>(4));
			system.attr("time") = t + c[4]*h;
			list_to_vector(k5, extract<list>(system()));

			//sixth slope
			vector_to_list(state, previous_state);
			scale_and_add(state, k1, h*a<6>(1));
			scale_and_add(state, k2, h*a<6>(2));
			scale_and_add(state, k3, h*a<6>(3));
			scale_and_add(state, k4, h*a<6>(4));
			scale_and_add(state, k5, h*a<6>(5));
			system.attr("time") = t + c[5]*h;
			list_to_vector(k6, extract<list>(system()));

			//state update
			vector_to_list(state, previous_state);
			scale_and_add(state, k1, h*b5[0]);
			scale_and_add(state, k2, h*b5[1]);
			scale_and_add(state, k3, h*b5[2]);
			scale_and_add(state, k4, h*b5[3]);
			scale_and_add(state, k5, h*b5[4]);
			scale_and_add(state, k6, h*b5[5]);
			system.attr("time") = t + h;

			//error update
			current_error = previous_state;
			scale_and_add(current_error, k1, h*b4[0]);
			scale_and_add(current_error, k2, h*b4[1]);
			scale_and_add(current_error, k3, h*b4[2]);
			scale_and_add(current_error, k4, h*b4[3]);
			scale_and_add(current_error, k5, h*b4[4]);
			scale_and_add(current_error, k6, h*b4[5]);
			for(unsigned int i=0; i<num_states; ++i) {
				current_error[i] -= extract<double>(state[i]);
				current_error[i] *= -1.0;
			}

			//flag prior state/error as saved
			saved = true;
		}
		void list_to_vector(std::vector<double>& y, boost::python::list x) const {
			for(unsigned int i=0; i<num_states; ++i) 
				y[i] = boost::python::extract<double>(x[i]);
		}
		void vector_to_list(boost::python::list y, std::vector<double> x) {
			for(unsigned int i=0; i<num_states; ++i)
				y[i] = x[i];
		}
		void euler_step(double h) {
			using namespace boost::python;
			//h is the step size
			list state = extract<list>(system.attr("state"));
			list_to_vector(previous_state, state);
			std::swap(previous_error, current_error);
			double t = extract<double>(system.attr("time"));

			list_to_vector(k1, extract<list>(system()));

			for(unsigned int i=0; i<num_states; ++i) 
				state[i] += k1[i] * h;
			system.attr("time") = t + h;
			saved = true;
		}
		bool revert() {
			using namespace boost::python;
			if(!saved) return false;
			list state = extract<list>(system.attr("state"));
			vector_to_list(state, previous_state);
			system.attr("time") = previous_time;
			std::swap(current_error, previous_error);
			saved = false;
			return true;
		}
		void save() {
			using namespace boost::python;
			previous_error = current_error;
			list state = extract<list>(system.attr("state"));
			list_to_vector(previous_state, state);
			previous_time = extract<double>(system.attr("time"));
			saved = true;
		}
		boost::python::list get_error() const { 
			boost::python::list error_list;
			for(auto x : current_error) error_list.append(x);
			return error_list; 
		}
		double get_step_size() const {
			if(saved) 
				return boost::python::extract<double>(system.attr("time")) - previous_time;
			else return 0.0;
		}
		unsigned int get_num_states() const { return num_states; }
	};

} //namespace lyapunov

#endif

