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

BOOST_PYTHON_MODULE(stepper) {
	using namespace boost::python;

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


	class Stepper {
		object system;
		bool saved;
		double previous_time;
		list previous_state, current_error, previous_error;

		list scale(list x, double factor) {
			auto n = len(x);
			for(unsigned int i=0; i<n; ++i) 
				x[i] += factor;
			return x;
		}

		list scale_and_add(list x, list y, double factor) {
			if(len(x) != len(y)) LengthError();
			auto n  = len(x);
			for(unsigned int i=0; i<n; ++i)
				x[i] += extract<double>(y[i]) * factor;
			return x;
		}

	public:
		Stepper(

		object get_system() const { return system; }
		void set_system(object new_system) { system = new_system; }
		void step(double h) {
			//h is the step size
			previous_state = extract<list>(system.attr("state"));
			previous_error = current_error;
			auto t = extract<double>(system.attr("time"));

			//slope calculations
			list k1 = scale(system(), h); 
			
			system.attr("state") = scale_and_add(total, k1, a<2>(1));
			system.attr("time") = t + c[1]*h;
			list k2 = scale(system(), h);
			
			list total(previous_state);
			scale_and_add(total, k1, a<3>(1));
			scale_and_add(total, k2, a<3>(2));
			system.attr("state") = total;
			system.attr("time") = t + c[2]*h;
			list k3 = scale(system(), h);
			
			total = list(previous_state);
			scale_and_add(total, k1, a<4>(1));
			scale_and_add(total, k2, a<4>(2));
			scale_and_add(total, k3, a<4>(3));
			system.attr("state") = total;
			system.attr("time") = t + c[3]*h;
			list k4 = scale(system(), h);
	
			total = list(previous_state);
			scale_and_add(total, k1, a<5>(1));
			scale_and_add(total, k2, a<5>(2));
			scale_and_add(total, k3, a<5>(3));
			scale_and_add(total, k4, a<5>(4));
			system.attr("state") = total;
			system.attr("time") = t + c[4]*h;
			list k5 = scale(system(), h);

			total = list(previous_state);
			scale_and_add(total, k1, a<6>(1));
			scale_and_add(total, k2, a<6>(2));
			scale_and_add(total, k3, a<6>(3));
			scale_and_add(total, k4, a<6>(4));
			scale_and_add(total, k5, a<6>(5));
			system.attr("state") = total;
			system.attr("time") = t + c[5]*h;
			list k6 = scale(system(), h);

			//state update
			total = list(previous_state);
			scale_and_add(total, k1, b5[0]);
			scale_and_add(total, k2, b5[1]);
			scale_and_add(total, k3, b5[2]);
			scale_and_add(total, k4, b5[3]);
			scale_and_add(total, k5, b5[4]);
			scale_and_add(total, k6, b5[5]);
			system.attr("state") = total;
			system.time("time") = t + h;

			//error
			current_error = list(previous_state);
			scale_and_add(current_error, k1, b4[0]);
			scale_and_add(current_error, k2, b4[1]);
			scale_and_add(current_error, k3, b4[2]);
			scale_and_add(current_error, k4, b4[3]);
			scale_and_add(current_error, k5, b4[4]);
			scale_and_add(current_error, k6, b4[5]);
			auto n = len(current_error);
			for(unsigned int i=0; i<n; ++i) {
				current_error[i] -= total[i];
				current_error[i] *= -1.0;
			}

			//flag prior state/error as saved
			saved = true;
		}

		bool revert() {
			if(!saved) return false;
			system.attr("state") = previous_state;
			current_error = previous_error;
			return true;
		}

		list get_error() const { return current_error; }
	};

	class_<Stepper>("Stepper")
		.def("step", &Stepper::step)
		.def("revert", &Stepper::revert)
		.add_property("error", &Stepper::get_error);
		.add_property("system", &Stepper::get_system, &Stepper::set_system);
}

