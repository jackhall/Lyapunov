#include <iostream>
#include <array>
#include <assert.h>
#include <type_traits>
//#include <boost/python.hpp>

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


int main() {
	std::cout << "a[4,2] = " << a<4>(2) << std::endl;
	std::cout << "a[6,5] = " << a<6>(5) << std::endl;
	std::cout << "a[2,1] = " << a<2>(1) << std::endl;
	return 0;
}

