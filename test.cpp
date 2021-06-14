#include <cmath>
#include <iostream>

#include "bubble.hpp"

using Real = double;
static bubble::Dim const DIM = 3;
static Real const PI = 3.1415926;

Real func_test(bubble::Point<DIM, Real> x) {
	Real result = 1;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		result *= std::sin(PI * x[dim]);
	}
	return result;
}

int main(int argc, char** argv) {
	auto builder = bubble::make_builder<DIM, Real>(func_test);
	std::cout << "Exploring." << std::endl;
	builder.explore();
	std::cout << "Tuning." << std::endl;
	builder.tune();
	std::cout << "Total prime: " << builder.prime() << std::endl;
	auto generator = bubble::make_generator(builder);
	std::cout << "Generating." << std::endl;
	for (std::size_t count = 0; count < 100; ++count) {
		Real weight;
		bubble::Point<DIM, Real> point;
		generator.generate(&weight, &point);
		std::cout << "point: ("
			<< point[0] << ", "
			<< point[1] << ", "
			<< point[2] << "), weight: "
			<< weight << std::endl;
	}
	return 0;
}

