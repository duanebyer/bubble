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
	std::vector<Real> weights;
	for (std::size_t count = 0; count < 100; ++count) {
		Real weight;
		bubble::Point<DIM, Real> point;
		generator.generate(&weight, &point);
		std::cout << "point: ("
			<< point[0] << ", "
			<< point[1] << ", "
			<< point[2] << "), weight: "
			<< weight << std::endl;
		weights.push_back(weight);
	}
	Real weight_mean;
	for (Real weight : weights) {
		weight_mean += weight / weights.size();
	}
	Real weight_var;
	for (Real weight : weights) {
		weight_var += (weight - weight_mean) * (weight - weight_mean) / weights.size();
	}
	std::cout << "Relative variance: " << weight_var / (weight_mean * weight_mean) << std::endl;
	return 0;
}

