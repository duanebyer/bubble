#include <cmath>
#include <iomanip>
#include <iostream>

#include "bubble.hpp"

using Real = double;
static bubble::Dim const DIM = 3;
static Real const PI = 3.1415926;

Real func_test(bubble::Point<DIM, Real> x);

int thread_count();

int main() {
	std::cout << std::setprecision(6) << std::fixed;
	std::cout << "Number of threads: " << thread_count() << std::endl;
	auto builder = bubble::make_builder<DIM, Real>(func_test);
	std::cout << "Exploring." << std::endl;
	std::size_t underexplored = builder.explore();
	if (underexplored != 0) {
		std::cout << "Underexplored " << underexplored << " cells." << std::endl;
	}
	std::cout << "Tuning." << std::endl;
	std::size_t undertuned = builder.tune();
	if (undertuned != 0) {
		std::cout << "Undertuned by " << undertuned << " samples." << std::endl;
	}
	Real est_eff = builder.mean() / builder.prime();
	Real est_rel_var = 1 / (est_eff * est_eff) - 1;
	std::cout << "Total prime:    " << builder.prime() << std::endl;
	std::cout << "Est. mean:      " << builder.mean() << std::endl;
	std::cout << "Est. eff.:      " << est_eff << std::endl;
	std::cout << "Est. rel. var.: " << est_rel_var << std::endl;
	std::cout << "Scaling exp.  : " << builder.scale_exp() << std::endl;
	auto generator = bubble::make_generator(builder);
	std::cout << "Generating." << std::endl;
	std::vector<Real> weights;
	std::cout << "Example points:" << std::endl;
	for (std::size_t count = 0; count < 10000; ++count) {
		Real weight;
		bubble::Point<DIM, Real> point;
		generator.generate(&weight, &point);
		if (count < 10) {
			std::cout
				<< "\tw: " << weight << ", x: ("
				<< point[0] << ", "
				<< point[1] << ", "
				<< point[2] << ")" << std::endl;
		}
		weights.push_back(weight);
	}
	Real weight_mean = 0;
	Real weight_sq = 0;
	for (Real weight : weights) {
		weight_mean += weight / weights.size();
		weight_sq += weight * weight / weights.size();
	}
	Real weight_var = 0;
	for (Real weight : weights) {
		weight_var += (weight - weight_mean) * (weight - weight_mean) / weights.size();
	}
	Real weight_eff = weight_mean / std::sqrt(weight_sq);
	Real weight_rel_var = weight_var / (weight_mean * weight_mean);
	std::cout << "Sample eff.:      " << weight_eff << std::endl;
	std::cout << "Sample rel. var.: " << weight_rel_var << std::endl;
	return 0;
}

Real func_test(bubble::Point<DIM, Real> x) {
	Real result = 1;
	Real r2 = 0;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		r2 += x[dim] * x[dim];
		//result *= std::sin(PI * x[dim]);
		//result += x[dim];
	}
	result = std::exp(-10 * std::sqrt(r2)) / std::sqrt(r2);
	return result;
}

int thread_count() {
	int n = 0;
	#pragma omp parallel reduction(+:n)
	n += 1;
	return n;
}

