#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include <bubble.hpp>

using Real = double;
static bubble::Dim const DIM = 3;
static Real const PI = 3.1415926;

using Func = Real (*)(bubble::Point<DIM, Real>);
Real func_constant(bubble::Point<DIM, Real> x);
Real func_ball(bubble::Point<DIM, Real> x);
Real func_sphere(bubble::Point<DIM, Real> x);
Real func_waves(bubble::Point<DIM, Real> x);
Real func_peak(bubble::Point<DIM, Real> x);
Real func_peak_sharp(bubble::Point<DIM, Real> x);

static Func const funcs[] = {
	&func_constant,
	&func_ball,
	&func_sphere,
	&func_waves,
	&func_peak,
	&func_peak_sharp,
};

int thread_count();

void build(int func_index, char const* file_name) {
	auto builder = bubble::make_builder<DIM, Real>(funcs[func_index]);
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
	std::cout << "Scaling exp.:   " << builder.scale_exp() << std::endl;
	std::cout << "Total prime:    " << builder.prime() << std::endl;
	std::cout << "Est. mean:      " << builder.mean() << std::endl;
	std::cout << "Est. eff.:      " << est_eff << std::endl;
	std::cout << "Est. rel. var.: " << est_rel_var << std::endl;
	std::cout << "Writing to file: " << file_name << std::endl;
	builder.write(file_name);
}

void generate(int func_index, char const* file_name, std::size_t samples) {
	auto generator = bubble::make_generator<DIM, Real>(funcs[func_index]);
	std::cout << "Reading from file: " << file_name << std::endl;
	generator.read(file_name);
	std::cout << "Generating." << std::endl;
	std::vector<Real> weights;
	std::cout << "Example points:" << std::endl;
	for (std::size_t sample = 0; sample < samples; ++sample) {
		Real weight;
		bubble::Point<DIM, Real> point;
		generator.generate(&weight, &point);
		if (sample < 10) {
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
}

int main(int argc, char** argv) {
	if (argc != 4) {
		std::cout << "Usage:" << std::endl
			<< "\ttest build <func-index> <file-name>" << std::endl
			<< "\ttest generate <func-index> <file-name>" << std::endl;
		return 0;
	}
	std::string command = argv[1];
	int func_index = std::stoi(std::string(argv[2]));
	std::string file = argv[3];
	if (command == "build") {
		std::cout << "Number of threads: " << thread_count() << std::endl;
		build(func_index, file.c_str());
	} else if (command == "generate") {
		generate(func_index, file.c_str(), 1000000);
	} else {
		std::cout << "Unrecognized command." << std::endl;
	}
	return 0;
}

int thread_count() {
	int n = 0;
	#ifdef BUBBLE_USE_OPENMP
	#pragma omp parallel reduction(+:n)
	#endif
	n += 1;
	return n;
}

Real func_constant(bubble::Point<DIM, Real> x) {
	static_cast<void>(x);
	return 1;
}
Real func_ball(bubble::Point<DIM, Real> x) {
	Real r2 = 0;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		r2 += (x[dim] - 0.5) * (x[dim] - 0.5);
	}
	if (r2 < 0.4 * 0.4) {
		return 1;
	} else {
		return 0;
	}
}
Real func_sphere(bubble::Point<DIM, Real> x) {
	Real r2 = 0;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		r2 += (x[dim] - 0.5) * (x[dim] - 0.5);
	}
	if (r2 < 0.5 * 0.5 && r2 > 0.4 * 0.4) {
		return 1;
	} else {
		return 0;
	}
}
Real func_waves(bubble::Point<DIM, Real> x) {
	Real result = 1;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		Real f = std::sin(4 * PI * x[dim]);
		result *= f * f;
	}
	return result;
}
Real func_peak(bubble::Point<DIM, Real> x) {
	Real r2 = 0;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		r2 += (x[dim] - 0.5) * (x[dim] - 0.5);
	}
	return std::exp(-r2 / (2 * 0.1 * 0.1));
}
Real func_peak_sharp(bubble::Point<DIM, Real> x) {
	Real r2 = 0;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		r2 += (x[dim] - 0.5) * (x[dim] - 0.5);
	}
	return std::exp(-r2 / (2 * 0.001 * 0.001));
}

