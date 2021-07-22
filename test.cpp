#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

#include <TFile.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TRandom3.h>

#include <bubble.hpp>

using Real = double;
static bubble::Dim const DIM = 3;
static Real const PI = 3.1415926;

using Func = Real (*)(bubble::Point<DIM, Real>) noexcept;
Real func_constant(bubble::Point<DIM, Real> x) noexcept;
Real func_linear(bubble::Point<DIM, Real> x) noexcept;
Real func_simple(bubble::Point<DIM, Real> x) noexcept;
Real func_quadratic(bubble::Point<DIM, Real> x) noexcept;
Real func_ball(bubble::Point<DIM, Real> x) noexcept;
Real func_sphere(bubble::Point<DIM, Real> x) noexcept;
Real func_waves(bubble::Point<DIM, Real> x) noexcept;
Real func_peak(bubble::Point<DIM, Real> x) noexcept;
Real func_peak_sharp(bubble::Point<DIM, Real> x) noexcept;
Real func_random(bubble::Point<DIM, Real> x) noexcept;

static Func const funcs[] = {
	&func_constant,
	&func_linear,
	&func_simple,
	&func_quadratic,
	&func_ball,
	&func_sphere,
	&func_waves,
	&func_peak,
	&func_peak_sharp,
	&func_random,
};

struct Functor : public TFoamIntegrand {
	Func func;
	Functor(int index) : func(funcs[index]) { }
	virtual Double_t Density(int dim, Double_t* vec_unit) override {
		if (dim != DIM) {
			return 0;
		} else {
			bubble::Point<DIM, Real> x;
			for (bubble::Dim dim = 0; dim < DIM; ++dim) {
				x[dim] = vec_unit[dim];
			}
			return operator()(x);
		}
	}
	Real operator()(bubble::Point<DIM, Real> x) const {
		return func(x);
	}
};

int thread_count();

void build_bubble(int func_index, char const* file_name) {
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
	Real prime_err;
	Real prime = builder.prime(&prime_err);
	Real mean_err;
	Real mean = builder.mean(&mean_err);
	Real rel_var_err;
	Real rel_var = builder.rel_var(&rel_var_err);
	std::cout << "# of cells:     " << builder.tree().size() << std::endl;
	std::cout << "Total prime:    " << prime
		<< " ± " << prime_err / prime * 100. << "%" << std::endl;
	std::cout << "Est. mean:      " << mean
		<< " ± " << mean_err / mean * 100. << "%" << std::endl;
	std::cout << "Est. rel. var.: " << rel_var
		<< " ± " << rel_var_err / rel_var * 100. << "%" << std::endl;
	std::cout << "Writing to file: " << file_name << std::endl;
	builder.write(file_name);
}

template<typename T>
void analyze_weights(Real prime, std::vector<T> const& weights) {
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
	std::cout << "Integral est.:    " << weight_mean * prime
		<< " ± " << std::sqrt(weight_sq / weights.size()) << std::endl;
	std::cout << "Sample eff.:      " << weight_eff << std::endl;
	std::cout << "Sample rel. var.: " << weight_rel_var << std::endl;
}

void generate_bubble(int func_index, char const* file_name, std::size_t samples) {
	auto generator = bubble::make_generator<DIM, Real>(funcs[func_index]);
	std::cout << "Reading from file: " << file_name << std::endl;
	generator.read(file_name);
	std::vector<Real> weights;
	std::cout << "Generating." << std::endl;
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
	analyze_weights(generator.prime(), weights);
}

void build_foam(int func_index, char const* file_name) {
	TFile file(file_name, "RECREATE");
	TFoam* foam = new TFoam("foam");
	TRandom3* random = new TRandom3();
	Functor func(func_index);
	//foam->SetChat(0);
	foam->SetkDim(DIM);
	foam->SetRho(&func);
	foam->SetPseRan(random);
	foam->SetnSampl(1024);
	foam->SetnCells(500000);
	foam->SetOptRej(0);
	foam->SetOptDrive(1);
	foam->Initialize();
	foam->Write("foam");
}

void generate_foam(int func_index, char const* file_name, std::size_t samples) {
	std::cout << "Reading from file: " << file_name << std::endl;
	TFile file(file_name);
	TFoam* foam = file.Get<TFoam>("foam");
	if (foam == nullptr) {
		std::cout << "Couldn't find foam." << std::endl;
		return;
	}
	TRandom3 random;
	Functor func(func_index);
	foam->ResetRho(&func);
	foam->SetPseRan(&random);
	std::vector<Double_t> weights;
	std::cout << "Generating." << std::endl;
	for (std::size_t sample = 0; sample < samples; ++sample) {
		Double_t weight;
		Double_t point[DIM];
		weight = foam->MCgenerate(&point[0]);
		if (sample < 10) {
			std::cout
				<< "\tw: " << weight << ", x: ("
				<< point[0] << ", "
				<< point[1] << ", "
				<< point[2] << ")" << std::endl;
		}
		weights.push_back(weight);
	}
	analyze_weights(foam->GetPrimary(), weights);
}

int main(int argc, char** argv) {
	if (argc != 5) {
		std::cout << "Usage:" << std::endl
			<< "\ttest (bubble,foam) build <file-name> <func-index>" << std::endl
			<< "\ttest (bubble,foam) generate <file-name> <func-index>" << std::endl;
		return 0;
	}
	std::string type = argv[1];
	std::string command = argv[2];
	std::string file = argv[3];
	int func_index = std::stoi(std::string(argv[4]));
	if (command == "build" && type == "bubble") {
		std::cout << "Number of threads: " << thread_count() << std::endl;
		build_bubble(func_index, file.c_str());
	} else if (command == "build" && type == "foam") {
		build_foam(func_index, file.c_str());
	} else if (command == "generate" && type == "bubble") {
		generate_bubble(func_index, file.c_str(), 1000000);
	} else if (command == "generate" && type == "foam") {
		generate_foam(func_index, file.c_str(), 1000000);
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

Real func_constant(bubble::Point<DIM, Real> x) noexcept {
	static_cast<void>(x);
	return 1;
}
Real func_simple(bubble::Point<DIM, Real> x) noexcept {
	Real result = 1;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		Real f = std::sin(PI * x[dim]);
		result *= f;
	}
	return result;
}
Real func_linear(bubble::Point<DIM, Real> x) noexcept {
	return x[0];
}
Real func_quadratic(bubble::Point<DIM, Real> x) noexcept {
	Real result = 0;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		result += x[dim] * x[dim];
	}
	return result;
}
Real func_ball(bubble::Point<DIM, Real> x) noexcept {
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
Real func_sphere(bubble::Point<DIM, Real> x) noexcept {
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
Real func_waves(bubble::Point<DIM, Real> x) noexcept {
	Real result = 1;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		Real f = std::sin(4 * PI * x[dim]);
		result *= f * f;
	}
	return result;
}
Real func_peak(bubble::Point<DIM, Real> x) noexcept {
	Real r2 = 0;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		r2 += (x[dim] - 0.5) * (x[dim] - 0.5);
	}
	return std::exp(-r2 / (2 * 0.1 * 0.1));
}
Real func_peak_sharp(bubble::Point<DIM, Real> x) noexcept {
	Real r2 = 0;
	for (bubble::Dim dim = 0; dim < DIM; ++dim) {
		r2 += (x[dim] - 0.5) * (x[dim] - 0.5);
	}
	return std::exp(-r2 / (2 * 0.001 * 0.001));
}
Real func_random(bubble::Point<DIM, Real> x) noexcept {
	// Use x as a seed for a random number generator.
	unsigned char const* seed = reinterpret_cast<unsigned char const*>(&x);
	std::size_t seed_len = sizeof(x);
	std::seed_seq seed_seq(seed, seed + seed_len);
	std::minstd_rand rnd(seed_seq);
	std::uniform_real_distribution<Real> dist(0., 1.);
	return dist(rnd);
}

