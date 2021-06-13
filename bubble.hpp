#ifndef BUBBLE_HPP__
#define BUBBLE_HPP__

#include <array>
#include <cmath>
#include <fstream>
#include <ios>
#include <limits>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace bubble {

// Some utility functions.
template<typename T>
static T sq(T x) {
	return x * x;
}

using Seed = unsigned long;
using Dim = unsigned char;
using CellHandle = std::size_t;

template<Dim D, typename R>
using Point = std::array<R, D>;

// A branch node for a k-d type tree.
template<Dim D, typename R>
struct KdBranch final {
	using Real = R;
	static Dim const DIM = D;

	static std::size_t const NUM_CHILDREN = 2;
	CellHandle children[NUM_CHILDREN];
	R split;
	Dim dim;

	KdBranch(Dim dim, R split) :
		children(),
		split(split),
		dim(dim) { }

	R volume(std::size_t child) const {
		if (child == 0) {
			return split;
		} else {
			return 1 - split;
		}
	}
	Point<D, R> extent(std::size_t child) const {
		Point<D, R> result;
		result.fill(1);
		if (child == 0) {
			result[dim] = split;
		} else {
			result[dim] = 1 - split;
		}
		return result;
	}
	Point<D, R> offset(std::size_t child) const {
		Point<D, R> result;
		result.fill(0);
		if (child == 0) {
		} else {
			result[dim] = split;
		}
		return result;
	}
};

// Cells can be either branches or leafs.
enum class CellType : unsigned char {
	BRANCH,
	LEAF,
};

// A spatial indexing tree of cells. Cells can be divided only, not merged, so
// the tree will only get larger with time.
template<typename B, typename DB, typename DL=DB>
class Tree final {
	template<typename B1, typename DB1, typename DL1>
	friend class Tree;

	// Represents a cell in a tree with branch type `B`, storing data of type
	// `DB` each branch node, and data of type `DL` on each leaf node.
	struct Cell final {
	private:
		CellType _type;
	public:
		union {
			B branch;
			char leaf;
		};
		union {
			DB branch_data;
			DL leaf_data;
		};

		CellType type() const {
			return _type;
		}

		// Basic constructors.
		template<typename... Ts>
		Cell(DB const& data, Ts... args) :
			_type(CellType::BRANCH),
			branch(args...),
			branch_data(data) { }
		Cell(DB const& data, B const& branch) :
			_type(CellType::BRANCH),
			branch(branch),
			branch_data(data) { }
		Cell(DL const& data) :
			_type(CellType::LEAF),
			leaf_data(data) { }
		// Copy constructor.
		Cell(Cell const& other) : _type(other._type) {
			if (_type == CellType::BRANCH) {
				new(&branch) B(other.branch);
				new(&branch_data) DB(other.branch_data);
			} else {
				new(&leaf_data) DL(other.leaf_data);
			}
		}
		// Destructor.
		~Cell() {
			if (_type == CellType::BRANCH) {
				branch.~B();
				branch_data.~DB();
			} else {
				leaf_data.~DL();
			}
		}
		// Assignment.
		Cell& operator=(Cell const& other) {
			if (_type == other._type) {
				if (_type == CellType::BRANCH) {
					branch = other.branch;
					branch_data = other.branch_data;
				} else {
					leaf_data = other.leaf_data;
				}
			} else {
				if (_type == CellType::BRANCH) {
					branch.~B();
					branch_data.~DB();
				} else {
					leaf_data.~DL();
				}
				if (other._type == CellType::BRANCH) {
					_type = CellType::BRANCH;
					new(&branch) B(other.branch);
					new(&branch_data) DB(other.branch_data);
				} else {
					_type = CellType::LEAF;
					new(&leaf_data) DL(other.leaf_data);
				}
			}
			return *this;
		}
	};

	std::vector<Cell> _cells;

	Tree(std::vector<Cell>&& cells) : _cells(cells) { }

public:
	using Branch = B;
	using BranchData = DB;
	using LeafData = DL;
	static Dim const DIM = Branch::DIM;
	using Real = typename Branch::Real;
	static std::size_t const NUM_CHILDREN = Branch::NUM_CHILDREN;

	static_assert(
		std::is_floating_point<Real>::value,
		"Branch::Real must be a floating point type.");

	Tree(LeafData root_data) : _cells{Cell(root_data)} { }

	CellHandle root() const {
		return 0;
	}
	std::size_t size() const {
		return _cells.size();
	}
	CellHandle get(std::size_t idx) const {
		return idx;
	}
	CellHandle child(CellHandle parent, std::size_t child_idx) const {
		return _cells[parent].branch.children[child_idx];
	}

	CellType type(CellHandle cell) const {
		return _cells[cell].type();
	}
	LeafData const& leaf_data(CellHandle cell) const {
		return _cells[cell].leaf_data;
	}
	LeafData& leaf_data(CellHandle cell) {
		return _cells[cell].leaf_data;
	}
	BranchData const& branch_data(CellHandle cell) const {
		return _cells[cell].branch_data;
	}
	BranchData& branch_data(CellHandle cell) {
		return _cells[cell].branch_data;
	}

	Point<DIM, Real> extent_local(CellHandle parent, std::size_t child) const {
		return _cells[parent].branch.extent(child);
	}
	Point<DIM, Real> offset_local(CellHandle parent, std::size_t child) const {
		return _cells[parent].branch.offset(child);
	}

	// Splits a leaf parent cell into new children.
	void split(
			CellHandle cell,
			Branch branch,
			BranchData branch_data,
			const LeafData (&leaf_data)[NUM_CHILDREN]) {
		// TODO: Get rid of check?
		if (_cells[cell].type() != CellType::LEAF) {
			throw std::runtime_error("Can only split a leaf cell.");
		}

		// Append the new cells to the tree.
		_cells.reserve(_cells.size() + NUM_CHILDREN);
		Cell parent(branch_data, branch);
		for (std::size_t child = 0; child < NUM_CHILDREN; ++child) {
			parent.branch.children[child] = _cells.size();
			_cells.push_back(Cell(leaf_data[child]));
		}
		_cells[cell] = parent;
	}

	template<typename DBP, typename DLP, typename FB, typename FL>
	Tree<B, DBP, DLP> transform(FB fb, FL fl) const {
		using CellOther = typename Tree<B, DBP, DLP>::Cell;
		std::vector<CellOther> cells_new;
		cells_new.reserve(_cells.size());
		for (Cell const& cell : _cells) {
			if (cell.type() == CellType::BRANCH) {
				cells_new.push_back(CellOther(fb(cell.branch_data), cell.branch));
			} else {
				cells_new.push_back(CellOther(fl(cell.leaf_data)));
			}
		}
		Tree<B, DBP, DLP> result(std::move(cells_new));
		return result;
	}
};

// Stores statistics about function evaluation in a region of space.
template<typename FI>
struct Stats final {
	FI mean_total;
	FI mom2_total;
	FI mom4_total;
	std::size_t count;

	Stats() :
		mean_total(0),
		mom2_total(0),
		mom4_total(0),
		count(0) { }
	Stats(
		FI mean_total,
		FI mom2_total, FI mom4_total,
		std::size_t count) :
		mean_total(mean_total),
		mom2_total(mom2_total),
		mom4_total(mom4_total),
		count(count) { }

	void update(FI measure) {
		FI measure_sq = sq(measure);
		mean_total += measure;
		mom2_total += measure_sq;
		mom4_total += sq(measure_sq);
		count += 1;
	}

	Stats<FI>& operator+=(Stats<FI> const& rhs) {
		mean_total += rhs.mean_total;
		mom2_total += rhs.mom2_total;
		mom4_total += rhs.mom4_total;
		count += rhs.count;
		return *this;
	}
	Stats<FI>& operator-=(Stats<FI> const& rhs) {
		mean_total -= rhs.mean_total;
		mom2_total -= rhs.mom2_total;
		mom4_total -= rhs.mom4_total;
		count -= rhs.count;
		return *this;
	}
	Stats<FI>& operator*=(FI scale) {
		FI scale_sq = sq(scale);
		mean_total *= scale;
		mom2_total *= scale_sq;
		mom4_total *= sq(scale_sq);
		return *this;
	}

	FI est_prime(FI* prime_var_out=nullptr) const {
		// Unbiased estimate of the second moment.
		FI mom2 = mom2_total / count;
		// Unbiased estimate of variance of square.
		FI var2 = (mom4_total - sq(mom2) * count) / (count - 1);
		// Variance in estimate of second moment.
		FI mom2_var = var2 / count;
		// Estimate the ideal prime value, which is given by:
		//     F_c = \sqrt{V_c \int dx f^2(x)}
		// with a small correction for bias added on.
		FI prime = std::sqrt(mom2) * (1 + mom2_var / (8 * mom2));
		// Variance in estimate of prime value.
		if (prime_var_out != nullptr) {
			*prime_var_out = mom2_var / (4 * mom2);
		}
		return prime;
	}
	FI est_eff(FI* eff_var_out=nullptr) const {
		// Unbiased estimate of mean.
		FI mean = mean_total / count;
		// Unbiased estimate of variance.
		FI var = (mom2_total - sq(mean) * count) / (count - 1);
		// Variance in estimate of mean.
		FI mean_var = var / count;
		// Prime estimate and variance.
		FI prime_var;
		FI prime = est_prime(&prime_var);
		// Estimate the efficiency.
		// TODO: Improve the quality of this estimator. Right now, it
		// may be quite biased in certain circumstances. Also improve
		// the quality of the variance estimate to go with it.
		FI eff = mean / prime;
		// Variance in estimate of efficiency.
		// TODO: Fix this calculation, since the two errors are not
		// independent. They are positively correlated, so this is an
		// overestimate of the true error (I think).
		if (eff_var_out != nullptr) {
			*eff_var_out = sq(eff) * (mean_var / sq(mean) + prime_var / sq(prime));
		}
		return eff;
	}
	FI est_vol(FI prime_tot) const {
		FI prime_var;
		FI prime = est_prime(&prime_var);
		FI mom2 = mom2_total / count;
		FI var2 = (mom4_total - sq(mom2) * count) / (count - 1);
		FI mom2_var = var2 / count;
		FI vol = 0.5 * std::sqrt((prime_tot - prime) / prime * (var2 / sq(mom2)));
		return vol;
	}
};
template<typename FI>
inline Stats<FI> operator+(Stats<FI> lhs, Stats<FI> const& rhs) {
	lhs += rhs;
	return lhs;
}
template<typename FI>
inline Stats<FI> operator-(Stats<FI> lhs, Stats<FI> const& rhs) {
	lhs -= rhs;
	return lhs;
}
template<typename FI>
inline Stats<FI> operator*(FI lhs, Stats<FI> rhs) {
	rhs *= lhs;
	return rhs;
}
template<typename FI>
inline Stats<FI> operator*(Stats<FI> lhs, FI rhs) {
	lhs *= rhs;
	return lhs;
}

template<Dim D, typename R, typename F>
class CellBuilder final {
	template<Dim D1, typename R1, typename F1>
	friend class CellGenerator;

	// Left and right random engines.
	using RndL = std::mt19937;
	using RndR = std::mt19937_64;

	// Return type of the function.
	using FR = decltype(std::declval<F>()(std::declval<Point<D, R> >()));
	// Type of integrals over function.
	using FI = decltype(std::declval<FR>() * std::declval<R>());
	static_assert(
		std::is_floating_point<R>::value,
		"R must be a floating point type.");
	static_assert(
		std::is_floating_point<FR>::value,
		"Function must return floating point type.");

	struct BranchData final {
		FI prime;
		FI vol;
		BranchData() : vol(0) { }
	};
	struct LeafData final {
		Stats<FI> stats;
		LeafData() : stats() { }
		LeafData(Stats<FI> const& stats) : stats(stats) { }
	};
	using TreeType = Tree<KdBranch<D, R>, BranchData, LeafData>;

	// Function.
	F _func;
	// Spatial tree storing the statistics.
	TreeType _tree;
	// Random number engine.
	RndL _rnd_left;

	// The target efficiency must be reached in all cells to completely
	// initialize the tree.
	FI _target_eff;
	FI _target_eff_rel_err;
	// Chi-squared needed for division.
	FI _chi_squared_for_divide;
	// Number of bins in edge histograms.
	std::size_t _hist_num_bins;
	// Maximum number of samples to be taken in a single cell in the exploration
	// phase.
	std::size_t _max_num_samples_explore;
	// The quality to tune the cell generator to.
	FI _tune_quality;
	// How many stages to tune. Usually a handful is good enough.
	std::size_t _tune_num_stages;
	// Maximum number of samples to be taken in a single cell in the tuning
	// phase.
	std::size_t _max_num_samples_tune;

	// Gets a seed from the left generator.
	std::vector<Seed> sample_rnd_left() {
		std::uniform_int_distribution<Seed> dist;
		std::vector<Seed> result {
			dist(_rnd_left),
			dist(_rnd_left),
			dist(_rnd_left),
			dist(_rnd_left),
		};
		return result;
	}

	// Fills in primes for a cell and all of its descendents.
	FI fill_primes(CellHandle cell) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).prime = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child);
				_tree.branch_data(parent).prime += fill_primes(child);
			}
			return _tree.branch_data(parent).prime;
		} else {
			return _tree.leaf_data(parent).stats.est_prime();
		}
	}
	// Fills in volatility for a cell and all of its descendents.
	FI fill_vols(CellHandle cell, FI prime_total) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).vol = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child);
				_tree.branch_data(parent).vol += fill_vols(child, prime_total);
			}
			return _tree.branch_data(parent).vol;
		} else {
			return _tree.leaf_data(parent).stats.est_vol(prime_total);
		}
	}

	// Generates points within a cell endlessly until one of three termination
	// conditions:
	//  * The efficiency reaches the target efficiency.
	//  * The minimum possible efficiency is determined to be worse than the
	//    target efficiency, so a cell division is requested.
	//  * The allotted number of function evaluations is exceeded.
	void explore_cell(
			CellHandle cell,
			Point<D, R> offset, Point<D, R> extent,
			std::mutex& tree_mutex,
			std::mutex& rnd_mutex) {
		// If the cell has children, then recursively explore those children.
		CellType type;
		{
			std::lock_guard<std::mutex> lock(tree_mutex);
			type = _tree.type(cell);
		}
		if (type == CellType::BRANCH) {
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				#pragma omp task
				{
					Point<D, R> offset_local;
					Point<D, R> extent_local;
					CellHandle child;
					{
						std::lock_guard<std::mutex> lock(tree_mutex);
						offset_local = _tree.offset_local(cell, child_idx);
						extent_local = _tree.extent_local(cell, child_idx);
						child = _tree.child(cell, child_idx);
					}
					Point<D, R> offset_child;
					Point<D, R> extent_child;
					for (Dim dim = 0; dim < D; ++dim) {
						offset_child[dim] = offset[dim] + extent[dim] * offset_local[dim];
						extent_child[dim] = extent[dim] * extent_local[dim];
					}
					explore_cell(
						child,
						offset_child, extent_child,
						tree_mutex, rnd_mutex);
				}
			}
			#pragma omp taskwait
			return;
		} else {
			R volume = 1;
			for (Dim dim = 0; dim < D; ++dim) {
				volume *= extent[dim];
			}
			// Statistics measures.
			Stats<FI> stats;
			Stats<FI> stats_init;
			{
				std::lock_guard<std::mutex> lock(tree_mutex);
				stats_init = _tree.leaf_data(cell).stats;
			}
			// Histograms, one for each axis of the hypercube.
			std::vector<Stats<FI> > hist_stats[D];
			for (Dim dim = 0; dim < D; ++dim) {
				hist_stats[dim].resize(_hist_num_bins);
			}

			// Create a random number generator (a "right" generator) which will
			// be used.
			std::vector<Seed> seed;
			{
				std::lock_guard<std::mutex> lock(rnd_mutex);
				seed = sample_rnd_left();
			}
			std::seed_seq seed_seq(seed.begin(), seed.end());
			RndR rnd(seed_seq);

			// The next number at which termination conditions will be checked.
			std::size_t next_count = 128;
			while (true) {
				for (std::size_t sample = stats.count; sample < next_count; ++sample) {
					// Generate a point randomly within the cell and evaluate
					// the function. First choose a histogram bin to generate
					// in, then produce the point itself.
					Point<D, R> point;
					std::size_t hist_idx[D];
					std::uniform_int_distribution<std::size_t> idx_dist(
						0,
						_hist_num_bins - 1);
					for (Dim dim = 0; dim < D; ++dim) {
						std::size_t idx = idx_dist(rnd);
						hist_idx[dim] = idx;
						R lower = offset[dim]
							+ extent[dim] * idx / _hist_num_bins;
						R upper = offset[dim]
							+ extent[dim] * (idx + 1) / _hist_num_bins;
						std::uniform_real_distribution<R> x_dist(lower, upper);
						point[dim] = x_dist(rnd);
					}
					// First choose a histogram index to generate in, then make
					// the point itself.
					// TODO: Consider computing these with a more numerically
					// accurate method, such as an averaging tree.
					FI f = volume * _func(point);
					stats.update(f);
					// Fill histograms.
					for (Dim dim = 0; dim < D; ++dim) {
						hist_stats[dim][hist_idx[dim]].update(f);
					}
				}

				// Check termination conditions. The termination conditions are
				// expensive, so scale by two so that they only get checked
				// every order of magnitude or so.
				next_count *= 2;

				// Get efficiency and prime.
				Stats<FI> stats_total = stats + stats_init;
				FI prime_var;
				FI prime = stats_total.est_prime(&prime_var);
				FI eff_var;
				FI eff = stats_total.est_eff(&eff_var);

				// Termination condition 1. Efficiency meets criteria.
				// There are two conditions for the efficiency to be good
				// enough:
				//  * Efficiency better than target efficiency.
				//  * Uncertainty in efficiency smaller than target uncertainty.
				FI eff_err = std::sqrt(eff_var);
				if (eff > _target_eff && eff_err / eff < _target_eff_rel_err) {
					std::lock_guard<std::mutex> lock(tree_mutex);
					_tree.leaf_data(cell).stats += stats;
					return;
				}

				// Termination condition 2. Cell division.
				// For cell division to occur, the edge histograms must be
				// different enough (by chi-squared test) from constant. Then,
				// the division site that minimizes the new sum of primes is
				// chosen.
				FI chi_squared = 0;
				FI prime_min = prime;
				Dim dim_min = 0;
				std::size_t bin_min = _hist_num_bins / 2;
				Stats<FI> stats_lower_min;
				Stats<FI> stats_upper_min;
				bool enough_samples = true;
				for (Dim dim = 0; dim < D; ++dim) {
					// Integrated statistics below and above the division site.
					Stats<FI> stats_lower = Stats<FI>();
					Stats<FI> stats_upper = stats;
					// Test division at each bin boundary. Note that the last
					// bin is left off: such a data point would have
					// `count_upper == 0`, which is difficult to deal with.
					for (std::size_t bin = 0; bin < _hist_num_bins - 1; ++bin) {
						// Update the integrated statistics.
						stats_lower += hist_stats[dim][bin];
						stats_upper -= hist_stats[dim][bin];
						// TODO: Have a better reason for the number of
						// statistics required per bin. Also, put a check
						// outside that the average number of counts per bin is
						// large enough, before starting to loop through bins.
						if (stats_lower.count <= 3 || stats_upper.count <= 3) {
							enough_samples = false;
							goto not_enough_samples;
						}
						// Estimate the total prime resulting from splitting the
						// cell at this bin boundary.
						// Need to scale by the volume fractions, to account for
						// volume being included in the statistics measures.
						R lambda_1 = R(bin + 1) / _hist_num_bins;
						R lambda_2 = R(_hist_num_bins - bin - 1) / _hist_num_bins;
						FI prime_var_lower;
						FI prime_lower = (lambda_1 * stats_lower).est_prime(
							&prime_var_lower);
						FI prime_var_upper;
						FI prime_upper = (lambda_2 * stats_upper).est_prime(
							&prime_var_upper);
						FI prime_new = prime_lower + prime_upper;
						// Because the two prime estimates are independent, the
						// variances can be summed.
						FI prime_var_new = prime_var_lower + prime_var_upper;
						// Chi squared test against the maximum possible value
						// for sum of primes, which is the current prime of the
						// cell.
						// TODO: The `prime_var_new` and `prime_var` are not
						// independent, and so technically shouldn't be added.
						// Not sure if there's enough dependence for it to
						// matter though. That would be good to work out.
						// TODO: Using chi-squared test here even though prime
						// is not normally distributed. At best it is a square
						// root of a normal variable (from central limit
						// theorem, which doesn't necessary work that well for
						// all distributions). It's not clear what a better test
						// would be though, so leave it for now.
						chi_squared += sq(prime_new - prime)
							/ (prime_var_new + prime_var);
						// Update minimum.
						if (prime_new < prime_min && prime_new < prime) {
							prime_min = prime_new;
							dim_min = dim;
							bin_min = bin;
							stats_lower_min = stats_lower;
							stats_upper_min = stats_upper;
						}
					}
				}
				not_enough_samples:
				if (enough_samples) {
					// Average the chi-squared over every degree of freedom.
					chi_squared /= D * _hist_num_bins;
					if (chi_squared > _chi_squared_for_divide) {
						std::lock_guard<std::mutex> lock(tree_mutex);
						R lambda_1 = R(bin_min + 1) / _hist_num_bins;
						R lambda_2 = R(_hist_num_bins - bin_min - 1) / _hist_num_bins;
						R split = offset[dim_min] + lambda_1 * extent[dim_min];
						LeafData leaf_data[2] = {
							{ lambda_1 * stats_lower_min },
							{ lambda_2 * stats_upper_min },
						};
						_tree.split(
							cell,
							typename TreeType::Branch(dim_min, split),
							BranchData(),
							leaf_data);
						// Recursively calling will cause the new children to be
						// explored.
						explore_cell(
							cell,
							offset, extent,
							tree_mutex, rnd_mutex);
						return;
					}
				}

				// Termination condition 3. Oversampled.
				// If too much time has been spent sampling from this cell and
				// a good division site has not been found, then just give up.
				// No more exploration will be done in this cell.
				if (stats.count >= _max_num_samples_explore) {
					std::lock_guard<std::mutex> lock(tree_mutex);
					_tree.leaf_data(cell).stats += stats;
					return;
				}
			}
		}
	}

	// Generates a number of points in a cell. The points are distributed among
	// the cell's children according to the "volatility" of the cell, which can
	// be estimated using `Stats::est_vol`.
	void tune_cell(
			CellHandle cell,
			Point<D, R> offset, Point<D, R> extent,
			FI prime_total,
			std::size_t count,
			std::mutex& rnd_mutex) {
		// No tree mutex is needed for tuning, since only one thread will ever
		// be writing to a single cell at a time.
		if (_tree.type(cell) == CellType::LEAF) {
			// For a leaf, there are no children, so all of the samples will be
			// uniformly distributed.
			R volume = 1;
			for (Dim dim = 0; dim < D; ++dim) {
				volume *= extent[dim];
			}
			Stats<FI> stats;

			// Create a random number generator.
			std::vector<Seed> seed;
			{
				std::lock_guard<std::mutex> lock(rnd_mutex);
				seed = sample_rnd_left();
			}
			std::seed_seq seed_seq(seed.begin(), seed.end());
			RndR rnd(seed_seq);

			for (std::size_t sample = 0; sample < count; ++sample) {
				Point<D, R> point;
				for (Dim dim = 0; dim < D; ++dim) {
					R lower = offset[dim];
					R upper = offset[dim] + extent[dim];
					std::uniform_real_distribution<R> x_dist(lower, upper);
					point[dim] = x_dist(rnd);
				}
				stats.update(volume * _func(point));
			}
			_tree.leaf_data(cell).stats += stats;
			return;
		} else {
			// For a branch, estimate the volatility of each child. The counts
			// will be divided up according to that.
			FI vols[TreeType::NUM_CHILDREN];
			FI vol_total = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(cell, child_idx);
				if (_tree.type(child) == CellType::BRANCH) {
					vols[child_idx] = _tree.branch_data(child).vol;
				} else {
					vols[child_idx] = _tree.leaf_data(child).stats.est_vol(prime_total);
				}
				vol_total += vols[child_idx];
			}
			std::size_t counts_child[TreeType::NUM_CHILDREN];
			std::size_t count_total = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				FI vol_frac = vols[child_idx] / vol_total;
				counts_child[child_idx] = std::size_t(vol_frac * count);
				count_total += counts_child[child_idx];
			}
			while (count_total < count) {
				counts_child[count_total % TreeType::NUM_CHILDREN] += 1;
				count_total += 1;
			}
			while (count_total > count) {
				counts_child[count_total % TreeType::NUM_CHILDREN] -= 1;
				count_total -= 1;
			}

			// Sample children according to volatility.
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				#pragma omp task
				{
					Point<D, R> offset_local;
					Point<D, R> extent_local;
					CellHandle child;
					offset_local = _tree.offset_local(cell, child_idx);
					extent_local = _tree.extent_local(cell, child_idx);
					child = _tree.child(cell, child_idx);
					Point<D, R> offset_child;
					Point<D, R> extent_child;
					for (Dim dim = 0; dim < D; ++dim) {
						offset_child[dim] = offset[dim] + extent[dim] * offset_local[dim];
						extent_child[dim] = extent[dim] * extent_local[dim];
					}
					tune_cell(
						child,
						offset_child, extent_child,
						prime_total,
						counts_child[child_idx],
						rnd_mutex);
				}
			}
			#pragma omp taskawait
			return;
		}
	}

public:
	CellBuilder(F func, Seed seed=std::random_device()()) :
			_func(func),
			_tree(LeafData()),
			_rnd_left(),
			// TODO: Properly initialize these parameters.
			_target_eff(0.9),
			_target_eff_rel_err(0.01),
			_chi_squared_for_divide(10.),
			_hist_num_bins(100),
			_max_num_samples_explore(1000000),
			_tune_quality(0.95),
			_tune_num_stages(4),
			_max_num_samples_tune(1000000) {
		std::seed_seq seed_seq { seed };
		_rnd_left.seed(seed_seq);
	}

	void explore() {
		Point<D, R> offset;
		Point<D, R> extent;
		offset.fill(0);
		extent.fill(1);
		CellHandle root = _tree.root();
		std::mutex tree_mutex;
		std::mutex rnd_mutex;
		#pragma omp parallel
		#pragma omp single
		explore_cell(root, offset, extent, tree_mutex, rnd_mutex);
		FI prime_total = fill_primes(root);
		fill_vols(root, prime_total);
	}

	// Tunes the cell generator to within some fraction of the ideal efficiency.
	void tune() {
		Point<D, R> offset;
		Point<D, R> extent;
		offset.fill(0);
		extent.fill(1);
		CellHandle root = _tree.root();
		// The tuning process is made somewhat complicated by the staged
		// process. Since the accuracy of the volatilities might not be good if
		// the cell generator is not well tuned, the tuning is done in stages
		// with the volatilities recalculated between each stage.
		for (std::size_t stage = 0; stage < _tune_num_stages; ++stage) {
			// Find total primes and volatilities.
			FI prime_total = fill_primes(root);
			FI vol_total = fill_vols(root, prime_total);
			// Find the number of samples needed to reach the desired accuracy.
			FI quality = 1 - (1 - _tune_quality) * stage / _tune_num_stages;
			std::size_t count =
				std::size_t(1 / (2 * (1 - quality)) * sq(vol_total));
			if (count > _max_num_samples_tune) {
				count = _max_num_samples_tune;
			}
			// Sample from the cell generator.
			std::mutex rnd_mutex;
			#pragma omp parallel
			#pragma omp single
			tune_cell(root, offset, extent, prime_total, count, rnd_mutex);
		}
		FI prime_total = fill_primes(root);
		fill_vols(root, prime_total);
	}

	FI prime() const {
		return _tree.branch_data(_tree.root()).prime;
	}
};

template<Dim D, typename R, typename F>
class CellGenerator final {
	using Rnd = std::mt19937_64;

	// Return type of the function.
	using FR = decltype(std::declval<F>()(std::declval<Point<D, R> >()));
	// Type of integrals over function.
	using FI = decltype(std::declval<FR>() * std::declval<R>());

	// The cells store the prime value.
	struct Data final {
		FI prime;
		Data() : prime(0) { }
		Data(FI prime) : prime(prime) { }
	};
	using TreeType = Tree<KdBranch<D, R>, Data>;
	using Builder = CellBuilder<D, R, F>;

	// Function.
	F _func;
	// Spatial tree.
	TreeType _tree;
	// Random number engine.
	Rnd _rnd;

	// Fills in the prime value for a cell and all of its descendents,
	// recursively. This is done on load so that parent prime values don't have
	// to be stored on disk. Returns the total prime value on completion.
	FI fill_primes(CellHandle cell) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).prime = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child);
				_tree.branch_data(parent).prime += fill_primes(child);
			}
			return _tree.branch_data(parent).prime;
		} else {
			return _tree.leaf_data(parent).prime;
		}
	}

	void select(
			CellHandle cell,
			R select_cell, Point<D, R> select_point,
			Point<D, R> offset, Point<D, R> extent,
			FR* weight_out, Point<D, R>* point_out) const {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				// Each child has a chance of being selected proportional to its
				// prime value.
				FI prime_ratio = _tree.leaf_data(child).prime / _tree.branch_data(parent).prime;
				if (select_cell <= prime_ratio || child_idx == TreeType::NUM_CHILDREN - 1) {
					// Scale `select_cell` to be between 0 and 1 again.
					select_cell /= prime_ratio;
					// Adjust offset and extent for the next level of recursion.
					Point<D, R> offset_local = _tree.offset_local(parent, child_idx);
					Point<D, R> extent_local = _tree.extent_local(parent, child_idx);
					for (Dim dim = 0; dim < D; ++dim) {
						offset[dim] += extent[dim] * offset_local[dim];
						extent[dim] *= extent_local[dim];
					}
					// Select the point from the child cell.
					select(
						child,
						select_cell, select_point,
						offset, extent,
						weight_out, point_out);
					return;
				} else {
					select_cell -= prime_ratio;
				}
			}
		} else if (_tree.type(parent) == CellType::LEAF) {
			// Rescale the function evaluation point.
			for (Dim dim = 0; dim < D; ++dim) {
				(*point_out)[dim] = offset[dim] + extent[dim] * select_point[dim];
			}
			*weight_out = _func(*point_out) / _tree.leaf_data(parent).prime;
		}
	}
	void select(
			CellHandle cell,
			R select_cell, Point<D, R> select_point,
			FR* weight_out, Point<D, R>* point_out) const {
		Point<D, R> offset;
		Point<D, R> extent;
		offset.fill(0);
		extent.fill(1);
		return select(
			cell,
			select_cell, select_point,
			offset, extent,
			weight_out, point_out);
	}

public:
	// Creates a `CellGenerator` directly from a `CellBuilder`. This should only
	// be used for testing purposes on small generators. For large ones, this
	// can lead to running out of memory very quickly.
	CellGenerator(
			Builder const& init,
			Seed seed=std::random_device()()) :
			_func(init._func),
			_tree(init._tree.template transform<Data, Data>(
				[](typename Builder::BranchData data) {
					return Data { data.prime };
				},
				[](typename Builder::LeafData data) {
					return Data { data.stats.est_prime() };
				})),
			_rnd() {
		std::seed_seq seed_seq { seed };
		_rnd.seed(seed_seq);
	}
	CellGenerator(
			char const* file_name,
			F func,
			Seed seed=std::random_device()()) :
			_func(func),
			_tree(Data{ 1 }),
			_rnd() {
		std::seed_seq seed_seq { seed };
		_rnd.seed(seed_seq);
		std::ifstream fin(file_name,
			std::ios_base::in
			| std::ios_base::binary);
		fill_primes(_tree.root());
	}

	void generate(FR* weight_out, Point<D, R>* point_out) {
		// Generate random number.
		std::uniform_real_distribution<R> dist(0, 1);
		R select_cell = dist(_rnd);
		Point<D, R> select_point;
		for (Dim dim = 0; dim < D; ++dim) {
			select_point[dim] = dist(_rnd);
		}
		// Draw from the root cell.
		select(_tree.root(), select_cell, select_point, weight_out, point_out);
	}
};

// Convenience functions to allow for template parameter inference.
template<Dim D, typename R, typename F>
inline CellBuilder<D, R, F> make_builder(
		F func,
		Seed seed=std::random_device()()) {
	return CellBuilder<D, R, F>(func, seed);
}
template<Dim D, typename R, typename F>
inline CellGenerator<D, R, F> make_generator(
		CellBuilder<D, R, F> const& builder,
		Seed seed=std::random_device()()) {
	return CellGenerator<D, R, F>(builder, seed);
}

}

#endif

