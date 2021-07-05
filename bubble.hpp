#ifndef BUBBLE_HPP__
#define BUBBLE_HPP__

#include <array>
#include <cmath>
#include <fstream>
#include <ios>
#include <limits>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace bubble {

// Identifier used in files.
static char const FILE_HEAD[] = "bubble file";

// Some utility functions.
template<typename T>
static T sq(T x) {
	return x * x;
}
template<typename T>
static T clamp(T x, T min, T max) {
	if (x < min) {
		return min;
	} else if (x < max) {
		return x;
	} else {
		return max;
	}
}
template<typename T>
static T clamp_above(T x, T min=0) {
	if (!(x >= min)) {
		return min;
	} else {
		return x;
	}
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
	static_assert(
		std::numeric_limits<Real>::is_iec559,
		"Branch::Real must satisfy IEC 559");

	Tree(LeafData root_data) : _cells{Cell(root_data)} { }

	void reserve(std::size_t cap) {
		_cells.reserve(cap);
	}
	void clear(LeafData root_data) {
		_cells.clear();
		_cells.push_back(Cell(root_data));
	}

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
	Branch const& branch(CellHandle cell) const {
		return _cells[cell].branch;
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
	Real volume_local(CellHandle parent, std::size_t child) const {
		return _cells[parent].branch.volume(child);
	}

	// Splits a leaf parent cell into new children.
	void split(
			CellHandle cell,
			Branch branch,
			BranchData branch_data,
			const LeafData (&leaf_data)[NUM_CHILDREN]) {
		// Append the new cells to the tree.
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
				cells_new.push_back(
					CellOther(fb(cell.branch_data), cell.branch));
			} else {
				cells_new.push_back(
					CellOther(fl(cell.leaf_data)));
			}
		}
		Tree<B, DBP, DLP> result(std::move(cells_new));
		return result;
	}
};

// Stores statistics about function evaluation in a region of space.
template<typename R>
struct Stats final {
	R mean_total;
	R mom2_total;
	R mom3_total;
	R mom4_total;
	std::size_t count;

	Stats() :
		mean_total(0),
		mom2_total(0),
		mom3_total(0),
		mom4_total(0),
		count(0) { }
	Stats(
		R mean_total,
		R mom2_total, R mom3_total, R mom4_total,
		std::size_t count) :
		mean_total(mean_total),
		mom2_total(mom2_total),
		mom3_total(mom3_total),
		mom4_total(mom4_total),
		count(count) { }

	void update(R measure) {
		R measure_sq = sq(measure);
		mean_total += measure;
		mom2_total += measure_sq;
		mom3_total += measure * measure_sq;
		mom4_total += sq(measure_sq);
		count += 1;
	}

	Stats<R>& operator+=(Stats<R> const& rhs) {
		mean_total += rhs.mean_total;
		mom2_total += rhs.mom2_total;
		mom3_total += rhs.mom3_total;
		mom4_total += rhs.mom4_total;
		count += rhs.count;
		return *this;
	}
	Stats<R>& operator-=(Stats<R> const& rhs) {
		mean_total -= rhs.mean_total;
		mom2_total -= rhs.mom2_total;
		mom3_total -= rhs.mom3_total;
		mom4_total -= rhs.mom4_total;
		count -= rhs.count;
		return *this;
	}
	Stats<R>& operator*=(R scale) {
		R scale_sq = sq(scale);
		mean_total *= scale;
		mom2_total *= scale_sq;
		mom3_total *= scale * scale_sq;
		mom4_total *= sq(scale_sq);
		return *this;
	}

	R est_mean(R* mean_err_out=nullptr) const {
		R mean = mean_total / count;
		if (mean_err_out != nullptr) {
			R var = clamp_above<R>((mom2_total - mean * mean_total) / (count - 1));
			*mean_err_out = std::sqrt(var / count);
		}
		return mean;
	}
	R est_var(R* var_err_out=nullptr) const {
		R mean = mean_total / count;
		R mean_sq = sq(mean);
		R var = clamp_above<R>((mom2_total - mean * mean_total) / (count - 1));
		if (var_err_out != nullptr) {
			// Central fourth moment (biased).
			R cmom4 = (
				mom4_total
				- 4 * mean * mom3_total
				+ 6 * mean_sq * mom2_total
				- 3 * mean_sq * mean * mean_total) / count;
			*var_err_out = std::sqrt(clamp_above<R>((cmom4 - sq(var)) / count));
		}
		return var;
	}
	R est_rel_var(R* rel_var_err_out=nullptr) const {
		R mean_err;
		R mean = est_mean(&mean_err);
		R var_err;
		R var = est_var(&var_err);
		// Small correction to relative variance for bias.
		R rel_var_0 = clamp<R>(var / sq(mean), 0, count);
		R correction = -3 * rel_var_0 / count;
		// Correction must be small.
		if (!(-0.5 < correction && correction < 0.5)) {
			correction = 0;
		}
		R rel_var = rel_var_0 * (1 + correction);
		if (rel_var_err_out != nullptr) {
			// TODO: Right now, the variance estimate is very bad! It assumes
			// that the mean and variance are independent and Gaussian, which is
			// most likely not true.
			R err = clamp<R>(
				1 / sq(mean) * std::sqrt(
					4 * var * rel_var_0 * sq(mean_err) + sq(var_err)),
				0, std::numeric_limits<R>::infinity());
			*rel_var_err_out = err;
		}
		return rel_var;
	}
	R est_prime(R* prime_err_out=nullptr) const {
		// Unbiased estimates of second moment with its variance.
		R mom2 = mom2_total / count;
		// Estimate the ideal prime value, which is given by:
		//     F_c = \sqrt{V_c \int dx f^2(x)}
		// with a small correction for bias added on.
		R prime_0 = std::sqrt(mom2);
		R mom4 = mom4_total / count;
		R ratio = clamp<R>(mom4 / sq(mom2), 0, count);
		R correction = clamp_above<R>((ratio - 1) / 8 / count);
		// The correction must be small.
		if (!(-0.5 < correction && correction < 0.5)) {
			correction = 0;
		}
		R prime = prime_0 * (1 + correction);
		// Error in estimate of prime value.
		if (prime_err_out != nullptr) {
			R ratio_mom2 = clamp<R>(mom4 / mom2 - mom2, 0, mom2_total);
			*prime_err_out = 0.5 * std::sqrt(ratio_mom2 / count);
		}
		return prime;
	}
	R est_explore_vol(R alpha) const {
		R prime = est_prime();
		R mean = est_mean();
		R explore_vol = std::pow(prime - mean, 1 / (alpha + 1));
		return explore_vol;
	}
	R est_tune_vol(R prime_tot) const {
		R prime = est_prime();
		R mom2 = mom2_total / count;
		R var2 = (mom4_total - mom2 * mom2_total) / (count - 1);
		R ratio = clamp<R>(var2 / sq(mom2), 0, count);
		R tune_vol = 0.5 * std::sqrt(
			(prime / prime_tot) * ((prime_tot - prime) / prime_tot) * ratio);
		return tune_vol;
	}
};
template<typename R>
inline Stats<R> operator+(Stats<R> lhs, Stats<R> const& rhs) {
	lhs += rhs;
	return lhs;
}
template<typename R>
inline Stats<R> operator-(Stats<R> lhs, Stats<R> const& rhs) {
	lhs -= rhs;
	return lhs;
}
template<typename R>
inline Stats<R> operator*(R lhs, Stats<R> rhs) {
	rhs *= lhs;
	return rhs;
}
template<typename R>
inline Stats<R> operator*(Stats<R> lhs, R rhs) {
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
	static_assert(
		std::is_floating_point<R>::value,
		"R must be a floating point type.");
	static_assert(
		std::numeric_limits<R>::is_iec559,
		"R must satisfy IEC 559");
	using FR = decltype(std::declval<F>()(std::declval<Point<D, R> >()));
	static_assert(
		std::is_same<R, FR>::value,
		"Function must return same floating point type as argument.");

	// The reason that we don't store a `Stats` and instead keep all of these
	// variables separate is because the `prime` and `volatility` measures are
	// not linear. So, the stats for different leafs can't just be added to get
	// the stats for a branch.
	struct BranchData final {
		// Estimate of mean within all children.
		R mean;
		// Total prime of all children.
		R prime;
		// Total exploration volatility of all children.
		R explore_vol;
		// Total tuning volatility of all children.
		R tune_vol;
		BranchData() : mean(0), prime(0), explore_vol(0), tune_vol(0) { }
	};
	struct LeafData final {
		Stats<R> stats;
		LeafData() : stats() { }
		LeafData(Stats<R> const& stats) : stats(stats) { }
	};
	using TreeType = Tree<KdBranch<D, R>, BranchData, LeafData>;

	// Function.
	F _func;
	// Spatial tree storing the statistics.
	TreeType _tree;
	// Random number engine.
	RndL _rnd_left;
	// The single-cell prime of the cell generator.
	R _prime_init;
	R _prime_init_err;
	R _mean_init;
	R _mean_init_err;
	R _scale_exp;

public:
	// Parameters that control various aspects of the exploration and tuning
	// process.

	// The target relative variance that must be reached in all cells.
	R target_rel_var = 0.01;
	// The precision (in standard deviations) to which the relative variance
	// must be reached.
	R target_rel_var_sigma = 1.5;
	// The minimum cell volume (with volume equal to length^d) at which
	// exploration terminates.
	R min_cell_length = 0.01;
	// An exponent related to how quickly the efficiency increases with more
	// cell divisions. It depends on the type of tree and the underlying
	// distribution. Empirically, lies between 0.5 and 2.0 for most functions.
	R scale_exp_init = 1.0;
	// Whether to adjust the scaling exponent in response to the behaviour of
	// the distribution at smaller and smaller scales.
	bool scale_exp_adjust = true;
	// Number of past points to consider when adjusting the scaling exponent.
	std::size_t scale_exp_points = 6;
	// The maximum cell volatility needed before exploration can ever terminate
	// in a cell. Generally, this should be between 0.1 and 1 (usually closely
	// to 1).
	R max_cell_explore_vol = 1;
	// The minimum cell volatility at which exploration terminates. A larger
	// value pushes for a more uniform efficiency over the entire generator
	// including in low-prime cells, where a smaller value will cause
	// exploration to focus on the high-prime cells.
	R min_cell_explore_vol = 1;
	// Chi-squared needed for a cell division to be triggered.
	R chi_squared_for_divide = 2;
	// Number of bins in edge histograms.
	std::size_t hist_num_bins = 128;
	// Minimum number of samples needed per bin before starting cell division.
	std::size_t hist_num_per_bin = 5;
	// Minimum and maximum number of samples to be taken in a single cell during
	// the exploration phase.
	std::size_t min_cell_explore_samples = 512;
	std::size_t max_cell_explore_samples = 524288;
	// Maximum number of cells to create during exploration.
	std::size_t max_explore_cells = 2*2*2*2*65536;
	// The accuracy to which to tune within.
	R tune_rel_accuracy = 0.01;
	// How many stages to tune. More stages gives slightly better tuning
	// performance.
	std::size_t tune_num_stages = 10;
	// Maximum number of samples to be taken total in the tuning phase.
	std::size_t max_tune_samples = 16777216;

private:
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

	// Fills in means for a cell and all of its descendents.
	R fill_means(CellHandle cell) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).mean = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				_tree.branch_data(parent).mean += fill_means(child);
			}
			return _tree.branch_data(parent).mean;
		} else {
			Stats<R> const& stats = _tree.leaf_data(parent).stats;
			if (stats.count == 0) {
				return 0;
			} else {
				return _tree.leaf_data(parent).stats.est_mean();
			}
		}
	}
	// Fills in primes for a cell and all of its descendents.
	R fill_primes(CellHandle cell) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).prime = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				_tree.branch_data(parent).prime += fill_primes(child);
			}
			return _tree.branch_data(parent).prime;
		} else {
			return _tree.leaf_data(parent).stats.est_prime();
		}
	}
	// Fills in exploration volatilities for a cell and all of its descendents.
	R fill_explore_vols(CellHandle cell) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).explore_vol = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				_tree.branch_data(parent).explore_vol += fill_explore_vols(child);
			}
			return _tree.branch_data(parent).explore_vol;
		} else {
			return _tree.leaf_data(parent).stats.est_explore_vol(_scale_exp);
		}
	}
	// Fills in tuning volatilities for a cell and all of its descendents.
	R fill_tune_vols(CellHandle cell, R prime_total) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).tune_vol = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				_tree.branch_data(parent).tune_vol += fill_tune_vols(child, prime_total);
			}
			return _tree.branch_data(parent).tune_vol;
		} else {
			return _tree.leaf_data(parent).stats.est_tune_vol(prime_total);
		}
	}

	// Determines the efficiency of the cell generator within a specific cell.
	// the efficiency is good enough, leaves the cell alone, otherwise, adds new
	// children to the cell to improve its efficiency.
	void explore_cell(
			// The cell to explore.
			CellHandle cell,
			// Dimensions of the cell.
			Point<D, R> offset, Point<D, R> extent,
			// Estimate of total integral of distribution.
			R mean_total,
			// Estimate of exploration volatility of tree.
			R vol_total,
			// Maximum number of layers of new cells to add on to the tree.
			std::size_t depth) {
		if (depth == 0) {
			return;
		}
		CellType type;
		Point<D, R> offset_local[TreeType::NUM_CHILDREN];
		Point<D, R> extent_local[TreeType::NUM_CHILDREN];
		CellHandle child[TreeType::NUM_CHILDREN];
		Stats<R> stats_init;
		#ifdef BUBBLE_USE_OPENMP
		#pragma omp critical (bubble_tree)
		#endif
		{
			type = _tree.type(cell);
			if (type == CellType::BRANCH) {
				for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
					offset_local[child_idx] = _tree.offset_local(cell, child_idx);
					extent_local[child_idx] = _tree.extent_local(cell, child_idx);
					child[child_idx] = _tree.child(cell, child_idx);
				}
			} else {
				stats_init = _tree.leaf_data(cell).stats;
			}
		}
		if (type == CellType::BRANCH) {
			// If the cell has children, then recursively explore them.
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				#ifdef BUBBLE_USE_OPENMP
				#pragma omp task
				#endif
				{
					Point<D, R> offset_child;
					Point<D, R> extent_child;
					for (Dim dim = 0; dim < D; ++dim) {
						offset_child[dim] = offset[dim];
						offset_child[dim] += extent[dim] * offset_local[child_idx][dim];
						extent_child[dim] = extent[dim] * extent_local[child_idx][dim];
					}
					explore_cell(
						child[child_idx],
						offset_child, extent_child,
						mean_total, vol_total,
						depth);
				}
			}
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp taskwait
			#endif
			return;
		} else {
			R volume = 1;
			for (Dim dim = 0; dim < D; ++dim) {
				volume *= extent[dim];
			}
			Stats<R> stats;
			// Histograms, one for each axis of the hypercube.
			std::vector<Stats<R> > hist_stats[D];
			for (Dim dim = 0; dim < D; ++dim) {
				hist_stats[dim].resize(hist_num_bins);
			}

			// Create a random number generator (a "right" generator) which will
			// be used.
			std::vector<Seed> seed;
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp critical (bubble_rand)
			#endif
			{
				seed = sample_rnd_left();
			}
			std::seed_seq seed_seq(seed.begin(), seed.end());
			RndR rnd(seed_seq);

			// The next number at which termination conditions will be checked.
			std::size_t next_samples = min_cell_explore_samples;
			while (true) {
				// Check termination conditions.
				Stats<R> stats_total = stats + stats_init;
				if (stats_total.count >= min_cell_explore_samples) {
					R prime_err;
					// We use `stats` instead of `stats_total` for the prime so
					// that it corresponds with the total statistics registered
					// in the histogram.
					R prime = stats.est_prime(&prime_err);
					R rel_var_err;
					R rel_var = stats_total.est_rel_var(&rel_var_err);

					// Termination condition 1. Efficiency meets criteria.
					// First check the cell division volatility. If the ratio is
					// greater than one, then check that the relative variance
					// has met the target accounting for uncertainty (and if so,
					// then terminate).
					R cell_vol =
						std::pow(
							vol_total / (0.5 * target_rel_var * mean_total),
							1 / _scale_exp)
						* stats_total.est_explore_vol(_scale_exp);
					// Termination will only be considered if the cell
					// volatility hits the threshold value. This ensures that no
					// matter what other termination conditions are considered,
					// the final cell generator will have good efficiency.
					if (cell_vol < max_cell_explore_vol) {
						R rel_var_max =
							rel_var + target_rel_var_sigma * rel_var_err;
						// Three possibilities for terminating:
						// * Relative variance is smaller than target.
						// * Cell is unimportant either because:
						//   * Cell is low volume.
						//   * Cell has low volatility.
						bool var_cond = (rel_var_max <= target_rel_var);
						bool volume_cond = (volume <= std::pow(min_cell_length, D));
						bool vol_cond = (cell_vol <= min_cell_explore_vol);
						if (var_cond || volume_cond || vol_cond) {
							#ifdef BUBBLE_USE_OPENMP
							#pragma omp critical (bubble_tree)
							#endif
							{
								_tree.leaf_data(cell).stats += stats;
							}
							return;
						}
					}

					// Termination condition 2. Cell division.
					// For cell division to occur, the edge histograms must be
					// different enough (by chi-squared test) from constant.
					// Then, the division site that minimizes the new sum of
					// primes is chosen.
					if (stats.count >= hist_num_per_bin * hist_num_bins) {
						R chi_squared = 0;
						R prime_min = std::numeric_limits<R>::max();
						Dim dim_min = 0;
						std::size_t bin_min = hist_num_bins / 2;
						Stats<R> stats_lower_min;
						Stats<R> stats_upper_min;
						for (Dim dim = 0; dim < D; ++dim) {
							// Integrated statistics below and above the
							// proposed division site.
							std::vector<Stats<R> > hist_stats_lower(hist_num_bins + 1);
							std::vector<Stats<R> > hist_stats_upper(hist_num_bins + 1);
							for (std::size_t bin = 1; bin < hist_num_bins + 1; ++bin) {
								std::size_t rbin = hist_num_bins - bin;
								hist_stats_lower[bin] = hist_stats_lower[bin - 1]
									+ hist_stats[dim][bin - 1];
								hist_stats_upper[rbin] = hist_stats_upper[rbin + 1]
									+ hist_stats[dim][rbin];
							}
							// Test division at each bin boundary. Note that the
							// last bin is left off: such a data point would
							// have `count_upper == 0`.
							for (std::size_t bin = 1; bin < hist_num_bins; ++bin) {
								Stats<R> const& stats_lower = hist_stats_lower[bin];
								Stats<R> const& stats_upper = hist_stats_upper[bin];
								// If a bin doesn't have the counts needed to
								// compute a standard deviation, just drop it.
								if (stats_lower.count <= 3 || stats_upper.count <= 3) {
									continue;
								}
								// Estimate the total prime resulting from
								// splitting the cell at this bin boundary. Need
								// to scale by the volume fractions, to account
								// for volume being included in the statistics
								// measures.
								R lambda_1 = R(bin)
									/ hist_num_bins;
								R lambda_2 = R(hist_num_bins - bin)
									/ hist_num_bins;
								R prime_1 = stats_lower.est_prime();
								R prime_2 = stats_upper.est_prime();
								R prime_diff = -prime
									+ lambda_1 * prime_1 + lambda_2 * prime_2;
								// The variance of the difference of primes is
								// given to leading order by this somewhat
								// complicated expression. Essentially, the
								// variance is reduced somewhat because of
								// positive correlations between the original
								// prime and the two new primes.
								R count = stats.count;
								R mom2 = stats.mom2_total / count;
								R mom4 = stats.mom4_total / count;
								R mom2_1 = lambda_1
									* stats_lower.mom2_total / count;
								R mom4_1 = sq(lambda_1) * lambda_1
									* stats_lower.mom4_total / count;
								R mom2_2 = lambda_2
									* stats_upper.mom2_total / count;
								R mom4_2 = sq(lambda_2) * lambda_2
									* stats_upper.mom4_total / count;
								R prime_var = 1. / (4. * count)
									* (mom4 - sq(mom2)) / mom2;
								R prime_var_1 = 1. / (4. * lambda_1 * count)
									* (mom4_1 - sq(mom2_1)) / mom2_1;
								R prime_var_2 = 1. / (4. * lambda_2 * count)
									* (mom4_2 - sq(mom2_2)) / mom2_2;
								R prime_cov_1 = 1. / (4. * count)
									* (mom4_1 - sq(mom2_1))
									/ std::sqrt(mom2_1 * mom2);
								R prime_cov_2 = 1. / (4. * count)
									* (mom4_2 - sq(mom2_2))
									/ std::sqrt(mom2_2 * mom2);
								R prime_diff_var = clamp_above<R>(
									prime_var + prime_var_1 + prime_var_2
									- 2. * prime_cov_1 - 2. * prime_cov_2);
								// Chi squared test against the maximum possible
								// value for sum of primes, which is the current
								// prime of the cell.
								// TODO: Using chi-squared test here even though
								// the prime difference is not normal.
								chi_squared += sq(prime_diff) / prime_diff_var;
								// Update minimum.
								if (prime_diff < prime_min) {
									prime_min = prime_diff;
									dim_min = dim;
									bin_min = bin;
									stats_lower_min = stats_lower;
									stats_upper_min = stats_upper;
								}
							}
						}
						// Average the chi-squared over every degree of freedom.
						chi_squared /= D * (hist_num_bins - 1);
						// If the chi-squared is large enough, or if just too
						// many samples have been requested, then make a
						// division.
						if (chi_squared > chi_squared_for_divide
								|| stats.count >= max_cell_explore_samples) {
							R lambda_1 = R(bin_min + 1)
								/ hist_num_bins;
							R lambda_2 = R(hist_num_bins - bin_min - 1)
								/ hist_num_bins;
							LeafData leaf_data[2] = {
								{ lambda_1 * stats_lower_min },
								{ lambda_2 * stats_upper_min },
							};
							#ifdef BUBBLE_USE_OPENMP
							#pragma omp critical (bubble_tree)
							#endif
							{
								_tree.split(
									cell,
									typename TreeType::Branch(dim_min, lambda_1),
									BranchData(),
									leaf_data);
							}
							// Recursively calling will cause the new children
							// to be explored. Reduce depth by 1 due to the
							// split.
							explore_cell(
								cell,
								offset, extent,
								mean_total, vol_total,
								depth - 1);
							return;
						}
					}
				}

				// If non of the termination conditions match, then produce more
				// samples.
				for (std::size_t sample = stats.count; sample < next_samples; ++sample) {
					// Generate a point randomly within the cell and evaluate
					// the function. First choose a histogram bin to generate
					// in, then produce the point itself.
					Point<D, R> point;
					std::size_t hist_idx[D];
					std::uniform_int_distribution<std::size_t> idx_dist(
						0,
						hist_num_bins - 1);
					for (Dim dim = 0; dim < D; ++dim) {
						std::size_t idx = idx_dist(rnd);
						hist_idx[dim] = idx;
						R lower = offset[dim]
							+ extent[dim] * idx / hist_num_bins;
						R upper = offset[dim]
							+ extent[dim] * (idx + 1) / hist_num_bins;
						std::uniform_real_distribution<R> x_dist(lower, upper);
						point[dim] = x_dist(rnd);
					}
					// First choose a histogram index to generate in, then make
					// the point itself.
					// TODO: Consider computing these with a more numerically
					// accurate method, such as an averaging tree.
					R f = volume * _func(point);
					stats.update(f);
					// Fill histograms.
					for (Dim dim = 0; dim < D; ++dim) {
						hist_stats[dim][hist_idx[dim]].update(f);
					}
				}
				// Double how many samples we produce each round.
				next_samples *= 2;
			}
		}
	}

	// Improves the estimate of the prime of a given cell.
	void tune_cell(
			// The cell to tune.
			CellHandle cell,
			// Dimensions of the cell.
			Point<D, R> offset, Point<D, R> extent,
			// Total prime of the entire cell generator.
			R prime_total,
			// Number of samples to use for estimating the prime.
			std::size_t samples) {
		if (_tree.type(cell) == CellType::BRANCH) {
			// For a branch, estimate the volatility of each child. The samples
			// will be divided up according to that.
			R vols[TreeType::NUM_CHILDREN];
			R vol_total = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(cell, child_idx);
				if (_tree.type(child) == CellType::BRANCH) {
					vols[child_idx] = _tree.branch_data(child).tune_vol;
				} else {
					vols[child_idx] = _tree.leaf_data(child).stats
						.est_tune_vol(prime_total);
				}
				vol_total += vols[child_idx];
			}
			std::size_t samples_child[TreeType::NUM_CHILDREN];
			std::size_t samples_total = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				R vol_frac = clamp<R>(vols[child_idx] / vol_total, 0, 1);
				samples_child[child_idx] = static_cast<std::size_t>(
					vol_frac * samples);
				samples_total += samples_child[child_idx];
			}
			while (samples_total < samples) {
				samples_child[samples_total % TreeType::NUM_CHILDREN] += 1;
				samples_total += 1;
			}
			while (samples_total > samples) {
				samples_child[samples_total % TreeType::NUM_CHILDREN] -= 1;
				samples_total -= 1;
			}

			// Sample children according to volatility.
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				#ifdef BUBBLE_USE_OPENMP
				#pragma omp task
				#endif
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
						offset_child[dim] = offset[dim];
						offset_child[dim] += extent[dim] * offset_local[dim];
						extent_child[dim] = extent[dim] * extent_local[dim];
					}
					tune_cell(
						child,
						offset_child, extent_child,
						prime_total,
						samples_child[child_idx]);
				}
			}
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp taskwait
			#endif
			return;
		} else {
			// For a leaf, there are no children, so all of the samples will be
			// uniformly distributed.
			R volume = 1;
			for (Dim dim = 0; dim < D; ++dim) {
				volume *= extent[dim];
			}
			Stats<R> stats;

			// Create a random number generator.
			std::vector<Seed> seed;
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp critical (bubble_rand)
			#endif
			{
				seed = sample_rnd_left();
			}
			std::seed_seq seed_seq(seed.begin(), seed.end());
			RndR rnd(seed_seq);

			for (std::size_t sample = 0; sample < samples; ++sample) {
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
		}
	}

public:
	CellBuilder(
			F func,
			std::size_t init_samples=16384,
			Seed seed=std::random_device()()) :
			_func(func),
			_tree(LeafData()),
			_rnd_left() {
		// Seed the left generator.
		std::seed_seq seed_seq { seed };
		_rnd_left.seed(seed_seq);
		// Sample a new generator for initialization.
		std::vector<Seed> seed_right = sample_rnd_left();
		std::seed_seq seed_seq_right(seed_right.begin(), seed_right.end());
		RndR rnd(seed_seq_right);
		// Do some sampling of the distribution before starting. Doing this
		// sampling gives an initial mean and prime estimate, and also lets us
		// catch some common errors.
		CellHandle root = _tree.root();
		Stats<R>& stats = _tree.leaf_data(root).stats;
		for (std::size_t sample = 0; sample < init_samples; ++sample) {
			// Assuming a volume of 1 for the unit hypercube.
			Point<D, R> point;
			for (Dim dim = 0; dim < D; ++dim) {
				std::uniform_real_distribution<R> x_dist(0, 1);
				point[dim] = x_dist(rnd);
			}
			R f = _func(point);
			if (!(f >= 0)) {
				throw std::runtime_error(
					"distribution is not non-negative everywhere");
			}
			stats.update(_func(point));
		}
		_prime_init = stats.est_prime(&_prime_init_err);
		_mean_init = stats.est_mean(&_mean_init_err);
		if (!std::isfinite(stats.mean_total)
				|| !std::isfinite(stats.mom2_total)
				|| !std::isfinite(stats.mom3_total)
				|| !std::isfinite(stats.mom4_total)) {
			throw std::runtime_error(
				"distribution statistics overflowed");
		}
		if (stats.mean_total <= 0
				|| stats.mom2_total <= 0
				|| stats.mom3_total <= 0
				|| stats.mom4_total <= 0) {
			throw std::runtime_error(
				"distribution cannot be distinguished from zero");
		}
	}

	// Explores the distribution to create a well-balanced cell generator that
	// can reach the desired performance. If unable to create a well-balanced
	// cell generator, returns estimated number of cells that the exploration
	// process fell short by.
	std::size_t explore() {
		Point<D, R> offset;
		Point<D, R> extent;
		offset.fill(0);
		extent.fill(1);
		CellHandle root = _tree.root();
		std::size_t cells;
		_scale_exp = scale_exp_init;
		std::vector<R> prev_primes;
		std::vector<std::size_t> prev_cells;
		do {
			cells = _tree.size();
			R mean_total = fill_means(root);
			R prime_total = fill_primes(root);

			prev_primes.push_back(prime_total);
			prev_cells.push_back(_tree.size());
			// Adjust the scaling exponent using a linear regression to some
			// number of previous steps.
			if (scale_exp_adjust && prev_primes.size() >= 2) {
				R mean_log_cells = 0;
				R mean_log_prime = 0;
				R mean_log_cells_sq = 0;
				R mean_log_prime_cells = 0;
				std::size_t idx_start = 0;
				if (prev_primes.size() >= scale_exp_points) {
					idx_start = prev_primes.size() - scale_exp_points;
				}
				std::size_t count = prev_primes.size() - idx_start;
				for (std::size_t idx = idx_start; idx < prev_primes.size(); ++idx) {
					auto log_cells = std::log(prev_cells[idx]);
					R log_prime = std::log(prev_primes[idx] - mean_total);
					mean_log_cells += log_cells / count;
					mean_log_prime += log_prime / count;
					mean_log_cells_sq += sq(log_cells) / count;
					mean_log_prime_cells += log_prime * log_cells / count;
				}
				R scale_exp_target =
					-(mean_log_prime_cells - mean_log_prime * mean_log_cells)
					/ (mean_log_cells_sq - sq(mean_log_cells));
				if (scale_exp_target < 0) {
					scale_exp_target = 0;
				}
				if (!std::isfinite(scale_exp_target)) {
					scale_exp_target = scale_exp_init;
				}
				_scale_exp = scale_exp_target;
			}
			R vol_total = fill_explore_vols(root);
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp parallel
			#pragma omp single
			#endif
			{
				explore_cell(
					root,
					offset, extent,
					mean_total, vol_total,
					1);
			}
		} while (cells != _tree.size() && _tree.size() < max_explore_cells);
		R mean_total = fill_means(root);
		R prime_total = fill_primes(root);
		R vol_total = fill_explore_vols(root);
		fill_tune_vols(root, prime_total);
		// Calculate ideal number of cells, and see if we surpassed that number.
		R ideal_cells = std::pow(
			std::pow(vol_total, _scale_exp + 1)
				/ (0.5 * target_rel_var * mean_total),
			1 / _scale_exp);
		if (ideal_cells > R(_tree.size())) {
			return ideal_cells - _tree.size();
		} else {
			return 0;
		}
	}

	// Tunes the cell generator to within some fraction of the ideal efficiency.
	// If unable to reach the desired efficiency, returns the estimated number
	// of samples that the tuning process fell short by.
	std::size_t tune() {
		Point<D, R> offset;
		Point<D, R> extent;
		offset.fill(0);
		extent.fill(1);
		CellHandle root = _tree.root();
		std::size_t undertuned_samples = 0;
		// The tuning process is made somewhat complicated by the staged
		// process. Since the accuracy of the volatilities might not be good if
		// the cell generator is not well tuned, the tuning is done in stages
		// with the volatilities recalculated between each stage.
		for (std::size_t stage = 0; stage < tune_num_stages; ++stage) {
			R prime_total = fill_primes(root);
			R vol_total = fill_tune_vols(root, prime_total);
			// Find the number of samples needed to reach the desired accuracy.
			R tune_accuracy = target_rel_var * tune_rel_accuracy;
			R tune_fraction = clamp<R>(R(stage + 1) / tune_num_stages, 0, 1);
			R samples_real = clamp<R>(
				tune_fraction * sq(vol_total) / tune_accuracy,
				// Multiply by 0.5 here so that the process of converting to `R`
				// won't cause `samples_real` to overflow the maximum possible
				// value.
				0, 0.5 * std::numeric_limits<std::size_t>::max());
			std::size_t samples = static_cast<std::size_t>(samples_real);
			std::size_t max_samples = static_cast<std::size_t>(
				tune_fraction * max_tune_samples);
			if (samples > max_samples) {
				undertuned_samples = samples - max_tune_samples;
				samples = max_tune_samples;
			}
			// Sample from the cell generator.
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp parallel
			#pragma omp single
			#endif
			{
				tune_cell(root, offset, extent, prime_total, samples);
			}
		}
		fill_means(root);
		R prime_total = fill_primes(root);
		fill_explore_vols(root);
		fill_tune_vols(root, prime_total);
		return undertuned_samples;
	}

	// Writes to file. A `CellBuilder` can be written to, but not read from, a
	// file. It reduces all of the statistics known about the distribution down
	// to a minimal set of data that can be constructed by a `CellGenerator`.
	void write(char const* file_name) const {
		std::ofstream file(file_name, std::ios::binary);
		file.exceptions(std::ios::failbit | std::ios::badbit);
		// Write identification string.
		file.write(FILE_HEAD, sizeof(FILE_HEAD));
		// Write sizes of important types, just as a quick validation.
		std::size_t size_size_t = sizeof(std::size_t);
		std::size_t size_tag = sizeof(unsigned char);
		std::size_t size_dim = sizeof(Dim);
		std::size_t size_split = sizeof(R);
		std::size_t size_prime = sizeof(R);
		file.write(reinterpret_cast<char const*>(&size_size_t), sizeof(std::size_t));
		file.write(reinterpret_cast<char const*>(&size_tag), sizeof(std::size_t));
		file.write(reinterpret_cast<char const*>(&size_dim), sizeof(std::size_t));
		file.write(reinterpret_cast<char const*>(&size_split), sizeof(std::size_t));
		file.write(reinterpret_cast<char const*>(&size_prime), sizeof(std::size_t));
		// Write total number of cells.
		std::size_t num_cells = _tree.size();
		file.write(reinterpret_cast<char const*>(&num_cells), sizeof(std::size_t));
		std::vector<CellHandle> cells;
		cells.push_back(_tree.root());
		while (!cells.empty()) {
			CellHandle cell = cells.back();
			cells.pop_back();
			if (_tree.type(cell) == CellType::BRANCH) {
				typename TreeType::Branch const& branch = _tree.branch(cell);
				// Branches store the following:
				// * Tag (unsigned char).
				// * Split dimension (unsigned char).
				// * Split location (real).
				unsigned char tag = 0;
				file.write(reinterpret_cast<char const*>(&tag), sizeof(unsigned char));
				file.write(reinterpret_cast<char const*>(&branch.dim), sizeof(Dim));
				file.write(reinterpret_cast<char const*>(&branch.split), sizeof(R));
				// Add branch children to be processed next, in order.
				for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
					cells.push_back(_tree.child(cell, child_idx));
				}
			} else {
				// Leafs store the following:
				// * Tag (unsigned char).
				// * Prime (integration real).
				unsigned char tag = 1;
				R prime = _tree.leaf_data(cell).stats.est_prime();
				file.write(reinterpret_cast<char const*>(&tag), sizeof(unsigned char));
				file.write(reinterpret_cast<char const*>(&prime), sizeof(R));
			}
		}
	}

	TreeType const& tree() const {
		return _tree;
	}

	R prime() const {
		CellHandle root = _tree.root();
		if (_tree.type(root) == CellType::BRANCH) {
			return _tree.branch_data(root).prime;
		} else {
			return _tree.leaf_data(root).stats.est_prime();
		}
	}
	R mean() const {
		CellHandle root = _tree.root();
		if (_tree.type(root) == CellType::BRANCH) {
			return _tree.branch_data(root).mean;
		} else {
			return _tree.leaf_data(root).stats.est_mean();
		}
	}
	R rel_var() const {
		return sq(prime() / mean()) - 1;
	}
	R scale_exp() const {
		return _scale_exp;
	}
};

template<Dim D, typename R, typename F>
class CellGenerator final {
	using Rnd = std::mt19937_64;

	// Return type of the function.
	static_assert(
		std::is_floating_point<R>::value,
		"R must be a floating point type.");
	static_assert(
		std::numeric_limits<R>::is_iec559,
		"R must satisfy IEC 559");
	using FR = decltype(std::declval<F>()(std::declval<Point<D, R> >()));
	static_assert(
		std::is_same<R, FR>::value,
		"Function must return same floating point type as argument.");

	// The cells store the prime value.
	struct Data final {
		R prime;
		Data() : prime(0) { }
		Data(R prime) : prime(prime) { }
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
	R fill_primes(CellHandle cell) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).prime = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				_tree.branch_data(parent).prime += fill_primes(child);
			}
			return _tree.branch_data(parent).prime;
		} else {
			return _tree.leaf_data(parent).prime;
		}
	}

	void select_cell(
			CellHandle cell,
			Point<D, R> offset, Point<D, R> extent,
			R choose_cell, Point<D, R> choose_point,
			R* weight_out, Point<D, R>* point_out) const {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				// Each child has a chance of being selected proportional to its
				// prime value.
				R prime_ratio = _tree.leaf_data(child).prime
					/ _tree.branch_data(parent).prime;
				if (choose_cell <= prime_ratio
						|| child_idx == TreeType::NUM_CHILDREN - 1) {
					// Scale `choose_cell` to be between 0 and 1 again.
					choose_cell /= prime_ratio;
					// Adjust offset and extent for the next level of recursion.
					Point<D, R> offset_local = _tree.offset_local(parent, child_idx);
					Point<D, R> extent_local = _tree.extent_local(parent, child_idx);
					for (Dim dim = 0; dim < D; ++dim) {
						offset[dim] += extent[dim] * offset_local[dim];
						extent[dim] *= extent_local[dim];
					}
					// Select the point from the child cell.
					select_cell(
						child,
						offset, extent,
						choose_cell, choose_point,
						weight_out, point_out);
					return;
				} else {
					choose_cell -= prime_ratio;
				}
			}
		} else if (_tree.type(parent) == CellType::LEAF) {
			// Rescale the function evaluation point.
			R volume = 1;
			for (Dim dim = 0; dim < D; ++dim) {
				(*point_out)[dim] = offset[dim];
				(*point_out)[dim] += extent[dim] * choose_point[dim];
				volume *= extent[dim];
			}
			R weight_norm = _tree.leaf_data(parent).prime;
			*weight_out = (volume * _func(*point_out)) / weight_norm;
		}
	}

public:
	// Creates a `CellGenerator` for use with a specific distribution.
	CellGenerator(
			F func,
			Seed seed=std::random_device()) :
			_func(func),
			_tree(Data{ 1 }),
			_rnd() {
		std::seed_seq seed_seq { seed };
		_rnd.seed(seed_seq);
		fill_primes(_tree.root());
	}
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
		fill_primes(_tree.root());
	}

	// Reads from a file.
	void read(char const* file_name) {
		std::ifstream file(file_name, std::ios::binary);
		file.exceptions(std::ios::failbit | std::ios::badbit);
		// Read identifier string.
		char head[sizeof(FILE_HEAD)];
		file.read(head, sizeof(FILE_HEAD));
		for(std::size_t idx = 0; idx < sizeof(FILE_HEAD); ++idx) {
			if (head[idx] != FILE_HEAD[idx]) {
				throw std::runtime_error("File header mismatch");
			}
		}
		// Read important data sizes.
		std::size_t size_size_t;
		std::size_t size_tag;
		std::size_t size_dim;
		std::size_t size_split;
		std::size_t size_prime;
		// `std::size_t` is especially important because the sizes themselves
		// are of that size.
		file.read(reinterpret_cast<char*>(&size_size_t), sizeof(std::size_t));
		if (size_size_t != sizeof(std::size_t)) {
			throw std::runtime_error("std::size_t size mismatch");
		}
		file.read(reinterpret_cast<char*>(&size_tag), sizeof(std::size_t));
		file.read(reinterpret_cast<char*>(&size_dim), sizeof(std::size_t));
		file.read(reinterpret_cast<char*>(&size_split), sizeof(std::size_t));
		file.read(reinterpret_cast<char*>(&size_prime), sizeof(std::size_t));
		if (size_tag != sizeof(unsigned char)) {
			throw std::runtime_error("unsigned char size mismatch");
		} else if (size_dim != sizeof(Dim)) {
			throw std::runtime_error("Dim size mismatch");
		} else if (size_split != sizeof(R)) {
			throw std::runtime_error("R size mismatch");
		} else if (size_prime != sizeof(R)) {
			throw std::runtime_error("R size mismatch");
		}
		// Read number of cells.
		std::size_t num_cells;
		file.read(reinterpret_cast<char*>(&num_cells), sizeof(std::size_t));
		// Read the cells themselves.
		_tree.clear(Data{ 0 });
		_tree.reserve(num_cells);
		std::vector<CellHandle> cells;
		cells.push_back(_tree.root());
		while (!cells.empty()) {
			// Read the next cell from the file.
			CellHandle cell = cells.back();
			cells.pop_back();
			unsigned char tag;
			file.read(reinterpret_cast<char*>(&tag), sizeof(unsigned char));
			if (tag == 0) {
				// Branch.
				Dim dim;
				R split;
				file.read(reinterpret_cast<char*>(&dim), sizeof(Dim));
				file.read(reinterpret_cast<char*>(&split), sizeof(R));
				_tree.split(
					cell,
					typename TreeType::Branch(dim, split),
					Data(0),
					{ Data(0), Data(0) });
				for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
					cells.push_back(_tree.child(cell, child_idx));
				}
			} else if (tag == 1) {
				// Leaf.
				R prime;
				file.read(reinterpret_cast<char*>(&prime), sizeof(R));
				_tree.leaf_data(cell).prime = prime;
			} else {
				throw std::runtime_error("invalid tag");
			}
		}
		if (_tree.size() > num_cells) {
			throw std::runtime_error("wrong number of cells");
		}
		/*
		if (file) {
			throw std::runtime_error("unexpected data");
		}
		*/
		fill_primes(_tree.root());
	}

	void generate(R* weight_out, Point<D, R>* point_out) {
		// Generate random number.
		std::uniform_real_distribution<R> dist(0, 1);
		R choose_cell = dist(_rnd);
		Point<D, R> choose_point;
		for (Dim dim = 0; dim < D; ++dim) {
			choose_point[dim] = dist(_rnd);
		}
		// Draw from the root cell.
		Point<D, R> offset;
		Point<D, R> extent;
		offset.fill(0);
		extent.fill(1);
		select_cell(
			_tree.root(),
			offset, extent,
			choose_cell, choose_point,
			weight_out, point_out);
	}
};

// Convenience functions to allow for template parameter inference.
template<Dim D, typename R, typename F>
inline CellBuilder<D, R, F> make_builder(
		F func,
		std::size_t init_samples=16384,
		Seed seed=std::random_device()()) {
	return CellBuilder<D, R, F>(func, init_samples, seed);
}
template<Dim D, typename R, typename F>
inline CellGenerator<D, R, F> make_generator(
		F func,
		Seed seed=std::random_device()()) {
	return CellGenerator<D, R, F>(func, seed);
}
template<Dim D, typename R, typename F>
inline CellGenerator<D, R, F> make_generator(
		CellBuilder<D, R, F> const& builder,
		Seed seed=std::random_device()()) {
	return CellGenerator<D, R, F>(builder, seed);
}

}

#endif

