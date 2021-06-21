#ifndef BUBBLE_HPP__
#define BUBBLE_HPP__

#include <array>
#include <cmath>
#include <fstream>
#include <ios>
#include <limits>
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
template<typename FI>
struct Stats final {
	FI mean_total;
	FI mom2_total;
	FI mom3_total;
	FI mom4_total;
	std::size_t count;

	Stats() :
		mean_total(0),
		mom2_total(0),
		mom3_total(0),
		mom4_total(0),
		count(0) { }
	Stats(
		FI mean_total,
		FI mom2_total, FI mom3_total, FI mom4_total,
		std::size_t count) :
		mean_total(mean_total),
		mom2_total(mom2_total),
		mom3_total(mom3_total),
		mom4_total(mom4_total),
		count(count) { }

	void update(FI measure) {
		FI measure_sq = sq(measure);
		mean_total += measure;
		mom2_total += measure_sq;
		mom3_total += measure_sq * sq(measure_sq);
		mom4_total += sq(measure_sq);
		count += 1;
	}

	Stats<FI>& operator+=(Stats<FI> const& rhs) {
		mean_total += rhs.mean_total;
		mom2_total += rhs.mom2_total;
		mom3_total += rhs.mom3_total;
		mom4_total += rhs.mom4_total;
		count += rhs.count;
		return *this;
	}
	Stats<FI>& operator-=(Stats<FI> const& rhs) {
		mean_total -= rhs.mean_total;
		mom2_total -= rhs.mom2_total;
		mom3_total -= rhs.mom3_total;
		mom4_total -= rhs.mom4_total;
		count -= rhs.count;
		return *this;
	}
	Stats<FI>& operator*=(FI scale) {
		FI scale_sq = sq(scale);
		mean_total *= scale;
		mom2_total *= scale_sq;
		mom3_total *= scale * sq(scale_sq);
		mom4_total *= sq(scale_sq);
		return *this;
	}

	FI est_mean(FI* mean_err_out=nullptr) const {
		// Unbiased estimate of mean.
		FI mean = mean_total / count;
		// Unbiased estimate of variance.
		FI var = (mom2_total - sq(mean) * count) / (count - 1);
		// Error in estimate of mean.
		if (mean_err_out != nullptr) {
			*mean_err_out = std::sqrt(var / count);
		}
		return mean;
	}
	FI est_var(FI* var_err_out=nullptr) const {
		// Unbiased estimate of mean.
		FI mean = mean_total / count;
		FI mean_sq = sq(mean);
		// Sample second and fourth central moments (biased).
		FI m2 = mom2_total / count - sq(mean);
		FI m4 = (
			mom4_total
			- 4 * mean * mom3_total
			+ 6 * mean_sq * mom2_total) / count
			- 3 * sq(mean_sq);
		// H-statistic for the fourth central moment.
		FI coeff_n_1 = FI(2) / (count - 2)
			+ 0.5 * FI(1) / (count - 1)
			-0.5 * FI(9) / (count - 3);
		FI coeff_n_2 = 1
			+ FI(1) / (count - 1)
			- FI(6) / (count - 2)
			+ FI(9) / (count - 3);
		FI h4 = 3 * coeff_n_1 * sq(m2) + coeff_n_2 * m4;
		// H-statistic for the second central moment (or, sample variance).
		FI h2 = m2 * count / (count - 1);
		if (var_err_out != nullptr) {
			// TODO: This error calculation is a little more complicated than it
			// needs to be. It's already quite biased, so the attempts to reduce
			// bias are likely ineffective.
			*var_err_out = std::sqrt(
				h4 / count - (FI(3) / count - FI(2) / (count - 1)) * sq(h2));
		}
		return h2;
	}
	FI est_rel_var(FI* rel_var_err_out=nullptr) const {
		// Mean estimate and variance.
		FI mean_err;
		FI mean = est_mean(&mean_err);
		// Variance estimate and variance.
		FI var_err;
		FI var = est_var(&var_err);
		// Relative variance, with a small correction for bias.
		FI rel_var = var / sq(mean) * (1 - 3 / count * var / sq(mean));
		if (rel_var_err_out != nullptr) {
			// TODO: Right now, the variance estimate is very bad! It assumes
			// that the mean and variance are independent and Gaussian, which is
			// most likely not true.
			*rel_var_err_out = rel_var
				* std::sqrt(4 * sq(mean_err / mean) + sq(var_err / var));
		}
		return rel_var;
	}
	FI est_prime(FI* prime_err_out=nullptr) const {
		// Unbiased estimate of the second moment.
		FI mom2 = mom2_total / count;
		// Unbiased estimate of variance of square.
		FI var2 = (mom4_total - sq(mom2) * count) / (count - 1);
		// Variance in estimate of second moment.
		FI mom2_var = var2 / count;
		// Estimate the ideal prime value, which is given by:
		//     F_c = \sqrt{V_c \int dx f^2(x)}
		// with a small correction for bias added on.
		FI prime = std::sqrt(mom2) * (1 + mom2_var / (8 * sq(mom2)));
		// Error in estimate of prime value.
		if (prime_err_out != nullptr) {
			*prime_err_out = 0.5 * std::sqrt(mom2_var / mom2);
		}
		return prime;
	}
	FI est_explore_vol(FI alpha) const {
		FI prime = est_prime();
		FI mean = est_mean();
		FI explore_vol = std::pow(prime - mean, 1 / (alpha + 1));
		return explore_vol;
	}
	FI est_tune_vol(FI prime_tot) const {
		FI prime = est_prime();
		FI mom2 = mom2_total / count;
		FI var2 = (mom4_total - sq(mom2) * count) / (count - 1);
		FI tune_vol = 0.5 * std::sqrt(
			(prime / prime_tot) * ((prime_tot - prime) / prime_tot)
			* (var2 / sq(mom2)));
		return tune_vol;
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

	// The reason that we don't store a `Stats` and instead keep all of these
	// variables separate is because the `prime` and `volatility` measures are
	// not linear. So, the stats for different leafs can't just be added to get
	// the stats for a branch.
	struct BranchData final {
		// Estimate of mean within all children.
		FI mean;
		// Total prime of all children.
		FI prime;
		// Total exploration volatility of all children.
		FI explore_vol;
		// Total tuning volatility of all children.
		FI tune_vol;
		BranchData() : mean(0), prime(0), explore_vol(0), tune_vol(0) { }
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

public:
	// Parameters that control various aspects of the exploration and tuning
	// process.

	// The target relative variance that must be reached in all cells.
	FI target_rel_var = 0.01;
	// The precision (in standard deviations) to which the relative variance
	// must be reached.
	FI target_rel_var_sigma = 1.5;
	// The minimum cell volume (with volume equal to length^d) at which
	// exploration terminates.
	FI min_cell_length = 0.01;
	// The exponent used when computing exploration volatility. This has been
	// empirically determined to be between 0.5 and 2.0 for k-d trees on smooth
	// functions.
	FI scale_exp = 1.0;
	// The maximum cell volatility needed before exploration can ever terminate
	// in a cell. Generally, this should be between 0.1 and 1 (usually closely
	// to 1).
	FI max_cell_explore_vol = 1;
	// The minimum cell volatility at which exploration terminates. A larger
	// value pushes for a more uniform efficiency over the entire generator
	// including in low-prime cells, where a smaller value will cause
	// exploration to focus on the high-prime cells.
	FI min_cell_explore_vol = 1;
	// Chi-squared needed for a cell division to be triggered.
	FI chi_squared_for_divide = 4;
	// Number of bins in edge histograms.
	std::size_t hist_num_bins = 100;
	// Minimum number of samples needed per bin before starting cell division.
	std::size_t hist_num_per_bin = 3;
	// Minimum and maximum number of samples to be taken in a single cell during
	// the exploration phase.
	std::size_t min_cell_explore_samples = 256;
	std::size_t max_cell_explore_samples = 524288;
	// Maximum number of cells to create during exploration.
	std::size_t max_explore_cells = 32768;
	// The accuracy to which to tune within.
	FI tune_rel_accuracy = 0.01;
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
	FI fill_means(CellHandle cell) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).mean = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				_tree.branch_data(parent).mean += fill_means(child);
			}
			return _tree.branch_data(parent).mean;
		} else {
			Stats<FI> const& stats = _tree.leaf_data(parent).stats;
			if (stats.count == 0) {
				return 0;
			} else {
				return _tree.leaf_data(parent).stats.est_mean();
			}
		}
	}
	// Fills in primes for a cell and all of its descendents.
	FI fill_primes(CellHandle cell) {
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
	FI fill_explore_vols(CellHandle cell) {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			_tree.branch_data(parent).explore_vol = 0;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				_tree.branch_data(parent).explore_vol += fill_explore_vols(child);
			}
			return _tree.branch_data(parent).explore_vol;
		} else {
			return _tree.leaf_data(parent).stats.est_explore_vol(scale_exp);
		}
	}
	// Fills in tuning volatilities for a cell and all of its descendents.
	FI fill_tune_vols(CellHandle cell, FI prime_total) {
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
			FI mean_total,
			// Estimate of exploration volatility of tree.
			FI vol_total,
			// Maximum number of layers of new cells to add on to the tree.
			std::size_t depth) {
		if (depth == 0) {
			return;
		}
		CellType type;
		#pragma omp critical (bubble_tree)
		{
			type = _tree.type(cell);
		}
		if (type == CellType::BRANCH) {
			// If the cell has children, then recursively explore them.
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				#pragma omp task
				{
					Point<D, R> offset_local;
					Point<D, R> extent_local;
					CellHandle child;
					#pragma omp critical (bubble_tree)
					{
						offset_local = _tree.offset_local(cell, child_idx);
						extent_local = _tree.extent_local(cell, child_idx);
						child = _tree.child(cell, child_idx);
					}
					Point<D, R> offset_child;
					Point<D, R> extent_child;
					for (Dim dim = 0; dim < D; ++dim) {
						offset_child[dim] = offset[dim];
						offset_child[dim] += extent[dim] * offset_local[dim];
						extent_child[dim] = extent[dim] * extent_local[dim];
					}
					explore_cell(
						child,
						offset_child, extent_child,
						mean_total, vol_total,
						depth);
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
			#pragma omp critical (bubble_tree)
			{
				stats_init = _tree.leaf_data(cell).stats;
			}
			// Histograms, one for each axis of the hypercube.
			std::vector<Stats<FI> > hist_stats[D];
			for (Dim dim = 0; dim < D; ++dim) {
				hist_stats[dim].resize(hist_num_bins);
			}

			// Create a random number generator (a "right" generator) which will
			// be used.
			std::vector<Seed> seed;
			#pragma omp critical (bubble_rand)
			{
				seed = sample_rnd_left();
			}
			std::seed_seq seed_seq(seed.begin(), seed.end());
			RndR rnd(seed_seq);

			// The next number at which termination conditions will be checked.
			std::size_t next_samples = min_cell_explore_samples;
			while (true) {
				// Check termination conditions.
				Stats<FI> stats_total = stats + stats_init;
				if (stats_total.count >= min_cell_explore_samples) {
					FI prime_err;
					FI prime = stats_total.est_prime(&prime_err);
					FI rel_var_err;
					FI rel_var = stats_total.est_rel_var(&rel_var_err);

					// Termination condition 1. Efficiency meets criteria.
					// First check the cell division volatility. If the ratio is
					// greater than one, then check that the relative variance
					// has met the target accounting for uncertainty (and if so,
					// then terminate).
					FI cell_vol =
						std::pow(
							vol_total / (0.5 * target_rel_var * mean_total),
							1 / scale_exp)
						* stats_total.est_explore_vol(scale_exp);
					// Termination will only be considered if the cell
					// volatility hits the threshold value. This ensures that no
					// matter what other termination conditions are considered,
					// the final cell generator will have good efficiency.
					if (cell_vol < max_cell_explore_vol) {
						FI rel_var_max =
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
							#pragma omp critical (bubble_tree)
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
						FI chi_squared = 0;
						FI prime_min = std::numeric_limits<FI>::max();
						Dim dim_min = 0;
						std::size_t bin_min = hist_num_bins / 2;
						Stats<FI> stats_lower_min;
						Stats<FI> stats_upper_min;
						for (Dim dim = 0; dim < D; ++dim) {
							// Integrated statistics below and above the
							// proposed division site.
							std::vector<Stats<FI> > hist_stats_lower(hist_num_bins + 1);
							std::vector<Stats<FI> > hist_stats_upper(hist_num_bins + 1);
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
								Stats<FI> const& stats_lower = hist_stats_lower[bin];
								Stats<FI> const& stats_upper = hist_stats_upper[bin];
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
								R lambda_1 = R(bin + 1)
									/ hist_num_bins;
								R lambda_2 = R(hist_num_bins - bin - 1)
									/ hist_num_bins;
								FI prime_var_lower;
								FI prime_lower = (lambda_1 * stats_lower)
									.est_prime(&prime_var_lower);
								FI prime_var_upper;
								FI prime_upper = (lambda_2 * stats_upper)
									.est_prime(&prime_var_upper);
								FI prime_new = prime_lower + prime_upper;
								// Because the two prime estimates are
								// independent, the variances can be summed.
								FI prime_var_new = prime_var_lower + prime_var_upper;
								// Chi squared test against the maximum possible
								// value for sum of primes, which is the current
								// prime of the cell.
								// TODO: The `prime_var_new` and `prime_err` are
								// not independent, and so technically shouldn't
								// be added. Not sure if there's enough
								// dependence for it to matter though.
								// TODO: Using chi-squared test here even though
								// prime is not normally distributed. At best it
								// is a square root of a normal variable (from
								// central limit theorem, which doesn't
								// necessary work that well for all
								// distributions). It's not clear what a better
								// test would be though.
								chi_squared += sq(prime_new - prime)
									/ (sq(prime_var_new) + sq(prime_err));
								// Update minimum.
								if (prime_new < prime_min) {
									prime_min = prime_new;
									dim_min = dim;
									bin_min = bin;
									stats_lower_min = stats_lower;
									stats_upper_min = stats_upper;
								}
							}
						}
						// Average the chi-squared over every degree of freedom.
						chi_squared /= D * hist_num_bins;
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
							#pragma omp critical (bubble_tree)
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
					FI f = volume * _func(point);
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
			FI prime_total,
			// Number of samples to use for estimating the prime.
			std::size_t samples) {
		if (_tree.type(cell) == CellType::BRANCH) {
			// For a branch, estimate the volatility of each child. The samples
			// will be divided up according to that.
			FI vols[TreeType::NUM_CHILDREN];
			FI vol_total = 0;
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
				FI vol_frac = vols[child_idx] / vol_total;
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
			#pragma omp taskwait
			return;
		} else {
			// For a leaf, there are no children, so all of the samples will be
			// uniformly distributed.
			R volume = 1;
			for (Dim dim = 0; dim < D; ++dim) {
				volume *= extent[dim];
			}
			Stats<FI> stats;

			// Create a random number generator.
			std::vector<Seed> seed;
			#pragma omp critical (bubble_rand)
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
	CellBuilder(F func, Seed seed=std::random_device()()) :
			_func(func),
			_tree(LeafData()),
			_rnd_left() {
		std::seed_seq seed_seq { seed };
		_rnd_left.seed(seed_seq);
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
		// TODO: Adjust `scale_exp` to fit the distribution as exploration
		// proceeds.
		std::size_t cells;
		do {
			cells = _tree.size();
			FI mean_total = fill_means(root);
			FI vol_total = fill_explore_vols(root);
			#pragma omp parallel
			#pragma omp single
			{
				explore_cell(
					root,
					offset, extent,
					mean_total, vol_total,
					1);
			}
			FI prime_total = fill_primes(root);
		} while (cells != _tree.size() && _tree.size() < max_explore_cells);
		FI mean_total = fill_means(root);
		FI prime_total = fill_primes(root);
		FI vol_total = fill_explore_vols(root);
		fill_tune_vols(root, prime_total);
		// Calculate ideal number of cells, and see if we surpassed that number.
		FI ideal_cells = static_cast<std::size_t>(
			std::pow(vol_total, scale_exp + 1)
			/ (0.5 * target_rel_var * mean_total));
		if (ideal_cells > _tree.size()) {
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
			FI prime_total = fill_primes(root);
			FI vol_total = fill_tune_vols(root, prime_total);
			// Find the number of samples needed to reach the desired accuracy.
			FI tune_accuracy = target_rel_var * tune_rel_accuracy;
			std::size_t samples = static_cast<std::size_t>(
				sq(vol_total) / tune_accuracy);
			samples = static_cast<std::size_t>(
				R(stage + 1) / tune_num_stages * samples);
			std::size_t max_samples = static_cast<std::size_t>(
				R(stage + 1) / tune_num_stages * max_tune_samples);
			if (samples > max_samples) {
				undertuned_samples = samples - max_tune_samples;
				samples = max_tune_samples;
			}
			// Sample from the cell generator.
			#pragma omp parallel
			#pragma omp single
			{
				tune_cell(root, offset, extent, prime_total, samples);
			}
		}
		fill_means(root);
		FI prime_total = fill_primes(root);
		fill_explore_vols(root);
		fill_tune_vols(root, prime_total);
		return undertuned_samples;
	}

	FI prime() const {
		CellHandle root = _tree.root();
		if (_tree.type(root) == CellType::BRANCH) {
			return _tree.branch_data(root).prime;
		} else {
			return _tree.leaf_data(root).stats.est_prime();
		}
	}
	FI mean() const {
		CellHandle root = _tree.root();
		if (_tree.type(root) == CellType::BRANCH) {
			return _tree.branch_data(root).mean;
		} else {
			return _tree.leaf_data(root).stats.est_mean();
		}
	}
	FI tune_vol() const {
		CellHandle root = _tree.root();
		if (_tree.type(root) == CellType::BRANCH) {
			return _tree.branch_data(root).tune_vol;
		} else {
			return _tree.leaf_data(root).stats.est_tune_vol(prime());
		}
	}
	FI explore_vol() const {
		CellHandle root = _tree.root();
		if (_tree.type(root) == CellType::BRANCH) {
			return _tree.branch_data(root).explore_vol;
		} else {
			return _tree.leaf_data(root).stats.est_explore_vol(scale_exp);
		}
	}

	FI rel_var() const {
		return sq(prime() / mean()) - 1;
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
			FI* weight_out, Point<D, R>* point_out) const {
		CellHandle parent = cell;
		if (_tree.type(parent) == CellType::BRANCH) {
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(parent, child_idx);
				// Each child has a chance of being selected proportional to its
				// prime value.
				FI prime_ratio = _tree.leaf_data(child).prime
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
			FI weight_norm = _tree.leaf_data(parent).prime;
			*weight_out = (volume * _func(*point_out)) / weight_norm;
		}
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
		fill_primes(_tree.root());
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

	void generate(FI* weight_out, Point<D, R>* point_out) {
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

