#ifndef BUBBLE_HPP__
#define BUBBLE_HPP__

#include <array>
#include <cmath>
#include <functional>
#include <ios>
#include <istream>
#include <limits>
#include <ostream>
#include <queue>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace bubble {

// Identifier used in streams.
constexpr char STREAM_HEAD[] = "bubble v. 0.1 stream\n";

// Some utility functions.
template<typename T>
static T sq(T x) {
	return x * x;
}
template<typename T>
static T sqrt1p_1m(T x) {
	// Numerically accurate computation of `sqrt(1 + x) - 1`.
	return std::expm1(0.5 * std::log1p(x));
}
template<typename T>
static T pow_inv(T x, T p) {
	if (p == 0.) {
		if (x > 1.) {
			return std::numeric_limits<T>::infinity();
		} else if (x == 1.) {
			return 1.;
		} else if (x >= 0.) {
			return 0.;
		} else {
			return std::numeric_limits<T>::quiet_NaN();
		}
	} else {
		return std::pow(x, 1. / p);
	}
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

template<Dim D, typename R, typename F, typename=void>
struct IsDistribution {
	static constexpr bool value = false;
};
template<Dim D, typename R, typename F>
struct IsDistribution<D, R, F,
	typename std::enable_if<
		// Requires that `F` can take a `Point<D, R>` to a type convertible to
		// `R` without exceptions.
		noexcept(R(std::declval<F>()(std::declval<Point<D, R> >())))>::type> {
	static constexpr bool value = true;
};

// A branch node for a k-d type tree.
template<Dim D, typename R>
struct KdBranch final {
	using Real = R;
	static constexpr Dim DIM = D;

	static constexpr std::size_t NUM_CHILDREN = 2;
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

	static_assert(
		std::is_trivially_copyable<B>::value,
		"Branch node type must be trivially copyable.");
	static_assert(
		std::is_trivially_copyable<DB>::value,
		"Branch data type must be trivially copyable.");
	static_assert(
		std::is_trivially_copyable<DL>::value,
		"Leaf data type must be trivially copyable.");

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
	};

	std::vector<Cell> _cells;

	Tree(std::vector<Cell>&& cells) : _cells(std::move(cells)) { }

public:
	using Branch = B;
	using BranchData = DB;
	using LeafData = DL;
	static constexpr Dim DIM = Branch::DIM;
	using Real = typename Branch::Real;
	static constexpr std::size_t NUM_CHILDREN = Branch::NUM_CHILDREN;

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
	std::size_t leaf_size() const {
		return (_cells.size() - 1) / NUM_CHILDREN * (NUM_CHILDREN - 1) + 1;
	}
	std::size_t branch_size() const {
		return (_cells.size() - 1) / NUM_CHILDREN;
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

// Stores statistics about independent samples from an arbitrary positive
// distribution.
template<typename R>
class Stats final {
	static_assert(
		std::is_floating_point<R>::value,
		"R must be a floating point type.");
	static_assert(
		std::numeric_limits<R>::is_iec559,
		"R must satisfy IEC 559.");
	// TODO: Assert that R can hold the maximum value of size_t.
	
	// The number of moments stored in a statistics.
	static constexpr std::size_t NUM_MOMS = 4;

	// Uncentered moments of the distribution normalized to be dimensionless by
	// the maximum of the distribtion. This normalization is used to prevent
	// floating point overflow/underflow, especially of higher moments.
	std::array<R, NUM_MOMS> _mom;
	// The maximum is the only dimensionful measure.
	R _max;
	// Number of events observed.
	std::size_t _count;

public:
	Stats() :
		_mom(),
		_max(0.),
		_count(0.) { }
	// Implicitly converts a single sample into statistics.
	Stats(R x) :
			_mom(),
			_max(x),
			_count(1) {
		_mom.fill(1.);
	}
	// Builds a statistics from aggregate information.
	Stats(std::array<R, NUM_MOMS> mom_to_max, R max, std::size_t count) :
		_mom(mom_to_max),
		_max(max),
		_count(count) { }

	Stats<R>& operator+=(Stats<R> const& rhs) {
		// Get relative contributions from each statistics.
		R count_tot = _count + rhs._count;
		R lambda_1 = 0.;
		R lambda_2 = 0.;
		if (!(count_tot == 0.)) {
			lambda_1 = _count / count_tot;
			lambda_2 = rhs._count / count_tot;
		}
		// Combine the maximums.
		R max_tot = std::max(_max, rhs._max);
		R max_scale_1 = 0.;
		R max_scale_2 = 0.;
		if (!(max_tot == 0.)) {
			max_scale_1 = _max / max_tot;
			max_scale_2 = rhs._max / max_tot;
		}
		// Adjust the moments based on the new maximum.
		_max = max_tot;
		for (std::size_t idx = 0; idx < NUM_MOMS; ++idx) {
			lambda_1 *= max_scale_1;
			lambda_2 *= max_scale_2;
			_mom[idx] *= lambda_1;
			_mom[idx] += lambda_2 * rhs._mom[idx];
		}
		_count += rhs._count;
		return *this;
	}
	Stats<R>& operator*=(R scale) {
		_max *= scale;
		return *this;
	}

	std::size_t count() const {
		return _count;
	}
	R mean() const {
		return _max * _mom[0];
	}
	R max() const {
		return _max;
	}

	// Rescales the number of events that have been recorded. If `new_count` is
	// greater than the current count, this is equivalent to adding a number of
	// zero events. If `new_count` is smaller, then it is equivalent to
	// neglecting a number of zero events.
	void rescale_count(std::size_t new_count) {
		R lambda = static_cast<R>(_count) / new_count;
		for (std::size_t idx = 0; idx < NUM_MOMS; ++idx) {
			_mom[idx] *= lambda;
		}
		_count = new_count;
	}

	// Returns whether the measures are valid. Non-valid measures could indicate
	// that negative numbers were entered into the statistics, or that overflow
	// occured.
	bool valid() const {
		if (!(std::isfinite(_max) && _max >= 0.)) {
			return false;
		}
		for (std::size_t idx = 0; idx < NUM_MOMS; ++idx) {
			if (!(std::isfinite(_mom[idx]) && _mom[idx] >= 0. && _mom[idx] <= 1.)) {
				return false;
			}
			if (idx > 0 && _mom[idx] > _mom[idx - 1]) {
				return false;
			}
		}
		return true;
	}
	// Returns whether the statistics have underflowed, so that certain measures
	// are zero when they should be small but non-zero.
	bool underflow() const {
		// If the first moment is zero, then all of the entries have been zero.
		if (_mom[0] == 0.) {
			return false;
		}
		for (std::size_t idx = 1; idx < NUM_MOMS; ++idx) {
			// If the first moment is non-zero, then every other entry must be
			// non-zero as well.
			if (_mom[idx] == 0.) {
				return true;
			}
		}
		return false;
	}

	// Raw statistical moments. Generally, the ratios should be used instead of
	// these, since there is a risk of overflow/underflow, especially with
	// higher moments.
	R m1() const {
		return _max * _mom[0];
	}
	R m2() const {
		return _max * _max * _mom[1];
	}
	R m3() const {
		return _max * _max * _max * _mom[2];
	}
	R m4() const {
		return _max * _max * _max * _max * _mom[3];
	}
	R var1() const {
		return ratio_var1_to_max2() * _max * _max;
	}
	R var2() const {
		return ratio_var2_to_max4() * _max * _max * _max * _max;
	}
	R skew1() const {
		return ratio_skew1_to_max3() * _max * _max * _max;
	}
	R kurt1() const {
		return ratio_kurt1_to_max4() * _max * _max * _max * _max;
	}
	R max1() const {
		return _max;
	}
	R max2() const {
		return _max * _max;
	}
	R max3() const {
		return _max * _max * _max;
	}
	R max4() const {
		return _max * _max * _max * _max;
	}

	// Ratios and functions of statistical moments. These can be computed with
	// low risk of overflow/underflow. Overflow/underflow will not be caused by
	// an intermediate step of the calculation, but can still happen if the
	// result just can't be stored in a floating point number. Note that while
	// many of these ratios are involve unbiased quantities (ex. `var1` refers
	// to the sample variance), they are not unbiased themselves.
	R ratio_m1_to_max1() const {
		return _mom[0];
	}
	R ratio_m2_to_max2() const {
		return _mom[1];
	}
	R ratio_m3_to_max3() const {
		return _mom[2];
	}
	R ratio_m4_to_max4() const {
		return _mom[3];
	}
	R ratio_var1_to_max2() const {
		R correction = 1. / (1. - 1. / _count);
		R sample = _mom[1] - _mom[0] * _mom[0];
		return correction * sample;
	}
	R ratio_var2_to_max4() const {
		R correction = 1. / (1. - 1. / _count);
		R sample = _mom[3] - _mom[1] * _mom[1];
		return correction * sample;
	}
	R ratio_var1_to_m1p2() const {
		R m1p2 = _mom[0] * _mom[0];
		if (!(m1p2 == 0.)) {
			return ratio_var1_to_max2() / m1p2;
		} else {
			return _count;
		}
	}
	R ratio_var2_to_m2p2() const {
		R m2p2 = _mom[1] * _mom[1];
		if (!(m2p2 == 0.)) {
			return ratio_var2_to_max4() / m2p2;
		} else {
			return _count;
		}
	}
	R ratio_var1_to_max1_m1() const {
		if (!(_mom[0] == 0.)) {
			return ratio_var1_to_max2() / _mom[0];
		} else {
			// TODO: Verify that returning 1 is correct in this case.
			return 1.;
		}
	}
	R ratio_var2_to_max2_m2() const {
		if (!(_mom[1] == 0.)) {
			return ratio_var2_to_max4() / _mom[1];
		} else {
			return 1.;
		}
	}
	R ratio_var2_to_var1_m2() const {
		R denom = _mom[1] * ratio_var1_to_max2();
		R var2 = ratio_var2_to_max4();
		if (!(denom == 0.)) {
			return var2 / denom;
		} else {
			return R(_count) * _count;
		}
	}
	R ratio_skew1_to_max3() const {
		// Unbiased estimate of third central moment (often called skewness).
		R correction = 1. / ((1. - 1. / _count) * (1. - 2. / _count));
		R sample = _mom[2]
			- 3. * _mom[1] * _mom[0]
			+ 2. * _mom[0] * _mom[0] * _mom[0];
		return correction * sample;
	}
	R ratio_skew1_to_var1_m1() const {
		R denom = _mom[0] * ratio_var1_to_max2();
		R skew1 = ratio_skew1_to_max3();
		if (!(denom == 0.)) {
			return skew1 / denom;
		} else {
			// TODO: Determine what to return in this case.
			return 0.;
		}
	}
	R ratio_kurt1_to_max4() const {
		// Unbiased estimate of fourth central moment (often called kurtosis).
		R correction_denom = (1. - 3. / _count + 3. / _count / _count)
			* (1. - 1. / _count) * (1. - 1. / _count) * (1. - 1. / _count);
		R correction_1 = (
				(1. - 1. / _count) * (1. - 1. / _count)
				+ (6. - 9. / _count) / _count / _count)
			/ correction_denom;
		R correction_2 = -(6. - 9. / _count) / _count / correction_denom;
		R sample_1 = _mom[3]
			- 4. * _mom[2] * _mom[0]
			+ 6. * _mom[1] * _mom[0] * _mom[0]
			- 3. * _mom[0] * _mom[0] * _mom[0] * _mom[0];
		R sample_2 = (_mom[1] - _mom[0] * _mom[0])
			* (_mom[1] - _mom[0] * _mom[0]);
		return correction_1 * sample_1 + correction_2 * sample_2;
	}
	R ratio_kurt1_to_var1p2() const {
		R var1 = ratio_var1_to_max2();
		R denom = var1 * var1;
		R kurt1 = ratio_kurt1_to_max4();
		if (!(denom == 0.)) {
			return kurt1 / denom;
		} else {
			// TODO: Determine what to return in this case.
			return 0.;
		}
	}

	// Estimators with errors. These estimators include an error estimate for
	// certain quantities useful for generators.
	R est_mean(R* err_out=nullptr) const {
		R est = mean();
		if (err_out != nullptr) {
			*err_out = max() * std::sqrt(clamp_above<R>(
				ratio_var1_to_max2() / count()));
		}
		return est;
	}
	// TODO: No standard deviation estimator yet.
	R est_rel_var(R* err_out=nullptr) const {
		R est_0 = ratio_var1_to_m1p2();
		// Small correction to relative variance for bias.
		R correction = (2. * ratio_skew1_to_var1_m1() - 3. * est_0) / count();
		// Correction must be small.
		if (!(-0.5 < correction && correction < 0.5)) {
			correction = 0.;
		}
		R est = est_0 * (1. + correction);
		if (err_out != nullptr) {
			R kurt_term = ratio_kurt1_to_var1p2();
			R skew_term = ratio_skew1_to_var1_m1();
			R var_term = ratio_var1_to_m1p2();
			R term = kurt_term - 4. * skew_term + 4. * var_term - 1.;
			*err_out = est_0 * std::sqrt(clamp_above<R>(term / count()));
		}
		return est;
	}
	R est_sqrt_m2(R* err_out=nullptr) const {
		R est_0 = max() * std::sqrt(clamp_above<R>(ratio_m2_to_max2()));
		R correction = (1. / 8) * ratio_var2_to_m2p2() / count();
		if (!(-0.5 < correction && correction < 0.5)) {
			correction = 0.;
		}
		R est = est_0 * (1. + correction);
		if (err_out != nullptr) {
			*err_out = max() * std::sqrt(clamp_above<R>(
				0.25 * ratio_var2_to_max2_m2() / count()));
		}
		return est;
	}
};
template<typename R>
inline Stats<R> operator+(Stats<R> lhs, Stats<R> const& rhs) {
	lhs += rhs;
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

// Precisely accumulates statistics by using a tree structure to add together
// large numbers of statistics in the most numerically precise way.
template<typename R>
class StatsAccum final {
	// The stack holds unsummed statistics, in decreasing order of count.
	// Generally, we try to sum statistics only when they have the same number
	// of samples, because this gives the best precision on average.
	std::vector<Stats<R> > _stack;
	// The total so far.
	Stats<R> _total;

public:
	StatsAccum() {
		// To avoid unnecessary future allocations, make a single allocation on
		// construction that stores all of the space that will ever be needed,
		// presuming that the number of statistics doesn't overflow a `size_t`.
		std::size_t capacity = std::numeric_limits<std::size_t>::digits;
		_stack.reserve(capacity);
	}

	// Adds a new entry to the accumulator.
	StatsAccum<R>& operator+=(R x) {
		Stats<R> next = x;
		while (!_stack.empty() && next.count() == _stack.back().count()) {
			next += _stack.back();
			_stack.pop_back();
		}
		// Reset the total to the stack everytime the stack empties out (which
		// happens on powers of 2). Otherwise, increment the total normally.
		if (_stack.empty()) {
			_total = next;
		} else {
			_total += x;
		}
		_stack.push_back(next);
		return *this;
	}

	// Returns the total statistics stored in this accumulator.
	Stats<R> total() const {
		Stats<R> stats_tot;
		for (Stats<R> const& stats : _stack) {
			stats_tot += stats;
		}
		return stats_tot;
	}

	// Returns a fast estimate of the total statistics in this accumulator.
	Stats<R> total_fast() const {
		return _total;
	}

	// Erases the stored sum.
	void reset() const {
		_stack.clear();
	}
};

// Types associated with `CellBuilder` that are used as callbacks for status
// updates during long-running `explore` and `tune` methods.
template<typename R>
struct ExploreProgress final {
	R progress;
};
template<typename R>
struct ExploreProgressReporter {
	virtual void operator()(ExploreProgress<R> progress) = 0;
};

template<typename R>
struct TuneProgress final {
	R progress;
};
template<typename R>
struct TuneProgressReporter {
	virtual void operator()(TuneProgress<R> progress) = 0;
};

// Parameters that control various aspects of the exploration and tuning
// processes.
template<typename R>
struct CellBuilderParams {
	// The target relative variance that must be reached in all cells.
	R target_rel_var = 0.05;
	// The precision (in standard deviations) to which the relative variance
	// must be reached.
	R target_rel_var_sigma = 1.;
	// The minimum cell volume (with volume equal to length^d) at which
	// exploration terminates.
	R min_cell_length = 0.01;
	// An exponent related to how quickly the efficiency increases with more
	// cell divisions. It depends on the type of tree and the underlying
	// distribution. Empirically, lies between 0.5 and 2.0 for most functions.
	R scale_exp_est = 1.5;
	// Whether to adjust the scaling exponent in response to the behaviour of
	// the distribution at smaller and smaller scales.
	bool scale_exp_adjust = false;
	// Number of past points (in logarithmic space) to consider when adjusting
	// the scaling exponent. For example, a value of 1 means to consider since
	// there were half as many cells, 2 means to consider since there were a
	// quarter as many cells, and so on.
	R scale_exp_history = 1.;
	// Minimum number of points needed to adjust the scaling exponent.
	std::size_t scale_exp_history_count = 8;
	// The maximum cell volatility needed for considering terminating cell
	// exploration in a cell. Generally, this should be around 1. Values larger
	// than 1 prioritize fewer cell divisions over reaching the target
	// efficiency (good for well-behaved distributions), and values smaller than
	// 1 do the opposite.
	R max_cell_split_vol = 1.;
	// The minimum cell volatility at which exploration terminates. A larger
	// value pushes for a more uniform efficiency over the entire generator
	// including in low-prime cells, where a smaller value will cause
	// exploration to focus on the high-prime cells. Should be smaller than the
	// `max_cell_split_vol` to have any effect.
	R min_cell_split_vol = 1.;
	// Relative error in prime difference for a cell division to be triggered.
	R prime_diff_rel_err_for_split = 0.4;
	// A factor related to how much of an improvement to the prime is worth
	// finding the best possible spot to split a cell. Typically should be
	// smaller than one. Smaller values mean better quality splits are chosen,
	// at the cost of more samples per cell.
	R flat_prime_diff_err_factor = 0.5;
	// Number of bins in edge histograms.
	std::size_t hist_num_bins = 128;
	// Minimum number of samples needed per bin before starting cell division.
	std::size_t hist_num_per_bin = 4;
	// Minimum and maximum number of samples to be taken in a single cell during
	// the exploration phase.
	std::size_t min_cell_explore_samples = 512;
	std::size_t max_cell_explore_samples = 1048576;
	// Maximum number of cells to create during exploration.
	std::size_t max_explore_cells = 1048576;
	// The accuracy to which to tune within.
	R tune_rel_accuracy = 0.01;
	// How many stages to tune. More stages gives slightly better tuning
	// performance.
	std::size_t tune_num_stages = 64;
	// Maximum number of samples to be taken total in the tuning phase.
	std::size_t max_tune_samples = 268435456;
	// Number of samples when doing the initial check of the distribution.
	std::size_t check_samples = 16384;
};

template<Dim D, typename R>
class CellGenerator;

template<Dim D, typename R, typename F>
class CellBuilder final {
	template<Dim D1, typename R1, typename F1>
	friend CellBuilder<D1, R1, F1> make_builder(F1 func);
	template<Dim D1, typename R1, typename F1>
	friend CellBuilder<D1, R1, F1> make_builder(
		F1 func,
		ExploreProgressReporter<R>* explore_progress_reporter,
		TuneProgressReporter<R>* tune_progress_reporter);
	template<Dim D1, typename R1, typename F1>
	friend inline CellGenerator<D1, R1> make_generator(
		CellBuilder<D1, R1, F1> const& builder);

	// Left and right random engines.
	using RndL = std::mt19937;
	using RndR = std::mt19937_64;

	static_assert(
		std::is_floating_point<R>::value,
		"R must be a floating point type.");
	static_assert(
		std::numeric_limits<R>::is_iec559,
		"R must satisfy IEC 559.");
	static_assert(
		IsDistribution<D, R, F>::value,
		"F must be a valid distribution type (able to take a Point<D, R> to a "
		"type convertible to R with noexcept).");

public:
	struct BranchData final { };
	struct LeafData final {
		Stats<R> stats;
		LeafData() : stats() { }
		LeafData(Stats<R> const& stats) : stats(stats) { }
	};
	using TreeType = Tree<KdBranch<D, R>, BranchData, LeafData>;

private:
	// Function.
	F _func;
	// Spatial tree storing the statistics.
	TreeType _tree;
	// Random number engine.
	RndL _rnd_left;
	// Whether the distribution function has been checked.
	bool _checked;
	// Persistent state for estimating scaling exponent while building.
	R _scale_exp;
	std::vector<R> _prev_primes;
	std::vector<std::size_t> _prev_leafs;
	// Global statistics.
	R _mean;
	R _mean_err;
	R _prime;
	R _prime_err;
	// TODO: Explain volatilities somewhere in a comment.
	R _split_vol;
	R _tune_vol;

public:
	// Parameters.
	CellBuilderParams<R> par;

	// Progress reporters.
	ExploreProgressReporter<R>* explore_progress_reporter;
	TuneProgressReporter<R>* tune_progress_reporter;

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

	// Estimates the different in the prime value from splitting a parent region
	// into two child regions.
	static R est_sqrt_m2_diff(
			Stats<R> const& stats,
			R lambda_1, Stats<R> const& stats_1,
			R lambda_2, Stats<R> const& stats_2,
			R* err_out=nullptr) {
		R prime = stats.est_sqrt_m2();
		R prime_1 = lambda_1 * stats_1.est_sqrt_m2();
		R prime_2 = lambda_2 * stats_2.est_sqrt_m2();
		R prime_diff = prime_1 + prime_2 - prime;
		if (err_out != nullptr) {
			// The variance of the difference of primes is given to subleading
			// order by this somewhat complicated expression. Essentially, the
			// variance is reduced somewhat because of positive covariance
			// between the original prime and the two new primes. The subleading
			// contribution is needed because the coefficient for the leading
			// contribution goes to zero as the upper and lower distributions
			// become more similar (share a mean/variance).
			R max = stats.max();
			R max_1 = stats_1.max();
			R max_2 = stats_2.max();
			R max_ratio_1, max_ratio_2;
			if (!(max == 0.)) {
				max_ratio_1 = max_1 / max;
				max_ratio_2 = max_2 / max;
			} else {
				max_ratio_1 = 0.5;
				max_ratio_2 = 0.5;
			}

			// Compute the leading order (in `1 / count`) contribution.
			R term = stats.ratio_var2_to_max2_m2();
			R term_1 = max_ratio_1 * max_ratio_1 * stats_1.ratio_var2_to_max2_m2();
			R term_2 = max_ratio_2 * max_ratio_2 * stats_2.ratio_var2_to_max2_m2();
			R m2 = stats.ratio_m2_to_max2();
			R m2_1 = max_ratio_1 * max_ratio_1 * stats_1.ratio_m2_to_max2();
			R m2_2 = max_ratio_2 * max_ratio_2 * stats_2.ratio_m2_to_max2();
			R ratio_1, ratio_2;
			if (!(m2 == 0.)) {
				ratio_1 = m2_1 / m2;
				ratio_2 = m2_2 / m2;
			} else {
				ratio_1 = m2_1 * stats.count();
				ratio_2 = m2_2 * stats.count();
			}
			R factor_1 = lambda_1 * (1. - 2. * std::sqrt(clamp_above<R>(ratio_1)));
			R factor_2 = lambda_2 * (1. - 2. * std::sqrt(clamp_above<R>(ratio_2)));
			R leading = term + factor_1 * term_1 + factor_2 * term_2;

			// Compute the subleading order contribution. This is included
			// because when the means of `stats_1` and `stats_2` are the same,
			// then the leading order term vanishes. The subleading order term
			// also mostly vanishes, but a small contribution remains, which is
			// included here.
			R mixed_var = max_ratio_1 * max_ratio_1 * max_ratio_1 * max_ratio_1
					* lambda_2 * stats_1.ratio_var2_to_max4()
				+ max_ratio_2 * max_ratio_2 * max_ratio_2 * max_ratio_2
					* lambda_1 * stats_2.ratio_var2_to_max4();
			R sum_mean = max_ratio_1 * max_ratio_1 * stats_1.ratio_m2_to_max2()
				+ max_ratio_2 * max_ratio_2 * stats_2.ratio_m2_to_max2();
			R subleading;
			if (!(sum_mean == 0.)) {
				subleading = (mixed_var / sum_mean)
					* (2. - mixed_var / sum_mean / sum_mean);
			} else {
				subleading = 0.;
			}

			// TODO: Consider adding in other contributions to the variance to
			// improve accuracy, especially in the identical means case. For
			// now, we make the questionable decision to limit leading and
			// subleading terms to be positive. This overestimates the variance
			// in extreme cases, but also ensures it is never negative.

			// This variance is calculated relative to the maximum of the
			// original statistics.
			R prime_diff_var =
				0.25 * clamp_above<R>(leading) / stats.count()
				+ 0.25 * clamp_above<R>(subleading) / stats.count() / stats.count();
			*err_out = stats.max() * std::sqrt(clamp_above<R>(prime_diff_var));
		}
		return prime_diff;
	}

	// Calculate the split volatility, a number indicating how much of an impact
	// splitting a cell has on the overall efficiency of the generator.
	static R split_vol(Stats<R> const& stats, R scale_exp_est) {
		R prime = stats.est_sqrt_m2();
		R mean = stats.est_mean();
		R split_vol = pow_inv(prime - mean, scale_exp_est + 1.);
		return split_vol;
	}

	// Calculate the tuning volatility, a number indicating how much of an
	// impact sampling from a cell has on the overall efficiency of the
	// generator.
	static R tune_vol(Stats<R> const& stats, R prime_tot) {
		R prime = stats.est_sqrt_m2();
		R ratio = stats.ratio_var2_to_m2p2();
		R tune_vol = 0.5 * std::sqrt(clamp_above<R>(
			(prime / prime_tot) * ((prime_tot - prime) / prime_tot) * ratio));
		return tune_vol;
	}

	// Finds the mean (with error) of a cell.
	R mean(CellHandle cell, R* err_out=nullptr) {
		if (_tree.type(cell) == CellType::BRANCH) {
			// Store the means and errors from all children.
			R means[TreeType::NUM_CHILDREN];
			R mean_errs[TreeType::NUM_CHILDREN];
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(cell, child_idx);
				means[child_idx] = mean(child, &mean_errs[child_idx]);
			}
			// Combine the means and errors.
			R mean_tot = 0.;
			for (R mean : means) {
				mean_tot += mean;
			}
			if (err_out != nullptr) {
				// Prevent overflow when combining errors in quadrature.
				R mean_err_sq_tot = 0.;
				R norm_err = 0.;
				for (R mean_err : mean_errs) {
					norm_err = std::max(norm_err, mean_err);
				}
				if (norm_err != 0.) {
					for (R mean_err : mean_errs) {
						mean_err_sq_tot += sq(mean_err / norm_err);
					}
				}
				R mean_err_tot = norm_err * std::sqrt(mean_err_sq_tot);
				*err_out = mean_err_tot;
			}
			return mean_tot;
		} else {
			return _tree.leaf_data(cell).stats.est_mean(err_out);
		}
	}
	// Finds the prime (with error) for a cell.
	R prime(CellHandle cell, R* err_out=nullptr) {
		if (_tree.type(cell) == CellType::BRANCH) {
			R primes[TreeType::NUM_CHILDREN];
			R prime_errs[TreeType::NUM_CHILDREN];
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(cell, child_idx);
				primes[child_idx] = prime(child, &prime_errs[child_idx]);
			}
			R prime_tot = 0.;
			for (R prime : primes) {
				prime_tot += prime;
			}
			if (err_out != nullptr) {
				// Prevent overflow when combining errors in quadrature.
				R prime_err_sq_tot = 0.;
				R norm_err = 0.;
				for (R prime_err : prime_errs) {
					norm_err = std::max(norm_err, prime_err);
				}
				if (norm_err != 0.) {
					for (R prime_err : prime_errs) {
						prime_err_sq_tot += sq(prime_err / norm_err);
					}
				}
				R prime_err_tot = norm_err * std::sqrt(prime_err_sq_tot);
				*err_out = prime_err_tot;
			}
			return prime_tot;
		} else {
			return _tree.leaf_data(cell).stats.est_sqrt_m2(err_out);
		}
	}
	// Finds split volatility for a cell.
	R split_vol(CellHandle cell, R scale_exp) {
		if (_tree.type(cell) == CellType::BRANCH) {
			R split_vol_tot = 0.;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(cell, child_idx);
				split_vol_tot += split_vol(child, scale_exp);
			}
			return split_vol_tot;
		} else {
			return split_vol(_tree.leaf_data(cell).stats, scale_exp);
		}
	}
	// Finds tuning volatility for a cell.
	R tune_vol(CellHandle cell, R prime_tot) {
		if (_tree.type(cell) == CellType::BRANCH) {
			R tune_vol_tot = 0.;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(cell, child_idx);
				tune_vol_tot += tune_vol(child, prime_tot);
			}
			return tune_vol_tot;
		} else {
			return tune_vol(_tree.leaf_data(cell).stats, prime_tot);
		}
	}

	// Fills in the global statistics measures.
	void update_total_stats() {
		CellHandle root = _tree.root();
		_mean = mean(root, &_mean_err);
		_prime = prime(root, &_prime_err);
		_split_vol = split_vol(root, _scale_exp);
		_tune_vol = tune_vol(root, _prime);
	}

	void explore_leaf_flat(
			CellHandle cell,
			Point<D, R> offset, Point<D, R> extent) noexcept {
		// Create a random number generator (a "right" generator) which will be
		// used.
		std::vector<Seed> seed;
		#ifdef BUBBLE_USE_OPENMP
		#pragma omp critical (bubble_rand)
		#endif
		{
			seed = sample_rnd_left();
		}
		std::seed_seq seed_seq(seed.begin(), seed.end());
		RndR rnd(seed_seq);

		R volume = 1.;
		for (Dim dim = 0; dim < D; ++dim) {
			volume *= extent[dim];
		}
		StatsAccum<R> stats_accum;
		std::array<std::vector<StatsAccum<R> >, D> hist_stats_accum;
		for (Dim dim = 0; dim < D; ++dim) {
			hist_stats_accum[dim].resize(par.hist_num_bins);
		}
		for (std::size_t sample = 0; sample < par.min_cell_explore_samples; ++sample) {
			Point<D, R> point;
			std::size_t hist_idx[D];
			std::uniform_int_distribution<std::size_t> idx_dist(
				0,
				par.hist_num_bins - 1);
			for (Dim dim = 0; dim < D; ++dim) {
				std::size_t idx = idx_dist(rnd);
				hist_idx[dim] = idx;
				R lower = offset[dim]
					+ extent[dim] * idx / par.hist_num_bins;
				R upper = offset[dim]
					+ extent[dim] * (idx + 1) / par.hist_num_bins;
				std::uniform_real_distribution<R> x_dist(lower, upper);
				point[dim] = x_dist(rnd);
			}
			// First choose a histogram index to generate in, then make the
			// point itself.
			R f = volume * _func(point);
			stats_accum += f;
			// Fill histograms.
			for (Dim dim = 0; dim < D; ++dim) {
				hist_stats_accum[dim][hist_idx[dim]] += f;
			}
		}
		Stats<R> stats = stats_accum.total();
		// Store the stats for a half-way division along the
		// dimension with greatest extent, just in case the cell is
		// too flat to choose a division site in any better way.
		Stats<R> stats_lower_half;
		Stats<R> stats_upper_half;
		// Find the dimension with the greatest extent.
		Dim dim_extent_max = 0;
		R extent_max = 0.;
		for (Dim dim = 0; dim < D; ++dim) {
			if (extent[dim] > extent_max) {
				extent_max = extent[dim];
				dim_extent_max = dim;
			}
		}
		// Find the division site with the greatest reduction in
		// overall prime value.
		R prime_diff_min = std::numeric_limits<R>::max();
		R prime_diff_err_min = 0.;
		Dim dim_min = dim_extent_max;
		std::size_t bin_min = par.hist_num_bins / 2;
		Stats<R> stats_lower_min;
		Stats<R> stats_upper_min;
		for (Dim dim = 0; dim < D; ++dim) {
			// Total up the statistics in the histograms.
			std::vector<Stats<R> > hist_stats(par.hist_num_bins);
			for (std::size_t bin = 0; bin < par.hist_num_bins; ++bin) {
				hist_stats[bin] = hist_stats_accum[dim][bin].total();
			}
			// Integrated statistics below and above the proposed
			// division site.
			std::vector<Stats<R> > hist_stats_lower;
			std::vector<Stats<R> > hist_stats_upper;
			hist_stats_lower.reserve(par.hist_num_bins - 1);
			hist_stats_upper.reserve(par.hist_num_bins - 1);
			// TODO: Use accumulators? The difficulty is that the
			// total is needed every iteration to push onto the
			// integrated histograms.
			Stats<R> stats_lower_tot;
			Stats<R> stats_upper_tot;
			for (std::size_t bin = 0; bin < par.hist_num_bins - 1; ++bin) {
				std::size_t rbin = par.hist_num_bins - bin - 1;
				stats_lower_tot += hist_stats[bin];
				stats_upper_tot += hist_stats[rbin];
				hist_stats_lower.push_back(stats_lower_tot);
				hist_stats_upper.push_back(stats_upper_tot);
			}
			// Store the statistics for a half-way division.
			if (dim == dim_extent_max) {
				std::size_t bin = par.hist_num_bins / 2;
				std::size_t bin_lower = bin - 1;
				std::size_t bin_upper = par.hist_num_bins - bin - 1;
				stats_lower_half = hist_stats_lower[bin_lower];
				stats_upper_half = hist_stats_upper[bin_upper];
			}
			// Test division at each internal bin boundary.
			for (std::size_t bin = 1; bin < par.hist_num_bins; ++bin) {
				std::size_t bin_lower = bin - 1;
				std::size_t bin_upper = par.hist_num_bins - bin - 1;
				Stats<R> stats_lower = hist_stats_lower[bin_lower];
				Stats<R> stats_upper = hist_stats_upper[bin_upper];
				// If a bin doesn't have the counts needed to
				// compute statistics, just drop it.
				if (stats_lower.count() <= 3
						|| stats_upper.count() <= 3) {
					continue;
				}
				// Estimate the total prime resulting from splitting
				// the cell at this bin boundary. Need to scale by
				// the volume fractions, to account for volume being
				// included in the statistics measures.
				R lambda_lower = R(bin) / par.hist_num_bins;
				R lambda_upper = R(par.hist_num_bins - bin) / par.hist_num_bins;
				R prime_diff_err;
				R prime_diff = est_sqrt_m2_diff(
					stats,
					lambda_lower, stats_lower,
					lambda_upper, stats_upper,
					&prime_diff_err);
				// Update minimum.
				if (prime_diff < prime_diff_min) {
					prime_diff_min = prime_diff;
					prime_diff_err_min = prime_diff_err;
					dim_min = dim;
					bin_min = bin;
					stats_lower_min = stats_lower;
					stats_upper_min = stats_upper;
				}
			}
		}
		// If the relative error is small enough, or if just too
		// many samples have been requested, then make a division.
		// TODO: Include a user-defined factor here.
		if (std::abs(prime_diff_min) <= prime_diff_err_min) {
			dim_min = dim_extent_max;
			bin_min = par.hist_num_bins / 2;
			stats_lower_min = stats_lower_half;
			stats_upper_min = stats_upper_half;
		}
		R lambda_lower = R(bin_min) / par.hist_num_bins;
		R lambda_upper = R(par.hist_num_bins - bin_min) / par.hist_num_bins;
		LeafData leaf_data[2] = {
			{ lambda_lower * stats_lower_min },
			{ lambda_upper * stats_upper_min },
		};
		#ifdef BUBBLE_USE_OPENMP
		#pragma omp critical (bubble_tree)
		#endif
		{
			_tree.split(
				cell,
				typename TreeType::Branch(
					dim_min, lambda_lower),
				BranchData(),
				leaf_data);
		}
	}
	// Determines the efficiency of the cell generator within a specific cell.
	// If the efficiency is good enough, leaves the cell alone, otherwise, adds
	// new children to the cell to improve its efficiency.
	void explore_leaf(
			// The cell to explore.
			CellHandle cell,
			// The initial statistics stored in the cell.
			Stats<R> stats_init,
			// Dimensions of the cell.
			Point<D, R> offset, Point<D, R> extent,
			// Scaling exponent.
			R scale_exp,
			// Number of leafs.
			std::size_t leafs,
			// Estimate of total integral of distribution.
			R mean_tot,
			// Estimate of split volatility of tree.
			R split_vol_tot) noexcept {
		R volume = 1.;
		for (Dim dim = 0; dim < D; ++dim) {
			volume *= extent[dim];
		}
		StatsAccum<R> stats_accum;
		// Histograms, one for each axis of the hypercube.
		std::array<std::vector<StatsAccum<R> >, D> hist_stats_accum;
		for (Dim dim = 0; dim < D; ++dim) {
			hist_stats_accum[dim].resize(par.hist_num_bins);
		}
		// Find the dimension with the greatest extent.
		Dim dim_extent_max = 0;
		R extent_max = 0.;
		for (Dim dim = 0; dim < D; ++dim) {
			if (extent[dim] > extent_max) {
				extent_max = extent[dim];
				dim_extent_max = dim;
			}
		}

		// Create a random number generator (a "right" generator) which will be
		// used.
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
		std::size_t next_samples = par.min_cell_explore_samples;
		while (true) {
			// Check termination conditions.
			Stats<R> stats = stats_accum.total();
			Stats<R> stats_tot = stats + stats_init;
			if (stats_tot.count() >= par.min_cell_explore_samples) {
				R rel_var_err;
				R rel_var = stats_tot.est_rel_var(&rel_var_err);
				// The target prime (given directly by the target relative
				// variance), with the mean subtracted from it.
				R target_prime_diff = mean_tot * sqrt1p_1m(par.target_rel_var);

				// Termination condition 1. Efficiency meets criteria.
				// First check the cell division volatility. If the ratio is
				// greater than one, then check that the relative variance has
				// met the target accounting for uncertainty (and if so, then
				// terminate).
				R cell_split_vol =
					pow_inv(
						clamp_above<R>(split_vol_tot / target_prime_diff),
						scale_exp)
					* split_vol(stats_tot, scale_exp);
				// Termination will only be considered if the cell volatility
				// is below the threshold value. This ensures that no matter
				// what other termination conditions are considered, the final
				// cell generator will have good efficiency. Note that we
				// additionally require that the total mean is non-zero. If the
				// total mean is zero, then all cells are required to divide in
				// hopes of measuring a non-zero mean eventually.
				if (mean_tot != 0. && cell_split_vol < par.max_cell_split_vol) {
					R rel_var_max =
						rel_var + par.target_rel_var_sigma * rel_var_err;
					// Three possibilities for terminating:
					// * Relative variance is smaller than target.
					// * Cell is unimportant either because:
					//   * Cell is low volume.
					//   * Cell has low volatility.
					bool var_cond = (rel_var_max <= par.target_rel_var);
					bool volume_cond = (volume <= std::pow(par.min_cell_length, D));
					bool vol_cond = (cell_split_vol <= par.min_cell_split_vol);
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
				if (stats.count() >= par.hist_num_per_bin * par.hist_num_bins) {
					// Store the stats for a half-way division along the
					// dimension with greatest extent, just in case the cell is
					// too flat to choose a division site in any better way.
					Stats<R> stats_lower_half;
					Stats<R> stats_upper_half;
					// Find the division site with the greatest reduction in
					// overall prime value.
					R prime_diff_min = std::numeric_limits<R>::max();
					R prime_diff_err_min = 0.;
					Dim dim_min = dim_extent_max;
					std::size_t bin_min = par.hist_num_bins / 2;
					Stats<R> stats_lower_min;
					Stats<R> stats_upper_min;
					for (Dim dim = 0; dim < D; ++dim) {
						// Total up the statistics in the histograms.
						std::vector<Stats<R> > hist_stats(par.hist_num_bins);
						for (std::size_t bin = 0; bin < par.hist_num_bins; ++bin) {
							hist_stats[bin] = hist_stats_accum[dim][bin].total();
						}
						// Integrated statistics below and above the proposed
						// division site.
						std::vector<Stats<R> > hist_stats_lower;
						std::vector<Stats<R> > hist_stats_upper;
						hist_stats_lower.reserve(par.hist_num_bins - 1);
						hist_stats_upper.reserve(par.hist_num_bins - 1);
						// TODO: Use accumulators? The difficulty is that the
						// total is needed every iteration to push onto the
						// integrated histograms.
						Stats<R> stats_lower_tot;
						Stats<R> stats_upper_tot;
						for (std::size_t bin = 0; bin < par.hist_num_bins - 1; ++bin) {
							std::size_t rbin = par.hist_num_bins - bin - 1;
							stats_lower_tot += hist_stats[bin];
							stats_upper_tot += hist_stats[rbin];
							hist_stats_lower.push_back(stats_lower_tot);
							hist_stats_upper.push_back(stats_upper_tot);
						}
						// Store the statistics for a half-way division.
						if (dim == dim_extent_max) {
							std::size_t bin = par.hist_num_bins / 2;
							std::size_t bin_lower = bin - 1;
							std::size_t bin_upper = par.hist_num_bins - bin - 1;
							stats_lower_half = hist_stats_lower[bin_lower];
							stats_upper_half = hist_stats_upper[bin_upper];
						}
						// Test division at each internal bin boundary.
						for (std::size_t bin = 1; bin < par.hist_num_bins; ++bin) {
							std::size_t bin_lower = bin - 1;
							std::size_t bin_upper = par.hist_num_bins - bin - 1;
							Stats<R> stats_lower = hist_stats_lower[bin_lower];
							Stats<R> stats_upper = hist_stats_upper[bin_upper];
							// If a bin doesn't have the counts needed to
							// compute statistics, just drop it.
							if (stats_lower.count() <= 3
									|| stats_upper.count() <= 3) {
								continue;
							}
							// Estimate the total prime resulting from splitting
							// the cell at this bin boundary. Need to scale by
							// the volume fractions, to account for volume being
							// included in the statistics measures.
							R lambda_lower =
								R(bin) / par.hist_num_bins;
							R lambda_upper =
								R(par.hist_num_bins - bin) / par.hist_num_bins;
							R prime_diff_err;
							R prime_diff = est_sqrt_m2_diff(
								stats,
								lambda_lower, stats_lower,
								lambda_upper, stats_upper,
								&prime_diff_err);
							// Update minimum.
							if (prime_diff < prime_diff_min) {
								prime_diff_min = prime_diff;
								prime_diff_err_min = prime_diff_err;
								dim_min = dim;
								bin_min = bin;
								stats_lower_min = stats_lower;
								stats_upper_min = stats_upper;
							}
						}
					}
					// If the relative error is small enough, or if just too
					// many samples have been requested, then make a division.
					R prime_diff_err_req = std::abs(
						par.prime_diff_rel_err_for_split * prime_diff_min);
					R flat_prime_diff_err_req = par.flat_prime_diff_err_factor
						* target_prime_diff / std::sqrt(leafs);
					bool err_cond = (prime_diff_err_min <= prime_diff_err_req);
					bool sample_cond = (stats.count() >= par.max_cell_explore_samples);
					bool flat_cond = (flat_prime_diff_err_req >= -prime_diff_min);
					if (err_cond || sample_cond || flat_cond) {
						// TODO: Include a user-defined factor here.
						if (std::abs(prime_diff_min) <= prime_diff_err_min) {
							dim_min = dim_extent_max;
							bin_min = par.hist_num_bins / 2;
							stats_lower_min = stats_lower_half;
							stats_upper_min = stats_upper_half;
						}
						R lambda_lower =
							R(bin_min) / par.hist_num_bins;
						R lambda_upper =
							R(par.hist_num_bins - bin_min) / par.hist_num_bins;
						LeafData leaf_data[2] = {
							{ lambda_lower * stats_lower_min },
							{ lambda_upper * stats_upper_min },
						};
						#ifdef BUBBLE_USE_OPENMP
						#pragma omp critical (bubble_tree)
						#endif
						{
							_tree.split(
								cell,
								typename TreeType::Branch(
									dim_min, lambda_lower),
								BranchData(),
								leaf_data);
						}
						return;
					}
				}
			}

			// If non of the termination conditions match, then produce more
			// samples.
			for (std::size_t sample = stats.count(); sample < next_samples; ++sample) {
				// Generate a point randomly within the cell and evaluate the
				// function. First choose a histogram bin to generate in, then
				// produce the point itself.
				Point<D, R> point;
				std::size_t hist_idx[D];
				std::uniform_int_distribution<std::size_t> idx_dist(
					0,
					par.hist_num_bins - 1);
				for (Dim dim = 0; dim < D; ++dim) {
					std::size_t idx = idx_dist(rnd);
					hist_idx[dim] = idx;
					R lower = offset[dim]
						+ extent[dim] * idx / par.hist_num_bins;
					R upper = offset[dim]
						+ extent[dim] * (idx + 1) / par.hist_num_bins;
					std::uniform_real_distribution<R> x_dist(lower, upper);
					point[dim] = x_dist(rnd);
				}
				// First choose a histogram index to generate in, then make the
				// point itself.
				R f = volume * _func(point);
				stats_accum += f;
				// Fill histograms.
				for (Dim dim = 0; dim < D; ++dim) {
					hist_stats_accum[dim][hist_idx[dim]] += f;
				}
			}
			// Double how many samples we produce each round.
			next_samples *= 2;
		}
	}
	// Recursively applies `explore_leaf` to any leaf cells descended from the
	// provided parent cell.
	void explore_cell(
			CellHandle cell,
			Point<D, R> offset, Point<D, R> extent,
			R scale_exp,
			std::size_t leafs,
			R mean_tot,
			R split_vol_tot) noexcept {
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
					scale_exp, leafs, mean_tot, split_vol_tot);
			}
		} else {
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp task
			#endif
			{
				explore_leaf(
					cell,
					stats_init,
					offset, extent,
					scale_exp, leafs, mean_tot, split_vol_tot);
			}
		}
	}

	// Improves the estimate of the prime of a given cell.
	void tune_leaf(
			// The cell to tune.
			CellHandle cell,
			// Dimensions of the cell.
			Point<D, R> offset, Point<D, R> extent,
			// Total prime of the generator.
			R prime_tot,
			// Total tuning volatility of the generator.
			R tune_vol_tot,
			// Number of samples to use for estimating the prime.
			std::size_t samples) noexcept {
		R volume = 1.;
		for (Dim dim = 0; dim < D; ++dim) {
			volume *= extent[dim];
		}
		Stats<R>& stats = _tree.leaf_data(cell).stats;
		// Sample by the fraction of the total volatility stored in this cell.
		R tune_vol_fraction = clamp<R>(
			tune_vol(stats, prime_tot) / tune_vol_tot,
			0., 1.);
		samples *= tune_vol_fraction;

		if (samples <= stats.count()) {
			return;
		}
		samples -= stats.count();

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

		StatsAccum<R> stats_accum;
		for (std::size_t sample = 0; sample < samples; ++sample) {
			Point<D, R> point;
			for (Dim dim = 0; dim < D; ++dim) {
				R lower = offset[dim];
				R upper = offset[dim] + extent[dim];
				std::uniform_real_distribution<R> x_dist(lower, upper);
				point[dim] = x_dist(rnd);
			}
			stats_accum += volume * _func(point);
		}
		stats += stats_accum.total();
	}
	// Recursively applies `tune_leaf` to any leaf cells descended from the
	// provided parent cell.
	void tune_cell(
			CellHandle cell,
			Point<D, R> offset, Point<D, R> extent,
			R prime_tot,
			R tune_vol_tot,
			std::size_t samples) noexcept {
		CellType type;
		Point<D, R> offset_local[TreeType::NUM_CHILDREN];
		Point<D, R> extent_local[TreeType::NUM_CHILDREN];
		CellHandle child[TreeType::NUM_CHILDREN];
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
			}
		}
		if (type == CellType::BRANCH) {
			// Sample children according to volatility.
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				Point<D, R> offset_child;
				Point<D, R> extent_child;
				for (Dim dim = 0; dim < D; ++dim) {
					offset_child[dim] = offset[dim];
					offset_child[dim] += extent[dim] * offset_local[child_idx][dim];
					extent_child[dim] = extent[dim] * extent_local[child_idx][dim];
				}
				tune_cell(
					child[child_idx],
					offset_child, extent_child,
					prime_tot, tune_vol_tot,
					samples);
			}
		} else {
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp task
			#endif
			{
				tune_leaf(
					cell,
					offset, extent,
					prime_tot, tune_vol_tot,
					samples);
			}
		}
	}

public:
	CellBuilder(F func) :
		CellBuilder(func, nullptr, nullptr) { }
	CellBuilder(
			F func,
			ExploreProgressReporter<R>* explore_progress_reporter,
			TuneProgressReporter<R>* tune_progress_reporter) :
			_func(func),
			_tree(LeafData()),
			_rnd_left(),
			_checked(false),
			explore_progress_reporter(explore_progress_reporter),
			tune_progress_reporter(tune_progress_reporter) {
		// Seed the left generator.
		std::seed_seq seed_seq { std::random_device()() };
		_rnd_left.seed(seed_seq);
	}

	F const& func() const {
		return _func;
	}

	TreeType const& tree() const {
		return _tree;
	}

	// Validates the distribution function by sampling some number of points
	// uniformly over the entire unit hypercube.
	void check() {
		if (_checked) {
			return;
		}
		_checked = true;
		// Sample a new generator for initialization.
		std::vector<Seed> seed_right = sample_rnd_left();
		std::seed_seq seed_seq_right(seed_right.begin(), seed_right.end());
		RndR rnd(seed_seq_right);
		// Do some sampling of the distribution before starting. Doing this
		// sampling gives an initial mean and prime estimate, and also lets us
		// catch some common errors.
		CellHandle root = _tree.root();
		StatsAccum<R> stats_accum;
		for (std::size_t sample = 0; sample < par.check_samples; ++sample) {
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
			stats_accum += _func(point);
		}
		Stats<R> stats = stats_accum.total();
		if (!stats.valid()) {
			throw std::runtime_error(
				"distribution invalid due to potential overflow or invalid "
				"floating point values");
		}
		if (stats.underflow()) {
			throw std::runtime_error(
				"distribution invalid due to underflow");
		}
		if (stats.mean() <= 0) {
			throw std::runtime_error(
				"distribution invalid because cannot be distinguished from "
				"zero");
		}
		// If the root cell is the only cell, then insert the stats into it.
		// This should always be the case, since the `explore` function calls
		// `check` before starting to divide up the root cell, but this check is
		// here just in case.
		if (_tree.type(root) == CellType::LEAF) {
			_tree.leaf_data(root).stats = stats;
			update_total_stats();
		}
	}

	void explore_flat() {
		if (!_checked) {
			check();
		}
		struct CellHandleSized {
			CellHandle handle;
			Point<D, R> offset;
			Point<D, R> extent;
		};
		auto compare = [&](CellHandleSized a, CellHandleSized b) {
			Stats<R> stats_a = _tree.leaf_data(a.handle).stats;
			Stats<R> stats_b = _tree.leaf_data(b.handle).stats;
			R prime_a = stats_a.est_sqrt_m2();
			R mean_a = stats_a.est_mean();
			R prime_b = stats_b.est_sqrt_m2();
			R mean_b = stats_b.est_mean();
			return prime_a - mean_a < prime_b - mean_b;
		};
		std::priority_queue<
			CellHandleSized,
			std::vector<CellHandleSized>,
			decltype(compare)> queue(compare);
		Point<D, R> offset;
		Point<D, R> extent;
		offset.fill(0.);
		extent.fill(1.);
		queue.push({ _tree.root(), offset, extent });
		std::vector<CellHandleSized> next_cells;
		while (_tree.size() < par.max_explore_cells) {
			next_cells.clear();
			for (std::size_t idx = 0; idx < 32; ++idx) {
				if (queue.empty()) {
					break;
				}
				CellHandleSized cell = queue.top();
				queue.pop();
				next_cells.push_back(cell);
			}
			if (explore_progress_reporter != nullptr) {
				R ratio = static_cast<R>(_tree.size()) / par.max_explore_cells;
				R progress = clamp<R>(ratio, 0., 1.);
				(*explore_progress_reporter)({ progress });
			}
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp parallel for shared(queue)
			#endif
			for (std::size_t idx = 0; idx < next_cells.size(); ++idx) {
				CellHandleSized cell = next_cells[idx];
				explore_leaf_flat(cell.handle, cell.offset, cell.extent);
				#ifdef BUBBLE_USE_OPENMP
				#pragma omp critical (bubble_queue)
				#pragma omp critical (bubble_tree)
				#endif
				{
					for (std::size_t child_idx = 0; child_idx < 2; ++child_idx) {
						Point<D, R> offset_child = cell.offset;
						Point<D, R> extent_child = cell.extent;
						Point<D, R> offset_local = _tree.offset_local(cell.handle, child_idx);
						Point<D, R> extent_local = _tree.extent_local(cell.handle, child_idx);
						for (Dim dim = 0; dim < D; ++dim) {
							offset_child[dim] += extent[dim] * offset_local[dim];
							extent_child[dim] *= extent_local[dim];
						}
						queue.push({ _tree.child(cell.handle, child_idx), offset_child, extent_child });
					}
				}
			}
		}
	}

	// Explores the distribution to create a well-balanced cell generator that
	// can reach the desired performance. If unable to create a well-balanced
	// cell generator, returns estimated number of cells that the exploration
	// process fell short by.
	std::size_t explore() {
		// Make sure that we have some basic validation of the distribution.
		if (!_checked) {
			check();
		}
		Point<D, R> offset;
		Point<D, R> extent;
		offset.fill(0.);
		extent.fill(1.);
		CellHandle root = _tree.root();
		std::size_t cells;
		_scale_exp = par.scale_exp_est;
		do {
			cells = _tree.size();
			R mean_tot = mean(root);
			R prime_tot = prime(root);
			R split_vol_tot = split_vol(root, _scale_exp);
			R target_prime_diff = mean_tot * sqrt1p_1m(par.target_rel_var);
			std::size_t leafs = _tree.leaf_size();
			// TODO: I've been using this to log behaviour of the exploration
			// process, so I'll shamelessly leave it in here until proper
			// progress evaluation is in place.
			//std::cout << "Status: " << leafs << "\t" << prime_tot << "\t" << mean_tot << std::endl;
			_prev_primes.push_back(prime_tot);
			_prev_leafs.push_back(leafs);
			// Estimate the true scaling exponent using a linear regression to
			// some number of previous steps.
			R scale_exp_regression = std::numeric_limits<R>::quiet_NaN();
			if (_prev_primes.size() >= par.scale_exp_history_count) {
				R mean_log_cells = 0.;
				R mean_log_prime = 0.;
				R mean_log_cells_sq = 0.;
				R mean_log_prime_cells = 0.;
				std::size_t idx_start = 0;
				for (std::size_t idx = _prev_leafs.size(); idx-- > 0; ) {
					R history = std::log2(R(leafs) / _prev_leafs[idx]);
					std::size_t history_count = _prev_leafs.size() - idx;
					if (history > par.scale_exp_history
							&& history_count >= par.scale_exp_history_count) {
						idx_start = idx;
						break;
					}
				}
				std::size_t count = _prev_primes.size() - idx_start;
				for (std::size_t idx = idx_start; idx < _prev_primes.size(); ++idx) {
					auto log_leafs = std::log(_prev_leafs[idx]);
					R log_prime = std::log(_prev_primes[idx] - mean_tot);
					mean_log_cells += log_leafs / count;
					mean_log_prime += log_prime / count;
					mean_log_cells_sq += sq(log_leafs) / count;
					mean_log_prime_cells += log_prime * log_leafs / count;
				}
				scale_exp_regression =
					-(mean_log_prime_cells - mean_log_prime * mean_log_cells)
					/ (mean_log_cells_sq - sq(mean_log_cells));
			}
			if (par.scale_exp_adjust) {
				if (!std::isfinite(scale_exp_regression)) {
					_scale_exp = par.scale_exp_est;
				} else if (!(scale_exp_regression >= 0.)) {
					_scale_exp = 0.;
				} else {
					_scale_exp = scale_exp_regression;
				}
			} else {
				_scale_exp = par.scale_exp_est;
			}
			// Estimate progress using regression scale exponent.
			if (explore_progress_reporter != nullptr) {
				R prime_diff = prime_tot - mean_tot;
				R progress = clamp<R>(
					pow_inv(target_prime_diff / prime_diff, scale_exp_regression),
					0., 1.);
				(*explore_progress_reporter)({ progress });
			}
			// TODO: I've been using this to log behaviour of the exploration
			// process, so I'll shamelessly leave it in here until proper
			// progress evaluation is in place.
			//std::cout << "Status: " << leafs << ", " << prime_tot << ", " << mean_tot << ", " << scale_exp_regression << std::endl;
			// Explore the cell generator.
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp parallel
			#pragma omp single
			#endif
			{
				explore_cell(
					root,
					offset, extent,
					_scale_exp, leafs, mean_tot, split_vol_tot);
			}
			// TODO: Termination conditions should be revised.
		} while (cells != _tree.size() && _tree.size() < par.max_explore_cells);
		update_total_stats();
		// Calculate ideal number of leaf cells, and see if we surpassed it.
		R target_prime_diff = _mean * sqrt1p_1m(par.target_rel_var);
		R ideal_leafs = pow_inv(
			std::pow(_split_vol, _scale_exp + 1.) / target_prime_diff,
			_scale_exp);
		if (ideal_leafs > R(_tree.leaf_size())) {
			return ideal_leafs - _tree.leaf_size();
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
		for (std::size_t stage = 0; stage < par.tune_num_stages; ++stage) {
			R mean_tot = mean(root);
			R prime_tot = prime(root);
			R tune_vol_tot = tune_vol(root, prime_tot);
			// Find the number of samples needed to reach the desired accuracy.
			R rel_var = sq(prime_tot / mean_tot) - 1.;
			R tune_accuracy = rel_var * par.tune_rel_accuracy;
			R tune_fraction = clamp<R>(R(stage + 1) / par.tune_num_stages, 0., 1.);
			R samples_real = clamp<R>(
				0.5 * tune_fraction * sq(tune_vol_tot) / tune_accuracy,
				// Multiply by 0.5 here so that the process of converting to `R`
				// won't cause `samples_real` to overflow the maximum possible
				// value.
				0.,
				0.5 * static_cast<R>(std::numeric_limits<std::size_t>::max()));
			std::size_t samples = static_cast<std::size_t>(samples_real);
			std::size_t max_samples = static_cast<std::size_t>(
				tune_fraction * par.max_tune_samples);
			if (samples > max_samples) {
				undertuned_samples = samples - par.max_tune_samples;
				samples = par.max_tune_samples;
			}
			// Estimate progress using the tuning stage.
			if (tune_progress_reporter != nullptr) {
				R progress = clamp<R>(R(stage) / par.tune_num_stages, 0., 1.);
				(*explore_progress_reporter)({ progress });
			}
			// Sample from the cell generator.
			#ifdef BUBBLE_USE_OPENMP
			#pragma omp parallel
			#pragma omp single
			#endif
			{
				tune_cell(
					root,
					offset, extent,
					prime_tot, tune_vol_tot,
					samples);
			}
		}
		update_total_stats();
		return undertuned_samples;
	}

	R mean(R* err_out=nullptr) const {
		if (err_out != nullptr) {
			*err_out = _mean_err;
		}
		return _mean;
	}
	R prime(R* err_out=nullptr) const {
		if (err_out != nullptr) {
			*err_out = _prime_err;
		}
		return _prime;
	}
	R rel_var(R* err_out=nullptr) const {
		R rel_var_tot = sq(_prime / _mean) - 1.;
		if (err_out != nullptr) {
			// TODO: This error estimate isn't very good because the total prime
			// and total mean values are strongly correlated with each other. In
			// general, though, it should overestimate the error which is
			// acceptable.
			*err_out = rel_var_tot * std::sqrt(
				sq(_mean_err / _mean) + sq(_prime_err / _prime));
		}
		return rel_var_tot;
	}
};

// Convenience functions to allow for template parameter inference.
template<Dim D, typename R, typename F>
inline CellBuilder<D, R, F> make_builder(F func) {
	return CellBuilder<D, R, F>(func);
}
template<Dim D, typename R, typename F>
inline CellBuilder<D, R, F> make_builder(
		F func,
		ExploreProgressReporter<R>* explore_progress_reporter,
		TuneProgressReporter<R>* tune_progress_reporter) {
	return CellBuilder<D, R, F>(
		func,
		explore_progress_reporter,
		tune_progress_reporter);
}

template<Dim D, typename R>
class CellGenerator final {
	template<Dim D1, typename R1>
	friend inline CellGenerator<D1, R1> make_generator();
	template<Dim D1, typename R1, typename F1>
	friend inline CellGenerator<D1, R1> make_generator(
		CellBuilder<D1, R1, F1> const& builder);

	static_assert(
		std::is_floating_point<R>::value,
		"R must be a floating point type.");
	static_assert(
		std::numeric_limits<R>::is_iec559,
		"R must satisfy IEC 559.");

public:
	// The cells store the prime value.
	struct Data final {
		R prime;
		Data() : prime(0.) { }
		Data(R prime) : prime(prime) { }
	};
	using TreeType = Tree<KdBranch<D, R>, Data>;

private:
	// Spatial tree.
	TreeType _tree;

	// Fills in the prime value for a cell and all of its descendents,
	// recursively. This is done on load so that parent prime values don't have
	// to be stored on disk. Returns the total prime value on completion.
	R fill_primes(CellHandle cell) {
		if (_tree.type(cell) == CellType::BRANCH) {
			R prime_tot = 0.;
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(cell, child_idx);
				prime_tot += fill_primes(child);
			}
			_tree.branch_data(cell).prime = prime_tot;
			return prime_tot;
		} else {
			return _tree.leaf_data(cell).prime;
		}
	}

	void select_cell(
			CellHandle cell,
			Point<D, R> offset, Point<D, R> extent,
			R choose_cell, Point<D, R> choose_point,
			R* weight_out, Point<D, R>* point_out) const {
		if (_tree.type(cell) == CellType::BRANCH) {
			for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
				CellHandle child = _tree.child(cell, child_idx);
				// Each child has a chance of being selected proportional to its
				// prime value.
				R prime_ratio = _tree.leaf_data(child).prime
					/ _tree.branch_data(cell).prime;
				if (choose_cell <= prime_ratio
						|| child_idx == TreeType::NUM_CHILDREN - 1) {
					// Scale `choose_cell` to be between 0 and 1 again.
					choose_cell /= prime_ratio;
					// Adjust offset and extent for the next level of recursion.
					Point<D, R> offset_local = _tree.offset_local(cell, child_idx);
					Point<D, R> extent_local = _tree.extent_local(cell, child_idx);
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
		} else if (_tree.type(cell) == CellType::LEAF) {
			// Rescale by the jacobian.
			R volume = 1;
			for (Dim dim = 0; dim < D; ++dim) {
				(*point_out)[dim] = offset[dim];
				(*point_out)[dim] += extent[dim] * choose_point[dim];
				volume *= extent[dim];
			}
			R weight_norm = _tree.leaf_data(cell).prime;
			*weight_out = volume / weight_norm;
		}
	}

	// Utility stream read/write functions.
	template<typename T>
	static std::ostream& write(std::ostream& os, T const& val) {
		return os.write(reinterpret_cast<char const*>(&val), sizeof(T));
	}
	template<typename T>
	static std::istream& read(std::istream& is, T& val) {
		return is.read(reinterpret_cast<char*>(&val), sizeof(T));
	}

public:
	// Creates an empty `CellGenerator`.
	CellGenerator() : _tree(Data{ 1. }) {
		fill_primes(_tree.root());
	}

	TreeType const& tree() const {
		return _tree;
	}

	// Returns the prime value. This can be combined with the weights to
	// integrate the distribution.
	R prime() const {
		CellHandle root = _tree.root();
		if (_tree.type(root) == CellType::BRANCH) {
			return _tree.branch_data(root).prime;
		} else {
			return _tree.leaf_data(root).prime;
		}
	}

	// Samples.
	template<typename Rnd>
	void generate(Rnd& rnd, R* weight_out, Point<D, R>* point_out) const {
		// Generate random number.
		std::uniform_real_distribution<R> dist(0, 1);
		R choose_cell = dist(rnd);
		Point<D, R> choose_point;
		for (Dim dim = 0; dim < D; ++dim) {
			choose_point[dim] = dist(rnd);
		}
		// Draw from the root cell.
		Point<D, R> offset;
		Point<D, R> extent;
		offset.fill(0.);
		extent.fill(1.);
		select_cell(
			_tree.root(),
			offset, extent,
			choose_cell, choose_point,
			weight_out, point_out);
	}

	// Samples with re-weighting by distribution.
	template<typename Rnd, typename F>
	void generate(Rnd& rnd, F const& func, R* weight_out, Point<D, R>* point_out) const {
		static_assert(
			IsDistribution<D, R, F>::value,
			"F must be a valid distribution type (able to take a Point<D, R> "
			"to a type convertible to R with noexcept).");
		generate(rnd, weight_out, point_out);
		*weight_out *= func(*point_out);
	}

	// Writes to an output stream using a binary format.
	std::ostream& write(std::ostream& os) const {
		std::vector<CellHandle> cells;
		std::size_t size_size_t, size_tag, size_dim, size_split, size_prime;
		// Write identification string.
		if (!os.write(STREAM_HEAD, sizeof(STREAM_HEAD))) { goto error; }
		// Write sizes of important types, just as a quick validation.
		size_size_t = sizeof(std::size_t);
		size_tag = sizeof(unsigned char);
		size_dim = sizeof(Dim);
		size_split = sizeof(R);
		size_prime = sizeof(R);
		if (!write(os, size_size_t)) { goto error; };
		if (!write(os, size_tag)) { goto error; };
		if (!write(os, size_dim)) { goto error; };
		if (!write(os, size_split)) { goto error; };
		if (!write(os, size_prime)) { goto error; };
		// Write total number of cells.
		if (!write(os, _tree.size())) { goto error; }
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
				if (!write(os, tag)) { goto error; }
				if (!write(os, branch.dim)) { goto error; }
				if (!write(os, branch.split)) { goto error; }
				// Add branch children to be processed next, in order.
				for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
					cells.push_back(_tree.child(cell, child_idx));
				}
			} else {
				// Leafs store the following:
				// * Tag (unsigned char).
				// * Prime (integration real).
				unsigned char tag = 1;
				R prime = _tree.leaf_data(cell).prime;
				if (!write(os, tag)) { goto error; }
				if (!write(os, prime)) { goto error; }
			}
		}
		return os;
	error:
		os.setstate(std::ios_base::failbit);
		return os;
	}

	// Reads from an input stream using a binary format.
	std::istream& read(std::istream& is) {
		std::vector<CellHandle> cells;
		std::size_t size_size_t, size_tag, size_dim, size_split, size_prime;
		// Read identifier string.
		char head[sizeof(STREAM_HEAD)];
		if (!is.read(head, sizeof(STREAM_HEAD))) { goto error; }
		for(std::size_t idx = 0; idx < sizeof(STREAM_HEAD); ++idx) {
			if (head[idx] != STREAM_HEAD[idx]) {
				goto error;
			}
		}
		// Read important data sizes.
		if (!read(is, size_size_t)) { goto error; }
		if (!read(is, size_tag)) { goto error; }
		if (!read(is, size_dim)) { goto error; }
		if (!read(is, size_split)) { goto error; }
		if (!read(is, size_prime)) { goto error; }
		if (size_size_t != sizeof(std::size_t)
				|| size_tag != sizeof(unsigned char)
				|| size_dim != sizeof(Dim)
				|| size_split != sizeof(R)
				|| size_prime != sizeof(R)) {
			goto error;
		}
		// Read number of cells.
		std::size_t num_cells;
		if (!read(is, num_cells)) { goto error; }
		// Read the cells themselves.
		_tree.clear(Data{ 0. });
		_tree.reserve(num_cells);
		cells.push_back(_tree.root());
		while (!cells.empty()) {
			// Read the next cell from the stream.
			CellHandle cell = cells.back();
			cells.pop_back();
			unsigned char tag;
			if (!read(is, tag)) { goto error; }
			if (tag == 0) {
				// Branch.
				Dim dim;
				R split;
				if (!read(is, dim)) { goto error; }
				if (!read(is, split)) { goto error; }
				_tree.split(
					cell,
					typename TreeType::Branch(dim, split),
					Data(0.),
					{ Data(0.), Data(0.) });
				for (std::size_t child_idx = 0; child_idx < TreeType::NUM_CHILDREN; ++child_idx) {
					cells.push_back(_tree.child(cell, child_idx));
				}
			} else if (tag == 1) {
				// Leaf.
				R prime;
				if (!read(is, prime)) { goto error; }
				_tree.leaf_data(cell).prime = prime;
			} else {
				goto error;
			}
		}
		if (_tree.size() > num_cells) {
			goto error;
		}
		fill_primes(_tree.root());
		return is;
	error:
		// On error, set a default value.
		_tree.clear(Data{ 0. });
		is.setstate(std::ios_base::failbit);
		return is;
	}
};

// Convenience functions to allow for template parameter inference.
template<Dim D, typename R>
inline CellGenerator<D, R> make_generator() {
	return CellGenerator<D, R>();
}
template<Dim D, typename R, typename F>
inline CellGenerator<D, R> make_generator(
		CellBuilder<D, R, F> const& builder) {
	using Generator = CellGenerator<D, R>;
	using Builder = CellBuilder<D, R, F>;
	CellGenerator<D, R> result;
	// Creates a `CellGenerator` directly from a `CellBuilder`.
	result._tree = builder._tree.template transform<
			typename Generator::Data, typename Generator::Data>(
		[](typename Builder::BranchData) {
			return typename Generator::Data { 0. };
		},
		[](typename Builder::LeafData data) {
			return typename Generator::Data { data.stats.est_sqrt_m2() };
		});
	result.fill_primes(result._tree.root());
	return result;
}

}

#endif

