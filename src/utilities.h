#ifndef __UTILITIES_H__
#define __UTILITIES_H__

// [[Rcpp::plugins(cpp17)]]

struct FlatMatrix;
struct IntMatrix;
struct BoolMatrix;
struct FlatArray;
struct DataFrameCpp;
struct ListCpp;

#include <algorithm>     // find, sort, unique,
#include <cmath>         // fabs, isnan, sqrt,
#include <cstddef>       // size_t
#include <cstdint>       // uint8_t
#include <cstring>       // memcpy, memmove
#include <functional>    // function
#include <iomanip>       // fixed, setprecision
#include <iostream>      // cout, ostream
#include <iterator>      // distance
#include <limits>        // numeric_limits
#include <sstream>       // ostringstream
#include <stdexcept>     // out_of_range
#include <string>        // string
#include <type_traits>   // is_convertible
#include <utility>       // declval
#include <vector>        // vector

#include <Rcpp.h>


inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
inline constexpr double POS_INF = std::numeric_limits<double>::infinity();

// Constants for numerical stability
// Machine epsilon for double precision
inline constexpr double EPSILON = 2.2204460492503131e-16;
inline constexpr double MIN_PROB = EPSILON;           // ~2.22e-16
inline constexpr double MAX_PROB = 1.0 - EPSILON;     // ~1.0 - 2.22e-16

// Maximum safe quantile value (corresponds to p ≈ 1 - 1e-16)
// qnorm(1 - 2.22e-16) ≈ 8.1258906647
inline constexpr double MAX_NORMAL_QUANTILE = 8.125890664701906;
inline constexpr double MIN_NORMAL_QUANTILE = -8.125890664701906;

// Extreme z-score threshold for pnorm
// For |z| > EXTREME_Z, pnorm returns 0 or 1 within machine precision
inline constexpr double EXTREME_Z = 37.5;

// --------------------------- Distribution helpers --------------------------
double boost_pnorm(double q, double mean = 0.0, double sd = 1.0,
                   bool lower_tail = true);
double boost_qnorm(double p, double mean = 0.0, double sd = 1.0,
                   bool lower_tail = true);
double boost_dnorm(double x, double mean = 0.0, double sd = 1.0);

double boost_plogis(double q, double location = 0.0, double scale = 1.0,
                    bool lower_tail = true);
double boost_qlogis(double p, double location = 0.0, double scale = 1.0,
                    bool lower_tail = true);
double boost_dlogis(double x, double location = 0.0, double scale = 1.0);

double boost_pextreme(double q, double location = 0.0, double scale = 1.0,
                      bool lower_tail = true);
double boost_qextreme(double p, double location = 0.0, double scale = 1.0,
                      bool lower_tail = true);
double boost_dextreme(double x, double location = 0.0, double scale = 1.0);

double boost_pchisq(double q, double df, bool lower_tail = true);
double boost_qchisq(double p, double df, bool lower_tail = true);

double boost_pt(double q, double df, bool lower_tail = true);
double boost_qt(double p, double df, bool lower_tail = true);

double pnorm_fast(double x);
double qnorm_acklam(double p);


struct DoubleView {
  const double* p = nullptr;
  std::size_t n = 0;

  const double& operator[](std::size_t i) const { return p[i]; }
  std::size_t size() const { return n; }
};

// --------------------------- Small utilities --------------------------------
inline double sq(double x) noexcept { return x * x; }


// seqcpp: inclusive sequence; inputs are int
std::vector<std::size_t> seqcpp(std::size_t start, std::size_t end);

// convertLogicalVector: convert Rcpp LogicalVector to std::vector<unsigned char>
std::vector<unsigned char> convertLogicalVector(const Rcpp::LogicalVector& vec);

// which: return indices of true values
std::vector<std::size_t> which(const std::vector<unsigned char>& vec);

// expand1: expand vector to full length
std::vector<double> expand1(
    const std::vector<double>& v,
    const std::size_t nintervals,
    const char* name);

// expand_stratified: expand to a matrix with nintv rows and nstrata columns
FlatMatrix expand_stratified(
    const std::vector<double>& v,
    const std::size_t nstrata,
    const std::size_t nintv,
    const char* name);

void expand_stratified_to_slice(
    const std::vector<double>& v,
    FlatArray& out,
    std::size_t slice_index,
    std::size_t nstrata,
    std::size_t nintv,
    const char* name);

// findInterval: adapted helper (return indices following R-like convention)
std::size_t findInterval1(const double x,
                          const std::vector<double>& v,
                          bool rightmost_closed = false,
                          bool all_inside = false,
                          bool left_open = false);

std::vector<std::size_t> findInterval3(const std::vector<double>& x,
                                       const std::vector<double>& v,
                                       bool rightmost_closed = false,
                                       bool all_inside = false,
                                       bool left_open = false);

// all_equal: check if all elements in v equal target within tolerance tol
inline bool all_equal(const std::vector<double>& v, double target, double tol = 0.0) {
  if (v.empty()) return true;              // mimic R: all(logical(0)) == TRUE
  if (tol == 0.0) {
    for (double x : v) if (!(x == target)) return false;
  } else {
    for (double x : v) if (std::fabs(x - target) > tol) return false;
  }
  return true;
}

// mean using Kahan summation for improved numerical stability
inline double mean_kahan(const std::vector<double>& v) {
  const std::size_t n = v.size();
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  double sum = 0.0;
  double c = 0.0; // compensation
  for (std::size_t i = 0; i < n; ++i) {
    double y = v[i] - c;        // corrected addend
    double t = sum + y;         // provisional sum
    c = (t - sum) - y;          // new compensation
    sum = t;
  }
  return sum / static_cast<double>(n);
}

// mean and sd using Welford's method
inline void mean_sd(const double* data, std::size_t n, double &omean, double &osd) {
  if (n == 0) {
    omean = std::numeric_limits<double>::quiet_NaN();
    osd = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  double mean = 0.0;
  double M2 = 0.0;     // sum of squares of differences
  double count = 0.0;

  for (std::size_t i = 0; i < n; ++i) {
    ++count;
    double x = data[i];
    double delta = x - mean;
    mean += delta / count;
    double delta2 = x - mean;
    M2 += delta * delta2;
  }

  omean = mean;
  osd = (count > 1) ? std::sqrt(M2 / (count - 1)) : 0.0;
}

inline void mean_sd(const double* data, int n, double &omean, double &osd) {
  mean_sd(data, static_cast<std::size_t>(n), omean, osd);
}

double extract_sum(const DataFrameCpp& df, const char* name);


// --------------------------- Root finders -----------------------------------
template <class F>
double brent(F&& f, double x1, double x2, double tol = 1e-8, int maxiter = 100) {
  constexpr double EPS = 3.0e-8;

  double a = x1, b = x2, c = x2;
  double fa = f(a), fb = f(b), fc = fb;
  double d = 0.0, d1 = 0.0;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    throw std::invalid_argument("Root must be bracketed in brent");

  for (int iter = 1; iter <= maxiter; ++iter) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a; fc = fa; d = b - a; d1 = d;
    }

    if (std::fabs(fc) < std::fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }

    double tol1 = 2.0 * EPS * std::fabs(b) + 0.5 * tol;
    double xm = 0.5 * (c - b);
    if (std::fabs(xm) <= tol1 || fb == 0.0) return b;

    double p, q, r, s;
    if (std::fabs(d1) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
      s = fb / fa;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) q = -q;
      p = std::fabs(p);
      double min1 = 3.0 * xm * q - std::fabs(tol1 * q);
      double min2 = std::fabs(d1 * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2)) {
        d1 = d; d = p / q;
      } else {
        d = xm; d1 = d;
      }
    } else {
      d = xm; d1 = d;
    }

    a = b; fa = fb;
    if (std::fabs(d) > tol1) b += d; else b += std::copysign(tol1, xm);
    fb = f(b);
  }
  throw std::runtime_error("Maximum iterations exceeded in brent");
}


// check if no elements are missing
inline bool none_na(const std::vector<double>& v) {
  return !v.empty() &&
    std::none_of(v.begin(), v.end(), [](double x){ return std::isnan(x); });
}

inline bool none_na(const std::vector<int>& v) {
  return !v.empty() &&
    std::none_of(v.begin(), v.end(), [](double x){ return x == INT_MIN; });
}

inline bool none_na(const std::vector<unsigned char>& v) {
  return !v.empty() &&
    std::none_of(v.begin(), v.end(), [](double x){ return x == 255; });
}

// check if any element is non-increasing compared to previous
template<typename T>
bool any_nonincreasing(const std::vector<T>& I) {
  if (I.size() < 2) return false;
  for (std::size_t i = 1; i < I.size(); ++i) {
    if (I[i] <= I[i-1]) return true;
  }
  return false;
}


// subset: return a subset of v according to 'order' (indices)
template <typename T>
std::vector<T> subset(const std::vector<T>& v, const std::vector<std::size_t>& order) {
  std::vector<T> result(order.size());
  std::size_t n = order.size();
  std::size_t nv = v.size();
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t index = order[i];
    if (index < 0 || index >= nv) {
      throw std::out_of_range(
          "Index in 'order' is out of bounds for the source vector.");
    }
    result[i] = v[index];
  }
  return result;
}

// subset_in_place: reorder/keep elements of v according to 'order' (indices)
template <typename T>
void subset_in_place(std::vector<T>& v, const std::vector<std::size_t>& order) {
  std::vector<T> temp_subset(order.size());
  std::size_t n = order.size();
  std::size_t nv = v.size();
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t index = order[i];
    if (index < 0 || index >= nv) {
      throw std::out_of_range(
          "Index in 'order' is out of bounds for the source vector.");
    }
    temp_subset[i] = v[index];
  }
  v = std::move(temp_subset);
}

// Return a new vector containing elements v[start, end).
// Preconditions required by you: 0 <= start < end (and end <= v.size()).
template <typename T>
std::vector<T> subset(const std::vector<T>& v, std::size_t start, std::size_t end) {
  if (start < 0) throw std::out_of_range("subset: start < 0");
  if (end < 0) throw std::out_of_range("subset: end < 0");
  const std::size_t vsz = v.size();
  if (static_cast<std::size_t>(end) > vsz)
    throw std::out_of_range("subset: end > v.size()");
  if (!(start < end)) throw std::invalid_argument("subset: require start < end");

  const std::size_t s = static_cast<std::size_t>(start);
  const std::size_t e = static_cast<std::size_t>(end);
  const std::size_t n = e - s;

  if constexpr (std::is_trivially_copyable_v<T>) {
    std::vector<T> out;
    out.resize(n); // allocate contiguous buffer
    if (n > 0) {
      std::memcpy(static_cast<void*>(out.data()),
                  static_cast<const void*>(v.data() + s),
                  n * sizeof(T));
    }
    return out;
  } else {
    // non-trivial types: element-wise copy constructor
    return std::vector<T>(v.begin() + static_cast<std::ptrdiff_t>(s),
                          v.begin() + static_cast<std::ptrdiff_t>(e));
  }
}

// In-place subset: keep elements [start, end) and discard the rest.
// Preconditions required by you: 0 <= start < end (and end <= v.size()).
template <typename T>
void subset_in_place(std::vector<T>& v, std::size_t start, std::size_t end) {
  if (start < 0) throw std::out_of_range("subset_in_place: start < 0");
  if (end < 0) throw std::out_of_range("subset_in_place: end < 0");
  const std::size_t vsz = v.size();
  if (static_cast<std::size_t>(end) > vsz)
    throw std::out_of_range("subset_in_place: end > v.size()");
  if (!(start < end))
    throw std::invalid_argument("subset_in_place: require start < end");

  const std::size_t s = static_cast<std::size_t>(start);
  const std::size_t e = static_cast<std::size_t>(end);
  const std::size_t n = e - s; // number of elements to keep

  if (s == 0) {
    // already at beginning; just resize down to requested length
    v.resize(n);
    return;
  }

  if constexpr (std::is_trivially_copyable_v<T>) {
    // overlapping move: use memmove (safe for overlapping ranges)
    std::memmove(static_cast<void*>(v.data()),
                 static_cast<const void*>(v.data() + s),
                 n * sizeof(T));
    v.resize(n);
  } else {
    // non-trivial types: use std::move for element-wise move construction
    std::move(v.begin() + static_cast<std::ptrdiff_t>(s),
              v.begin() + static_cast<std::ptrdiff_t>(e),
              v.begin());
    v.resize(n);
  }
}

// unique_sorted: return sorted unique values
template <typename T>
std::vector<T> unique_sorted(const std::vector<T>& v) {
  std::vector<T> w = v;
  std::sort(w.begin(), w.end());
  w.erase(std::unique(w.begin(), w.end()), w.end());
  return w;
}

// matchcpp: for each element of x find its index in table or -1 if not found
template <typename T>
std::vector<int> matchcpp(const std::vector<T>& x, const std::vector<T>& table,
                          const std::size_t start_index = 0) {
  std::vector<int> result(x.size());
  std::size_t n = x.size();
  for (std::size_t i = 0; i < n; ++i) {
    auto it = std::find(table.begin(), table.end(), x[i]);
    if (it != table.end()) {
      result[i] = static_cast<int>(std::distance(table.begin(), it)) + start_index;
    } else {
      result[i] = -1;
    }
  }
  return result;
}


double dtpwexpcpp1(
    const double q,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound = 0.0,
    const bool logd = false);

double ptpwexpcpp1(
    const double q,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound = 0.0,
    const bool lowertail = true,
    const bool logp = false);


template <class VTIME, class VLam>
inline double qtpwexpcpp1(
    const double p,
    const VTIME& piecewiseSurvivalTime,
    const VLam& lambda,
    const double lowerBound = 0.0,
    const bool lowertail = true,
    const bool logp = false) {

  std::size_t m = piecewiseSurvivalTime.size();
  double u = logp ? std::exp(p) : p;
  if (!lowertail) u = 1.0 - u;
  if (u <= 0.0) return lowerBound;
  if (u >= 1.0) return std::numeric_limits<double>::infinity();
  double v1 = -log1p(-u);
  std::size_t j = 0;
  while (j < m && piecewiseSurvivalTime[j] <= lowerBound) ++j;
  std::size_t j1 = (j == 0) ? 0 : (j - 1);
  double v = 0.0;
  if (j1 == m - 1) {
    double lj = lambda[j1];
    if (lj <= 0.0) return std::numeric_limits<double>::infinity();
    return lowerBound + v1 / lj;
  }
  for (j = j1; j < m - 1; ++j) {
    double dt = (j == j1) ? piecewiseSurvivalTime[j + 1] - lowerBound :
    piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j];
    double lj = lambda[j];
    if (lj > 0.0) v += lj * dt;
    if (v >= v1) break;
  }
  double lj = lambda[j];
  if (lj <= 0.0) return std::numeric_limits<double>::infinity();
  if (j == m - 1) {
    double dt = (v1 - v) / lj;
    return piecewiseSurvivalTime[j] + dt;
  }
  double dt = (v - v1) / lj;
  return piecewiseSurvivalTime[j + 1] - dt;
}


double intnorm(const std::function<double(double)>& f,
               double mu, double sigma, double a, double b);

std::pair<double, double> mini(const std::function<double(double)>& f,
                               double x1, double x2);

double integrate3(const std::function<double(double)>& f,
                  const std::vector<double>& breaks, double tol = 1e-8);

double quad2d(const std::function<double(double,double)>& f,
              double ax, double bx, double ay, double by,
              double tol = 1.0e-5);

double pbvnormcpp(const std::vector<double>& lower,
                  const std::vector<double>& upper,
                  const double rho);

ListCpp hazard_pdcpp(const std::vector<double>& piecewiseSurvivalTime,
                     const std::vector<double>& hazard_pfs,
                     const std::vector<double>& hazard_os,
                     const double rho_pd_os);

double corr_pfs_oscpp(const std::vector<double>& piecewiseSurvivalTime,
                      const std::vector<double>& hazard_pfs,
                      const std::vector<double>& hazard_os,
                      const double rho_pd_os);

ListCpp hazard_subcpp(const std::vector<double>& piecewiseSurvivalTime,
                      const std::vector<double>& hazard_itt,
                      const std::vector<double>& hazard_pos,
                      const double p_pos);


ListCpp bygroup(const DataFrameCpp& data,
                const std::vector<std::string>& variables);

int cholesky2(FlatMatrix& matrix, std::size_t n, double toler = 1e-12);
void chsolve2(FlatMatrix& matrix, std::size_t n, double* y);
FlatMatrix invsympd(const FlatMatrix& matrix, std::size_t n, double toler = 1e-12);
FlatMatrix invchol(FlatMatrix& matrix, std::size_t n);

std::vector<double> mat_vec_mult(const FlatMatrix& A, const std::vector<double>& x);
FlatMatrix mat_mat_mult(const FlatMatrix& A, const FlatMatrix& B);

FlatMatrix transpose(const FlatMatrix& M);
IntMatrix transpose(const IntMatrix& M);
BoolMatrix transpose(const BoolMatrix& M);

double quadsym(const std::vector<double>& u, const FlatMatrix& v);

// Print a std::vector<T> to std::cout.
// Requirements: T must be streamable via operator<< to std::ostream.
//
// Parameters:
//  - v: vector to print
//  - label: optional prefix printed before the vector
//  - precision: if >= 0, sets std::fixed and std::setprecision(precision)
//    for floating values
//  - head: number of leading elements to show when truncated
//  - tail: number of trailing elements to show when truncated
//  - sep: separator between elements (default ", ")
//  - show_indices: if true prints each element as "idx: value"
//  - endline: whether to append a newline at the end (true by default)
template <typename T>
void print_vector(const std::vector<T>& v,
                  const std::string& label = "",
                  int precision = -1,
                  std::size_t head = 5,
                  std::size_t tail = 5,
                  const std::string& sep = ", ",
                  bool show_indices = false,
                  bool endline = true) {
  static_assert(
    std::is_convertible<decltype(
      std::declval<std::ostream&>() << std::declval<T>()), std::ostream&>::value,
                   "Type T must be streamable to std::ostream (operator<<)");
  std::ostringstream ss;
  if (!label.empty()) ss << label << ": ";

  std::size_t n = v.size();
  if (n == 0) {
    ss << "[]";
    if (endline) ss << '\n';
    std::cout << ss.str();
    return;
  }

  // Configure precision only if requested
  bool use_precision = (precision >= 0);
  if (use_precision) ss << std::fixed << std::setprecision(precision);

  ss << "[";
  auto print_elem = [&](std::size_t i) {
    if (show_indices) ss << i << ": ";
    if constexpr (std::is_same_v<T, unsigned char> ||
                  std::is_same_v<T, std::uint8_t>) {
      // print unsigned char / uint8_t as integer 0/1 (not as a character)
      ss << static_cast<int>(v[i]);
    } else {
      ss << v[i];
    }
  };

  if (n <= head + tail || head + tail == 0) {
    for (std::size_t i = 0; i < n; ++i) {
      if (i) ss << sep;
      print_elem(i);
    }
  } else {
    // print head
    for (std::size_t i = 0; i < head; ++i) {
      if (i) ss << sep;
      print_elem(i);
    }
    ss << sep << "..." << sep;
    // print tail
    for (std::size_t j = n - tail; j < n; ++j) {
      if (j != n - tail) ss << sep;
      print_elem(j);
    }
  }
  ss << "]";
  if (endline) ss << '\n';

  std::cout << ss.str();
}

#endif // __UTILITIES__
