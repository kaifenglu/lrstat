#ifndef __UTILITIES_H__
#define __UTILITIES_H__

// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>

struct FlatMatrix;
struct IntMatrix;
struct BoolMatrix;
struct DataFrameCpp;
struct ListCpp;

#include <algorithm>   // copy, find, sort, unique,
#include <cmath>       // sqrt, isnan
#include <cstddef>     // size_t
#include <cstdint>     // uint64_t
#include <cstring>     // memcpy, memmove
#include <functional>  // function
#include <iomanip>     // fixed, setprecision
#include <iostream>    // cout, ostream
#include <iterator>    // distance
#include <limits>      // numeric_limits
#include <numeric>     // accumulate
#include <sstream>     // ostringstream
#include <stdexcept>   // out_of_range
#include <string>      // string
#include <type_traits> // is_convertible
#include <utility>     // declval
#include <vector>      // vector

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

// --------------------------- Small utilities --------------------------------
inline double sq(double x) noexcept { return x * x; }

std::vector<double> stl_sort(const std::vector<double>& x);

// seqcpp: inclusive sequence; inputs are int
std::vector<int> seqcpp(int start, int end);

// convertLogicalVector: convert Rcpp LogicalVector to std::vector<unsigned char>
std::vector<unsigned char> convertLogicalVector(const Rcpp::LogicalVector& vec);

// which: return indices of true values
std::vector<int> which(const std::vector<unsigned char>& vec);

// findInterval3: adapted helper (return indices following R-like convention)
std::vector<int> findInterval3(const std::vector<double>& x,
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

// --------------------------- Root finders -----------------------------------
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol = 1e-8, int maxiter = 100);
double bisect(const std::function<double(double)>& f,
              double x1, double x2, double tol = 1e-8, int maxiter = 100);

// check if no elements are missing
inline bool none_na(const std::vector<double>& v) {
  return std::none_of(v.begin(), v.end(), [](double x){ return std::isnan(x); });
}

inline bool none_na(const std::vector<int>& v) {
  return std::none_of(v.begin(), v.end(), [](double x){ return x == INT_MIN; });
}

inline bool none_na(const std::vector<unsigned char>& v) {
  return std::none_of(v.begin(), v.end(), [](double x){ return x == 255; });
}

// check if any element is non-increasing compared to previous
template<typename T>
bool any_nonincreasing(const std::vector<T>& I) {
  if (I.size() < 2) return false;
  for (size_t i = 1; i < I.size(); ++i) {
    if (I[i] <= I[i-1]) return true;
  }
  return false;
}


// subset: return a subset of v according to 'order' (indices)
template <typename T>
std::vector<T> subset(const std::vector<T>& v, const std::vector<int>& order) {
  std::vector<T> result(order.size());
  int n = order.size();
  int nv = v.size();
  for (int i = 0; i < n; ++i) {
    int index = order[i];
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
void subset_in_place(std::vector<T>& v, const std::vector<int>& order) {
  std::vector<T> temp_subset(order.size());
  int n = order.size();
  int nv = v.size();
  for (int i = 0; i < n; ++i) {
    int index = order[i];
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
std::vector<T> subset(const std::vector<T>& v, int start, int end) {
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
void subset_in_place(std::vector<T>& v, int start, int end) {
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
                          const int start_index = 0) {
  std::vector<int> result(x.size());
  int n = x.size();
  for (int i = 0; i < n; ++i) {
    auto it = std::find(table.begin(), table.end(), x[i]);
    if (it != table.end()) {
      result[i] = static_cast<int>(std::distance(table.begin(), it)) + start_index;
    } else {
      result[i] = -1;
    }
  }
  return result;
}


double errorSpentcpp(const double t,
                     const double error = 0.025,
                     const std::string& sf = "sfOF",
                     const double sfpar = 0.0);

ListCpp exitprobcpp(const std::vector<double>& b,
                    const std::vector<double>& a,
                    const std::vector<double>& theta,
                    const std::vector<double>& I);

double dtpwexpcpp1(const double q,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound = 0.0,
                   const bool logd = false);

std::vector<double> dtpwexpcpp(
    const std::vector<double>& q,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound = 0.0,
    const bool logd = false);

double ptpwexpcpp1(const double q,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound = 0.0,
                   const bool lowertail = true,
                   const bool logp = false);

std::vector<double> ptpwexpcpp(
    const std::vector<double>& q,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound = 0.0,
    const bool lowertail = true,
    const bool logp = false);

double qtpwexpcpp1(const double p,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound = 0.0,
                   const bool lowertail = true,
                   const bool logp = false);

std::vector<double> qtpwexpcpp(
    const std::vector<double>& p,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound = 0.0,
    const bool lowertail = true,
    const bool logp = false);

ListCpp mtpwexpcpp(const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound = 0.0);

std::vector<double> getBoundcpp(
    const int k,
    const std::vector<double>& informationRates,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const std::vector<unsigned char>& efficacyStopping);

ListCpp getPower(
    const double alpha,
    const int kMax,
    const std::vector<double>& b,
    const std::vector<double>& theta,
    const std::vector<double>& I,
    const std::string& bsf,
    const double bsfpar,
    const std::vector<double>& st,
    const std::vector<unsigned char>& futilityStopping,
    const std::vector<double>& w);

double intnorm(const std::function<double(double)>& f,
               double mu, double sigma, double a, double b);

std::pair<double, double> mini(
    const std::function<double(double)>& f, double x1, double x2);

double quad(const std::function<double(double)>& f,
            double lower, double upper, double tol = 1e-8, unsigned maxiter = 1000);

double quad2d(const std::function<double(double,double)>& f,
              double ax, double bx, double ay, double by, double tol = 1.0e-5);

double pbvnormcpp(std::vector<double> lower, std::vector<double> upper, double rho);

ListCpp hazard_pdcpp(const std::vector<double>& piecewiseSurvivalTime,
                     const std::vector<double>& hazard_pfs,
                     const std::vector<double>& hazard_os,
                     const double rho_pd_os);

ListCpp hazard_subcpp(const std::vector<double>& piecewiseSurvivalTime,
                      const std::vector<double>& hazard_itt,
                      const std::vector<double>& hazard_pos,
                      const double p_pos);

std::vector<double> accrual(const std::vector<double>& time,
                            const std::vector<double>& accrualTime,
                            const std::vector<double>& accrualIntensity,
                            const double accrualDuration);

std::vector<double> getAccrualDurationFromN(
    const std::vector<double>& nsubjects,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity);

std::vector<double> patrisk(const std::vector<double>& time,
                            const std::vector<double>& piecewiseSurvivalTime,
                            const std::vector<double>& lambda,
                            const std::vector<double>& gamma);

std::vector<double> pevent(const std::vector<double>& time,
                           const std::vector<double>& piecewiseSurvivalTime,
                           const std::vector<double>& lambda,
                           const std::vector<double>& gamma);

FlatMatrix natrisk(const std::vector<double>& time,
                   const double allocationRatioPlanned,
                   const std::vector<double>& accrualTime,
                   const std::vector<double>& accrualIntensity,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda1,
                   const std::vector<double>& lambda2,
                   const std::vector<double>& gamma1,
                   const std::vector<double>& gamma2,
                   const double accrualDuration,
                   const double minFollowupTime,
                   const double maxFollowupTime);

FlatMatrix nevent(const std::vector<double>& time,
                  const double allocationRatioPlanned,
                  const std::vector<double>& accrualTime,
                  const std::vector<double>& accrualIntensity,
                  const std::vector<double>& piecewiseSurvivalTime,
                  const std::vector<double>& lambda1,
                  const std::vector<double>& lambda2,
                  const std::vector<double>& gamma1,
                  const std::vector<double>& gamma2,
                  const double accrualDuration,
                  const double minFollowupTime,
                  const double maxFollowupTime);

FlatMatrix nevent2(const std::vector<double>& time,
                   const double allocationRatioPlanned,
                   const std::vector<double>& accrualTime,
                   const std::vector<double>& accrualIntensity,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda1,
                   const std::vector<double>& lambda2,
                   const std::vector<double>& gamma1,
                   const std::vector<double>& gamma2,
                   const double accrualDuration,
                   const double minFollowupTime,
                   const double maxFollowupTime);

ListCpp getDesigncpp(const double beta,
                     const double IMax,
                     const double theta,
                     const int kMax,
                     const std::vector<double>& informationRates,
                     const std::vector<unsigned char>& efficacyStopping,
                     const std::vector<unsigned char>& futilityStopping,
                     const std::vector<double>& criticalValues,
                     const double alpha,
                     const std::string& typeAlphaSpending,
                     const double parameterAlphaSpending,
                     const std::vector<double>& userAlphaSpending,
                     const std::vector<double>& futilityBounds,
                     const std::string& typeBetaSpending,
                     const double parameterBetaSpending,
                     const std::vector<double>& userBetaSpending,
                     const std::vector<double>& spendingTime,
                     const double varianceRatio = 1.0);

ListCpp getDesignEquivcpp(const double beta,
                          const double IMax,
                          const double thetaLower,
                          const double thetaUpper,
                          const double theta,
                          const int kMax,
                          const std::vector<double>& informationRates,
                          const std::vector<double>& criticalValues,
                          const double alpha,
                          const std::string& typeAlphaSpending,
                          const double parameterAlphaSpending,
                          const std::vector<double>& userAlphaSpending,
                          const std::vector<double>& spendingTime);

ListCpp adaptDesigncpp(double betaNew,
                       double INew,
                       const int L,
                       const double zL,
                       const double theta,
                       const double IMax,
                       const int kMax,
                       const std::vector<double>& informationRates,
                       const std::vector<unsigned char>& efficacyStopping,
                       const std::vector<unsigned char>& futilityStopping,
                       const std::vector<double>& criticalValues,
                       const double alpha,
                       const std::string& typeAlphaSpending,
                       const double parameterAlphaSpending,
                       const std::vector<double>& userAlphaSpending,
                       const std::vector<double>& futilityBounds,
                       const std::string& typeBetaSpending,
                       const double parameterBetaSpending,
                       const std::vector<double>& spendingTime,
                       const bool MullerSchafer,
                       const int kNew,
                       const std::vector<double>& informationRatesNew,
                       const std::vector<unsigned char>& efficacyStoppingNew,
                       const std::vector<unsigned char>& futilityStoppingNew,
                       const std::string& typeAlphaSpendingNew,
                       const double parameterAlphaSpendingNew,
                       const std::string& typeBetaSpendingNew,
                       const double parameterBetaSpendingNew,
                       const std::vector<double>& userBetaSpendingNew,
                       const std::vector<double>& spendingTimeNew,
                       const double varianceRatio);

ListCpp bygroup(const DataFrameCpp& data, const std::vector<std::string>& variables);

int cholesky2(FlatMatrix& matrix, int n, double toler = 1e-12);
void chsolve2(FlatMatrix& matrix, int n, double* y);
FlatMatrix invsympd(const FlatMatrix& matrix, int n, double toler = 1e-12);

FlatMatrix transpose(const FlatMatrix& M);
IntMatrix transpose(const IntMatrix& M);
BoolMatrix transpose(const BoolMatrix& M);

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
