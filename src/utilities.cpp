#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <algorithm>  // lower_bound, sort, upper_bound
#include <cmath>      // copysign, exp, fabs, isinf, isnan, log, sqrt
#include <functional> // function
#include <limits>     // numeric_limits
#include <numeric>    // inner_product, iota
#include <queue>      // priority_queue
#include <stdexcept>  // invalid_argument, runtime_error
#include <string>     // string
#include <utility>    // pair
#include <vector>     // vector

#include <Rcpp.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/tools/minima.hpp>

using std::size_t;


double boost_pnorm(double q, double mean, double sd, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (sd <= 0) throw std::invalid_argument("Standard deviation must be positive.");

  double z = (q - mean) / sd;
  if (lower_tail) {
    if (z <= -EXTREME_Z) return 0.0;
    if (z >= EXTREME_Z) return 1.0;
  } else {
    if (z >= EXTREME_Z) return 0.0;
    if (z <= -EXTREME_Z) return 1.0;
  }

  boost::math::normal_distribution<> dist(mean, sd);
  if (lower_tail) return boost::math::cdf(dist, q);
  else return boost::math::cdf(boost::math::complement(dist, q));
}

double boost_qnorm(double p, double mean, double sd, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (sd <= 0) throw std::invalid_argument("Standard deviation must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");

  // Clamp extreme probabilities to avoid overflow in boost::math::quantile
  bool at_extreme = false;
  double extreme_quantile = 0.0;

  if (lower_tail) {
    if (p <= MIN_PROB) {
      at_extreme = true;
      extreme_quantile = MIN_NORMAL_QUANTILE;
    } else if (p >= MAX_PROB) {
      at_extreme = true;
      extreme_quantile = MAX_NORMAL_QUANTILE;
    }
  } else {
    // When lower_tail = false, we compute quantile(1 - p)
    if (p <= MIN_PROB) {
      at_extreme = true;
      extreme_quantile = MAX_NORMAL_QUANTILE;
    } else if (p >= MAX_PROB) {
      at_extreme = true;
      extreme_quantile = MIN_NORMAL_QUANTILE;
    }
  }

  // If at extreme, return scaled value directly
  if (at_extreme) {
    return mean + sd * extreme_quantile;
  }

  // Safe to call boost quantile
  boost::math::normal_distribution<> dist(mean, sd);
  return lower_tail ? boost::math::quantile(dist, p) :
    boost::math::quantile(dist, 1.0 - p);
}

double boost_dnorm(double x, double mean, double sd) {
  if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
  if (sd <= 0) throw std::invalid_argument("Standard deviation must be positive.");
  boost::math::normal_distribution<> dist(mean, sd);
  return boost::math::pdf(dist, x);
}

double boost_plogis(double q, double location, double scale, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::logistic_distribution<> dist(location, scale);
  if (lower_tail) return boost::math::cdf(dist, q);
  else return boost::math::cdf(boost::math::complement(dist, q));
}

double boost_qlogis(double p, double location, double scale, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");
  boost::math::logistic_distribution<> dist(location, scale);
  return lower_tail ? boost::math::quantile(dist, p) :
    boost::math::quantile(dist, 1.0 - p);
}

double boost_dlogis(double x, double location, double scale) {
  if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::logistic_distribution<> dist(location, scale);
  return boost::math::pdf(dist, x);
}

double boost_pextreme(double q, double location, double scale, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::extreme_value_distribution<> dist(location, scale);
  // keep semantics consistent with complementary log-log link
  if (lower_tail) return boost::math::cdf(complement(dist, 2.0 * location - q));
  else return boost::math::cdf(dist, 2.0 * location - q);
}

double boost_qextreme(double p, double location, double scale, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");
  boost::math::extreme_value_distribution<> dist(location, scale);
  return lower_tail? -boost::math::quantile(complement(dist, p)) :
    -boost::math::quantile(complement(dist, 1.0 - p));
}

double boost_dextreme(double x, double location, double scale) {
  if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::extreme_value_distribution<> dist(location, scale);
  return boost::math::pdf(dist, -x);
}

double boost_pchisq(double q, double df, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  if (std::isinf(q)) {
    if (q > 0.0) return lower_tail ? 1.0 : 0.0;
    else return lower_tail ? 0.0 : 1.0;
  }
  boost::math::chi_squared_distribution<> dist(df);
  if (lower_tail) return boost::math::cdf(dist, q);
  else return boost::math::cdf(boost::math::complement(dist, q));
}

double boost_qchisq(double p, double df, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");
  boost::math::chi_squared_distribution<> dist(df);
  return lower_tail ? boost::math::quantile(dist, p) :
    boost::math::quantile(dist, 1.0 - p);
}

double boost_pt(double q, double df, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  boost::math::students_t_distribution<> dist(df);
  if (lower_tail) return boost::math::cdf(dist, q);
  else return boost::math::cdf(boost::math::complement(dist, q));
}

double boost_qt(double p, double df, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");
  boost::math::students_t_distribution<> dist(df);
  return lower_tail ? boost::math::quantile(dist, p) :
    boost::math::quantile(dist, 1.0 - p);
}


// Fast normal CDF approximation
// see https://doi.org/10.2139/ssrn.2842681
inline constexpr double g2 = -0.0150234471495426236132;
inline constexpr double g4 = 0.000666098511701018747289;
inline constexpr double g6 = 5.07937324518981103694e-06;
inline constexpr double g8 = -2.92345273673194627762e-06;
inline constexpr double g10 = 1.34797733516989204361e-07;
inline constexpr double m2dpi = -0.6366197723675813824329; // -2.0 / pi

// [[Rcpp::export]]
double pnorm_fast(double x) {
  if (!std::isfinite(x)) return (x > 0 ? 1.0 : 0.0);

  const double x2  = x * x;
  const double x4  = x2 * x2;
  const double x6  = x4 * x2;
  const double x8  = x6 * x2;
  const double x10 = x8 * x2;

  double tmp = 1.0 + g2 * x2 + g4 * x4 + g6 * x6 + g8 * x8 + g10 * x10;
  double t = 1.0 - std::exp(tmp * m2dpi * x2);
  return 0.5 + 0.5 * ((x > 0) - (x < 0)) * std::sqrt(t);
}


// Approximate inverse standard normal CDF (qnorm).
// Based on Peter John Acklam's rational approximation.
// Good accuracy for double precision when p in (0,1).
// [[Rcpp::export]]
double qnorm_acklam(double p) {
  // Coefficients in rational approximations
  static constexpr double a[] = {
    -3.969683028665376e+01,
    2.209460984245205e+02,
    -2.759285104469687e+02,
    1.383577518672690e+02,
    -3.066479806614716e+01,
    2.506628277459239e+00
  };
  static constexpr double b[] = {
    -5.447609879822406e+01,
    1.615858368580409e+02,
    -1.556989798598866e+02,
    6.680131188771972e+01,
    -1.328068155288572e+01
  };
  static constexpr double c[] = {
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
    4.374664141464968e+00,
    2.938163982698783e+00
  };
  static constexpr double d[] = {
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e+00,
    3.754408661907416e+00
  };

  // Protect against invalid inputs; caller should clamp too.
  p = std::min(std::max(p, 1e-300), 1.0 - 1e-16);

  // Define break-points.
  constexpr double plow  = 0.02425;
  constexpr double phigh = 1.0 - plow;

  double q, r;
  if (p < plow) {
    // Rational approximation for lower region.
    q = std::sqrt(-2.0 * std::log(p));
    return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
      ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
  } else if (p > phigh) {
    // Rational approximation for upper region.
    q = std::sqrt(-2.0 * std::log(1.0 - p));
    return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
      ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
  } else {
    // Rational approximation for central region.
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
      (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
  }
}


std::vector<size_t> seqcpp(size_t start, size_t end) {
  if (start > end) throw std::invalid_argument(
      "start must be less than or equal to end for the sequence function.");
  size_t size = end - start + 1;
  std::vector<size_t> result(size);
  std::iota(result.begin(), result.end(), start);
  return result;
}

std::vector<unsigned char> convertLogicalVector(const Rcpp::LogicalVector& vec) {
  size_t n = vec.size();
  std::vector<unsigned char> result;
  result.resize(n);
  for (size_t i = 0; i < n; ++i) {
    int v = vec[i];
    if (v == NA_LOGICAL) result[i] = 255; // NA representation
    else result[i] = v ? 1 : 0; // TRUE -> 1, FALSE -> 0
  }
  return result;
}

std::vector<size_t> which(const std::vector<unsigned char>& vec) {
  std::vector<size_t> indices;
  indices.reserve(vec.size());
  size_t n = vec.size();
  for (size_t i = 0; i < n; ++i) {
    if (vec[i] != 0 && vec[i] != 255) indices.push_back(i);
  }
  return indices;
}

std::vector<double> expand1(
    const std::vector<double>& v,
    const size_t nintervals,
    const char* name) {
  if (v.size() == 1) {
    return std::vector<double>(nintervals, v[0]);
  } else if (v.size() == nintervals) {
    return v;
  } else {
    throw std::invalid_argument(std::string("Invalid length for ") + name);
  }
}

FlatMatrix expand_stratified(
    const std::vector<double>& v,
    const size_t nstrata,
    const size_t nintv,
    const char* name) {

  FlatMatrix out(nintv, nstrata);
  if (v.size() == 1) {
    out.fill(v[0]);
  } else if (v.size() == nintv) {
    for (size_t s = 0; s < nstrata; ++s) {
      flatmatrix_set_column(out, s, v);
    }
  } else if (v.size() == nstrata * nintv) {
    out.data = v;
  } else {
    throw std::invalid_argument(std::string("Invalid length for ") + name);
  }
  return out;
}

void expand_stratified_to_slice(
    const std::vector<double>& v,
    FlatArray& out,
    size_t slice_index,
    size_t nstrata,
    size_t nintv,
    const char* name) {

  double* dst = out.slice_ptr(slice_index);
  if (v.size() == 1) {
    for (size_t i = 0; i < nintv * nstrata; ++i) {
      dst[i] = v[0];
    }
  } else if (v.size() == nintv) {
    const double* src = v.data();
    for (size_t s = 0; s < nstrata; ++s) {
      std::memcpy(dst + s * nintv, src, nintv * sizeof(double));
    }
  } else if (v.size() == nstrata * nintv) {
    const double* src = v.data();
    std::memcpy(dst, src, nintv * nstrata * sizeof(double));
  } else {
    throw std::invalid_argument(std::string("Invalid length for ") + name);
  }
}


size_t findInterval1(const double x,
                     const std::vector<double>& v,
                     bool rightmost_closed,
                     bool all_inside,
                     bool left_open) {

  const double* v_begin = v.data();
  const double* v_end   = v_begin + v.size();
  const size_t nv = v.size();

  const double* pos = left_open ? std::lower_bound(v_begin, v_end, x) :
    std::upper_bound(v_begin, v_end, x);
  size_t idx = static_cast<size_t>(pos - v_begin);
  if (rightmost_closed) {
    if (left_open) {
      if (x == v[0]) idx = 1;
    } else {
      if (x == v[nv - 1]) idx = nv - 1;
    }
  }
  if (all_inside) {
    if (idx == 0) idx = 1;
    else if (idx == nv) idx = nv - 1;
  }

  return idx;
}


std::vector<size_t> findInterval3(const std::vector<double>& x,
                                  const std::vector<double>& v,
                                  bool rightmost_closed,
                                  bool all_inside,
                                  bool left_open) {
  std::vector<size_t> out(x.size());
  const double* v_begin = v.data();
  const double* v_end   = v_begin + v.size();
  const size_t n = x.size();
  const size_t nv = v.size();

  for (size_t i = 0; i < n; ++i) {
    double xi = x[i];
    const double* pos = left_open ? std::lower_bound(v_begin, v_end, xi) :
      std::upper_bound(v_begin, v_end, xi);
    size_t idx = static_cast<size_t>(pos - v_begin);
    if (rightmost_closed) {
      if (left_open) {
        if (nv > 0 && xi == v[0]) idx = 1;
      } else {
        if (nv > 0 && xi == v[nv - 1]) idx = nv - 1;
      }
    }
    if (all_inside) {
      if (idx == 0) idx = 1;
      else if (idx == nv) idx = nv - 1;
    }
    out[i] = idx;
  }
  return out;
}


double extract_sum(const DataFrameCpp& df, const char* name) {
  auto vec = df.get<double>(name);
  return std::accumulate(vec.begin(), vec.end(), 0.0);
};



double dtpwexpcpp1(
    const double q,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound,
    const bool logd) {
  double d;
  if (q <= lowerBound) {
    d = 0.0;
  } else {
    double v;
    size_t i0 = findInterval1(lowerBound, piecewiseSurvivalTime) - 1;
    size_t i1 = findInterval1(q, piecewiseSurvivalTime) - 1;
    if (i0 == i1) {
      v = lambda[i0] * (q - lowerBound);
    } else {
      v = lambda[i0] * (piecewiseSurvivalTime[i0 + 1] - lowerBound);
      for (size_t j = i0 + 1; j < i1; ++j) {
        v += lambda[j] * (piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j]);
      }
      v += lambda[i1] * (q - piecewiseSurvivalTime[i1]);
    }
    d = lambda[i1] * std::exp(-v);
  }
  if (logd) d = std::log(d);
  return d;
}


// [[Rcpp::export]]
std::vector<double> dtpwexpcpp(
    const std::vector<double>& q,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound,
    const bool logd) {
  size_t m = piecewiseSurvivalTime.size();
  size_t i0 = findInterval1(lowerBound, piecewiseSurvivalTime) - 1;
  size_t m1 = m - i0;
  std::vector<double> ch(m1, 0.0);
  for (size_t j = 1; j < m1; ++j) {
    size_t jj = j + i0 - 1;
    double lb = (j == 1) ? lowerBound : piecewiseSurvivalTime[jj];
    ch[j] = ch[j-1] + lambda[jj] * (piecewiseSurvivalTime[jj+1] - lb);
  }

  size_t n = q.size();
  std::vector<double> d(n);

  for (size_t h = 0; h < n; ++h) {
    if (q[h] <= lowerBound) {
      d[h] = 0.0;
    } else {
      size_t i1 = findInterval1(q[h], piecewiseSurvivalTime) - 1;
      size_t j1 = i1 + i0;
      double lb = (j1 == i0) ? lowerBound : piecewiseSurvivalTime[j1];
      double ch1 = ch[i1] + lambda[i1] * (q[h] - lb);
      d[h] = lambda[i1] * std::exp(-ch1);
    }
  }

  if (logd) {
    for (size_t h = 0; h < n; ++h) {
      d[h] = std::log(d[h]);
    }
  }

  return d;
}


double ptpwexpcpp1(
    const double q,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound,
    const bool lowertail,
    const bool logp) {
  double p;
  if (q <= lowerBound) {
    p = 0.0;
  } else {
    size_t i0 = findInterval1(lowerBound, piecewiseSurvivalTime);
    size_t i1 = findInterval1(q, piecewiseSurvivalTime);
    double v;
    if (i0 == i1) {
      v = lambda[i0 - 1] * (q - lowerBound);
    } else {
      v = lambda[i0 - 1] * (piecewiseSurvivalTime[i0] - lowerBound);
      for (size_t j = i0; j < i1 - 1; ++j) {
        v += lambda[j] * (piecewiseSurvivalTime[j + 1] -
          piecewiseSurvivalTime[j]);
      }
      v += lambda[i1 - 1] * (q - piecewiseSurvivalTime[i1 - 1]);
    }
    p = 1.0 - std::exp(-v);
  }

  if (!lowertail) p = 1.0 - p;
  if (logp) p = std::log(p);

  return p;
}


// [[Rcpp::export]]
std::vector<double> ptpwexpcpp(
    const std::vector<double>& q,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound,
    const bool lowertail,
    const bool logp) {
  size_t m = piecewiseSurvivalTime.size();
  size_t i0 = findInterval1(lowerBound, piecewiseSurvivalTime) - 1;
  size_t m1 = m - i0;
  std::vector<double> ch(m1, 0.0);
  for (size_t j = 1; j < m1; ++j) {
    size_t jj = j + i0 - 1;
    double lb = (j == 1) ? lowerBound : piecewiseSurvivalTime[jj];
    ch[j] = ch[j-1] + lambda[jj] * (piecewiseSurvivalTime[jj+1] - lb);
  }

  size_t n = q.size();
  std::vector<double> p(n);

  for (size_t h = 0; h < n; ++h) {
    if (q[h] <= lowerBound) {
      p[h] = 0.0;
    } else {
      size_t i1 = findInterval1(q[h], piecewiseSurvivalTime) - 1;
      size_t j1 = i1 + i0;
      double lb = (j1 == i0) ? lowerBound : piecewiseSurvivalTime[j1];
      double ch1 = ch[i1] + lambda[i1] * (q[h] - lb);
      p[h] = 1.0 - std::exp(-ch1);
    }
  }

  if (!lowertail) {
    for (size_t h = 0; h < n; ++h) {
      p[h] = 1.0 - p[h];
    }
  }

  if (logp) {
    for (size_t h = 0; h < n; ++h) {
      p[h] = std::log(p[h]);
    }
  }

  return p;
}


// [[Rcpp::export]]
std::vector<double> qtpwexpcpp(
    const std::vector<double>& p,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound,
    const bool lowertail,
    const bool logp) {
  size_t n = p.size();
  std::vector<double> q(n);
  for (size_t h = 0; h < n; ++h) {
    q[h] = qtpwexpcpp1(p[h], piecewiseSurvivalTime, lambda,
                       lowerBound, lowertail, logp);
  }

  return q;
}


// mean and variance of a truncated piecewise exponential distribution
ListCpp mtpwexpcpp(
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound) {

  if (piecewiseSurvivalTime[0] != 0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");
  if (!none_na(lambda)) throw std::invalid_argument("lambda must be provided");
  if (lambda.size() != piecewiseSurvivalTime.size())
    throw std::invalid_argument("Invalid length for lambda");
  for (double v : lambda) {
    if (v < 0) throw std::invalid_argument("lambda must be nonnegative");
  }
  if (lowerBound < 0) throw std::invalid_argument("lowerBound must be nonnegative");

  size_t m = piecewiseSurvivalTime.size();
  if (lambda[m-1] == 0.0) {
    throw std::invalid_argument("The last hazard rate must be positive");
  }

  size_t i0 = findInterval1(lowerBound, piecewiseSurvivalTime) - 1;

  double s1 = 0.0, s2 = 0.0, v = 0.0;
  for (size_t j = i0; j < m; ++j) {
    double lb = (j == i0) ? lowerBound : piecewiseSurvivalTime[j];
    double ub = (j == m-1) ? POS_INF : piecewiseSurvivalTime[j+1];

    if (j > i0) {
      double lb0 = (j == i0 + 1) ? lowerBound : piecewiseSurvivalTime[j-1];
      v += lambda[j-1] * (piecewiseSurvivalTime[j] - lb0);
    }

    if (lambda[j] == 0.0) {
      s1 += std::exp(-v) * (ub - lb);
      s2 += std::exp(-v) * (ub*ub - lb*lb);
    } else {
      double a1, a2;
      if (j < m - 1) {
        double a0 = std::exp(-lambda[j] * (ub - lb));
        a1 = 1.0 - a0;
        a2 = (1.0 + lambda[j] * lb) - (1.0 + lambda[j] * ub) * a0;
      } else {
        a1 = 1.0;
        a2 = 1.0 + lambda[j] * lb;
      }

      double ilam = 1.0 / lambda[j];
      s1 += std::exp(-v) * ilam * a1;
      s2 += std::exp(-v) * 2.0 * ilam * ilam * a2;
    }
  }

  double m1 = s1 + lowerBound;
  double m2 = s2 + lowerBound * lowerBound;

  ListCpp result;
  result.push_back(m1, "mean");
  result.push_back(m2 - m1 * m1, "variance");
  return result;
}


//' @title Mean and Variance of Truncated Piecewise Exponential Distribution
//' @description Obtains the mean and variance from a truncated piecewise
//' exponential distribution.
//'
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @param lowerBound The left truncation time point for the survival time.
//'   Defaults to 0 for no truncation.
//'
//' @return A list with two components, one for the mean, and the other for
//' the variance of the truncated piecewise exponential distribution.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' mtpwexp(piecewiseSurvivalTime = c(0, 6, 9, 15),
//'         lambda = c(0.025, 0.04, 0.015, 0.007))
//'
//' @export
//  [[Rcpp::export]]
Rcpp::List mtpwexp(
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& lambda = NA_REAL,
    const double lowerBound = 0) {

  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lam = Rcpp::as<std::vector<double>>(lambda);

  auto out = mtpwexpcpp(pwSurvT, lam, lowerBound);
  return Rcpp::wrap(out);
}


// Integrate a function f(theta) with respect to a normal density of theta:
//   integrate[f(theta) * dnorm(theta, mu, sigma), a, b) /
//             (pnorm(b, mu, sigma) - pnorm(a, mu, sigma)), {theta, a, b}].
double intnorm(const std::function<double(double)>& f,
               double mu, double sigma, double a, double b) {

  size_t r = 18, r1 = 6 * r - 1;
  double a1 = (a - mu) / sigma , b1 = (b - mu) / sigma;

  std::vector<double> x1(r1);
  for (size_t i = 0; i < r1; ++i) {
    if (i < r - 1) {
      x1[i] = -3.0 - 4.0 * std::log(static_cast<double>(r) / (i + 1.0));
    } else if (i < 5 * r) {
      x1[i] = -3.0 + 3.0 * ( (i + 1.0 - r) / (2.0 * r) );
    } else {
      x1[i] = 3.0 + 4.0 * std::log(static_cast<double>(r) / (6.0 * r - i - 1.0));
    }
  }

  // trim off x values outside (a1, b1)
  // trim from below
  auto it1 = std::upper_bound(x1.begin(), x1.end(), a1);
  size_t i1 = static_cast<size_t>(it1 - x1.begin()); // x1[i1-1] <= a1 < x1[i1]

  // trim from above
  auto it2 = std::lower_bound(x1.begin(), x1.end(), b);
  size_t i2 = static_cast<size_t>(it2 - x1.begin()); // x1[i2-1] < b1 <= x1[i2]

  // save the trimmed portion to x
  size_t m1 = i2 - i1 + 2;
  std::vector<double> x(m1);
  x[0] = a1;
  x[m1 - 1] = b1;
  for (size_t i = 1; i < m1 - 1; ++i) x[i] = x1[i + i1 - 1];

  // derive the grid points for z
  size_t m = 2 * m1 - 1;
  std::vector<double> z(m), w(m);
  // odd points
  for (size_t i = 0; i < m1; ++i) z[2*i] = x[i];
  // even points (midpoints)
  for (size_t i = 0; i < m1 - 1; ++i) z[2*i + 1] = 0.5 * (z[2*i] + z[2*i + 2]);

  // weights w as Simpson-like composite rule (same formulas)
  // first weight
  w[0] = (z[2] - z[0]) / 6.0;

  // interior even indices (i = 2, 4, ..., 2*(m1-2))
  for (size_t i0 = 1; i0 < m1 - 1; ++i0) {
    size_t i = 2 * i0;
    w[i] = (z[i + 2] - z[i - 2]) / 6.0;
  }
  // interior odd indices (i = 1,3,...)
  for (size_t i0 = 1; i0 < m1; ++i0) {
    size_t i = 2 * i0 - 1;
    w[i] = 4.0 * (z[i + 1] - z[i - 1]) / 6.0;
  }
  // last weight
  w[m - 1] = (z[m - 1] - z[m - 3]) / 6.0;

  // integrate
  double aval = 0.0;
  for (size_t i = 0; i < m; ++i) {
    aval += w[i] * f(mu + sigma * z[i]) * boost_dnorm(z[i]);
  }

  double denom = boost_pnorm(b1) - boost_pnorm(a1);

  return aval / denom;
}


std::pair<double, double> mini(const std::function<double(double)>& f,
                               double x1, double x2) {

  // Validate inputs
  if (!(x1 < x2)) throw std::invalid_argument("mini: require x1 < x2");

  // Use Boost's brent_find_minima
  return boost::math::tools::brent_find_minima(
    f, x1, x2, std::numeric_limits<double>::digits10);
}


// Integrate f over multiple intervals defined by breaks, i.e. sum of integrals
double integrate3(const std::function<double(double)>& f,
                  const std::vector<double>& breaks, double tol) {

  boost::math::quadrature::tanh_sinh<double> integrator_ts;

  double sum = 0.0;
  for (size_t i = 0; i + 1 < breaks.size(); ++i) {
    double a = breaks[i];
    double b = breaks[i+1];
    if (b <= a) continue;
    double val = integrator_ts.integrate(f, a, b, tol);
    sum += val;
  }
  return sum;
}



// 2D adaptive quadrature using 3x3 and 5x5 tensor-product Gauss rules
// 1D Gauss nodes
static const double x3[3] = { -0.774596669241483,
                              0.0,
                              0.774596669241483 };

static const double w3[3] = { 0.555555555555556,
                              0.888888888888889,
                              0.555555555555556 };

static const double x5[5] = { -0.906179845938664,
                              -0.538469310105683,
                              0.0,
                              0.538469310105683,
                              0.906179845938664 };

static const double w5[5] = { 0.236926885056189,
                              0.478628670499366,
                              0.568888888888889,
                              0.478628670499366,
                              0.236926885056189 };

// structure to store a region
struct Region {
  double ax, bx, ay, by;
  double I_low;    // 3x3 rule
  double I_high;   // 5x5 rule
  double err;      // |I_high - I_low|
};

// compute 3x3 and 5x5 tensor-product Gauss rules on a region
static void eval_region(const std::function<double(double,double)>& f,
                        double ax, double bx, double ay, double by,
                        double &I3, double &I5) {
  double cx = 0.5*(ax + bx), dx = 0.5*(bx - ax);
  double cy = 0.5*(ay + by), dy = 0.5*(by - ay);

  I3 = 0.0;
  I5 = 0.0;

  // Store f at Gauss nodes to reuse in Kronrod
  double f3[3][3];

  // 3x3 Gauss
  for(size_t i=0;i<3;++i){
    double xi = cx + dx * x3[i];
    for(size_t j=0;j<3;++j){
      double yj = cy + dy * x3[j];
      f3[i][j] = f(xi,yj);
      I3 += w3[i]*w3[j]*f3[i][j];
    }
  }
  I3 *= dx*dy;

  // 5×5 rule
  for(size_t i=0;i<5;++i){
    double xi = cx + dx * x5[i];
    for(size_t j=0;j<5;++j){
      double yj = cy + dy * x5[j];
      if (i==2 && j==2) {
        I5 += w5[i]*w5[j] * f3[1][1];
      } else {
        I5 += w5[i]*w5[j] * f(xi,yj);
      }
    }
  }
  I5 *= dx * dy;
}

// priority comparison: larger error first
struct RegionCompare {
  bool operator()(const Region& a, const Region& b) const {
    return a.err < b.err;   // max-heap
  }
};

double quad2d(const std::function<double(double,double)>& f,
              double ax, double bx, double ay, double by, double tol) {
  size_t maxRegions = 1000000;

  double I3, I5;
  eval_region(f, ax, bx, ay, by, I3, I5);

  double global_value = I5;
  double global_err   = std::fabs(I5 - I3);

  std::priority_queue<Region, std::vector<Region>, RegionCompare> pq;

  pq.push({ax, bx, ay, by, I3, I5, global_err});

  while(global_err > tol && global_err > std::fabs(global_value) * tol &&
        pq.size() < maxRegions)
  {
    Region R = pq.top();
    pq.pop();

    // split in x direction (larger span)
    double dx = R.bx - R.ax;
    double dy = R.by - R.ay;
    bool splitX = (dx >= dy);

    Region R1 = R, R2 = R;
    if(splitX){
      double m = 0.5*(R.ax + R.bx);
      R1.bx = m; R2.ax = m;
    } else {
      double m = 0.5*(R.ay + R.by);
      R1.by = m; R2.ay = m;
    }

    eval_region(f, R1.ax, R1.bx, R1.ay, R1.by, R1.I_low, R1.I_high);
    eval_region(f, R2.ax, R2.bx, R2.ay, R2.by, R2.I_low, R2.I_high);

    R1.err = std::fabs(R1.I_high - R1.I_low);
    R2.err = std::fabs(R2.I_high - R2.I_low);

    // update global integral & error incrementally
    global_value += (R1.I_high + R2.I_high) - R.I_high;
    global_err   += (R1.err + R2.err) - R.err;

    pq.push(R1);
    pq.push(R2);
  }

  return global_value;
}


double pbvnormcpp(const std::vector<double>& lower,
                  const std::vector<double>& upper,
                  const double rho) {
  if (!none_na(lower)) throw std::invalid_argument("lower must be provided");
  if (!none_na(upper)) throw std::invalid_argument("upper must be provided");

  auto lowerv = expand1(lower, 2, "lower");
  auto upperv = expand1(upper, 2, "upper");
  for (size_t i = 0; i < 2; ++i) {
    if (lowerv[i] >= upperv[i])
      throw std::invalid_argument("lower must be less than upper");
  }

  if (rho >= 1 || rho <= -1)
    throw std::invalid_argument("corr must lie between -1 and 1");

  if (rho == 0.0) {
    double v1 = boost_pnorm(upperv[0]) - boost_pnorm(lowerv[0]);
    double v2 = boost_pnorm(upperv[1]) - boost_pnorm(lowerv[1]);
    return v1 * v2;
  }

  double a2 = lowerv[1];
  double b2 = upperv[1];
  double s = std::sqrt(1.0 - rho * rho);
  auto f = [&](double x)->double {
    double a = (a2 - rho * x) / s;
    double b = (b2 - rho * x) / s;
    double t1 = boost_dnorm(x);
    double t2 = boost_pnorm(b) - boost_pnorm(a);
    return t1 * t2;
  };
  std::vector<double> breaks = {lowerv[0], upperv[0]};
  return integrate3(f, breaks, 1e-8);
}



//' @title Distribution Function of the Standard Bivariate Normal
//' @description Computes the cumulative distribution function (CDF) of
//' the standard bivariate normal distribution with specified lower and
//' upper integration limits and correlation coefficient.
//'
//' @param lower A numeric vector of length 2 specifying the lower limits
//'   of integration.
//' @param upper A numeric vector of length 2 specifying the upper limits
//'   of integration.
//' @param rho A numeric value specifying the correlation coefficient of
//'   the standard bivariate normal distribution.
//'
//' @details This function evaluates the probability
//' \eqn{P(\code{lower[1]} < X < \code{upper[1]},
//' \code{lower[2]} < Y < \code{upper[2]})} where
//' \eqn{(X, Y)} follows a standard bivariate normal
//' distribution with correlation \code{corr}.
//'
//' @return A numeric value representing the probability that a standard
//' bivariate normal vector falls within the specified rectangular region.
//'
//' @examples
//' pbvnorm(c(-1, -1), c(1, 1), 0.5)
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @export
// [[Rcpp::export]]
double pbvnorm(
    const Rcpp::NumericVector& lower = NA_REAL,
    const Rcpp::NumericVector& upper = NA_REAL,
    const double rho = 0) {

  auto lowerv = Rcpp::as<std::vector<double>>(lower);
  auto upperv = Rcpp::as<std::vector<double>>(upper);
  return pbvnormcpp(lowerv, upperv, rho);
}




ListCpp hazard_pdcpp(const std::vector<double>& piecewiseSurvivalTime,
                     const std::vector<double>& hazard_pfs,
                     const std::vector<double>& hazard_os,
                     const double rho_pd_os) {

  if (piecewiseSurvivalTime[0] != 0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");
  size_t n = piecewiseSurvivalTime.size();

  if (!none_na(hazard_pfs)) throw std::invalid_argument("hazard_pfs must be provided");
  if (!none_na(hazard_os)) throw std::invalid_argument("hazard_os must be provided");
  for (double v : hazard_pfs) {
    if (v <= 0) throw std::invalid_argument("hazard_pfs must be positive");
  }
  for (double v : hazard_os) {
    if (v <= 0) throw std::invalid_argument("hazard_os must be positive");
  }
  auto hazard_pfsv = expand1(hazard_pfs, n, "hazard_pfs");
  auto hazard_osv = expand1(hazard_os, n, "hazard_os");
  for (size_t i = 0; i < n; ++i) {
    if (hazard_pfsv[i] <= hazard_osv[i]) {
      throw std::invalid_argument("hazard_pfs must be greater than hazard_os");
    }
  }
  if (rho_pd_os <= -1 || rho_pd_os >= 1)
    throw std::invalid_argument("corr_pd_os must lie between -1 and 1");

  // append additional time points for pfs quantiles
  std::vector<double> p(10);
  for (size_t i=0; i<9; ++i) p[i] = (i+1.0)/10.0;
  p[9] = 0.95;

  std::vector<double> u(n+10);
  for (size_t i=0; i<n-1; ++i) u[i] = piecewiseSurvivalTime[i+1];
  u[n-1] = piecewiseSurvivalTime[n-1] + std::log(2.0)/hazard_pfsv[n-1];
  for (size_t i=0; i<10; ++i) {
    u[n+i] = qtpwexpcpp1(p[i], piecewiseSurvivalTime, hazard_pfsv);
  }

  // obtain sorted and unique time points
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  size_t m = u.size();

  // shifted time points
  std::vector<double> u1(m);
  u1[0] = 0;
  for (size_t i=1; i<m; ++i) u1[i] = u[i-1];

  // get corresponding hazards
  std::vector<size_t> index = findInterval3(u1, piecewiseSurvivalTime);
  for (size_t i=0; i<m; ++i) index[i] = index[i] - 1;
  std::vector<double> hazard_pfs1 = subset(hazard_pfsv, index);
  std::vector<double> hazard_os1 = subset(hazard_osv, index);

  // solve for hazard_pd
  double t;
  std::vector<double> hazard_pd(m);
  std::vector<double> v(0), hazard(0), haz_pfs(0), haz_os(0);
  auto f = [&](double haz)->double {
    std::vector<double> haz_pd = hazard;
    haz_pd.push_back(haz);
    std::vector<double> lower(2);
    double a = ptpwexpcpp1(t, v, haz_pd, 0, 1, 0);
    double b = ptpwexpcpp1(t, v, haz_os, 0, 1, 0);
    lower[0] = boost_qnorm(a);
    lower[1] = boost_qnorm(b);
    std::vector<double> upper(2, POS_INF);
    double q = pbvnormcpp(lower, upper, rho_pd_os);
    return q - ptpwexpcpp1(t, v, haz_pfs, 0, 0, 0);
  };

  double tol = 1e-6;
  for (size_t i=0; i<m; ++i) {
    t = u[i];
    v.push_back(u1[i]);
    haz_pfs.push_back(hazard_pfs1[i]);
    haz_os.push_back(hazard_os1[i]);
    hazard_pd[i] = brent(f, 0.5 * (hazard_pfs1[i] - hazard_os1[i]),
                         2.0 * hazard_pfs1[i], tol);
    hazard.push_back(hazard_pd[i]);
  }

  ListCpp result;
  result.push_back(std::move(u1), "piecewiseSurvivalTime");
  result.push_back(std::move(hazard_pd), "hazard_pd");
  result.push_back(std::move(hazard_os1), "hazard_os");
  result.push_back(rho_pd_os, "rho_pd_os");
  return result;
}


//' @title Hazard Function for Progressive Disease (PD) Given Correlation
//' Between PD and OS
//'
//' @description
//' Computes the hazard function of a piecewise exponential
//' distribution for progressive disease (PD), such that the
//' resulting hazard function for progression-free survival (PFS)
//' closely matches a given piecewise hazard for PFS.
//'
//' @inheritParams param_piecewiseSurvivalTime
//' @param hazard_pfs A scalar or numeric vector specifying the
//'   hazard(s) for PFS based on a piecewise exponential distribution.
//' @param hazard_os A scalar or numeric vector specifying the
//'   hazard(s) for overall survival (OS) based on a piecewise
//'   exponential distribution.
//' @param rho_pd_os A numeric value specifying the correlation
//'   between PD and OS times.
//'
//' @details
//' This function determines the hazard vector \eqn{\lambda_{\text{pd}}}
//' for the piecewise exponential distribution of PD, so that the
//' implied survival function for PFS time,
//' \eqn{T_{\text{pfs}} = \min(T_{\text{pd}}, T_{\text{os}})}, closely
//' matches the specified piecewise exponential distribution for PFS
//' with hazard vector \eqn{\lambda_{\text{pfs}}}.
//'
//' To achieve this, we simulate
//' \eqn{(Z_{\text{pd}}, Z_{\text{os}})} from
//' a standard bivariate normal distribution with correlation
//' \eqn{\rho}. Then, \eqn{U_{\text{pd}} = \Phi(Z_{\text{pd}})}
//' and \eqn{U_{\text{os}} = \Phi(Z_{\text{os}})} are generated, where
//' \eqn{\Phi} denotes the standard normal CDF.
//'
//' The times to PD and OS are obtained via the inverse transform
//' method using quantile functions of the piecewise exponential distribution:
//' \deqn{T_{\text{pd}} = \text{qpwexp}(U_{\text{pd}},u,\lambda_{\text{pd}})}
//' \deqn{T_{\text{os}} = \text{qpwexp}(U_{\text{os}},u,\lambda_{\text{os}})}
//' where \code{u = piecewiseSurvivalTime}.
//'
//' The function solves for \eqn{\lambda_{\text{pd}}} such that
//' the survival function of \eqn{T_{\text{pfs}}} closely matches that
//' of a piecewise exponential distribution with hazard
//' \eqn{\lambda_{\text{pfs}}}:
//' \deqn{P(\min(T_{\text{pd}}, T_{\text{os}}) > t) = S_{\text{pfs}}(t)}
//' Since \deqn{Z_{\text{pd}} =
//'   \Phi^{-1}(\text{ppwexp}(T_\text{pd}, u, \lambda_{\text{pd}}))} and
//' \deqn{Z_{\text{os}} =
//'   \Phi^{-1}(\text{ppwexp}(T_\text{os}, u, \lambda_{\text{os}}))}
//' we have
//' \deqn{P(\min(T_{\text{pd}}, T_{\text{os}}) > t) =
//' P(Z_{\text{pd}} >
//'     \Phi^{-1}(\text{ppwexp}(t,u,\lambda_{\text{pd}})),
//'   Z_{\text{os}} >
//'     \Phi^{-1}(\text{ppwexp}(t,u,\lambda_{\text{os}})))}
//' while
//' \deqn{S_{\text{pfs}}(t) = 1 - \text{ppwexp}(t,u,\lambda_{\text{pfs}})}
//'
//' Matching is performed sequentially at the internal cut points
//' \eqn{u_2, ..., u_J} and at the point
//' \eqn{u_J + \log(2)/\lambda_{\text{pfs},J}} for the final interval,
//' as well as the percentile points at 10%, 20%, ..., 90%, and 95%
//' to solve for \eqn{\lambda_{\text{pd},1}, \ldots,
//' \lambda_{\text{pd},K}}, where \eqn{K} is the total number of
//' unique cut points.
//'
//' @return A list with the following components:
//'
//' * \code{piecewiseSurvivalTime}: A vector that specifies the starting time
//'   points of the intervals for the piecewise exponential distribution
//'   for PD.
//'
//' * \code{hazard_pd}: A numeric vector representing the calculated hazard
//'   rates for the piecewise exponential distribution of PD.
//'
//' * \code{hazard_os}: A numeric vector representing the hazard rates for
//'   the piecewise exponential distribution of OS at the same time points
//'   as PD.
//'
//' * \code{rho_pd_os}: The correlation between PD and OS times (as input).
//'
//' @author
//' Kaifeng Lu (\email{kaifenglu@gmail.com})
//'
//' @examples
//' u <- c(0, 1, 3, 4)
//' lambda1 <- c(0.0151, 0.0403, 0.0501, 0.0558)
//' lambda2 <- 0.0145
//' rho_pd_os <- 0.5
//' hazard_pd(u, lambda1, lambda2, rho_pd_os)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List hazard_pd(
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& hazard_pfs = NA_REAL,
    const Rcpp::NumericVector& hazard_os = NA_REAL,
    const double rho_pd_os = 0.5) {

  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto haz_pfs = Rcpp::as<std::vector<double>>(hazard_pfs);
  auto haz_os = Rcpp::as<std::vector<double>>(hazard_os);

  auto out = hazard_pdcpp(pwSurvT, haz_pfs, haz_os, rho_pd_os);
  return Rcpp::wrap(out);
}



double pdf_pfs(const double time,
               const std::vector<double>& piecewiseSurvivalTime,
               const std::vector<double>& hazard_pd,
               const std::vector<double>& hazard_os,
               const double rho_pd_os) {
  double s = std::sqrt(1.0 - rho_pd_os * rho_pd_os);
  double u1 = ptpwexpcpp1(time, piecewiseSurvivalTime, hazard_pd);
  double u2 = ptpwexpcpp1(time, piecewiseSurvivalTime, hazard_os);
  double d1 = dtpwexpcpp1(time, piecewiseSurvivalTime, hazard_pd);
  double d2 = dtpwexpcpp1(time, piecewiseSurvivalTime, hazard_os);
  double z1 = boost_qnorm(u1);
  double z2 = boost_qnorm(u2);
  double a1 = 1.0 - boost_pnorm((z2 - rho_pd_os * z1) / s);
  double a2 = 1.0 - boost_pnorm((z1 - rho_pd_os * z2) / s);
  double result = a1 * d1 + a2 * d2;
  return result;
};

double sdf_pfs(const double time,
               const std::vector<double>& piecewiseSurvivalTime,
               const std::vector<double>& hazard_pd,
               const std::vector<double>& hazard_os,
               const double rho_pd_os) {
  double u1 = ptpwexpcpp1(time, piecewiseSurvivalTime, hazard_pd);
  double u2 = ptpwexpcpp1(time, piecewiseSurvivalTime, hazard_os);
  double z1 = boost_qnorm(u1);
  double z2 = boost_qnorm(u2);
  std::vector<double> lower = {z1, z2};
  std::vector<double> upper = {POS_INF, POS_INF};
  double p = pbvnormcpp(lower, upper, rho_pd_os);
  return p;
};

double upper_pfs(const std::vector<double>& piecewiseSurvivalTime,
                 const std::vector<double>& hazard_pd,
                 const std::vector<double>& hazard_os,
                 const double rho_pd_os) {
  double tol = 1e-12;
  double upper = 1.0;
  double p = sdf_pfs(upper, piecewiseSurvivalTime, hazard_pd, hazard_os, rho_pd_os);
  while (p > tol) {
    upper *= 2.0;
    p = sdf_pfs(upper, piecewiseSurvivalTime, hazard_pd, hazard_os, rho_pd_os);
  }
  return upper;
}

ListCpp m_pfs(const std::vector<double>& piecewiseSurvivalTime,
              const std::vector<double>& hazard_pd,
              const std::vector<double>& hazard_os,
              const double rho_pd_os) {
  double upper = upper_pfs(piecewiseSurvivalTime, hazard_pd, hazard_os, rho_pd_os);
  auto fm_pfs = [&](double t)->double {
    return t * pdf_pfs(t, piecewiseSurvivalTime, hazard_pd, hazard_os, rho_pd_os);
  };
  auto fm2_pfs = [&](double t)->double {
    return t*t * pdf_pfs(t, piecewiseSurvivalTime, hazard_pd, hazard_os, rho_pd_os);
  };

  double tol = 1e-5;

  // Candidate set of break points for integration
  std::vector<double> breaks;
  breaks.reserve(piecewiseSurvivalTime.size() + 1);
  for (double t : piecewiseSurvivalTime) {
    if (t >= upper) continue;
    breaks.push_back(t);
  }
  breaks.push_back(upper);
  double m1 = integrate3(fm_pfs, breaks, tol);
  double m2 = integrate3(fm2_pfs, breaks, tol);

  ListCpp result;
  result.push_back(m1, "mean");
  result.push_back(m2 - m1 * m1, "variance");
  return result;
}

double cor_pfs_os(const std::vector<double>& piecewiseSurvivalTime,
                  const std::vector<double>& hazard_pd,
                  const std::vector<double>& hazard_os,
                  const double rho_pd_os) {
  ListCpp mv1 = m_pfs(piecewiseSurvivalTime, hazard_pd, hazard_os, rho_pd_os);
  double m1 = mv1.get<double>("mean");
  double v1 = mv1.get<double>("variance");

  ListCpp mv2 = mtpwexpcpp(piecewiseSurvivalTime, hazard_os, 0);
  double m2 = mv2.get<double>("mean");
  double v2 = mv2.get<double>("variance");

  double s = std::sqrt(1.0 - rho_pd_os * rho_pd_os);
  double c1 = 1.0 / (2.0 * M_PI * s);
  double c2 = 1.0 / (2.0 * s * s);

  // integrand for the joint moment
  auto f = [&](double u1, double u2)->double {
    double t1 = qtpwexpcpp1(u1, piecewiseSurvivalTime, hazard_pd);
    double t2 = qtpwexpcpp1(u2, piecewiseSurvivalTime, hazard_os);
    double z1 = boost_qnorm(u1);
    double z2 = boost_qnorm(u2);
    double a1 = std::min(t1, t2) * t2;

    // joint density of standard bivariate normal
    double a2 = c1 * exp(-c2 * (z1*z1 - 2.0*rho_pd_os*z1*z2 + z2*z2));
    double a3 = boost_dnorm(z1) * boost_dnorm(z2);
    return a1 * a2 / a3;
  };

  double tol = 1e-4;
  double m12 = quad2d(f, 0.0, 1.0, 0.0, 1.0, tol);

  double cov = m12 - m1 * m2;
  return cov / sqrt(v1 * v2);
}


double corr_pfs_oscpp(const std::vector<double>& piecewiseSurvivalTime,
                      const std::vector<double>& hazard_pfs,
                      const std::vector<double>& hazard_os,
                      const double rho_pd_os) {

  if (piecewiseSurvivalTime[0] != 0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");
  size_t n = piecewiseSurvivalTime.size();

  if (!none_na(hazard_pfs)) throw std::invalid_argument("hazard_pfs must be provided");
  if (!none_na(hazard_os)) throw std::invalid_argument("hazard_os must be provided");
  for (double v : hazard_pfs) {
    if (v <= 0) throw std::invalid_argument("hazard_pfs must be positive");
  }
  for (double v : hazard_os) {
    if (v <= 0) throw std::invalid_argument("hazard_os must be positive");
  }
  auto hazard_pfsv = expand1(hazard_pfs, n, "hazard_pfs");
  auto hazard_osv = expand1(hazard_os, n, "hazard_os");
  for (size_t i = 0; i < n; ++i) {
    if (hazard_pfsv[i] <= hazard_osv[i]) {
      throw std::invalid_argument("hazard_pfs must be greater than hazard_os");
    }
  }
  if (rho_pd_os <= -1 || rho_pd_os >= 1)
    throw std::invalid_argument("corr_pd_os must lie between -1 and 1");

  ListCpp a = hazard_pdcpp(piecewiseSurvivalTime, hazard_pfsv, hazard_osv, rho_pd_os);
  std::vector<double> u = a.get<std::vector<double>>("piecewiseSurvivalTime");
  std::vector<double> hazard_pd1 = a.get<std::vector<double>>("hazard_pd");
  std::vector<double> hazard_os1 = a.get<std::vector<double>>("hazard_os");

  return cor_pfs_os(u, hazard_pd1, hazard_os1, rho_pd_os);
}


//' @title Correlation Between PFS and OS Given Correlation Between PD and OS
//'
//' @description
//' Computes the correlation between PFS and OS given the correlation
//' between PD and OS.
//'
//' @inheritParams param_piecewiseSurvivalTime
//' @param hazard_pfs A scalar or numeric vector specifying the
//'   hazard(s) for PFS based on a piecewise exponential distribution.
//' @param hazard_os A scalar or numeric vector specifying the
//'   hazard(s) for overall survival (OS) based on a piecewise
//'   exponential distribution.
//' @param rho_pd_os A numeric value specifying the correlation
//'   between PD and OS times.
//'
//' @details
//' This function first determines the piecewise exponential distribution
//' for PD such that the implied survival function for PFS time,
//' \eqn{T_{\text{pfs}} = \min(T_{\text{pd}}, T_{\text{os}})}, closely
//' matches the specified piecewise exponential distribution for PFS
//' with hazard vector \eqn{\lambda_{\text{pfs}}}. Then, it calculates
//' the correlation between PFS and OS times based on the derived
//' piecewise exponential distribution for PD and the given piecewise
//' exponential distribution for OS.
//'
//' @return The estimated correlation between PFS and OS.
//'
//' @author
//' Kaifeng Lu (\email{kaifenglu@gmail.com})
//'
//' @examples
//' u <- c(0, 1, 3, 4)
//' lambda1 <- c(0.0151, 0.0403, 0.0501, 0.0558)
//' lambda2 <- 0.0145
//' rho_pd_os <- 0.5
//' corr_pfs_os(u, lambda1, lambda2, rho_pd_os)
//'
//' @export
// [[Rcpp::export]]
double corr_pfs_os(
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& hazard_pfs = NA_REAL,
    const Rcpp::NumericVector& hazard_os = NA_REAL,
    const double rho_pd_os = NA_REAL) {

  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto haz_pfs = Rcpp::as<std::vector<double>>(hazard_pfs);
  auto haz_os = Rcpp::as<std::vector<double>>(hazard_os);

  return corr_pfs_oscpp(pwSurvT, haz_pfs, haz_os, rho_pd_os);
}



ListCpp hazard_subcpp(const std::vector<double>& piecewiseSurvivalTime,
                      const std::vector<double>& hazard_itt,
                      const std::vector<double>& hazard_pos,
                      const double p_pos) {

  if (piecewiseSurvivalTime[0] != 0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");
  if (!none_na(hazard_itt)) throw std::invalid_argument("hazard_itt must be provided");
  if (!none_na(hazard_pos)) throw std::invalid_argument("hazard_pos must be provided");
  for (double v : hazard_itt) {
    if (v <= 0) throw std::invalid_argument("hazard_itt must be positive");
  }
  for (double v : hazard_pos) {
    if (v <= 0) throw std::invalid_argument("hazard_pos must be positive");
  }
  if (std::isnan(p_pos)) throw std::invalid_argument("p_pos must be provided");
  if (p_pos <= 0 || p_pos >= 1)
    throw std::invalid_argument("p_pos must lie between 0 and 1");

  size_t n = piecewiseSurvivalTime.size();
  auto hazard_ittv = expand1(hazard_itt, n, "hazard_itt");
  auto hazard_posv = expand1(hazard_pos, n, "hazard_pos");

  // append additional time points for pfs quantiles
  std::vector<double> p(10);
  for (size_t i=0; i<9; ++i) p[i] = (i+1.0)/10.0;
  p[9] = 0.95;

  std::vector<double> u(n+10);
  for (size_t i=0; i<n-1; ++i) u[i] = piecewiseSurvivalTime[i+1];
  u[n-1] = piecewiseSurvivalTime[n-1] + std::log(2.0)/hazard_ittv[n-1];
  for (size_t i=0; i<10; ++i) {
    u[n+i] = qtpwexpcpp1(p[i], piecewiseSurvivalTime, hazard_ittv);
  }

  // obtain sorted and unique time points
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  size_t m = u.size();

  // shifted time points
  std::vector<double> u1(m);
  u1[0] = 0;
  for (size_t i=1; i<m; ++i) u1[i] = u[i-1];

  // get corresponding hazards
  std::vector<size_t> index = findInterval3(u1, piecewiseSurvivalTime);
  for (size_t i=0; i<m; ++i) index[i] = index[i] - 1;
  std::vector<double> hazard_itt1 = subset(hazard_ittv, index);
  std::vector<double> hazard_pos1 = subset(hazard_posv, index);

  // solve for hazard_sub
  double t;
  std::vector<double> hazard_neg(m);
  std::vector<double> v(0), hazard(0), haz_itt(0), haz_pos(0);
  auto f = [&](double haz)->double {
    std::vector<double> haz_neg = hazard;
    haz_neg.push_back(haz);
    double a = ptpwexpcpp1(t, v, haz_pos, 0, 1, 0);
    double b = ptpwexpcpp1(t, v, haz_neg, 0, 1, 0);
    double q = p_pos * a + (1.0 - p_pos) * b;
    return q - ptpwexpcpp1(t, v, haz_itt, 0, 1, 0);
  };

  double tol = 1e-6;
  for (size_t i=0; i<m; ++i) {
    t = u[i];
    v.push_back(u1[i]);
    haz_itt.push_back(hazard_itt1[i]);
    haz_pos.push_back(hazard_pos1[i]);
    hazard_neg[i] = brent(f, 0.5 * (hazard_itt1[i] - p_pos * hazard_pos1[i]),
                          2.0 * hazard_itt1[i] / (1.0 - p_pos), tol);
    hazard.push_back(hazard_neg[i]);
  }

  ListCpp result;
  result.push_back(std::move(u1), "piecewiseSurvivalTime");
  result.push_back(std::move(hazard_pos1), "hazard_pos");
  result.push_back(std::move(hazard_neg), "hazard_neg");
  result.push_back(p_pos, "p_pos");
  return result;
}


//' @title Hazard Function for Sub Population
//'
//' @description
//' Computes the hazard function of a piecewise exponential
//' distribution for the biomarker negative sub population, such that the
//' resulting survival function for the ITT population
//' closely matches a given piecewise survival function.
//'
//' @inheritParams param_piecewiseSurvivalTime
//' @param hazard_itt A scalar or numeric vector specifying the
//'   hazard(s) for the ITT population based on a piecewise exponential
//'   distribution.
//' @param hazard_pos A scalar or numeric vector specifying the
//'   hazard(s) for the biomarker positive sub population
//'   based on a piecewise exponential distribution.
//' @param p_pos A numeric value specifying the prevalence of the
//'   biomarker positive sub population.
//'
//' @details
//' This function determines the hazard vector \eqn{\lambda_{\text{neg}}}
//' for the piecewise exponential distribution of the biomarker negative
//' sub population, so that the implied survival function for the ITT
//' population closely matches the specified piecewise exponential
//' distribution with hazard vector \eqn{\lambda_{\text{itt}}}.
//'
//' Let \eqn{p_{\text{pos}}} be the
//' prevalence of the biomarker positive sub population,
//' then the survival function for the ITT population is given by
//' \deqn{S_{\text{itt}}(t) = p_{\text{pos}} S_{\text{pos}}(t) +
//' (1 - p_{\text{pos}}) S_{\text{neg}}(t)}
//' where \eqn{S_{\text{pos}}(t)} and \eqn{S_{\text{neg}}(t)} are
//' the survival functions for the biomarker positive and
//' biomarker negative sub populations, respectively.
//'
//' Matching is performed sequentially at the internal cutpoints
//' \eqn{u_2, ..., u_J} and at the point
//' \eqn{u_J + \log(2)/\lambda_{\text{itt},J}} for the final interval,
//' as well as the percentile points at 10%, 20%, ..., 90%, and 95%,
//' to solve for \eqn{\lambda_{\text{neg},1}, \ldots,
//' \lambda_{\text{neg},K}}, where \eqn{K} is the total number of
//' unique cut points.
//'
//' @return A list with the following components:
//'
//' * \code{piecewiseSurvivalTime}: A vector that specifies the starting time
//'   points of the intervals for the piecewise exponential distribution
//'   for the biomarker negative sub population.
//'
//' * \code{hazard_pos}: A numeric vector representing the hazard rates for
//'   the piecewise exponential distribution of the biomarker positive
//'   sub population at the same time points as the biomarker negative
//'   sub population.
//'
//' * \code{hazard_neg}: A numeric vector representing the estimated hazard
//'   rates for the piecewise exponential distribution of the biomarker
//'   negative sub population.
//'
//' * \code{p_pos}: The prevalence of the biomarker positive sub population
//'  (as input).
//'
//' @author
//' Kaifeng Lu (\email{kaifenglu@gmail.com})
//'
//' @examples
//' u <- c(0, 1, 3, 4)
//' lambda_itt <- c(0.0151, 0.0403, 0.0501, 0.0558)
//' lambda_pos <- c(0.0115, 0.0302, 0.0351, 0.0404)
//' p_pos <- 0.3
//' hazard_sub(u, lambda_itt, lambda_pos, p_pos)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List hazard_sub(
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& hazard_itt = NA_REAL,
    const Rcpp::NumericVector& hazard_pos = NA_REAL,
    const double p_pos = NA_REAL) {

  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto haz_itt = Rcpp::as<std::vector<double>>(hazard_itt);
  auto haz_pos = Rcpp::as<std::vector<double>>(hazard_pos);

  auto out = hazard_subcpp(pwSurvT, haz_itt, haz_pos, p_pos);
  return Rcpp::wrap(out);
}


// bygroup: group-by helper that builds lookup tables and combined indices
ListCpp bygroup(const DataFrameCpp& data,
                const std::vector<std::string>& variables) {

  size_t n = data.nrows();
  size_t p = variables.size();
  ListCpp result;
  std::vector<int> nlevels(p);

  // IntMatrix for indices (n rows, p cols), column-major storage
  IntMatrix indices(n, p);

  // Flattened lookup buffers and per-variable metadata
  struct VarLookupInfo {
    int type; // 0=int, 1=double, 2=bool, 3=string
    size_t offset;
  };
  std::vector<VarLookupInfo> var_info(p);

  std::vector<int> int_flat;
  std::vector<double> dbl_flat;
  std::vector<unsigned char> bool_flat;
  std::vector<std::string> str_flat;

  ListCpp lookups_per_variable; // will contain a std::vector for each variable

  for (size_t i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    if (!data.containElementNamed(var))
      throw std::invalid_argument("Data must contain variable: " + var);

    if (data.int_cols.count(var)) {
      const auto& col = data.int_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w); // indices 0..(levels-1)

      // append w to flat buffer and record metadata
      size_t off = int_flat.size();
      int_flat.insert(int_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{0, off};

      // fill indices column i (column-major layout)
      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.numeric_cols.count(var)) {
      const auto& col = data.numeric_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);

      size_t off = dbl_flat.size();
      dbl_flat.insert(dbl_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{1, off};

      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.bool_cols.count(var)) {
      const auto& col = data.bool_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);

      size_t off = bool_flat.size();
      bool_flat.insert(bool_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{2, off};

      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.string_cols.count(var)) {
      const auto& col = data.string_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);

      size_t off = str_flat.size();
      str_flat.insert(str_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{3, off};

      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else {
      throw std::invalid_argument("Unsupported variable type in bygroup: " + var);
    }
  } // end for variables

  // compute combined index
  std::vector<int> combined_index(n, 0);
  int orep = 1;
  for (size_t i = 0; i < p; ++i) orep *= nlevels[i];
  size_t lookup_nrows = orep;

  for (size_t i = 0; i < p; ++i) {
    orep /= nlevels[i];
    const int* col_ptr = indices.data_ptr() + i * n;
    for (size_t j = 0; j < n; ++j) {
      combined_index[j] += col_ptr[j] * orep;
    }
  }

  // Build lookup_df with columns repeated in the same pattern as before.
  DataFrameCpp lookup_df;
  size_t repeat_each = lookup_nrows;
  for (size_t i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    size_t nlevels_i = nlevels[i];
    repeat_each /= nlevels_i;
    size_t times = lookup_nrows / ( nlevels_i * repeat_each );

    VarLookupInfo info = var_info[i];
    if (info.type == 0) {
      const int* base = int_flat.data() + info.offset;
      std::vector<int> col(lookup_nrows);
      size_t idxw = 0;
      for (size_t t = 0; t < times; ++t) {
        for (size_t level = 0; level < nlevels_i; ++level)
          for (size_t r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    } else if (info.type == 1) {
      const double* base = dbl_flat.data() + info.offset;
      std::vector<double> col(lookup_nrows);
      size_t idxw = 0;
      for (size_t t = 0; t < times; ++t) {
        for (size_t level = 0; level < nlevels_i; ++level)
          for (size_t r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    } else if (info.type == 2) {
      const unsigned char* base = bool_flat.data() + info.offset;
      std::vector<unsigned char> col(lookup_nrows);
      size_t idxw = 0;
      for (size_t t = 0; t < times; ++t) {
        for (size_t level = 0; level < nlevels_i; ++level) {
          for (size_t r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
        }
      }
      lookup_df.push_back(std::move(col), var);
    } else { // string
      const std::string* base = str_flat.data() + info.offset;
      std::vector<std::string> col(lookup_nrows);
      size_t idxw = 0;
      for (size_t t = 0; t < times; ++t) {
        for (size_t level = 0; level < nlevels_i; ++level)
          for (size_t r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    }
  }

  result.push_back(std::move(nlevels), "nlevels");
  result.push_back(std::move(indices), "indices");
  result.push_back(std::move(lookups_per_variable), "lookups_per_variable");
  result.push_back(std::move(combined_index), "index");
  result.push_back(std::move(lookup_df), "lookup");
  return result;
}


// --------------------------- Linear algebra helpers (FlatMatrix-backed) ----
// cholesky2: in-place working on FlatMatrix (n x n), returns rank * nonneg
int cholesky2(FlatMatrix& matrix, size_t n, double toler) {
  double* base = matrix.data_ptr();
  double eps = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double val = matrix(i, i);
    if (val > eps) eps = val;
  }
  if (eps == 0.0) eps = toler; else eps *= toler;
  int nonneg = 1;
  int rank = 0;

  for (size_t i = 0; i < n; ++i) {
    double* col_i = base + i * n;
    double pivot = col_i[i];
    if (std::isinf(pivot) || pivot < eps) {
      col_i[i] = 0.0;
      if (pivot < -8.0 * eps) nonneg = -1;
    } else {
      ++rank;
      for (size_t j = i + 1; j < n; ++j) {
        double* col_j = base + j * n;
        double temp = col_i[j] / pivot;
        col_i[j] = temp;
        col_j[j] -= temp * temp * pivot;
        for (size_t k = j + 1; k < n; ++k) {
          col_j[k] -= temp * col_i[k];
        }
      }
    }
  }
  return rank * nonneg;
}


// chsolve2 assumes matrix holds the representation produced by cholesky2
void chsolve2(FlatMatrix& matrix, size_t n, double* y) {
  // Forward substitution L * z = y
  double* base = matrix.data_ptr();
  for (size_t j = 0; j < n-1; ++j) {
    double yj = y[j];
    if (yj == 0.0) continue;
    double* col_j = base + j * n;
    for (size_t i = j + 1; i < n; ++i) {
      y[i] -= yj * col_j[i];
    }
  }
  // Now y holds z; solve L^T * x = z
  if (n == 0) return;
  for (int i = n; i-- > 0; ) {
    double* col_i = base + i * n;
    double diag = col_i[i];
    if (diag == 0.0) {
      y[i] = 0.0;
    } else {
      double temp = y[i] / diag;
      for (size_t j = i + 1; j < n; ++j) temp -= y[j] * col_i[j];
      y[i] = temp;
    }
  }
}


// invsympd: returns the inverse of a symmetric positive definite matrix
FlatMatrix invsympd(const FlatMatrix& matrix, size_t n, double toler) {
  FlatMatrix v = matrix; // copy
  cholesky2(v, n, toler);
  FlatMatrix iv(n, n);
  for (size_t i = 0; i < n; ++i) {
    iv(i,i) = 1.0;
    double* ycol = iv.data_ptr() + i * n;
    chsolve2(v, n, ycol);
  }
  return iv;
}


// invchol assumes matrix holds the representation produced by cholesky2
FlatMatrix invchol(FlatMatrix& matrix, size_t n) {
  // Forward substitution L * z = y
  double* base = matrix.data_ptr();
  FlatMatrix iv(n, n);
  for (size_t k = 0; k < n; ++k) {
    iv(k,k) = 1.0;
    double* y = iv.data_ptr() + k * n;
    for (size_t j = 0; j < n-1; ++j) {
      double yj = y[j];
      if (yj == 0.0) continue;
      double* col_j = base + j * n;
      for (size_t i = j + 1; i < n; ++i) {
        y[i] -= yj * col_j[i];
      }
    }
  }
  return iv;
}


std::vector<double> mat_vec_mult(const FlatMatrix& A, const std::vector<double>& x) {
  size_t m = A.nrow;
  size_t p = A.ncol;
  if (x.size() != p)
    throw std::invalid_argument("Vector size mismatch");
  std::vector<double> result(m, 0.0);
  for (size_t c = 0; c < p; ++c) {
    double xc = x[c];
    size_t offset = c * m;
    const double* colptr = A.data_ptr() + offset;
    for (size_t r = 0; r < m; ++r) result[r] += colptr[r] * xc;
  }
  return result;
}

FlatMatrix mat_mat_mult(const FlatMatrix& A, const FlatMatrix& B) {
  size_t m = A.nrow;
  size_t k = A.ncol;
  size_t k2 = B.nrow;
  size_t n = B.ncol;
  if (k != k2) throw std::invalid_argument("Matrix dimensions mismatch");
  if (m == 0 || k == 0 || n == 0) return FlatMatrix();
  FlatMatrix C(m, n);
  // Column-major: For each column j in B/C, compute
  // C[:,j] = sum_{t=0..k-1} A[:,t] * B[t,j]
  for (size_t j = 0; j < n; ++j) {
    const double* bcol = B.data_ptr() + j * k;
    double* ccol = C.data_ptr() + j * m;
    for (size_t t = 0; t < k; ++t) {
      const double* acol = A.data_ptr() + t * m;
      double scale = bcol[t];
      if (scale == 0.0) continue;
      for (size_t i = 0; i < m; ++i) {
        ccol[i] += acol[i] * scale;
      }
    }
  }
  return C;
}

// Transpose a FlatMatrix (double)
FlatMatrix transpose(const FlatMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return FlatMatrix();

  const size_t src_nrow = M.nrow;
  const size_t src_ncol = M.ncol;
  FlatMatrix out(src_ncol, src_nrow); // swapped dims

  const double* src = M.data_ptr();
  double* dst = out.data_ptr();

  for (size_t c = 0; c < src_ncol; ++c) {
    const double* src_col = src + c * src_nrow;
    for (size_t r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }

  return out;
}

// Transpose an IntMatrix (int)
IntMatrix transpose(const IntMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return IntMatrix();

  const size_t src_nrow = M.nrow;
  const size_t src_ncol = M.ncol;
  IntMatrix out(src_ncol, src_nrow); // swapped dims

  const int* src = M.data_ptr();
  int* dst = out.data_ptr();

  for (size_t c = 0; c < src_ncol; ++c) {
    const int* src_col = src + c * src_nrow;
    for (size_t r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }

  return out;
}

// Transpose a BoolMatrix (unsigned char)
BoolMatrix transpose(const BoolMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return BoolMatrix();

  const size_t src_nrow = M.nrow;
  const size_t src_ncol = M.ncol;
  BoolMatrix out(src_ncol, src_nrow); // swapped dims

  const unsigned char* src = M.data_ptr();
  unsigned char* dst = out.data_ptr();

  for (size_t c = 0; c < src_ncol; ++c) {
    const unsigned char* src_col = src + c * src_nrow;
    for (size_t r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }

  return out;
}

double quadsym(const std::vector<double>& u, const FlatMatrix& v) {
  size_t p = u.size();
  const double* vptr = v.data_ptr();
  const double* uptr = u.data();
  double sum = 0.0;

  for (size_t j = 0; j < p; ++j) {
    const double* col = vptr + j * p;
    // diagonal term
    sum += uptr[j] * uptr[j] * col[j]; // col[j] == v(j,j)
    // off-diagonals i < j. Access column j contiguous for i = 0..j-1
    double s = 0.0;
    for (size_t i = 0; i < j; ++i) s += col[i] * uptr[i];
    sum += 2.0 * uptr[j] * s; // account for symmetric pair (i,j) and (j,i)
  }
  return sum;
}


// Householder vector
// Given an n-vector x, this function computes an n-vector v with v(1) = 1
// such that (I - 2*v*t(v)/t(v)*v)*x is zero in all but the first component.
std::vector<double> house(const std::vector<double>& x) {
  size_t n = x.size();
  double sumxx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  double mu = std::sqrt(sumxx);
  std::vector<double> v = x;
  if (mu > 0.0) {
    double beta = x[0] + std::copysign(1.0, x[0])*mu;
    for (size_t i=1; i<n; ++i) {
      v[i] /= beta;
    }
  }
  v[0] = 1.0;
  return v;
}

// Householder pre-multiplication
// Given an m-by-n matrix A and a nonzero m-vector v with v(1) = 1,
// the following algorithm overwrites A with P*A where
// P = I - 2*v*t(v)/t(v)*v.
void row_house(FlatMatrix& A, const size_t i1, const size_t i2, const size_t j1,
               const size_t j2, const std::vector<double>& v) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow) {
    throw std::invalid_argument("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol) {
    throw std::invalid_argument("Invalid column indices j1 and j2");
  }

  size_t m = i2-i1+1, n = j2-j1+1;
  double sumvv = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double beta = -2.0 / sumvv;
  std::vector<double> w(n);
  for (size_t j=0; j<n; ++j) {
    for (size_t i=0; i<m; ++i) {
      w[j] += A(i+i1,j+j1)*v[i];
    }
    w[j] *= beta;
  }

  for (size_t j=0; j<n; ++j) {
    for (size_t i=0; i<m; ++i) {
      A(i+i1,j+j1) += v[i]*w[j];
    }
  }
}

// Householder post-multiplication
// Given an m-by-n matrix A and a nonzero n-vector v with v(1) = 1,
// the following algorithm overwrites A with A*P where
// P = I - 2*v*t(v)/t(v)*v.
void col_house(FlatMatrix& A, const size_t i1, const size_t i2, const size_t j1,
               const size_t j2, const std::vector<double>& v) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow) {
    throw std::invalid_argument("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol) {
    throw std::invalid_argument("Invalid column indices j1 and j2");
  }

  size_t m = i2-i1+1, n = j2-j1+1;
  double sumvv = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double beta = -2.0 / sumvv;
  std::vector<double> w(m);
  for (size_t j=0; j<n; ++j) {
    for (size_t i=0; i<m; ++i) {
      w[i] += A(i+i1,j+j1)*v[j];
    }
  }
  for (size_t i=0; i<m; ++i) {
    w[i] *= beta;
  }
  for (size_t j=0; j<n; ++j) {
    for (size_t i=0; i<m; ++i) {
      A(i+i1,j+j1) += w[i]*v[j];
    }
  }
}

// Given scalars a and b, this function computes
// c = cos(theta) and s = sin(theta) so that
//               |  c   s  |^T | a |  =  | r |
//               | -s   c  |   | b |     | 0 |
std::vector<double> givens(const double a, const double b) {
  double c, s, tau;

  if (b == 0.0) {
    c = 1.0; s = 0.0;
  } else {
    if (std::fabs(b) > std::fabs(a)) {
      double d = -std::copysign(1.0, b);
      tau = -a/b; s = d*1.0/std::sqrt(1.0 + tau*tau); c = s*tau;
    } else {
      double d = std::copysign(1.0, a);
      tau = -b/a; c = d*1.0/std::sqrt(1.0 + tau*tau); s = c*tau;
    }
  }

  std::vector<double> result = {c,s};
  return result;
}

// Given A in R^(2xq), c = cos(theta), and s = sin(theta),
// the following algorithm overwrites A with the matrix
//               |  c   s  |^T  A
//               | -s   c  |
void row_rot(FlatMatrix& A, const size_t i1, const size_t i2,
             const size_t j1, const size_t j2,
             const double c, const double s) {
  if (i1 < 0 || i1 >= i2 || i2 >= A.nrow) {
    throw std::invalid_argument("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol) {
    throw std::invalid_argument("Invalid column indices j1 and j2");
  }

  size_t q = j2-j1+1;
  for (size_t j=0; j<q; ++j) {
    size_t jj1 = j+j1;
    double tau1 = A(i1,jj1);
    double tau2 = A(i2,jj1);
    A(i1,jj1) = c*tau1 - s*tau2;
    A(i2,jj1) = s*tau1 + c*tau2;
  }
}

// Given A in R^(qx2), c = cos(theta), and s = sin(theta),
// the following algorithm overwrites A with the matrix
//               A  |  c   s  |
//                  | -s   c  |
void col_rot(FlatMatrix& A, const size_t i1, const size_t i2,
             const size_t j1, const size_t j2,
             const double c, const double s) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow) {
    throw std::invalid_argument("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 >= j2 || j2 >= A.ncol) {
    throw std::invalid_argument("Invalid column indices j1 and j2");
  }

  size_t q = i2-i1+1;
  for (size_t i=0; i<q; ++i) {
    size_t ii1 = i+i1;
    double tau1 = A(ii1,j1);
    double tau2 = A(ii1,j2);
    A(ii1,j1) = c*tau1 - s*tau2;
    A(ii1,j2) = s*tau1 + c*tau2;
  }
}

// Householder Bidiagonalization
// Given A in R^(mxn) with m>=n, the following algorithm overwrites
// the upper bidiagonal part of A with the upper bidiagonal part of
// t(U)*A*V = B, where B is upper bidiagonal and U = U_1 ... U_n and
// V = V_1 ... V_{n-2}. The essential part of U_j's Householder vector
// is stored in A((j+1):m, j), while the essential part of V_j's
// Householder vector is stored in A(j, (j+2):n).
ListCpp house_bidiag(FlatMatrix& A, const bool outtransform = true) {
  size_t m = A.nrow, n = A.ncol;
  if (m < n) {
    throw std::invalid_argument("The input matrix must have # rows >= # columns");
  }

  double tol = 1e-12;
  bool bidiag = true;
  for (size_t j=2; j<n; ++j) {
    for (size_t i=0; i<j-1; ++i) {
      if (std::fabs(A(i,j)) > tol) {
        bidiag = false; break;
      }
    }
  }
  for (size_t j=0; j<n-1; ++j) {
    for (size_t i=j+1; i<n; ++i) {
      if (std::fabs(A(i,j)) > tol) {
        bidiag = false; break;
      }
    }
  }
  for (size_t j=0; j<n; ++j) {
    for (size_t i=n; i<m-1; ++i) {
      if (std::fabs(A(i,j)) > tol) {
        bidiag = false; break;
      }
    }
  }

  FlatMatrix B(n,n);
  FlatMatrix U(m,m), V(n,n);
  for (size_t i=0; i<m; ++i) U(i,i) = 1.0;
  for (size_t i=0; i<n; ++i) V(i,i) = 1.0;

  if (bidiag) {
    B = A;
  } else {
    for (size_t j=0; j<n; ++j) {
      std::vector<double> v(m-j);
      for (size_t i=0; i<m-j; ++i) {
        v[i] = A(i+j,j);
      }
      v = house(v);

      row_house(A, j, m-1, j, n-1, v);

      // update the sub-diagonal elements of column j
      for (size_t i=1; i<m-j; ++i) {
        A(i+j,j) = v[i];
      }

      if (j < n-2) {
        std::vector<double> v(n-j-1);
        for (size_t i=0; i<n-j-1; ++i) {
          v[i] = A(j,i+j+1);
        }
        v = house(v);

        col_house(A, j, m-1, j+1, n-1, v);

        // update the elements of row j
        for (size_t i=1; i<n-j-1; ++i) {
          A(j,i+j+1) = v[i];
        }
      }
    }

    if (outtransform) {
      for (size_t j = n; j-- > 0; ) {
        std::vector<double> v(m-j);
        v[0] = 1.0;
        for (size_t i=1; i<m-j; ++i) {
          v[i] = A(i+j,j);
        }

        row_house(U, j, m-1, j, m-1, v);
      }

      for (size_t j = n-2; j-- > 0; ) {
        std::vector<double> v(n-j-1);
        v[0] = 1.0;
        for (size_t i = 1; i < n-j-1; ++i) {
          v[i] = A(j,i+j+1);
        }

        row_house(V, j+1, n-1, j+1, n-1, v);
      }
    }

    for (size_t j=0; j<n; ++j) {
      B(j,j) = A(j,j);
      if (j<n-1) {
        B(j,j+1) = A(j,j+1);
      }
    }
  }

  ListCpp result;
  if (outtransform) {
    result.push_back(std::move(B), "B");
    result.push_back(std::move(U), "U");
    result.push_back(std::move(V), "V");
  } else {
    result.push_back(std::move(B), "B");
  }
  return result;
}

// Given a bidiagonal matrix with a zero diagonal, premultiplication
// by a sequence of Givens transformations to zero the entire row
ListCpp zero_diagonal(FlatMatrix& B, const size_t k, const bool outtransform = true) {
  size_t n = B.nrow;
  if (B.ncol != n) {
    throw std::invalid_argument("The input matrix must be a square matrix");
  }
  if (k < 0 || k >= n-1) {
    throw std::invalid_argument("Invalid value for index k");
  }
  FlatMatrix U(n,n);
  for (size_t i=0; i<n; ++i) U(i,i) = 1.0;

  for (size_t j=k+1; j<n; ++j) {
    std::vector<double> v = givens(B(k,j), B(j,j));
    double w = v[0];
    v[0] = -v[1]; v[1] = w;
    size_t j1 = j < n-1 ? j+1 : n-1;
    row_rot(B, k, j, j, j1, v[0], v[1]);
    if (outtransform) col_rot(U, k, j, k, j, v[0], v[1]);
  }

  ListCpp result;
  if (outtransform) {
    result.push_back(B, "B");
    result.push_back(std::move(U), "U");
  } else {
    result.push_back(B, "B");
  }
  return result;
}

// Golub-Kahan SVD Step
// Given a bidiagonal matrix B having no zeros on its diagonal or
// superdiagonal, the following algorithm overwrites B with the
// bidiagonal matrix t(U)*B*V, where U and V are orthogonal and V
// is essentially the orthogonal matrix that would be obtained by
// applying Algorithm 8.2.2 in Golub and Van Loan (1989) to T = t(B)*B.
ListCpp svd_step(FlatMatrix& B, const bool outtransform = true) {
  size_t n = B.ncol;
  FlatMatrix U(n,n), V(n,n);
  for (size_t i=0; i<n; ++i) {
    U(i,i) = 1.0; V(i,i) = 1.0;
  }

  if (n == 1) {
    ListCpp result;
    if (outtransform) {
      result.push_back(B, "B");
      result.push_back(U, "U");
      result.push_back(V, "V");
    } else {
      result.push_back(B, "B");
    }
    return result;
  }

  double f1 = (n > 2) ? B(n-3,n-2) : 0.0;
  double f2 = B(n-2,n-1);
  double d1 = B(n-2,n-2), d2 = B(n-1,n-1);
  double a1 = f1*f1 + d1*d1, a2 = f2*f2 + d2*d2, b1 = f2*d1;
  double d = 0.5*(a1-a2);
  double mu = a2 + d - std::copysign(1.0, d)*std::sqrt(d*d + b1*b1);
  double y = B(0,0)*B(0,0) - mu;
  double z = B(0,0)*B(0,1);
  std::vector<double> v(2);
  for (size_t k=0; k<n-1; ++k) {
    v = givens(y,z);
    size_t k1 = k > 0 ? k-1 : 0;
    col_rot(B, k1, k+1, k, k+1, v[0], v[1]);
    if (outtransform) col_rot(V, 0, k+1, k, k+1, v[0], v[1]);

    y = B(k,k);
    z = B(k+1,k);
    v = givens(y,z);
    size_t k2 = k < n-2 ? k+2 : n-1;
    row_rot(B, k, k+1, k, k2, v[0], v[1]);
    if (outtransform) col_rot(U, 0, k+1, k, k+1, v[0], v[1]);

    if (k < n-2) {
      y = B(k,k+1);
      z = B(k,k+2);
    }
  }

  ListCpp result;
  if (outtransform) {
    result.push_back(B, "B");
    result.push_back(std::move(U), "U");
    result.push_back(std::move(V), "V");
  } else {
    result.push_back(B, "B");
  }
  return result;
}

ListCpp svdcpp1(const FlatMatrix& X, const bool outtransform,
                const bool decreasing) {
  size_t m1 = X.nrow, n1 = X.ncol, m, n;
  if (m1 >= n1) {
    m = m1; n = n1;
  } else {
    m = n1; n = m1;
  }

  FlatMatrix Y(m,n);
  if (m1 >= n1) {
    Y = X;
  } else {
    Y = transpose(X);
  }

  ListCpp a = house_bidiag(Y, outtransform);
  FlatMatrix B = a.get<FlatMatrix>("B");

  FlatMatrix U(m,m), V(n,n);
  for (size_t i=0; i<m; ++i) U(i,i) = 1.0;
  for (size_t i=0; i<n; ++i) V(i,i) = 1.0;
  if (outtransform) {
    U = a.get<FlatMatrix>("U");
    V = a.get<FlatMatrix>("V");
  }

  double tol = 1e-12;
  size_t p, q = 0;
  while (q < n) {
    for (size_t i=1; i<n; ++i) {
      if (std::fabs(B(i-1,i)) <= tol*(std::fabs(B(i-1,i-1)) + std::fabs(B(i,i)))) {
        B(i-1,i) = 0.0;
      }
    }

    // find the largest non-negative q and the smallest non-negative p
    // such that
    //               |  B11   0    0   |   p
    //           B = |   0   B22   0   |   n-p-q
    //               |   0    0   B33  |   q
    //                   p  n-p-q  q
    // where B33 is diagonal and B22 has nonzero superdiagonal

    q = n;

    // Find last nonzero superdiagonal entry B(i-1, i) scanning i = n-1 ... 1
    if (n >= 2) {
      for (size_t i = n; i-- > 1; ) {          // i: n-1, ..., 1
        if (B(i-1, i) != 0.0) {
          q = n - i - 1;
          break;
        }
      }
    }

    p = 0;

    // We search for a zero superdiagonal in the active part.
    const size_t active = n - q;
    if (active >= 2) {
      for (size_t i = active - 1; i-- > 1; ) { // i: active-2, ..., 1
        if (B(i-1, i) == 0.0) {
          p = i;
          break;
        }
      }
    }


    if (q < n) {
      // if any diagonal entry in B22 is zero, then zero the superdiagonal
      // entry in the same row
      FlatMatrix B22 = subset_flatmatrix(B, p, n-q, p, n-q);
      for (size_t i=0; i<n-p-q-1; ++i) {
        if (std::fabs(B22(i,i)) < tol) {
          ListCpp b = zero_diagonal(B22, i, outtransform);
          if (outtransform) {
            FlatMatrix Z = b.get<FlatMatrix>("U");
            FlatMatrix W = subset_flatmatrix(U, 0, m, p, n-q);
            for (size_t k=0; k<n-p-q; ++k) {
              for (size_t j=0; j<m; ++j) {
                U(j,k+p) = 0.0;
              }
            }
            for (size_t k=0; k<n-p-q; ++k) {
              for (size_t l=0; l<n-p-q; ++l) {
                for (size_t j=0; j<m; ++j) {
                  U(j,k+p) += W(j,l)*Z(l,k);
                }
              }
            }
          }
        }
      }

      // apply Algorithm 8.3.1 to B22
      ListCpp c = svd_step(B22, outtransform);

      // update B22
      for (size_t j=0; j<n-p-q; ++j) {
        for (size_t i=0; i<n-p-q; ++i) {
          B(i+p,j+p) = B22(i,j);
        }
      }

      if (outtransform) {
        FlatMatrix Z1 = c.get<FlatMatrix>("U");
        FlatMatrix W1 = subset_flatmatrix(U, 0, m, p, n-q);
        for (size_t j=0; j<n-p-q; ++j) {
          for (size_t i=0; i<m; ++i) {
            U(i,j+p) = 0.0;
          }
        }
        for (size_t j=0; j<n-p-q; ++j) {
          for (size_t k=0; k<n-p-q; ++k) {
            for (size_t i=0; i<m; ++i) {
              U(i,j+p) += W1(i,k)*Z1(k,j);
            }
          }
        }

        FlatMatrix Z2 = c.get<FlatMatrix>("V");
        FlatMatrix W2 = subset_flatmatrix(V, 0, n, p, n-q);
        for (size_t j=0; j<n-p-q; ++j) {
          for (size_t i=0; i<n; ++i) {
            V(i,j+p) = 0.0;
          }
        }
        for (size_t j=0; j<n-p-q; ++j) {
          for (size_t k=0; k<n-p-q; ++k) {
            for (size_t i=0; i<n; ++i) {
              V(i,j+p) += W2(i,k)*Z2(k,j);
            }
          }
        }
      }
    }
  }

  std::vector<double> d(n);
  for (size_t i=0; i<n; ++i) d[i] = B(i,i);

  // ensure the singular values are positive
  for (size_t i=0; i<n; ++i) {
    if (d[i] < 0.0) {
      d[i] = -d[i];
      for (size_t j=0; j<n; ++j) V(j,i) = -V(j,i);
    }
  }

  if (decreasing) {
    // order the singular values from the largest to the smallest
    // and the arrange the associated vectors accordingly
    std::vector<size_t> order = seqcpp(0, n-1);
    std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
      return d[i] > d[j];
    });

    subset_in_place(d, order);
    if (outtransform) {
      FlatMatrix Z = U;
      FlatMatrix W = V;
      for (size_t i=0; i<n; ++i) {
        size_t k = order[i];
        for (size_t j=0; j<m; ++j) U(j,i) = Z(j,k);
        for (size_t j=0; j<n; ++j) V(j,i) = W(j,k);
      }
    }
  }

  // switch U and V if m1 < n1
  FlatMatrix U1(m1,m1), V1(n1,n1);
  if (m1 >= n1) {
    U1 = U; V1 = V;
  } else {
    U1 = V; V1 = U;
  }

  ListCpp result;
  if (outtransform) {
    result.push_back(std::move(d), "d");
    result.push_back(std::move(U1), "U");
    result.push_back(std::move(V1), "V");
  } else {
    result.push_back(std::move(d), "d");
  }
  return result;
}

//' @title Singular Value Decomposition of a Matrix
//' @description Computes the singular-value decomposition of a
//' rectangular matrix.
//'
//' @param X A numeric matrix whose SVD decomposition is to be computed.
//' @param outtransform Whether the orthogonal matrices composing of the
//'   left and right singular vectors are to be computed.
//' @param decreasing Whether the singular values should be sorted in
//'   decreasing order and the corresponding singular vectors rearranged
//'   accordingly.
//'
//' @details
//' Given \eqn{A \in R^{m\times n} (m \geq n)}, the following algorithm
//' overwrites \eqn{A} with \eqn{U^T A V = D}, where
//' \eqn{U\in R^{m\times m}} is orthogonal, \eqn{V \in R^{n\times n}} is
//' orthogonal, and \eqn{D \in R^{m\times n}} is diagonal.
//'
//' @return A list with the following components:
//'
//' * \code{d}: A vector containing the singular values of \eqn{X}.
//'
//' * \code{U}: A matrix whose columns contain the left singular vectors
//'   of \eqn{X}.
//'
//' * \code{V}: A matrix whose columns contain the right singular vectors
//'   of \eqn{X}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Gene N. Golub and Charles F. Van Loan.
//' Matrix Computations, second edition. Baltimore, Maryland:
//' The John Hopkins University Press, 1989, p.434.
//'
//' @examples
//'
//' A <- matrix(c(1,0,0,0, 1,2,0,0, 0,1,3,0, 0,0,1,4), 4, 4)
//' svdcpp(A)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List svdcpp(const Rcpp::NumericMatrix& X,
                  const bool outtransform = true,
                  const bool decreasing = true) {
  auto fm = flatmatrix_from_Rmatrix(X);
  auto cpp_result = svdcpp1(fm, outtransform, decreasing);
  return Rcpp::wrap(cpp_result);
}


//' @title Converting a decimal to a fraction
//' @description Converts a decimal to a fraction based on the algorithm
//' from http://stackoverflow.com/a/5128558/221955.
//'
//' @param x The fraction in decimal form.
//' @param tol The tolerance level for the conversion error.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' float_to_fraction(5/3)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> float_to_fraction(const double x,
                                      const double tol = 0.000001) {
  std::vector<double> v(2);
  double x1 = x;
  double n = std::floor(x1);
  x1 = x1 - n;
  if (x1 < tol) {
    v[0] = n;
    v[1] = 1;
  } else if (1 - tol < x1) {
    v[0] = n+1;
    v[1] = 1;
  } else {
    // The lower fraction is 0/1
    double lower_n = 0;
    double lower_d = 1;

    // The upper fraction is 1/1
    double upper_n = 1;
    double upper_d = 1;

    bool cond = true;
    while (cond) {
      // The middle fraction is (lower_n + upper_n) / (lower_d + upper_d)
      double middle_n = lower_n + upper_n;
      double middle_d = lower_d + upper_d;

      // If x + tol < middle
      if (middle_d * (x1 + tol) < middle_n) {
        // middle is our new upper
        upper_n = middle_n;
        upper_d = middle_d;
      } else if (middle_n < (x1 - tol) * middle_d) {
        // Else If middle < x - tol
        // middle is our new lower
        lower_n = middle_n;
        lower_d = middle_d;
      } else {
        // Else middle is our best fraction
        v[0] = n * middle_d + middle_n;
        v[1] = middle_d;
        cond = false;
      }
    }
  }

  return v;
}

