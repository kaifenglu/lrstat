#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <Rcpp.h>

#include <algorithm>  // fill, sort, lower_bound, upper_bound, max_element,
// swap, for_each, none_of, all_of
#include <cmath>      // fabs, isnan, isinf, exp, log, pow, sqrt, copysign
#include <cstddef>    // size_t
#include <cstring>    // memcpy
#include <functional> // function
#include <limits>     // numeric_limits
#include <memory>     // make_shared, shared_ptr
#include <numeric>    // accumulate, inner_product, iota
#include <queue>      // priority_queue
#include <sstream>    // ostringstream
#include <stdexcept>  // invalid_argument, runtime_error, out_of_range
#include <string>     // string, tolower
#include <utility>    // pair, swap
#include <vector>     // vector

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/tools/minima.hpp>

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


std::vector<double> stl_sort(const std::vector<double>& x) {
  std::vector<double> y = x;
  std::sort(y.begin(), y.end());
  return y;
}

std::vector<int> seqcpp(int start, int end) {
  if (start > end) throw std::invalid_argument(
      "start must be less than or equal to end for the sequence function.");
  int size = end - start + 1;
  std::vector<int> result(size);
  std::iota(result.begin(), result.end(), start);
  return result;
}

std::vector<unsigned char> convertLogicalVector(const Rcpp::LogicalVector& vec) {
  int n = vec.size();
  std::vector<unsigned char> result;
  result.resize(n);
  for (int i = 0; i < n; ++i) {
    int v = vec[i];
    if (v == NA_LOGICAL) result[i] = 255; // NA representation
    else result[i] = v ? 1 : 0; // TRUE -> 1, FALSE -> 0
  }
  return result;
}

std::vector<int> which(const std::vector<unsigned char>& vec) {
  std::vector<int> indices;
  indices.reserve(vec.size());
  int n = vec.size();
  for (int i = 0; i < n; ++i) {
    if (vec[i] != 0 && vec[i] != 255) indices.push_back(i);
  }
  return indices;
}

std::vector<int> findInterval3(const std::vector<double>& x,
                               const std::vector<double>& v,
                               bool rightmost_closed,
                               bool all_inside,
                               bool left_open) {
  std::vector<int> out(x.size());
  const double* v_begin = v.data();
  const double* v_end   = v_begin + v.size();
  const int n = x.size();
  const int nv = v.size();

  for (int i = 0; i < n; ++i) {
    double xi = x[i];
    if (std::isnan(xi)) { out[i] = -1; continue; }
    const double* pos = left_open ? std::lower_bound(v_begin, v_end, xi) :
      std::upper_bound(v_begin, v_end, xi);
    int idx = static_cast<int>(pos - v_begin);
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

// --------------------------- Root finders -----------------------------------
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol, int maxiter) {
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

double bisect(const std::function<double(double)>& f,
              double x1, double x2, double tol, int maxiter) {
  double a = x1, b = x2;
  double fa = f(a), fb = f(b);
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    throw std::invalid_argument("Root must be bracketed in bisect");
  if (std::fabs(fa) < tol) return a;
  if (std::fabs(fb) < tol) return b;
  double xmid, fmid;
  for (int j = 1; j <= maxiter; ++j) {
    xmid = a + 0.5 * (b - a);
    fmid = f(xmid);
    if (std::fabs(fmid) < tol || (b - a) < tol) return xmid;
    if ((fa > 0.0 && fmid < 0.0) || (fa < 0.0 && fmid > 0.0)) {
      b = xmid; fb = fmid; }
    else { a = xmid; fa = fmid; }
  }
  throw std::runtime_error("Maximum number of iterations exceeded in bisect");
}


// [[Rcpp::export]]
double errorSpentcpp(const double t,
                     const double error,
                     const std::string& sf,
                     const double sfpar) {
  if (error <= 0 || error >= 1) {
    throw std::invalid_argument("error must be a number between 0 and 1");
  }
  if (t <= 0 || t > 1) {
    throw std::invalid_argument("t must be a number between 0 and 1");
  }

  std::string asf = sf;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  double aval;
  if (asf == "sfp") {
    aval = error * std::log(1.0 + (std::exp(1.0) - 1.0) * t);
  } else if (asf == "sfof") {
    aval = boost_qnorm(1.0 - error / 2.0);
    aval = 2.0 * (1.0 - boost_pnorm(aval / std::sqrt(t)));
  } else if (asf == "sfkd") {
    if (std::isnan(sfpar)) {
      throw std::invalid_argument("sfpar is missing for sfKD");
    } else if (sfpar <= 0) {
      throw std::invalid_argument("sfpar must be positive for sfKD");
    } else {
      aval = error * std::pow(t, sfpar);
    }
  } else if (asf == "sfhsd") {
    if (std::isnan(sfpar)) {
      throw std::invalid_argument("sfpar is missing for sfHSD");
    } else if (sfpar == 0) {
      aval = error * t;
    } else {
      aval = error * (1.0 - std::exp(-sfpar * t)) / (1.0 - std::exp(-sfpar));
    }
  } else {
    throw std::invalid_argument("Invalid spending function");
  }
  return aval;
}


// [[Rcpp::export]]
ListCpp exitprobcpp(const std::vector<double>& b,
                    const std::vector<double>& a,
                    const std::vector<double>& theta,
                    const std::vector<double>& I) {

  // kMax is the total number of stages
  int kMax = static_cast<int>(b.size());

  // Integer value controlling grid for numerical integration as in
  // Jennison and Turnbull (2000)
  const int r = 18;
  const int r1 = 6 * r - 1;   // size of x1 and shift
  const int r2 = 12 * r - 3;  // size of z, h vectors when fully expanded

  // Prepare a1, theta1, I1 only if defaults / expansion necessary.
  std::vector<double> a1; a1.reserve(kMax);
  if (none_na(a)) {
    if (static_cast<int>(a.size()) < kMax)
      throw std::invalid_argument("Insufficient length for a");
    a1 = a; // copy once
  } else {
    a1.assign(kMax, -6.0);
    a1[kMax - 1] = b[kMax - 1]; // set last element to b[kMax-1]
  }

  // check lower < upper
  for (int i = 0; i < kMax; ++i) {
    if (a1[i] > b[i]) throw std::invalid_argument(
        "Lower bounds (a) must be less than upper bounds (b)");
  }

  // theta expansion
  std::vector<double> theta1;
  if (none_na(theta)) {
    if (static_cast<int>(theta.size()) == 1) {
      theta1.assign(kMax, theta[0]);
    } else if (static_cast<int>(theta.size()) >= kMax) {
      theta1 = theta;
    } else {
      throw std::invalid_argument("Insufficient length for theta");
    }
  } else {
    theta1.assign(kMax, 0.0);
  }


  // information times expansion / validation
  std::vector<double> I1;
  if (none_na(I)) {
    if (static_cast<int>(I.size()) < kMax)
      throw std::invalid_argument("Insufficient length for I");
    if (I[0] <= 0.0)
      throw std::invalid_argument("I must be positive");
    if (any_nonincreasing(I))
      throw std::invalid_argument("I must be increasing");
    I1 = I;
  } else {
    I1.resize(kMax);
    std::iota(I1.begin(), I1.end(), 1.0);
  }

  // Precompute shifts (constant across stages)
  std::vector<double> shift(r1);
  for (int i = 0; i < r1; ++i) {
    if (i < r - 1) {
      shift[i] = -3.0 - 4.0 * std::log(static_cast<double>(r) / (i + 1.0));
    } else if (i < 5 * r) {
      shift[i] = -3.0 + 3.0 * ( (i + 1.0 - r) / (2.0 * r) );
    } else {
      shift[i] = 3.0 + 4.0 * std::log(static_cast<double>(r) / (6.0 * r - i - 1.0));
    }
  }


  // Precompute sqrt and theta*sqrt/I combos
  std::vector<double> sqrtI(kMax), thetaSqrtI(kMax), thetaI(kMax);
  for (int j = 0; j < kMax; ++j) {
    sqrtI[j] = std::sqrt(I1[j]);
    thetaSqrtI[j] = theta1[j] * sqrtI[j];
    thetaI[j] = theta1[j] * I1[j];
  }


  // dI and dThetaI
  std::vector<double> dI(kMax), dThetaI(kMax);
  dI[0] = I1[0];
  dThetaI[0] = thetaI[0];
  for (int j = 1; j < kMax; ++j) {
    dI[j] = I1[j] - I1[j-1];
    dThetaI[j] = thetaI[j] - thetaI[j-1];
  }

  // pre-allocate buffers once, reuse across stages
  std::vector<double> exitProbUpper(kMax), exitProbLower(kMax);

  // allocate working arrays at max sizes once
  std::vector<double> x1(r1), x(r1);  // x1 is untrimmed, x is trimmed
  std::vector<double> w(r2), z(r2), h(r2); // z, h for the current stage
  std::vector<double> z0, h0; // z0, h0 for the previous stage
  z0.reserve(r2); h0.reserve(r2);

  int m0 = 0; // size of previous z0/h0
  // z0 and h0 are represented by the first m0 entries of z and h vectors.

  for (int j = 0; j < kMax; ++j) {
    const double thetaSqrtIj = thetaSqrtI[j];
    const double sqrtIj = sqrtI[j];
    const double sqrtIjm1 = sqrtI[j-1];
    const double a1j = a1[j];
    const double bj = b[j];
    const double dThetaIj = dThetaI[j];
    const double sqrtdIj = std::sqrt(dI[j]);
    const double sqrtI1dIj = std::sqrt(I1[j] / dI[j]);

    // initialize x1 = thetaSqrtI[j] + shift
    for (int i = 0; i < r1; ++i) x1[i] = thetaSqrtIj + shift[i];

    // find trimming indices using binary search (x1 is sorted b/c shift is monotone)
    int i1 = 0;
    if (a1[j] >= x1[0]) {
      // first index i1 s.t. x1[i1] > a1[j]; then use i1-1 as lower trimming index
      auto it = std::upper_bound(x1.begin(), x1.end(), a1[j]);
      i1 = static_cast<int>(it - x1.begin()) - 1;
    }

    int i2 = r1 - 1;
    if (b[j] <= x1[r1 - 1]) {
      // find last index i2 such that x1[i2] < b[j]; then i2+1
      auto it2 = std::lower_bound(x1.begin(), x1.end(), b[j]);
      i2 = static_cast<int>(it2 - x1.begin());
    }


    // m1 is number of retained x nodes after trimming
    const int m1 = i2 - i1 + 1;

    // build x (trimmed) and set x[0]=xlower, x[m1-1]=xupper
    x[0] = std::max(a1[j], x1[0]);
    x[m1 - 1] = std::min(b[j], x1[r1 - 1]);

    // copy interior trimmed x1 values
    for (int i = 1; i < m1 - 1; ++i) x[i] = x1[i + i1];

    // derive z grid (odd + even interleaving)
    const int m = 2 * m1 - 1;
    // odd points
    for (int i = 0; i < m1; ++i) z[2*i] = x[i];
    // even points (midpoints)
    for (int i = 0; i < m1 - 1; ++i) z[2*i + 1] = 0.5 * (z[2*i] + z[2*i + 2]);


    // weights w as Simpson-like composite rule (same formulas)
    // w[0]
    w[0] = (z[2] - z[0]) / 6.0;

    // interior even indices (i = 2, 4, ..., 2*(m1-2))
    for (int i0 = 1; i0 <= m1 - 2; ++i0) {
      int i = 2 * i0;
      w[i] = (z[i + 2] - z[i - 2]) / 6.0;
    }
    // interior odd indices (i = 1,3,...)
    for (int i0 = 1; i0 <= m1 - 1; ++i0) {
      int i = 2 * i0 - 1;
      w[i] = 4.0 * (z[i + 1] - z[i - 1]) / 6.0;
    }
    // last weight
    w[m - 1] = (z[m - 1] - z[m - 3]) / 6.0;


    // first stage is easy
    if (j == 0) {
      // exit probabilities
      exitProbUpper[j] = boost_pnorm(-bj + thetaSqrtIj);
      exitProbLower[j] = boost_pnorm(a1j - thetaSqrtIj);

      // prepare h0, m0, z0 for the next stage
      if (kMax > 1) {
        m0 = m;
        h0.resize(m0);
        z0.resize(m0);
        for (int i = 0; i < m0; ++i) {
          h0[i] = w[i] * boost_dnorm(z[i] - thetaSqrtIj);
          z0[i] = z[i];
        }
      }
    } else {
      // calculate exit probabilities using h0 from the previous stage
      double sumUpper = 0.0, sumLower = 0.0;
      for (int i0 = 0; i0 < m0; ++i0) {
        double tupper = (z0[i0] * sqrtIjm1 - bj * sqrtIj + dThetaIj) / sqrtdIj;
        double tlower = (-z0[i0] * sqrtIjm1 + a1j * sqrtIj - dThetaIj) / sqrtdIj;
        sumUpper += h0[i0] * boost_pnorm(tupper);
        sumLower += h0[i0] * boost_pnorm(tlower);
      }
      exitProbUpper[j] = sumUpper;
      exitProbLower[j] = sumLower;

      // prepare h0, m0, z0 for the next stage
      if (j < kMax-1) {
        for (int i = 0; i < m; ++i) {
          double sum = 0.0;
          for (int i0 = 0; i0 < m0; ++i0) {
            double t = (z[i] * sqrtIj - z0[i0] * sqrtIjm1 - dThetaIj) / sqrtdIj;
            sum += h0[i0] * boost_dnorm(t);
          }
          h[i] = sum * w[i] * sqrtI1dIj; // factors invariant to i0
        }

        m0 = m;
        h0.resize(m0);
        z0.resize(m0);
        for (int i = 0; i < m0; ++i) {
          h0[i] = h[i];
          z0[i] = z[i];
        }
      }
    }
  }

  // return a list of stagewise exit probabilities
  ListCpp exitProb;
  exitProb.push_back(std::move(exitProbUpper), "exitProbUpper");
  exitProb.push_back(std::move(exitProbLower), "exitProbLower");
  return exitProb;
}


double dtpwexpcpp1(const double q,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound,
                   const bool logd) {
  double d;
  if (q <= lowerBound) {
    d = 0.0;
  } else {
    std::vector<double> y = {lowerBound, q};
    std::vector<int> i = findInterval3(y, piecewiseSurvivalTime);
    double v;
    int i0 = i[0] - 1;
    int i1 = i[1] - 1;
    if (i0 == i1) {
      v = lambda[i0] * (q - lowerBound);
    } else {
      v = lambda[i0] * (piecewiseSurvivalTime[i0 + 1] - lowerBound);
      for (int j = i0 + 1; j < i1; ++j) {
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
  int m = static_cast<int>(piecewiseSurvivalTime.size());
  std::vector t0 = {lowerBound};
  int i0 = findInterval3(t0, piecewiseSurvivalTime)[0] - 1;
  int m1 = m - i0;
  std::vector<double> ch(m1, 0.0);
  for (int j = 1; j < m1; ++j) {
    int jj = j + i0 - 1;
    double lb = (j == 1) ? lowerBound : piecewiseSurvivalTime[jj];
    ch[j] = ch[j-1] + lambda[jj] * (piecewiseSurvivalTime[jj+1] - lb);
  }

  std::vector<int> idx = findInterval3(q, piecewiseSurvivalTime);
  int n = static_cast<int>(q.size());
  std::vector<double> d(n);

  for (int h = 0; h < n; ++h) {
    if (q[h] <= lowerBound) {
      d[h] = 0.0;
    } else {
      int i1 = idx[h] - 1;
      int j1 = i1 + i0;
      double lb = (j1 == i0) ? lowerBound : piecewiseSurvivalTime[j1];
      double ch1 = ch[i1] + lambda[i1] * (q[h] - lb);
      d[h] = lambda[i1] * std::exp(-ch1);
    }
  }

  if (logd) {
    for (int h = 0; h < n; ++h) {
      d[h] = std::log(d[h]);
    }
  }

  return d;
}


double ptpwexpcpp1(const double q,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound,
                   const bool lowertail,
                   const bool logp) {

  double p;
  if (q <= lowerBound) {
    p = 0.0;
  } else {
    std::vector<double> y = {lowerBound, q};
    std::vector<int> i = findInterval3(y, piecewiseSurvivalTime);
    double v;
    if (i[0] == i[1]) {
      v = lambda[i[0] - 1] * (q - lowerBound);
    } else {
      v = lambda[i[0] - 1] * (piecewiseSurvivalTime[i[0]] - lowerBound);
      for (int j = i[0]; j < i[1] - 1; ++j) {
        v += lambda[j] * (piecewiseSurvivalTime[j + 1] -
          piecewiseSurvivalTime[j]);
      }
      v += lambda[i[1] - 1] * (q - piecewiseSurvivalTime[i[1] - 1]);
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
  int m = static_cast<int>(piecewiseSurvivalTime.size());
  std::vector t0 = {lowerBound};
  int i0 = findInterval3(t0, piecewiseSurvivalTime)[0] - 1;
  int m1 = m - i0;
  std::vector<double> ch(m1, 0.0);
  for (int j = 1; j < m1; ++j) {
    int jj = j + i0 - 1;
    double lb = (j == 1) ? lowerBound : piecewiseSurvivalTime[jj];
    ch[j] = ch[j-1] + lambda[jj] * (piecewiseSurvivalTime[jj+1] - lb);
  }

  std::vector<int> idx = findInterval3(q, piecewiseSurvivalTime);
  int n = static_cast<int>(q.size());
  std::vector<double> p(n);

  for (int h = 0; h < n; ++h) {
    if (q[h] <= lowerBound) {
      p[h] = 0.0;
    } else {
      int i1 = idx[h] - 1;
      int j1 = i1 + i0;
      double lb = (j1 == i0) ? lowerBound : piecewiseSurvivalTime[j1];
      double ch1 = ch[i1] + lambda[i1] * (q[h] - lb);
      p[h] = 1.0 - std::exp(-ch1);
    }
  }

  if (!lowertail) {
    for (int h = 0; h < n; ++h) {
      p[h] = 1.0 - p[h];
    }
  }

  if (logp) {
    for (int h = 0; h < n; ++h) {
      p[h] = std::log(p[h]);
    }
  }

  return p;
}


double qtpwexpcpp1(const double p,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound,
                   const bool lowertail,
                   const bool logp) {
  int m = static_cast<int>(piecewiseSurvivalTime.size());
  double u = logp ? std::exp(p) : p;
  if (!lowertail) u = 1.0 - u;
  if (u <= 0.0) return lowerBound;
  if (u >= 1.0) return std::numeric_limits<double>::infinity();
  double v1 = -log1p(-u);
  int j = 0;
  while (j < m && piecewiseSurvivalTime[j] <= lowerBound) ++j;
  int j1 = (j == 0) ? 0 : (j - 1);
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

// [[Rcpp::export]]
std::vector<double> qtpwexpcpp(
    const std::vector<double>& p,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const double lowerBound,
    const bool lowertail,
    const bool logp) {
  int n = static_cast<int>(p.size());
  std::vector<double> q(n);
  for (int h = 0; h < n; ++h) {
    q[h] = qtpwexpcpp1(p[h], piecewiseSurvivalTime, lambda,
                       lowerBound, lowertail, logp);
  }

  return q;
}


// mean and variance of a truncated piecewise exponential distribution
// [[Rcpp::export]]
ListCpp mtpwexpcpp(const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound) {
  int m = static_cast<int>(piecewiseSurvivalTime.size());
  if (lambda[m-1] == 0.0) {
    throw std::invalid_argument("The last hazard rate must be positive");
  }

  std::vector t0 = {lowerBound};
  int i0 = findInterval3(t0, piecewiseSurvivalTime)[0] - 1;

  double s1 = 0.0, s2 = 0.0, v = 0.0;
  for (int j = i0; j < m; ++j) {
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


std::vector<double> getBoundcpp(
    const int k,
    const std::vector<double>& informationRates,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const std::vector<unsigned char>& efficacyStopping) {

  if (k <= 0) throw std::invalid_argument("k must be provided and positive");

  // infoRates: if missing create 1/k, 2/k, ..., k/k
  std::vector<double> infoRates(k);
  if (none_na(informationRates)) {
    if (static_cast<int>(informationRates.size()) < k)
      throw std::invalid_argument("Insufficient length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[k-1] > 1.0)
      throw std::invalid_argument("informationRates must not exceed 1");
    infoRates = informationRates; // copy
  } else {
    for (int i = 0; i < k; ++i) {
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(k);
    }
  }

  // spendTime: default to infoRates if missing
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (static_cast<int>(spendingTime.size()) < k)
      throw std::invalid_argument("Insufficient length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[k-1] > 1.0)
      throw std::invalid_argument("spendingTime must not exceed 1");
    spendTime = spendingTime;
  } else {
    spendTime = infoRates;
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (static_cast<int>(efficacyStopping.size()) < k)
      throw std::invalid_argument("Insufficient length for efficacyStopping");
    if (efficacyStopping[k-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(k, 1);
  }

  // asf (alpha spending function) to lower-case
  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  // userAlphaSpending checks when asf == "user"
  if (asf == "user") {
    if (!none_na(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be specified");
    if (static_cast<int>(userAlphaSpending.size()) < k)
      throw std::invalid_argument("Insufficient length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[k-1] > alpha)
      throw std::invalid_argument("userAlphaSpending must not exceed alpha");
  }

  if (asf == "of" || asf == "p" || asf == "wt") {
    bool equal_spacing = true;
    for (int i = 0; i < k; ++i) {
      double expected = static_cast<double>(i+1) / static_cast<double>(k);
      if (std::fabs(infoRates[i] - expected) > 1e-6) {
        equal_spacing = false; break;
      }
    }
    if (!equal_spacing) {
      thread_utils::push_thread_warning(
        "Equal spacing is used for OF, P, and WT boundaries");
    }
  }

  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }

  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  std::vector<double> criticalValues(k);

  if (asf == "none") {
    for (int i = 0; i < k-1; ++i) criticalValues[i] = 6.0;
    criticalValues[k-1] = boost_qnorm(1.0 - alpha);
    return criticalValues;
  }

  if (asf == "of" || asf == "p" || asf == "wt") {
    double Delta;
    if (asf == "of") Delta = 0.0;
    else if (asf == "p") Delta = 0.5;
    else Delta = parameterAlphaSpending; // parameterAlphaSpending holds delta for WT

    // for a given multiplier, compute cumulative upper exit probability - alpha
    std::vector<double> u(k);
    std::vector<double> l(k, -6.0);
    std::vector<double> theta(k, 0.0);
    std::vector<double> I(k), u0(k);
    for (int i = 0; i < k; ++i) {
      I[i] = static_cast<double>(i+1) / static_cast<double>(k);
      u0[i] = std::pow(I[i], Delta - 0.5);
    }

    auto f = [&](double aval)->double {
      for (int i = 0; i < k; ++i) {
        u[i] = aval * u0[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      ListCpp probs = exitprobcpp(u, l, theta, I);
      auto v = probs.get<std::vector<double>>("exitProbUpper");
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    double cwt = brent(f, 0.0, 10.0, 1e-6);
    for (int i = 0; i < k; ++i) {
      criticalValues[i] = cwt * u0[i];
      if (!effStopping[i]) criticalValues[i] = 6.0;
    }
    return criticalValues;
  }

  if (asf == "sfof" || asf == "sfp" || asf == "sfkd" ||
      asf == "sfhsd" || asf == "user") {
    // stage 1
    double cumAlpha;
    if (asf == "user") cumAlpha = userAlphaSpending[0];
    else cumAlpha = errorSpentcpp(spendTime[0], alpha, asf, parameterAlphaSpending);

    if (!effStopping[0]) criticalValues[0] = 6.0;
    else criticalValues[0] = boost_qnorm(1.0 - cumAlpha);

    // Preallocate reusable buffers used by the root-finding lambda
    std::vector<double> u_vec; u_vec.reserve(k);
    std::vector<double> l_vec(k, -6.0);
    std::vector<double> theta_vec(k, 0.0);

    // subsequent stages
    for (int k1 = 1; k1 < k; ++k1) {
      // determine cumulative alpha at this stage
      if (asf == "user") cumAlpha = userAlphaSpending[k1];
      else cumAlpha = errorSpentcpp(spendTime[k1], alpha, asf, parameterAlphaSpending);

      if (!effStopping[k1]) {
        criticalValues[k1] = 6.0;
        continue;
      }

      // Ensure reusable buffers have size k1+1 and capacity >= k
      u_vec.resize(k1 + 1);

      // - copy already computed criticalValues[0..k1-1] into u_vec[0..k1-1]
      // the last entry (u_vec[k1]) will be set by the lambda
      std::memcpy(u_vec.data(), criticalValues.data(), k1 * sizeof(double));

      // Define lambda that only sets the last element of u_vec
      auto f = [&](double aval)->double {
        // set the last element to the current candidate critical value
        u_vec[k1] = aval;

        // exitprobcpp expects exact-sized vectors
        ListCpp probs = exitprobcpp(u_vec, l_vec, theta_vec, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - cumAlpha;
      };


      double f_6 = f(6.0);
      if (f_6 > 0.0) { // no alpha spent at current visit
        criticalValues[k1] = 6.0;
      } else {
        auto f_for_brent = [&](double x)->double {
          if (x == 6.0) return f_6; // avoid recomputation at 6.0
          return f(x);
        };
        criticalValues[k1] = brent(f_for_brent, -5.0, 6.0, 1e-6);
      }
    }

    return criticalValues;
  }

  throw std::invalid_argument("Invalid value for typeAlphaSpending");
}

//' @title Efficacy Boundaries for Group Sequential Design
//' @description Obtains the efficacy stopping boundaries for a group
//' sequential design.
//'
//' @param k Look number for the current analysis.
//' @param informationRates Information rates up to the current look. Must be
//'   increasing and less than or equal to 1.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A vector of length \code{k} for the error spending
//'   time at each analysis. Must be increasing and less than or equal to 1.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @inheritParams param_efficacyStopping
//'
//' @details
//' If \code{typeAlphaSpending} is "OF", "P", or "WT", then the boundaries
//' will be based on equally spaced looks.
//'
//' @return A numeric vector of critical values up to the current look.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' getBound(k = 2, informationRates = c(0.5,1),
//'          alpha = 0.025, typeAlphaSpending = "sfOF")
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector getBound(
    const int k = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  auto result = getBoundcpp(
    k, infoRates, alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlpha, spendTime, effStopping
  );

  return Rcpp::wrap(result);
}


BoundCacheAlpha::BoundCacheAlpha(int k,
                                 const std::vector<double>& infoRates,
                                 const std::string& asf,
                                 double asfpar,
                                 const std::vector<double>& userAlphaSpending,
                                 const std::vector<double>& spendTime,
                                 const std::vector<unsigned char>& effStopping,
                                 std::size_t maxEntries,
                                 int alphaPrecision)
  : k_(k),
    infoRates_(infoRates),
    asf_(asf),
    asfpar_(asfpar),
    userAlphaSpending_(userAlphaSpending),
    spendTime_(spendTime),
    effStopping_(effStopping),
    maxEntries_(maxEntries),
    alphaPrecision_(alphaPrecision) {}

int64_t BoundCacheAlpha::discretize(double alpha) const {
  const double scale = std::pow(10.0, alphaPrecision_);
  return static_cast<int64_t>(std::llround(alpha * scale));
}

std::vector<double> BoundCacheAlpha::get(double alpha) {
  int64_t key = discretize(alpha);
  { std::lock_guard<std::mutex> lg(mu_);
    auto it = map_.find(key);
    if (it != map_.end()) {
      usage_.splice(usage_.begin(), usage_, it->second.lruIt);
      return it->second.value;
    }
  }

  // Compute without holding the lock
  std::vector<double> computed =
    getBoundcpp(k_, infoRates_, alpha, asf_, asfpar_, userAlphaSpending_,
                spendTime_, effStopping_);

  std::lock_guard<std::mutex> lg(mu_);
  auto it2 = map_.find(key);
  if (it2 != map_.end()) {
    usage_.splice(usage_.begin(), usage_, it2->second.lruIt);
    return it2->second.value;
  }

  if (map_.size() >= maxEntries_) {
    int64_t lastkey = usage_.back();
    usage_.pop_back();
    map_.erase(lastkey);
  }
  usage_.push_front(key);
  CacheEntry e;
  e.value = std::move(computed);
  e.lruIt = usage_.begin();
  map_.emplace(key, std::move(e));
  return map_[key].value;
}



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
    const std::vector<double>& w) { // w is the sqrt of variance ratio

  double thetaSqrtI0 = theta[0] * std::sqrt(I[0]);
  std::vector<double> a(kMax, 0.0);
  std::vector<double> u(kMax), l(kMax, 0.0);
  for (int i = 0; i < kMax; ++i) u[i] = b[i] * w[i];

  // reusable buffers for prefixes
  std::vector<double> u1; u1.reserve(kMax);
  std::vector<double> l1; l1.reserve(kMax);

  auto f = [&](double x) -> double {
    // reset futility bounds
    std::fill(a.begin(), a.end(), 0.0);
    double eps = 0.0;

    // first stage
    double cb = errorSpentcpp(st[0], x, bsf, bsfpar); // cumulative beta spent
    if (!futilityStopping[0]) {
      a[0] = -6.0;
    } else {
      eps = boost_pnorm(u[0] - thetaSqrtI0) - cb;
      if (eps < 0.0) return -1.0; // to decrease beta
      a[0] = (boost_qnorm(cb) + thetaSqrtI0) / w[0];
    }

    // subsequent stages
    for (int k = 1; k < kMax; ++k) {
      l[k-1] = a[k-1] * w[k-1];
      cb = errorSpentcpp(st[k], x, bsf, bsfpar);
      if (!futilityStopping[k]) {
        a[k] = -6.0;
      } else {
        u1.resize(k + 1);
        l1.resize(k + 1);
        std::memcpy(u1.data(), u.data(), k * sizeof(double));
        std::memcpy(l1.data(), l.data(), k * sizeof(double));
        u1[k] = 6.0;

        // lambda expression for finding futility bound at stage k
        auto g = [&](double aval) -> double {
          l1[k] = aval * w[k];
          ListCpp probs = exitprobcpp(u1, l1, theta, I);
          auto v = probs.get<std::vector<double>>("exitProbLower");
          double cpl = std::accumulate(v.begin(), v.end(), 0.0);
          return cpl - cb;
        };

        double bk = b[k];
        eps = g(bk);
        double g_minus6 = g(-6.0);
        if (g_minus6 > 0.0) { // no beta spent at current visit
          a[k] = -6.0;
        } else if (eps > 0.0) {
          auto g_for_brent = [&](double x)->double {
            if (x == -6.0) return g_minus6;  // avoid recomputation at 6.0
            if (x == bk) return eps;         // avoid recomputation at b[k]
            return g(x);
          };
          a[k] = brent(g_for_brent, -6.0, bk, 1e-6);
        } else if (k < kMax-1) {
          return -1.0;
        }
      }
    }

    return eps;
  };

  double v1 = f(0.0001), v2 = f(1.0 - alpha);
  double beta = 0.0;
  ListCpp probs;
  if (v1 == -1.0 || (v1 < 0.0 && a[kMax-1] == 0.0)) {
    throw std::invalid_argument("Power must be less than 0.9999");
  } else if (v2 > 0.0) {
    throw std::invalid_argument("Power must be greater than alpha");
  } else {
    beta = brent(f, 0.0001, 1.0 - alpha, 1e-6);
    a[kMax-1] = b[kMax-1];
    l[kMax-1] = a[kMax-1] * w[kMax-1];
    probs = exitprobcpp(u, l, theta, I);
  }

  ListCpp result;
  result.push_back(1.0 - beta, "power");
  result.push_back(std::move(a), "futilityBounds");
  result.push_back(std::move(probs), "probs");
  return result;
}


// Integrate a function f(theta) with respect to a normal density of theta:
//   integrate[f(theta) * dnorm(theta, mu, sigma), a, b) /
//             (pnorm(b, mu, sigma) - pnorm(a, mu, sigma)), {theta, a, b}].
double intnorm(const std::function<double(double)>& f,
               double mu, double sigma, double a, double b) {

  int r = 18, r1 = 6 * r - 1;
  double a1 = (a - mu) / sigma , b1 = (b - mu) / sigma;

  std::vector<double> x1(r1);
  for (int i = 0; i < r1; ++i) {
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
  int i1 = 0;
  if (a1 >= x1[0]) {
    // first index i1 such that x1[i1] > a1[j]; then use i1-1 as lower trimming index
    auto it = std::upper_bound(x1.begin(), x1.end(), a1);
    i1 = static_cast<int>(it - x1.begin()) - 1;
  }

  // trim from above
  int i2 = r1 - 1;
  if (b1 <= x1[r1 - 1]) {
    // find last index i2 such that x1[i2] < b[j]; then i2+1
    auto it2 = std::lower_bound(x1.begin(), x1.end(), b1);
    i2 = static_cast<int>(it2 - x1.begin());
  }

  // save the trimmed portion to x
  const int m1 = i2 - i1 + 1;
  std::vector<double> x(m1);
  x[0] = std::max(a1, x1[0]);
  x[m1 - 1] = std::min(b1, x1[r1 - 1]);
  for (int i = 1; i < m1 - 1; ++i) x[i] = x1[i + i1];

  // derive the grid points for z
  const int m = 2 * m1 - 1;
  std::vector<double> z(m), w(m);
  // odd points
  for (int i = 0; i < m1; ++i) z[2*i] = x[i];
  // even points (midpoints)
  for (int i = 0; i < m1 - 1; ++i) z[2*i + 1] = 0.5 * (z[2*i] + z[2*i + 2]);

  // weights w as Simpson-like composite rule (same formulas)
  // w[0]
  w[0] = (z[2] - z[0]) / 6.0;

  // interior even indices (i = 2, 4, ..., 2*(m1-2))
  for (int i0 = 1; i0 <= m1 - 2; ++i0) {
    int i = 2 * i0;
    w[i] = (z[i + 2] - z[i - 2]) / 6.0;
  }
  // interior odd indices (i = 1,3,...)
  for (int i0 = 1; i0 <= m1 - 1; ++i0) {
    int i = 2 * i0 - 1;
    w[i] = 4.0 * (z[i + 1] - z[i - 1]) / 6.0;
  }
  // last weight
  w[m - 1] = (z[m - 1] - z[m - 3]) / 6.0;

  // integrate
  double aval = 0.0;
  for (int i = 0; i < m; ++i) {
    aval += w[i] * f(mu + sigma * z[i]) * boost_dnorm(z[i]);
  }

  double denom = boost_pnorm(b1) - boost_pnorm(a1);

  return aval / denom;
}


std::pair<double, double> mini(
    const std::function<double(double)>& f, double x1, double x2) {

  // Validate inputs
  if (!(x1 < x2)) throw std::invalid_argument("mini: require x1 < x2");

  // Use Boost's brent_find_minima
  return boost::math::tools::brent_find_minima(
    f, x1, x2, std::numeric_limits<double>::digits10);
}


// Numerical integration of f over [lower, upper] with specified tolerance.
// - tol is absolute tolerance; relative behavior depends on integrator.
// - maxiter is max subdivisions/recursions for the GK integrator.
double quad(const std::function<double(double)>& f,
            double lower, double upper, double tol, unsigned maxiter) {

  // Both finite -> use Gauss-Kronrod (good default for finite intervals)
  if (!std::isinf(lower) && !std::isinf(upper)) {
    // Use 15-point GK rule (common choice)
    boost::math::quadrature::gauss_kronrod<double, 15> integrator;
    return integrator.integrate(f, lower, upper, maxiter, tol);
  }

  // Any endpoint infinite -> use tanh-sinh which handles infinite endpoints well
  // tanh_sinh::integrate accepts infinite endpoints (pass +/- inf)
  boost::math::quadrature::tanh_sinh<double> integrator;
  return integrator.integrate(f, lower, upper, tol);
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
static void eval_region(
    const std::function<double(double,double)>& f,
    double ax, double bx,
    double ay, double by,
    double &I3, double &I5) {
  double cx = 0.5*(ax + bx), dx = 0.5*(bx - ax);
  double cy = 0.5*(ay + by), dy = 0.5*(by - ay);

  I3 = 0.0;
  I5 = 0.0;

  // Store f at Gauss nodes to reuse in Kronrod
  double f3[3][3];

  // 3x3 Gauss
  for(int i=0;i<3;++i){
    double xi = cx + dx * x3[i];
    for(int j=0;j<3;++j){
      double yj = cy + dy * x3[j];
      f3[i][j] = f(xi,yj);
      I3 += w3[i]*w3[j]*f3[i][j];
    }
  }
  I3 *= dx*dy;

  // 5×5 rule
  for(int i=0;i<5;++i){
    double xi = cx + dx * x5[i];
    for(int j=0;j<5;++j){
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
  int maxRegions = 1000000;

  double I3, I5;
  eval_region(f, ax, bx, ay, by, I3, I5);

  double global_value = I5;
  double global_err   = std::fabs(I5 - I3);

  std::priority_queue<Region, std::vector<Region>, RegionCompare> pq;

  pq.push({ax, bx, ay, by, I3, I5, global_err});

  while(global_err > tol && global_err > std::fabs(global_value) * tol &&
        (int)pq.size() < maxRegions)
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


// [[Rcpp::export]]
double pbvnormcpp(std::vector<double> lower, std::vector<double> upper, double rho) {
  if (rho == 0.0) {
    double v1 = boost_pnorm(upper[0]) - boost_pnorm(lower[0]);
    double v2 = boost_pnorm(upper[1]) - boost_pnorm(lower[1]);
    return v1 * v2;
  }

  double a2 = lower[1];
  double b2 = upper[1];
  double s = std::sqrt(1.0 - rho * rho);
  auto f = [&](double x)->double {
    double a = (a2 - rho * x) / s;
    double b = (b2 - rho * x) / s;
    double t1 = boost_dnorm(x);
    double t2 = boost_pnorm(b) - boost_pnorm(a);
    return t1 * t2;
  };
  return quad(f, lower[0], upper[0]);
}


// [[Rcpp::export]]
ListCpp hazard_pdcpp(const std::vector<double>& piecewiseSurvivalTime,
                     const std::vector<double>& hazard_pfs,
                     const std::vector<double>& hazard_os,
                     const double rho_pd_os) {
  int n = static_cast<int>(piecewiseSurvivalTime.size());

  // append additional time points for pfs quantiles
  std::vector<double> p(10);
  for (int i=0; i<9; ++i) p[i] = (i+1.0)/10.0;
  p[9] = 0.95;

  std::vector<double> u0 = qtpwexpcpp(p, piecewiseSurvivalTime, hazard_pfs);
  std::vector<double> u(n+10);
  for (int i=0; i<n-1; ++i) u[i] = piecewiseSurvivalTime[i+1];
  u[n-1] = piecewiseSurvivalTime[n-1] + std::log(2.0)/hazard_pfs[n-1];
  for (int i=0; i<10; ++i) u[n+i] = u0[i];

  // obtain sorted and unique time points
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  int m = static_cast<int>(u.size());

  // shifted time points
  std::vector<double> u1(m);
  u1[0] = 0;
  for (int i=1; i<m; ++i) u1[i] = u[i-1];

  // get corresponding hazards
  std::vector<int> index = findInterval3(u1, piecewiseSurvivalTime);
  for (int i=0; i<m; ++i) index[i] = index[i] - 1;
  std::vector<double> hazard_pfs1 = subset(hazard_pfs, index);
  std::vector<double> hazard_os1 = subset(hazard_os, index);

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
  for (int i=0; i<m; ++i) {
    t = u[i];
    v.push_back(u1[i]);
    haz_pfs.push_back(hazard_pfs1[i]);
    haz_os.push_back(hazard_os1[i]);
    hazard_pd[i] = brent(f, 0.5 * (hazard_pfs1[i] - hazard_os1[i]),
                         2.0 * hazard_pfs1[i], tol);
    hazard.push_back(hazard_pd[i]);
  }

  ListCpp result;
  result.push_back(u1, "piecewiseSurvivalTime");
  result.push_back(hazard_pd, "hazard_pd");
  result.push_back(hazard_os1, "hazard_os");
  result.push_back(rho_pd_os, "rho_pd_os");
  return result;
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
  double m1 = quad(fm_pfs, 0.0, upper, tol);
  double m2 = quad(fm2_pfs, 0.0, upper, tol);

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

  ListCpp mv2 = mtpwexpcpp(piecewiseSurvivalTime, hazard_os);
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

// [[Rcpp::export]]
double corr_pfs_oscpp(const std::vector<double>& piecewiseSurvivalTime,
                      const std::vector<double>& hazard_pfs,
                      const std::vector<double>& hazard_os,
                      const double rho_pd_os) {
  ListCpp a = hazard_pdcpp(piecewiseSurvivalTime, hazard_pfs, hazard_os, rho_pd_os);
  std::vector<double> u = a.get<std::vector<double>>("piecewiseSurvivalTime");
  std::vector<double> hazard_pd1 = a.get<std::vector<double>>("hazard_pd");
  std::vector<double> hazard_os1 = a.get<std::vector<double>>("hazard_os");
  return cor_pfs_os(u, hazard_pd1, hazard_os1, rho_pd_os);
}


// [[Rcpp::export]]
ListCpp hazard_subcpp(const std::vector<double>& piecewiseSurvivalTime,
                      const std::vector<double>& hazard_itt,
                      const std::vector<double>& hazard_pos,
                      const double p_pos) {
  int n = static_cast<int>(piecewiseSurvivalTime.size());

  // append additional time points for pfs quantiles
  std::vector<double> p(10);
  for (int i=0; i<9; ++i) p[i] = (i+1.0)/10.0;
  p[9] = 0.95;

  std::vector<double> u0 = qtpwexpcpp(p, piecewiseSurvivalTime, hazard_itt);
  std::vector<double> u(n+10);
  for (int i=0; i<n-1; ++i) u[i] = piecewiseSurvivalTime[i+1];
  u[n-1] = piecewiseSurvivalTime[n-1] + std::log(2.0)/hazard_itt[n-1];
  for (int i=0; i<10; ++i) u[n+i] = u0[i];

  // obtain sorted and unique time points
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  int m = static_cast<int>(u.size());

  // shifted time points
  std::vector<double> u1(m);
  u1[0] = 0;
  for (int i=1; i<m; ++i) u1[i] = u[i-1];

  // get corresponding hazards
  std::vector<int> index = findInterval3(u1, piecewiseSurvivalTime);
  for (int i=0; i<m; ++i) index[i] = index[i] - 1;
  std::vector<double> hazard_itt1 = subset(hazard_itt, index);
  std::vector<double> hazard_pos1 = subset(hazard_pos, index);

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
  for (int i=0; i<m; ++i) {
    t = u[i];
    v.push_back(u1[i]);
    haz_itt.push_back(hazard_itt1[i]);
    haz_pos.push_back(hazard_pos1[i]);
    hazard_neg[i] = brent(f, 0.5 * (hazard_itt1[i] - p_pos * hazard_pos1[i]),
                          2.0 * hazard_itt1[i] / (1.0 - p_pos), tol);
    hazard.push_back(hazard_neg[i]);
  }

  ListCpp result;
  result.push_back(u1, "piecewiseSurvivalTime");
  result.push_back(hazard_pos1, "hazard_pos");
  result.push_back(hazard_neg, "hazard_neg");
  result.push_back(p_pos, "p_pos");
  return result;
}


//' @title Number of Enrolled Subjects
//' @description Obtains the number of subjects enrolled by given calendar
//' times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of enrolled subjects.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_accrualDuration
//'
//' @return A vector of total number of subjects enrolled by the
//' specified calendar times.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Example 1: Uniform enrollment with 20 patients per month for 12 months.
//'
//' accrual(time = 3, accrualTime = 0, accrualIntensity = 20,
//'         accrualDuration = 12)
//'
//'
//' # Example 2: Piecewise accrual, 10 patients per month for the first
//' # 3 months, and 20 patients per month thereafter. Patient recruitment
//' # ends at 12 months for the study.
//'
//' accrual(time = c(2, 9), accrualTime = c(0, 3),
//'         accrualIntensity = c(10, 20), accrualDuration = 12)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> accrual(const std::vector<double>& time,
                            const std::vector<double>& accrualTime,
                            const std::vector<double>& accrualIntensity,
                            const double accrualDuration) {

  int k = static_cast<int>(time.size());
  std::vector<double> n(k);

  // up to end of enrollment
  std::vector<double> t(k);
  for (int i = 0; i < k; ++i) {
    t[i] = std::max(std::min(time[i], accrualDuration), 0.0);
  }

  // identify the time interval containing t
  std::vector<int> m = findInterval3(t, accrualTime);

  // sum up patients enrolled in each interval up to t
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < m[i]; ++j) {
      if (j < m[i] - 1) {
        n[i] += accrualIntensity[j] * (accrualTime[j + 1] - accrualTime[j]);
      } else {
        n[i] += accrualIntensity[j] * (t[i] - accrualTime[j]);
      }
    }
  }

  return n;
}

//' @title Accrual Duration to Enroll Target Number of Subjects
//' @description Obtains the accrual duration to enroll the target number
//' of subjects.
//'
//' @param nsubjects The vector of target number of subjects.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//'
//' @return A vector of accrual durations.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' getAccrualDurationFromN(nsubjects = c(20, 150), accrualTime = c(0, 3),
//'                         accrualIntensity = c(10, 20))
//'
//' @export
// [[Rcpp::export]]
std::vector<double> getAccrualDurationFromN(
    const std::vector<double>& nsubjects,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity) {
  int I = static_cast<int>(nsubjects.size());
  int J = static_cast<int>(accrualTime.size());
  std::vector<double> t(I), p(J);

  p[0] = 0;
  for (int j = 0; j < J - 1; ++j) {
    p[j+1] = p[j] + accrualIntensity[j] * (accrualTime[j+1] - accrualTime[j]);
  }

  std::vector<int> m = findInterval3(nsubjects, p);

  for (int i = 0; i < I; ++i) {
    int j = m[i] - 1;
    t[i] = accrualTime[j] + (nsubjects[i] - p[j]) / accrualIntensity[j];
  }

  return t;
}

//' @title Probability of Being at Risk
//' @description Obtains the probability of being at risk at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//'   probability of being at risk.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of being at risk at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' patrisk(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'         lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> patrisk(const std::vector<double>& time,
                            const std::vector<double>& piecewiseSurvivalTime,
                            const std::vector<double>& lambda,
                            const std::vector<double>& gamma) {
  std::size_t k = time.size();
  std::size_t J = piecewiseSurvivalTime.size();

  // Validate and replicate lambda and gamma
  std::vector<double> lambdax(J), gammax(J);

  if (lambda.size() == 1) {
    lambdax = std::vector<double>(J, lambda[0]);
  } else if (lambda.size() == J) {
    lambdax = lambda;
  } else {
    throw std::invalid_argument("Invalid length for lambda");
  }

  if (gamma.size() == 1) {
    gammax = std::vector<double>(J, gamma[0]);
  } else if (gamma.size() == J) {
    gammax = gamma;
  } else {
    throw std::invalid_argument("Invalid length for gamma");
  }

  // Cumulative hazards for lambda + gamma
  std::vector<double> lamgam(J);
  for (std::size_t j = 0; j < J; ++j) {
    lamgam[j] = lambdax[j] + gammax[j];
  }

  std::vector<double> cumulativeRisk(k, 0.0);

  // Find intervals containing specified analysis time
  std::vector<int> m = findInterval3(time, piecewiseSurvivalTime);

  // Compute cumulative hazard for each time point
  for (std::size_t i = 0; i < k; ++i) {
    double a = 0.0;  // Hazard accumulator for this time point
    for (int j = 0; j < m[i]; ++j) {
      if (j < m[i] - 1) {
        // Contribution from intervals fully covered by time[i]
        a += lamgam[j] * (piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j]);
      } else {
        // Contribution from the remaining portion of the last interval
        a += lamgam[j] * (time[i] - piecewiseSurvivalTime[j]);
      }
    }
    cumulativeRisk[i] = std::exp(-a);  // Apply exponential decay
  }

  return cumulativeRisk;
}

//' @title Probability of Having an Event
//' @description Obtains the probability of having an event at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//'   probability of having an event.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of having an event at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' pevent(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'        lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> pevent(const std::vector<double>& time,
                           const std::vector<double>& piecewiseSurvivalTime,
                           const std::vector<double>& lambda,
                           const std::vector<double>& gamma) {
  std::size_t k = time.size();
  std::size_t J = piecewiseSurvivalTime.size();

  // Validate and replicate lambda and gamma
  std::vector<double> lambdax(J), gammax(J);

  if (lambda.size() == 1) {
    lambdax = std::vector<double>(J, lambda[0]);
  } else if (lambda.size() == J) {
    lambdax = lambda;
  } else {
    throw std::invalid_argument("Invalid length for lambda");
  }

  if (gamma.size() == 1) {
    gammax = std::vector<double>(J, gamma[0]);
  } else if (gamma.size() == J) {
    gammax = gamma;
  } else {
    throw std::invalid_argument("Invalid length for gamma");
  }

  // Compute lambda + gamma
  std::vector<double> lamgam(J);
  for (std::size_t j = 0; j < J; ++j) {
    lamgam[j] = lambdax[j] + gammax[j];
  }

  // Get risk of patients up to each time interval
  std::vector<double> n = patrisk(piecewiseSurvivalTime, piecewiseSurvivalTime,
                                  lambda, gamma);

  std::vector<int> m = findInterval3(time, piecewiseSurvivalTime);
  std::vector<double> a(k, 0.0);

  // Compute cumulative hazard contributions for each time point
  for (std::size_t i = 0; i < k; ++i) {
    double ai = 0.0;  // Accumulator for this time point
    for (int j = 0; j < m[i]; ++j) {
      double p;
      if (j < m[i] - 1) {
        // Full interval is covered
        p = lambdax[j] / lamgam[j] * (1.0 - std::exp(-lamgam[j] *
          (piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j])));
      } else {
        // Partial interval is covered
        p = lambdax[j] / lamgam[j] * (1.0 - std::exp(-lamgam[j] *
          (time[i] - piecewiseSurvivalTime[j])));
      }
      ai += n[j] * p;  // Add risk-weighted probability
    }
    a[i] = ai;
  }

  return a;
}

//' @title Integrated Event Probability Over an Interval With Constant Hazard
//' @description Obtains the integrated probability of having an event
//' during an interval with constant hazard.
//'
//' @param j The analysis time interval with constant hazard.
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A value for the integrated probability of having an event
//' during an interval with constant hazard for a treatment
//' group with specified piecewise exponential survival and dropout
//' distributions.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' hd(j = 1, t1 = 1, t2 = 3, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double hd(const int j,
          const double t1,
          const double t2,
          const std::vector<double>& piecewiseSurvivalTime,
          const std::vector<double>& lambda,
          const std::vector<double>& gamma) {

  int j1 = j-1;

  // lower bound of time interval j for piecewise exponential distribution
  double t0 = piecewiseSurvivalTime[j1];
  std::vector<double> t0_vec = {t0};

  // probability of being at risk at the start of interval j
  double n0 = patrisk(t0_vec, piecewiseSurvivalTime, lambda, gamma)[0];

  // Compute probability of having an event at the start of interval j
  double d0 = pevent(t0_vec, piecewiseSurvivalTime, lambda, gamma)[0];

  // Prepare lambda and gamma vectors
  std::size_t J = piecewiseSurvivalTime.size();
  std::vector<double> lambdax(J), gammax(J);

  if (lambda.size() == 1) {
    lambdax = std::vector<double>(J, lambda[0]);
  } else if (lambda.size() == J) {
    lambdax = lambda;
  } else {
    throw std::invalid_argument("Invalid length for lambda");
  }

  if (gamma.size() == 1) {
    gammax = std::vector<double>(J, gamma[0]);
  } else if (gamma.size() == J) {
    gammax = gamma;
  } else {
    throw std::invalid_argument("Invalid length for gamma");
  }

  // Compute total hazard (lambda + gamma)
  std::vector<double> lamgam(J);
  for (std::size_t i = 0; i < J; ++i) {
    lamgam[i] = lambdax[i] + gammax[i];
  }

  // Integration for conditional probability over (t1, t2)
  double lamgam_j1 = lamgam[j1];
  double exp1 = std::exp(-lamgam_j1 * (t1 - t0));
  double exp2 = std::exp(-lamgam_j1 * (t2 - t0));
  double q1 = (exp1 - exp2) / lamgam_j1;
  double q = lambdax[j1] / lamgam_j1 * (t2 - t1 - q1);

  // Sum up the integration for already failed and to-be-failed
  return d0 * (t2 - t1) + n0 * q;
}

//' @title Integrated Event Probability Over an Interval
//' @description Obtains the integration of the probability of having an
//' event during an interval. The specified analysis time interval can span
//' more than one analysis time interval with constant hazard.
//'
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A value for the integrated probability of having an event
//' during an interval for a treatment group with specified
//' piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' pd(t1 = 1, t2 = 8, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double pd(const double t1,
          const double t2,
          const std::vector<double>& piecewiseSurvivalTime,
          const std::vector<double>& lambda,
          const std::vector<double>& gamma) {

  // Identify analysis time intervals containing t1 and t2
  std::vector<double> t12 = {t1, t2};
  std::vector<int> j12 = findInterval3(t12, piecewiseSurvivalTime);

  int j1 = std::max(j12[0] - 1, 0);  // Ensure index is not less than 0
  int j2 = std::max(j12[1] - 1, 0);

  double a = 0.0;

  // Sum up the integrated event probabilities across analysis time intervals
  for (int j = j1; j <= j2; ++j) {
    double x = 0.0;
    if (j1 == j2) {
      // Both t1 and t2 are in the same interval
      x = hd(j + 1, t1, t2, piecewiseSurvivalTime, lambda, gamma);
    } else if (j == j1) {
      // First interval
      x = hd(j + 1, t1, piecewiseSurvivalTime[j + 1],
             piecewiseSurvivalTime, lambda, gamma);
    } else if (j == j2) {
      // Last interval
      x = hd(j + 1, piecewiseSurvivalTime[j], t2,
             piecewiseSurvivalTime, lambda, gamma);
    } else {
      // Intermediate intervals
      x = hd(j + 1, piecewiseSurvivalTime[j], piecewiseSurvivalTime[j + 1],
             piecewiseSurvivalTime, lambda, gamma);
    }
    a += x;
  }

  return a;
}

//' @title Number of Patients Enrolled During an Interval and Having an Event
//' by Specified Calendar Times
//' @description Obtains the number of patients who are enrolled during a
//' specified enrollment time interval and have an event by the specified
//' calendar times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of patients having an event.
//' @param u1 Lower bound of the accrual time interval.
//' @param u2 Upper bound of the accrual time interval.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of number of patients who are enrolled during a
//' specified enrollment time interval and have an event by the specified
//' calendar times for a given treatment group had the enrollment being
//' restricted to the treatment group. By definition, we must have
//' \code{time >= u2}.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, 10 patients per month for the first 3 months, and
//' # 20 patients per month thereafter. Piecewise exponential survival with
//' # hazard 0.0533 in the first 6 months, and hazard 0.0309 thereafter,
//' # and 5% dropout by the end of 1 year.
//'
//' ad(time = c(9, 15), u1 = 1, u2 = 8, accrualTime = c(0, 3),
//'    accrualIntensity = c(10, 20), piecewiseSurvivalTime=c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> ad(const std::vector<double>& time,
                       const double u1,
                       const double u2,
                       const std::vector<double>& accrualTime,
                       const std::vector<double>& accrualIntensity,
                       const std::vector<double>& piecewiseSurvivalTime,
                       const std::vector<double>& lambda,
                       const std::vector<double>& gamma) {

  // Identify accrual time intervals containing u1 and u2
  std::vector<double> u12 = {u1, u2};
  std::vector<int> j12 = findInterval3(u12, accrualTime);
  int j1 = std::max(j12[0] - 1, 0);  // 0-based index for j1
  int j2 = std::max(j12[1] - 1, 0);  // 0-based index for j2

  std::size_t k = time.size();
  std::vector<double> a(k, 0.0);  // Initialize the result vector with zeroes

  // Sum up the number of patients with an event across accrual time intervals
  for (std::size_t i = 0; i < k; ++i) {
    double t = time[i];  // Current time
    for (int j = j1; j <= j2; ++j) {
      double x = 0.0;
      // Check intervals
      if (j1 == j2) {
        // Both u1 and u2 are in the same interval
        x = pd(t - u2, t - u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j == j1) {
        // First interval
        x = pd(t - accrualTime[j + 1], t - u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j == j2) {
        // Last interval
        x = pd(t - u2, t - accrualTime[j], piecewiseSurvivalTime, lambda, gamma);
      } else {
        // Intermediate intervals
        x = pd(t - accrualTime[j + 1], t - accrualTime[j],
               piecewiseSurvivalTime, lambda, gamma);
      }
      // Add the contribution from this interval
      a[i] += accrualIntensity[j] * x;
    }
  }

  return a;
}


FlatMatrix natriskcpp(const std::vector<double>& time,
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
                      const double maxFollowupTime) {

  std::size_t k = time.size();

  // truncate the analysis time by the maximum follow-up
  std::vector<double> t(k), u(k);
  for (std::size_t i = 0; i < k; ++i) {
    t[i] = std::min(time[i], maxFollowupTime);
    u[i] = std::min(accrualDuration + minFollowupTime - t[i], accrualDuration);
  }

  // Number of patients enrolled
  std::vector<double> a = accrual(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  FlatMatrix n(k, 2);  // FlatMatrix with column-major storage

  // Compute probabilities for the active and control groups
  std::vector<double> patrisk1 = patrisk(t, piecewiseSurvivalTime, lambda1, gamma1);
  std::vector<double> patrisk2 = patrisk(t, piecewiseSurvivalTime, lambda2, gamma2);

  // Compute values for FlatMatrix directly
  for (std::size_t i = 0; i < k; ++i) {
    n(i, 0) = phi * a[i] * patrisk1[i];  // Patients at risk in active treatment
    n(i, 1) = (1.0 - phi) * a[i] * patrisk2[i];  // Patients at risk in control
  }

  return n;
}

//' @title Number of Subjects at Risk
//' @description Obtains the number of subjects at risk at given analysis
//' times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
//'   of patients at risk.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients at risk at the specified
//' analysis times (row) for each treatment group (column).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' natrisk(time = c(9, 24), allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, minFollowupTime = 18,
//'         maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix natrisk(
    const Rcpp::NumericVector& time = NA_REAL,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double minFollowupTime = NA_REAL,
    const double maxFollowupTime = NA_REAL) {

  auto timecpp = Rcpp::as<std::vector<double>>(time);
  auto accrualTimecpp = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualIntensitycpp = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTimecpp = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lambda1cpp = Rcpp::as<std::vector<double>>(lambda1);
  auto lambda2cpp = Rcpp::as<std::vector<double>>(lambda2);
  auto gamma1cpp = Rcpp::as<std::vector<double>>(gamma1);
  auto gamma2cpp = Rcpp::as<std::vector<double>>(gamma2);

  auto cpp_result = natriskcpp(
    timecpp, allocationRatioPlanned, accrualTimecpp, accrualIntensitycpp,
    piecewiseTimecpp, lambda1cpp, lambda2cpp, gamma1cpp, gamma2cpp,
    accrualDuration, minFollowupTime, maxFollowupTime
  );

  return Rcpp::wrap(cpp_result);
}


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
                  const double maxFollowupTime) {

  std::size_t k = time.size();

  // truncate the analysis time by the maximum follow-up
  std::vector<double> t(k), u(k);
  for (std::size_t i = 0; i < k; ++i) {
    t[i] = std::min(time[i], maxFollowupTime);
    u[i] = std::min(accrualDuration + minFollowupTime - t[i], accrualDuration);
  }

  // Number of patients enrolled
  std::vector<double> a = accrual(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Prepare FlatMatrix for results (k rows, 2 columns)
  FlatMatrix d(k, 2);  // FlatMatrix with column-major storage

  // Compute the probabilities of having events in both treatment groups
  std::vector<double> pevent1 = pevent(t, piecewiseSurvivalTime, lambda1, gamma1);
  std::vector<double> pevent2 = pevent(t, piecewiseSurvivalTime, lambda2, gamma2);

  // Constant for ad() calculations
  std::vector<double> u1(1, accrualDuration + minFollowupTime);

  // Compute the number of patients having an event in each group
  for (std::size_t i = 0; i < k; ++i) {
    double ad1 = ad(u1, u[i], accrualDuration, accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, lambda1, gamma1)[0];
    double ad2 = ad(u1, u[i], accrualDuration, accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, lambda2, gamma2)[0];

    double d1 = a[i] * pevent1[i];
    double d2 = a[i] * pevent2[i];

    d(i, 0) = phi * (d1 + ad1);  // Active treatment group
    d(i, 1) = (1.0 - phi) * (d2 + ad2);  // Control group
  }

  return d;
}

//' @title Number of Subjects Having an Event
//' @description Obtains the number of subjects having an event by given
//' analysis times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
//'   of patients having an event.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients having an event at the
//' specified analysis times (row) for each treatment group (column).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' nevent(time = c(9, 24), allocationRatioPlanned = 1,
//'        accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'        piecewiseSurvivalTime = c(0, 6),
//'        lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'        gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'        accrualDuration = 12, minFollowupTime = 18,
//'        maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix nevent(
    const Rcpp::NumericVector& time = NA_REAL,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double minFollowupTime = NA_REAL,
    const double maxFollowupTime = NA_REAL) {

  auto timecpp = Rcpp::as<std::vector<double>>(time);
  auto accrualTimecpp = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualIntensitycpp = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTimecpp = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lambda1cpp = Rcpp::as<std::vector<double>>(lambda1);
  auto lambda2cpp = Rcpp::as<std::vector<double>>(lambda2);
  auto gamma1cpp = Rcpp::as<std::vector<double>>(gamma1);
  auto gamma2cpp = Rcpp::as<std::vector<double>>(gamma2);

  auto cpp_result = nevent(
    timecpp, allocationRatioPlanned, accrualTimecpp, accrualIntensitycpp,
    piecewiseTimecpp, lambda1cpp, lambda2cpp, gamma1cpp, gamma2cpp,
    accrualDuration, minFollowupTime, maxFollowupTime
  );

  return Rcpp::wrap(cpp_result);
}


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
                   const double maxFollowupTime) {

  std::size_t k = time.size();

  // truncate the analysis time by the maximum follow-up
  std::vector<double> t(k), u(k), w(k);
  for (std::size_t i = 0; i < k; ++i) {
    t[i] = std::min(time[i], accrualDuration + minFollowupTime);
    u[i] = std::min(std::max(t[i] - maxFollowupTime, 0.0), accrualDuration);
    w[i] = std::min(t[i], accrualDuration);
  }

  // Number of patients enrolled
  std::vector<double> a = accrual(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Prepare FlatMatrix for results (k rows, 2 columns)
  FlatMatrix d(k, 2);  // FlatMatrix with column-major storage

  // Precompute probabilities using pevent
  std::vector<double> s(1, maxFollowupTime); // s contains maxFollowupTime
  std::vector<double> d1 = a; // Copy of a
  std::vector<double> d2 = a; // Copy of a

  double pevent1 = pevent(s, piecewiseSurvivalTime, lambda1, gamma1)[0];
  double pevent2 = pevent(s, piecewiseSurvivalTime, lambda2, gamma2)[0];

  for (std::size_t i = 0; i < k; ++i) {
    d1[i] *= pevent1;
    d2[i] *= pevent2;
  }

  // Compute the number of patients experiencing events in each group
  for (std::size_t i = 0; i < k; ++i) {
    std::vector<double> v(1, t[i]);
    double ad1 = ad(v, u[i], w[i], accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, lambda1, gamma1)[0];
    double ad2 = ad(v, u[i], w[i], accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, lambda2, gamma2)[0];

    d(i, 0) = phi * (d1[i] + ad1);  // Active treatment group
    d(i, 1) = (1.0 - phi) * (d2[i] + ad2);  // Control group
  }

  return d;
}

//' @title Number of Subjects Having an Event by Calendar Time
//' @description Obtains the number of subjects having an event by given
//' calendar times for each treatment group.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of patients having an event.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients having an event at the
//' specified calendar times (row) for each treatment group (column).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//' nevent2(time = c(9, 24), allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, minFollowupTime = 18,
//'         maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix nevent2(
    const Rcpp::NumericVector& time = NA_REAL,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double minFollowupTime = NA_REAL,
    const double maxFollowupTime = NA_REAL) {

  auto timecpp = Rcpp::as<std::vector<double>>(time);
  auto accrualTimecpp = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualIntensitycpp = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTimecpp = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lambda1cpp = Rcpp::as<std::vector<double>>(lambda1);
  auto lambda2cpp = Rcpp::as<std::vector<double>>(lambda2);
  auto gamma1cpp = Rcpp::as<std::vector<double>>(gamma1);
  auto gamma2cpp = Rcpp::as<std::vector<double>>(gamma2);

  auto cpp_result = nevent2(
    timecpp, allocationRatioPlanned, accrualTimecpp, accrualIntensitycpp,
    piecewiseTimecpp, lambda1cpp, lambda2cpp, gamma1cpp, gamma2cpp,
    accrualDuration, minFollowupTime, maxFollowupTime
  );

  return Rcpp::wrap(cpp_result);
}


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
                     const double varianceRatio) {

  // ----------- Input Validation ----------- //
  if (std::isnan(beta) && std::isnan(IMax)) {
    throw std::invalid_argument("beta and IMax cannot be missing simultaneously");
  }
  if (!std::isnan(beta) && !std::isnan(IMax)) {
    throw std::invalid_argument("Only one of beta and IMax should be provided");
  }
  if (!std::isnan(IMax) && IMax <= 0) {
    throw std::invalid_argument("IMax must be positive");
  }
  if (std::isnan(theta)) {
    throw std::invalid_argument("theta must be provided");
  }
  if (kMax < 1) {
    throw std::invalid_argument("kMax must be a positive integer");
  }

  // Alpha and Beta must be within valid ranges
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(beta) && (beta >= 1 - alpha || beta < 0.0001)) {
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");
  }

  std::string unknown = std::isnan(beta) ? "beta" : "IMax";

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (static_cast<int>(informationRates.size()) != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[kMax-1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (int i = 0; i < kMax; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(kMax);
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (static_cast<int>(efficacyStopping.size()) != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[kMax-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (static_cast<int>(futilityStopping.size()) != kMax)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping[kMax-1] != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds);

  if (!missingCriticalValues) {
    if (static_cast<int>(criticalValues.size()) != kMax) {
      throw std::invalid_argument("Invalid length for criticalValues");
    }
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "user" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  if (missingCriticalValues && asf == "user") {
    if (!none_na(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be specified");
    if (static_cast<int>(userAlphaSpending.size()) != kMax)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[kMax-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (!(static_cast<int>(futilityBounds.size()) == kMax - 1 ||
        static_cast<int>(futilityBounds.size()) == kMax)) {
      throw std::invalid_argument("Invalid length for futilityBounds");
    }
  }
  if (!missingCriticalValues && !missingFutilityBounds) {
    for (int i = 0; i < kMax - 1; ++i) {
      if (futilityBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
    if (static_cast<int>(futilityBounds.size()) == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      throw std::invalid_argument(
          "futilityBounds must meet criticalValues at the final look");
    }
  }

  std::string bsf = typeBetaSpending;
  for (char &c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (unknown == "IMax") {
    if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
        bsf == "sfkd" || bsf == "sfhsd" || bsf == "user" || bsf == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpending");
    }
  } else {
    if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
        bsf == "sfkd" || bsf == "sfhsd" || bsf == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpending");
    }
  }

  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument ("parameterBetaSpending must be positive for sfKD");
  }

  if (unknown == "IMax" && bsf == "user") {
    if (!none_na(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be specified");
    if (static_cast<int>(userBetaSpending.size()) != kMax)
      throw std::invalid_argument("Invalid length of userBetaSpending");
    if (userBetaSpending[0] < 0.0)
      throw std::invalid_argument("userBetaSpending must be nonnegative");
    if (any_nonincreasing(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be nondecreasing");
    if (userBetaSpending[kMax-1] != beta)
      throw std::invalid_argument("userBetaSpending must end with specified beta");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (static_cast<int>(spendingTime.size()) != kMax)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[kMax-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  if (varianceRatio <= 0.0) {
    throw std::invalid_argument("varianceRatio must be positive");
  }
  // ----------- End of Input Validation ----------- //

  // set up efficacy bounds
  std::vector<double> l(kMax, -6.0), zero(kMax, 0.0);
  std::vector<double> critValues = criticalValues;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (kMax > 1 && static_cast<int>(criticalValues.size()) == kMax) {
      bool hasNaN = false;
      for (int i = 0; i < kMax - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[kMax-1])) {
        haybittle = true;
      }
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(kMax);
      for (int i = 0; i < kMax - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = cumAlphaSpent[kMax-1];

  // set up futility bounds
  std::vector<double> futBounds = futilityBounds;
  if (kMax > 1) {
    if (missingFutilityBounds && bsf == "none") {
      futBounds = std::vector<double>(kMax, -6.0);
      futBounds[kMax-1] = critValues[kMax-1];
    } else if (!missingFutilityBounds &&
      static_cast<int>(futBounds.size()) == kMax-1) {
      futBounds.push_back(critValues[kMax-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // multiplier for the boundaries under the alternative hypothesis
  std::vector<double> w(kMax, std::sqrt(varianceRatio));
  std::vector<double> u(kMax);
  for (int i = 0; i < kMax; ++i) {
    u[i] = critValues[i] * w[i];
  }

  double beta1 = beta;
  double IMax1 = IMax;
  double drift;
  if (unknown == "IMax") {
    std::vector<double> u1; u1.reserve(kMax);
    std::vector<double> l1; l1.reserve(kMax);
    double sqrtt0 = std::sqrt(infoRates[0]);

    auto f = [&](double aval)->double {
      std::vector<double> delta = std::vector<double>(kMax, aval);

      // compute stagewise exit probabilities
      if (!missingFutilityBounds || bsf == "none" || kMax == 1) {
        for (int i = 0; i < kMax; ++i) {
          l[i] = futBounds[i] * w[i];
        }
        ListCpp probs = exitprobcpp(u, l, delta, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
        return overallReject - (1.0 - beta);
      } else {
        // initialize futility bound to be updated
        futBounds = std::vector<double>(kMax);
        double eps = 0.0;

        // first stage
        double cb = (bsf == "user") ? userBetaSpending[0] :
          errorSpentcpp(spendTime[0], beta, bsf, parameterBetaSpending);

        if (!futStopping[0]) {
          futBounds[0] = -6.0;
        } else {
          double dt0 = delta[0] * sqrtt0;
          eps = boost_pnorm(u[0] - dt0) - cb;
          if (eps < 0.0) return -1.0; // to decrease drift
          futBounds[0] = (boost_qnorm(cb) + dt0) / w[0];
        }

        // subsequent stages
        for (int k = 1; k < kMax; ++k) {
          l[k-1] = futBounds[k-1] * w[k-1];
          cb = (bsf == "user") ? userBetaSpending[k] :
            errorSpentcpp(spendTime[k], beta, bsf, parameterBetaSpending);

          if (!futStopping[k]) {
            futBounds[k] = -6.0;
          } else {
            u1.resize(k + 1);
            l1.resize(k + 1);

            std::memcpy(u1.data(), u.data(), k * sizeof(double));
            u1[k] = 6.0;
            std::memcpy(l1.data(), l.data(), k * sizeof(double));

            // lambda expression for finding futility bound at stage k
            auto g = [&](double aval)->double {
              l1[k] = aval * w[k];
              ListCpp probs = exitprobcpp(u1, l1, delta, infoRates);
              auto v = probs.get<std::vector<double>>("exitProbLower");
              double cpl = std::accumulate(v.begin(), v.end(), 0.0);
              return cpl - cb;
            };

            double bk = critValues[k];
            eps = g(bk);
            double g_minus6 = g(-6.0);

            if (g_minus6 > 0.0) { // no beta spent at current visit
              futBounds[k] = -6.0;
            } else if (eps > 0.0) {
              auto g_for_brent = [&](double x)->double {
                if (x == -6.0) return g_minus6;  // avoid recomputation at 6.0
                if (x == bk) return eps;         // avoid recomputation at b[k]
                return g(x);
              };

              futBounds[k] = brent(g_for_brent, -6.0, bk, 1e-6);
            } else if (k < kMax-1) {
              return -1.0;
            }
          }
        }

        return eps;
      }
    };

    drift = brent(f, 0.0, 6.0, 1e-6);
    IMax1 = sq(drift / theta);
    futBounds[kMax-1] = critValues[kMax-1];
    l[kMax-1] = futBounds[kMax-1] * w[kMax-1];
    std::vector<double> delta = std::vector<double>(kMax, drift);
    probs = exitprobcpp(u, l, delta, infoRates);
  } else {
    drift = theta * std::sqrt(IMax1);
    std::vector<double> delta = std::vector<double>(kMax, drift);
    if (!missingFutilityBounds || bsf=="none" || kMax==1) {
      for (int i = 0; i < kMax; ++i) {
        l[i] = futBounds[i] * w[i];
      }
      probs = exitprobcpp(u, l, delta, infoRates);
      auto v = probs.get<std::vector<double>>("exitProbUpper");
      double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
      beta1 = 1.0 - overallReject;
    } else {
      ListCpp out = getPower(alpha1, kMax, critValues, delta, infoRates, bsf,
                             parameterBetaSpending, spendTime, futStopping, w);
      double overallReject = out.get<double>("power");
      beta1 = 1.0 - overallReject;
      futBounds = out.get<std::vector<double>>("futilityBounds");
      for (int i = 0; i < kMax; ++i) {
        l[i] = futBounds[i] * w[i];
      }
      probs = out.get_list("probs");
    }
  }

  double driftf = boost_qnorm(1.0 - alpha1) * w[0] + boost_qnorm(1.0 - beta1);
  double inflationFactor = sq(drift / driftf);

  // output the results
  std::vector<double> information(kMax);
  std::vector<double> efficacyTheta(kMax);
  std::vector<double> futilityTheta(kMax);
  std::vector<double> efficacyP(kMax);
  std::vector<double> futilityP(kMax);
  for (int i = 0; i < kMax; ++i) {
    information[i] = IMax1 * infoRates[i];
    efficacyTheta[i] = u[i] / std::sqrt(information[i]);
    futilityTheta[i] = l[i] / std::sqrt(information[i]);
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
    futilityP[i] = 1.0 - boost_pnorm(futBounds[i]);
  }

  // stagewise exit probabilities under H1
  auto pu = probs.get<std::vector<double>>("exitProbUpper");
  auto pl = probs.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpu(kMax), cpl(kMax);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());
  double overallReject = cpu[kMax-1];
  std::vector<double> ptotal(kMax);
  for (int i = 0; i < kMax; ++i) ptotal[i] = pu[i] + pl[i];
  double expectedInformationH1 = std::inner_product(
    ptotal.begin(), ptotal.end(), information.begin(), 0.0);

  // stagewise exit probabilities under H0 with binding futility
  ListCpp probsH0 = exitprobcpp(critValues, futBounds, zero, infoRates);
  auto puH0 = probsH0.get<std::vector<double>>("exitProbUpper");
  auto plH0 = probsH0.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpuH0(kMax), cplH0(kMax);
  std::partial_sum(puH0.begin(), puH0.end(), cpuH0.begin());
  std::partial_sum(plH0.begin(), plH0.end(), cplH0.begin());
  double overallRejectH0 = cpuH0[kMax-1];
  std::vector<double> ptotalH0(kMax);
  for (int i = 0; i < kMax; ++i) ptotalH0[i] = puH0[i] + plH0[i];
  double expectedInformationH0 = std::inner_product(
    ptotalH0.begin(), ptotalH0.end(), information.begin(), 0.0);

  for (int i = 0; i < kMax; ++i) {
    if (critValues[i] == 6) effStopping[i] = 0;
    if (futBounds[i] == -6) futStopping[i] = 0;
  }

  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(futBounds), "futilityBounds");
  byStageResults.push_back(std::move(pu), "rejectPerStage");
  byStageResults.push_back(std::move(pl), "futilityPerStage");
  byStageResults.push_back(std::move(cpu), "cumulativeRejection");
  byStageResults.push_back(std::move(cpl), "cumulativeFutility");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(efficacyTheta), "efficacyTheta");
  byStageResults.push_back(std::move(futilityTheta), "futilityTheta");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(futilityP), "futilityP");
  byStageResults.push_back(std::move(information), "information");
  byStageResults.push_back(std::move(effStopping), "efficacyStopping");
  byStageResults.push_back(std::move(futStopping), "futilityStopping");
  byStageResults.push_back(std::move(puH0), "rejectPerStageH0");
  byStageResults.push_back(std::move(plH0), "futilityPerStageH0");
  byStageResults.push_back(std::move(cpuH0), "cumulativeRejectionH0");
  byStageResults.push_back(std::move(cplH0), "cumulativeFutilityH0");

  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha1, "alpha");
  overallResults.push_back(overallRejectH0, "attainedAlpha");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(theta, "theta");
  overallResults.push_back(IMax1, "information");
  overallResults.push_back(expectedInformationH1, "expectedInformationH1");
  overallResults.push_back(expectedInformationH0, "expectedInformationH0");
  overallResults.push_back(drift, "drift");
  overallResults.push_back(inflationFactor, "inflationFactor");

  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(typeBetaSpending, "typeBetaSpending");
  settings.push_back(parameterBetaSpending, "parameterBetaSpending");
  settings.push_back(userBetaSpending, "userBetaSpending");
  settings.push_back(spendingTime, "spendingTime");
  settings.push_back(varianceRatio, "varianceRatio");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");
  return result;
}

//' @title Power and Sample Size for a Generic Group Sequential Design
//' @description Obtains the maximum information and stopping boundaries
//' for a generic group sequential design assuming a constant treatment
//' effect, or obtains the power given the maximum information and
//' stopping boundaries.
//'
//' @param beta The type II error.
//' @param IMax The maximum information. Either \code{beta} or \code{IMax}
//'   should be provided while the other one should be missing.
//' @param theta The parameter value.
//' @inheritParams param_kMax
//' @param informationRates The information rates. Fixed prior to the trial.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param varianceRatio The ratio of the variance under H0 to the
//'   variance under H1.
//'
//' @return An S3 class \code{design} object with three components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{attainedAlpha}: The attained significance level, which is
//'       different from the overall significance level in the presence of
//'       futility stopping.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{theta}: The parameter value.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedInformationH1}: The expected information under H1.
//'
//'     - \code{expectedInformationH0}: The expected information under H0.
//'
//'     - \code{drift}: The drift parameter, equal to
//'       \code{theta*sqrt(information)}.
//'
//'     - \code{inflationFactor}: The inflation factor (relative to the
//'       fixed design).
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'
//'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{futilityPerStage}: The probability for futility stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeFutility}: The cumulative probability for futility
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'
//'     - \code{efficacyTheta}: The efficacy boundaries on the parameter
//'       scale.
//'
//'     - \code{futilityTheta}: The futility boundaries on the parameter
//'       scale.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{futilityP}: The futility boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'
//'     - \code{futilityStopping}: Whether to allow futility stopping.
//'
//'     - \code{rejectPerStageH0}: The probability for efficacy stopping
//'       under H0.
//'
//'     - \code{futilityPerStageH0}: The probability for futility stopping
//'       under H0.
//'
//'     - \code{cumulativeRejectionH0}: The cumulative probability for
//'       efficacy stopping under H0.
//'
//'     - \code{cumulativeFutilityH0}: The cumulative probability for
//'       futility stopping under H0.
//'
//' * \code{settings}: A list containing the following input parameters:
//'
//'     - \code{typeAlphaSpending}: The type of alpha spending.
//'
//'     - \code{parameterAlphaSpending}: The parameter value for alpha
//'       spending.
//'
//'     - \code{userAlphaSpending}: The user defined alpha spending.
//'
//'     - \code{typeBetaSpending}: The type of beta spending.
//'
//'     - \code{parameterBetaSpending}: The parameter value for beta
//'       spending.
//'
//'     - \code{userBetaSpending}: The user defined beta spending.
//'
//'     - \code{spendingTime}: The error spending time at each analysis.
//'
//'     - \code{varianceRatio}: The ratio of the variance under H0
//'       to the variance under H1.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Christopher Jennison, Bruce W. Turnbull.
//' Group Sequential Methods with Applications to Clinical Trials.
//' Chapman & Hall/CRC: Boca Raton, 2000, ISBN:0849303168
//'
//' @examples
//'
//' # Example 1: obtain the maximum information given power
//' (design1 <- getDesign(
//'   beta = 0.2, theta = -log(0.7),
//'   kMax = 2, informationRates = c(0.5,1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   typeBetaSpending = "sfP"))
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesign(
//'   IMax = 72.5, theta = -log(0.7),
//'   kMax = 3, informationRates = c(0.5, 0.75, 1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   typeBetaSpending = "sfP"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getDesign(
    const double beta = NA_REAL,
    const double IMax = NA_REAL,
    const double theta = NA_REAL,
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& futilityBounds = NA_REAL,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const Rcpp::NumericVector& userBetaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const double varianceRatio = 1) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto userBeta = Rcpp::as<std::vector<double>>(userBetaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  auto cpp_result = getDesigncpp(
    beta, IMax, theta, kMax, infoRates, effStopping, futStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, futBounds,
    typeBetaSpending, parameterBetaSpending, userBeta, spendTime, varianceRatio
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "design";
  return result;
}


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
                          const std::vector<double>& spendingTime) {

  // ----------- Input Validation ----------- //
  if (std::isnan(beta) && std::isnan(IMax)) {
    throw std::invalid_argument("beta and IMax cannot be missing simultaneously");
  }
  if (!std::isnan(beta) && !std::isnan(IMax)) {
    throw std::invalid_argument("Only one of beta and IMax should be provided");
  }
  if (!std::isnan(IMax) && IMax <= 0) {
    throw std::invalid_argument("IMax must be positive");
  }
  if (std::isnan(theta)) {
    throw std::invalid_argument("theta must be provided");
  }
  if (std::isnan(thetaLower)) {
    throw std::invalid_argument("thetaLower must be provided");
  }
  if (std::isnan(thetaUpper)) {
    throw std::invalid_argument("thetaUpper must be provided");
  }
  if (thetaLower >= theta) {
    throw std::invalid_argument("thetaLower must be less than theta");
  }
  if (thetaUpper <= theta) {
    throw std::invalid_argument("thetaUpper must be greater than theta");
  }
  if (kMax < 1) {
    throw std::invalid_argument("kMax must be a positive integer");
  }

  // Alpha and Beta must be within valid ranges
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(beta) && (beta >= 1 - alpha || beta < 0.0001)) {
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");
  }

  std::string unknown = std::isnan(beta) ? "beta" : "IMax";


  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (static_cast<int>(informationRates.size()) != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[kMax-1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (int i = 0; i < kMax; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(kMax);
  }


  bool missingCriticalValues = !none_na(criticalValues);

  if (!missingCriticalValues) {
    if (static_cast<int>(criticalValues.size()) != kMax) {
      throw std::invalid_argument("Invalid length for criticalValues");
    }
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "user" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  if (missingCriticalValues && asf == "user") {
    if (!none_na(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be specified");
    if (static_cast<int>(userAlphaSpending.size()) != kMax)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[kMax-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (static_cast<int>(spendingTime.size()) != kMax)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[kMax-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // ----------- End of Input Validation ----------- //

  std::vector<double> u(kMax), l(kMax, -6.0), zero(kMax, 0.0);
  std::vector<double> critValues = criticalValues;

  // obtain criticalValues
  if (missingCriticalValues) {
    bool haybittle = false;
    if (kMax > 1 && static_cast<int>(criticalValues.size()) == kMax) {
      bool hasNaN = false;
      for (int i = 0; i < kMax - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[kMax-1])) {
        haybittle = true;
      }
    }

    if (haybittle) { // Haybittle & Peto
      for (int i = 0; i < kMax - 1; ++i) {
        u[i] = criticalValues[i];
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      std::vector<unsigned char> effStopping(kMax, 1);
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  std::vector<double> li(kMax, -6.0), ui(kMax, 6.0);
  ListCpp probs = exitprobcpp(critValues, li, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());

  std::vector<double> efficacyP(kMax);
  for (int i = 0; i < kMax; ++i) {
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
  }

  // calculate cumulative rejection probability under H1
  double deltaLower = thetaLower - theta;
  double deltaUpper = thetaUpper - theta;

  // obtain IMax if needed
  double IMax1 = IMax;
  std::vector<double> I(kMax);
  std::vector<double> b(kMax);
  std::vector<double> a(kMax);
  std::vector<double> cpl(kMax);
  std::vector<double> cpu(kMax);

  if (unknown == "IMax") {
    auto f = [&](double aval)->double {
      for (int i = 0; i < kMax; ++i) {
        I[i] = infoRates[i] * aval;
        double sqrtIi = std::sqrt(I[i]);
        l[i] = critValues[i] + deltaLower * sqrtIi;
        u[i] = -critValues[i] + deltaUpper * sqrtIi;
        b[i] = std::max(l[i], li[i]);
        a[i] = std::min(u[i], ui[i]);
      }

      ListCpp probs1 = exitprobcpp(b, li, zero, I);
      ListCpp probs2 = exitprobcpp(ui, a, zero, I);
      auto v1 = probs1.get<std::vector<double>>("exitProbUpper");
      auto v2 = probs2.get<std::vector<double>>("exitProbLower");
      std::partial_sum(v1.begin(), v1.end(), cpl.begin());
      std::partial_sum(v2.begin(), v2.end(), cpu.begin());
      double p1 = cpl[kMax-1];
      double p2 = cpu[kMax-1];

      bool cross = false;
      for (int i = 0; i < kMax; ++i) {
        if (l[i] <= u[i]) { cross = true; break; }
      }

      double power;
      if (cross) {
        power = p1 + p2 - 1.0;
      } else {
        ListCpp probs = exitprobcpp(l, u, zero, I);
        auto v1x = probs.get<std::vector<double>>("exitProbUpper");
        auto v2x = probs.get<std::vector<double>>("exitProbLower");
        double p1x = std::accumulate(v1x.begin(), v1x.end(), 0.0);
        double p2x = std::accumulate(v2x.begin(), v2x.end(), 0.0);
        power = p1 + p2 - p1x - p2x;
      }

      return power - (1.0 - beta);
    };

    double z0 = boost_qnorm(1.0 - alpha);
    double z1 = boost_qnorm(1.0 - beta);
    double IMax0 = sq((z0 + z1) / deltaLower);
    double IMaxLower = 0.5 * IMax0;
    double IMaxUpper = 1.5 * IMax0;
    IMax1 = brent(f, IMaxLower, IMaxUpper, 1e-6);
  }

  // cumulative rejection probability under H1
  for (int i = 0; i < kMax; ++i) {
    I[i] = infoRates[i] * IMax1;
    double sqrtIi = std::sqrt(I[i]);
    l[i] = critValues[i] + deltaLower * sqrtIi;
    u[i] = -critValues[i] + deltaUpper * sqrtIi;
    b[i] = std::max(l[i], li[i]);
    a[i] = std::min(u[i], ui[i]);
  }

  ListCpp probs1 = exitprobcpp(b, li, zero, I);
  ListCpp probs2 = exitprobcpp(ui, a, zero, I);
  auto v1 = probs1.get<std::vector<double>>("exitProbUpper");
  auto v2 = probs2.get<std::vector<double>>("exitProbLower");
  std::partial_sum(v1.begin(), v1.end(), cpl.begin());
  std::partial_sum(v2.begin(), v2.end(), cpu.begin());

  std::vector<unsigned char> nocross(kMax);
  for (int i = 0; i < kMax; ++i) {
    nocross[i] = (l[i] >= u[i]) ? 1 : 0;
  }
  std::vector<int> k = which(nocross);

  std::vector<double> cp(kMax);
  if (k.empty()) {
    for (int i = 0; i < kMax; ++i) {
      cp[i] = cpl[i] + cpu[i] - 1.0;
    }
  } else {
    int K = *std::max_element(k.begin(), k.end());
    std::vector l1 = subset(l, 0, K+1);
    std::vector u1 = subset(u, 0, K+1);
    std::vector d1 = subset(zero, 0, K+1);
    std::vector I1 = subset(I, 0, K+1);
    ListCpp probs = exitprobcpp(l1, u1, d1, I1);
    auto v1x = probs.get<std::vector<double>>("exitProbUpper");
    auto v2x = probs.get<std::vector<double>>("exitProbLower");
    std::vector<double> cplx(kMax), cpux(kMax);
    std::partial_sum(v1x.begin(), v1x.end(), cplx.begin());
    std::partial_sum(v2x.begin(), v2x.end(), cpux.begin());
    for (int i = 0; i <= K; ++i) {
      cp[i] = cpl[i] + cpu[i] - cplx[i] - cpux[i];
    }
    for (int i = K+1; i < kMax; ++i) {
      cp[i] = cpl[i] + cpu[i] - 1.0;
    }
  }

  // incremental exit probabilities under H1
  std::vector<double> q(kMax);
  q[0] = cp[0];
  for (int i = 1; i < kMax - 1; ++i) q[i] = cp[i] - cp[i-1];
  q[kMax-1] = 1.0 - cp[kMax-2];

  std::vector<double> rejectPerStage(kMax);
  rejectPerStage[0] = cp[0];
  for (int i = 1; i < kMax; ++i) {
    rejectPerStage[i] = cp[i] - cp[i-1];
  }

  double overallReject = cp[kMax-1];
  double expectedInformationH1 = std::inner_product(
    q.begin(), q.end(), I.begin(), 0.0);

  std::vector<double> efficacyThetaLower(kMax), efficacyThetaUpper(kMax);
  for (int i = 0; i < kMax; ++i) {
    double thetaBound = critValues[i] / std::sqrt(I[i]);
    efficacyThetaLower[i] = thetaBound + thetaLower;
    efficacyThetaUpper[i] = -thetaBound + thetaUpper;
  }

  std::vector<double> theta10(kMax, deltaLower);
  std::vector<double> theta20(kMax, deltaUpper);

  // cumulative rejection probability under H10
  ListCpp probsH10 = exitprobcpp(ui, a, theta10, I);
  auto vH10 = probsH10.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpuH10(kMax);
  std::partial_sum(vH10.begin(), vH10.end(), cpuH10.begin());
  std::vector<double> cplH10 = cumAlphaSpent;

  std::vector<double> cpH10(kMax);
  if (k.empty()) {
    for (int i = 0; i < kMax; ++i) {
      cpH10[i] = cplH10[i] + cpuH10[i] - 1.0;
    }
  } else {
    int K = *std::max_element(k.begin(), k.end());
    std::vector l1 = subset(l, 0, K+1);
    std::vector u1 = subset(u, 0, K+1);
    std::vector d1 = subset(theta10, 0, K+1);
    std::vector I1 = subset(I, 0, K+1);
    ListCpp probs = exitprobcpp(l1, u1, d1, I1);
    auto v1x = probs.get<std::vector<double>>("exitProbUpper");
    auto v2x = probs.get<std::vector<double>>("exitProbLower");
    std::vector<double> cplH10x(kMax), cpuH10x(kMax);
    std::partial_sum(v1x.begin(), v1x.end(), cplH10x.begin());
    std::partial_sum(v2x.begin(), v2x.end(), cpuH10x.begin());
    for (int i = 0; i <= K; ++i) {
      cpH10[i] = cplH10[i] + cpuH10[i] - cplH10x[i] - cpuH10x[i];
    }
    for (int i = K+1; i < kMax; ++i) {
      cpH10[i] = cplH10[i] + cpuH10[i] - 1.0;
    }
  }

  // incremental exit probabilities under H10
  std::vector<double> qH10(kMax);
  qH10[0] = cpH10[0];
  for (int i = 1; i < kMax - 1; ++i) qH10[i] = cpH10[i] - cpH10[i-1];
  qH10[kMax-1] = 1.0 - cpH10[kMax-2];

  double attainedAlphaH10 = cpH10[kMax-1];
  double expectedInformationH10 = std::inner_product(
    qH10.begin(), qH10.end(), I.begin(), 0.0);

  // cumulative rejection probability under H20
  ListCpp probsH20 = exitprobcpp(b, li, theta20, I);
  auto vH20 = probsH20.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cplH20(kMax);
  std::partial_sum(vH20.begin(), vH20.end(), cplH20.begin());
  std::vector<double> cpuH20 = cumAlphaSpent;

  std::vector<double> cpH20(kMax);
  if (k.empty()) {
    for (int i = 0; i < kMax; ++i) {
      cpH20[i] = cplH20[i] + cpuH20[i] - 1.0;
    }
  } else {
    int K = *std::max_element(k.begin(), k.end());
    std::vector l1 = subset(l, 0, K+1);
    std::vector u1 = subset(u, 0, K+1);
    std::vector d1 = subset(theta20, 0, K+1);
    std::vector I1 = subset(I, 0, K+1);
    ListCpp probs = exitprobcpp(l1, u1, d1, I1);
    auto v1x = probs.get<std::vector<double>>("exitProbUpper");
    auto v2x = probs.get<std::vector<double>>("exitProbLower");
    std::vector<double> cplH20x(kMax), cpuH20x(kMax);
    std::partial_sum(v1x.begin(), v1x.end(), cplH20x.begin());
    std::partial_sum(v2x.begin(), v2x.end(), cpuH20x.begin());
    for (int i = 0; i <= K; ++i) {
      cpH20[i] = cplH20[i] + cpuH20[i] - cplH20x[i] - cpuH20x[i];
    }
    for (int i = K+1; i < kMax; ++i) {
      cpH20[i] = cplH20[i] + cpuH20[i] - 1.0;
    }
  }

  // incremental exit probabilities under H20
  std::vector<double> qH20(kMax);
  qH20[0] = cpH20[0];
  for (int i = 1; i < kMax - 1; ++i) qH20[i] = cpH20[i] - cpH20[i-1];
  qH20[kMax-1] = 1.0 - cpH20[kMax-2];

  double attainedAlphaH20 = cpH20[kMax-1];
  double expectedInformationH20 = std::inner_product(
    qH20.begin(), qH20.end(), I.begin(), 0.0);

  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(rejectPerStage), "rejectPerStage");
  byStageResults.push_back(std::move(cp), "cumulativeRejection");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(cpH10), "cumulativeAttainedAlphaH10");
  byStageResults.push_back(std::move(cpH20), "cumulativeAttainedAlphaH20");
  byStageResults.push_back(std::move(efficacyThetaLower), "efficacyThetaLower");
  byStageResults.push_back(std::move(efficacyThetaUpper), "efficacyThetaUpper");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(I), "information");

  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha, "alpha");
  overallResults.push_back(attainedAlphaH10, "attainedAlphaH10");
  overallResults.push_back(attainedAlphaH20, "attainedAlphaH20");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(thetaLower, "thetaLower");
  overallResults.push_back(thetaUpper, "thetaUpper");
  overallResults.push_back(theta, "theta");
  overallResults.push_back(IMax1, "information");
  overallResults.push_back(expectedInformationH1, "expectedInformationH1");
  overallResults.push_back(expectedInformationH10, "expectedInformationH10");
  overallResults.push_back(expectedInformationH20, "expectedInformationH20");

  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(spendingTime, "spendingTime");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");
  return result;
}

//' @title Power and Sample Size for a Generic Group Sequential Equivalence
//' Design
//'
//' @description Obtains the maximum information and stopping boundaries
//' for a generic group sequential equivalence design assuming a constant
//' treatment effect, or obtains the power given the maximum information
//' and stopping boundaries.
//'
//' @param beta The type II error.
//' @param IMax The maximum information. Either \code{beta} or \code{IMax}
//'   should be provided while the other one should be missing.
//' @param thetaLower The parameter value at the lower equivalence limit.
//' @param thetaUpper The parameter value at the upper equivalence limit.
//' @param theta The parameter value under the alternative hypothesis.
//' @inheritParams param_kMax
//' @param informationRates The information rates. Fixed prior to the trial.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_criticalValues
//' @param alpha The significance level for each of the two one-sided
//'   tests, e.g., 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//'
//' @details
//' Consider the equivalence design with two one-sided hypotheses:
//' \deqn{H_{10}: \theta \leq \theta_{10},}
//' \deqn{H_{20}: \theta \geq \theta_{20}.}
//' We reject \eqn{H_{10}} at or before look \eqn{k} if
//' \deqn{Z_{1j} = (\hat{\theta}_j - \theta_{10})\sqrt{I_j}
//' \geq b_j}
//' for some \eqn{j=1,\ldots,k}, where \eqn{\{b_j:j=1,\ldots,K\}} are the
//' critical values associated with the specified alpha-spending function,
//' and \eqn{I_j} is the information for \eqn{\theta} (inverse variance of
//' \eqn{\hat{\theta}}) at the
//' \eqn{j}th look. For example,
//' for estimating the risk difference \eqn{\theta = \pi_1 - \pi_2},
//' \deqn{I_j = \left\{\frac{\pi_1 (1-\pi_1)}{n_{1j}} +
//' \frac{\pi_2(1-\pi_2)}{n_{2j}}\right\}^{-1}.}
//' It follows that
//' \deqn{(Z_{1j} \geq b_j) = (Z_j \geq b_j +
//' \theta_{10}\sqrt{I_j}),}
//' where \eqn{Z_j = \hat{\theta}_j \sqrt{I_j}}.
//'
//' Similarly, we reject \eqn{H_{20}} at or before look \eqn{k} if
//' \deqn{Z_{2j} = (\hat{\theta}_j - \theta_{20})\sqrt{I_j}
//' \leq -b_j} for some \eqn{j=1,\ldots,k}. We have
//' \deqn{(Z_{2j} \leq -b_j) = (Z_j \leq - b_j +
//' \theta_{20}\sqrt{I_j}).}
//'
//' Let \eqn{l_j = b_j + \theta_{10}\sqrt{I_j}},
//' and \eqn{u_j = -b_j + \theta_{20}\sqrt{I_j}}.
//' The cumulative probability to reject \eqn{H_0 = H_{10} \cup H_{20}} at
//' or before look \eqn{k} under the alternative hypothesis \eqn{H_1} is
//' given by
//' \deqn{P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j) \cap
//' \cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right) = p_1 + p_2 - p_{12},}
//' where
//' \deqn{p_1 = P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \geq l_j)\right),}
//' \deqn{p_2 = P_\theta\left(\cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \leq u_j)\right),}
//' and
//' \deqn{p_{12} = P_\theta\left(\cup_{j=1}^{k} (Z_j \geq l_j) \cup
//' (Z_j \leq u_j)\right).}
//' Of note, both \eqn{p_1} and \eqn{p_2} can be evaluated using
//' one-sided exit probabilities for group sequential designs.
//' If there exists \eqn{j\leq k} such that \eqn{l_j \leq u_j}, then
//' \eqn{p_{12} = 1}. Otherwise, \eqn{p_{12}} can be evaluated using
//' two-sided exit probabilities for group sequential designs.
//'
//' Since the equivalent hypothesis is tested using two one-sided tests,
//' the type I error is controlled. To evaluate the attained type I error
//' of the equivalence trial under \eqn{H_{10}} (or \eqn{H_{20}}),
//' we simply fix the control group parameters, update the active
//' treatment group parameters according to the null hypothesis, and
//' use the parameters in the power calculation outlined above.
//'
//' @return An S3 class \code{designEquiv} object with three components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{attainedAlphaH10}: The attained significance level under H10.
//'
//'     - \code{attainedAlphaH20}: The attained significance level under H20.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{thetaLower}: The parameter value at the lower equivalence
//'       limit.
//'
//'     - \code{thetaUpper}: The parameter value at the upper equivalence
//'       limit.
//'
//'     - \code{theta}: The parameter value under the alternative hypothesis.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedInformationH1}: The expected information under H1.
//'
//'     - \code{expectedInformationH10}: The expected information under H10.
//'
//'     - \code{expectedInformationH20}: The expected information under H20.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
//'       each of the two one-sided tests.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
//'       the two one-sided tests.
//'
//'     - \code{cumulativeAttainedAlphaH10}: The cumulative probability for
//'       efficacy stopping under H10.
//'
//'     - \code{cumulativeAttainedAlphaH20}: The cumulative probability for
//'       efficacy stopping under H20.
//'
//'     - \code{efficacyThetaLower}: The efficacy boundaries on the
//'       parameter scale for the one-sided null hypothesis at the
//'       lower equivalence limit.
//'
//'     - \code{efficacyThetaUpper}: The efficacy boundaries on the
//'       parameter scale for the one-sided null hypothesis at the
//'       upper equivalence limit.
//'
//'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
//'       each of the two one-sided tests.
//'
//'     - \code{information}: The cumulative information.
//'
//' * \code{settings}: A list containing the following components:
//'
//'     - \code{typeAlphaSpending}: The type of alpha spending.
//'
//'     - \code{parameterAlphaSpending}: The parameter value for alpha
//'       spending.
//'
//'     - \code{userAlphaSpending}: The user defined alpha spending.
//'
//'     - \code{spendingTime}: The error spending time at each analysis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' # Example 1: obtain the maximum information given power
//' (design1 <- getDesignEquiv(
//'   beta = 0.2, thetaLower = log(0.8), thetaUpper = log(1.25),
//'   kMax = 2, informationRates = c(0.5, 1),
//'   alpha = 0.05, typeAlphaSpending = "sfOF"))
//'
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesignEquiv(
//'   IMax = 72.5, thetaLower = log(0.7), thetaUpper = -log(0.7),
//'   kMax = 3, informationRates = c(0.5, 0.75, 1),
//'   alpha = 0.05, typeAlphaSpending = "sfOF"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getDesignEquiv(
    const double beta = NA_REAL,
    const double IMax = NA_REAL,
    const double thetaLower = NA_REAL,
    const double thetaUpper = NA_REAL,
    const double theta = 0,
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.05,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  auto cpp_result = getDesignEquivcpp(
    beta, IMax, thetaLower, thetaUpper, theta, kMax, infoRates, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, spendTime
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "designEquiv";
  return result;
}


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
                       const double varianceRatio) {

  // ----------- Input Validation ----------- //
  if (std::isnan(betaNew) && std::isnan(INew)) {
    throw std::invalid_argument("betaNew and INew cannot be missing simultaneously");
  }
  if (!std::isnan(betaNew) && !std::isnan(INew)) {
    throw std::invalid_argument("Only one of betaNew and INew should be provided");
  }
  if (!std::isnan(INew) && INew <= 0.0) {
    throw std::invalid_argument("INew must be positive");
  }
  if (std::isnan(theta)) throw std::invalid_argument("theta must be provided");
  if (L < 1) throw std::invalid_argument("L must be a positive integer");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  // Alpha and Beta must be within valid ranges
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(betaNew) && (betaNew < 0.0001 || betaNew >= 1)) {
    throw std::invalid_argument("betaNew must lie in [0.0001, 1)");
  }

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (static_cast<int>(informationRates.size()) != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[kMax-1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (int i = 0; i < kMax; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(kMax);
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (static_cast<int>(efficacyStopping.size()) != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[kMax-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (static_cast<int>(futilityStopping.size()) != kMax)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping[kMax-1] != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds);

  if (!missingCriticalValues) {
    if (static_cast<int>(criticalValues.size()) != kMax) {
      throw std::invalid_argument("Invalid length for criticalValues");
    }
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "user" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  if (missingCriticalValues && asf == "user") {
    if (!none_na(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be specified");
    if (static_cast<int>(userAlphaSpending.size()) != kMax)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[kMax-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (!(static_cast<int>(futilityBounds.size()) == kMax - 1 ||
        static_cast<int>(futilityBounds.size()) == kMax)) {
      throw std::invalid_argument("Invalid length for futilityBounds");
    }
  }
  if (!missingCriticalValues && !missingFutilityBounds) {
    for (int i = 0; i < kMax - 1; ++i) {
      if (futilityBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
    if (static_cast<int>(futilityBounds.size()) == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      throw std::invalid_argument(
          "futilityBounds must meet criticalValues at the final look");
    }
  }

  std::string bsf = typeBetaSpending;
  for (char &c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
      bsf == "sfkd" || bsf == "sfhsd" || bsf == "none")) {
    throw std::invalid_argument("Invalid value for typeBetaSpending");
  }

  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument ("parameterBetaSpending must be positive for sfKD");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (static_cast<int>(spendingTime.size()) != kMax)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[kMax-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // ----------- New Design Input Validation ----------- //
  std::vector<double> infoRatesNew;
  std::vector<unsigned char> effStoppingNew;
  std::vector<unsigned char> futStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  std::string bsfNew = typeBetaSpendingNew;
  for (char &c : bsfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (MullerSchafer) {
    if (kNew < 1) {
      throw std::invalid_argument("kNew must be a positive integer");
    }

    // informationRatesNew: default to (1:kNew)/kNew if missing
    infoRatesNew.resize(kNew);
    if (none_na(informationRatesNew)) {
      if (static_cast<int>(informationRatesNew.size()) != kNew)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew[kNew-1] != 1.0)
        throw std::invalid_argument("informationRatesNew must end with 1");
      infoRatesNew = informationRatesNew; // copy
    } else {
      for (int i = 0; i < kNew; ++i)
        infoRatesNew[i] = static_cast<double>(i+1) / static_cast<double>(kNew);
    }

    // effStoppingNew: default to all 1s if missing
    if (none_na(efficacyStoppingNew)) {
      if (static_cast<int>(efficacyStoppingNew.size()) != kNew)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew[kNew-1] != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew; // copy
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    // futStoppingNew: default to all 1s if missing
    if (none_na(futilityStoppingNew)) {
      if (static_cast<int>(futilityStoppingNew.size()) != kNew)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew[kNew-1] != 1)
        throw std::invalid_argument("futilityStoppingNew must end with 1");
      futStoppingNew = futilityStoppingNew; // copy
    } else {
      futStoppingNew.assign(kNew, 1);
    }

    if (!(asfNew == "of" || asfNew == "p" || asfNew == "wt" ||
        asfNew == "sfof" || asfNew == "sfp" ||
        asfNew == "sfkd" || asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }
    if ((asfNew == "wt" || asfNew == "sfkd" || asfNew == "sfhsd") &&
        std::isnan(parameterAlphaSpendingNew)) {
      throw std::invalid_argument("Missing value for parameterAlphaSpendingNew");
    }
    if (asfNew == "sfkd" && parameterAlphaSpendingNew <= 0.0) {
      throw std::invalid_argument (
          "parameterAlphaSpendingNew must be positive for sfKD");
    }

    if (std::isnan(INew)) {
      if (!(bsfNew == "sfof" || bsfNew == "sfp" ||
          bsfNew == "sfkd" || bsfNew == "sfhsd" || bsfNew == "user" ||
          bsfNew == "none")) {
        throw std::invalid_argument("Invalid value for typeBetaSpendingNew");
      }
    } else {
      if (!(bsfNew == "sfof" || bsfNew == "sfp" ||
          bsfNew == "sfkd" || bsfNew == "sfhsd" || bsfNew == "none")) {
        throw std::invalid_argument("Invalid value for typeBetaSpendingNew");
      }
    }

    if ((bsfNew == "sfkd" || bsfNew == "sfhsd") &&
        std::isnan(parameterBetaSpendingNew)) {
      throw std::invalid_argument("Missing value for parameterBetaSpendingNew");
    }
    if (bsfNew == "sfkd" && parameterBetaSpendingNew <= 0.0) {
      throw std::invalid_argument(
          "parameterBetaSpendingNew must be positive for sfKD");
    }

    if (std::isnan(INew) && bsfNew == "user") {
      if (!none_na(userBetaSpendingNew))
        throw std::invalid_argument("userBetaSpendingNew must be specified");
      if (static_cast<int>(userBetaSpendingNew.size()) != kNew)
        throw std::invalid_argument("Invalid length of userBetaSpendingNew");
      if (userBetaSpendingNew[0] < 0.0)
        throw std::invalid_argument("userBetaSpendingNew must be nonnegative");
      if (any_nonincreasing(userBetaSpendingNew))
        throw std::invalid_argument("userBetaSpendingNew must be nondecreasing");
      if (userBetaSpendingNew[kNew] != betaNew)
        throw std::invalid_argument(
            "userBetaSpendingNew must end with specified betaNew");
    }

    // spendingTimeNew: default to informationRatesNew if missing
    if (none_na(spendingTimeNew)) {
      if (static_cast<int>(spendingTimeNew.size()) != kNew)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew[kNew-1] != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = infoRatesNew;
    }
  }

  if (varianceRatio <= 0.0) {
    throw std::invalid_argument("varianceRatio must be positive");
  }
  // ----------- End of Input Validation ----------- //

  // obtain critical values for the primary trial
  std::vector<double> l(kMax, -6.0), zero(kMax, 0.0);
  std::vector<double> critValues = criticalValues;
  double alpha1 = alpha;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (kMax > 1 && static_cast<int>(criticalValues.size()) == kMax) {
      bool hasNaN = false;
      for (int i = 0; i < kMax - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[kMax-1])) {
        haybittle = true;
      }
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(kMax);
      for (int i = 0; i < kMax - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  } else {
    for (int i = 0; i < kMax; ++i) {
      if (!effStopping[i]) critValues[i] = 6.0;
    }
    ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    alpha1 = std::accumulate(v.begin(), v.end(), 0.0);
  }

  // obtain futility bounds for the primary trial
  std::vector<double> futBounds = futilityBounds;
  if (kMax > 1) {
    if (missingFutilityBounds && bsf == "none") {
      futBounds = std::vector<double>(kMax, -6.0);
      futBounds[kMax-1] = critValues[kMax-1];
    } else if (!missingFutilityBounds &&
      static_cast<int>(futBounds.size()) == kMax-1) {
      futBounds.push_back(critValues[kMax-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // multiplier for the boundaries under the alternative hypothesis
  std::vector<double> w(kMax, std::sqrt(varianceRatio));
  if (!none_na(futBounds)) {
    if (std::isnan(IMax)) {
      throw std::invalid_argument("IMax must be provided");
    }
    if (IMax <= 0.0) {
      throw std::invalid_argument("IMax must be positive");
    }

    std::vector<double> delta(kMax, theta);
    std::vector<double> information(kMax);
    for (int i = 0; i < kMax; ++i) {
      information[i] = IMax * infoRates[i];
    }

    ListCpp out = getPower(alpha1, kMax, critValues, delta, information, bsf,
                           parameterBetaSpending, spendTime, futStopping, w);
    futBounds = out.get<std::vector<double>>("futilityBounds");
  }

  // compute conditional alpha, conditional power, and predictive power
  int k1 = kMax - L;

  std::vector<double> t1(k1), r1(k1), b1(k1), a1(k1, -6.0), theta1(k1, 0.0);
  for (int l = 0; l < k1; ++l) {
    t1[l] = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
    r1[l] = infoRates[L - 1] / infoRates[l + L];
    b1[l] = (critValues[l + L] - std::sqrt(r1[l]) * zL) / std::sqrt(1.0 - r1[l]);
    if (!effStopping[l + L]) b1[l] = 6.0;
  }

  // conditional type I error
  ListCpp probs = exitprobcpp(b1, a1, theta1, t1);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  double alphaNew = std::accumulate(v.begin(), v.end(), 0.0);

  // conditional power and predictive power
  double conditionalPower, predictivePower;
  for (int l = 0; l < k1; ++l) {
    a1[l] = (futBounds[l + L] - std::sqrt(r1[l]) * zL) / std::sqrt(1.0 - r1[l]);
    if (!futStopping[l + L]) a1[l] = -6.0;
  }

  if (!std::isnan(IMax)) {
    double sigma = 1.0 / std::sqrt(IMax * infoRates[L - 1]);
    double mu = zL * sigma;
    std::vector<double> theta1(k1, mu);

    std::vector<double> I1(k1);
    for (int l = 0; l < k1; ++l) {
      I1[l] = IMax * (infoRates[l + L] - infoRates[L - 1]);
    }

    ListCpp probs = exitprobcpp(b1, a1, theta1, I1);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    conditionalPower = std::accumulate(v.begin(), v.end(), 0.0);

    // predictive power
    auto f = [&](double theta)->double {
      std::vector<double> theta1(k1, theta);
      ListCpp probs = exitprobcpp(b1, a1, theta1, I1);
      auto v = probs.get<std::vector<double>>("exitProbUpper");
      return std::accumulate(v.begin(), v.end(), 0.0);
    };

    double lower = mu - 6.0*sigma, upper = mu + 6.0*sigma;
    predictivePower = intnorm(f, mu, sigma, lower, upper);
  } else {
    conditionalPower = NaN;
    predictivePower = NaN;
  }

  ListCpp des1;
  des1.push_back(L, "L");
  des1.push_back(zL, "zL");
  des1.push_back(theta, "theta");
  des1.push_back(kMax, "kMax");
  des1.push_back(std::move(infoRates), "informationRates");
  des1.push_back(std::move(critValues), "efficacyBounds");
  des1.push_back(std::move(futBounds), "futilityBounds");
  des1.push_back(alphaNew, "conditionalAlpha");
  des1.push_back(conditionalPower, "conditionalPower");
  des1.push_back(predictivePower, "predictivePower");
  des1.push_back(MullerSchafer, "MullerSchafer");

  ListCpp des2;
  if (!MullerSchafer) {
    effStoppingNew = subset(effStopping, L, kMax);
    futStoppingNew = subset(futStopping, L, kMax);
    des2 = getDesigncpp(betaNew, INew, theta, k1, t1, effStoppingNew,
                        futStoppingNew, b1, NaN, typeAlphaSpendingNew,
                        parameterAlphaSpendingNew, {NaN}, a1,
                        typeBetaSpendingNew, parameterBetaSpendingNew,
                        userBetaSpendingNew, spendTimeNew, varianceRatio);
  } else {
    if (!std::isnan(betaNew) && betaNew >= 1.0 - alphaNew) {
      throw std::invalid_argument(
          "betaNew must be less than 1 minus the conditional type I error");
    }

    std::vector<double> b1New(kNew, NaN), a1New(kNew, NaN);
    des2 = getDesigncpp(betaNew, INew, theta, kNew, infoRatesNew, effStoppingNew,
                        futStoppingNew, b1New, alphaNew, typeAlphaSpendingNew,
                        parameterAlphaSpendingNew, {NaN}, a1New,
                        typeBetaSpendingNew, parameterBetaSpendingNew,
                        userBetaSpendingNew, spendTimeNew, varianceRatio);
  }

  ListCpp result;
  result.push_back(std::move(des1), "primaryTrial");
  result.push_back(std::move(des2), "secondaryTrial");
  return result;
}

//' @title Adaptive Design at an Interim Look
//' @description
//' Calculates the conditional power for specified incremental
//' information, given the interim results, parameter value,
//' data-dependent changes in the error spending function, and
//' the number and spacing of interim looks. Conversely,
//' calculates the incremental information required to attain
//' a specified conditional power, given the interim results,
//' parameter value, data-dependent changes in the error
//' spending function, and the number and spacing of interim looks.
//'
//' @param betaNew The type II error for the secondary trial.
//' @param INew The maximum information of the secondary trial. Either
//'   \code{betaNew} or \code{INew} should be provided, while the other
//'   must be missing.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of
//'   the primary trial.
//' @param theta The parameter value.
//' @param IMax The maximum information of the primary trial. Must be
//'   provided if \code{futilityBounds} is missing and
//'   \code{typeBetaSpending} is not equal to "none", or
//'   if conditional power calculation is desired.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param futilityStopping Indicators of whether futility stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping for the primary trial.
//' @param alpha The significance level of the primary trial.
//'   Defaults to 0.025.
//' @param typeAlphaSpending The type of alpha spending for the primary
//'   trial. One of the following:
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function,
//'   "user" for user defined spending, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to Delta for "WT", rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param userAlphaSpending The user-defined alpha spending for the
//'   primary trial. Represents the cumulative alpha spent up to each stage.
//' @param futilityBounds The lower boundaries on the z-test statistic scale
//'   for futility stopping for the primary trial. Defaults to
//'   \code{rep(-6, kMax-1)} if left unspecified.
//' @param typeBetaSpending The type of beta spending for the primary trial.
//'   One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early futility stopping.
//'   Defaults to "none".
//' @param parameterBetaSpending The parameter value of beta spending
//'   for the primary trial. Corresponds to rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param futilityStoppingNew The indicators of whether futility stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param typeAlphaSpendingNew The type of alpha spending for the secondary
//'   trial. One of the following:
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to Delta for "WT", rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param typeBetaSpendingNew The type of beta spending for the secondary
//'   trial. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function,
//'   "user" for user defined spending, and
//'   "none" for no early futility stopping.
//'   Defaults to "none".
//' @param parameterBetaSpendingNew The parameter value of beta spending
//'   for the secondary trial. Corresponds to rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param userBetaSpendingNew The user-defined cumulative beta spending.
//'   Represents the cumulative beta spent up to each stage of the
//'   secondary trial.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRatesNew}.
//' @param varianceRatio The ratio of the variance under H0 to the
//'   variance under H1.
//'
//' @return An \code{adaptDesign} object with two list components:
//'
//' * \code{primaryTrial}: A list of selected information for the primary
//'   trial, including \code{L}, \code{zL}, \code{theta}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
//'   and \code{MullerSchafer}.
//'
//' * \code{secondaryTrial}: A \code{design} object for the secondary trial.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Lu Chi, H. M. James Hung, and Sue-Jane Wang.
//' Modification of sample size in group sequential clinical trials.
//' Biometrics 1999;55:853-857.
//'
//' Hans-Helge Muller and Helmut Schafer.
//' Adaptive group sequential designs for clinical trials:
//' Combining the advantages of adaptive and of
//' classical group sequential approaches.
//' Biometrics 2001;57:886-891.
//'
//' @seealso \code{\link{getDesign}}
//'
//' @examples
//'
//' # original group sequential design with 90% power to detect delta = 6
//' delta = 6
//' sigma = 17
//' n = 282
//' (des1 = getDesign(IMax = n/(4*sigma^2), theta = delta, kMax = 3,
//'                   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'                   parameterAlphaSpending = -4))
//'
//' # interim look results
//' L = 1
//' n1 = n/3
//' delta1 = 4.5
//' sigma1 = 20
//' zL = delta1/sqrt(4/n1*sigma1^2)
//'
//' t = des1$byStageResults$informationRates
//'
//' # conditional power with sample size increase
//' (des2 = adaptDesign(
//'   betaNew = NA, INew = 420/(4*sigma1^2),
//'   L = L, zL = zL, theta = delta1,
//'   IMax = n/(4*sigma1^2), kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4))
//'
//' # Muller & Schafer (2001) method to design the secondary trial:
//' # 3-look gamma(-2) spending with 84% power at delta = 4.5 and sigma = 20
//' (des2 = adaptDesign(
//'   betaNew = 0.16, INew = NA,
//'   L = L, zL = zL, theta = delta1,
//'   IMax = n/(4*sigma1^2), kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpendingNew = "sfHSD",
//'   parameterAlphaSpendingNew = -2))
//'
//' # incremental sample size for sigma = 20
//' (nNew = 4*sigma1^2*des2$secondaryTrial$overallResults$information)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List adaptDesign(
    double betaNew = NA_REAL,
    double INew = NA_REAL,
    const int L = NA_INTEGER,
    const double zL = NA_REAL,
    const double theta = NA_REAL,
    const double IMax = NA_REAL,
    const int kMax = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& futilityBounds = NA_REAL,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int kNew = NA_INTEGER,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const std::string& typeBetaSpendingNew = "none",
    const double parameterBetaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& userBetaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL,
    const double varianceRatio = 1.0) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto futStoppingNew = convertLogicalVector(futilityStoppingNew);
  auto userBetaNew = Rcpp::as<std::vector<double>>(userBetaSpendingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  auto cpp_result = adaptDesigncpp(
    betaNew, INew, L, zL, theta, IMax, kMax, infoRates, effStopping,
    futStopping, critValues, alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlpha, futBounds, typeBetaSpending, parameterBetaSpending, spendTime,
    MullerSchafer, kNew, infoRatesNew, effStoppingNew, futStoppingNew,
    typeAlphaSpendingNew, parameterAlphaSpendingNew, typeBetaSpendingNew,
    parameterBetaSpendingNew, userBetaNew, spendTimeNew, varianceRatio
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "adaptDesign";
  return result;
}


// bygroup: group-by helper that builds lookup tables and combined indices
ListCpp bygroup(const DataFrameCpp& data,
                const std::vector<std::string>& variables) {
  int n = data.nrows();
  int p = variables.size();
  ListCpp result;
  std::vector<int> nlevels(p);

  // IntMatrix for indices (n rows, p cols), column-major storage
  IntMatrix indices(n, p);

  // Flattened lookup buffers and per-variable metadata
  struct VarLookupInfo {
    int type; // 0=int, 1=double, 2=bool, 3=string
    int offset;
  };
  std::vector<VarLookupInfo> var_info(p);

  std::vector<int> int_flat;
  std::vector<double> dbl_flat;
  std::vector<unsigned char> bool_flat;
  std::vector<std::string> str_flat;

  ListCpp lookups_per_variable; // will contain a std::vector for each variable

  for (int i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    if (!data.containElementNamed(var))
      throw std::invalid_argument("Data must contain variable: " + var);

    if (data.int_cols.count(var)) {
      const auto& col = data.int_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w); // indices 0..(levels-1)

      // append w to flat buffer and record metadata
      int off = int_flat.size();
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

      int off = dbl_flat.size();
      dbl_flat.insert(dbl_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{1, off};

      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.bool_cols.count(var)) {
      const auto& col = data.bool_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);

      int off = bool_flat.size();
      bool_flat.insert(bool_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{2, off};

      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.string_cols.count(var)) {
      const auto& col = data.string_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);

      int off = str_flat.size();
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
  for (int i = 0; i < p; ++i) orep *= nlevels[i];
  int lookup_nrows = orep;

  for (int i = 0; i < p; ++i) {
    orep /= nlevels[i];
    const int* col_ptr = indices.data_ptr() + i * n;
    for (int j = 0; j < n; ++j) {
      combined_index[j] += col_ptr[j] * orep;
    }
  }

  // Build lookup_df with columns repeated in the same pattern as before.
  DataFrameCpp lookup_df;
  int repeat_each = lookup_nrows;
  for (int i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    int nlevels_i = nlevels[i];
    repeat_each /= nlevels_i;
    int times = lookup_nrows / ( nlevels_i * repeat_each );

    VarLookupInfo info = var_info[i];
    if (info.type == 0) {
      const int* base = int_flat.data() + info.offset;
      std::vector<int> col(lookup_nrows);
      int idxw = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level)
          for (int r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    } else if (info.type == 1) {
      const double* base = dbl_flat.data() + info.offset;
      std::vector<double> col(lookup_nrows);
      int idxw = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level)
          for (int r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    } else if (info.type == 2) {
      const unsigned char* base = bool_flat.data() + info.offset;
      std::vector<unsigned char> col(lookup_nrows);
      int idxw = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) {
          for (int r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
        }
      }
      lookup_df.push_back(std::move(col), var);
    } else { // string
      const std::string* base = str_flat.data() + info.offset;
      std::vector<std::string> col(lookup_nrows);
      int idxw = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level)
          for (int r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
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
int cholesky2(FlatMatrix& matrix, int n, double toler) {
  double* base = matrix.data_ptr();
  double eps = 0.0;
  for (int i = 0; i < n; ++i) {
    double val = matrix(i, i);
    if (val > eps) eps = val;
  }
  if (eps == 0.0) eps = toler; else eps *= toler;
  int nonneg = 1;
  int rank = 0;

  for (int i = 0; i < n; ++i) {
    double* col_i = base + i * n;
    double pivot = col_i[i];
    if (std::isinf(pivot) || pivot < eps) {
      col_i[i] = 0.0;
      if (pivot < -8.0 * eps) nonneg = -1;
    } else {
      ++rank;
      for (int j = i + 1; j < n; ++j) {
        double* col_j = base + j * n;
        double temp = col_i[j] / pivot;
        col_i[j] = temp;
        col_j[j] -= temp * temp * pivot;
        for (int k = j + 1; k < n; ++k) {
          col_j[k] -= temp * col_i[k];
        }
      }
    }
  }
  return rank * nonneg;
}

// chsolve2 assumes matrix holds the representation produced by cholesky2
void chsolve2(FlatMatrix& matrix, int n, double* y) {
  // Forward substitution L * z = y
  double* base = matrix.data_ptr();
  for (int j = 0; j < n-1; ++j) {
    double yj = y[j];
    if (yj == 0.0) continue;
    double* col_j = base + j * n;
    for (int i = j + 1; i < n; ++i) {
      y[i] -= yj * col_j[i];
    }
  }
  // Now y holds z; solve L^T * x = z
  if (n == 0) return;
  for (int i = n - 1; i >= 0; --i) {
    double* col_i = base + i * n;
    double diag = col_i[i];
    if (diag == 0.0) {
      y[i] = 0.0;
    } else {
      double temp = y[i] / diag;
      for (int j = i + 1; j < n; ++j) temp -= y[j] * col_i[j];
      y[i] = temp;
    }
  }
}

// invsympd: returns the inverse of a symmetric positive definite matrix
FlatMatrix invsympd(const FlatMatrix& matrix, int n, double toler) {
  FlatMatrix v = matrix; // copy
  cholesky2(v, n, toler);
  FlatMatrix iv(n, n);
  for (int i = 0; i < n; ++i) {
    iv(i,i) = 1.0;
    double* ycol = iv.data_ptr() + i * n;
    chsolve2(v, n, ycol);
  }
  return iv;
}


// Transpose a FlatMatrix (double)
FlatMatrix transpose(const FlatMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return FlatMatrix();

  const int src_nrow = M.nrow;
  const int src_ncol = M.ncol;
  FlatMatrix out(src_ncol, src_nrow); // swapped dims

  const double* src = M.data_ptr();
  double* dst = out.data_ptr();

  for (int c = 0; c < src_ncol; ++c) {
    const double* src_col = src + c * src_nrow;
    for (int r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }

  return out;
}

// Transpose an IntMatrix (int)
IntMatrix transpose(const IntMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return IntMatrix();

  const int src_nrow = M.nrow;
  const int src_ncol = M.ncol;
  IntMatrix out(src_ncol, src_nrow); // swapped dims

  const int* src = M.data_ptr();
  int* dst = out.data_ptr();

  for (int c = 0; c < src_ncol; ++c) {
    const int* src_col = src + c * src_nrow;
    for (int r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }

  return out;
}

// Transpose a BoolMatrix (unsigned char)
BoolMatrix transpose(const BoolMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return BoolMatrix();

  const int src_nrow = M.nrow;
  const int src_ncol = M.ncol;
  BoolMatrix out(src_ncol, src_nrow); // swapped dims

  const unsigned char* src = M.data_ptr();
  unsigned char* dst = out.data_ptr();

  for (int c = 0; c < src_ncol; ++c) {
    const unsigned char* src_col = src + c * src_nrow;
    for (int r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }

  return out;
}


// Householder vector
// Given an n-vector x, this function computes an n-vector v with v(1) = 1
// such that (I - 2*v*t(v)/t(v)*v)*x is zero in all but the first component.
std::vector<double> house(const std::vector<double>& x) {
  int n = static_cast<int>(x.size());
  double sumxx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  double mu = std::sqrt(sumxx);
  std::vector<double> v = x;
  if (mu > 0.0) {
    double beta = x[0] + std::copysign(1.0, x[0])*mu;
    for (int i=1; i<n; ++i) {
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
void row_house(FlatMatrix& A, const int i1, const int i2, const int j1,
               const int j2, const std::vector<double>& v) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow) {
    throw std::invalid_argument("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol) {
    throw std::invalid_argument("Invalid column indices j1 and j2");
  }

  int m = i2-i1+1, n = j2-j1+1;
  double sumvv = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double beta = -2.0 / sumvv;
  std::vector<double> w(n);
  for (int j=0; j<n; ++j) {
    for (int i=0; i<m; ++i) {
      w[j] += A(i+i1,j+j1)*v[i];
    }
    w[j] *= beta;
  }

  for (int j=0; j<n; ++j) {
    for (int i=0; i<m; ++i) {
      A(i+i1,j+j1) += v[i]*w[j];
    }
  }
}

// Householder post-multiplication
// Given an m-by-n matrix A and a nonzero n-vector v with v(1) = 1,
// the following algorithm overwrites A with A*P where
// P = I - 2*v*t(v)/t(v)*v.
void col_house(FlatMatrix& A, const int i1, const int i2, const int j1,
               const int j2, const std::vector<double>& v) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow) {
    throw std::invalid_argument("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol) {
    throw std::invalid_argument("Invalid column indices j1 and j2");
  }

  int m = i2-i1+1, n = j2-j1+1;
  double sumvv = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double beta = -2.0 / sumvv;
  std::vector<double> w(m);
  for (int j=0; j<n; ++j) {
    for (int i=0; i<m; ++i) {
      w[i] += A(i+i1,j+j1)*v[j];
    }
  }
  for (int i=0; i<m; ++i) {
    w[i] *= beta;
  }
  for (int j=0; j<n; ++j) {
    for (int i=0; i<m; ++i) {
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
void row_rot(FlatMatrix& A, const int i1, const int i2, const int j1, const int j2,
             const double c, const double s) {
  if (i1 < 0 || i1 >= i2 || i2 >= A.nrow) {
    throw std::invalid_argument("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol) {
    throw std::invalid_argument("Invalid column indices j1 and j2");
  }

  int q = j2-j1+1;
  for (int j=0; j<q; ++j) {
    int jj1 = j+j1;
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
void col_rot(FlatMatrix& A, const int i1, const int i2, const int j1, const int j2,
             const double c, const double s) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow) {
    throw std::invalid_argument("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 >= j2 || j2 >= A.ncol) {
    throw std::invalid_argument("Invalid column indices j1 and j2");
  }

  int q = i2-i1+1;
  for (int i=0; i<q; ++i) {
    int ii1 = i+i1;
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
  int m = A.nrow, n = A.ncol;
  if (m < n) {
    throw std::invalid_argument("The input matrix must have # rows >= # columns");
  }

  double tol = 1e-12;
  bool bidiag = true;
  for (int j=2; j<n; ++j) {
    for (int i=0; i<j-1; ++i) {
      if (std::fabs(A(i,j)) > tol) {
        bidiag = false; break;
      }
    }
  }
  for (int j=0; j<n-1; ++j) {
    for (int i=j+1; i<n; ++i) {
      if (std::fabs(A(i,j)) > tol) {
        bidiag = false; break;
      }
    }
  }
  for (int j=0; j<n; ++j) {
    for (int i=n; i<m-1; ++i) {
      if (std::fabs(A(i,j)) > tol) {
        bidiag = false; break;
      }
    }
  }

  FlatMatrix B(n,n);
  FlatMatrix U(m,m), V(n,n);
  for (int i=0; i<m; ++i) U(i,i) = 1.0;
  for (int i=0; i<n; ++i) V(i,i) = 1.0;

  if (bidiag) {
    B = A;
  } else {
    for (int j=0; j<n; ++j) {
      std::vector<double> v(m-j);
      for (int i=0; i<m-j; ++i) {
        v[i] = A(i+j,j);
      }
      v = house(v);

      row_house(A, j, m-1, j, n-1, v);

      // update the sub-diagonal elements of column j
      for (int i=1; i<m-j; ++i) {
        A(i+j,j) = v[i];
      }

      if (j < n-2) {
        std::vector<double> v(n-j-1);
        for (int i=0; i<n-j-1; ++i) {
          v[i] = A(j,i+j+1);
        }
        v = house(v);

        col_house(A, j, m-1, j+1, n-1, v);

        // update the elements of row j
        for (int i=1; i<n-j-1; ++i) {
          A(j,i+j+1) = v[i];
        }
      }
    }

    if (outtransform) {
      for (int j=n-1; j>=0; --j) {
        std::vector<double> v(m-j);
        v[0] = 1.0;
        for (int i=1; i<m-j; ++i) {
          v[i] = A(i+j,j);
        }

        row_house(U, j, m-1, j, m-1, v);
      }

      for (int j=n-3; j>=0; --j) {
        std::vector<double> v(n-j-1);
        v[0] = 1.0;
        for (int i=1; i<n-j-1; ++i) {
          v[i] = A(j,i+j+1);
        }

        row_house(V, j+1, n-1, j+1, n-1, v);
      }
    }

    for (int j=0; j<n; ++j) {
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
ListCpp zero_diagonal(FlatMatrix& B, const int k, const bool outtransform = true) {
  int n = B.nrow;
  if (B.ncol != n) {
    throw std::invalid_argument("The input matrix must be a square matrix");
  }
  if (k < 0 || k >= n-1) {
    throw std::invalid_argument("Invalid value for index k");
  }
  FlatMatrix U(n,n);
  for (int i=0; i<n; ++i) U(i,i) = 1.0;

  for (int j=k+1; j<n; ++j) {
    std::vector<double> v = givens(B(k,j), B(j,j));
    double w = v[0];
    v[0] = -v[1]; v[1] = w;
    int j1 = j < n-1 ? j+1 : n-1;
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
  int n = B.ncol;
  FlatMatrix U(n,n), V(n,n);
  for (int i=0; i<n; ++i) {
    U(i,i) = 1.0; V(i,i) = 1.0;
  }

  double f1 = B(n-3,n-2), f2 = B(n-2,n-1);
  double d1 = B(n-2,n-2), d2 = B(n-1,n-1);
  double a1 = f1*f1 + d1*d1, a2 = f2*f2 + d2*d2, b1 = f2*d1;
  double d = 0.5*(a1-a2);
  double mu = a2 + d - std::copysign(1.0, d)*std::sqrt(d*d + b1*b1);
  double y = B(0,0)*B(0,0) - mu;
  double z = B(0,0)*B(0,1);
  std::vector<double> v(2);
  for (int k=0; k<n-1; ++k) {
    v = givens(y,z);
    int k1 = k > 0 ? k-1 : 0;
    col_rot(B, k1, k+1, k, k+1, v[0], v[1]);
    if (outtransform) col_rot(V, 0, k+1, k, k+1, v[0], v[1]);

    y = B(k,k);
    z = B(k+1,k);
    v = givens(y,z);
    int k2 = k < n-2 ? k+2 : n-1;
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
  int m1 = X.nrow, n1 = X.ncol, m, n;
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
  for (int i=0; i<m; ++i) U(i,i) = 1.0;
  for (int i=0; i<n; ++i) V(i,i) = 1.0;
  if (outtransform) {
    U = a.get<FlatMatrix>("U");
    V = a.get<FlatMatrix>("V");
  }

  double tol = 1e-12;
  int p, q = 0;
  while (q < n) {
    for (int i=1; i<n; ++i) {
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
    for (int i=n-1; i>=1; --i) {
      if (B(i-1,i) != 0.0) {
        q = n-i-1; break;
      }
    }

    p = 0;
    for (int i=n-q-2; i>=1; --i) {
      if (B(i-1,i) == 0.0) {
        p = i; break;
      }
    }

    if (q < n) {
      // if any diagonal entry in B22 is zero, then zero the superdiagonal
      // entry in the same row
      FlatMatrix B22 = subset_flatmatrix(B, p, n-q, p, n-q);
      for (int i=0; i<n-p-q-1; ++i) {
        if (std::fabs(B22(i,i)) < tol) {
          ListCpp b = zero_diagonal(B22, i, outtransform);
          if (outtransform) {
            FlatMatrix Z = b.get<FlatMatrix>("U");
            FlatMatrix W = subset_flatmatrix(U, 0, m, p, n-q);
            for (int k=0; k<n-p-q; ++k) {
              for (int j=0; j<m; ++j) {
                U(j,k+p) = 0.0;
              }
            }
            for (int k=0; k<n-p-q; ++k) {
              for (int l=0; l<n-p-q; ++l) {
                for (int j=0; j<m; ++j) {
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
      for (int j=0; j<n-p-q; ++j) {
        for (int i=0; i<n-p-q; ++i) {
          B(i+p,j+p) = B22(i,j);
        }
      }

      if (outtransform) {
        FlatMatrix Z1 = c.get<FlatMatrix>("U");
        FlatMatrix W1 = subset_flatmatrix(U, 0, m, p, n-q);
        for (int j=0; j<n-p-q; ++j) {
          for (int i=0; i<m; ++i) {
            U(i,j+p) = 0.0;
          }
        }
        for (int j=0; j<n-p-q; ++j) {
          for (int k=0; k<n-p-q; ++k) {
            for (int i=0; i<m; ++i) {
              U(i,j+p) += W1(i,k)*Z1(k,j);
            }
          }
        }

        FlatMatrix Z2 = c.get<FlatMatrix>("V");
        FlatMatrix W2 = subset_flatmatrix(V, 0, n, p, n-q);
        for (int j=0; j<n-p-q; ++j) {
          for (int i=0; i<n; ++i) {
            V(i,j+p) = 0.0;
          }
        }
        for (int j=0; j<n-p-q; ++j) {
          for (int k=0; k<n-p-q; ++k) {
            for (int i=0; i<n; ++i) {
              V(i,j+p) += W2(i,k)*Z2(k,j);
            }
          }
        }
      }
    }
  }

  std::vector<double> d(n);
  for (int i=0; i<n; ++i) d[i] = B(i,i);

  // ensure the singular values are positive
  for (int i=0; i<n; ++i) {
    if (d[i] < 0.0) {
      d[i] = -d[i];
      for (int j=0; j<n; ++j) V(j,i) = -V(j,i);
    }
  }

  if (decreasing) {
    // order the singular values from the largest to the smallest
    // and the arrange the associated vectors accordingly
    std::vector<int> order = seqcpp(0, n-1);
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return d[i] > d[j];
    });

    subset_in_place(d, order);
    if (outtransform) {
      FlatMatrix Z = U;
      FlatMatrix W = V;
      for (int i=0; i<n; ++i) {
        int k = order[i];
        for (int j=0; j<m; ++j) U(j,i) = Z(j,k);
        for (int j=0; j<n; ++j) V(j,i) = W(j,k);
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
                                      const double tol=0.000001) {
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
