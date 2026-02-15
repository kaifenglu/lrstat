#include "utilities.h"
#include "dataframe_list.h"
#include "miettinen_nurminen.h"

#include <algorithm>     // any_of, distance, fill, min, min_element
#include <cctype>        // tolower
#include <cmath>         // fabs, isnan
#include <cstddef>       // size_t
#include <limits>        // numeric_limits
#include <numeric>       // accumulate
#include <stdexcept>     // invalid_argument
#include <string>        // string
#include <vector>        // vector
#include <unordered_map> // unordered_map
#include <utility>       // make_pair, pair

#include <Rcpp.h>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/random.hpp>


struct OneSampleExactResult {
  double attainedAlpha;
  double power;
  int n;
  int r;
};

// Exact one-proportion power/alpha
OneSampleExactResult powerOnePropExactcpp(
    const int n,
    const double piH0,
    const double pi,
    const double alpha) {

  using boost::math::binomial_distribution;

  OneSampleExactResult out;
  out.n = n;
  out.r = 0;
  out.attainedAlpha = 0.0;
  out.power = 0.0;

  bool directionUpper = (pi > piH0);

  if (directionUpper) {
    // Upper-tail test: reject for X >= r
    binomial_distribution<double> d0(n, piH0);
    // quantile gives smallest k with P[X <= k] >= p
    double qd = boost::math::quantile(d0, 1.0 - alpha);
    int q = static_cast<int>(std::floor(qd + 1e-12));
    int r = q + 1;
    out.r = r;

    // attained alpha = P_H0[X >= r] = 1 - P_H0[X <= r-1]
    out.attainedAlpha = 1.0 - boost::math::cdf(d0, r - 1);

    binomial_distribution<double> d1(n, pi);
    out.power = 1.0 - boost::math::cdf(d1, r - 1);
  } else {
    // Lower-tail test: reject for X <= r
    binomial_distribution<double> d0(n, piH0);
    double qd = boost::math::quantile(d0, alpha);
    int rstar = static_cast<int>(std::floor(qd + 1e-12));
    double c_at_rstar = boost::math::cdf(d0, rstar);
    int r = (c_at_rstar <= alpha) ? rstar : (rstar - 1);
    out.r = r;

    // attained alpha = P_H0[X <= r]
    out.attainedAlpha = boost::math::cdf(d0, r);

    binomial_distribution<double> d1(n, pi);
    out.power = boost::math::cdf(d1, r);
  }

  return out;
}


//' @title Power for Binomial One-Sample Exact Test
//' @description Obtains the power for binomial one-sample exact test.
//'
//' @param n The sample size.
//' @param piH0 The response probability under the null hypothesis.
//' @param pi The response probability under the alternative hypothesis.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//'
//' @return A data frame containing the critical value of the number of
//' responses for rejecting the null hypothesis, the attained type I
//' error, the power for the exact test, the sample size, the
//' response probabilities under the null and alternative hypotheses,
//' and the direction of the alternative.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' powerOnePropExact(n = 110, piH0 = 0.15, pi = 0.25, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame powerOnePropExact(
    const int n = NA_INTEGER,
    const double piH0 = NA_REAL,
    const double pi = NA_REAL,
    const double alpha = 0.025) {

  auto out = powerOnePropExactcpp(n, piH0, pi, alpha);

  DataFrameCpp df;
  df.push_back(alpha, "alpha");
  df.push_back(out.attainedAlpha, "attainedAlpha");
  df.push_back(out.power, "power");
  df.push_back(out.n, "n");
  df.push_back(piH0, "piH0");
  df.push_back(pi, "pi");
  df.push_back(out.r, "r");

  return Rcpp::wrap(df);
}



OneSampleExactResult samplesizeOnePropExactcpp(
    const double beta,
    const double piH0,
    const double pi,
    const double alpha,
    const int max_n_search,
    const int window) {

  // Normal quantiles
  double z0 = boost_qnorm(1.0 - alpha);
  double z1 = boost_qnorm(1.0 - beta);
  double tau = sq((z0 + z1) / (pi - piH0));

  double v0 = piH0 * (1.0 - piH0);
  double v1 = pi * (1.0 - pi);
  double vmin = std::min(v0, v1);

  double n0 = tau * vmin;

  if (n0 > 0.5 * max_n_search) throw std::invalid_argument(
      std::string("Initial sample size estimate (") +
        std::to_string(static_cast<int>(std::floor(n0))) +
        ") is too large for exact test; consider normal approximation.");

  // Choose starting integer n
  int start_n = static_cast<int>(std::floor(n0));
  if (start_n < 1) start_n = 1;

  const double target = 1.0 - beta;

  // cache power evaluations
  std::unordered_map<int, double> power_cache;
  power_cache.reserve(1024);

  auto power_at = [&](int n)->double {
    if (n <= 0) return 0.0;
    auto it = power_cache.find(n);
    if (it != power_cache.end()) return it->second;
    auto res = powerOnePropExactcpp(n, piH0, pi, alpha);
    double p = res.power;
    power_cache.emplace(n, p);
    return p;
  };

  // Step-down: if start_n's power >= target, decrease until power < target
  int n = start_n;
  double p_n = power_at(n);
  while (p_n >= target && n > 1) {
    --n;
    p_n = power_at(n);
  }

  // Linear upward scan from candidate to max_n_search, enforcing window criterion
  int found_n = -1;
  for (int nn = n + 1; nn <= max_n_search; ++nn) {
    double p_nn = power_at(nn);
    if (p_nn < target) continue;

    // verify next window-1 sizes are also >= target
    bool ok = true;
    int fail_at = -1;
    for (int j = 1; j < window; ++j) {
      double pj = power_at(nn + j);
      if (pj < target) { ok = false; fail_at = j; break; }
    }
    if (ok) { found_n = nn; break; }
    n += fail_at; // skip ahead by fail_at for next iteration
  }

  if (found_n == -1) {
    throw std::runtime_error("No sample size <= max_n_search satisfies the "
                               "windowed power criterion");
  }

  auto out = powerOnePropExactcpp(found_n, piH0, pi, alpha);
  return out;
}


//' @title Sample Size for Binomial One-Sample Exact Test
//' @description Obtains the sample size for binomial one-sample exact test.
//'
//' @param beta The type II error.
//' @param piH0 The response probability under the null hypothesis.
//' @param pi The response probability under the alternative hypothesis.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//' @param max_n_search The maximum sample size to search up to. If no
//'   sample size up to this value satisfies the windowed power criterion,
//'   an error is thrown.
//' @param window The number of consecutive sample sizes that must all
//'   satisfy the power criterion to confirm the found sample size.
//'   This is to mitigate non-monotonicity of power in sample size for
//'   the exact test.
//'
//' @return A data frame containing the critical value of the number of
//' responses for rejecting the null hypothesis, the attained type I
//' error, the power for the exact test, the sample size, the
//' response probabilities under the null and alternative hypotheses,
//' and the direction of the alternative.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' samplesizeOnePropExact(beta = 0.2, piH0 = 0.15, pi = 0.25, alpha = 0.025)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame samplesizeOnePropExact(
    const double beta = 0.2,
    const double piH0 = NA_REAL,
    const double pi = NA_REAL,
    const double alpha = 0.025,
    const int max_n_search = 1000,
    const int window = 10) {

  auto out = samplesizeOnePropExactcpp(
    beta, piH0, pi, alpha, max_n_search, window);

  DataFrameCpp df;
  df.push_back(alpha, "alpha");
  df.push_back(out.attainedAlpha, "attainedAlpha");
  df.push_back(out.power, "power");
  df.push_back(out.n, "n");
  df.push_back(piH0, "piH0");
  df.push_back(pi, "pi");
  df.push_back(out.r, "r");

  return Rcpp::wrap(df);
}



OneSampleExactResult powerOneRateExactcpp(
    const int n,
    const double lambdaH0,
    const double lambda,
    const double D,
    const double alpha) {

  using boost::math::poisson_distribution;

  OneSampleExactResult out;
  out.n = n;
  out.r = 0;
  out.attainedAlpha = 0.0;
  out.power = 0.0;

  // Poisson mean under H0
  const double mean0 = static_cast<double>(n) * lambdaH0 * D;
  const double mean1 = static_cast<double>(n) * lambda * D;

  bool directionUpper = (lambda > lambdaH0);

  if (directionUpper) {
    // Upper-tail: reject for X >= r
    poisson_distribution<double> d0(mean0);
    double qd = boost::math::quantile(d0, 1.0 - alpha);
    int q = static_cast<int>(std::floor(qd + 1e-12));
    int r = q + 1;
    out.r = r;

    // attainedAlpha = P_H0[X >= r] = 1 - P_H0[X <= r-1]
    out.attainedAlpha = 1.0 - boost::math::cdf(d0, r - 1);

    poisson_distribution<double> d1(mean1);
    out.power = 1.0 - boost::math::cdf(d1, r - 1);
  } else {
    // Lower-tail: reject for X <= r
    poisson_distribution<double> d0(mean0);
    double qd = boost::math::quantile(d0, alpha);
    int rstar = static_cast<int>(std::floor(qd + 1e-12));
    double c_at_rstar = boost::math::cdf(d0, rstar);
    int r = (c_at_rstar <= alpha) ? rstar : (rstar - 1);

    out.attainedAlpha = boost::math::cdf(d0, r);

    poisson_distribution<double> d1(mean1);
    out.power = boost::math::cdf(d1, r);
  }

  return out;
}


//' @title Power for Poisson One-Sample Exact Test
//' @description Obtains the power for Poisson one-sample exact test.
//'
//' @param n The sample size.
//' @param lambdaH0 The Poisson rate under the null hypothesis.
//' @param lambda The Poisson rate under the alternative hypothesis.
//' @param D The average exposure per subject.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//'
//' @return A data frame containing the critical value of the number of
//' events for rejecting the null hypothesis, the attained type I
//' error, the power for the exact test, the sample size,
//' the Poisson rates under the null and alternative hypotheses,
//' the average exposure, and the direction of the alternative.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' powerOneRateExact(n = 525, lambdaH0 = 0.049, lambda = 0.012,
//'                   D = 0.5, alpha = 0.025)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame powerOneRateExact(
    const int n = NA_INTEGER,
    const double lambdaH0 = NA_REAL,
    const double lambda = NA_REAL,
    const double D = 1,
    const double alpha = 0.025) {

  auto out = powerOneRateExactcpp(n, lambdaH0, lambda, D, alpha);

  DataFrameCpp df;
  df.push_back(alpha, "alpha");
  df.push_back(out.attainedAlpha, "attainedAlpha");
  df.push_back(out.power, "power");
  df.push_back(out.n, "n");
  df.push_back(lambdaH0, "lambdaH0");
  df.push_back(lambda, "lambda");
  df.push_back(D, "D");
  df.push_back(out.r, "r");

  return Rcpp::wrap(df);
}



OneSampleExactResult samplesizeOneRateExactcpp(
    const double beta,
    const double lambdaH0,
    const double lambda,
    const double D,
    const double alpha,
    const int max_n_search,
    const int window) {

  // Normal quantiles
  double z0 = boost_qnorm(1.0 - alpha);
  double z1 = boost_qnorm(1.0 - beta);
  double tau = sq((z0 + z1) / log(lambda / lambdaH0));

  // per-unit variance of log(lambda_hat) under H0
  double v0 = 1.0 / (lambdaH0 * D);
  double v1 = 1.0 / (lambda * D);
  double vmin = std::min(v0, v1);

  double n0 = tau * vmin;

  if (n0 > 0.5 * max_n_search) throw std::invalid_argument(
      std::string("Initial sample size estimate (") +
        std::to_string(static_cast<int>(std::floor(n0))) +
        ") is too large for exact test; consider normal approximation.");

  // Choose starting integer n
  int start_n = static_cast<int>(std::floor(n0));
  if (start_n < 1) start_n = 1;

  const double target = 1.0 - beta;

  // cache power evaluations
  std::unordered_map<int, double> power_cache;
  power_cache.reserve(1024);

  auto power_at = [&](int n)->double {
    if (n <= 0) return 0.0;
    auto it = power_cache.find(n);
    if (it != power_cache.end()) return it->second;
    auto res = powerOneRateExactcpp(n, lambdaH0, lambda, D, alpha);
    double p = res.power;
    power_cache.emplace(n, p);
    return p;
  };

  // Step-down: if start_n's power >= target, decrease until power < target
  int n = start_n;
  double p_n = power_at(n);
  while (p_n >= target && n > 1) {
    --n;
    p_n = power_at(n);
  }

  // Linear upward scan from candidate to max_n_search, enforcing window criterion
  int found_n = -1;
  for (int nn = n + 1; nn <= max_n_search; ++nn) {
    double p_nn = power_at(nn);
    if (p_nn < target) continue;

    // verify next window-1 sizes are also >= target
    bool ok = true;
    int fail_at = -1;
    for (int j = 1; j < window; ++j) {
      double pj = power_at(nn + j);
      if (pj < target) { ok = false; fail_at = j; break; }
    }
    if (ok) { found_n = nn; break; }
    n += fail_at; // skip ahead by fail_at for next iteration
  }

  if (found_n == -1) {
    throw std::runtime_error("No sample size <= max_n_search satisfies the "
                               "windowed power criterion");
  }

  auto out = powerOneRateExactcpp(found_n, lambdaH0, lambda, D, alpha);
  return out;
}


//' @title Sample Size for Poisson One-Sample Exact Test
//' @description Obtains the sample size for Poisson one-sample exact test.
//'
//' @param beta The type II error.
//' @param lambdaH0 The Poisson rate under the null hypothesis.
//' @param lambda The Poisson rate under the alternative hypothesis.
//' @param D The average exposure per subject.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//' @param max_n_search The maximum sample size to search up to. If no
//'   sample size up to this value satisfies the windowed power criterion,
//'   an error is thrown.
//' @param window The number of consecutive sample sizes that must all
//'   satisfy the power criterion to confirm the found sample size.
//'   This is to mitigate non-monotonicity of power in sample size for
//'   the exact test.
//'
//' @return A data frame containing the critical value of the number of
//' events for rejecting the null hypothesis, the attained type I
//' error, the power for the exact test, the sample size,
//' the Poisson rates under the null and alternative hypotheses,
//' the average exposure, and the direction of the alternative.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' samplesizeOneRateExact(beta = 0.2, lambdaH0 = 0.2, lambda = 0.3,
//'                        D = 1, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame samplesizeOneRateExact(
    const double beta = 0.2,
    const double lambdaH0 = NA_REAL,
    const double lambda = NA_REAL,
    const double D = 1,
    const double alpha = 0.025,
    const int max_n_search = 1000,
    const int window = 10) {

  auto out = samplesizeOneRateExactcpp(
    beta, lambdaH0, lambda, D, alpha, max_n_search, window);

  DataFrameCpp df;
  df.push_back(alpha, "alpha");
  df.push_back(out.attainedAlpha, "attainedAlpha");
  df.push_back(out.power, "power");
  df.push_back(out.n, "n");
  df.push_back(lambdaH0, "lambdaH0");
  df.push_back(lambda, "lambda");
  df.push_back(D, "D");
  df.push_back(out.r, "r");

  return Rcpp::wrap(df);
}



// log choose(n,k) using lgamma for stability
double log_choose(int n, int k) {
  if (k < 0 || k > n) return -std::numeric_limits<double>::infinity();
  return std::lgamma(static_cast<double>(n) + 1.0)
    - std::lgamma(static_cast<double>(k) + 1.0)
    - std::lgamma(static_cast<double>(n - k) + 1.0);
}

// hypergeometric pmf: P(X = x) where X ~ Hypergeom(n1, n2, m) but expressed
// via choose(n1, x) * choose(n2, m-x) / choose(n1+n2, m)
double hypergeom_pmf(int x, int n1, int n2, int m) {
  // bounds check
  if (x < 0 || x > n1) return 0.0;
  int y = m - x;
  if (y < 0 || y > n2) return 0.0;
  double ln = log_choose(n1, x) + log_choose(n2, y) - log_choose(n1 + n2, m);
  return std::exp(ln);
}

// binomial pmf: choose(n,k) p^k (1-p)^(n-k)
double binomial_pmf(int k, int n, double p) {
  if (k < 0 || k > n) return 0.0;
  if (p <= 0.0) return (k == 0) ? 1.0 : 0.0;
  if (p >= 1.0) return (k == n) ? 1.0 : 0.0;
  double ln = log_choose(n, k) + k * std::log(p) + (n - k) * std::log(1.0 - p);
  return std::exp(ln);
}


// computes Fisher-exact power for two proportions (conditional test)
double powerFisherExactcpp(
    const int n,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha) {

  // determine sample split: n1 = round(n * r), n2 = n - n1
  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);
  int n1 = static_cast<int>(std::round(n * r));
  int n2 = n - n1;

  double power = 0.0; // unconditional power

  // iterate over all possible column sums m = total successes
  // for each m compute conditional distribution of x ~ Hypergeom(n1, n2, m)
  for (int m = 0; m <= n; ++m) { // compute conditional power given m
    // feasible x range given margins
    int lower = std::max(0, m - n2);
    int upper = std::min(n1, m);
    int k = upper - lower + 1;

    // compute p0 for x = lower..upper (hypergeometric PMF under H0)
    std::vector<double> p0(k);
    std::vector<int> x_vec(k);
    for (int x = lower; x <= upper; ++x) {
      int i = x - lower;
      p0[i] = hypergeom_pmf(x, n1, n2, m);
      x_vec[i] = x;
    }

    // produce order (indices) that sort p0 ascending
    std::vector<int> order = seqcpp(0, k - 1);
    std::sort(order.begin(), order.end(),
              [&p0](int i, int j) { return p0[i] < p0[j]; });

    // build sorted p0 using order
    std::vector<double> p0sorted = subset(p0, order);

    // group equal p0sorted values
    // collect indices of group boundaries
    std::vector<int> idxs;
    idxs.reserve(k + 1);
    idxs.push_back(0);
    for (int i = 1; i < k; ++i) {
      if (p0sorted[i] != p0sorted[i - 1]) idxs.push_back(i);
    }
    int k1 = static_cast<int>(idxs.size());
    idxs.push_back(k); // sentinel for end

    // compute cumulative probabilities for groups
    std::vector<double> cp0(k1, 0.0);
    double s = 0.0;
    for (int gi = 0; gi < k1; ++gi) {
      int start = idxs[gi];
      int end = idxs[gi + 1]; // exclusive
      for (int j = start; j < end; ++j) s += p0sorted[j];
      cp0[gi] = s;
    }

    // reject only if the smallest table-probability cumulative is <= alpha
    if (cp0[0] > alpha) continue;

    // find group index i where cp0[i] > alpha (first group exceeding alpha)
    int i_grp = std::upper_bound(cp0.begin() + 1, cp0.end(), alpha) - cp0.begin();
    // i_grp is the index of first cp0 > alpha, or k1 if none

    // compute unconditional alternative probabilities for corresponding x:
    // p1_x = Binom(n1, pi1)[x] * Binom(n2, pi2)[m-x]
    // then sort p1 according to the same order and sum up to idxs[i_grp]-1

    // number of smallest-table-probability cells included in rejection region
    for (int i = 0; i < i_grp; ++i) {
      for (int j = idxs[i]; j < idxs[i + 1]; ++j) {
        int orig_idx = order[j];
        int x = x_vec[orig_idx];
        double p_x = binomial_pmf(x, n1, pi1) * binomial_pmf(m - x, n2, pi2);
        power += p_x;
      }
    }
  } // end m loop

  return power;
}


//' @title Power for Fisher's Exact Test for Two Proportions
//' @description Obtains the power given sample size for Fisher's exact
//' test for two proportions.
//'
//' @param n The total sample size.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The two-sided significance level. Defaults to 0.05.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The two-sided significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @keywords internal
//'
//' @examples
//'
//' (design1 <- powerFisherExact(
//'   n = 136, pi1 = 0.25, pi2 = 0.05, alpha = 0.05))
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame powerFisherExact(
    const int n = NA_INTEGER,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05) {

  double power = powerFisherExactcpp(
    n, pi1, pi2, allocationRatioPlanned, alpha);

  DataFrameCpp df;
  df.push_back(alpha, "alpha");
  df.push_back(power, "power");
  df.push_back(n, "n");
  df.push_back(pi1, "pi1");
  df.push_back(pi2, "pi2");
  df.push_back(allocationRatioPlanned, "allocationRatioPlanned");

  return df;
}


DataFrameCpp samplesizeFisherExactcpp(
    const double beta,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha,
    const int max_n_search,
    const int window) {

  // allocation ratio r for group1 (fraction of total sample)
  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);
  double v1 = 1.0 / (4.0 * r * (1.0 - r));

  // normal approx estimate n0
  double z0 = boost_qnorm(1.0 - alpha);
  double z1  = boost_qnorm(1.0 - beta);

  double delta0 = std::asin(std::sqrt(pi1)) - std::asin(std::sqrt(pi2));
  if (delta0 == 0.0) throw std::invalid_argument(
    "pi1 and pi2 are equal; sample size is not finite");

  double n0 = sq(z0 + z1) * v1 / (delta0 * delta0);

  if (n0 > 0.5 * max_n_search) throw std::invalid_argument(
      std::string("Initial sample size estimate (") +
        std::to_string(static_cast<int>(std::floor(n0))) +
        ") is too large for exact test; consider normal approximation.");

  // Choose starting integer n
  int start_n = static_cast<int>(std::floor(n0));
  if (start_n < 1) start_n = 1;

  const double target = 1.0 - beta;

  // caching of power evaluations
  std::unordered_map<int, double> power_cache;
  power_cache.reserve(4096);

  auto power_at = [&](int n)->double {
    if (n <= 0) return 0.0;
    auto it = power_cache.find(n);
    if (it != power_cache.end()) return it->second;
    double p = powerFisherExactcpp(n, pi1, pi2, allocationRatioPlanned, alpha);
    power_cache.emplace(n, p);
    return p;
  };

  // Step-down: if start_n's power >= target, decrease until power < target
  int n = start_n;
  double p_n = power_at(n);
  while (p_n >= target && n > 1) {
    --n;
    p_n = power_at(n);
  }

  // Linear upward scan from candidate to max_n_search, enforcing window criterion
  int found_n = -1;
  for (int nn = n + 1; nn <= max_n_search; ++nn) {
    double p_nn = power_at(nn);
    if (p_nn < target) continue;

    // verify next window-1 sizes are also >= target
    bool ok = true;
    int fail_at = -1;
    for (int j = 1; j < window; ++j) {
      double pj = power_at(nn + j);
      if (pj < target) { ok = false; fail_at = j; break; }
    }
    if (ok) { found_n = nn; break; }
    n += fail_at; // skip ahead by fail_at for next iteration
  }

  if (found_n == -1) {
    throw std::runtime_error("No sample size <= max_n_search satisfies the "
                               "windowed power criterion");
  }

  // construct DataFrameCpp result using cached power
  double best_power = power_at(found_n);

  DataFrameCpp df;
  df.push_back(alpha, "alpha");
  df.push_back(best_power, "power");
  df.push_back(found_n, "n");
  df.push_back(pi1, "pi1");
  df.push_back(pi2, "pi2");
  df.push_back(allocationRatioPlanned, "allocationRatioPlanned");
  return df;
}


//' @title Sample Size for Fisher's Exact Test for Two Proportions
//' @description Obtains the sample size given power for Fisher's exact
//' test for two proportions.
//'
//' @param beta The type II error.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The two-sided significance level. Defaults to 0.05.
//' @param max_n_search The maximum sample size to search up to. If no
//'   sample size up to this value satisfies the windowed power criterion,
//'   an error is thrown.
//' @param window The number of consecutive sample sizes that must all
//'   satisfy the power criterion to confirm the found sample size.
//'   This is to mitigate non-monotonicity of power in sample size for
//'   the exact test.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The two-sided significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @keywords internal
//'
//' @examples
//'
//' (design1 <- samplesizeFisherExact(
//'   beta = 0.1, pi1 = 0.25, pi2 = 0.05, alpha = 0.05))
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame samplesizeFisherExact(
    const double beta = NA_REAL,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05,
    const int max_n_search = 1000,
    const int window = 10) {

  auto df = samplesizeFisherExactcpp(
    beta, pi1, pi2, allocationRatioPlanned, alpha, max_n_search, window);

  return Rcpp::wrap(df);
}


// -------------------- small helpers --------------------

std::vector<double> make_log_choose(int n) {
  std::vector<double> lg(n + 1);
  double lnNfact = std::lgamma(static_cast<double>(n) + 1.0);
  for (int k = 0; k <= n; ++k) {
    lg[k] = lnNfact - std::lgamma(static_cast<double>(k) + 1.0)
    - std::lgamma(static_cast<double>(n - k) + 1.0);
  }
  return lg;
}

double binomial_pmf_from_logchoose(
    int k, int n, double p, const std::vector<double>& logchoose) {
  if (k < 0 || k > n) return 0.0;
  if (p <= 0.0) return (k == 0) ? 1.0 : 0.0;
  if (p >= 1.0) return (k == n) ? 1.0 : 0.0;
  return std::exp(logchoose[k] + k * std::log(p) + (n - k) * std::log(1.0 - p));
}


DataFrameCpp powerRiskDiffExactcpp(
    const int n,
    const double riskDiffH0,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha,
    const bool calculateAttainedAlpha) {

  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);
  int n1 = static_cast<int>(std::round(n * r));
  int n2 = n - n1;

  const int k = (n1 + 1) * (n2 + 1);

  // Precompute T and flattened mapping (flat_y1, flat_y2)
  std::vector<double> T; T.reserve(k);
  std::vector<int> flat_y1; flat_y1.reserve(k);
  std::vector<int> flat_y2; flat_y2.reserve(k);

  for (int y1 = 0; y1 <= n1; ++y1) {
    for (int y2 = 0; y2 <= n2; ++y2) {
      auto pp = remlRiskDiff(n1, y1, n2, y2, riskDiffH0);
      double p1_hat = static_cast<double>(y1) / static_cast<double>(n1);
      double p2_hat = static_cast<double>(y2) / static_cast<double>(n2);
      double md = p1_hat - p2_hat - riskDiffH0;
      double mv = pp[0] * (1.0 - pp[0]) / static_cast<double>(n1) +
        pp[1] * (1.0 - pp[1]) / static_cast<double>(n2);
      if (mv < 1.0e-8) mv = 1.0e-8;
      T.push_back(md / std::sqrt(mv));
      flat_y1.push_back(y1);
      flat_y2.push_back(y2);
    }
  }

  // Order indices by ascending T
  std::vector<int> order(k);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int i, int j) { return T[i] < T[j]; });

  // Build Tsorted and find unique group boundaries (idx)
  std::vector<double> Tsorted(k);
  for (int i = 0; i < k; ++i) Tsorted[i] = T[order[i]];

  std::vector<int> idx; idx.reserve(k / 4 + 4);
  idx.push_back(0);
  for (int i = 1; i < k; ++i) if (Tsorted[i] != Tsorted[i - 1]) idx.push_back(i);
  int k1 = static_cast<int>(idx.size());
  std::vector<double> Tunique; Tunique.reserve(k1);
  for (int i = 0; i < k1; ++i) Tunique.push_back(Tsorted[idx[i]]);
  idx.push_back(k); // sentinel

  // Determine direction: higher T means better response if pi1 - pi2 > riskDiffH0
  bool directionUpper = (pi1 - pi2) > riskDiffH0;

  // Precompute log-choose tables for n1 and n2 to speed repeated
  // binomial pmf calculations
  std::vector<double> logchoose_n1 = make_log_choose(n1);
  std::vector<double> logchoose_n2 = make_log_choose(n2);

  // Flatten mapping from order positions back to flattened index -> y1,y2
  // but we already have flat_y1/flat_y2; we'll use order[j] to identify which (y1,y2)
  // that corresponds to j-th smallest T.

  // f(p2) computes the critical value (signed)
  auto f = [&](double p2)->double {
    double p1 = p2 + riskDiffH0;

    // compute q1 and q2 arrays
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int y = 0; y <= n1; ++y)
      q1[y] = binomial_pmf_from_logchoose(y, n1, p1, logchoose_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2, logchoose_n2);

    // sum probabilities group-wise until exceed alpha
    if (directionUpper) {
      double s = 0.0;
      int i_group;
      for (i_group = k1 - 1; i_group >= 0; --i_group) {
        // sum all entries in group i_group
        for (int j = idx[i_group]; j < idx[i_group + 1]; ++j) {
          int flatIndex = order[j];
          s += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
        }
        if (s > alpha) break;
      }
      double aval;
      if (i_group == k1 - 1) aval = Tunique[k1 - 1] + 1.0; // impossible to reject
      else aval = Tunique[i_group + 1];
      return -aval; // negate to make it a minimization problem
    } else {
      double s = 0.0;
      int i_group;
      for (i_group = 0; i_group < k1; ++i_group) {
        for (int j = idx[i_group]; j < idx[i_group + 1]; ++j) {
          int flatIndex = order[j];
          s += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
        }
        if (s > alpha) break;
      }
      double aval;
      if (i_group == 0) aval = Tunique[0] - 1.0; // impossible to reject
      else aval = Tunique[i_group - 1];
      return aval;
    }
  };

  // Domain for p2
  double pi2lower = std::max(0.0, -riskDiffH0);
  double pi2upper = std::min(1.0, 1.0 - riskDiffH0);
  if (pi2upper < pi2lower) std::swap(pi2upper, pi2lower);

  // Partition domain into K subintervals, run local minimizer on each,
  // pick global best
  const int K = 100;
  const double delta = (pi2upper - pi2lower) / static_cast<double>(K);
  std::vector<double> a(K), b(K);

  for (int ii = 0; ii < K; ++ii) {
    double lo = pi2lower + ii * delta;
    double hi = lo + delta;
    // minimize f on [lo,hi]
    auto res = mini(f, lo, hi);
    a[ii] = res.first;
    b[ii] = res.second;
  }

  // pick interval with minimal b
  int best_i = 0;
  for (int ii = 1; ii < K; ++ii) if (b[ii] < b[best_i]) best_i = ii;
  double pi2star = a[best_i];
  double t = directionUpper ? -b[best_i] : b[best_i];

  // g(p1, p2) computes the rejection probability (power) for given probabilities
  auto g = [&](double p1_val, double p2_val)->double {
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int y = 0; y <= n1; ++y)
      q1[y] = binomial_pmf_from_logchoose(y, n1, p1_val, logchoose_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2_val, logchoose_n2);

    double preject = 0.0;
    for (int flatIndex = 0; flatIndex < k; ++flatIndex) {
      if ((directionUpper && (T[flatIndex] >= t)) ||
          ((!directionUpper) && (T[flatIndex] <= t))) {
        preject += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
      }
    }
    return preject;
  };

  double power = g(pi1, pi2);

  DataFrameCpp result;

  if (calculateAttainedAlpha) {
    // h(p2) = -preject under H0 -> minimize h to find maximal preject (attained alpha)
    auto h = [&](double p2)->double {
      double p1 = p2 + riskDiffH0;
      return -g(p1, p2);
    };

    // search K intervals again for h minimization
    for (int ii = 0; ii < K; ++ii) {
      double lo = pi2lower + ii * delta;
      double hi = lo + delta;
      auto res = mini(h, lo, hi);
      a[ii] = res.first;
      b[ii] = res.second;
    }
    int best2 = 0;
    for (int ii = 1; ii < K; ++ii) if (b[ii] < b[best2]) best2 = ii;
    double attainedAlpha = -b[best2];

    result.push_back(alpha, "alpha");
    result.push_back(attainedAlpha, "attainedAlpha");
    result.push_back(power, "power");
    result.push_back(n, "n");
    result.push_back(riskDiffH0, "riskDiffH0");
    result.push_back(pi1, "pi1");
    result.push_back(pi2, "pi2");
    result.push_back(allocationRatioPlanned, "allocationRatioPlanned");
    result.push_back(t, "zstatRiskDiffBound");
    result.push_back(pi2star, "pi2star");
  } else {
    result.push_back(alpha, "alpha");
    result.push_back(power, "power");
    result.push_back(n, "n");
    result.push_back(riskDiffH0, "riskDiffH0");
    result.push_back(pi1, "pi1");
    result.push_back(pi2, "pi2");
    result.push_back(allocationRatioPlanned, "allocationRatioPlanned");
    result.push_back(t, "zstatRiskDiffBound");
    result.push_back(pi2star, "pi2star");
  }

  return result;
}


//' @title Power for Exact Unconditional Test of Risk Difference
//' @description Obtains the power given sample size for exact unconditional
//' test of risk difference.
//'
//' @param n The total sample size.
//' @param riskDiffH0 The risk difference under the null hypothesis.
//'   Defaults to 0.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified one-sided significance level.
//'
//' * \code{attainedAlpha}: The attained one-sided significance level if
//'   requested.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskDiffH0}: The risk difference under the null hypothesis.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskDiffBound}: The critical value on the scale of
//'   score test statistic for risk difference.
//'
//' * \code{pi2star}: The response probability in the control group
//'   at which the critical value of the test statistic is attained.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' powerRiskDiffExact(n = 68, pi1 = 0.6, pi2 = 0.25, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame powerRiskDiffExact(
    const int n = NA_INTEGER,
    const double riskDiffH0 = 0,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.025,
    const bool calculateAttainedAlpha = true) {

  auto df = powerRiskDiffExactcpp(
    n, riskDiffH0, pi1, pi2, allocationRatioPlanned, alpha,
    calculateAttainedAlpha);

  return Rcpp::wrap(df);
}


DataFrameCpp samplesizeRiskDiffExactcpp(
    const double beta,
    const double riskDiffH0,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha,
    const bool calculateAttainedAlpha,
    const int max_n_search,
    const int window) {

  // allocation ratio r for group1 (fraction of total sample)
  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Use remlRiskDiff with fractional "n" to obtain p1, p2 for v0
  auto mr = remlRiskDiff(r, r * pi1, 1.0 - r, (1.0 - r) * pi2, riskDiffH0);
  double p1 = mr[0], p2 = mr[1];

  // variance under "H0-ish" and under alternative (normal approx)
  double v0 = p1 * (1.0 - p1) / r + p2 * (1.0 - p2) / (1.0 - r);
  double v1 = pi1 * (1.0 - pi1) / r + pi2 * (1.0 - pi2) / (1.0 - r);

  double z0 = boost_qnorm(1.0 - alpha);
  double z1 = boost_qnorm(1.0 - beta);

  double theta = std::fabs(pi1 - pi2 - riskDiffH0);
  if (theta <= 0.0) throw std::invalid_argument(
    "pi1 - pi2 equals riskDiffH0; sample size is not finite");

  // Normal-approximate continuous sample size
  double n0 = sq(z0 * std::sqrt(v0) + z1 * std::sqrt(v1)) / (theta * theta);

  if (n0 > 0.5 * max_n_search) throw std::invalid_argument(
      std::string("Initial sample size estimate (") +
        std::to_string(static_cast<int>(std::floor(n0))) +
        ") is too large for exact test; consider normal approximation.");

  // Choose starting integer n
  int start_n = static_cast<int>(std::floor(n0));
  if (start_n < 1) start_n = 1;

  const double target = 1.0 - beta;

  // cache power evaluations
  std::unordered_map<int, double> power_cache;
  power_cache.reserve(1024);

  auto power_at = [&](int n)->double {
    if (n <= 0) return 0.0;
    auto it = power_cache.find(n);
    if (it != power_cache.end()) return it->second;
    DataFrameCpp df = powerRiskDiffExactcpp(
      n, riskDiffH0, pi1, pi2, allocationRatioPlanned, alpha, false);
    double p = df.get<double>("power")[0];
    power_cache.emplace(n, p);
    return p;
  };

  // Step-down: if start_n's power >= target, decrease until power < target
  int n = start_n;
  double p_n = power_at(n);
  while (p_n >= target && n > 1) {
    --n;
    p_n = power_at(n);
  }

  // Linear upward scan from candidate to max_n_search, enforcing window criterion
  int found_n = -1;
  for (int nn = n + 1; nn <= max_n_search; ++nn) {
    double p_nn = power_at(nn);
    if (p_nn < target) continue;

    // verify next window-1 sizes are also >= target
    bool ok = true;
    int fail_at = -1;
    for (int j = 1; j < window; ++j) {
      double pj = power_at(nn + j);
      if (pj < target) { ok = false; fail_at = j; break; }
    }
    if (ok) { found_n = nn; break; }
    n += fail_at; // skip ahead by fail_at for next iteration
  }

  if (found_n == -1) {
    throw std::runtime_error("No sample size <= max_n_search satisfies the "
                               "windowed power criterion");
  }

  // Return final DataFrame
  DataFrameCpp final_df = powerRiskDiffExactcpp(
    found_n, riskDiffH0, pi1, pi2, allocationRatioPlanned, alpha,
    calculateAttainedAlpha);
  return final_df;
}

//' @title Sample Size for Exact Unconditional Test of Risk Difference
//' @description Obtains the sample size given power for exact unconditional
//' test of risk difference.
//'
//' @param beta The type II error.
//' @param riskDiffH0 The risk difference under the null hypothesis.
//'   Defaults to 0.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The one-sided significance level.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//' @param max_n_search The maximum sample size to search up to. If no
//'   sample size up to this value satisfies the windowed power criterion,
//'   an error is thrown.
//' @param window The number of consecutive sample sizes that must all
//'   satisfy the power criterion to confirm the found sample size.
//'   This is to mitigate non-monotonicity of power in sample size for
//'   the exact test.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified one-sided significance level.
//'
//' * \code{attainedAlpha}: The attained one-sided significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskDiffH0}: The risk difference under the null hypothesis.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskDiffBound}: The critical value on the scale of
//'   score test statistic for risk difference.
//'
//' * \code{pi2star}: The response probability in the control group
//'   at which the critical value of the test statistic is attained.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' samplesizeRiskDiffExact(beta = 0.2, riskDiffH0 = -0.3,
//'                         pi1 = 0.8, pi2 = 0.8, alpha = 0.025)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame samplesizeRiskDiffExact(
    const double beta = NA_REAL,
    const double riskDiffH0 = 0,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.025,
    const bool calculateAttainedAlpha = true,
    const int max_n_search = 1000,
    const int window = 10) {

  auto df = samplesizeRiskDiffExactcpp(
    beta, riskDiffH0, pi1, pi2, allocationRatioPlanned, alpha,
    calculateAttainedAlpha, max_n_search, window);

  return Rcpp::wrap(df);
}


DataFrameCpp powerRiskRatioExactcpp(
    const int n,
    const double riskRatioH0,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha,
    const bool calculateAttainedAlpha) {

  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);
  int n1 = static_cast<int>(std::round(n * r));
  int n2 = n - n1;

  const int k = (n1 + 1) * (n2 + 1);

  // Precompute T and flattened mapping (flat_y1, flat_y2)
  std::vector<double> T; T.reserve(k);
  std::vector<int> flat_y1; flat_y1.reserve(k);
  std::vector<int> flat_y2; flat_y2.reserve(k);

  for (int y1 = 0; y1 <= n1; ++y1) {
    for (int y2 = 0; y2 <= n2; ++y2) {
      auto pp = remlRiskRatio(n1, y1, n2, y2, riskRatioH0);
      double p1 = pp[0], p2 = pp[1];
      double p1hat = static_cast<double>(y1) / static_cast<double>(n1);
      double p2hat = static_cast<double>(y2) / static_cast<double>(n2);
      double md = p1hat - p2hat * riskRatioH0;
      double mv = p1 * (1.0 - p1) / static_cast<double>(n1) +
        (riskRatioH0 * riskRatioH0) * p2 * (1.0 - p2) / static_cast<double>(n2);
      if (mv < 1e-8) mv = 1e-8;
      T.push_back(md / std::sqrt(mv));
      flat_y1.push_back(y1);
      flat_y2.push_back(y2);
    }
  }

  // order indices by ascending T
  std::vector<int> order(k);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int i, int j) { return T[i] < T[j]; });

  // Build Tsorted and find unique group boundaries (idx)
  std::vector<double> Tsorted(k);
  for (int i = 0; i < k; ++i) Tsorted[i] = T[order[i]];

  std::vector<int> idx; idx.reserve(k / 4 + 4);
  idx.push_back(0);
  for (int i = 1; i < k; ++i) if (Tsorted[i] != Tsorted[i - 1]) idx.push_back(i);
  int k1 = static_cast<int>(idx.size());
  std::vector<double> Tunique; Tunique.reserve(k1);
  for (int i = 0; i < k1; ++i) Tunique.push_back(Tsorted[idx[i]]);
  idx.push_back(k); // sentinel

  // determine direction: higher T means better response if pi1 / pi2 > riskRatioH0
  bool directionUpper = (pi1 / pi2) > riskRatioH0;

  // Precompute log-choose tables for n1 and n2 to speed repeated
  // binomial pmf calculations
  std::vector<double> logchoose_n1 = make_log_choose(n1);
  std::vector<double> logchoose_n2 = make_log_choose(n2);

  // Flatten mapping from order positions back to flattened index -> y1,y2
  // but we already have flat_y1/flat_y2; we'll use order[j] to identify which (y1,y2)
  // that corresponds to j-th smallest T.

  // f(p2) computes the critical value (signed)
  auto f = [&](double p2)->double {
    double p1 = p2 * riskRatioH0;

    // compute q1 and q2 arrays
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int y = 0; y <= n1; ++y)
      q1[y] = binomial_pmf_from_logchoose(y, n1, p1, logchoose_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2, logchoose_n2);

    // sum probabilities group-wise until exceed alpha
    if (directionUpper) {
      double s = 0.0;
      int i_group;
      for (i_group = k1 - 1; i_group >= 0; --i_group) {
        // sum all entries in group i_group
        for (int j = idx[i_group]; j < idx[i_group + 1]; ++j) {
          int flatIndex = order[j];
          s += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
        }
        if (s > alpha) break;
      }
      double aval;
      if (i_group == k1 - 1) aval = Tunique[k1 - 1] + 1.0; // impossible to reject
      else aval = Tunique[i_group + 1];
      return -aval; // negate to make it a minimization problem
    } else {
      double s = 0.0;
      int i_group;
      for (i_group = 0; i_group < k1; ++i_group) {
        for (int j = idx[i_group]; j < idx[i_group + 1]; ++j) {
          int flatIndex = order[j];
          s += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
        }
        if (s > alpha) break;
      }
      double aval;
      if (i_group == 0) aval = Tunique[0] - 1.0; // impossible to reject
      else aval = Tunique[i_group - 1];
      return aval;
    }
  };

  // domain for p2
  double pi2lower = 0.0;
  double pi2upper = std::min(1.0, 1.0 / riskRatioH0);
  if (pi2upper < pi2lower) std::swap(pi2upper, pi2lower);

  // Partition domain into K subintervals, run local minimizer on each,
  // pick global best
  const int K = 100;
  const double delta = (pi2upper - pi2lower) / static_cast<double>(K);
  std::vector<double> a(K), b(K);

  for (int ii = 0; ii < K; ++ii) {
    double lo = pi2lower + ii * delta;
    double hi = lo + delta;
    // minimize f on [lo,hi]
    auto res = mini(f, lo, hi);
    a[ii] = res.first;
    b[ii] = res.second;
  }

  // pick interval with minimal b
  int best_i = 0;
  for (int ii = 1; ii < K; ++ii) if (b[ii] < b[best_i]) best_i = ii;
  double pi2star = a[best_i];
  double t = directionUpper ? -b[best_i] : b[best_i];

  // g(p1,p2) computes the rejection probability (power) for given probabilities
  auto g = [&](double p1_val, double p2_val)->double {
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int y = 0; y <= n1; ++y)
      q1[y] = binomial_pmf_from_logchoose(y, n1, p1_val, logchoose_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2_val, logchoose_n2);

    double preject = 0.0;
    for (int flatIndex = 0; flatIndex < k; ++flatIndex) {
      if ((directionUpper && (T[flatIndex] >= t)) ||
          ((!directionUpper) && (T[flatIndex] <= t))) {
        preject += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
      }
    }
    return preject;
  };

  double power = g(pi1, pi2);

  DataFrameCpp result;

  if (calculateAttainedAlpha) {
    // h(p2) = -preject under H0 -> minimize h to find maximal preject (attained alpha)
    auto h = [&](double p2)->double {
      double p1 = p2 * riskRatioH0;
      return -g(p1, p2);
    };

    for (int ii = 0; ii < K; ++ii) {
      double lo = pi2lower + ii * delta;
      double hi = lo + delta;
      auto res = mini(h, lo, hi);
      a[ii] = res.first;
      b[ii] = res.second;
    }
    int best2 = 0;
    for (int ii = 1; ii < K; ++ii) if (b[ii] < b[best2]) best2 = ii;
    double attainedAlpha = -b[best2];

    result.push_back(alpha, "alpha");
    result.push_back(attainedAlpha, "attainedAlpha");
    result.push_back(power, "power");
    result.push_back(n, "n");
    result.push_back(riskRatioH0, "riskRatioH0");
    result.push_back(pi1, "pi1");
    result.push_back(pi2, "pi2");
    result.push_back(allocationRatioPlanned, "allocationRatioPlanned");
    result.push_back(t, "zstatRiskRatioBound");
    result.push_back(pi2star, "pi2star");
  } else {
    result.push_back(alpha, "alpha");
    result.push_back(power, "power");
    result.push_back(n, "n");
    result.push_back(riskRatioH0, "riskRatioH0");
    result.push_back(pi1, "pi1");
    result.push_back(pi2, "pi2");
    result.push_back(allocationRatioPlanned, "allocationRatioPlanned");
    result.push_back(t, "zstatRiskRatioBound");
    result.push_back(pi2star, "pi2star");
  }

  return result;
}


//' @title Power for Exact Unconditional Test of Risk Ratio
//' @description Obtains the power given sample size for exact unconditional
//' test of risk ratio.
//'
//' @param n The total sample size.
//' @param riskRatioH0 The risk ratio under the null hypothesis.
//'   Defaults to 1.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified one-sided significance level.
//'
//' * \code{attainedAlpha}: The attained one-sided significance level if
//'   requested.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskRatioH0}: The risk ratio under the null hypothesis.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskRatioBound}: The critical value on the scale of
//'   score test statistic for risk ratio.
//'
//' * \code{pi2star}: The response probability in the control group
//'   at which the critical value of the test statistic is attained.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' powerRiskRatioExact(n = 68, pi1 = 0.6, pi2 = 0.25, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame powerRiskRatioExact(
    const int n = NA_INTEGER,
    const double riskRatioH0 = 1,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.025,
    const bool calculateAttainedAlpha = true) {

  auto df = powerRiskRatioExactcpp(
    n, riskRatioH0, pi1, pi2, allocationRatioPlanned, alpha,
    calculateAttainedAlpha);

  return Rcpp::wrap(df);
}


DataFrameCpp samplesizeRiskRatioExactcpp(
    const double beta,
    const double riskRatioH0,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha,
    const bool calculateAttainedAlpha,
    const int max_n_search,
    const int window) {

  // allocation ratio r for group1 (fraction of total sample)
  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Use remlRiskRatio with fractional "n" to obtain p1, p2 for v0
  auto mr = remlRiskRatio(r, r * pi1, 1.0 - r, (1.0 - r) * pi2, riskRatioH0);
  double p1 = mr[0], p2 = mr[1];

  // variance under "H0-ish" and under alternative (normal approx)
  double v0 = (1.0 - p1) / (r * p1) + (1.0 - p2) / ((1.0 - r) * p2);
  double v1 = (1.0 - pi1) / (r * pi1) + (1.0 - pi2) / ((1.0 - r) * pi2);

  double z0 = boost_qnorm(1.0 - alpha);
  double z1 = boost_qnorm(1.0 - beta);

  double theta = std::fabs(std::log(pi1 / pi2) - std::log(riskRatioH0));
  if (theta <= 0.0) throw std::invalid_argument(
    "pi1 / pi2 equals riskRatioH0; sample size is not finite");

  // Normal-approximate continuous sample size
  double n0 = sq(z0 * std::sqrt(v0) + z1 * std::sqrt(v1)) / (theta * theta);

  if (n0 > 0.5 * max_n_search) throw std::invalid_argument(
      std::string("Initial sample size estimate (") +
        std::to_string(static_cast<int>(std::floor(n0))) +
        ") is too large for exact test; consider normal approximation.");

  // Choose starting integer n
  int start_n = static_cast<int>(std::floor(n0));
  if (start_n < 1) start_n = 1;

  const double target = 1.0 - beta;

  // cache power evaluations
  std::unordered_map<int, double> power_cache;
  power_cache.reserve(1024);

  auto power_at = [&](int n)->double {
    if (n <= 0) return 0.0;
    auto it = power_cache.find(n);
    if (it != power_cache.end()) return it->second;
    DataFrameCpp df = powerRiskRatioExactcpp(
      n, riskRatioH0, pi1, pi2, allocationRatioPlanned, alpha, false);
    double p = df.get<double>("power")[0];
    power_cache.emplace(n, p);
    return p;
  };

  // Step-down: if start_n's power >= target, decrease until power < target
  int n = start_n;
  double p_n = power_at(n);
  while (p_n >= target && n > 1) {
    --n;
    p_n = power_at(n);
  }

  // Linear upward scan from candidate to max_n_search, enforcing window criterion
  int found_n = -1;
  for (int nn = n + 1; nn <= max_n_search; ++nn) {
    double p_nn = power_at(nn);
    if (p_nn < target) continue;

    // verify next window-1 sizes are also >= target
    bool ok = true;
    int fail_at = -1;
    for (int j = 1; j < window; ++j) {
      double pj = power_at(nn + j);
      if (pj < target) { ok = false; fail_at = j; break; }
    }
    if (ok) { found_n = nn; break; }
    n += fail_at; // skip ahead by fail_at for next iteration
  }

  if (found_n == -1) {
    throw std::runtime_error("No sample size <= max_n_search satisfies the "
                               "windowed power criterion");
  }

  // Return final DataFrame
  DataFrameCpp final_df = powerRiskRatioExactcpp(
    found_n, riskRatioH0, pi1, pi2, allocationRatioPlanned, alpha,
    calculateAttainedAlpha);
  return final_df;
}


//' @title Sample Size for Exact Unconditional Test of Risk Ratio
//' @description Obtains the sample size given power for exact unconditional
//' test of risk ratio.
//'
//' @param beta The type II error.
//' @param riskRatioH0 The risk ratio under the null hypothesis.
//'   Defaults to 1.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The one-sided significance level.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//' @param max_n_search The maximum sample size to search up to. If no
//'   sample size up to this value satisfies the windowed power criterion,
//'   an error is thrown.
//' @param window The number of consecutive sample sizes that must all
//'   satisfy the power criterion to confirm the found sample size.
//'   This is to mitigate non-monotonicity of power in sample size for
//'   the exact test.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified one-sided significance level.
//'
//' * \code{attainedAlpha}: The attained one-sided significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskRatioH0}: The risk ratio under the null hypothesis.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskRatioBound}: The critical value on the scale of
//'   score test statistic for risk ratio.
//'
//' * \code{pi2star}: The response probability in the control group
//'   at which the critical value of the test statistic is attained.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' samplesizeRiskRatioExact(beta = 0.2, riskRatioH0 = 0.8,
//'                          pi1 = 0.95, pi2 = 0.95, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame samplesizeRiskRatioExact(
    const double beta = NA_REAL,
    const double riskRatioH0 = 1,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.025,
    const bool calculateAttainedAlpha = true,
    const int max_n_search = 1000,
    const int window = 10) {

  auto df = samplesizeRiskRatioExactcpp(
    beta, riskRatioH0, pi1, pi2, allocationRatioPlanned, alpha,
    calculateAttainedAlpha, max_n_search, window);

  return Rcpp::wrap(df);
}




DataFrameCpp powerRiskDiffExactEquivcpp(
    const int n,
    const double riskDiffLower,
    const double riskDiffUpper,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha,
    const bool calculateAttainedAlpha) {

  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);
  int n1 = static_cast<int>(std::round(n * r));
  int n2 = n - n1;

  const int k = (n1 + 1) * (n2 + 1);

  // Precompute log-choose tables for binomial pmfs
  std::vector<double> logc_n1 = make_log_choose(n1);
  std::vector<double> logc_n2 = make_log_choose(n2);

  // Build T1 (for riskDiffLower) and flat y maps
  std::vector<double> T1; T1.reserve(k);
  std::vector<int> flat_y1; flat_y1.reserve(k);
  std::vector<int> flat_y2; flat_y2.reserve(k);

  for (int y1 = 0; y1 <= n1; ++y1) {
    for (int y2 = 0; y2 <= n2; ++y2) {
      auto pp = remlRiskDiff(n1, y1, n2, y2, riskDiffLower);
      double p1hat = static_cast<double>(y1) / static_cast<double>(n1);
      double p2hat = static_cast<double>(y2) / static_cast<double>(n2);
      double md = p1hat - p2hat - riskDiffLower;
      double mv = pp[0] * (1.0 - pp[0]) / static_cast<double>(n1) +
        pp[1] * (1.0 - pp[1]) / static_cast<double>(n2);
      if (mv < 1e-8) mv = 1e-8;
      T1.push_back(md / std::sqrt(mv));
      flat_y1.push_back(y1);
      flat_y2.push_back(y2);
    }
  }

  // order1 and unique breakpoints for T1
  std::vector<int> order1(k);
  std::iota(order1.begin(), order1.end(), 0);
  std::sort(order1.begin(), order1.end(), [&](int a, int b){ return T1[a] < T1[b]; });

  std::vector<double> T1sorted(k);
  for (int i = 0; i < k; ++i) T1sorted[i] = T1[order1[i]];

  std::vector<int> idx1; idx1.reserve(k/4 + 4);
  idx1.push_back(0);
  for (int i = 1; i < k; ++i) if (T1sorted[i] != T1sorted[i-1]) idx1.push_back(i);
  int k1 = static_cast<int>(idx1.size());
  std::vector<double> T1unique; T1unique.reserve(k1);
  for (int i = 0; i < k1; ++i) T1unique.push_back(T1sorted[idx1[i]]);
  idx1.push_back(k); // sentinel

  // f1: given p2 under H10, compute signed critical value (-aval)
  auto f1 = [&](double p2)->double {
    double p1 = p2 + riskDiffLower;

    // compute binomial pmfs for rows and cols
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int x = 0; x <= n1; ++x)
      q1[x] = binomial_pmf_from_logchoose(x, n1, p1, logc_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2, logc_n2);

    // iterate groups from largest T downwards, accumulating group probabilities
    double s = 0.0;
    int i_group;
    for (i_group = k1 - 1; i_group >= 0; --i_group) {
      for (int j = idx1[i_group]; j < idx1[i_group + 1]; ++j) {
        int flatIndex = order1[j];
        s += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
      }
      if (s > alpha) break;
    }
    double aval;
    if (i_group == k1 - 1) aval = T1unique[k1 - 1] + 1.0; // cannot reject
    else aval = T1unique[i_group + 1];
    return -aval;
  };

  // find critical value for H10 by partitioned minimization over p2 in valid domain
  double pi2lower1 = std::max(0.0, -riskDiffLower);
  double pi2upper1 = std::min(1.0, 1.0 - riskDiffLower);
  const int K = 100;
  double delta1 = (pi2upper1 - pi2lower1) / static_cast<double>(K);
  std::vector<double> a1(K), b1(K);

  for (int ii = 0; ii < K; ++ii) {
    double lo = pi2lower1 + ii * delta1;
    double hi = lo + delta1;
    auto res = mini(f1, lo, hi);
    a1[ii] = res.first;
    b1[ii] = res.second;
  }
  double t1 = - *std::min_element(b1.begin(), b1.end());

  // -------------------- T2 (for riskDiffUpper) --------------------
  std::vector<double> T2; T2.reserve(k);
  // Reuse flat_y1, flat_y2 ordering (same layout)
  for (int y1 = 0; y1 <= n1; ++y1) {
    for (int y2 = 0; y2 <= n2; ++y2) {
      auto pp = remlRiskDiff(n1, y1, n2, y2, riskDiffUpper);
      double p1hat = static_cast<double>(y1) / static_cast<double>(n1);
      double p2hat = static_cast<double>(y2) / static_cast<double>(n2);
      double md = p1hat - p2hat - riskDiffUpper;
      double mv = pp[0] * (1.0 - pp[0]) / static_cast<double>(n1) +
        pp[1] * (1.0 - pp[1]) / static_cast<double>(n2);
      if (mv < 1e-8) mv = 1e-8;
      T2.push_back(md / std::sqrt(mv));
    }
  }

  double t2;
  double pi2lower2 = std::max(0.0, -riskDiffUpper);
  double pi2upper2 = std::min(1.0, 1.0 - riskDiffUpper);
  double delta2 = (pi2upper2 - pi2lower2) / static_cast<double>(K);

  if (riskDiffLower == -riskDiffUpper) {
    t2 = -t1;
  } else {
    // build order2, unique breakpoints
    std::vector<int> order2(k);
    std::iota(order2.begin(), order2.end(), 0);
    std::sort(order2.begin(), order2.end(),
              [&](int a, int b){ return T2[a] < T2[b]; });

    std::vector<double> T2sorted(k);
    for (int i = 0; i < k; ++i) T2sorted[i] = T2[order2[i]];

    std::vector<int> idx2; idx2.reserve(k/4 + 4);
    idx2.push_back(0);
    for (int i = 1; i < k; ++i) if (T2sorted[i] != T2sorted[i-1]) idx2.push_back(i);
    int k2 = static_cast<int>(idx2.size());
    std::vector<double> T2unique; T2unique.reserve(k2);
    for (int i = 0; i < k2; ++i) T2unique.push_back(T2sorted[idx2[i]]);
    idx2.push_back(k);

    auto f2 = [&](double p2)->double {
      double p1 = p2 + riskDiffUpper;

      std::vector<double> q1(n1 + 1), q2(n2 + 1);
      for (int x = 0; x <= n1; ++x)
        q1[x] = binomial_pmf_from_logchoose(x, n1, p1, logc_n1);
      for (int y = 0; y <= n2; ++y)
        q2[y] = binomial_pmf_from_logchoose(y, n2, p2, logc_n2);

      double s = 0.0;
      int i_group;
      for (i_group = 0; i_group < k2; ++i_group) {
        for (int j = idx2[i_group]; j < idx2[i_group + 1]; ++j) {
          int flatIndex = order2[j];
          s += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
        }
        if (s > alpha) break;
      }

      double aval;
      if (i_group == 0) aval = T2unique[0] - 1.0;
      else aval = T2unique[i_group - 1];
      return aval;
    };

    // partition minimization
    std::vector<double> a2(K), b2(K);
    for (int ii = 0; ii < K; ++ii) {
      double lo = pi2lower2 + ii * delta2;
      double hi = lo + delta2;
      auto res = mini(f2, lo, hi);
      a2[ii] = res.first;
      b2[ii] = res.second;
    }
    t2 = *std::min_element(b2.begin(), b2.end());
  }

  // -------------------- compute power under (pi1,pi2) --------------------
  // g(p1, p2) computes the rejection probability (power) for given probabilities
  auto g = [&](double p1_val, double p2_val)->double {
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int x = 0; x <= n1; ++x)
      q1[x] = binomial_pmf_from_logchoose(x, n1, p1_val, logc_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2_val, logc_n2);

    double preject = 0.0;
    for (int flatIndex = 0; flatIndex < k; ++flatIndex) {
      if ((T1[flatIndex] >= t1) && (T2[flatIndex] <= t2)) {
        preject += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
      }
    }
    return preject;
  };

  double power = g(pi1, pi2);


  DataFrameCpp out;

  if (calculateAttainedAlpha) {
    // h10 function: for a given p2, compute -preject under H10 (for use with mini)
    auto h10 = [&](double p2)->double {
      double p1 = p2 + riskDiffLower;
      return -g(p1, p2);
    };

    // partition and minimize h10 over pi2lower1..pi2upper1
    for (int ii = 0; ii < K; ++ii) {
      double lo = pi2lower1 + ii * delta1;
      double hi = lo + delta1;
      auto res = mini(h10, lo, hi);
      b1[ii] = res.second;
    }
    double attainedAlphaH10 = - *std::min_element(b1.begin(), b1.end());

    double attainedAlphaH20;
    if (riskDiffLower == -riskDiffUpper) {
      attainedAlphaH20 = attainedAlphaH10;
    } else {
      // h20 function: for a given p2, compute -preject under H20 (for use with mini)
      auto h20 = [&](double p2)->double {
        double p1 = p2 + riskDiffUpper;
        return -g(p1, p2);
      };

      // partition and minimize h20 over across pi2lower2..pi2upper2
      std::vector<double> b2(K);
      for (int ii = 0; ii < K; ++ii) {
        double lo = pi2lower2 + ii * delta2;
        double hi = lo + delta2;
        auto res = mini(h20, lo, hi);
        b2[ii] = res.second;
      }

      attainedAlphaH20 = - *std::min_element(b2.begin(), b2.end());
    }

    out.push_back(alpha, "alpha");
    out.push_back(attainedAlphaH10, "attainedAlphaH10");
    out.push_back(attainedAlphaH20, "attainedAlphaH20");
    out.push_back(power, "power");
    out.push_back(n, "n");
    out.push_back(riskDiffLower, "riskDiffLower");
    out.push_back(riskDiffUpper, "riskDiffUpper");
    out.push_back(pi1, "pi1");
    out.push_back(pi2, "pi2");
    out.push_back(pi1 - pi2, "riskDiff");
    out.push_back(allocationRatioPlanned, "allocationRatioPlanned");
    out.push_back(t1, "zstatRiskDiffLower");
    out.push_back(t2, "zstatRiskDiffUpper");
  } else {
    out.push_back(alpha, "alpha");
    out.push_back(power, "power");
    out.push_back(n, "n");
    out.push_back(riskDiffLower, "riskDiffLower");
    out.push_back(riskDiffUpper, "riskDiffUpper");
    out.push_back(pi1, "pi1");
    out.push_back(pi2, "pi2");
    out.push_back(pi1 - pi2, "riskDiff");
    out.push_back(allocationRatioPlanned, "allocationRatioPlanned");
    out.push_back(t1, "zstatRiskDiffLower");
    out.push_back(t2, "zstatRiskDiffUpper");
  }

  return out;
}


//' @title Power for Exact Unconditional Test of Equivalence in Risk
//' Difference
//' @description Obtains the power given sample size for exact unconditional
//' test of equivalence in risk difference.
//'
//' @param n The total sample size.
//' @param riskDiffLower The lower equivalence limit of risk difference.
//' @param riskDiffUpper The upper equivalence limit of risk difference.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified significance level for each of the two
//'   one-sided tests.
//'
//' * \code{attainedAlphaH10}: The attained significance level under H10
//'   if requested.
//'
//' * \code{attainedAlphaH20}: The attained significance level under H20
//'   if requested.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskDiffLower}: The lower equivalence limit of risk difference.
//'
//' * \code{riskDiffUpper}: The upper equivalence limit of risk difference.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{riskDiff}: The risk difference.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskDiffLower}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   lower equivalence limit.
//'
//' * \code{zstatRiskDiffUpper}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   upper equivalence limit.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' powerRiskDiffExactEquiv(
//'   n = 200, riskDiffLower = -0.2, riskDiffUpper = 0.2,
//'   pi1 = 0.775, pi2 = 0.775, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame powerRiskDiffExactEquiv(
    const int n = NA_INTEGER,
    const double riskDiffLower = NA_REAL,
    const double riskDiffUpper = NA_REAL,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05,
    const bool calculateAttainedAlpha = true) {

  auto df = powerRiskDiffExactEquivcpp(
    n, riskDiffLower, riskDiffUpper, pi1, pi2, allocationRatioPlanned,
    alpha, calculateAttainedAlpha);

  return Rcpp::wrap(df);
}


DataFrameCpp samplesizeRiskDiffExactEquivcpp(
    const double beta,
    const double riskDiffLower,
    const double riskDiffUpper,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha,
    const bool calculateAttainedAlpha,
    const int max_n_search,
    const int window) {

  // allocation ratio r for group1 (fraction of total sample)
  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  double diff = pi1 - pi2;
  double theta = std::min(diff - riskDiffLower, riskDiffUpper - diff);
  if (theta <= 0.0) throw std::invalid_argument(
    "pi1 - pi2 is outside or on the equivalence margins; sample size is not finite");

  double v1 = pi1 * (1.0 - pi1) / r + pi2 * (1.0 - pi2) / (1.0 - r);
  double z0 = boost_qnorm(1.0 - alpha);
  double z1 = boost_qnorm(1.0 - beta);

  // Normal-approximate continuous sample size
  double n0 = sq(z0 + z1) * v1 / (theta * theta);

  if (n0 > 0.5 * max_n_search) throw std::invalid_argument(
      std::string("Initial sample size estimate (") +
        std::to_string(static_cast<int>(std::floor(n0))) +
        ") is too large for exact test; consider normal approximation.");

  // Choose starting integer n
  int start_n = static_cast<int>(std::floor(n0));
  if (start_n < 1) start_n = 1;

  const double target = 1.0 - beta;

  // cache power evaluations
  std::unordered_map<int, double> power_cache;
  power_cache.reserve(1024);

  auto power_at = [&](int n)->double {
    if (n <= 0) return 0.0;
    auto it = power_cache.find(n);
    if (it != power_cache.end()) return it->second;
    DataFrameCpp df = powerRiskDiffExactEquivcpp(
      n, riskDiffLower, riskDiffUpper, pi1, pi2, allocationRatioPlanned,
      alpha, false);
    double p = df.get<double>("power")[0];
    power_cache.emplace(n, p);
    return p;
  };

  // Step-down: if start_n's power >= target, decrease until power < target
  int n = start_n;
  double p_n = power_at(n);
  while (p_n >= target && n > 1) {
    --n;
    p_n = power_at(n);
  }

  // Linear upward scan from candidate to max_n_search, enforcing window criterion
  int found_n = -1;
  for (int nn = n + 1; nn <= max_n_search; ++nn) {
    double p_nn = power_at(nn);
    if (p_nn < target) continue;

    // verify next window-1 sizes are also >= target
    bool ok = true;
    int fail_at = -1;
    for (int j = 1; j < window; ++j) {
      double pj = power_at(nn + j);
      if (pj < target) { ok = false; fail_at = j; break; }
    }
    if (ok) { found_n = nn; break; }
    n += fail_at; // skip ahead by fail_at for next iteration
  }

  if (found_n == -1) {
    throw std::runtime_error("No sample size <= max_n_search satisfies the "
                               "windowed power criterion");
  }

  // Return final DataFrame
  DataFrameCpp final_df = powerRiskDiffExactEquivcpp(
    found_n, riskDiffLower, riskDiffUpper, pi1, pi2, allocationRatioPlanned,
    alpha, calculateAttainedAlpha);
  return final_df;
}


//' @title Sample Size for Exact Unconditional Test of Equivalence in Risk
//' Difference
//' @description Obtains the sample size given power for exact unconditional
//' test of equivalence in risk difference.
//'
//' @param beta The type II error.
//' @param riskDiffLower The lower equivalence limit of risk difference.
//' @param riskDiffUpper The upper equivalence limit of risk difference.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//' @param max_n_search The maximum sample size to search up to. If no
//'   sample size up to this value satisfies the windowed power criterion,
//'   an error is thrown.
//' @param window The number of consecutive sample sizes that must all
//'   satisfy the power criterion to confirm the found sample size.
//'   This is to mitigate non-monotonicity of power in sample size for
//'   the exact test.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified significance level for each of the two
//'   one-sided tests.
//'
//' * \code{attainedAlpha}: The attained significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskDiffLower}: The lower equivalence limit of risk difference.
//'
//' * \code{riskDiffUpper}: The upper equivalence limit of risk difference.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{riskDiff}: The risk difference.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskDiffLower}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   lower equivalence limit.
//'
//' * \code{zstatRiskDiffUpper}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   upper equivalence limit.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' samplesizeRiskDiffExactEquiv(
//'   beta = 0.2, riskDiffLower = -0.3, riskDiffUpper = 0.3,
//'   pi1 = 0.9, pi2 = 0.9, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame samplesizeRiskDiffExactEquiv(
    const double beta = NA_REAL,
    const double riskDiffLower = NA_REAL,
    const double riskDiffUpper = NA_REAL,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05,
    const bool calculateAttainedAlpha = true,
    const int max_n_search = 1000,
    const int window = 10) {

  auto df = samplesizeRiskDiffExactEquivcpp(
    beta, riskDiffLower, riskDiffUpper, pi1, pi2, allocationRatioPlanned,
    alpha, calculateAttainedAlpha, max_n_search, window);

  return Rcpp::wrap(df);
}



DataFrameCpp powerRiskRatioExactEquivcpp(
    const int n,
    const double riskRatioLower,
    const double riskRatioUpper,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha,
    const bool calculateAttainedAlpha) {

  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);
  int n1 = static_cast<int>(std::round(n * r));
  int n2 = n - n1;

  const int k = (n1 + 1) * (n2 + 1);

  // Precompute log-choose tables for binomial pmfs
  std::vector<double> logc_n1 = make_log_choose(n1);
  std::vector<double> logc_n2 = make_log_choose(n2);

  // Build T1 (for riskRatioLower) and flat y maps
  std::vector<double> T1; T1.reserve(k);
  std::vector<int> flat_y1; flat_y1.reserve(k);
  std::vector<int> flat_y2; flat_y2.reserve(k);

  for (int y1 = 0; y1 <= n1; ++y1) {
    for (int y2 = 0; y2 <= n2; ++y2) {
      auto pp = remlRiskRatio(n1, y1, n2, y2, riskRatioLower);
      double p1hat = static_cast<double>(y1) / static_cast<double>(n1);
      double p2hat = static_cast<double>(y2) / static_cast<double>(n2);
      double md = p1hat - p2hat * riskRatioLower;
      double mv = pp[0] * (1.0 - pp[0]) / static_cast<double>(n1) +
        sq(riskRatioLower) * pp[1] * (1.0 - pp[1]) / static_cast<double>(n2);
      if (mv < 1e-8) mv = 1e-8;
      T1.push_back(md / std::sqrt(mv));
      flat_y1.push_back(y1);
      flat_y2.push_back(y2);
    }
  }

  // order1 and unique breakpoints for T1
  std::vector<int> order1(k);
  std::iota(order1.begin(), order1.end(), 0);
  std::sort(order1.begin(), order1.end(), [&](int a, int b){ return T1[a] < T1[b]; });

  std::vector<double> T1sorted(k);
  for (int i = 0; i < k; ++i) T1sorted[i] = T1[order1[i]];

  std::vector<int> idx1; idx1.reserve(k/4 + 4);
  idx1.push_back(0);
  for (int i = 1; i < k; ++i) if (T1sorted[i] != T1sorted[i-1]) idx1.push_back(i);
  int k1 = static_cast<int>(idx1.size());
  std::vector<double> T1unique; T1unique.reserve(k1);
  for (int i = 0; i < k1; ++i) T1unique.push_back(T1sorted[idx1[i]]);
  idx1.push_back(k); // sentinel

  // f1: given p2 under H10, compute signed critical value (-aval)
  auto f1 = [&](double p2)->double {
    double p1 = p2 * riskRatioLower;

    // compute binomial pmfs for rows and cols
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int x = 0; x <= n1; ++x)
      q1[x] = binomial_pmf_from_logchoose(x, n1, p1, logc_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2, logc_n2);

    // iterate groups from largest T downwards, accumulating group probabilities
    double s = 0.0;
    int i_group;
    for (i_group = k1 - 1; i_group >= 0; --i_group) {
      for (int j = idx1[i_group]; j < idx1[i_group + 1]; ++j) {
        int flatIndex = order1[j];
        s += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
      }
      if (s > alpha) break;
    }
    double aval;
    if (i_group == k1 - 1) aval = T1unique[k1 - 1] + 1.0; // cannot reject
    else aval = T1unique[i_group + 1];
    return -aval;
  };

  // find critical value for H10 by partitioned minimization over p2 in valid domain
  double pi2lower1 = 0;
  double pi2upper1 = std::min(1.0, 1.0 / riskRatioLower);
  const int K = 100;
  double delta1 = (pi2upper1 - pi2lower1) / static_cast<double>(K);
  std::vector<double> a1(K), b1(K);

  for (int ii = 0; ii < K; ++ii) {
    double lo = pi2lower1 + ii * delta1;
    double hi = lo + delta1;
    auto res = mini(f1, lo, hi);
    a1[ii] = res.first;
    b1[ii] = res.second;
  }
  double t1 = - *std::min_element(b1.begin(), b1.end());

  // -------------------- T2 (for riskRatioUpper) --------------------
  std::vector<double> T2; T2.reserve(k);
  // Reuse flat_y1, flat_y2 ordering (same layout)
  for (int y1 = 0; y1 <= n1; ++y1) {
    for (int y2 = 0; y2 <= n2; ++y2) {
      auto pp = remlRiskRatio(n1, y1, n2, y2, riskRatioUpper);
      double p1hat = static_cast<double>(y1) / static_cast<double>(n1);
      double p2hat = static_cast<double>(y2) / static_cast<double>(n2);
      double md = p1hat - p2hat * riskRatioUpper;
      double mv = pp[0] * (1.0 - pp[0]) / static_cast<double>(n1) +
        sq(riskRatioUpper) * pp[1] * (1.0 - pp[1]) / static_cast<double>(n2);
      if (mv < 1e-8) mv = 1e-8;
      T2.push_back(md / std::sqrt(mv));
    }
  }

  double t2;
  double pi2lower2 = 0.0;
  double pi2upper2 = std::min(1.0, 1.0 / riskRatioUpper);
  double delta2 = (pi2upper2 - pi2lower2) / static_cast<double>(K);

  if (std::fabs(riskRatioLower * riskRatioUpper - 1.0) < 1e-8) {
    t2 = -t1;
  } else {
    // build order2, unique breakpoints
    std::vector<int> order2(k);
    std::iota(order2.begin(), order2.end(), 0);
    std::sort(order2.begin(), order2.end(),
              [&](int a, int b){ return T2[a] < T2[b]; });

    std::vector<double> T2sorted(k);
    for (int i = 0; i < k; ++i) T2sorted[i] = T2[order2[i]];

    std::vector<int> idx2; idx2.reserve(k/4 + 4);
    idx2.push_back(0);
    for (int i = 1; i < k; ++i) if (T2sorted[i] != T2sorted[i-1]) idx2.push_back(i);
    int k2 = static_cast<int>(idx2.size());
    std::vector<double> T2unique; T2unique.reserve(k2);
    for (int i = 0; i < k2; ++i) T2unique.push_back(T2sorted[idx2[i]]);
    idx2.push_back(k);

    auto f2 = [&](double p2)->double {
      double p1 = p2 * riskRatioUpper;

      std::vector<double> q1(n1 + 1), q2(n2 + 1);
      for (int x = 0; x <= n1; ++x)
        q1[x] = binomial_pmf_from_logchoose(x, n1, p1, logc_n1);
      for (int y = 0; y <= n2; ++y)
        q2[y] = binomial_pmf_from_logchoose(y, n2, p2, logc_n2);

      double s = 0.0;
      int i_group;
      for (i_group = 0; i_group < k2; ++i_group) {
        for (int j = idx2[i_group]; j < idx2[i_group + 1]; ++j) {
          int flatIndex = order2[j];
          s += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
        }
        if (s > alpha) break;
      }

      double aval;
      if (i_group == 0) aval = T2unique[0] - 1.0;
      else aval = T2unique[i_group - 1];
      return aval;
    };

    // partition minimization
    std::vector<double> a2(K), b2(K);
    for (int ii = 0; ii < K; ++ii) {
      double lo = pi2lower2 + ii * delta2;
      double hi = lo + delta2;
      auto res = mini(f2, lo, hi);
      a2[ii] = res.first;
      b2[ii] = res.second;
    }
    t2 = *std::min_element(b2.begin(), b2.end());
  }

  // -------------------- compute power under (pi1,pi2) --------------------
  // g(p1, p2) computes the rejection probability (power) for given probabilities
  auto g = [&](double p1_val, double p2_val)->double {
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int x = 0; x <= n1; ++x)
      q1[x] = binomial_pmf_from_logchoose(x, n1, p1_val, logc_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2_val, logc_n2);

    double preject = 0.0;
    for (int flatIndex = 0; flatIndex < k; ++flatIndex) {
      if ((T1[flatIndex] >= t1) && (T2[flatIndex] <= t2)) {
        preject += q1[flat_y1[flatIndex]] * q2[flat_y2[flatIndex]];
      }
    }
    return preject;
  };

  double power = g(pi1, pi2);


  DataFrameCpp out;

  if (calculateAttainedAlpha) {
    // h10 function: for a given p2, compute -preject under H10 (for use with mini)
    auto h10 = [&](double p2)->double {
      double p1 = p2 * riskRatioLower;
      return -g(p1, p2);
    };

    // partition and minimize h10 over pi2lower1..pi2upper1
    for (int ii = 0; ii < K; ++ii) {
      double lo = pi2lower1 + ii * delta1;
      double hi = lo + delta1;
      auto res = mini(h10, lo, hi);
      b1[ii] = res.second;
    }
    double attainedAlphaH10 = - *std::min_element(b1.begin(), b1.end());

    double attainedAlphaH20;
    if (std::fabs(riskRatioLower * riskRatioUpper - 1.0) < 1e-8) {
      attainedAlphaH20 = attainedAlphaH10;
    } else {
      // h20 function: for a given p2, compute -preject under H20 (for use with mini)
      auto h20 = [&](double p2)->double {
        double p1 = p2 * riskRatioUpper;
        return -g(p1, p2);
      };

      // partition and minimize h20 over across pi2lower2..pi2upper2
      std::vector<double> b2(K);
      for (int ii = 0; ii < K; ++ii) {
        double lo = pi2lower2 + ii * delta2;
        double hi = lo + delta2;
        auto res = mini(h20, lo, hi);
        b2[ii] = res.second;
      }

      attainedAlphaH20 = - *std::min_element(b2.begin(), b2.end());
    }

    out.push_back(alpha, "alpha");
    out.push_back(attainedAlphaH10, "attainedAlphaH10");
    out.push_back(attainedAlphaH20, "attainedAlphaH20");
    out.push_back(power, "power");
    out.push_back(n, "n");
    out.push_back(riskRatioLower, "riskRatioLower");
    out.push_back(riskRatioUpper, "riskRatioUpper");
    out.push_back(pi1, "pi1");
    out.push_back(pi2, "pi2");
    out.push_back(pi1 / pi2, "riskRatio");
    out.push_back(allocationRatioPlanned, "allocationRatioPlanned");
    out.push_back(t1, "zstatRiskRatioLower");
    out.push_back(t2, "zstatRiskRatioUpper");
  } else {
    out.push_back(alpha, "alpha");
    out.push_back(power, "power");
    out.push_back(n, "n");
    out.push_back(riskRatioLower, "riskRatioLower");
    out.push_back(riskRatioUpper, "riskRatioUpper");
    out.push_back(pi1, "pi1");
    out.push_back(pi2, "pi2");
    out.push_back(pi1 / pi2, "riskRatio");
    out.push_back(allocationRatioPlanned, "allocationRatioPlanned");
    out.push_back(t1, "zstatRiskRatioLower");
    out.push_back(t2, "zstatRiskRatioUpper");
  }

  return out;
}


//' @title Power for Exact Unconditional Test of Equivalence in Risk
//' Ratio
//' @description Obtains the power given sample size for exact unconditional
//' test of equivalence in risk ratio.
//'
//' @param n The total sample size.
//' @param riskRatioLower The lower equivalence limit of risk ratio.
//' @param riskRatioUpper The upper equivalence limit of risk ratio.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified significance level for each of the two
//'   one-sided tests.
//'
//' * \code{attainedAlphaH10}: The attained significance level under H10
//'   if requested.
//'
//' * \code{attainedAlphaH20}: The attained significance level under H20
//'   if requested.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskRatioLower}: The lower equivalence limit of risk ratio.
//'
//' * \code{riskRatioUpper}: The upper equivalence limit of risk ratio.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{riskRatio}: The risk ratio.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskRatioLower}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   lower equivalence limit.
//'
//' * \code{zstatRiskRatioUpper}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   upper equivalence limit.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' powerRiskRatioExactEquiv(
//'   n = 200, riskRatioLower = 0.8, riskRatioUpper = 1.25,
//'   pi1 = 0.775, pi2 = 0.775, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame powerRiskRatioExactEquiv(
    const int n = NA_INTEGER,
    const double riskRatioLower = NA_REAL,
    const double riskRatioUpper = NA_REAL,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05,
    const bool calculateAttainedAlpha = true) {

  auto df = powerRiskRatioExactEquivcpp(
    n, riskRatioLower, riskRatioUpper, pi1, pi2, allocationRatioPlanned,
    alpha, calculateAttainedAlpha);

  return Rcpp::wrap(df);
}


DataFrameCpp samplesizeRiskRatioExactEquivcpp(
    const double beta,
    const double riskRatioLower,
    const double riskRatioUpper,
    const double pi1,
    const double pi2,
    const double allocationRatioPlanned,
    const double alpha,
    const bool calculateAttainedAlpha,
    const int max_n_search,
    const int window) {

  // allocation ratio r for group1 (fraction of total sample)
  double r = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  double ratio = pi1 / pi2;
  double theta = std::min(std::log(ratio) - std::log(riskRatioLower),
                          std::log(riskRatioUpper) - std::log(ratio));
  if (theta <= 0.0) throw std::invalid_argument(
    "pi1 / pi2 is outside or on the equivalence margins; sample size is not finite");

  double v1 = (1.0 - pi1) / (r * pi1) + (1.0 - pi2) / ((1.0 - r) * pi2);
  double z0 = boost_qnorm(1.0 - alpha);
  double z1 = boost_qnorm(1.0 - beta);

  // Normal-approximate continuous sample size
  double n0 = sq(z0 + z1) * v1 / (theta * theta);

  if (n0 > 0.5 * max_n_search) throw std::invalid_argument(
      std::string("Initial sample size estimate (") +
        std::to_string(static_cast<int>(std::floor(n0))) +
        ") is too large for exact test; consider normal approximation.");

  // Choose starting integer n
  int start_n = static_cast<int>(std::floor(n0));
  if (start_n < 1) start_n = 1;

  const double target = 1.0 - beta;

  // cache power evaluations
  std::unordered_map<int, double> power_cache;
  power_cache.reserve(1024);

  auto power_at = [&](int n)->double {
    if (n <= 0) return 0.0;
    auto it = power_cache.find(n);
    if (it != power_cache.end()) return it->second;
    DataFrameCpp df = powerRiskRatioExactEquivcpp(
      n, riskRatioLower, riskRatioUpper, pi1, pi2, allocationRatioPlanned,
      alpha, false);
    double p = df.get<double>("power")[0];
    power_cache.emplace(n, p);
    return p;
  };

  // Step-down: if start_n's power >= target, decrease until power < target
  int n = start_n;
  double p_n = power_at(n);
  while (p_n >= target && n > 1) {
    --n;
    p_n = power_at(n);
  }

  // Linear upward scan from candidate to max_n_search, enforcing window criterion
  int found_n = -1;
  for (int nn = n + 1; nn <= max_n_search; ++nn) {
    double p_nn = power_at(nn);
    if (p_nn < target) continue;

    // verify next window-1 sizes are also >= target
    bool ok = true;
    int fail_at = -1;
    for (int j = 1; j < window; ++j) {
      double pj = power_at(nn + j);
      if (pj < target) { ok = false; fail_at = j; break; }
    }
    if (ok) { found_n = nn; break; }
    n += fail_at; // skip ahead by fail_at for next iteration
  }

  if (found_n == -1) {
    throw std::runtime_error("No sample size <= max_n_search satisfies the "
                               "windowed power criterion");
  }

  // Return final DataFrame
  DataFrameCpp final_df = powerRiskRatioExactEquivcpp(
    found_n, riskRatioLower, riskRatioUpper, pi1, pi2, allocationRatioPlanned,
    alpha, calculateAttainedAlpha);
  return final_df;
}


//' @title Sample Size for Exact Unconditional Test of Equivalence in Risk
//' Ratio
//' @description Obtains the sample size given power for exact unconditional
//' test of equivalence in risk ratio.
//'
//' @param beta The type II error.
//' @param riskRatioLower The lower equivalence limit of risk ratio.
//' @param riskRatioUpper The upper equivalence limit of risk ratio.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//' @param max_n_search The maximum sample size to search up to. If no
//'   sample size up to this value satisfies the windowed power criterion,
//'   an error is thrown.
//' @param window The number of consecutive sample sizes that must all
//'   satisfy the power criterion to confirm the found sample size.
//'   This is to mitigate non-monotonicity of power in sample size for
//'   the exact test.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified significance level for each of the two
//'   one-sided tests.
//'
//' * \code{attainedAlpha}: The attained significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskRatioLower}: The lower equivalence limit of risk ratio.
//'
//' * \code{riskRatioUpper}: The upper equivalence limit of risk ratio.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{riskRatio}: The risk ratio.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskRatioLower}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   lower equivalence limit.
//'
//' * \code{zstatRiskRatioUpper}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   upper equivalence limit.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' samplesizeRiskRatioExactEquiv(
//'   beta = 0.2, riskRatioLower = 0.7, riskRatioUpper = 1/0.7,
//'   pi1 = 0.95, pi2 = 0.95, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame samplesizeRiskRatioExactEquiv(
    const double beta = NA_REAL,
    const double riskRatioLower = NA_REAL,
    const double riskRatioUpper = NA_REAL,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05,
    const bool calculateAttainedAlpha = true,
    const int max_n_search = 1000,
    const int window = 10) {

  auto df = samplesizeRiskRatioExactEquivcpp(
    beta, riskRatioLower, riskRatioUpper, pi1, pi2, allocationRatioPlanned,
    alpha, calculateAttainedAlpha, max_n_search, window);

  return Rcpp::wrap(df);
}


DataFrameCpp riskDiffExactPValuecpp(
    const int n1,
    const int y1,
    const int n2,
    const int y2,
    const double riskDiffH0,
    bool directionUpper) {

  const int k = (n1 + 1) * (n2 + 1);

  // Precompute T for all tables (row-major: index = y1*(n2+1) + y2)
  std::vector<double> T; T.reserve(k);
  for (int x = 0; x <= n1; ++x) {
    for (int y = 0; y <= n2; ++y) {
      auto pp = remlRiskDiff(n1, x, n2, y, riskDiffH0);
      double p1hat = static_cast<double>(x) / static_cast<double>(n1);
      double p2hat = static_cast<double>(y) / static_cast<double>(n2);
      double md = p1hat - p2hat - riskDiffH0;
      double mv = pp[0] * (1.0 - pp[0]) / static_cast<double>(n1) +
        pp[1] * (1.0 - pp[1]) / static_cast<double>(n2);
      if (mv < 1e-8) mv = 1e-8;
      T.push_back(md / std::sqrt(mv));
    }
  }

  // observed z-stat and observed difference
  double Tobs = T[y1 * (n2 + 1) + y2];
  double riskDiff = static_cast<double>(y1) / static_cast<double>(n1) -
    static_cast<double>(y2) / static_cast<double>(n2);

  // precompute log-choose tables for speed
  std::vector<double> logc_n1 = make_log_choose(n1);
  std::vector<double> logc_n2 = make_log_choose(n2);

  // f(p2) returns -s where s is the total probability under
  // H0 (p1 = p2 + riskDiffH0) of all tables with statistic as or
  // more extreme than observed according to directionUpper.
  auto f = [&](double p2)->double {
    double p1 = p2 + riskDiffH0;

    // compute binomial pmfs for margins
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int x = 0; x <= n1; ++x)
      q1[x] = binomial_pmf_from_logchoose(x, n1, p1, logc_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2, logc_n2);

    double s = 0.0;
    const int sign = (directionUpper ? +1 : -1);
    // iterate all flat indices and sum q1[y1] * q2[y2] where
    // (T[i] - Tobs) has same sign or zero
    int idx = 0;
    for (int yy1 = 0; yy1 <= n1; ++yy1) {
      double q1val = q1[yy1];
      for (int yy2 = 0; yy2 <= n2; ++yy2) {
        if (sign * (T[idx] - Tobs) >= 0.0) {
          s += q1val * q2[yy2];
        }
        ++idx;
      }
    }
    return -s;
  };

  // search domain for p2
  double pi2lower = std::max(0.0, -riskDiffH0);
  double pi2upper = std::min(1.0, 1.0 - riskDiffH0);

  if (pi2upper - pi2lower <= 0.0) {
    // degenerate case: p2 is fixed under H0; just compute p-value directly
    double pvalue = -f(pi2lower);
    DataFrameCpp out;
    out.push_back(riskDiffH0, "riskDiffH0");
    out.push_back(directionUpper, "directionUpper");
    out.push_back(riskDiff, "riskDiff");
    out.push_back(Tobs, "zstat");
    out.push_back(pvalue, "pvalue");
    out.push_back(pi2lower, "pi2star");
    return out;
  }

  const int K = 100;
  std::vector<double> a(K), b(K);

  // Partition domain into K subintervals; run local minimizer on each
  double delta = (pi2upper - pi2lower) / static_cast<double>(K);

  for (int ii = 0; ii < K; ++ii) {
    double lo = pi2lower + ii * delta;
    double hi = lo + delta;
    auto res = mini(f, lo, hi);
    a[ii] = res.first;
    b[ii] = res.second;
  }

  // pick best (minimum b)
  int best = 0;
  for (int ii = 1; ii < K; ++ii) if (b[ii] < b[best]) best = ii;

  double pi2star = a[best];
  double pvalue = -b[best];

  DataFrameCpp out;
  out.push_back(riskDiffH0, "riskDiffH0");
  out.push_back(directionUpper, "directionUpper");
  out.push_back(riskDiff, "riskDiff");
  out.push_back(Tobs, "zstat");
  out.push_back(pvalue, "pvalue");
  out.push_back(pi2star, "pi2star");
  return out;
}


//' @title P-Value for Exact Unconditional Test of Risk Difference
//' @description Obtains the p-value for exact unconditional
//' test of risk difference.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param riskDiffH0 The risk difference under the null hypothesis.
//'   Defaults to 0.
//' @param directionUpper Whether larger values represent better
//'   responses.
//'
//' @return A data frame containing the following variables:
//'
//' * \code{riskDiffH0}: The risk difference under the null hypothesis.
//'
//' * \code{directionUpper}: Whether larger values represent better
//'   responses.
//'
//' * \code{riskDiff}: The observed risk difference.
//'
//' * \code{zstat}: The observed value of the Z test statistic.
//'
//' * \code{pvalue}: The one-sided p-value for the unconditional exact test.
//'
//' * \code{pi2star}: The value of pi2 that yields the p-value.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' riskDiffExactPValue(riskDiffH0 = 0, directionUpper = 1,
//'                     n1 = 68, y1 = 2, n2 = 65, y2 = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame riskDiffExactPValue(
    const int n1 = NA_INTEGER,
    const int y1 = NA_INTEGER,
    const int n2 = NA_INTEGER,
    const int y2 = NA_INTEGER,
    const double riskDiffH0 = 0,
    bool directionUpper = true) {

  auto df = riskDiffExactPValuecpp(n1, y1, n2, y2, riskDiffH0, directionUpper);
  return Rcpp::wrap(df);
}


DataFrameCpp riskDiffExactCIcpp(
    const int n1,
    const int y1,
    const int n2,
    const int y2,
    const double cilevel) {

  // estimate (risk difference)
  double p1 = static_cast<double>(y1) / static_cast<double>(n1);
  double p2 = static_cast<double>(y2) / static_cast<double>(n2);
  double estimate = p1 - p2;

  double alpha = 1.0 - cilevel; // two-sided alpha
  const double target = alpha / 2.0;

  auto f1 = [&](double riskDiff)->double {
    auto df = riskDiffExactPValuecpp(n1, y1, n2, y2, riskDiff, true);
    double pvalue = df.get<double>("pvalue")[0];
    return pvalue - target;
  };

  auto f2 = [&](double riskDiff)->double {
    auto df = riskDiffExactPValuecpp(n1, y1, n2, y2, riskDiff, false);
    double pvalue = df.get<double>("pvalue")[0];
    return pvalue - target;
  };

  // Compute lower and upper with brent
  double lower = brent(f1, -1.0, estimate, 1e-6);
  double upper = brent(f2, estimate, 1.0, 1e-6);

  DataFrameCpp result;
  result.push_back(std::string("risk difference"), "scale");
  result.push_back(estimate, "estimate");
  result.push_back(lower, "lower");
  result.push_back(upper, "upper");
  result.push_back(cilevel, "cilevel");
  return result;
}


//' @title Exact Unconditional Confidence Interval for Risk Difference
//' @description Obtains the exact unconditional confidence interval for
//' risk difference based on the standardized score statistic.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @return A data frame containing the following variables:
//'
//' * \code{scale}: The scale of treatment effect.
//'
//' * \code{estimate}: The point estimate.
//'
//' * \code{lower}: The lower limit of the confidence interval.
//'
//' * \code{upper}: The upper limit of the confidence interval.
//'
//' * \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' riskDiffExactCI(n1 = 30, y1 = 2, n2 = 30, y2 = 1, cilevel = 0.95)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame riskDiffExactCI(
    const int n1 = NA_INTEGER,
    const int y1 = NA_INTEGER,
    const int n2 = NA_INTEGER,
    const int y2 = NA_INTEGER,
    const double cilevel = 0.95) {

  auto df = riskDiffExactCIcpp(n1, y1, n2, y2, cilevel);
  return Rcpp::wrap(df);
}

DataFrameCpp riskRatioExactPValuecpp(
    const int n1,
    const int y1,
    const int n2,
    const int y2,
    const double riskRatioH0,
    bool directionUpper) {

  const int k = (n1 + 1) * (n2 + 1);

  // Precompute T for all tables (row-major: index = y1*(n2+1) + y2)
  std::vector<double> T; T.reserve(k);
  for (int x = 0; x <= n1; ++x) {
    for (int y = 0; y <= n2; ++y) {
      auto pp = remlRiskRatio(n1, x, n2, y, riskRatioH0);
      double p1hat = static_cast<double>(x) / static_cast<double>(n1);
      double p2hat = static_cast<double>(y) / static_cast<double>(n2);
      double md = p1hat - p2hat * riskRatioH0;
      double mv = pp[0] * (1.0 - pp[0]) / static_cast<double>(n1) +
        sq(riskRatioH0) * pp[1] * (1.0 - pp[1]) / static_cast<double>(n2);
      if (mv < 1e-8) mv = 1e-8;
      T.push_back(md / std::sqrt(mv));
    }
  }

  // observed z-stat and observed ratio
  double Tobs = T[y1 * (n2 + 1) + y2];
  double p1hat = static_cast<double>(y1) / static_cast<double>(n1);
  double p2hat = static_cast<double>(y2) / static_cast<double>(n2);
  double riskRatio = p2hat > 0.0 ? (p1hat / p2hat) : POS_INF;

  // precompute log-choose tables for speed
  std::vector<double> logc_n1 = make_log_choose(n1);
  std::vector<double> logc_n2 = make_log_choose(n2);

  // f(p2) returns -s where s is the total probability under
  // H0 (p1 = p2 * riskRatioH0) of all tables with statistic as or
  // more extreme than observed according to directionUpper.
  auto f = [&](double p2)->double {
    double p1 = p2 * riskRatioH0;

    // compute binomial pmfs for margins
    std::vector<double> q1(n1 + 1), q2(n2 + 1);
    for (int x = 0; x <= n1; ++x)
      q1[x] = binomial_pmf_from_logchoose(x, n1, p1, logc_n1);
    for (int y = 0; y <= n2; ++y)
      q2[y] = binomial_pmf_from_logchoose(y, n2, p2, logc_n2);

    double s = 0.0;
    const int sign = (directionUpper ? +1 : -1);
    // iterate all flat indices and sum q1[y1] * q2[y2] where
    // (T[i] - Tobs) has same sign or zero
    int idx = 0;
    for (int yy1 = 0; yy1 <= n1; ++yy1) {
      double q1val = q1[yy1];
      for (int yy2 = 0; yy2 <= n2; ++yy2) {
        if (sign * (T[idx] - Tobs) >= 0.0) {
          s += q1val * q2[yy2];
        }
        ++idx;
      }
    }
    return -s;
  };

  // search domain for p2
  double pi2lower = 0.0;
  double pi2upper = std::min(1.0, 1.0 / riskRatioH0);

  if (pi2upper - pi2lower <= 0.0) {
    // degenerate case: p2 is fixed under H0; just compute p-value directly
    double pvalue = -f(pi2lower);
    DataFrameCpp out;
    out.push_back(riskRatioH0, "riskRatioH0");
    out.push_back(directionUpper, "directionUpper");
    out.push_back(riskRatio, "riskRatio");
    out.push_back(Tobs, "zstat");
    out.push_back(pvalue, "pvalue");
    out.push_back(pi2lower, "pi2star");
    return out;
  }

  const int K = 100;
  std::vector<double> a(K), b(K);

  // Partition domain into K subintervals; run local minimizer on each
  double delta = (pi2upper - pi2lower) / static_cast<double>(K);

  for (int ii = 0; ii < K; ++ii) {
    double lo = pi2lower + ii * delta;
    double hi = lo + delta;
    auto res = mini(f, lo, hi);
    a[ii] = res.first;
    b[ii] = res.second;
  }

  // pick best (minimum b)
  int best = 0;
  for (int ii = 1; ii < K; ++ii) if (b[ii] < b[best]) best = ii;

  double pi2star = a[best];
  double pvalue = -b[best];

  DataFrameCpp out;
  out.push_back(riskRatioH0, "riskRatioH0");
  out.push_back(directionUpper, "directionUpper");
  out.push_back(riskRatio, "riskRatio");
  out.push_back(Tobs, "zstat");
  out.push_back(pvalue, "pvalue");
  out.push_back(pi2star, "pi2star");
  return out;
}


//' @title P-Value for Exact Unconditional Test of Risk Ratio
//' @description Obtains the p-value for exact unconditional
//' test of risk ratio.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param riskRatioH0 The risk ratio under the null hypothesis.
//'   Defaults to 1.
//' @param directionUpper Whether larger values represent better
//'   responses.
//'
//' @return A data frame containing the following variables:
//'
//' * \code{riskRatioH0}: The risk ratio under the null hypothesis.
//'
//' * \code{directionUpper}: Whether larger values represent better
//'   responses.
//'
//' * \code{riskRatio}: The observed risk ratio.
//'
//' * \code{zstat}: The observed value of the Z test statistic.
//'
//' * \code{pvalue}: The one-sided p-value for the unconditional exact test.
//'
//' * \code{pi2star}: The value of pi2 that yields the p-value.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' riskRatioExactPValue(riskRatioH0 = 1, directionUpper = 1,
//'                      n1 = 68, y1 = 2, n2 = 65, y2 = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame riskRatioExactPValue(
    const int n1 = NA_INTEGER,
    const int y1 = NA_INTEGER,
    const int n2 = NA_INTEGER,
    const int y2 = NA_INTEGER,
    const double riskRatioH0 = 1,
    bool directionUpper = true) {

  auto df = riskRatioExactPValuecpp(n1, y1, n2, y2, riskRatioH0, directionUpper);
  return Rcpp::wrap(df);
}


DataFrameCpp riskRatioExactCIcpp(
    const int n1,
    const int y1,
    const int n2,
    const int y2,
    const double cilevel) {

  // estimate (risk difference)
  double p1 = static_cast<double>(y1) / static_cast<double>(n1);
  double p2 = static_cast<double>(y2) / static_cast<double>(n2);
  double estimate = p2 > 0.0 ? p1 / p2 : POS_INF;

  double alpha = 1.0 - cilevel; // two-sided alpha
  const double target = alpha / 2.0;

  auto f1 = [&](double riskRatio)->double {
    auto df = riskRatioExactPValuecpp(n1, y1, n2, y2, riskRatio, true);
    double pvalue = df.get<double>("pvalue")[0];
    return pvalue - target;
  };

  auto f2 = [&](double riskRatio)->double {
    auto df = riskRatioExactPValuecpp(n1, y1, n2, y2, riskRatio, false);
    double pvalue = df.get<double>("pvalue")[0];
    return pvalue - target;
  };

  // Compute lower and upper with brent
  double lower = brent(f1, 0.0001, estimate, 1e-6);
  double upper = brent(f2, estimate, 1000.0, 1e-6);

  DataFrameCpp result;
  result.push_back(std::string("risk ratio"), "scale");
  result.push_back(estimate, "estimate");
  result.push_back(lower, "lower");
  result.push_back(upper, "upper");
  result.push_back(cilevel, "cilevel");
  return result;
}


//' @title Exact Unconditional Confidence Interval for Risk Ratio
//' @description Obtains the exact unconditional confidence interval for
//' risk ratio based on the standardized score statistic.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @return A data frame containing the following variables:
//'
//' * \code{scale}: The scale of treatment effect.
//'
//' * \code{estimate}: The point estimate.
//'
//' * \code{lower}: The lower limit of the confidence interval.
//'
//' * \code{upper}: The upper limit of the confidence interval.
//'
//' * \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' riskRatioExactCI(n1 = 30, y1 = 2, n2 = 30, y2 = 1, cilevel = 0.95)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame riskRatioExactCI(
    const int n1 = NA_INTEGER,
    const int y1 = NA_INTEGER,
    const int n2 = NA_INTEGER,
    const int y2 = NA_INTEGER,
    const double cilevel = 0.95) {

  auto df = riskRatioExactCIcpp(n1, y1, n2, y2, cilevel);
  return Rcpp::wrap(df);
}
