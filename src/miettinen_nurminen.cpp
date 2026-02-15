#include "miettinen_nurminen.h"
#include "utilities.h"
#include "dataframe_list.h"

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


//' @title REML Estimates of Individual Proportions With Specified Risk
//' difference
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified risk difference.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param riskDiffH0 The specified risk difference.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the response probabilities for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlRiskDiff(n1 = 10, y1 = 4, n2 = 20, y2 = 0, riskDiffH0 = 0.1)
//'
//' @export
//'
// [[Rcpp::export]]
std::vector<double> remlRiskDiff(const double n1,
                                 const double y1,
                                 const double n2,
                                 const double y2,
                                 const double riskDiffH0 = 0.0) {
  double n = n1 + n2;
  double y = y1 + y2;

  constexpr double EPS = 1e-12;
  double p1, p2;

  // If H0 difference is (effectively) zero, pooled estimate
  if (std::fabs(riskDiffH0) < 1e-8) {
    p1 = y / n;
    p2 = p1;
  } else {
    // Coefficients used in the cubic-solving approach
    double L3 = n;
    double L2 = (n1 + 2.0 * n2) * riskDiffH0 - n - y;
    double L1 = (n2 * riskDiffH0 - n - 2.0 * y2) * riskDiffH0 + y;
    double L0 = y2 * riskDiffH0 * (1.0 - riskDiffH0);

    // Compute q and p (Cardano-style)
    const double denom3 = 3.0 * L3;
    const double denom3_sq = denom3 * denom3;
    const double denom3_cu = denom3_sq * denom3;

    double q = (L2 * L2 * L2) / denom3_cu -
      (L1 * L2) / (6.0 * (L3 * L3)) + L0 / (2.0 * L3);

    double inside = (L2 * L2) / denom3_sq - L1 / denom3;

    // If inside is negative (numerical error), clamp to zero
    if (inside < 0.0) inside = 0.0;
    double sign = (q > 0.0) ? 1.0 : -1.0;
    double p = sign * std::sqrt(inside);

    if (std::fabs(p) < EPS) {
      // If p is effectively zero, we have a repeated root
      p2 = -L2 / denom3;
    } else {
      // ratio passed to acos
      double ratio = q / (p * p * p);

      // clamp ratio to [-1,1] to avoid NaN from acos due to small numerical overshoot
      ratio = std::max(-1.0, std::min(1.0, ratio));

      const double PI = std::acos(-1.0);
      double a = (PI + std::acos(ratio)) / 3.0;

      p2 = 2.0 * p * std::cos(a) - L2 / denom3;
    }

    p1 = p2 + riskDiffH0;
  }

  std::vector<double> result = {p1, p2};
  return result;
}


//' @title Miettinen-Nurminen Score Test Statistic for Two-Sample Risk
//' difference
//' @description Obtains the Miettinen-Nurminen score test statistic
//' for two-sample risk difference possibly with stratification.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param riskDiffH0 The risk difference under the null hypothesis.
//'   Defaults to 0.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatRiskDiff(n1 = c(10, 10), y1 = c(4, 3),
//'               n2 = c(20, 10), y2 = c(2, 0), riskDiffH0 = 0)
//'
//' @export
//'
// [[Rcpp::export]]
double zstatRiskDiff(const std::vector<double>& n1,
                     const std::vector<double>& y1,
                     const std::vector<double>& n2,
                     const std::vector<double>& y2,
                     const double riskDiffH0 = 0.0) {

  double num = 0.0, den = 0.0;
  const std::size_t k = n1.size();
  for (std::size_t i = 0; i < k; ++i) {
    double ni = n1[i] + n2[i];
    double wi = n1[i] * n2[i] / ni;
    double mdi = y1[i] / n1[i] - y2[i] / n2[i] - riskDiffH0;
    auto ppi = remlRiskDiff(n1[i], y1[i], n2[i], y2[i], riskDiffH0);
    double p1i = ppi[0];
    double p2i = ppi[1];
    double mvi = p1i * (1.0 - p1i) / n1[i] + p2i * (1.0 - p2i) / n2[i];
    mvi = std::max(mvi * ni / (ni - 1.0), 1e-8);
    num += wi * mdi;
    den += wi * wi * mvi;
  }

  return num / std::sqrt(den);
}


DataFrameCpp mnRiskDiffCIcpp(const std::vector<double>& n1,
                             const std::vector<double>& y1,
                             const std::vector<double>& n2,
                             const std::vector<double>& y2,
                             const double cilevel = 0.95) {
  // Input validation
  const std::size_t k = n1.size();
  if (y1.size() != k || n2.size() != k || y2.size() != k) {
    throw std::invalid_argument("All input vectors must have the same length");
  }
  for (std::size_t i = 0; i < k; ++i) {
    if (!(n1[i] > 0.0)) throw std::invalid_argument("n1 must be positive");
    if (!(n2[i] > 0.0)) throw std::invalid_argument("n2 must be positive");
    if (!(y1[i] >= 0.0 && y1[i] <= n1[i]))
      throw std::invalid_argument("y1 must be between 0 and n1");
    if (!(y2[i] >= 0.0 && y2[i] <= n2[i]))
      throw std::invalid_argument("y2 must be between 0 and n2");
  }
  if (!(cilevel > 0.0 && cilevel < 1.0))
    throw std::invalid_argument("cilevel must lie between 0 and 1");


  double sum_w = 0.0, estimate = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    double ni = n1[i] + n2[i];
    double wi = n1[i] * n2[i] / ni;
    double mdi = y1[i] / n1[i] - y2[i] / n2[i];
    sum_w += wi;
    estimate += wi * mdi;
  }
  estimate /= sum_w;

  // Normal quantile for the confidence level
  double b = boost_qnorm((1.0 + cilevel) / 2.0);

  auto f1 = [&](double d) { return zstatRiskDiff(n1, y1, n2, y2, d) - b; };
  auto f2 = [&](double d) { return zstatRiskDiff(n1, y1, n2, y2, d) + b; };
  double lower = brent(f1, -1.0, estimate, 1e-6);
  double upper = brent(f2, estimate, 1.0, 1e-6);

  DataFrameCpp df;
  df.push_back(std::string("risk difference"), "scale");
  df.push_back(estimate, "estimate");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");
  df.push_back(cilevel, "cilevel");

  return df;
}



//' @title Miettinen-Nurminen Score Confidence Interval for
//' Two-Sample Risk Difference
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample risk difference possibly with
//' stratification.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input sample size
//'   and number of responses for each treatment group.
//'   It has the following variables:
//'
//'     - \code{n1}: The sample size for the active treatment group.
//'
//'     - \code{y1}: The number of responses for the active treatment group.
//'
//'     - \code{n2}: The sample size for the control group.
//'
//'     - \code{y2}: The number of responses for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for risk difference. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnRiskDiffCI(n1 = c(10, 10), y1 = c(4, 3),
//'              n2 = c(20, 10), y2 = c(2, 0))
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::DataFrame mnRiskDiffCI(const std::vector<double>& n1,
                             const std::vector<double>& y1,
                             const std::vector<double>& n2,
                             const std::vector<double>& y2,
                             const double cilevel = 0.95) {
  auto df = mnRiskDiffCIcpp(n1, y1, n2, y2, cilevel);
  return Rcpp::wrap(df);
}



//' @title REML Estimates of Individual Proportions With Specified Risk
//' Ratio
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified risk ratio.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param riskRatioH0 The specified risk ratio.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the response probabilities for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlRiskRatio(n1 = 10, y1 = 4, n2 = 20, y2 = 2, riskRatioH0 = 1.2)
//'
//' @export
//'
// [[Rcpp::export]]
std::vector<double> remlRiskRatio(const double n1,
                                  const double y1,
                                  const double n2,
                                  const double y2,
                                  const double riskRatioH0 = 1.0) {
  double n = n1 + n2;
  double y = y1 + y2;

  double p1, p2;
  if (std::fabs(riskRatioH0 - 1.0) < 1e-8) {
    p1 = y / n;
    p2 = p1;
  } else {
    double a = n * riskRatioH0;
    double b = -(n1 * riskRatioH0 + y1 + n2 + y2 * riskRatioH0);
    double c = y;

    double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
      // If the discriminant is negative (numerical error), clamp to zero
      disc = 0.0;
    }

    p2 = (-b - std::sqrt(disc)) / (2.0 * a);
    p1 = p2 * riskRatioH0;
  }

  std::vector<double> result = {p1, p2};
  return result;
}


//' @title Miettinen-Nurminen Score Test Statistic for Two-Sample Risk Ratio
//' @description Obtains the Miettinen-Nurminen score test statistic for
//' two-sample risk ratio possibly with stratification.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param riskRatioH0 The risk ratio under the null hypothesis.
//'   Defaults to 1.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatRiskRatio(n1 = c(10, 10), y1 = c(4, 3),
//'                n2 = c(20, 10), y2 = c(2, 0), riskRatioH0 = 1)
//'
//' @export
//'
// [[Rcpp::export]]
double zstatRiskRatio(const std::vector<double>& n1,
                      const std::vector<double>& y1,
                      const std::vector<double>& n2,
                      const std::vector<double>& y2,
                      const double riskRatioH0 = 1.0) {

  double num = 0.0, den = 0.0;
  const std::size_t k = n1.size();
  for (std::size_t i = 0; i < k; ++i) {
    double ni = n1[i] + n2[i];
    double wi = n1[i] * n2[i] / ni;
    double mdi = y1[i] / n1[i] - y2[i] / n2[i] * riskRatioH0;
    auto ppi = remlRiskRatio(n1[i], y1[i], n2[i], y2[i], riskRatioH0);
    double p1i = ppi[0];
    double p2i = ppi[1];
    double mvi = p1i * (1.0 - p1i) / n1[i] +
      p2i * (1.0 - p2i) / n2[i] * riskRatioH0 * riskRatioH0;
    mvi = std::max(mvi * ni / (ni - 1.0), 1e-8);
    num += wi * mdi;
    den += wi * wi * mvi;
  }

  return num / std::sqrt(den);
}


DataFrameCpp mnRiskRatioCIcpp(const std::vector<double>& n1,
                              const std::vector<double>& y1,
                              const std::vector<double>& n2,
                              const std::vector<double>& y2,
                              const double cilevel = 0.95) {
  // Input validation
  const std::size_t k = n1.size();
  if (y1.size() != k || n2.size() != k || y2.size() != k) {
    throw std::invalid_argument("All input vectors must have the same length");
  }
  for (std::size_t i = 0; i < k; ++i) {
    if (!(n1[i] > 0.0)) throw std::invalid_argument("n1 must be positive");
    if (!(n2[i] > 0.0)) throw std::invalid_argument("n2 must be positive");
    if (!(y1[i] >= 0.0 && y1[i] <= n1[i]))
      throw std::invalid_argument("y1 must be between 0 and n1");
    if (!(y2[i] >= 0.0 && y2[i] <= n2[i]))
      throw std::invalid_argument("y2 must be between 0 and n2");
  }
  if (!(cilevel > 0.0 && cilevel < 1.0))
    throw std::invalid_argument("cilevel must lie between 0 and 1");

  double sum_w = 0.0, p1 = 0.0, p2 = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    double ni = n1[i] + n2[i];
    double wi = (n1[i] * n2[i]) / ni;
    sum_w += wi;
    p1 += wi * (y1[i] / n1[i]);
    p2 += wi * (y2[i] / n2[i]);
  }
  p1 /= sum_w;
  p2 /= sum_w;

  double b = boost_qnorm((1.0 + cilevel) / 2.0);

  // special cases
  bool any_both_zero = false;
  for (std::size_t i = 0; i < k; ++i) {
    if (y1[i] == 0 && y2[i] == 0) {
      any_both_zero = true;
      break;
    }
  }

  bool all_y1_zero = true;
  for (std::size_t i = 0; i < k; ++i) {
    if (y1[i] != 0) {
      all_y1_zero = false;
      break;
    }
  }

  bool all_y2_zero = true;
  for (std::size_t i = 0; i < k; ++i) {
    if (y2[i] != 0) {
      all_y2_zero = false;
      break;
    }
  }

  double estimate, lower, upper;
  if (any_both_zero) {
    // If any stratum has zero events in both groups, the risk ratio is undefined
    estimate = NaN;
    lower = NaN;
    upper = NaN;
  } else if (all_y1_zero) {
    // If all events in the active group are zero, the risk ratio is zero and
    // the upper limit is finite
    estimate = 0;
    lower = 0;
    auto f2 = [&](double r) { return zstatRiskRatio(n1, y1, n2, y2, r) + b; };
    upper = brent(f2, 0.001, 1000.0, 1e-6);
  } else if (all_y2_zero) {
    // If all events in the control group are zero, the risk ratio is infinite and
    // the lower limit is finite
    estimate = POS_INF;
    upper = POS_INF;
    auto f1 = [&](double r) { return zstatRiskRatio(n1, y1, n2, y2, r) - b; };
    lower = brent(f1, 0.001, 1000.0, 1e-6);
  } else {
    estimate = p1 / p2;

    auto f1 = [&](double r) { return zstatRiskRatio(n1, y1, n2, y2, r) - b; };
    auto f2 = [&](double r) { return zstatRiskRatio(n1, y1, n2, y2, r) + b; };

    lower = brent(f1, 0.001, estimate, 1e-6);
    upper = brent(f2, estimate, 1000.0, 1e-6);
  }

  DataFrameCpp df;
  df.push_back(std::string("risk ratio"), "scale");
  df.push_back(estimate, "estimate");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");
  df.push_back(cilevel, "cilevel");

  return df;
}


//' @title Miettinen-Nurminen Score Confidence Interval for
//' Two-Sample Risk Ratio
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample risk ratio possibly with
//' stratification.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input sample size
//'   and number of responses for each treatment group.
//'   It has the following variables:
//'
//'     - \code{n1}: The sample size for the active treatment group.
//'
//'     - \code{y1}: The number of responses for the active treatment group.
//'
//'     - \code{n2}: The sample size for the control group.
//'
//'     - \code{y2}: The number of responses for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for risk ratio. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnRiskRatioCI(n1 = c(10, 10), y1 = c(4, 3),
//'               n2 = c(20, 10), y2 = c(2, 0))
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::DataFrame mnRiskRatioCI(const std::vector<double>& n1,
                              const std::vector<double>& y1,
                              const std::vector<double>& n2,
                              const std::vector<double>& y2,
                              const double cilevel = 0.95) {

  auto df = mnRiskRatioCIcpp(n1, y1, n2, y2, cilevel);
  return Rcpp::wrap(df);
}


//' @title REML Estimates of Individual Proportions With Specified Odds
//' Ratio
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified odds ratio.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param oddsRatioH0 The specified odds ratio.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the response probabilities for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlOddsRatio(n1 = 10, y1 = 4, n2 = 20, y2 = 2, oddsRatioH0 = 1.25)
//'
//' @export
//'
// [[Rcpp::export]]
std::vector<double> remlOddsRatio(const double n1,
                                  const double y1,
                                  const double n2,
                                  const double y2,
                                  const double oddsRatioH0 = 1.0) {
  double n = n1 + n2;
  double y = y1 + y2;

  double p1, p2;
  if (std::fabs(oddsRatioH0 - 1.0) < 1e-8) {
    p1 = y / n;
    p2 = p1;
  } else {
    double a = n2 * (oddsRatioH0 - 1.0);
    double b = n1 * oddsRatioH0 + n2 - y * (oddsRatioH0 - 1.0);
    double c = -y;

    double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
      // If the discriminant is negative (numerical error), clamp to zero
      disc = 0.0;
    }

    p2 = (-b + std::sqrt(disc)) / (2.0 * a);
    p1 = p2 * oddsRatioH0 / (1.0 + p2 * (oddsRatioH0 - 1.0));
  }

  std::vector<double> result = {p1, p2};
  return result;
}



//' @title Miettinen-Nurminen Score Test Statistic for Two-Sample Odds Ratio
//' @description Obtains the Miettinen-Nurminen score test statistic for
//' two-sample odds ratio possibly with stratification.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param oddsRatioH0 The odds ratio under the null hypothesis.
//'   Defaults to 1.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatOddsRatio(n1 = c(10, 10), y1 = c(4, 3),
//'                n2 = c(20, 10), y2 = c(2, 0), oddsRatioH0 = 1)
//'
//' @export
//'
// [[Rcpp::export]]
double zstatOddsRatio(const std::vector<double>& n1,
                      const std::vector<double>& y1,
                      const std::vector<double>& n2,
                      const std::vector<double>& y2,
                      const double oddsRatioH0 = 1.0) {

  double num = 0.0, den = 0.0;
  const std::size_t k = n1.size();
  for (std::size_t i = 0; i < k; ++i) {
    double ni = n1[i] + n2[i];
    double wi = n1[i] * n2[i] / ni;
    auto ppi = remlOddsRatio(n1[i], y1[i], n2[i], y2[i], oddsRatioH0);
    double p1i = ppi[0];
    double p2i = ppi[1];
    double v1i = p1i * (1.0 - p1i);
    double v2i = p2i * (1.0 - p2i);
    double mdi = (y1[i] / n1[i] - p1i) / v1i - (y2[i] / n2[i] - p2i) / v2i;
    double mvi = 1.0 / (n1[i] * v1i) + 1.0 / (n2[i] * v2i);
    mvi = std::max(mvi * ni / (ni - 1.0), 1e-8);
    num += wi * mdi;
    den += wi * wi * mvi;
  }

  return num / std::sqrt(den);
}


DataFrameCpp mnOddsRatioCIcpp(const std::vector<double>& n1,
                              const std::vector<double>& y1,
                              const std::vector<double>& n2,
                              const std::vector<double>& y2,
                              const double cilevel = 0.95) {

  // Input validation
  const std::size_t k = n1.size();
  if (y1.size() != k || n2.size() != k || y2.size() != k) {
    throw std::invalid_argument("All input vectors must have the same length");
  }
  for (std::size_t i = 0; i < k; ++i) {
    if (!(n1[i] > 0.0)) throw std::invalid_argument("n1 must be positive");
    if (!(n2[i] > 0.0)) throw std::invalid_argument("n2 must be positive");
    if (!(y1[i] >= 0.0 && y1[i] <= n1[i]))
      throw std::invalid_argument("y1 must be between 0 and n1");
    if (!(y2[i] >= 0.0 && y2[i] <= n2[i]))
      throw std::invalid_argument("y2 must be between 0 and n2");
  }
  if (!(cilevel > 0.0 && cilevel < 1.0))
    throw std::invalid_argument("cilevel must lie between 0 and 1");

  double b = boost_qnorm((1.0 + cilevel) / 2.0);

  // special cases
  bool any_both_zero = false;
  for (std::size_t i = 0; i < k; ++i) {
    if (y1[i] == 0 && y2[i] == 0) {
      any_both_zero = true;
      break;
    }
  }

  bool any_both_n = false;
  for (std::size_t i = 0; i < k; ++i) {
    if (y1[i] == n1[i] && y2[i] == n2[i]) {
      any_both_n = true;
      break;
    }
  }

  bool all_y1_zero = true;
  for (std::size_t i = 0; i < k; ++i) {
    if (y1[i] != 0) {
      all_y1_zero = false;
      break;
    }
  }

  bool all_y1_n = true;
  for (std::size_t i = 0; i < k; ++i) {
    if (y1[i] != n1[i]) {
      all_y1_n = false;
      break;
    }
  }

  bool all_y2_zero = true;
  for (std::size_t i = 0; i < k; ++i) {
    if (y2[i] != 0) {
      all_y2_zero = false;
      break;
    }
  }

  bool all_y2_n = true;
  for (std::size_t i = 0; i < k; ++i) {
    if (y2[i] != n2[i]) {
      all_y2_n = false;
      break;
    }
  }

  double estimate, lower, upper;
  if (any_both_zero || any_both_n) {
    estimate = NaN;
    lower = NaN;
    upper = NaN;
  } else if (all_y1_zero || all_y2_n) {
    estimate = 0;
    lower = 0;
    auto f2 = [&](double r) { return zstatOddsRatio(n1, y1, n2, y2, r) + b; };
    upper = brent(f2, 0.001, 1000.0, 1e-6);
  } else if (all_y2_zero || all_y1_n) {
    estimate = POS_INF;
    upper = POS_INF;
    auto f1 = [&](double r) { return zstatOddsRatio(n1, y1, n2, y2, r) - b; };
    lower = brent(f1, 0.001, 1000.0, 1e-6);
  } else {
    auto f0 = [&](double r) { return zstatOddsRatio(n1, y1, n2, y2, r); };
    auto f1 = [&](double r) { return zstatOddsRatio(n1, y1, n2, y2, r) - b; };
    auto f2 = [&](double r) { return zstatOddsRatio(n1, y1, n2, y2, r) + b; };
    estimate = brent(f0, 0.001, 1000.0, 1e-6);
    lower = brent(f1, 0.001, estimate, 1e-6);
    upper = brent(f2, estimate, 1000.0, 1e-6);
  }

  DataFrameCpp df;
  df.push_back(std::string("odds ratio"), "scale");
  df.push_back(estimate, "estimate");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");
  df.push_back(cilevel, "cilevel");

  return df;
}


//' @title Miettinen-Nurminen Score Confidence Interval for
//' Two-Sample Odds Ratio
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample odds ratio possibly with
//' stratification.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input sample size
//'   and number of responses for each treatment group.
//'   It has the following variables:
//'
//'     - \code{n1}: The sample size for the active treatment group.
//'
//'     - \code{y1}: The number of responses for the active treatment group.
//'
//'     - \code{n2}: The sample size for the control group.
//'
//'     - \code{y2}: The number of responses for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for odds ratio. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnOddsRatioCI(n1 = c(10,10), y1 = c(4,3), n2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::DataFrame mnOddsRatioCI(const std::vector<double>& n1,
                              const std::vector<double>& y1,
                              const std::vector<double>& n2,
                              const std::vector<double>& y2,
                              const double cilevel = 0.95) {

  auto df = mnOddsRatioCIcpp(n1, y1, n2, y2, cilevel);
  return Rcpp::wrap(df);
}


//' @title REML Estimates of Individual Rates With Specified Rate
//' Difference
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified rate difference.
//'
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//' @param rateDiffH0 The specified rate difference.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the incidence rates for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlRateDiff(t1 = 10, y1 = 4, t2 = 20, y2 = 2, rateDiffH0 = 0.1)
//'
//' @export
//'
// [[Rcpp::export]]
std::vector<double> remlRateDiff(const double t1,
                                 const double y1,
                                 const double t2,
                                 const double y2,
                                 const double rateDiffH0 = 0.0) {
  double t = t1 + t2;
  double y = y1 + y2;

  double r1, r2;
  if (std::fabs(rateDiffH0) < 1e-8) {
    r1 = y / t;
    r2 = r1;
  } else {
    double a = t;
    double b = t * rateDiffH0 - y;
    double c = -y2 * rateDiffH0;

    double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
      // If the discriminant is negative (numerical error), clamp to zero
      disc = 0.0;
    }

    r2 = (-b + std::sqrt(disc)) / (2.0 * a);
    r1 = r2 + rateDiffH0;
  }

  std::vector<double> result = {r1, r2};
  return result;
}


//' @title Miettinen-Nurminen Score Test Statistic for Two-Sample Rate
//' Difference
//' @description Obtains the Miettinen-Nurminen score test statistic for
//' two-sample rate difference possibly with stratification.
//'
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//' @param rateDiffH0 The rate difference under the null hypothesis.
//'   Defaults to 0.
//'
//' @details
//' The Mantel-Haenszel weights are used for stratified samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatRateDiff(t1 = c(10, 10), y1 = c(4, 3),
//'               t2 = c(20, 10), y2 = c(2, 0), rateDiffH0 = 0)
//'
//' @export
//'
// [[Rcpp::export]]
double zstatRateDiff(const std::vector<double>& t1,
                     const std::vector<double>& y1,
                     const std::vector<double>& t2,
                     const std::vector<double>& y2,
                     const double rateDiffH0 = 0.0) {

  double num = 0.0, den = 0.0;
  const std::size_t k = t1.size();
  for (std::size_t i = 0; i < k; ++i) {
    double ti = t1[i] + t2[i];
    double wi = t1[i] * t2[i] / ti;
    double mdi = y1[i] / t1[i] - y2[i] / t2[i] - rateDiffH0;
    auto rri = remlRateDiff(t1[i], y1[i], t2[i], y2[i], rateDiffH0);
    double r1i = rri[0];
    double r2i = rri[1];
    double mvi = r1i / t1[i] + r2i / t2[i];
    mvi = std::max(mvi, 1e-8);
    num += wi * mdi;
    den += wi * wi * mvi;
  }

  return num / std::sqrt(den);
}


DataFrameCpp mnRateDiffCIcpp(const std::vector<double>& t1,
                             const std::vector<double>& y1,
                             const std::vector<double>& t2,
                             const std::vector<double>& y2,
                             const double cilevel = 0.95) {

  // Input validation
  const std::size_t k = t1.size();
  if (y1.size() != k || t2.size() != k || y2.size() != k) {
    throw std::invalid_argument("All input vectors must have the same length");
  }
  for (std::size_t i = 0; i < k; ++i) {
    if (!(t1[i] > 0.0)) throw std::invalid_argument("t1 must be positive");
    if (!(t2[i] > 0.0)) throw std::invalid_argument("t2 must be positive");
    if (!(y1[i] >= 0.0))
      throw std::invalid_argument("y1 must be nonnegative");
    if (!(y2[i] >= 0.0))
      throw std::invalid_argument("y2 must be nonnegative");
  }
  if (!(cilevel > 0.0 && cilevel < 1.0))
    throw std::invalid_argument("cilevel must lie between 0 and 1");

  double sum_w = 0.0, estimate = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    double ti = t1[i] + t2[i];
    double wi = t1[i] * t2[i] / ti;
    double mdi = y1[i] / t1[i] - y2[i] / t2[i];
    sum_w += wi;
    estimate += wi * mdi;
  }
  estimate /= sum_w;

  double b = boost_qnorm((1.0 + cilevel) / 2.0);

  auto f1 = [&](double d) { return zstatRateDiff(t1, y1, t2, y2, d) - b; };
  auto f2 = [&](double d) { return zstatRateDiff(t1, y1, t2, y2, d) + b; };

  double lower = brent(f1, -1000.0, estimate, 1e-6);
  double upper = brent(f2, estimate, 1000.0, 1e-6);

  DataFrameCpp df;
  df.push_back(std::string("rate difference"), "scale");
  df.push_back(estimate, "estimate");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");
  df.push_back(cilevel, "cilevel");

  return df;
}


//' @title Miettinen-Nurminen Score Confidence Interval for
//' Two-Sample Rate Difference
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample rate difference possibly with
//' stratification.
//'
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel weights are used for stratified samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input exposure
//'   and number of events for each treatment group.
//'   It has the following variables:
//'
//'     - \code{t1}: The exposure for the active treatment group.
//'
//'     - \code{y1}: The number of events for the active treatment group.
//'
//'     - \code{t2}: The exposure for the control group.
//'
//'     - \code{y2}: The number of events for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for rate difference. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnRateDiffCI(t1 = c(10,10), y1 = c(4,3), t2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::DataFrame mnRateDiffCI(const std::vector<double>& t1,
                             const std::vector<double>& y1,
                             const std::vector<double>& t2,
                             const std::vector<double>& y2,
                             const double cilevel = 0.95) {

  auto df = mnRateDiffCIcpp(t1, y1, t2, y2, cilevel);
  return Rcpp::wrap(df);
}



//' @title REML Estimates of Individual Rates With Specified Rate Ratio
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified rate ratio.
//'
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//' @param rateRatioH0 The specified rate ratio.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the incidence rates for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlRateRatio(t1 = 10, y1 = 4, t2 = 20, y2 = 2, rateRatioH0 = 1.1)
//'
//' @export
//'
// [[Rcpp::export]]
std::vector<double> remlRateRatio(const double t1,
                                  const double y1,
                                  const double t2,
                                  const double y2,
                                  const double rateRatioH0 = 1.0) {

  double r2 = (y1 + y2) / (t1 * rateRatioH0 + t2);
  double r1 = r2 * rateRatioH0;

  std::vector<double> result = {r1, r2};
  return result;
}


//' @title Miettinen-Nurminen Score Test Statistic for Two-Sample Rate Ratio
//' @description Obtains the Miettinen-Nurminen score test statistic for
//' two-sample rate ratio possibly with stratification.
//'
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//' @param rateRatioH0 The rate ratio under the null hypothesis.
//'   Defaults to 1.
//'
//' @details
//' The Mantel-Haenszel weights are used for stratified samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatRateRatio(t1 = c(10, 10), y1 = c(4, 3),
//'                t2 = c(20, 10), y2 = c(2, 0), rateRatioH0 = 1)
//'
//' @export
//'
// [[Rcpp::export]]
double zstatRateRatio(const std::vector<double>& t1,
                      const std::vector<double>& y1,
                      const std::vector<double>& t2,
                      const std::vector<double>& y2,
                      const double rateRatioH0 = 1.0) {

  double num = 0.0, den = 0.0;
  const std::size_t k = t1.size();
  for (std::size_t i = 0; i < k; ++i) {
    double ti = t1[i] + t2[i];
    double wi = t1[i] * t2[i] / ti;
    double mdi = y1[i] / t1[i] - y2[i] / t2[i] * rateRatioH0;
    auto rri = remlRateRatio(t1[i], y1[i], t2[i], y2[i], rateRatioH0);
    double r1i = rri[0];
    double r2i = rri[1];
    double mvi = r1i / t1[i] + r2i / t2[i] * rateRatioH0 * rateRatioH0;
    mvi = std::max(mvi, 1e-8);
    num += wi * mdi;
    den += wi * wi * mvi;
  }

  return num / std::sqrt(den);
}


DataFrameCpp mnRateRatioCIcpp(const std::vector<double>& t1,
                              const std::vector<double>& y1,
                              const std::vector<double>& t2,
                              const std::vector<double>& y2,
                              const double cilevel = 0.95) {

  // Input validation
  const std::size_t k = t1.size();
  if (y1.size() != k || t2.size() != k || y2.size() != k) {
    throw std::invalid_argument("All input vectors must have the same length");
  }
  for (std::size_t i = 0; i < k; ++i) {
    if (!(t1[i] > 0.0)) throw std::invalid_argument("t1 must be positive");
    if (!(t2[i] > 0.0)) throw std::invalid_argument("t2 must be positive");
    if (!(y1[i] >= 0.0))
      throw std::invalid_argument("y1 must be nonnegative");
    if (!(y2[i] >= 0.0))
      throw std::invalid_argument("y2 must be nonnegative");
  }
  if (!(cilevel > 0.0 && cilevel < 1.0))
    throw std::invalid_argument("cilevel must lie between 0 and 1");

  double sum_w = 0.0, r1 = 0.0, r2 = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    double ti = t1[i] + t2[i];
    double wi = (t1[i] * t2[i]) / ti;
    sum_w += wi;
    r1 += wi * (y1[i] / t1[i]);
    r2 += wi * (y2[i] / t2[i]);
  }
  r1 /= sum_w;
  r2 /= sum_w;

  double b = boost_qnorm((1.0 + cilevel) / 2.0);

  // special cases
  bool any_both_zero = false;
  for (std::size_t i = 0; i < k; ++i) {
    if (y1[i] == 0 && y2[i] == 0) {
      any_both_zero = true;
      break;
    }
  }

  bool all_y1_zero = true;
  for (std::size_t i = 0; i < k; ++i) {
    if (y1[i] != 0) {
      all_y1_zero = false;
      break;
    }
  }

  bool all_y2_zero = true;
  for (std::size_t i = 0; i < k; ++i) {
    if (y2[i] != 0) {
      all_y2_zero = false;
      break;
    }
  }

  double estimate, lower, upper;
  if (any_both_zero) {
    estimate = NaN;
    lower = NaN;
    upper = NaN;
  } else if (all_y1_zero) {
    estimate = 0;
    lower = 0;
    auto f2 = [&](double r) { return zstatRateRatio(t1, y1, t2, y2, r) + b; };
    upper = brent(f2, 0.001, 1000.0, 1e-6);
  } else if (all_y2_zero) {
    estimate = POS_INF;
    upper = POS_INF;
    auto f1 = [&](double r) { return zstatRateRatio(t1, y1, t2, y2, r) - b; };
    lower = brent(f1, 0.001, 1000.0, 1e-6);
  } else {
    estimate = r1 / r2;
    auto f1 = [&](double r) { return zstatRateRatio(t1, y1, t2, y2, r) - b; };
    auto f2 = [&](double r) { return zstatRateRatio(t1, y1, t2, y2, r) + b; };
    lower = brent(f1, 0.001, estimate, 1e-6);
    upper = brent(f2, estimate, 1000.0, 1e-6);
  }

  DataFrameCpp df;
  df.push_back(std::string("rate ratio"), "scale");
  df.push_back(estimate, "estimate");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");
  df.push_back(cilevel, "cilevel");

  return df;
}


//' @title Miettinen-Nurminen Score Confidence Interval for
//' Two-Sample Rate Ratio
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample rate ratio possibly with
//' stratification.
//'
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel weights are used for stratified samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input exposure
//'   and number of events for each treatment group.
//'   It has the following variables:
//'
//'     - \code{t1}: The exposure for the active treatment group.
//'
//'     - \code{y1}: The number of events for the active treatment group.
//'
//'     - \code{t2}: The exposure for the control group.
//'
//'     - \code{y2}: The number of events for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for rate ratio. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnRateRatioCI(t1 = c(10,10), y1 = c(4,3), t2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::DataFrame mnRateRatioCI(const std::vector<double>& t1,
                              const std::vector<double>& y1,
                              const std::vector<double>& t2,
                              const std::vector<double>& y2,
                              const double cilevel = 0.95) {

  auto df = mnRateRatioCIcpp(t1, y1, t2, y2, cilevel);
  return Rcpp::wrap(df);
}
