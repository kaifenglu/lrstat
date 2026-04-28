#include "generic_design.h"
#include "utilities.h"
#include "dataframe_list.h"

#include <algorithm>     // any_of, fill
#include <cctype>        // tolower
#include <cmath>         // isnan
#include <cstring>       // memcpy
#include <numeric>       // accumulate
#include <stdexcept>     // invalid_argument
#include <string>        // string
#include <vector>        // vector
#include <utility>       // pair, make_pair

#include <Rcpp.h>

using std::size_t;


// Helper to compute the confidence interval at the end of a group sequential trial
DataFrameCpp getCIcpp(const size_t L,
                      const double zL,
                      const double IMax,
                      const std::vector<double>& informationRates,
                      const std::vector<unsigned char>& efficacyStopping,
                      const std::vector<double>& criticalValues,
                      const double alpha,
                      const std::string& typeAlphaSpending,
                      const double parameterAlphaSpending,
                      const std::vector<double>& spendingTime) {

  // Basic argument checks
  if (L <= 0) throw std::invalid_argument("L must be a positive integer");

  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (!none_na(informationRates))
    throw std::invalid_argument("informationRates must be provided");

  size_t kMax = informationRates.size();
  if (kMax < L)
    throw std::invalid_argument("Insufficient length for informationRates");
  if (informationRates[0] <= 0.0)
    throw std::invalid_argument("informationRates must be positive");
  if (any_nonincreasing(informationRates))
    throw std::invalid_argument("informationRates must be increasing");
  if (informationRates.back() > 1.0)
    throw std::invalid_argument("informationRates must not exceed 1");

  // efficacyStopping: if provided, validate, otherwise default to all ones
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() < kMax)
      throw std::invalid_argument("Insufficient length for efficacyStopping");
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping;
  } else {
    effStopping.assign(kMax, 1);
  }

  // spendingTime: if provided validate, otherwise use informationRates
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() < kMax)
      throw std::invalid_argument("Insufficient length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime.back() > 1.0)
      throw std::invalid_argument("spendingTime must not exceed 1");
    spendTime = spendingTime;
  } else {
    spendTime = informationRates;
  }

  // alpha checks
  if (std::isnan(alpha)) throw std::invalid_argument("alpha must be provided");
  if (alpha < 0.00001 || alpha >= 0.5)
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");

  // typeAlphaSpending to lower-case
  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (!none_na(criticalValues) && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument("parameterAlphaSpending must be positive for sfKD");
  }

  if (asf == "of" || asf == "p" || asf == "wt" || asf == "none") {
    if (informationRates.back() != 1.0) {
      throw std::invalid_argument(
          "informationRates must end with 1 for OF, P, WT, or NONE");
    }
    if (spendTime.back() != 1.0) {
      throw std::invalid_argument(
          "spendingTime must end with 1 for OF, P, WT, or NONE");
    }
  }

  // critical values: if not provided, compute using getBoundcpp
  std::vector<double> b;
  if (none_na(criticalValues)) {
    if (criticalValues.size() < L)
      throw std::invalid_argument("Insufficient length for criticalValues");
    b = criticalValues;
  } else {
    b = getBoundcpp(L, informationRates, alpha, asf, parameterAlphaSpending,
                    std::vector<double>{}, spendTime, effStopping);
  }

  // Build full information vector I = IMax * informationRates
  std::vector<double> I(L);
  for (size_t i = 0; i < L; ++i) I[i] = IMax * informationRates[i];

  // p-value at theta = 0
  double pvalue = f_pvalue(0.0, L, zL, b, I);

  double cilevel = 1.0 - 2.0 * alpha;

  // initial bracketing interval for theta: (zL +/- 6) / sqrt(I_L)
  double sqrtIL = std::sqrt(I[L - 1]);
  double left = (zL - 6.0) / sqrtIL;
  double right = (zL + 6.0) / sqrtIL;
  double tol = 1.0e-6;

  // median-unbiased estimate thetahat: solve f_pvalue(theta) - 0.5 = 0
  auto f_med = [&](double theta)->double {
    return f_pvalue(theta, L, zL, b, I) - 0.5;
  };
  double thetahat = brent(f_med, left, right, tol);

  // lower bound: solve f_pvalue(theta) - (1-cilevel)/2 = 0, in [left, thetahat]
  double target_lower = (1.0 - cilevel) / 2.0;
  auto f_lower = [&](double theta)->double {
    return f_pvalue(theta, L, zL, b, I) - target_lower;
  };
  double lower = brent(f_lower, left, thetahat, tol);

  // upper bound: solve f_pvalue(theta) - (1+cilevel)/2 = 0, in [thetahat, right]
  double target_upper = (1.0 + cilevel) / 2.0;
  auto f_upper = [&](double theta)->double {
    return f_pvalue(theta, L, zL, b, I) - target_upper;
  };
  double upper = brent(f_upper, thetahat, right, tol);

  // Build DataFrameCpp result
  DataFrameCpp df;
  df.push_back(pvalue, "pvalue");
  df.push_back(thetahat, "thetahat");
  df.push_back(cilevel, "cilevel");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");

  return df;
}



//' @title Confidence Interval After Trial Termination
//' @description Obtains the p-value, median unbiased point estimate, and
//' confidence interval after the end of a group sequential trial.
//'
//' @param L The termination look.
//' @param zL The z-test statistic at the termination look.
//' @param IMax The maximum information of the trial.
//' @param informationRates The information rates up to look \code{L}.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage up to look \code{L}.
//'   Defaults to true if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping up to look \code{L}.
//' @inheritParams param_alpha
//' @param typeAlphaSpending The type of alpha spending for the trial.
//'   One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @inheritParams param_parameterAlphaSpending
//' @param spendingTime The error spending time up to look \code{L}.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//'
//' @return A data frame with the following components:
//'
//' * \code{pvalue}: p-value for rejecting the null hypothesis.
//'
//' * \code{thetahat}: Median unbiased point estimate of the parameter.
//'
//' * \code{cilevel}: Confidence interval level.
//'
//' * \code{lower}: Lower bound of confidence interval.
//'
//' * \code{upper}: Upper bound of confidence interval.
//'
//' @details
//' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
//' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
//' and \code{spendingTime} must be of full length \code{kMax}, and
//' \code{informationRates} and \code{spendingTime} must end with 1.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Anastasios A. Tsiatis, Gary L. Rosner and Cyrus R. Mehta.
//' Exact confidence intervals following a group sequential test.
//' Biometrics 1984;40:797-803.
//'
//' @examples
//'
//' # group sequential design with 90% power to detect delta = 6
//' delta = 6
//' sigma = 17
//' n = 282
//' (des1 = getDesign(IMax = n/(4*sigma^2), theta = delta, kMax = 3,
//'                   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'                   parameterAlphaSpending = -4))
//'
//' # crossed the boundary at the second look
//' L = 2
//' n1 = n*2/3
//' delta1 = 7
//' sigma1 = 20
//' zL = delta1/sqrt(4/n1*sigma1^2)
//'
//' # confidence interval
//' getCI(L = L, zL = zL, IMax = n/(4*sigma1^2),
//'       informationRates = c(1/3, 2/3), alpha = 0.05,
//'       typeAlphaSpending = "sfHSD", parameterAlphaSpending = -4)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame getCI(
    const int L = NA_INTEGER,
    const double zL = NA_REAL,
    const double IMax = NA_REAL,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {
  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto result = getCIcpp(static_cast<size_t>(L), zL, IMax, infoRates,
                         effStopping, critValues, alpha, typeAlphaSpending,
                         parameterAlphaSpending, spendTime);
  return Rcpp::wrap(result);
}


// Helper to compute repeated confidence interval
DataFrameCpp getRCIcpp(
    const size_t L,
    const double zL,
    const double IMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& spendingTime) {

  // Basic argument checks
  if (L <= 0) throw std::invalid_argument("L must be a positive integer");

  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (!none_na(informationRates))
    throw std::invalid_argument("informationRates must be provided");

  size_t kMax = informationRates.size();
  if (kMax < L)
    throw std::invalid_argument("Insufficient length for informationRates");
  if (informationRates[0] <= 0.0)
    throw std::invalid_argument("informationRates must be positive");
  if (any_nonincreasing(informationRates))
    throw std::invalid_argument("informationRates must be increasing");
  if (informationRates.back() > 1.0)
    throw std::invalid_argument("informationRates must not exceed 1");

  // efficacyStopping: if provided, validate, otherwise default to all ones
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() < kMax)
      throw std::invalid_argument("Insufficient length for efficacyStopping");
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping;
  } else {
    effStopping.assign(kMax, 1);
  }

  // spendingTime: if provided validate, otherwise use informationRates
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() < kMax)
      throw std::invalid_argument("Insufficient length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime.back() > 1.0)
      throw std::invalid_argument("spendingTime must not exceed 1");
    spendTime = spendingTime;
  } else {
    spendTime = informationRates;
  }

  // alpha checks
  if (std::isnan(alpha)) throw std::invalid_argument("alpha must be provided");
  if (alpha < 0.00001 || alpha >= 0.5)
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");

  // typeAlphaSpending to lower-case
  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (!none_na(criticalValues) && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument("parameterAlphaSpending must be positive for sfKD");
  }

  if (asf == "of" || asf == "p" || asf == "wt" || asf == "none") {
    if (informationRates.back() != 1.0) {
      throw std::invalid_argument(
          "informationRates must end with 1 for OF, P, WT, or NONE");
    }
    if (spendTime.back() != 1.0) {
      throw std::invalid_argument(
          "spendingTime must end with 1 for OF, P, WT, or NONE");
    }
  }

  // critical values: if not provided, compute using getBoundcpp with caching
  BoundCacheAlpha cache(L, informationRates, asf, parameterAlphaSpending,
                        std::vector<double>{}, spendTime, effStopping, 64, 12);
  std::vector<double> b;
  if (none_na(criticalValues)) {
    if (criticalValues.size() < L)
      throw std::invalid_argument("Insufficient length for criticalValues");
    b = criticalValues;
  } else {
    b = cache.get(alpha);
  }

  // Build full information vector I = IMax * informationRates
  std::vector<double> I(L);
  for (size_t i = 0; i < L; ++i) I[i] = IMax * informationRates[i];

  // repeated confidence interval
  double sqrtIL = std::sqrt(I[L-1]);
  double thetahat = zL / sqrtIL;
  double lower = (zL - b[L-1]) / sqrtIL;
  double upper = (zL + b[L-1]) / sqrtIL;

  // repeated p-value: alpha for which the lower bound of theta is zero:
  // Solve f(aval) = zL - u[ L-1 ] = 0 for aval in (1e-6, 0.999999)
  auto f = [&](double aval)->double {
    std::vector<double> u_local = cache.get(aval);
    return zL - u_local[L-1];
  };

  double pvalue;
  double lo_p = 0.000001;
  double hi_p = 0.999999;
  double f_lo = f(lo_p);
  if (f_lo > 0.0) {
    pvalue = lo_p;
  } else {
    double f_hi = f(hi_p);
    if (f_hi < 0.0) {
      pvalue = hi_p;
    } else {
      pvalue = brent(f, lo_p, hi_p, 1e-6);
    }
  }

  DataFrameCpp df;
  df.push_back(pvalue, "pvalue");
  df.push_back(thetahat, "thetahat");
  df.push_back(1.0 - 2.0 * alpha, "cilevel");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");

  return df;
}


//' @title Repeated Confidence Interval for Group Sequential Design
//' @description Obtains the repeated confidence interval
//' for a group sequential trial.
//'
//' @param L The look of interest.
//' @param zL The z-test statistic at the look.
//' @param IMax The maximum information of the trial.
//' @param informationRates The information rates up to look \code{L}.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage up to look \code{L}. Defaults to true
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping up to look \code{L}.
//' @inheritParams param_alpha
//' @param typeAlphaSpending The type of alpha spending for the trial.
//'   One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @param spendingTime The error spending time up to look \code{L}.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//'
//' @return A data frame with the following components:
//'
//' * \code{pvalue}: Repeated p-value for rejecting the null hypothesis.
//'
//' * \code{thetahat}: Point estimate of the parameter.
//'
//' * \code{cilevel}: Confidence interval level.
//'
//' * \code{lower}: Lower bound of repeated confidence interval.
//'
//' * \code{upper}: Upper bound of repeated confidence interval.
//'
//' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
//' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
//' and \code{spendingTime} must be of full length \code{kMax}, and
//' \code{informationRates} and \code{spendingTime} must end with 1.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Christopher Jennison and Bruce W. Turnbull.
//' Interim analyses: the repeated confidence interval approach
//' (with discussion).
//' J R Stat Soc Series B. 1989;51:305-361.
//'
//' @examples
//'
//' # group sequential design with 90% power to detect delta = 6
//' delta = 6
//' sigma = 17
//' n = 282
//' (des1 = getDesign(IMax = n/(4*sigma^2), theta = delta, kMax = 3,
//'                   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'                   parameterAlphaSpending = -4))
//'
//' # results at the second look
//' L = 2
//' n1 = n*2/3
//' delta1 = 7
//' sigma1 = 20
//' zL = delta1/sqrt(4/n1*sigma1^2)
//'
//' # repeated confidence interval
//' getRCI(L = L, zL = zL, IMax = n/(4*sigma1^2),
//'        informationRates = c(1/3, 2/3), alpha = 0.05,
//'        typeAlphaSpending = "sfHSD", parameterAlphaSpending = -4)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame getRCI(
    const int L = NA_INTEGER,
    const double zL = NA_REAL,
    const double IMax = NA_REAL,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {
  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto result = getRCIcpp(static_cast<size_t>(L), zL, IMax, infoRates,
                          effStopping, critValues, alpha, typeAlphaSpending,
                          parameterAlphaSpending, spendTime);
  return Rcpp::wrap(result);
}


// Compute the backward image (J, zJ)
std::pair<size_t, double> f_bwimage(const double theta,
                                    const size_t kMax,
                                    const size_t L,
                                    const double zL,
                                    const std::vector<double>& b,
                                    const std::vector<double>& I,
                                    const size_t L2,
                                    const double zL2,
                                    const std::vector<double>& b2,
                                    const std::vector<double>& I2) {

  // compute astar for the adapted secondary trial
  double astar = f_pvalue(theta, L2, zL2, b2, I2);

  // prepare b1, a1, mu, I1 for the original secondary trial
  size_t k1 = kMax - L;
  std::vector<double> I1(k1);
  for (size_t l = 0; l < k1; ++l) {
    I1[l] = I[l + L] - I[L - 1];
  }

  std::vector<double> b1(k1);
  std::vector<double> a1(k1, -6.0);
  std::vector<double> mu(k1, theta);
  for (size_t l = 0; l < k1; ++l) {
    double r1 = I[L - 1] / I[l + L];
    b1[l] = (b[l + L] - zL * std::sqrt(r1)) / std::sqrt(1.0 - r1);
  }

  // compute exit probabilities for b1 / a1
  ListCpp probs = exitprobcpp(b1, a1, mu, I1);
  auto pu = probs.get<std::vector<double>>("exitProbUpper");

  // find interval containing astar
  std::vector<double> cpu(k1);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  size_t j = std::min(findInterval1(astar, cpu) + 1, k1);
  size_t J = L + j; // combined stage index in primary trial numbering

  // find zJ
  double zJ;
  if (j == 1) {
    double r1 = I[L - 1] / I[L];
    zJ = boost_qnorm(1.0 - astar) * std::sqrt(1.0 - r1) + zL * std::sqrt(r1);
  } else {
    // root find for z in stagewise exit probability difference
    double r1 = I[L - 1] / I[L + j - 1];
    auto f = [&](double z)->double {
      double zj = (z - zL * std::sqrt(r1)) / std::sqrt(1.0 - r1);
      return f_pvalue(theta, j, zj, b1, I1) - astar;
    };

    if (j < k1) {
      zJ = brent(f, b[L + j - 1], 6.0, 1e-6);
    } else {
      double r1 = I[L - 1] / I[L + j - 1];
      double lo = -6.0 * std::sqrt(1.0 - r1) + zL * std::sqrt(r1);
      zJ = brent(f, lo, 6.0, 1e-6);
    }
  }

  return std::make_pair(J, zJ);
}


// compute backward p-value for adapted trial
double f_bwpvalue(const double theta,
                  const size_t kMax,
                  const size_t L,
                  const double zL,
                  const std::vector<double>& b,
                  const std::vector<double>& I,
                  const size_t L2,
                  const double zL2,
                  const std::vector<double>& b2,
                  const std::vector<double>& I2) {
  auto bw = f_bwimage(theta, kMax, L, zL, b, I, L2, zL2, b2, I2);
  return f_pvalue(theta, bw.first, bw.second, b, I);
}


// Helper to compute confidence interval after the end of an adaptive trial
DataFrameCpp getADCIcpp(
    const size_t L,
    const double zL,
    const double IMax,
    const size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
    const size_t Lc,
    const double zLc,
    const double INew,
    const std::vector<double>& informationRatesNew,
    const std::vector<unsigned char>& efficacyStoppingNew,
    const std::string& typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double>& spendingTimeNew) {

  // Input validation and defaults
  if (L <= 0) throw std::invalid_argument("L must be provided and positive");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (informationRates.size() != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[kMax - 1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy if provided
  } else {
    infoRates.resize(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      infoRates[i] = static_cast<double>(i + 1) / static_cast<double>(kMax);
    }
  }

  // efficacyStopping: default to all ones if not provided
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[kMax - 1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy if provided
  } else {
    effStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  if (!missingCriticalValues && criticalValues.size() != kMax) {
    throw std::invalid_argument("Invalid length for criticalValues");
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 0.5)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax)
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


  // Now handle new trial inputs
  if (Lc <= L) throw std::invalid_argument("Lc must be greater than L");
  if (std::isnan(zLc)) throw std::invalid_argument("zLc must be provided");
  if (std::isnan(INew)) throw std::invalid_argument("INew must be provided");
  if (INew <= 0.0) throw std::invalid_argument("INew must be positive");

  size_t L2 = Lc - L;

  std::vector<unsigned char> effStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  size_t kNew = L2;
  if (MullerSchafer) {
    if (none_na(informationRatesNew)) {
      kNew = informationRatesNew.size();
      if (kNew < L2)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew.back() > 1.0)
        throw std::invalid_argument("informationRatesNew must not exceed 1");
    } else {
      throw std::invalid_argument(
          "informationRatesNew must be provided for MullerSchafer");
    }

    if (none_na(efficacyStoppingNew)) {
      if (efficacyStoppingNew.size() < kNew)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew.back() != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew;
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    if (!(asfNew == "of" || asfNew == "p" || asfNew == "wt" ||
        asfNew == "sfof" || asfNew == "sfp" ||
        asfNew == "sfkd" || asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }

    if ((asfNew == "wt" || asfNew == "sfkd" || asfNew == "sfhsd") &&
        std::isnan(parameterAlphaSpendingNew))
      throw std::invalid_argument("Missing value for parameterAlphaSpendingNew");

    if (asfNew == "sfkd" && parameterAlphaSpendingNew <= 0.0)
      throw std::invalid_argument(
          "parameterAlphaSpendingNew must be positive for sfKD");

    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() < kNew)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew.back() > 1.0)
        throw std::invalid_argument("spendingTimeNew must not exceed 1");
      spendTimeNew = spendingTimeNew;
    } else {
      spendTimeNew = informationRatesNew;
    }

    if (asfNew == "of" || asfNew == "p" || asfNew == "wt" || asfNew == "none") {
      if (informationRatesNew.back() != 1.0) {
        throw std::invalid_argument(
            "informationRatesNew must end with 1 for OF, P, WT, or NONE");
      }
      if (spendTimeNew.back() != 1.0) {
        throw std::invalid_argument(
            "spendingTimeNew must end with 1 for OF, P, WT, or NONE");
      }
    }
  }


  // obtain critical values for the primary trial
  std::vector<double> l(kMax, -6.0), zero(kMax, 0.0);
  std::vector<double> b = criticalValues;
  double alpha1 = alpha;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (kMax > 1 && criticalValues.size() == kMax) {
      bool hasNaN = false;
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[kMax-1])) haybittle = true;
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(kMax);
      for (size_t i = 0; i < kMax - 1; ++i) {
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

      b[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      b = getBoundcpp(kMax, infoRates, alpha, asf, parameterAlphaSpending,
                      std::vector<double>{}, spendTime, effStopping);
    }
  } else {
    for (size_t i = 0; i < kMax; ++i) {
      if (!effStopping[i]) b[i] = 6.0;
    }
    ListCpp probs = exitprobcpp(b, l, zero, infoRates);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    alpha1 = std::accumulate(v.begin(), v.end(), 0.0);
  }


  // Primary information vector
  std::vector<double> I(kMax);
  for (size_t i = 0; i < kMax; ++i) I[i] = IMax * informationRates[i];

  // compute b2 and I2 for secondary trial depending on MullerSchafer
  std::vector<double> b2(L2), I2(L2);
  if (!MullerSchafer) {
    for (size_t l = 0; l < L2; ++l) {
      double t1 = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      double r1 = infoRates[L - 1] / infoRates[l + L];
      b2[l] = (b[l + L] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
      if (!effStopping[l + L]) b2[l] = 6.0;
      I2[l] = INew * t1;
    }
  } else { // conditional type I error
    size_t k1 = kMax - L;
    std::vector<double> t1(k1), b1(k1), a1(k1, -6.0);
    for (size_t l = 0; l < k1; ++l) {
      t1[l] = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      double r1 = infoRates[L - 1] / infoRates[l + L];
      b1[l] = (b[l + L] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
      if (!effStopping[l + L]) b1[l] = 6.0;
    }
    ListCpp probs = exitprobcpp(b1, a1, zero, t1);
    auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
    double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);

    b2 = getBoundcpp(L2, informationRatesNew, alphaNew, asfNew,
                     parameterAlphaSpendingNew, std::vector<double>{},
                     spendTimeNew, effStoppingNew);

    for (size_t l = 0; l < L2; ++l) I2[l] = INew * informationRatesNew[l];
  }

  // confidence level
  double cilevel = 1.0 - 2.0 * alpha1;

  // compute pvalue under theta=0 using f_bwpvalue
  double zL2 = (zLc * std::sqrt(I[L - 1] + I2[L2 - 1]) -
                zL * std::sqrt(I[L - 1])) / std::sqrt(I2[L2 - 1]);
  double pvalue = f_bwpvalue(0.0, kMax, L, zL, b, I, L2, zL2, b2, I2);

  // interval brackets and root-finding to obtain thetahat, lower, upper
  double sqrtIL = std::sqrt(I[L - 1]);
  double left = (zL - b[L - 1]) / sqrtIL;
  double right = (zL + b[L - 1]) / sqrtIL;
  double tol = 1.0e-6;

  auto f_med = [&](double theta)->double {
    return f_bwpvalue(theta, kMax, L, zL, b, I, L2, zL2, b2, I2) - 0.5;
  };
  double thetahat = brent(f_med, left, right, tol);

  double target_lower = (1.0 - cilevel) / 2.0;
  auto f_low = [&](double theta)->double {
    return f_bwpvalue(theta, kMax, L, zL, b, I, L2, zL2, b2, I2) - target_lower;
  };
  double lower = brent(f_low, left, thetahat, tol);

  double target_upper = (1.0 + cilevel) / 2.0;
  auto f_high = [&](double theta)->double {
    return f_bwpvalue(theta, kMax, L, zL, b, I, L2, zL2, b2, I2) - target_upper;
  };
  double upper = brent(f_high, thetahat, right, tol);

  DataFrameCpp df;
  df.push_back(pvalue, "pvalue");
  df.push_back(thetahat, "thetahat");
  df.push_back(cilevel, "cilevel");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");
  return df;
}


//' @title Confidence Interval After Adaptation
//' @description Obtains the p-value, median unbiased point estimate, and
//' confidence interval after the end of an adaptive trial.
//'
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of
//'   the primary trial.
//' @param IMax The maximum information of the primary trial.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to true
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping for the primary trial.
//' @param alpha The significance level of the primary trial.
//'   Defaults to 0.025.
//' @param typeAlphaSpending The type of alpha spending for the primary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param Lc The termination look of the integrated trial.
//' @param zLc The z-test statistic at the termination look of the
//'   integrated trial.
//' @param INew The maximum information of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial
//'   up to look \code{L2}.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial up to look \code{L2}.
//'   Defaults to true if left unspecified.
//' @param typeAlphaSpendingNew The type of alpha spending for the secondary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param spendingTimeNew The error spending time of the secondary trial
//'   up to look \code{L2}. Defaults to missing, in which case, it is
//'   the same as \code{informationRatesNew}.
//'
//' @return A data frame with the following variables:
//'
//' * \code{pvalue}: p-value for rejecting the null hypothesis.
//'
//' * \code{thetahat}: Median unbiased point estimate of the parameter.
//'
//' * \code{cilevel}: Confidence interval level.
//'
//' * \code{lower}: Lower bound of confidence interval.
//'
//' * \code{upper}: Upper bound of confidence interval.
//'
//' @details
//' If \code{typeAlphaSpendingNew} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
//' \code{"none"}, then \code{informationRatesNew}, \code{efficacyStoppingNew},
//' and \code{spendingTimeNew} must be of full length \code{kNew}, and
//' \code{informationRatesNew} and \code{spendingTimeNew} must end with 1.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Lingyun Liu and Cyrus Mehta.
//' Exact inference for adaptive group sequential designs.
//' Stat Med. 2013;32(23):3991-4005.
//'
//' @seealso \code{\link{adaptDesign}}
//'
//' @examples
//' # two-arm randomized clinical trial with a normally distributed endpoint
//' # 90% power to detect mean difference of 15 with a standard deviation of 50
//' # Design the Stage I Trial with 3 looks and Lan-DeMets O'Brien-Fleming type
//' # spending function
//' delta <- 15
//' sigma <- 50
//'
//' (des1 <- getDesignMeanDiff(
//'   beta = 0.1, meanDiff = delta, stDev = sigma,
//'   kMax = 3, alpha = 0.025, typeAlphaSpending = "sfOF"
//' ))
//'
//' s1 <- des1$byStageResults$informationRates
//' b1 <- des1$byStageResults$efficacyBounds
//' n <- des1$overallResults$numberOfSubjects
//'
//' # Monitoring the Stage I Trial
//' L <- 1
//' nL <- des1$byStageResults$numberOfSubjects[L]
//' deltahat <- 8
//' sigmahat <- 55
//' sedeltahat <- sigmahat * sqrt( 4 / nL)
//' zL <- deltahat / sedeltahat
//'
//' # Making an Adaptive Change: Stage I to Stage II
//' # revised clinically meaningful difference downward to 10 power the study
//' # retain the standard deviation at the design stage
//' # Muller & Schafer (2001) method to design the secondary trial
//' # with 2 looks and Lan-DeMets Pocock type spending function
//' # re-estimate sample size to reach 90% conditional power
//' deltaNew <- 10
//'
//' (des2 <- adaptDesign(
//'   betaNew = 0.1, L = L, zL = zL, theta = deltaNew,
//'   IMax = n / (4 * sigma^2), kMax = 3, informationRates = s1,
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   MullerSchafer = TRUE, kNew = 2, typeAlphaSpendingNew = "sfP"
//' ))
//'
//' INew <- des2$secondaryTrial$overallResults$information
//' (nNew <- ceiling(INew * 4 * sigma^2))
//' (nTotal <- nL + nNew)
//'
//' # Monitoring the Integrated Trial
//' s2 <- des2$secondaryTrial$byStageResults$informationRates
//'
//' Lc <- 2
//' deltahatc <- 9.5
//' sigmahatc <- 52.759
//' L2 <- Lc - L
//' nL2 <-  nNew * s2[L2]
//' nc <- nL + nL2
//' sedeltahatc <- sigmahatc * sqrt(4 / nc)
//' zLc <- deltahatc / sedeltahatc
//' zL2 <- (zLc * sqrt(nc) - zL * sqrt(nL)) / sqrt(nL2)
//'
//' getADCI(
//'   L = L, zL = zL, IMax = n / (4 * sigmahatc^2), kMax = 3,
//'   informationRates = s1, alpha = 0.025, typeAlphaSpending = "sfOF",
//'   MullerSchafer = TRUE, Lc = Lc, zLc = zLc,
//'   INew = nNew / (4 * sigmahatc^2), informationRatesNew = s2,
//'   typeAlphaSpendingNew = "sfP")
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame getADCI(
    const int L = NA_INTEGER,
    const double zL = NA_REAL,
    const double IMax = NA_REAL,
    const int kMax = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.25,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int Lc = NA_INTEGER,
    const double zLc = NA_REAL,
    const double INew = NA_REAL,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL) {
  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);
  auto result = getADCIcpp(
    static_cast<size_t>(L), zL, IMax, static_cast<size_t>(kMax), infoRates,
    effStopping, critValues, alpha, typeAlphaSpending, parameterAlphaSpending,
    spendTime, MullerSchafer, static_cast<size_t>(Lc), zLc, INew,
    infoRatesNew, effStoppingNew, typeAlphaSpendingNew,
    parameterAlphaSpendingNew, spendTimeNew);
  return Rcpp::wrap(result);
}


// Helper to calculate repeated confidence interval after adaptation
DataFrameCpp getADRCIcpp(
    const size_t L,
    const double zL,
    const double IMax,
    const size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
    const size_t Lc,
    const double zLc,
    const double INew,
    const std::vector<double>& informationRatesNew,
    const std::vector<unsigned char>& efficacyStoppingNew,
    const std::string& typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double>& spendingTimeNew) {

  // Input validation and defaults
  if (L <= 0) throw std::invalid_argument("L must be provided and positive");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (informationRates.size() != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[kMax - 1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy if provided
  } else {
    infoRates.resize(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      infoRates[i] = static_cast<double>(i + 1) / static_cast<double>(kMax);
    }
  }

  // efficacyStopping: default to all ones if not provided
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[kMax - 1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy if provided
  } else {
    effStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  if (!missingCriticalValues && criticalValues.size() != kMax) {
    throw std::invalid_argument("Invalid length for criticalValues");
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 0.5)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax)
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


  // Now handle new trial inputs
  if (Lc <= L) throw std::invalid_argument("Lc must be greater than L");
  if (std::isnan(zLc)) throw std::invalid_argument("zLc must be provided");
  if (std::isnan(INew)) throw std::invalid_argument("INew must be provided");
  if (INew <= 0.0) throw std::invalid_argument("INew must be positive");

  size_t L2 = Lc - L;

  std::vector<unsigned char> effStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  size_t kNew = L2;
  if (MullerSchafer) {
    if (none_na(informationRatesNew)) {
      kNew = informationRatesNew.size();
      if (kNew < L2)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew.back() > 1.0)
        throw std::invalid_argument("informationRatesNew must not exceed 1");
    } else {
      throw std::invalid_argument(
          "informationRatesNew must be provided for MullerSchafer");
    }

    if (none_na(efficacyStoppingNew)) {
      if (efficacyStoppingNew.size() < kNew)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew.back() != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew;
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    if (!(asfNew == "of" || asfNew == "p" || asfNew == "wt" ||
        asfNew == "sfof" || asfNew == "sfp" ||
        asfNew == "sfkd" || asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }

    if ((asfNew == "wt" || asfNew == "sfkd" || asfNew == "sfhsd") &&
        std::isnan(parameterAlphaSpendingNew))
      throw std::invalid_argument("Missing value for parameterAlphaSpendingNew");

    if (asfNew == "sfkd" && parameterAlphaSpendingNew <= 0.0)
      throw std::invalid_argument(
          "parameterAlphaSpendingNew must be positive for sfKD");

    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() < kNew)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew.back() > 1.0)
        throw std::invalid_argument("spendingTimeNew must not exceed 1");
      spendTimeNew = spendingTimeNew;
    } else {
      spendTimeNew = informationRatesNew;
    }

    if (asfNew == "of" || asfNew == "p" || asfNew == "wt" || asfNew == "none") {
      if (informationRatesNew.back() != 1.0) {
        throw std::invalid_argument(
            "informationRatesNew must end with 1 for OF, P, WT, or NONE");
      }
      if (spendTimeNew.back() != 1.0) {
        throw std::invalid_argument(
            "spendingTimeNew must end with 1 for OF, P, WT, or NONE");
      }
    }
  }


  // obtain critical values for the primary trial
  std::vector<double> l(kMax, -6.0), zero(kMax, 0.0);
  std::vector<double> b = criticalValues;
  double alpha1 = alpha;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (kMax > 1 && criticalValues.size() == kMax) {
      bool hasNaN = false;
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[kMax-1])) haybittle = true;
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(kMax);
      for (size_t i = 0; i < kMax - 1; ++i) {
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

      b[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      b = getBoundcpp(kMax, infoRates, alpha, asf, parameterAlphaSpending,
                      std::vector<double>{}, spendTime, effStopping);
    }
  } else {
    for (size_t i = 0; i < kMax; ++i) {
      if (!effStopping[i]) b[i] = 6.0;
    }
    ListCpp probs = exitprobcpp(b, l, zero, infoRates);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    alpha1 = std::accumulate(v.begin(), v.end(), 0.0);
  }


  // Primary information vector
  std::vector<double> I(kMax);
  for (size_t i = 0; i < kMax; ++i) I[i] = IMax * informationRates[i];

  // compute b2 and I2 for secondary trial depending on MullerSchafer
  std::vector<double> b2(L2), I2(L2);
  if (!MullerSchafer) {
    for (size_t l = 0; l < L2; ++l) {
      double t1 = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      double r1 = infoRates[L - 1] / infoRates[l + L];
      b2[l] = (b[l + L] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
      if (!effStopping[l + L]) b2[l] = 6.0;
      I2[l] = INew * t1;
    }
  } else { // conditional type I error
    size_t k1 = kMax - L;
    std::vector<double> t1(k1), b1(k1), a1(k1, -6.0);
    for (size_t l = 0; l < k1; ++l) {
      t1[l] = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      double r1 = infoRates[L - 1] / infoRates[l + L];
      b1[l] = (b[l + L] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
      if (!effStopping[l + L]) b1[l] = 6.0;
    }
    ListCpp probs = exitprobcpp(b1, a1, zero, t1);
    auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
    double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);

    b2 = getBoundcpp(L2, informationRatesNew, alphaNew, asfNew,
                     parameterAlphaSpendingNew, std::vector<double>{},
                     spendTimeNew, effStoppingNew);

    for (size_t l = 0; l < L2; ++l) I2[l] = INew * informationRatesNew[l];
  }

  // confidence level
  double cilevel = 1.0 - 2.0 * alpha1;

  double zL2 = (zLc * std::sqrt(I[L - 1] + I2[L2 - 1]) -
                zL * std::sqrt(I[L - 1])) / std::sqrt(I2[L2 - 1]);

  // compute repeated pvalue, thetahat, lower, upper
  double lower = NAN, upper = NAN, thetahat = NAN, pvalue = NAN;

  if (!MullerSchafer) {
    // simple combination formula
    double I1 = IMax * infoRates[L - 1];
    double I2 = INew * (infoRates[L + L2 - 1] - infoRates[L - 1]) /
      (1.0 - infoRates[L - 1]);
    double r1 = infoRates[L - 1] / infoRates[L + L2 - 1];
    double w1 = std::sqrt(r1);
    double w2 = std::sqrt(1.0 - r1);
    double c1 = w1 * zL + w2 * zL2;
    double c2 = w1 * std::sqrt(I1) + w2 * std::sqrt(I2);

    thetahat = c1 / c2;
    lower = (c1 - b[L + L2 - 1]) / c2;
    upper = (c1 + b[L + L2 - 1]) / c2;

    size_t J = L + L2;

    // create cache_J for J
    std::vector<double> t_prefix(J);
    for (size_t i = 0; i < J; ++i) t_prefix[i] = infoRates[i];
    BoundCacheAlpha cache_J(J, t_prefix, asf, parameterAlphaSpending,
                            std::vector<double>{}, spendTime, effStopping, 64, 12);

    // repeated p-value: solve f(aval) = c1 - u[J-1] = 0 for aval in (1e-6, 0.999999)
    auto f_alpha = [&](double aval)->double {
      std::vector<double> u_local = cache_J.get(aval);
      return c1 - u_local[J - 1];
    };

    double p_lo = 0.000001;
    double p_hi = 0.999999;
    double f_lo = f_alpha(p_lo);
    if (f_lo > 0.0) {
      pvalue = p_lo;
    } else {
      double f_hi = f_alpha(p_hi);
      if (f_hi < 0.0) {
        pvalue = p_hi;
      } else {
        pvalue = brent(f_alpha, p_lo, p_hi, 1.0e-6);
      }
    }
  } else {
    // -------------------- Setup small alpha-only caches --------------------
    // cache for kMax
    BoundCacheAlpha cache_kMax(kMax, infoRates, asf, parameterAlphaSpending,
                               std::vector<double>{}, spendTime, effStopping,
                               64, 12);
    // cache for L2 (used inside g)
    BoundCacheAlpha cache_L2(L2, informationRatesNew, asfNew,
                             parameterAlphaSpendingNew, std::vector<double>{},
                             spendTimeNew, effStoppingNew, 64, 12);

    // Muller-Schafer branch: more complex iterative root-finding
    double I1 = IMax * infoRates[L - 1];
    double I2 = INew * informationRatesNew[L2 - 1];
    double sqrtI1 = std::sqrt(I1);
    double sqrtI2 = std::sqrt(I2);

    size_t k1 = kMax - L;

    std::vector<double> t1(k1), w1(k1), w2(k1);
    for (size_t l = 0; l < k1; ++l) {
      t1[l] = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      double r1 = infoRates[L - 1] / infoRates[l + L];
      w1[l] = std::sqrt(r1);
      w2[l] = std::sqrt(1.0 - r1);
    }

    // interval for root-finding of theta
    double left = (zL - b[L - 1]) / sqrtI1;
    double right = (zL + b[L - 1]) / sqrtI1;
    double tol = 1.0e-6;

    std::vector<double> b1(k1);
    std::vector<double> a1(k1, -6.0);
    std::vector<double> b2_local(L2);

    thetahat = zL2 / sqrtI2;

    // lower: root of f1(theta) = 0 on [left, thetahat]
    auto f1 = [&](double theta)->double {
      double zL1 = zL - theta * sqrtI1;
      for (size_t l = 0; l < k1; ++l) {
        b1[l] = (b[l + L] - w1[l] * zL1) / w2[l];
        if (!effStopping[l + L]) b1[l] = 6.0;
      }
      ListCpp probs = exitprobcpp(b1, a1, zero, t1);
      auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
      double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);
      b2_local = cache_L2.get(alphaNew);
      return zL2 - theta * sqrtI2 - b2_local[L2 - 1];
    };
    lower = brent(f1, left, thetahat, tol);

    // upper: root of f2(theta) = 0 on [thetahat, right]
    auto f2 = [&](double theta)->double {
      double zL1 = -zL + theta * sqrtI1;
      for (size_t l = 0; l < k1; ++l) {
        b1[l] = (b[l + L] - w1[l] * zL1) / w2[l];
        if (!effStopping[l + L]) b1[l] = 6.0;
      }
      ListCpp probs = exitprobcpp(b1, a1, zero, t1);
      auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
      double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);
      b2_local = cache_L2.get(alphaNew);
      return -zL2 + theta * sqrtI2 - b2_local[L2 - 1];
    };
    upper = brent(f2, thetahat, right, tol);

    // repeated p-value: more complex, uses nested root-find (brent)
    auto f = [&](double aval)->double {
      std::vector<double> u_local = cache_kMax.get(aval);

      auto g = [&](double theta)->double {
        double zL1 = zL - theta * sqrtI1;
        for (size_t l = 0; l < k1; ++l) {
          b1[l] = (u_local[l + L] - w1[l] * zL1) / w2[l];
          if (!effStopping[l + L]) b1[l] = 6.0;
        }
        ListCpp probs = exitprobcpp(b1, a1, zero, t1);
        auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
        double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);
        b2_local = cache_L2.get(alphaNew);
        return zL2 - theta * sqrtI2 - b2_local[L2 - 1];
      };

      double theta_root = brent(g, left, right, tol);
      return theta_root;
    };

    // Now locate pvalue
    if (f(0.000001) >= 0.0) {
      pvalue = 0.000001;
    } else {
      double left_p = 0.000001;
      double right_p = 0.5;
      int count = 0;
      while (true) {
        double f_right = f(right_p);
        if (f_right > 0.0) break;
        left_p = right_p;
        right_p = (left_p + 1.0) / 2.0;
        ++count;
        if (count > 18) break;
      }
      if (count <= 18) {
        pvalue = brent(f, left_p, right_p, 1.0e-6);
      } else {
        pvalue = right_p;
      }
    }
  } // end MullerSchafer branch

  DataFrameCpp df;
  df.push_back(pvalue, "pvalue");
  df.push_back(thetahat, "thetahat");
  df.push_back(cilevel, "cilevel");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");
  return df;
}


//' @title Repeated Confidence Interval After Adaptation
//' @description Obtains the repeated p-value, conservative point estimate,
//' and repeated confidence interval for an adaptive group sequential trial.
//'
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of
//'   the primary trial.
//' @param IMax The maximum information of the primary trial.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to true
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping for the primary trial.
//' @param alpha The significance level of the primary trial.
//'   Defaults to 0.025.
//' @param typeAlphaSpending The type of alpha spending for the primary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param Lc The look of interest in the integrated trial.
//' @param zLc The z-test statistic at the look of the integrated trial.
//' @param INew The maximum information of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial up to look \code{L2}.
//'   Defaults to true if left unspecified.
//' @param typeAlphaSpendingNew The type of alpha spending for the secondary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   up to look \code{L2}. Defaults to missing, in which case, it is
//'   the same as \code{informationRatesNew}.
//'
//' @return A data frame with the following variables:
//'
//' * \code{pvalue}: Repeated p-value for rejecting the null hypothesis.
//'
//' * \code{thetahat}: Point estimate of the parameter.
//'
//' * \code{cilevel}: Confidence interval level.
//'
//' * \code{lower}: Lower bound of repeated confidence interval.
//'
//' * \code{upper}: Upper bound of repeated confidence interval.
//'
//' @details
//' If \code{typeAlphaSpendingNew} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
//' \code{"none"}, then \code{informationRatesNew}, \code{efficacyStoppingNew},
//' and \code{spendingTimeNew} must be of full length \code{kNew}, and
//' \code{informationRatesNew} and \code{spendingTimeNew} must end with 1.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Cyrus R. Mehta, Peter Bauer, Martin Posch and Werner Brannath.
//' Repeated confidence intervals for adaptive group sequential trials.
//' Stat Med. 2007;26:5422–5433.
//'
//' @seealso \code{\link{adaptDesign}}
//'
//' @examples
//' # two-arm randomized clinical trial with a normally distributed endpoint
//' # 90% power to detect mean difference of 15 with a standard deviation of 50
//' # Design the Stage I Trial with 3 looks and Lan-DeMets O'Brien-Fleming type
//' # spending function
//' delta <- 15
//' sigma <- 50
//'
//' (des1 <- getDesignMeanDiff(
//'   beta = 0.1, meanDiff = delta, stDev = sigma,
//'   kMax = 3, alpha = 0.025, typeAlphaSpending = "sfOF"
//' ))
//'
//' s1 <- des1$byStageResults$informationRates
//' b1 <- des1$byStageResults$efficacyBounds
//' n <- des1$overallResults$numberOfSubjects
//'
//' # Monitoring the Stage I Trial
//' L <- 1
//' nL <- des1$byStageResults$numberOfSubjects[L]
//' deltahat <- 8
//' sigmahat <- 55
//' sedeltahat <- sigmahat * sqrt( 4 / nL)
//' zL <- deltahat / sedeltahat
//'
//' # Making an Adaptive Change: Stage I to Stage II
//' # revised clinically meaningful difference downward to 10 power the study
//' # retain the standard deviation at the design stage
//' # Muller & Schafer (2001) method to design the secondary trial
//' # with 2 looks and Lan-DeMets Pocock type spending function
//' # re-estimate sample size to reach 90% conditional power
//' deltaNew <- 10
//'
//' (des2 <- adaptDesign(
//'   betaNew = 0.1, L = L, zL = zL, theta = deltaNew,
//'   IMax = n / (4 * sigma^2), kMax = 3, informationRates = s1,
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   MullerSchafer = TRUE, kNew = 2, typeAlphaSpendingNew = "sfP"
//' ))
//'
//' INew <- des2$secondaryTrial$overallResults$information
//' (nNew <- ceiling(INew * 4 * sigma^2))
//' (nTotal <- nL + nNew)
//'
//' # Monitoring the Integrated Trial
//' s2 <- des2$secondaryTrial$byStageResults$informationRates
//'
//' Lc <- 2
//' deltahatc <- 9.5
//' sigmahatc <- 52.759
//' L2 <- Lc - L
//' nL2 <-  nNew * s2[L2]
//' nc <- nL + nL2
//' sedeltahatc <- sigmahatc * sqrt(4 / nc)
//' zLc <- deltahatc / sedeltahatc
//' zL2 <- (zLc * sqrt(nc) - zL * sqrt(nL)) / sqrt(nL2)
//'
//' getADRCI(
//'   L = L, zL = zL, IMax = n / (4 * sigmahatc^2), kMax = 3,
//'   informationRates = s1, alpha = 0.025, typeAlphaSpending = "sfOF",
//'   MullerSchafer = TRUE, Lc = Lc, zLc = zLc,
//'   INew = nNew / (4 * sigmahatc^2), informationRatesNew = s2,
//'   typeAlphaSpendingNew = "sfP")
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame getADRCI(
    const int L = NA_INTEGER,
    const double zL = NA_REAL,
    const double IMax = NA_REAL,
    const int kMax = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int Lc = NA_INTEGER,
    const double zLc = NA_REAL,
    const double INew = NA_REAL,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  auto result = getADRCIcpp(
    static_cast<size_t>(L), zL, IMax, static_cast<size_t>(kMax), infoRates,
    effStopping, critValues, alpha, typeAlphaSpending, parameterAlphaSpending,
    spendTime, MullerSchafer, static_cast<size_t>(Lc), zLc, INew,
    infoRatesNew, effStoppingNew, typeAlphaSpendingNew,
    parameterAlphaSpendingNew, spendTimeNew);
  return Rcpp::wrap(result);
}
