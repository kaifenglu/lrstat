#include "utilities.h"
#include "dataframe_list.h"

#include <Rcpp.h>

#include <algorithm>     // any_of, fill
#include <cctype>        // tolower
#include <cmath>         // isnan
#include <cstring>       // memcpy
#include <numeric>       // accumulate
#include <stdexcept>     // invalid_argument
#include <string>        // string
#include <vector>        // vector
#include <utility>       // pair, make_pair


// Compute the p-value given theta, look L, observed z at look L (zL),
// critical values vector b (length L) and information vector I (length L).
double f_pvalue(const double theta,
                const int L,
                const double zL,
                const std::vector<double>& b,
                const std::vector<double>& I) {
  // Build the vectors required by exitprobcpp:
  // upper: first L-1 from b, last = zL
  // lower: all -6.0
  // mu: all = theta
  std::vector<double> upper(L);
  if (L > 1) {
    std::memcpy(upper.data(), b.data(), (L - 1) * sizeof(double));
  }
  upper[L - 1] = zL;

  std::vector<double> lower(L, -6.0);
  std::vector<double> mu(L, theta);

  ListCpp probs = exitprobcpp(upper, lower, mu, I);
  auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
  double sum_up = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);

  return sum_up;
}


// Helper to compute the confidence interval at the end of a group sequential trial
DataFrameCpp getCIcpp(const int L,
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
  if (static_cast<int>(informationRates.size()) < L)
    throw std::invalid_argument("Insufficient length for informationRates");
  if (informationRates[0] <= 0.0)
    throw std::invalid_argument("informationRates must be positive");
  if (any_nonincreasing(informationRates))
    throw std::invalid_argument("informationRates must be increasing");
  if (informationRates[L - 1] > 1.0)
    throw std::invalid_argument("informationRates must not exceed 1");

  // efficacyStopping: if provided, validate, otherwise default to all ones
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (static_cast<int>(efficacyStopping.size()) < L)
      throw std::invalid_argument("Insufficient length for efficacyStopping");
    if (efficacyStopping[L - 1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping;
  } else {
    effStopping.assign(L, 1);
  }

  // spendingTime: if provided validate, otherwise use informationRates
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (static_cast<int>(spendingTime.size()) < L)
      throw std::invalid_argument("Insufficient length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[L - 1] > 1.0)
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

  // critical values: if not provided, compute using getBoundcpp
  std::vector<double> b;
  if (none_na(criticalValues)) {
    if (static_cast<int>(criticalValues.size()) < L)
      throw std::invalid_argument("Insufficient length for criticalValues");
    b = criticalValues;
  } else {
    b = getBoundcpp(L, informationRates, alpha, asf, parameterAlphaSpending,
                    std::vector<double>{}, spendTime, effStopping);
  }

  // Build full information vector I = IMax * informationRates
  std::vector<double> I(L);
  for (int i = 0; i < L; ++i) I[i] = IMax * informationRates[i];

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

  // lower bound: solve f_pvalue(theta) - (1-cilevel)/2 = 0, bracket [left, thetahat]
  double target_lower = (1.0 - cilevel) / 2.0;
  auto f_lower = [&](double theta)->double {
    return f_pvalue(theta, L, zL, b, I) - target_lower;
  };
  double lower = brent(f_lower, left, thetahat, tol);

  // upper bound: solve f_pvalue(theta) - (1+cilevel)/2 = 0, bracket [thetahat, right]
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
//' @param typeAlphaSpending The type of alpha spending.
//'   One of the following:
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpending The parameter value of alpha spending.
//'   Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".
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
  auto result = getCIcpp(L, zL, IMax, infoRates, effStopping,
                         critValues, alpha, typeAlphaSpending,
                         parameterAlphaSpending, spendTime);
  return Rcpp::wrap(result);
}


// Helper to compute repeated confidence interval
DataFrameCpp getRCIcpp(
    const int L,
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
  if (static_cast<int>(informationRates.size()) < L)
    throw std::invalid_argument("Insufficient length for informationRates");
  if (informationRates[0] <= 0.0)
    throw std::invalid_argument("informationRates must be positive");
  if (any_nonincreasing(informationRates))
    throw std::invalid_argument("informationRates must be increasing");
  if (informationRates[L - 1] > 1.0)
    throw std::invalid_argument("informationRates must not exceed 1");

  // efficacyStopping: if provided, validate, otherwise default to all ones
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (static_cast<int>(efficacyStopping.size()) < L)
      throw std::invalid_argument("Insufficient length for efficacyStopping");
    if (efficacyStopping[L - 1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping;
  } else {
    effStopping.assign(L, 1);
  }

  // spendingTime: if provided validate, otherwise use informationRates
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (static_cast<int>(spendingTime.size()) < L)
      throw std::invalid_argument("Insufficient length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[L - 1] > 1.0)
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

  // critical values: if not provided, compute using getBoundcpp with caching
  BoundCacheAlpha cache(L, informationRates, asf, parameterAlphaSpending,
                        std::vector<double>{}, spendTime, effStopping, 64, 12);
  std::vector<double> b;
  if (none_na(criticalValues)) {
    if (static_cast<int>(criticalValues.size()) < L)
      throw std::invalid_argument("Insufficient length for criticalValues");
    b = criticalValues;
  } else {
    b = cache.get(alpha);
  }

  // Build full information vector I = IMax * informationRates
  std::vector<double> I(L);
  for (int i = 0; i < L; ++i) I[i] = IMax * informationRates[i];

  // repeated confidence interval
  double sqrtIL = std::sqrt(I[L-1]);
  double lower = (zL - b[L-1]) / sqrtIL;
  double upper = (zL + b[L-1]) / sqrtIL;

  // point estimate is the lower bound for alpha = 0.5
  std::vector<double> u = cache.get(0.5);
  double thetahat = (zL - u[L-1]) / sqrtIL;

  // repeated p-value: alpha for which the lower bound of theta is zero:
  // Solve f(aval) = zL - u[ L-1 ] = 0 for aval in (1e-6, 0.999999)
  auto f = [&](double aval)->double {
    std::vector<double> u_local = cache.get(aval);
    return zL - u_local[L-1];
  };

  double pvalue;
  const double lo_p = 0.000001;
  const double hi_p = 0.999999;
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
//' @param typeAlphaSpending The type of alpha spending.
//'   One of the following:
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpending The parameter value of alpha spending.
//'   Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".
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
  auto result = getRCIcpp(L, zL, IMax, infoRates, effStopping,
                          critValues, alpha, typeAlphaSpending,
                          parameterAlphaSpending, spendTime);
  return Rcpp::wrap(result);
}



// Compute astar = sum(exitProbUpper) for secondary trial parameters.
double f_astar(const double theta,
               const int L2,
               const double zL2,
               const std::vector<double>& b2,
               const std::vector<double>& I2) {
  std::vector<double> upper(L2);
  std::vector<double> lower(L2, -6.0);
  std::vector<double> mu(L2, theta);

  if (L2 > 1) {
    std::memcpy(upper.data(), b2.data(), (L2 - 1) * sizeof(double));
  }
  upper[L2 - 1] = zL2;

  ListCpp probs = exitprobcpp(upper, lower, mu, I2);
  auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
  double sum_up = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);
  return sum_up;
}


// Compute the backward image (J, zJ)
std::pair<int,double> f_bwimage(const double theta,
                                const int kMax,
                                const int L,
                                const double zL,
                                const std::vector<double>& b,
                                const std::vector<double>& I,
                                const int L2,
                                const double zL2,
                                const std::vector<double>& b2,
                                const std::vector<double>& I2) {
  // compute astar for the secondary trial under shifted null
  double astar = f_astar(theta, L2, zL2, b2, I2);

  int k1 = kMax - L;

  // prepare b1, a1, mu, I1 for the "shifted" secondary trial
  std::vector<double> b1(k1);
  std::vector<double> a1(k1, -6.0);
  std::vector<double> mu(k1, theta);
  std::vector<double> I1(k1);

  for (int l = 0; l < k1; ++l) {
    const double ratio = I[L-1] / I[l + L];
    b1[l] = (b[l + L] - std::sqrt(ratio) * zL) / std::sqrt(1.0 - ratio);
    I1[l] = I[l + L] - I[L - 1];
  }

  // compute exit probabilities for b1 / a1
  ListCpp probs = exitprobcpp(b1, a1, mu, I1);
  auto pu = probs.get<std::vector<double>>("exitProbUpper");

  // cumulative probabilities p[0]=0, p[1]=pu[0], ...
  std::vector<double> p(k1 + 1, 0.0);
  for (int l = 0; l < k1; ++l) p[l + 1] = p[l] + pu[l];

  // find interval containing astar (j such that p[j-1] <= astar < p[j])
  int j = findInterval3({astar}, p)[0];

  // find z1j
  double z1j;
  if (j == 1) {
    z1j = boost_qnorm(1.0 - astar);
  } else {
    // root find for z in stagewise exit probability difference
    std::vector<double> upper(j);
    std::vector<double> lower(j, -6.0);
    std::vector<double> mu(j, theta);
    if (j > 1) std::memcpy(upper.data(), b1.data(), (j - 1) * sizeof(double));

    auto f = [&](double z)->double {
      upper[j - 1] = z;
      ListCpp probs = exitprobcpp(upper, lower, mu, I1);
      auto v = probs.get<std::vector<double>>("exitProbUpper");
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - astar;
    };

    z1j = brent(f, -6.0, 6.0, 1e-6);
  }

  int J = L + j; // combined stage index in primary trial numbering
  double ratio = I[L-1] / I[J-1];
  double zJ = std::sqrt(ratio) * zL + std::sqrt(1.0 - ratio) * z1j;

  return std::make_pair(J, zJ);
}


// compute backward p-value for adapted trial (wrapper)
double f_bwpvalue(const double theta,
                  const int kMax,
                  const int L,
                  const double zL,
                  const std::vector<double>& b,
                  const std::vector<double>& I,
                  const int L2,
                  const double zL2,
                  const std::vector<double>& b2,
                  const std::vector<double>& I2) {
  auto bw = f_bwimage(theta, kMax, L, zL, b, I, L2, zL2, b2, I2);
  int J = bw.first;
  double zJ = bw.second;

  std::vector<double> upper(J);
  std::vector<double> lower(J, -6.0);
  std::vector<double> mu(J, theta);

  if (J > 1) std::memcpy(upper.data(), b.data(), (J - 1) * sizeof(double));
  upper[J - 1] = zJ;

  ListCpp probs = exitprobcpp(upper, lower, mu, I);
  auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
  double sum_up = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);
  return sum_up;
}


// Helper to compute confidence interval after the end of an adaptive trial
DataFrameCpp getADCIcpp(
    const int L,
    const double zL,
    const double IMax,
    const int kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& spendingTime,
    const int L2,
    const double zL2,
    const double INew,
    const bool MullerSchafer,
    const std::vector<double>& informationRatesNew,
    const std::vector<unsigned char>& efficacyStoppingNew,
    const std::string& typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double>& spendingTimeNew) {

  // -----------------------------------------------------------------------
  // Input validation and defaults
  if (L <= 0) throw std::invalid_argument("L must be provided and positive");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (static_cast<int>(informationRates.size()) != kMax)
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
    for (int i = 0; i < kMax; ++i) {
      infoRates[i] = static_cast<double>(i + 1) / static_cast<double>(kMax);
    }
  }

  // efficacyStopping: default to all ones if not provided
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (static_cast<int>(efficacyStopping.size()) != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[kMax - 1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy if provided
  } else {
    effStopping.assign(kMax, 1);
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
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 0.5)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");
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


  // Now handle new trial inputs
  if (L2 <= 0) throw std::invalid_argument("L2 must be provided and positive");
  if (std::isnan(zL2)) throw std::invalid_argument("zL2 must be provided");
  if (std::isnan(INew)) throw std::invalid_argument("INew must be provided");
  if (INew <= 0.0) throw std::invalid_argument("INew must be positive");

  std::vector<unsigned char> effStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (MullerSchafer) {
    if (none_na(informationRatesNew)) {
      if (static_cast<int>(informationRatesNew.size()) != L2)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew[L2 - 1] > 1.0)
        throw std::invalid_argument("informationRatesNew must not exceed 1");
    } else {
      throw std::invalid_argument(
          "informationRatesNew must be provided for MullerSchafer");
    }

    if (none_na(efficacyStoppingNew)) {
      if (static_cast<int>(efficacyStoppingNew.size()) != L2)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew[L2 - 1] != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew;
    } else {
      effStoppingNew.assign(L2, 1);
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
      if (static_cast<int>(spendingTimeNew.size()) != L2)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew[L2 - 1] > 1.0)
        throw std::invalid_argument("spendingTimeNew must not exceed 1");
      spendTimeNew = spendingTimeNew;
    } else {
      spendTimeNew = informationRatesNew;
    }
  }


  // -----------------------------------------------------------------------
  // obtain critical values for the primary trial
  std::vector<double> l(kMax, -6.0), zero(kMax, 0.0);
  std::vector<double> b = criticalValues;
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

      b[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      b = getBoundcpp(kMax, infoRates, alpha, asf, parameterAlphaSpending,
                      std::vector<double>{}, spendTime, effStopping);
    }
  } else {
    for (int i = 0; i < kMax; ++i) {
      if (!effStopping[i]) b[i] = 6.0;
    }
    ListCpp probs = exitprobcpp(b, l, zero, infoRates);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    alpha1 = std::accumulate(v.begin(), v.end(), 0.0);
  }


  // Primary information vector
  std::vector<double> I(kMax);
  for (int i = 0; i < kMax; ++i) I[i] = IMax * informationRates[i];

  // compute b2 and I2 for secondary trial depending on MullerSchafer
  std::vector<double> b2(L2), I2(L2);
  if (!MullerSchafer) {
    for (int l = 0; l < L2; ++l) {
      double t1 = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      double r1 = infoRates[L - 1] / infoRates[l + L];
      b2[l] = (b[l + L] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
      if (!effStopping[l + L]) b2[l] = 6.0;
      I2[l] = INew * t1;
    }
  } else { // conditional type I error
    int k1 = kMax - L;
    std::vector<double> t1(k1), r1(k1), b1(k1), a1(k1, -6.0);
    for (int l = 0; l < k1; ++l) {
      t1[l] = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      r1[l] = infoRates[L - 1] / infoRates[l + L];
      b1[l] = (b[l + L] - std::sqrt(r1[l]) * zL) / std::sqrt(1.0 - r1[l]);
      if (!effStopping[l + L]) b1[l] = 6.0;
    }
    ListCpp probs = exitprobcpp(b1, a1, zero, t1);
    auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
    double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);

    b2 = getBoundcpp(L2, informationRatesNew, alphaNew, asfNew,
                     parameterAlphaSpendingNew, std::vector<double>{},
                     spendTimeNew, effStoppingNew);

    for (int l = 0; l < L2; ++l) I2[l] = INew * informationRatesNew[l];
  }

  // confidence level
  double cilevel = 1.0 - 2.0 * alpha1;

  // compute pvalue under theta=0 using f_bwpvalue
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
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to Delta for "WT", rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @param L2 The termination look of the secondary trial.
//' @param zL2 The z-test statistic at the termination look of the
//'   secondary trial.
//' @param INew The maximum information of the secondary trial.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param informationRatesNew The spacing of looks of the secondary trial
//'   up to look \code{L2}.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial up to look \code{L2}.
//'   Defaults to true if left unspecified.
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
//'   for the secondary trial. Corresponds to Delta for "WT",
//'   rho for "sfKD", and gamma for "sfHSD".
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
//' # Muller & Schafer (2001) method to design the secondary trial:
//' des2 = adaptDesign(
//'   betaNew = 0.2, L = L, zL = zL, theta = 5,
//'   kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpendingNew = "sfHSD",
//'   parameterAlphaSpendingNew = -2)
//'
//' n2 = ceiling(des2$secondaryTrial$overallResults$information*4*20^2)
//' ns = round(n2*(1:3)/3)
//'  (des2 = adaptDesign(
//'    INew = n2/(4*20^2), L = L, zL = zL, theta = 5,
//'    kMax = 3, informationRates = t,
//'    alpha = 0.05, typeAlphaSpending = "sfHSD",
//'    parameterAlphaSpending = -4,
//'    MullerSchafer = TRUE,
//'    kNew = 3, informationRatesNew = ns/n2,
//'    typeAlphaSpendingNew = "sfHSD",
//'    parameterAlphaSpendingNew = -2))
//'
//' # termination at the second look of the secondary trial
//' L2 = 2
//' delta2 = 6.86
//' sigma2 = 21.77
//' zL2 = delta2/sqrt(4/197*sigma2^2)
//'
//' t2 = des2$secondaryTrial$byStageResults$informationRates[1:L2]
//'
//' # confidence interval
//' getADCI(L = L, zL = zL,
//'         IMax = n/(4*sigma1^2), kMax = 3,
//'         informationRates = t,
//'         alpha = 0.05, typeAlphaSpending = "sfHSD",
//'         parameterAlphaSpending = -4,
//'         L2 = L2, zL2 = zL2,
//'         INew = n2/(4*sigma2^2),
//'         MullerSchafer = TRUE,
//'         informationRatesNew = t2,
//'         typeAlphaSpendingNew = "sfHSD",
//'         parameterAlphaSpendingNew = -2)
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
    const int L2 = NA_INTEGER,
    const double zL2 = NA_REAL,
    const double INew = NA_REAL,
    const bool MullerSchafer = 0,
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
  auto result = getADCIcpp(L, zL, IMax, kMax, infoRates, effStopping,
                           critValues, alpha, typeAlphaSpending,
                           parameterAlphaSpending, spendTime,
                           L2, zL2, INew, MullerSchafer, infoRatesNew,
                           effStoppingNew, typeAlphaSpendingNew,
                           parameterAlphaSpendingNew, spendTimeNew);
  return Rcpp::wrap(result);
}


// Helper to calculate repeated confidence interval after adaptation
DataFrameCpp getADRCIcpp(
    const int L,
    const double zL,
    const double IMax,
    const int kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& spendingTime,
    const int L2,
    const double zL2,
    const double INew,
    const bool MullerSchafer,
    const std::vector<double>& informationRatesNew,
    const std::vector<unsigned char>& efficacyStoppingNew,
    const std::string& typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double>& spendingTimeNew) {

  // -----------------------------------------------------------------------
  // Input validation and defaults
  if (L <= 0) throw std::invalid_argument("L must be provided and positive");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (static_cast<int>(informationRates.size()) != kMax)
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
    for (int i = 0; i < kMax; ++i) {
      infoRates[i] = static_cast<double>(i + 1) / static_cast<double>(kMax);
    }
  }

  // efficacyStopping: default to all ones if not provided
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (static_cast<int>(efficacyStopping.size()) != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[kMax - 1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy if provided
  } else {
    effStopping.assign(kMax, 1);
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
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 0.5)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");
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


  // Now handle new trial inputs
  if (L2 <= 0) throw std::invalid_argument("L2 must be provided and positive");
  if (std::isnan(zL2)) throw std::invalid_argument("zL2 must be provided");
  if (std::isnan(INew)) throw std::invalid_argument("INew must be provided");
  if (INew <= 0.0) throw std::invalid_argument("INew must be positive");

  std::vector<unsigned char> effStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (MullerSchafer) {
    if (none_na(informationRatesNew)) {
      if (static_cast<int>(informationRatesNew.size()) != L2)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew[L2 - 1] > 1.0)
        throw std::invalid_argument("informationRatesNew must not exceed 1");
    } else {
      throw std::invalid_argument(
          "informationRatesNew must be provided for MullerSchafer");
    }

    if (none_na(efficacyStoppingNew)) {
      if (static_cast<int>(efficacyStoppingNew.size()) != L2)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew[L2 - 1] != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew;
    } else {
      effStoppingNew.assign(L2, 1);
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
      if (static_cast<int>(spendingTimeNew.size()) != L2)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew[L2 - 1] > 1.0)
        throw std::invalid_argument("spendingTimeNew must not exceed 1");
      spendTimeNew = spendingTimeNew;
    } else {
      spendTimeNew = informationRatesNew;
    }
  }


  // -----------------------------------------------------------------------
  // obtain critical values for the primary trial
  std::vector<double> l(kMax, -6.0), zero(kMax, 0.0);
  std::vector<double> b = criticalValues;
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

      b[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      b = getBoundcpp(kMax, infoRates, alpha, asf, parameterAlphaSpending,
                      std::vector<double>{}, spendTime, effStopping);
    }
  } else {
    for (int i = 0; i < kMax; ++i) {
      if (!effStopping[i]) b[i] = 6.0;
    }
    ListCpp probs = exitprobcpp(b, l, zero, infoRates);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    alpha1 = std::accumulate(v.begin(), v.end(), 0.0);
  }


  // Primary information vector
  std::vector<double> I(kMax);
  for (int i = 0; i < kMax; ++i) I[i] = IMax * informationRates[i];

  // compute b2 and I2 for secondary trial depending on MullerSchafer
  std::vector<double> b2(L2), I2(L2);
  if (!MullerSchafer) {
    for (int l = 0; l < L2; ++l) {
      double t1 = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      double r1 = infoRates[L - 1] / infoRates[l + L];
      b2[l] = (b[l + L] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
      if (!effStopping[l + L]) b2[l] = 6.0;
      I2[l] = INew * t1;
    }
  } else { // conditional type I error
    int k1 = kMax - L;
    std::vector<double> t1(k1), r1(k1), b1(k1), a1(k1, -6.0);
    for (int l = 0; l < k1; ++l) {
      t1[l] = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      r1[l] = infoRates[L - 1] / infoRates[l + L];
      b1[l] = (b[l + L] - std::sqrt(r1[l]) * zL) / std::sqrt(1.0 - r1[l]);
      if (!effStopping[l + L]) b1[l] = 6.0;
    }
    ListCpp probs = exitprobcpp(b1, a1, zero, t1);
    auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
    double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);

    b2 = getBoundcpp(L2, informationRatesNew, alphaNew, asfNew,
                     parameterAlphaSpendingNew, std::vector<double>{},
                     spendTimeNew, effStoppingNew);

    for (int l = 0; l < L2; ++l) I2[l] = INew * informationRatesNew[l];
  }


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

    lower = (c1 - b[L + L2 - 1]) / c2;
    upper = (c1 + b[L + L2 - 1]) / c2;

    // point estimate is lower bound with alpha=0.5
    int J = L + L2;
    // create cache_J for J
    std::vector<double> t_prefix(J);
    for (int i = 0; i < J; ++i) t_prefix[i] = infoRates[i];
    BoundCacheAlpha cache_J(J, t_prefix, asf, parameterAlphaSpending,
                            std::vector<double>{}, spendTime, effStopping, 64, 12);
    std::vector<double> u = cache_J.get(0.5);

    thetahat = (c1 - u[J - 1]) / c2;

    // repeated p-value: solve f(aval) = c1 - u[J-1] = 0 for aval in (1e-6, 0.999999)
    auto f_alpha = [&](double aval)->double {
      std::vector<double> u_local = cache_J.get(aval);
      return c1 - u_local[J - 1];
    };

    const double p_lo = 0.000001;
    const double p_hi = 0.999999;
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
                               std::vector<double>{}, spendTime, effStopping, 64, 12);
    // cache for L2 (used inside g)
    BoundCacheAlpha cache_L2(L2, informationRatesNew, asfNew,
                             parameterAlphaSpendingNew, std::vector<double>{},
                             spendTimeNew, effStoppingNew, 64, 12);

    // Muller-Schafer branch: more complex iterative root-finding
    double I1 = IMax * infoRates[L - 1];
    double I2 = INew * informationRatesNew[L2 - 1];
    double sqrtI1 = std::sqrt(I1);
    double sqrtI2 = std::sqrt(I2);

    int k1 = kMax - L;

    std::vector<double> t1(k1), r1(k1), w1(k1), w2(k1);
    for (int l = 0; l < k1; ++l) {
      t1[l] = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
      r1[l] = infoRates[L - 1] / infoRates[l + L];
      w1[l] = std::sqrt(r1[l]);
      w2[l] = std::sqrt(1.0 - r1[l]);
    }

    // interval for root-finding of theta
    double left = (zL - b[L - 1]) / sqrtI1;
    double right = (zL + b[L - 1]) / sqrtI1;
    double tol = 1.0e-6;

    // point estimate is the lower bound for alpha = 0.5
    std::vector<double> u = cache_kMax.get(0.5);

    std::vector<double> b1(k1);
    std::vector<double> a1(k1, -6.0);
    std::vector<double> b2_local(L2);

    // thetahat: root of f0(theta) = 0 on [left, right]
    auto f0 = [&](double theta)->double {
      double zL1 = zL - theta * sqrtI1;
      for (int l = 0; l < k1; ++l) {
        b1[l] = (u[l + L] - w1[l] * zL1) / w2[l];
        if (!effStopping[l + L]) b1[l] = 6.0;
      }
      ListCpp probs = exitprobcpp(b1, a1, zero, t1);
      auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
      double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);
      b2_local = cache_L2.get(alphaNew);
      return zL2 - theta * sqrtI2 - b2_local[L2 - 1];
    };
    thetahat = brent(f0, left, right, tol);

    // lower: root of f1(theta) = 0 on [left, thetahat]
    auto f1 = [&](double theta)->double {
      double zL1 = zL - theta * sqrtI1;
      for (int l = 0; l < k1; ++l) {
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
      for (int l = 0; l < k1; ++l) {
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
        for (int l = 0; l < k1; ++l) {
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
  df.push_back(1.0 - 2.0 * alpha1, "cilevel");
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
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to Delta for "WT", rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @param L2 The look of interest in the secondary trial.
//' @param zL2 The z-test statistic at the look of the secondary trial.
//' @param INew The maximum information of the secondary trial.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial up to look \code{L2}.
//'   Defaults to true if left unspecified.
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
//'   for the secondary trial. Corresponds to Delta for "WT",
//'   rho for "sfKD", and gamma for "sfHSD".
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
//' # Muller & Schafer (2001) method to design the secondary trial:
//' des2 = adaptDesign(
//'   betaNew = 0.2, L = L, zL = zL, theta = 5,
//'   kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpendingNew = "sfHSD",
//'   parameterAlphaSpendingNew = -2)
//'
//' n2 = ceiling(des2$secondaryTrial$overallResults$information*4*20^2)
//' ns = round(n2*(1:3)/3)
//' (des2 = adaptDesign(
//'   INew = n2/(4*20^2), L = L, zL = zL, theta = 5,
//'   kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, informationRatesNew = ns/n2,
//'   typeAlphaSpendingNew = "sfHSD",
//'   parameterAlphaSpendingNew = -2))
//'
//' # termination at the second look of the secondary trial
//' L2 = 2
//' delta2 = 6.86
//' sigma2 = 21.77
//' zL2 = delta2/sqrt(4/197*sigma2^2)
//'
//' t2 = des2$secondaryTrial$byStageResults$informationRates[1:L2]
//'
//' # repeated confidence interval
//' getADRCI(L = L, zL = zL,
//'          IMax = n/(4*sigma1^2), kMax = 3,
//'          informationRates = t,
//'          alpha = 0.05, typeAlphaSpending = "sfHSD",
//'          parameterAlphaSpending = -4,
//'          L2 = L2, zL2 = zL2,
//'          INew = n2/(4*sigma2^2),
//'          MullerSchafer = TRUE,
//'          informationRatesNew = t2,
//'          typeAlphaSpendingNew = "sfHSD",
//'          parameterAlphaSpendingNew = -2)
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
    const int L2 = NA_INTEGER,
    const double zL2 = NA_REAL,
    const double INew = NA_REAL,
    const bool MullerSchafer = 0,
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
  auto result = getADRCIcpp(L, zL, IMax, kMax, infoRates, effStopping,
                            critValues, alpha, typeAlphaSpending,
                            parameterAlphaSpending, spendTime,
                            L2, zL2, INew, MullerSchafer, infoRatesNew,
                            effStoppingNew, typeAlphaSpendingNew,
                            parameterAlphaSpendingNew, spendTimeNew);
  return Rcpp::wrap(result);
}
