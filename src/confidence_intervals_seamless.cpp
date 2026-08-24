#include "dataframe_list.h"
#include "generic_design.h"
#include "seamless_design.h"
#include "utilities.h"

#include <algorithm> // any_of, fill
#include <cctype>    // tolower
#include <cmath>     // isnan
#include <cstring>   // memcpy
#include <numeric>   // accumulate
#include <stdexcept> // invalid_argument
#include <string>    // string
#include <utility>   // pair, make_pair
#include <vector>    // vector

#include <Rcpp.h>

using std::size_t;

// Compute the p-value for two-stage seamless sequential design
// given theta, look L, observed z at look L (zL),
// number of active arms M, allocation ratio to common control r,
// whether the correlation is known,
// critical values vector b (length L), and information vector I (length L).
double f_pvalue_seamless(const double theta, const size_t M, const double r,
                         const bool corr_known, const size_t L, const double zL,
                         const std::vector<double> &b,
                         const std::vector<double> &I, const size_t rankp0) {

  // Build the components required by exitprobcpp:
  // upper vector: first L components from b, last component = zL
  // theta vector: all = theta scalar
  std::vector<double> upper(L + 1);
  std::memcpy(upper.data(), b.data(), L * sizeof(double));
  upper[L] = zL;

  std::vector<double> mu(M, theta);
  auto probs = exitprob_seamless_cpp(M, r, mu, corr_known, L, upper, I, rankp0);
  double sum_up = std::accumulate(probs.exitProbUpper.begin(),
                                  probs.exitProbUpper.end(), 0.0);
  return sum_up;
}

// Helper to compute the confidence interval at the end of a group sequential
// trial
DataFrameCpp getCI_seamless_cpp(
    const size_t M, const double r, const bool corr_known, const size_t L,
    const double zL, const double IMax,
    const std::vector<double> &informationRates,
    const std::vector<unsigned char> &efficacyStopping,
    const std::vector<double> &criticalValues, const double alpha,
    const std::string &typeAlphaSpending, const double parameterAlphaSpending,
    const std::vector<double> &spendingTime, const size_t rankp0) {

  // Basic argument checks
  if (M < 1)
    throw std::invalid_argument("M should be at least 1");
  if (rankp0 < 1 || rankp0 > M) {
    throw std::invalid_argument("rankp0 must be an integer between 1 and M");
  }
  if (rankp0 > 1 && !corr_known) {
    throw std::invalid_argument("corr_known must be true when rankp0 > 1");
  }
  if (r <= 0.0)
    throw std::invalid_argument("r should be positive");
  if (L <= 0)
    throw std::invalid_argument("L must be a positive integer");
  if (std::isnan(zL))
    throw std::invalid_argument("zL must be provided");
  if (std::isnan(IMax))
    throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0)
    throw std::invalid_argument("IMax must be positive");
  if (!none_na(informationRates))
    throw std::invalid_argument("informationRates must be provided");

  size_t kMax = informationRates.size();

  if (kMax < L + 1)
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
  if (std::isnan(alpha))
    throw std::invalid_argument("alpha must be provided");
  if (alpha < 0.00001 || alpha >= 0.5)
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");

  // typeAlphaSpending to lower-case
  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (!none_na(criticalValues) &&
      !(asf == "of" || asf == "p" || asf == "wt" || asf == "sfof" ||
        asf == "sfp" || asf == "sfkd" || asf == "sfhsd" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument(
        "parameterAlphaSpending must be positive for sfKD");
  }

  if (asf == "of" || asf == "p" || asf == "wt" || asf == "none") {
    if (informationRates.back() != 1.0) {
      throw std::invalid_argument(
          "informationRates must end with 1 for OF, P, WT, or none");
    }
    if (spendTime.back() != 1.0) {
      throw std::invalid_argument(
          "spendingTime must end with 1 for OF, P, WT, or none");
    }
  }

  // set up efficacy bounds
  std::vector<double> b;
  if (none_na(criticalValues)) {
    if (criticalValues.size() < L + 1)
      throw std::invalid_argument("Insufficient length for criticalValues");
    b = criticalValues;
  } else {
    b = getBound_seamless_cpp(M, r, corr_known, L, informationRates, alpha, asf,
                              parameterAlphaSpending, std::vector<double>{},
                              spendTime, effStopping, rankp0);
  }

  // Build full information vector I = IMax * informationRates
  std::vector<double> I(L + 1);
  for (size_t i = 0; i <= L; ++i)
    I[i] = IMax * informationRates[i];

  // p-value at theta = 0
  double pvalue = f_pvalue_seamless(0.0, M, r, corr_known, L, zL, b, I, rankp0);

  double cilevel = 1.0 - 2.0 * alpha;

  // initial bracketing interval for theta: (zL +/- 8) / sqrt(I_L)
  double sqrtIL = std::sqrt(I[L]);
  double left = (zL - 8.0) / sqrtIL;
  double right = (zL + 8.0) / sqrtIL;
  double tol = 1.0e-6;

  // median-unbiased estimate thetahat: solve f_pvalue(theta) - 0.5 = 0
  auto f_med = [&](double theta) -> double {
    return f_pvalue_seamless(theta, M, r, corr_known, L, zL, b, I, rankp0) -
           0.5;
  };
  double thetahat = brent(f_med, left, right, tol);

  // lower bound: solve f_pvalue(theta) - (1-cilevel)/2 = 0, in [left, thetahat]
  double target_lower = (1.0 - cilevel) / 2.0;
  auto f_lower = [&](double theta) -> double {
    return f_pvalue_seamless(theta, M, r, corr_known, L, zL, b, I, rankp0) -
           target_lower;
  };
  double lower = brent(f_lower, left, thetahat, tol);

  // upper bound: solve f_pvalue(theta) - (1+cilevel)/2 = 0, in [thetahat,
  // right]
  double target_upper = (1.0 + cilevel) / 2.0;
  auto f_upper = [&](double theta) -> double {
    return f_pvalue_seamless(theta, 1, r, corr_known, L, zL, b, I, 1) -
           target_upper;
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

// [[Rcpp::export]]
Rcpp::DataFrame getCI_seamless_Rcpp(
    const int M = NA_INTEGER, const double r = 1, const bool corr_known = true,
    const int L = NA_INTEGER, const double zL = NA_REAL,
    const double IMax = NA_REAL,
    const Rcpp::NumericVector &informationRates = NA_REAL,
    const Rcpp::LogicalVector &efficacyStopping = NA_LOGICAL,
    const Rcpp::NumericVector &criticalValues = NA_REAL,
    const double alpha = 0.025, const std::string &typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector &spendingTime = NA_REAL, const int rankp0 = 1) {

  std::vector<double> infoRates(informationRates.begin(),
                                informationRates.end());
  auto effStopping = convertLogicalVector(efficacyStopping);
  std::vector<double> critValues(criticalValues.begin(), criticalValues.end());
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());

  auto result = getCI_seamless_cpp(
      static_cast<size_t>(M), r, corr_known, static_cast<size_t>(L), zL, IMax,
      infoRates, effStopping, critValues, alpha, typeAlphaSpending,
      parameterAlphaSpending, spendTime, static_cast<size_t>(rankp0));
  return Rcpp::wrap(result);
}

// Compute the backward image (J, zJ)
std::pair<size_t, double> f_bwimage_seamless(const double theta, const size_t K,
                                             const size_t L, const double zL,
                                             const std::vector<double> &b,
                                             const std::vector<double> &I,
                                             const size_t L2, const double zL2,
                                             const std::vector<double> &b2,
                                             const std::vector<double> &I2) {

  // compute astar for the adapted secondary trial
  double astar = f_pvalue(theta, L2, zL2, b2, I2);

  // prepare b1, mu, I1 for the original secondary trial
  size_t k1 = K - L;
  std::vector<double> I1(k1);
  for (size_t l = 0; l < k1; ++l) {
    I1[l] = I[l + L + 1] - I[L];
  }

  std::vector<double> b1(k1);
  std::vector<double> a1(k1, -8.0);
  std::vector<double> mu(k1, theta);
  for (size_t l = 0; l < k1; ++l) {
    double r1 = I[L] / I[l + L + 1];
    b1[l] = (b[l + L + 1] - zL * std::sqrt(r1)) / std::sqrt(1.0 - r1);
  }

  // compute exit probabilities for b1
  auto probs = exitprobcpp(b1, a1, mu, I1);

  // find interval containing astar
  std::vector<double> cpu(k1);
  std::partial_sum(probs.exitProbUpper.begin(), probs.exitProbUpper.end(),
                   cpu.begin());
  size_t j = std::min(findInterval1(astar, cpu) + 1, k1);
  size_t J = L + j; // combined stage index in primary trial numbering

  // find zJ
  double zJ;
  if (j == 1) {
    double r1 = I[L] / I[L + 1];
    zJ = boost_qnorm(1.0 - astar) * std::sqrt(1.0 - r1) + zL * std::sqrt(r1);
  } else {
    double r1 = I[L] / I[L + j];
    auto f = [&](double z) -> double {
      double zj = (z - zL * std::sqrt(r1)) / std::sqrt(1.0 - r1);
      return f_pvalue(theta, j, zj, b1, I1) - astar;
    };

    if (j < k1) {
      zJ = brent(f, b[L + j], 8.0, 1e-6);
    } else {
      double r1 = I[L] / I[L + j];
      double lo = -8.0 * std::sqrt(1.0 - r1) + zL * std::sqrt(r1);
      zJ = brent(f, lo, 8.0, 1e-6);
    }
  }

  return std::make_pair(J, zJ);
}

// compute backward p-value for adapted trial
double f_bwpvalue_seamless(const double theta, const size_t M, const double r,
                           const bool corr_known, const size_t K,
                           const size_t L, const double zL,
                           const std::vector<double> &b,
                           const std::vector<double> &I, const size_t L2,
                           const double zL2, const std::vector<double> &b2,
                           const std::vector<double> &I2, const size_t rankp0) {

  auto bw = f_bwimage_seamless(theta, K, L, zL, b, I, L2, zL2, b2, I2);
  return f_pvalue_seamless(theta, M, r, corr_known, bw.first, bw.second, b, I,
                           rankp0);
}

// Helper to compute confidence interval after the end of an adaptive trial
DataFrameCpp getADCI_seamless_cpp(
    const size_t M, const double r, const bool corr_known, const size_t L,
    const double zL, const double IMax, const size_t K,
    const std::vector<double> &informationRates,
    const std::vector<unsigned char> &efficacyStopping,
    const std::vector<double> &criticalValues, const double alpha,
    const std::string &typeAlphaSpending, const double parameterAlphaSpending,
    const std::vector<double> &spendingTime, const bool MullerSchafer,
    const size_t Lc, const double zLc, const double INew,
    const std::vector<double> &informationRatesNew,
    const std::vector<unsigned char> &efficacyStoppingNew,
    const std::string &typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double> &spendingTimeNew, const size_t rankp0) {

  // Input validation and defaults
  if (M < 1)
    throw std::invalid_argument("M should be at least 1");
  if (rankp0 < 1 || rankp0 > M) {
    throw std::invalid_argument("rankp0 must be an integer between 1 and M");
  }
  if (rankp0 > 1 && !corr_known) {
    throw std::invalid_argument("corr_known must be true when rankp0 > 1");
  }
  if (r <= 0.0)
    throw std::invalid_argument("r should be positive");
  if (L <= 0)
    throw std::invalid_argument("L must be a positive integer");
  if (std::isnan(zL))
    throw std::invalid_argument("zL must be provided");
  if (std::isnan(IMax))
    throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0)
    throw std::invalid_argument("IMax must be positive");
  if (K <= L)
    throw std::invalid_argument("K must be greater than L");

  size_t kMax = K + 1;

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (informationRates.size() != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates.back() != 1.0)
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
    if (efficacyStopping.back() != 1)
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
    throw std::invalid_argument(
        "alpha must be provided for missing criticalValues");
  }
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 0.5)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues &&
      !(asf == "of" || asf == "p" || asf == "wt" || asf == "sfof" ||
        asf == "sfp" || asf == "sfkd" || asf == "sfhsd" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument(
        "parameterAlphaSpending must be positive for sfKD");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime.back() != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // Now handle new trial inputs
  std::vector<double> infoRatesNew;
  std::vector<unsigned char> effStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (Lc <= L)
    throw std::invalid_argument("Lc must be greater than L");
  if (std::isnan(zLc))
    throw std::invalid_argument("zLc must be provided");
  if (std::isnan(INew))
    throw std::invalid_argument("INew must be provided");
  if (INew <= 0.0)
    throw std::invalid_argument("INew must be positive");

  size_t L2 = Lc - L;

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
          asfNew == "sfof" || asfNew == "sfp" || asfNew == "sfkd" ||
          asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }

    if ((asfNew == "wt" || asfNew == "sfkd" || asfNew == "sfhsd") &&
        std::isnan(parameterAlphaSpendingNew))
      throw std::invalid_argument(
          "Missing value for parameterAlphaSpendingNew");

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
            "informationRatesNew must end with 1 for OF, P, WT, or none");
      }
      if (spendTimeNew.back() != 1.0) {
        throw std::invalid_argument(
            "spendingTimeNew must end with 1 for OF, P, WT, or none");
      }
    }
  }

  // set up efficacy bounds
  std::vector<double> b;
  double alpha1;
  if (!missingCriticalValues) {
    if (criticalValues.size() < K + 1)
      throw std::invalid_argument("Insufficient length for criticalValues");
    b = criticalValues;
    std::vector<double> zero(M, 0.0);
    auto probs = exitprob_seamless_cpp(M, r, zero, corr_known, K, b,
                                       informationRates, rankp0);
    alpha1 = std::accumulate(probs.exitProbUpper.begin(),
                             probs.exitProbUpper.end(), 0.0);
  } else {
    b = getBound_seamless_cpp(M, r, corr_known, K, infoRates, alpha, asf,
                              parameterAlphaSpending, std::vector<double>{},
                              spendTime, effStopping, rankp0);
    alpha1 = alpha;
  }

  // Primary information vector
  std::vector<double> I(kMax);
  for (size_t i = 0; i < kMax; ++i)
    I[i] = IMax * informationRates[i];

  // compute b2 and I2 for secondary trial depending on MullerSchafer
  std::vector<double> b2(L2), I2(L2);
  if (!MullerSchafer) {
    for (size_t l = 0; l < L2; ++l) {
      double t1 = (infoRates[l + L + 1] - infoRates[L]) / (1.0 - infoRates[L]);
      double r1 = infoRates[L] / infoRates[l + L + 1];
      b2[l] = (b[l + L + 1] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
      if (!effStopping[l + L + 1])
        b2[l] = 8.0;
      I2[l] = INew * t1;
    }
  } else { // conditional type I error
    size_t k1 = K - L;
    std::vector<double> t1(k1), b1(k1), a1(k1, -8.0);
    for (size_t l = 0; l < k1; ++l) {
      t1[l] = (infoRates[l + L + 1] - infoRates[L]) / (1.0 - infoRates[L]);
      double r1 = infoRates[L] / infoRates[l + L + 1];
      b1[l] = (b[l + L + 1] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
      if (!effStopping[l + L + 1])
        b1[l] = 8.0;
    }

    std::vector<double> zero1(k1, 0.0);
    auto probs = exitprobcpp(b1, a1, zero1, t1);
    double alphaNew = std::accumulate(probs.exitProbUpper.begin(),
                                      probs.exitProbUpper.end(), 0.0);

    b2 = getBoundcpp(L2, informationRatesNew, alphaNew, asfNew,
                     parameterAlphaSpendingNew, std::vector<double>{},
                     spendTimeNew, effStoppingNew);

    for (size_t l = 0; l < L2; ++l)
      I2[l] = INew * informationRatesNew[l];
  }

  // confidence level
  double cilevel = 1.0 - 2.0 * alpha1;

  // compute pvalue under theta=0 using f_bwpvalue
  double zL2 = (zLc * std::sqrt(I[L] + I2[L2 - 1]) - zL * std::sqrt(I[L])) /
               std::sqrt(I2[L2 - 1]);
  double pvalue = f_bwpvalue_seamless(0.0, M, r, corr_known, K, L, zL, b, I, L2,
                                      zL2, b2, I2, rankp0);

  // interval brackets and root-finding to obtain thetahat, lower, upper
  double sqrtIL = std::sqrt(I[L]);
  double left = (zL - b[L]) / sqrtIL;
  double right = (zL + b[L]) / sqrtIL;
  double tol = 1.0e-6;

  auto f_med = [&](double theta) -> double {
    return f_bwpvalue_seamless(theta, M, r, corr_known, K, L, zL, b, I, L2, zL2,
                               b2, I2, rankp0) -
           0.5;
  };
  double thetahat = brent(f_med, left, right, tol);

  double target_lower = (1.0 - cilevel) / 2.0;
  auto f_low = [&](double theta) -> double {
    return f_bwpvalue_seamless(theta, M, r, corr_known, K, L, zL, b, I, L2, zL2,
                               b2, I2, rankp0) -
           target_lower;
  };
  double lower = brent(f_low, left, thetahat, tol);

  double target_upper = (1.0 + cilevel) / 2.0;
  auto f_high = [&](double theta) -> double {
    return f_bwpvalue_seamless(theta, M, r, corr_known, K, L, zL, b, I, L2, zL2,
                               b2, I2, rankp0) -
           target_upper;
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

// [[Rcpp::export]]
Rcpp::DataFrame getADCI_seamless_Rcpp(
    const int M = NA_INTEGER, const double r = 1, const bool corr_known = true,
    const int L = NA_INTEGER, const double zL = NA_REAL,
    const double IMax = NA_REAL, const int K = NA_INTEGER,
    const Rcpp::NumericVector &informationRates = NA_REAL,
    const Rcpp::LogicalVector &efficacyStopping = NA_LOGICAL,
    const Rcpp::NumericVector &criticalValues = NA_REAL,
    const double alpha = 0.25, const std::string &typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector &spendingTime = NA_REAL,
    const bool MullerSchafer = false, const int Lc = NA_INTEGER,
    const double zLc = NA_REAL, const double INew = NA_REAL,
    const Rcpp::NumericVector &informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector &efficacyStoppingNew = NA_LOGICAL,
    const std::string &typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const Rcpp::NumericVector &spendingTimeNew = NA_REAL,
    const int rankp0 = 1) {

  std::vector<double> infoRates(informationRates.begin(),
                                informationRates.end());
  auto effStopping = convertLogicalVector(efficacyStopping);
  std::vector<double> critValues(criticalValues.begin(), criticalValues.end());
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());
  std::vector<double> infoRatesNew(informationRatesNew.begin(),
                                   informationRatesNew.end());
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  std::vector<double> spendTimeNew(spendingTimeNew.begin(),
                                   spendingTimeNew.end());

  auto result = getADCI_seamless_cpp(
      static_cast<size_t>(M), r, corr_known, static_cast<size_t>(L), zL, IMax,
      static_cast<size_t>(K), infoRates, effStopping, critValues, alpha,
      typeAlphaSpending, parameterAlphaSpending, spendTime, MullerSchafer,
      static_cast<size_t>(Lc), zLc, INew, infoRatesNew, effStoppingNew,
      typeAlphaSpendingNew, parameterAlphaSpendingNew, spendTimeNew,
      static_cast<size_t>(rankp0));
  return Rcpp::wrap(result);
}
