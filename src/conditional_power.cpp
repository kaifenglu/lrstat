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

#include <Rcpp.h>


// Helper function to compute conditional power
double getCPcpp(
    const double INew,
    const int L,
    const double zL,
    const std::vector<double>& theta,
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
    const std::vector<double>& spendingTimeNew,
    const double varianceRatio) {

  // Basic validations
  if (std::isnan(INew)) throw std::invalid_argument("INew must be provided");
  if (INew <= 0.0) throw std::invalid_argument("INew must be positive");
  if (L <= 0) throw std::invalid_argument("L must be a positive integer");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::any_of(theta.begin(), theta.end(),
                  [](double v){ return std::isnan(v); })) {
    throw std::invalid_argument("theta must be provided");
  }
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }

  std::size_t K = static_cast<std::size_t>(kMax);
  std::size_t KNew = static_cast<std::size_t>(kNew);

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(K);
  if (none_na(informationRates)) {
    if (informationRates.size() != K)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[K-1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (std::size_t i = 0; i < K; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(K);
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != K)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[K-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(K, 1);
  }

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != K)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping[K-1] != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(K, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds);

  if (!missingCriticalValues && criticalValues.size() != K) {
    throw std::invalid_argument("Invalid length for criticalValues");
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
    if (userAlphaSpending.size() != K)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[K-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (!(futilityBounds.size() == K - 1 || futilityBounds.size() == K)) {
      throw std::invalid_argument("Invalid length for futilityBounds");
    }
  }
  if (!missingCriticalValues && !missingFutilityBounds) {
    for (std::size_t i = 0; i < K - 1; ++i) {
      if (futilityBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
    if (futilityBounds.size() == K && futilityBounds[K-1] != criticalValues[K-1]) {
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
    if (spendingTime.size() != K)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[K-1] != 1.0)
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
    if (KNew < 1) {
      throw std::invalid_argument("kNew must be a positive integer");
    }

    // informationRatesNew: default to (1:kNew)/kNew if missing
    infoRatesNew.resize(KNew);
    if (none_na(informationRatesNew)) {
      if (informationRatesNew.size() != KNew)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew[KNew-1] != 1.0)
        throw std::invalid_argument("informationRatesNew must end with 1");
      infoRatesNew = informationRatesNew; // copy
    } else {
      for (std::size_t i = 0; i < KNew; ++i)
        infoRatesNew[i] = static_cast<double>(i+1) / static_cast<double>(KNew);
    }

    // effStoppingNew: default to all 1s if missing
    if (none_na(efficacyStoppingNew)) {
      if (efficacyStoppingNew.size() != KNew)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew[KNew-1] != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew; // copy
    } else {
      effStoppingNew.assign(KNew, 1);
    }

    // futStoppingNew: default to all 1s if missing
    if (none_na(futilityStoppingNew)) {
      if (futilityStoppingNew.size() != KNew)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew[KNew-1] != 1)
        throw std::invalid_argument("futilityStoppingNew must end with 1");
      futStoppingNew = futilityStoppingNew; // copy
    } else {
      futStoppingNew.assign(KNew, 1);
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

    if (!(bsfNew == "sfof" || bsfNew == "sfp" ||
        bsfNew == "sfkd" || bsfNew == "sfhsd" || bsfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpendingNew");
    }
    if ((bsfNew == "sfkd" || bsfNew == "sfhsd") &&
        std::isnan(parameterBetaSpendingNew)) {
      throw std::invalid_argument("Missing value for parameterBetaSpendingNew");
    }
    if (bsfNew == "sfkd" && parameterBetaSpendingNew <= 0.0) {
      throw std::invalid_argument(
          "parameterBetaSpendingNew must be positive for sfKD");
    }

    // spendingTimeNew: default to informationRatesNew if missing
    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() != KNew)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew[KNew-1] != 1.0)
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
  std::vector<double> l(K, -6.0), zero(K, 0.0);
  std::vector<double> critValues = criticalValues;
  double alpha1 = alpha;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (K > 1 && criticalValues.size() == K) {
      bool hasNaN = false;
      for (std::size_t i = 0; i < K - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[K-1])) {
        haybittle = true;
      }
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(K);
      for (std::size_t i = 0; i < K - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[K-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[K-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  } else {
    for (std::size_t i = 0; i < K; ++i) {
      if (!effStopping[i]) critValues[i] = 6.0;
    }
    ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    alpha1 = std::accumulate(v.begin(), v.end(), 0.0);
  }

  // obtain futility bounds for the primary trial
  std::vector<double> futBounds = futilityBounds;
  if (K > 1) {
    if (missingFutilityBounds && bsf == "none") {
      futBounds = std::vector<double>(K, -6.0);
      futBounds[K-1] = critValues[K-1];
    } else if (!missingFutilityBounds && futBounds.size() == K-1) {
      futBounds.push_back(critValues[K-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // multiplier for the boundaries under the alternative hypothesis
  std::vector<double> w(K, std::sqrt(varianceRatio));
  if (!none_na(futBounds)) {
    // build delta prefix of length kMax
    std::vector<double> delta(K);
    if (theta.size() == 1) {
      std::fill(delta.begin(), delta.end(), theta[0]);
    } else if (theta.size() >= K) {
      std::memcpy(delta.data(), theta.data(), K * sizeof(double));
    } else {
      throw std::invalid_argument("Invalid length for theta");
    }

    std::vector<double> information(K);
    for (std::size_t i = 0; i < K; ++i) {
      information[i] = IMax * infoRates[i];
    }

    ListCpp out = getPower(alpha1, kMax, critValues, delta, information, bsf,
                           parameterBetaSpending, spendTime, futStopping, w);
    futBounds = out.get<std::vector<double>>("futilityBounds");
  }

  // compute transformed quantities for adaptation
  std::size_t K1 = K - L;
  std::vector<double> t1(K1), r1(K1), b1(K1), a1(K1, -6.0);
  for (std::size_t l = 0; l < K1; ++l) {
    t1[l] = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
    r1[l] = infoRates[L - 1] / infoRates[l + L];
    b1[l] = (critValues[l + L] - std::sqrt(r1[l]) * zL) / std::sqrt(1.0 - r1[l]);
    if (!effStopping[l+L]) b1[l] = 6.0;
  }

  double result = 0.0;
  if (!MullerSchafer) {
    // construct theta1 of length K1 + 1
    std::vector<double> theta1(K1 + 1);
    if (theta.size() == 1) {
      std::fill(theta1.begin(), theta1.end(), theta[0]);
    } else if (theta.size() == K + K1) {
      theta1[0] = theta[L - 1];
      for (std::size_t l = 0; l < K1; ++l) theta1[l + 1] = theta[K + l];
    } else {
      throw std::invalid_argument("Invalid length for theta");
    }

    // compute a1 (futility bounds for secondary trial)
    for (std::size_t l = 0; l < K1; ++l) {
      a1[l] = (futBounds[l + L] - std::sqrt(r1[l]) * zL) / std::sqrt(1.0 - r1[l]);
      if (!futStopping[l + L]) a1[l] = -6.0;
    }

    std::vector<double> mu(K1), I2(K1);
    for (std::size_t l = 0; l < K1; ++l) {
      double r = (IMax * infoRates[L - 1]) / (IMax * infoRates[L - 1] + INew * t1[l]);
      mu[l] = (theta1[l + 1] - r * theta1[0]) / (1.0 - r);
      I2[l] = INew * t1[l];
    }

    ListCpp probs = exitprobcpp(b1, a1, mu, I2);
    auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
    result = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);
  } else {
    // Muller-Schafer branch
    std::vector<double> theta1(KNew + 1);
    if (theta.size() == 1) {
      std::fill(theta1.begin(), theta1.end(), theta[0]);
    } else if (theta.size() == K + KNew) {
      theta1[0] = theta[L - 1];
      for (std::size_t l = 0; l < KNew; ++l) theta1[l + 1] = theta[K + l];
    } else {
      throw std::invalid_argument("Invalid length for theta");
    }

    // conditional type I error
    ListCpp probs = exitprobcpp(b1, a1, zero, t1);
    auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
    double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);

    // efficacy boundaries for secondary trial
    auto b2 = getBoundcpp(kNew, informationRatesNew, alphaNew, asfNew,
                          parameterAlphaSpendingNew, std::vector<double>{},
                          spendTimeNew, effStoppingNew);

    // conditional power
    std::vector<double> mu(KNew), I2(KNew);
    for (std::size_t l = 0; l < KNew; ++l) {
      double r = (IMax * infoRates[L - 1]) /
        (IMax * infoRates[L - 1] + INew * informationRatesNew[l]);
      mu[l] = (theta1[l + 1] - r * theta1[0]) / (1.0 - r);
      I2[l] = INew * informationRatesNew[l];
    }

    ListCpp out = getPower(alphaNew, kNew, b2, mu, I2, bsfNew,
                           parameterBetaSpendingNew, spendTimeNew,
                           futStoppingNew, w);
    result = out.get<double>("power");
  }

  return result;
}


//' @title Conditional Power Allowing for Varying Parameter Values
//' @description Obtains the conditional power for specified incremental
//' information given the interim results, parameter values, and
//' data-dependent changes in the error spending function, as well as the
//' number and spacing of interim looks.
//'
//' @param INew The maximum information of the secondary trial.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of
//'   the primary trial.
//' @param theta A scalar or a vector of parameter values of
//'   length \code{kMax + kMax - L} if \code{MullerSchafer = FALSE} or
//'   length \code{kMax + kNew} if \code{MullerSchafer = TRUE}.
//' @param IMax The maximum information of the primary trial.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to true
//'   if left unspecified.
//' @param futilityStopping Indicators of whether futility stopping is
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
//'   "sfHSD" for Hwang, Shi & DeCani spending function,
//'   "user" for user defined spending, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to Delta for "WT", rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param userAlphaSpending The user defined alpha spending for the primary
//'   trial. Cumulative alpha spent up to each stage.
//' @param futilityBounds	The lower boundaries on the z-test statistic scale
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
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to true
//'   if left unspecified.
//' @param futilityStoppingNew The indicators of whether futility stopping is
//'   allowed at each look of the secondary trial. Defaults to true
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
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early futility stopping.
//'   Defaults to "none".
//' @param parameterBetaSpendingNew The parameter value of beta spending
//'   for the secondary trial. Corresponds to rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRatesNew}.
//' @param varianceRatio The ratio of the variance under H0 to the variance
//'   under H1.
//'
//' @return The conditional power given the interim results, parameter
//' values, and data-dependent design changes.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Cyrus R. Mehta and Stuart J. Pocock.
//' Adaptive increase in sample size when interim results are promising:
//' A practical guide with examples.
//' Stat Med. 2011;30:3267–3284.
//'
//' @seealso \code{\link{getDesign}}
//'
//' @examples
//'
//' # Conditional power calculation with delayed treatment effect
//'
//' # Two interim analyses have occurred with 179 and 266 events,
//' # respectively. The observed hazard ratio at the second interim
//' # look is 0.81.
//'
//' trialsdt = as.Date("2020-03-04")                       # trial start date
//' iadt = c(as.Date("2022-02-01"), as.Date("2022-11-01")) # interim dates
//' mo1 = as.numeric(iadt - trialsdt + 1)/30.4375          # interim months
//'
//' # Assume a piecewise Poisson enrollment process with a 8-month ramp-up
//' # and 521 patients were enrolled after 17.94 months
//' N = 521                   # total number of patients
//' Ta = 17.94                # enrollment duration
//' Ta1 = 8                   # assumed end of enrollment ramp-up
//' enrate = N / (Ta - Ta1/2) # enrollment rate after ramp-up
//'
//' # Assume a median survival of 16.7 months for the control group, a
//' # 5-month delay in treatment effect, and a hazard ratio of 0.7 after
//' # the delay
//' lam1 = log(2)/16.7  # control group hazard of exponential distribution
//' t1 = 5              # months of delay in treatment effect
//' hr = 0.7            # hazard ratio after delay
//' lam2 = hr*lam1      # treatment group hazard after delay
//'
//' # Assume an annual dropout rate of 5%
//' gam = -log(1-0.05)/12  # hazard for dropout
//'
//' # The original target number of events was 298 and the new target is 335
//' mo2 <- caltime(
//'   nevents = c(298, 335),
//'   allocationRatioPlanned = 1,
//'   accrualTime = seq(0, Ta1),
//'   accrualIntensity = enrate*seq(1, Ta1+1)/(Ta1+1),
//'   piecewiseSurvivalTime = c(0, t1),
//'   lambda1 = c(lam1, lam2),
//'   lambda2 = c(lam1, lam1),
//'   gamma1 = gam,
//'   gamma2 = gam,
//'   accrualDuration = Ta,
//'   followupTime = 1000)
//'
//' # expected number of events and average hazard ratios
//' (lr1 <- lrstat(
//'   time = c(mo1, mo2),
//'   accrualTime = seq(0, Ta1),
//'   accrualIntensity = enrate*seq(1, Ta1+1)/(Ta1+1),
//'   piecewiseSurvivalTime = c(0, t1),
//'   lambda1 = c(lam1, lam2),
//'   lambda2 = c(lam1, lam1),
//'   gamma1 = gam,
//'   gamma2 = gam,
//'   accrualDuration = Ta,
//'   followupTime = 1000,
//'   predictTarget = 3))
//'
//'
//' hr2 = 0.81                    # observed hazard ratio at interim 2
//' z2 = (-log(hr2))*sqrt(266/4)  # corresponding z-test statistic value
//'
//' # expected mean of -log(HR) at the original looks and the new final look
//' theta = -log(lr1$HR[c(1,2,3,4)])
//'
//' # conditional power with sample size increase
//' getCP(INew = (335 - 266)/4,
//'       L = 2, zL = z2, theta = theta,
//'       IMax = 298/4, kMax = 3,
//'       informationRates = c(179, 266, 298)/298,
//'       alpha = 0.025, typeAlphaSpending = "sfOF")
//'
//' @export
// [[Rcpp::export]]
double getCP(double INew = NA_REAL,
             const int L = NA_INTEGER,
             const double zL = NA_REAL,
             const Rcpp::NumericVector& theta = NA_REAL,
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
             const bool MullerSchafer = 0,
             const int kNew = NA_INTEGER,
             const Rcpp::NumericVector& informationRatesNew = NA_REAL,
             const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
             const Rcpp::LogicalVector& futilityStoppingNew = NA_LOGICAL,
             const std::string& typeAlphaSpendingNew = "sfOF",
             const double parameterAlphaSpendingNew = NA_REAL,
             const std::string& typeBetaSpendingNew = "none",
             const double parameterBetaSpendingNew = NA_REAL,
             const Rcpp::NumericVector& spendingTimeNew = NA_REAL,
             const double varianceRatio = 1) {

  auto theta1 = Rcpp::as<std::vector<double>>(theta);
  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto futStoppingNew = convertLogicalVector(futilityStoppingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);
  auto result = getCPcpp(
    INew, L, zL, theta1, IMax, kMax, infoRates, effStopping,
    futStopping, critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha, futBounds,
    typeBetaSpending, parameterBetaSpending, spendTime,
    MullerSchafer, kNew, infoRatesNew, effStoppingNew,
    futStoppingNew, typeAlphaSpendingNew,
    parameterAlphaSpendingNew, typeBetaSpendingNew,
    parameterBetaSpendingNew, spendTimeNew, varianceRatio);
  return result;
}
