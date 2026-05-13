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

using std::size_t;


// Helper function to compute conditional power
std::vector<double> getCPcpp(
    const double INew,
    const size_t L,
    const double zL,
    const std::vector<double>& theta,
    const double IMax,
    const size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<unsigned char>& futilityStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& futilityBounds,
    const std::vector<double>& futilityCP,
    const std::vector<double>& futilityTheta,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
    const size_t kNew,
    const std::vector<double>& informationRatesNew,
    const std::vector<unsigned char>& efficacyStoppingNew,
    const std::vector<unsigned char>& futilityStoppingNew,
    const std::string& typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double>& futilityBoundsInt,
    const std::vector<double>& futilityCPInt,
    const std::vector<double>& futilityThetaInt,
    const std::string& typeBetaSpendingNew,
    const double parameterBetaSpendingNew,
    const std::vector<double>& spendingTimeNew,
    const double varianceRatio) {

  // Basic validations
  if (std::isnan(INew)) throw std::invalid_argument("INew must be provided");
  if (INew <= 0.0) throw std::invalid_argument("INew must be positive");
  if (L <= 0) throw std::invalid_argument("L must be a positive integer");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (!none_na(theta)) throw std::invalid_argument("theta must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }


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
    infoRates = informationRates; // copy
  } else {
    for (size_t i = 0; i < kMax; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(kMax);
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping.back() != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds)
    && !none_na(futilityCP) && !none_na(futilityTheta);

  if (!missingCriticalValues && criticalValues.size() != kMax) {
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
    if (userAlphaSpending.size() != kMax)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending.back() != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (none_na(futilityBounds)) {
      if (futilityBounds.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityBounds");
      }
    } else if (none_na(futilityCP)) {
      if (futilityCP.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityCP");
      }
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (futilityCP[i] < 0.0 || futilityCP[i] > 1.0) {
          throw std::invalid_argument("futilityCP must lie in [0, 1]");
        }
      }
    } else {
      if (futilityTheta.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityTheta");
      }
    }
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

  size_t k1 = kMax - L;
  if (!MullerSchafer) {
    infoRatesNew.resize(k1);
    for (size_t i = 0; i < k1; ++i) {
      infoRatesNew[i] =
        (infoRates[L + i] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
    }

    effStoppingNew.resize(k1);
    std::memcpy(effStoppingNew.data(), effStopping.data() + L,
                k1 * sizeof(unsigned char));

    if (none_na(futilityStoppingNew)) {
      if (futilityStoppingNew.size() != k1)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew.back() != 1)
        throw std::invalid_argument("futilityStoppingNew must end with 1");
      futStoppingNew = futilityStoppingNew; // copy
    } else {
      futStoppingNew.resize(k1);
      std::memcpy(futStoppingNew.data(), futStopping.data() + L,
                  k1 * sizeof(unsigned char));
    }

    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() != k1)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew.back() != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = infoRatesNew;
    }
  } else {
    if (kNew < 1) {
      throw std::invalid_argument("kNew must be a positive integer");
    }

    // informationRatesNew: default to (1:kNew)/kNew if missing
    infoRatesNew.resize(kNew);
    if (none_na(informationRatesNew)) {
      if (informationRatesNew.size() != kNew)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew.back() != 1.0)
        throw std::invalid_argument("informationRatesNew must end with 1");
      infoRatesNew = informationRatesNew; // copy
    } else {
      for (size_t i = 0; i < kNew; ++i)
        infoRatesNew[i] = static_cast<double>(i+1) / static_cast<double>(kNew);
    }

    // effStoppingNew: default to all 1s if missing
    if (none_na(efficacyStoppingNew)) {
      if (efficacyStoppingNew.size() != kNew)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew.back() != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew; // copy
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    // futStoppingNew: default to all 1s if missing
    if (none_na(futilityStoppingNew)) {
      if (futilityStoppingNew.size() != kNew)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew.back() != 1)
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

    // spendingTimeNew: default to informationRatesNew if missing
    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() != kNew)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew.back() != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = infoRatesNew;
    }
  }

  size_t k2 = MullerSchafer ? kNew : k1;
  std::vector<double> thetav = expand1(theta, kMax + k2, "theta");

  bool missingFutilityBoundsInt = !none_na(futilityBoundsInt)
    && !none_na(futilityCPInt) && !none_na(futilityThetaInt);

  if (!missingFutilityBoundsInt) {
    if (none_na(futilityBoundsInt)) {
      if (futilityBoundsInt.size() < k2 - 1) {
        throw std::invalid_argument("Insufficient length for futilityBoundsInt");
      }
    } else if (none_na(futilityCPInt)) {
      if (futilityCPInt.size() < k2 - 1) {
        throw std::invalid_argument("Insufficient length for futilityCPInt");
      }
      for (size_t i = 0; i < k2 - 1; ++i) {
        if (futilityCPInt[i] < 0.0 || futilityCPInt[i] > 1.0) {
          throw std::invalid_argument("futilityCPInt must lie in [0, 1]");
        }
      }
    } else {
      if (futilityThetaInt.size() < k2 - 1) {
        throw std::invalid_argument("Insufficient length for futilityThetaInt");
      }
    }
  }

  if (missingFutilityBoundsInt && !(bsfNew == "sfof" || bsfNew == "sfp" ||
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

  if (varianceRatio <= 0.0) {
    throw std::invalid_argument("varianceRatio must be positive");
  }
  // ----------- End of Input Validation ----------- //

  // scratch space for probabilities and other vectors
  ExitProbResult probs;

  // obtain critical values for the primary trial
  std::vector<double> l(kMax, -8.0), zero(kMax, 0.0);
  std::vector<double> critValues = criticalValues;
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
        if (!effStopping[i]) u[i] = 8.0;
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        probs = exitprobcpp(u, l, zero, infoRates);
        double cpu = std::accumulate(probs.exitProbUpper.begin(),
                                     probs.exitProbUpper.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, -5.0, 8.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  } else {
    for (size_t i = 0; i < kMax; ++i) {
      if (!effStopping[i]) critValues[i] = 8.0;
    }
  }

  // multiplier for the boundaries under the alternative hypothesis
  std::vector<double> w(kMax, std::sqrt(varianceRatio));

  // obtain futility bounds for the primary trial
  std::vector<double> futBounds(kMax);
  if (kMax > 1) {
    if (missingFutilityBounds) {
      std::fill_n(futBounds.begin(), kMax-1, -8.0);
      futBounds[kMax-1] = critValues[kMax-1];
    } else if (!missingFutilityBounds) {
      if (none_na(futilityBounds)) {
        for (size_t i = 0; i < kMax - 1; ++i) {
          if (futilityBounds[i] > critValues[i]) {
            throw std::invalid_argument(
                "futilityBounds must lie below criticalValues");
          }
        }
        std::copy_n(futilityBounds.begin(), kMax-1, futBounds.begin());
        futBounds[kMax-1] = critValues[kMax-1];
      } else if (none_na(futilityCP)) {
        double c2 = critValues[kMax - 1] * w[kMax - 1];
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i]) *
            (c2 - std::sqrt(1 - infoRates[i]) * boost_qnorm(1 - futilityCP[i]));
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityCP values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax-1] = critValues[kMax-1];
      } else {
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i] * IMax) * futilityTheta[i] / w[i];
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityTheta values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax-1] = critValues[kMax-1];
      }
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // information for the primary trial
  std::vector<double> information1(kMax);
  for (size_t i = 0; i < kMax; ++i) information1[i] = IMax * infoRates[i];

  // compute transformed quantities for adaptation
  std::vector<double> r1(k1), b1(k1), a1(k1, -8.0), zero1(k1, 0.0);
  for (size_t i = 0; i < k1; ++i) {
    r1[i] = infoRates[L - 1] / infoRates[i + L];
    b1[i] = (critValues[i + L] - std::sqrt(r1[i]) * zL) / std::sqrt(1.0 - r1[i]);
    if (!effStoppingNew[i]) b1[i] = 8.0;
  }

  probs = exitprobcpp(b1, a1, zero1, infoRatesNew);
  auto v0 = probs.exitProbUpper;
  double alphaNew = std::accumulate(v0.begin(), v0.end(), 0.0);

  for (size_t i = 0; i < k1; ++i) {
    a1[i] = (futBounds[i + L] - std::sqrt(r1[i]) * zL) / std::sqrt(1.0 - r1[i]);
    if (!futStoppingNew[i]) a1[i] = -8.0;
  }

  std::vector<double> I1(k1);
  for (size_t i = 0; i < k1; ++i) {
    I1[i] = information1[i + L] - information1[L - 1];
  }

  std::vector<double> theta1(k1);
  for (size_t i = 0; i < k1; ++i) {
    theta1[i] = (thetav[i + L] * information1[i + L] -
      thetav[L - 1] * information1[L - 1]) / I1[i];
  }

  probs = exitprobcpp(b1, a1, theta1, I1);
  double conditionalPower = std::accumulate(probs.exitProbUpper.begin(),
                                            probs.exitProbUpper.end(), 0.0);


  // critical values for the secondary trial
  std::vector<double> b2;
  std::string asf2;
  double asfpar2;
  std::vector<double> cpu0(k2);
  if (!MullerSchafer) {
    asf2 = "user";
    asfpar2 = NaN;
    std::partial_sum(v0.begin(), v0.end(), cpu0.begin());
    b2 = b1;
  } else {
    asf2 = asfNew;
    asfpar2 = parameterAlphaSpendingNew;
    if (asf2 != "none" && asf2 != "of" && asf2 != "p" && asf2 != "wt") {
      for (size_t i = 0; i < k2; ++i) {
        cpu0[i] = errorSpentcpp(spendTimeNew[i], alphaNew, asf2, asfpar2);
      }
    }
    b2 = getBoundcpp(k2, infoRatesNew, alphaNew, asfNew,
                     parameterAlphaSpendingNew, {NaN},
                     spendTimeNew, effStoppingNew);
  }

  // futility boundaries for the secondary trial
  std::vector<double> a2(k2, -8.0);

  // critical values and futility boundaries for the integrated trial
  std::vector<double> wc(k2, std::sqrt(varianceRatio));

  std::vector<double> critValues2(k2);
  std::vector<double> futBounds2(k2, -8.0);
  std::vector<double> theta2(k2);

  double thetaL = thetav[L - 1];
  double IL = information1[L - 1];
  double zscaled = zL * std::sqrt(IL);

  std::vector<double> I2(k2), Ic(k2), sqrtI2(k2), sqrtIc(k2), wsqrtIc(k2);
  for (size_t i = 0; i < k2; ++i) {
    I2[i] = INew * infoRatesNew[i];
    Ic[i] = I2[i] + IL;
    sqrtI2[i] = std::sqrt(I2[i]);
    sqrtIc[i] = std::sqrt(Ic[i]);
    wsqrtIc[i] = wc[i] * sqrtIc[i];
    theta2[i] = (thetav[kMax + i] * Ic[i] - thetav[L - 1] * IL) / I2[i];
  }


  for (size_t i = 0; i < k2; ++i) {
    critValues2[i] = (b2[i] * sqrtI2[i] + zscaled) / wsqrtIc[i];
  }

  // now compute futility bounds if needed
  if (k2 > 1) {
    if (missingFutilityBoundsInt && bsfNew == "none") {
      std::fill_n(futBounds2.begin(), k2 - 1, -8.0);
      futBounds2[k2-1] = critValues2[k2-1];
    } else if (!missingFutilityBoundsInt) {
      if (none_na(futilityBoundsInt)) {
        for (size_t i = 0; i < k2 - 1; ++i) {
          if (futilityBoundsInt[i] > critValues2[i]) {
            throw std::invalid_argument(
                "futilityBoundsInt must lie below critical values "
                "for the integrated trial");
          }
        }
        std::copy_n(futilityBoundsInt.begin(), k2-1, futBounds2.begin());
        futBounds2[k2-1] = critValues2[k2-1];
      } else if (none_na(futilityCPInt)) {
        double c2 = critValues2[k2 - 1] * wc[k2 - 1];
        for (size_t i = 0; i < k2 - 1; ++i) {
          double q = boost_qnorm(1 - futilityCPInt[i]);
          double sc = Ic[i] / Ic[k2 - 1];
          futBounds2[i] = std::sqrt(sc) * (c2 - std::sqrt(1 - sc) * q);
          if (futBounds2[i] > critValues2[i]) {
            throw std::invalid_argument(
                "futilityCPInt values are too large to be compatible with "
                "critical values for the integrated trial");
          }
        }
        futBounds2[k2-1] = critValues2[k2-1];
      } else {
        for (size_t i = 0; i < k2 - 1; ++i) {
          futBounds2[i] = std::sqrt(Ic[i]) * futilityThetaInt[i];
          if (futBounds2[i] > critValues2[i]) {
            throw std::invalid_argument(
                "futilityThetaInt values are too large to be compatible with "
                "critical values for the integrated trial");
          }
        }
        futBounds2[k2-1] = critValues2[k2-1];
      }
    }
  } else {
    if (missingFutilityBoundsInt) {
      futBounds2 = critValues2;
    }
  }

  if (missingFutilityBoundsInt && bsfNew != "none" && k2 > 1) { // beta-spending
    auto out = getPower(
      alphaNew, k2, critValues2, theta2, Ic, bsfNew,
      parameterBetaSpendingNew, spendTimeNew, futStoppingNew,
      wc, IL, thetaL, zL);
    futBounds2 = out.futilityBounds;
  }

  for (size_t i = 0; i < k2 - 1; ++i) {
    a2[i] = (futBounds2[i] * wsqrtIc[i] - zscaled) / sqrtI2[i];
  }
  a2[k2 - 1] = b2[k2 - 1];

  probs = exitprobcpp(b2, a2, theta2, I2);
  double conditionalPowerNew = std::accumulate(probs.exitProbUpper.begin(),
                                               probs.exitProbUpper.end(), 0.0);

  std::vector<double> result = {conditionalPower, conditionalPowerNew};
  return result;
}


//' @title Conditional Power for Generic Group Sequential Design
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
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function,
//'   \code{"user"} for user defined spending, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param userAlphaSpending The user defined alpha spending for the primary
//'   trial. Cumulative alpha spent up to each stage.
//' @param futilityBounds	The lower boundaries on the z-test statistic scale
//'   for futility stopping for the primary trial. Defaults to
//'   \code{rep(-8, kMax-1)} if left unspecified.
//' @param futilityCP The conditional power-based futility bounds for the
//'   primary trial.
//' @param futilityTheta The parameter value-based futility bounds for the
//'   primary trial.
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
//' @param futilityBoundsInt The futility boundaries on the z statistic
//'   scale for new stages of the integrated trial.
//' @param futilityCPInt The conditional power-based futility bounds for
//'   new stages of the integrated trial.
//' @param futilityThetaInt The parameter value-based futility bounds for the
//'   new stages of the integrated trial.
//' @param typeBetaSpendingNew The type of beta spending for the secondary
//'   trial. One of the following:
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early futility stopping.
//'   Defaults to \code{"none"}.
//' @param parameterBetaSpendingNew The parameter value of beta spending
//'   for the secondary trial. Corresponds to \eqn{\rho} for \code{"sfKD"},
//'   and \eqn{\gamma} for \code{"sfHSD"}.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRatesNew}.
//' @param varianceRatio The ratio of the variance under H0 to the variance
//'   under H1.
//'
//' @return A vector of two conditional powers given the interim results and
//' parameter values, one without design change and the other with
//' data-dependent design changes.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Cyrus R. Mehta and Stuart J. Pocock.
//' Adaptive increase in sample size when interim results are promising:
//' A practical guide with examples.
//' Stat Med. 2011;30:3267–3284.
//'
//' @seealso \code{\link{adaptDesign}}
//'
//' @examples
//'
//' # Conditional power calculation with delayed treatment effect
//'
//' # Two interim analyses have occurred with 179 and 266 events,
//' # respectively. The observed hazard ratio at the second interim
//' # look is 0.81.
//'
//' trialsdt <- as.Date("2020-03-04")                       # trial start date
//' iadt <- c(as.Date("2022-02-01"), as.Date("2022-11-01")) # interim dates
//' mo1 <- as.numeric(iadt - trialsdt + 1)/30.4375          # interim months
//'
//' # Assume a piecewise Poisson enrollment process with a 8-month ramp-up
//' # and 521 patients were enrolled after 17.94 months
//' N <- 521                   # total number of patients
//' Ta <- 17.94                # enrollment duration
//' Ta1 <- 8                   # assumed end of enrollment ramp-up
//' enrate <- N / (Ta - Ta1/2) # enrollment rate after ramp-up
//'
//' # Assume a median survival of 16.7 months for the control group, a
//' # 5-month delay in treatment effect, and a hazard ratio of 0.7 after
//' # the delay
//' lam1 <- log(2)/16.7  # control group hazard of exponential distribution
//' t1 <- 5              # months of delay in treatment effect
//' hr <- 0.7            # hazard ratio after delay
//' lam2 <- hr*lam1      # treatment group hazard after delay
//'
//' # Assume an annual dropout rate of 5%
//' gam <- -log(1-0.05)/12  # hazard for dropout
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
//' hr2 <- 0.81                    # observed hazard ratio at interim 2
//' z2 <- (-log(hr2))*sqrt(266/4)  # corresponding z-test statistic value
//'
//' # expected mean of -log(HR) at the original looks and the new final look
//' theta <- -log(lr1$HR[c(1,2,3,4)])
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
Rcpp::NumericVector getCP(
    const double INew = NA_REAL,
    const int L = NA_INTEGER,
    const double zL = NA_REAL,
    const Rcpp::NumericVector& theta = NA_REAL,
    const double IMax = NA_REAL,
    const int kMax = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::Nullable<Rcpp::NumericVector> criticalValues = R_NilValue,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBounds = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityCP = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityTheta = R_NilValue,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int kNew = NA_INTEGER,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBoundsInt = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityCPInt = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityThetaInt = R_NilValue,
    const std::string& typeBetaSpendingNew = "none",
    const double parameterBetaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL,
    const double varianceRatio = 1) {

  auto theta1 = Rcpp::as<std::vector<double>>(theta);
  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto futStoppingNew = convertLogicalVector(futilityStoppingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  std::vector<double> critValues, futBounds, futCP, futTheta;
  if (criticalValues.isNotNull()) {
    critValues = Rcpp::as<std::vector<double>>(criticalValues);
  } else {
    critValues = std::vector<double>(1, NaN);
  }

  if (futilityBounds.isNotNull()) {
    futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  } else {
    futBounds = std::vector<double>(1, NaN);
  }

  if (futilityCP.isNotNull()) {
    futCP = Rcpp::as<std::vector<double>>(futilityCP);
  } else {
    futCP = std::vector<double>(1, NaN);
  }

  if (futilityTheta.isNotNull()) {
    futTheta = Rcpp::as<std::vector<double>>(futilityTheta);
  } else {
    futTheta = std::vector<double>(1, NaN);
  }

  std::vector<double> futBoundsInt, futCPInt, futThetaInt;
  if (futilityBoundsInt.isNotNull()) {
    futBoundsInt = Rcpp::as<std::vector<double>>(futilityBoundsInt);
  } else {
    futBoundsInt = std::vector<double>(1, NaN);
  }

  if (futilityCPInt.isNotNull()) {
    futCPInt = Rcpp::as<std::vector<double>>(futilityCPInt);
  } else {
    futCPInt = std::vector<double>(1, NaN);
  }

  if (futilityThetaInt.isNotNull()) {
    futThetaInt = Rcpp::as<std::vector<double>>(futilityThetaInt);
  } else {
    futThetaInt = std::vector<double>(1, NaN);
  }

  auto result = getCPcpp(
    INew, static_cast<size_t>(L), zL, theta1, IMax,
    static_cast<size_t>(kMax), infoRates, effStopping,
    futStopping, critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds,futCP, futTheta, spendTime,
    MullerSchafer, static_cast<size_t>(kNew), infoRatesNew,
    effStoppingNew, futStoppingNew, typeAlphaSpendingNew,
    parameterAlphaSpendingNew, futBoundsInt, futCPInt, futThetaInt,
    typeBetaSpendingNew, parameterBetaSpendingNew,
    spendTimeNew, varianceRatio);

  return Rcpp::wrap(result);
}
