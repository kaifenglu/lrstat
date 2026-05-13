#pragma once

#include <cstddef>       // size_t
#include <cstdint>       // std::int64_t
#include <list>          // std::list
#include <mutex>         // std::mutex
#include <string>        // std::string
#include <unordered_map> // std::unordered_map
#include <vector>        // std::vector

#include "utilities.h"
#include "dataframe_list.h"

struct ExitProbSeamless {
  std::vector<double> exitProbUpper;
  std::vector<double> exitProbLower;
  FlatMatrix exitProbByArmUpper;
  FlatMatrix exitProbByArmLower;
  std::vector<double> selectAsBest;
};

ExitProbSeamless exitprob_seamless_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t K,
    const std::vector<double>& b,
    const std::vector<double>& a,
    const std::vector<double>& I);

ExitProbSeamless exitprob_seamless_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t K,
    const std::vector<double>& b,
    const std::vector<double>& I);

std::vector<double> getBound_seamless_cpp(
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t k,
    const std::vector<double>& informationRates,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const std::vector<unsigned char>& efficacyStopping);

struct GetPowerSeamless {
  double power;
  std::vector<double> futilityBounds;
  ExitProbSeamless probs;
};

GetPowerSeamless getPower_seamless(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const double alpha,
    const size_t K,
    const std::vector<double>& critValues,
    const std::vector<double>& I,
    const std::string& bsf,
    const double bsfpar,
    const std::vector<double>& st,
    const std::vector<unsigned char>& futStopping);

ListCpp getDesign_seamless_cpp(
    const double beta,
    const double IMax,
    const std::vector<double>& theta,
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t K,
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
    const std::string& typeBetaSpending,
    const double parameterBetaSpending,
    const std::vector<double>& userBetaSpending,
    const std::vector<double>& spendingTime);


ListCpp adaptDesign_seamless_cpp(
    double betaNew,
    double INew,
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t L,
    const double zL,
    const double theta,
    const double IMax,
    const size_t K,
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
    const std::vector<double>& userBetaSpendingNew,
    const std::vector<double>& spendingTimeNew);
