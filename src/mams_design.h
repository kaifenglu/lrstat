#pragma once

#include <cstddef>       // size_t
#include <cstdint>       // std::int64_t
#include <list>          // std::list
#include <mutex>         // std::mutex
#include <string>        // std::string
#include <unordered_map> // std::unordered_map
#include <vector>        // std::vector

#include "utilities.h"

struct FlatMatrix;
struct ListCpp;

ListCpp exitprob_mams_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t kMax,
    const FlatMatrix& b,
    const FlatMatrix& a,
    const std::vector<double>& I);

ListCpp exitprob_mams_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t kMax,
    const FlatMatrix& b,
    const std::vector<double>& I);

std::vector<double> getBound_mams_cpp(
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

ListCpp getPower_mams(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const double alpha,
    const size_t kMax,
    const std::vector<double>& critValues,
    const std::vector<double>& I,
    const std::string& bsf,
    const double bsfpar,
    const std::vector<double>& st,
    const std::vector<unsigned char>& futStopping,
    const double IL,
    const std::vector<double>& zL);

ListCpp getDesign_mams_cpp(
    const double beta,
    const double IMax,
    const std::vector<double>& theta,
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<unsigned char>& futilityStopping,
    const FlatMatrix& criticalValues,
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

ListCpp adaptDesign_mams_cpp(
    double betaNew,
    double INew,
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t L,
    const std::vector<double>& zL,
    const std::vector<double>& theta,
    const double IMax,
    const size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<unsigned char>& futilityStopping,
    const FlatMatrix& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& futilityBounds,
    const std::vector<double>& futilityCP,
    const std::vector<double>& futilityTheta,
    const std::string& typeBetaSpending,
    const double parameterBetaSpending,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
    const size_t MNew,
    const std::vector<int>& selected,
    const double rNew,
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
