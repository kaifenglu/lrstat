#pragma once

#include <cstddef>       // size_t
#include <cstdint>       // std::int64_t
#include <list>          // std::list
#include <mutex>         // std::mutex
#include <string>        // std::string
#include <unordered_map> // std::unordered_map
#include <vector>        // std::vector

#include "utilities.h"

struct ListCpp;


double errorSpentcpp(
    const double t,
    const double error = 0.025,
    const std::string& sf = "sfOF",
    const double sfpar = 0.0);


ListCpp exitprobcpp(
    const std::vector<double>& b,
    const std::vector<double>& a,
    const std::vector<double>& theta,
    const std::vector<double>& I);

double f_pvalue(
    const double theta,
    const size_t L,
    const double zL,
    const std::vector<double>& b,
    const std::vector<double>& I);

std::vector<double> getBoundcpp(
    const std::size_t k,
    const std::vector<double>& informationRates,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const std::vector<unsigned char>& efficacyStopping);


class BoundCacheAlpha {
public:
  BoundCacheAlpha(
    std::size_t k,
    const std::vector<double>& infoRates,
    const std::string& asf,
    double asfpar,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendTime,
    const std::vector<unsigned char>& effStopping,
    std::size_t maxEntries = 64,
    int alphaPrecision = 12);

  std::vector<double> get(double alpha);

private:
  struct CacheEntry {
    std::vector<double> value;
    std::list<int64_t>::iterator lruIt;
  };

  int64_t discretize(double alpha) const;

  std::size_t k_;
  std::vector<double> infoRates_;
  std::string asf_;
  double asfpar_;
  std::vector<double> userAlphaSpending_;
  std::vector<double> spendTime_;
  std::vector<unsigned char> effStopping_;

  std::size_t maxEntries_;
  int alphaPrecision_;
  std::unordered_map<int64_t, CacheEntry> map_;
  std::list<int64_t> usage_;
  std::mutex mu_;
};


ListCpp getPower(
    const double alpha2,
    const size_t k2,
    const std::vector<double>& critValues2,
    const std::vector<double>& theta,
    const std::vector<double>& Ic,
    const std::string& bsf2,
    const double bsfpar2,
    const std::vector<double>& st2,
    const std::vector<unsigned char>& futStopping2,
    const std::vector<double>& wc,
    const double IL = 0,
    const double thetaL = 0,
    const double zL = 0);


ListCpp getDesigncpp(
    const double beta,
    const double IMax,
    const double theta,
    const std::size_t kMax,
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
    const std::vector<double>& spendingTime,
    const double varianceRatio = 1.0);


ListCpp getDesignEquivcpp(
    const double beta,
    const double IMax,
    const double thetaLower,
    const double thetaUpper,
    const double theta,
    const std::size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime);


ListCpp adaptDesigncpp(
    double betaNew,
    double INew,
    const size_t L,
    const double zL,
    const double theta,
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
    const std::vector<double>& userBetaSpendingNew,
    const std::vector<double>& spendingTimeNew,
    const double varianceRatio = 1.0);
