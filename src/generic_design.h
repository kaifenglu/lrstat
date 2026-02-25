#ifndef __GENERIC_DESIGN__
#define __GENERIC_DESIGN__

struct ListCpp;

#include <cstddef>       // size_t
#include <cstdint>       // std::int64_t
#include <list>          // std::list
#include <mutex>         // std::mutex
#include <string>        // std::string
#include <unordered_map> // std::unordered_map
#include <vector>        // std::vector

using std::size_t;


double errorSpentcpp(const double t,
                     const double error = 0.025,
                     const std::string& sf = "sfOF",
                     const double sfpar = 0.0);


ListCpp exitprobcpp(const std::vector<double>& b,
                    const std::vector<double>& a,
                    const std::vector<double>& theta,
                    const std::vector<double>& I);


std::vector<double> getBoundcpp(
    const int k,
    const std::vector<double>& informationRates,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const std::vector<unsigned char>& efficacyStopping);


class BoundCacheAlpha {
public:
  BoundCacheAlpha(int k,
                  const std::vector<double>& infoRates,
                  const std::string& asf,
                  double asfpar,
                  const std::vector<double>& userAlphaSpending,
                  const std::vector<double>& spendTime,
                  const std::vector<unsigned char>& effStopping,
                  size_t maxEntries = 64,
                  int alphaPrecision = 12);

  std::vector<double> get(double alpha);

private:
  struct CacheEntry {
    std::vector<double> value;
    std::list<int64_t>::iterator lruIt;
  };

  int64_t discretize(double alpha) const;

  int k_;
  std::vector<double> infoRates_;
  std::string asf_;
  double asfpar_;
  std::vector<double> userAlphaSpending_;
  std::vector<double> spendTime_;
  std::vector<unsigned char> effStopping_;

  size_t maxEntries_;
  int alphaPrecision_;
  std::unordered_map<int64_t, CacheEntry> map_;
  std::list<int64_t> usage_;
  std::mutex mu_;
};


ListCpp getPower(
    const double alpha,
    const int kMax,
    const std::vector<double>& b,
    const std::vector<double>& theta,
    const std::vector<double>& I,
    const std::string& bsf,
    const double bsfpar,
    const std::vector<double>& st,
    const std::vector<unsigned char>& futilityStopping,
    const std::vector<double>& w);


ListCpp getDesigncpp(
    const double beta,
    const double IMax,
    const double theta,
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
    const std::vector<double>& userBetaSpending,
    const std::vector<double>& spendingTime,
    const double varianceRatio = 1.0);


ListCpp getDesignEquivcpp(
    const double beta,
    const double IMax,
    const double thetaLower,
    const double thetaUpper,
    const double theta,
    const int kMax,
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
    const int L,
    const double zL,
    const double theta,
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
    const std::vector<double>& userBetaSpendingNew,
    const std::vector<double>& spendingTimeNew,
    const double varianceRatio = 1.0);


#endif // __GENERIC_DESIGN__
