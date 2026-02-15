#ifndef __MIETTINEN_NURMINEN__
#define __MIETTINEN_NURMINEN__

struct DataFrameCpp;

#include <vector>        // vector

std::vector<double> remlRiskDiff(const double n1,
                                 const double y1,
                                 const double n2,
                                 const double y2,
                                 const double riskDiffH0);

double zstatRiskDiff(const std::vector<double>& n1,
                     const std::vector<double>& y1,
                     const std::vector<double>& n2,
                     const std::vector<double>& y2,
                     const double riskDiffH0);

DataFrameCpp mnRiskDiffCIcpp(const std::vector<double>& n1,
                             const std::vector<double>& y1,
                             const std::vector<double>& n2,
                             const std::vector<double>& y2,
                             const double cilevel);

std::vector<double> remlRiskRatio(const double n1,
                                  const double y1,
                                  const double n2,
                                  const double y2,
                                  const double riskRatioH0);

double zstatRiskRatio(const std::vector<double>& n1,
                      const std::vector<double>& y1,
                      const std::vector<double>& n2,
                      const std::vector<double>& y2,
                      const double riskRatioH0);

DataFrameCpp mnRiskRatioCIcpp(const std::vector<double>& n1,
                              const std::vector<double>& y1,
                              const std::vector<double>& n2,
                              const std::vector<double>& y2,
                              const double cilevel);

std::vector<double> remlOddsRatio(const double n1,
                                  const double y1,
                                  const double n2,
                                  const double y2,
                                  const double oddsRatioH0);

double zstatOddsRatio(const std::vector<double>& n1,
                      const std::vector<double>& y1,
                      const std::vector<double>& n2,
                      const std::vector<double>& y2,
                      const double oddsRatioH0);

DataFrameCpp mnOddsRatioCIcpp(const std::vector<double>& n1,
                              const std::vector<double>& y1,
                              const std::vector<double>& n2,
                              const std::vector<double>& y2,
                              const double cilevel);

std::vector<double> remlRateDiff(const double t1,
                                 const double y1,
                                 const double t2,
                                 const double y2,
                                 const double rateDiffH0);

double zstatRateDiff(const std::vector<double>& t1,
                     const std::vector<double>& y1,
                     const std::vector<double>& t2,
                     const std::vector<double>& y2,
                     const double rateDiffH0);

DataFrameCpp mnRateDiffCIcpp(const std::vector<double>& t1,
                             const std::vector<double>& y1,
                             const std::vector<double>& t2,
                             const std::vector<double>& y2,
                             const double cilevel);

std::vector<double> remlRateRatio(const double t1,
                                  const double y1,
                                  const double t2,
                                  const double y2,
                                  const double rateRatioH0);

double zstatRateRatio(const std::vector<double>& t1,
                      const std::vector<double>& y1,
                      const std::vector<double>& t2,
                      const std::vector<double>& y2,
                      const double rateRatioH0);

DataFrameCpp mnRateRatioCIcpp(const std::vector<double>& t1,
                              const std::vector<double>& y1,
                              const std::vector<double>& t2,
                              const std::vector<double>& y2,
                              const double cilevel);

#endif // __MIETTINEN_NURMINEN__
