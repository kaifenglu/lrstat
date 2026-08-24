#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "dataframe_list.h"
#include "utilities.h"

struct Graph {
  std::vector<double> w;
  FlatMatrix G;
  std::vector<std::size_t> I;
};

struct WeightMatrix {
  IntMatrix inthyp;
  FlatMatrix wgtmat;
};

struct AdjustedPValues {
  IntMatrix inthyp;
  FlatMatrix pinter;
  FlatMatrix padj;
};

inline std::vector<size_t> validateConvertIdx(const std::vector<int> &indices,
                                              int max_val, const char *name) {
  std::vector<size_t> result;
  result.reserve(indices.size());
  for (int idx : indices) {
    if (idx <= 0 || idx > max_val) {
      throw std::invalid_argument(std::string(name) +
                                  " contains invalid indices");
    }
    result.push_back(static_cast<size_t>(idx - 1));
  }
  return result;
}

Graph updateGraphcpp(const std::vector<double> &w, const FlatMatrix &G,
                     const std::vector<size_t> &I, const size_t j);

WeightMatrix fDefaultWgtmatcpp(size_t m);

BoolMatrix fDefaultFamilycpp(size_t m);

FlatMatrix fDefaultCorrcpp(const BoolMatrix &family);

WeightMatrix fwgtmatcpp(const std::vector<double> &w, const FlatMatrix &G);

AdjustedPValues fadjpboncpp(const FlatMatrix &p, const WeightMatrix &wgtmat);

AdjustedPValues fadjpboncpp(const FlatMatrix &p);

AdjustedPValues fadjpsimcpp(const FlatMatrix &p, const WeightMatrix &wgtmat,
                            const BoolMatrix &family);

AdjustedPValues fadjpsimcpp(const FlatMatrix &p, const BoolMatrix &family);

AdjustedPValues fadjpsimcpp(const FlatMatrix &p, const WeightMatrix &wgtmat);

AdjustedPValues fadjpduncpp(const FlatMatrix &p, const WeightMatrix &wgtmat,
                            const BoolMatrix &family, const FlatMatrix &corr);

AdjustedPValues fadjpduncpp(const FlatMatrix &p, const BoolMatrix &family,
                            const FlatMatrix &corr);

AdjustedPValues fadjpduncpp(const FlatMatrix &p, const WeightMatrix &wgtmat,
                            const BoolMatrix &family);

AdjustedPValues fadjpduncpp(const FlatMatrix &p, const BoolMatrix &family);

AdjustedPValues fadjpduncpp(const FlatMatrix &p, const WeightMatrix &wgtmat);

AdjustedPValues fadjpduncpp(const FlatMatrix &p);

FlatMatrix repeatedPValuecpp(const size_t kMax,
                             const std::string &typeAlphaSpending,
                             const double parameterAlphaSpending,
                             const double maxInformation, const FlatMatrix &p,
                             const FlatMatrix &information,
                             const FlatMatrix &spendingTime);

IntMatrix fseqboncpp(const std::vector<double> &w, const FlatMatrix &G,
                     const double alpha, const size_t kMax,
                     const std::vector<std::string> &typeAlphaSpending,
                     const std::vector<double> &parameterAlphaSpending,
                     const std::vector<double> &maxInformation,
                     const BoolMatrix &incidenceMatrix, const size_t k1,
                     const FlatMatrix &p, const FlatMatrix &information,
                     const FlatMatrix &spendingTime, const bool lookback);

FlatMatrix fstp2seqcpp(const FlatMatrix &p, const std::vector<double> &gamma,
                       const std::string &test = "hochberg",
                       const bool retest = true);

AdjustedPValues fstdmixcpp(const FlatMatrix &p, const BoolMatrix &family,
                           const BoolMatrix &serial, const BoolMatrix &parallel,
                           const std::vector<double> &gamma,
                           const std::string &test = "hommel",
                           const bool exhaust = true);

AdjustedPValues fmodmixcpp(const FlatMatrix &p, const BoolMatrix &family,
                           const BoolMatrix &serial, const BoolMatrix &parallel,
                           const std::vector<double> &gamma,
                           const std::string &test = "hommel",
                           const bool exhaust = true);

AdjustedPValues ftrunccpp(const FlatMatrix &p,
                          const std::string &test = "hommel",
                          const double gamma = 1.0);
