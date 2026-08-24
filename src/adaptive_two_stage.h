#pragma once

#include <cstddef>       // size_t
#include <cstdint>       // std::int64_t
#include <list>          // std::list
#include <mutex>         // std::mutex
#include <string>        // std::string
#include <unordered_map> // std::unordered_map
#include <vector>        // std::vector

#include "dataframe_list.h"
#include "multiplicity.h"
#include "utilities.h"

struct LocalPValues {
  std::vector<size_t> inthyp_idx;
  IntMatrix inthyp;
  std::vector<double> pinter;
};

struct PCStage1Result {
  std::vector<size_t> stg1_elemhyp_r_idx;
  std::vector<size_t> stg1_inthyp_nr_idx;
  IntMatrix inthyp;
};

struct PCStage2Result {
  std::vector<size_t> stg1_inthyp_nr_idx;
  std::vector<size_t> stg2_elemhyp_idx;
  IntMatrix inthyp;
  std::vector<double> stg1_pinter;
  std::vector<double> stg2_pinter;
  std::vector<double> comb_pinter;
  std::vector<int> rej_elem;
};

struct StageBoundaries {
  IntMatrix inthyp;
  std::vector<double> stg1_coef;
  std::vector<double> stg2_coef;
  FlatMatrix stg1_bnd;
  FlatMatrix stg2_bnd;
};

struct CER {
  std::vector<size_t> stg1_elemhyp_r_idx;
  std::vector<size_t> stg1_inthyp_nr_idx;
  IntMatrix inthyp;
  std::vector<double> CER;
};

struct AdjustedBoundaries {
  IntMatrix inthyp;
  std::vector<double> stg2_coef_new;
  FlatMatrix stg2_bnd_new;
};

LocalPValues fPCStagewiseCpp(const std::vector<double> &stg2_p,
                             const WeightMatrix &wgtmat,
                             const BoolMatrix &family, const FlatMatrix &corr,
                             const std::vector<size_t> &stg1_inthyp_nr,
                             const std::vector<size_t> &stg2_elemhyp,
                             const WeightMatrix &stg2_wgtmat,
                             const std::string &test);

LocalPValues fPCStagewiseCpp(const std::vector<double> &stg2_p,
                             const BoolMatrix &family, const FlatMatrix &corr,
                             const std::vector<size_t> &stg1_inthyp_nr,
                             const std::vector<size_t> &stg2_elemhyp,
                             const std::string &test);

LocalPValues fPCStagewiseCpp(const std::vector<double> &stg2_p,
                             const WeightMatrix &wgtmat,
                             const std::vector<size_t> &stg1_inthyp_nr,
                             const std::vector<size_t> &stg2_elemhyp,
                             const WeightMatrix &stg2_wgtmat,
                             const std::string &test);

LocalPValues fPCStagewiseCpp(const std::vector<double> &stg2_p,
                             const WeightMatrix &wgtmat,
                             const BoolMatrix &family,
                             const std::vector<size_t> &stg1_inthyp_nr,
                             const std::vector<size_t> &stg2_elemhyp,
                             const WeightMatrix &stg2_wgtmat,
                             const std::string &test);

PCStage1Result fPCStage1cpp(const LocalPValues &stg1_loc_p,
                            const double alpha1);

PCStage2Result fPCRejCpp(const LocalPValues &stg1_loc_p,
                         const LocalPValues &stg2_loc_p,
                         const std::vector<size_t> &stg1_elemhyp_r_idx,
                         const std::vector<size_t> &stg2_elemhyp_idx,
                         const double alpha, const double info_frac);

StageBoundaries fCERStageBoundCpp(const WeightMatrix &wgtmat,
                                  const BoolMatrix &family,
                                  const FlatMatrix &corr, const double alpha,
                                  const double alpha1, const double info_frac);

CER fCERCerCpp(const std::vector<double> &stg1_p, const WeightMatrix &wgtmat,
               const BoolMatrix &family, const FlatMatrix &corr,
               const double info_frac, const FlatMatrix &stg1_bnd,
               const FlatMatrix &stg2_bnd);

AdjustedBoundaries
fCERNewBoundCpp(const std::vector<double> &stg1_p, const WeightMatrix &wgtmat,
                const BoolMatrix &family, const FlatMatrix &corr,
                const std::vector<size_t> &stg1_inthyp_nr_idx,
                const std::vector<double> &CER,
                const std::vector<size_t> &stg2_elemhyp_idx,
                const WeightMatrix &stg2_wgtmat, const double info_frac_new);

std::vector<int> fCERRejCpp(const std::vector<double> &cum_p,
                            const std::vector<size_t> &stg1_elemhyp_r_idx,
                            const std::vector<size_t> &stg2_elemhyp_idx,
                            const IntMatrix &stg2_inthyp,
                            const FlatMatrix &stg2_bnd_new);
