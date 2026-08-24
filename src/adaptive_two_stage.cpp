#include "adaptive_two_stage.h"
#include "dataframe_list.h"
#include "generic_design.h"
#include "multiplicity.h"
#include "mvnormr.h"
#include "utilities.h"

#include <algorithm> // find
#include <cmath>     // sqrt
#include <stdexcept> // invalid_argument
#include <string>    // string
#include <utility>   // pair, make_pair
#include <vector>    // vector

#include <Rcpp.h>

using std::size_t;

// Helper to compute adjusted p-values for graphical approaches
LocalPValues fPCStagewiseCpp(const std::vector<double> &stg2_p,
                             const WeightMatrix &wgtmat,
                             const BoolMatrix &family, const FlatMatrix &corr,
                             const std::vector<size_t> &stg1_inthyp_nr,
                             const std::vector<size_t> &stg2_elemhyp,
                             const WeightMatrix &stg2_wgtmat,
                             const std::string &test) {

  // Normalize test string
  std::string test1 = test;
  for (char &c : test1) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  size_t nfams = family.nrow;
  size_t m = family.ncol;
  for (size_t i : stg2_elemhyp) {
    if (i >= m) {
      throw std::invalid_argument("stg2_elemhyp contains invalid indices");
    }
  }

  size_t m2 = stg2_elemhyp.size();
  size_t ntests2 = (1 << m2) - 1;

  std::vector<size_t> active_fams;
  active_fams.reserve(nfams);
  for (size_t h = 0; h < nfams; ++h) {
    bool has_hypothesis = false;
    for (size_t j2 = 0; j2 < m2; ++j2) {
      if (family(h, stg2_elemhyp[j2]) != 0) {
        has_hypothesis = true;
        break;
      }
    }
    if (has_hypothesis)
      active_fams.push_back(h);
  }

  size_t nfams2 = active_fams.size();
  BoolMatrix family2(nfams2, m2);
  for (size_t j2 = 0; j2 < m2; ++j2) {
    for (size_t h2 = 0; h2 < nfams2; ++h2) {
      family2(h2, j2) = family(active_fams[h2], stg2_elemhyp[j2]);
    }
  }

  FlatMatrix corr2(m2, m2);
  for (size_t j = 0; j < m2; ++j) {
    for (size_t i = 0; i < m2; ++i) {
      corr2(i, j) = corr(stg2_elemhyp[i], stg2_elemhyp[j]);
    }
  }

  // convert the vector p into a matrix with one row;
  FlatMatrix stg2_pmat(1, m2);
  for (size_t j = 0; j < m2; ++j) {
    stg2_pmat(0, j) = stg2_p[j];
  }

  AdjustedPValues adj_p;
  if (test1.compare(0, 3, "bon") == 0) { // Bonferroni
    adj_p = fadjpboncpp(stg2_pmat, stg2_wgtmat);
  } else if (test1.compare(0, 3, "sim") == 0) { // Simes
    adj_p = fadjpsimcpp(stg2_pmat, stg2_wgtmat, family2);
  } else if (test1.compare(0, 3, "dun") == 0) { // Dunnett
    adj_p = fadjpduncpp(stg2_pmat, stg2_wgtmat, family2, corr2);
  } else {
    throw std::invalid_argument(
        "test must be 'bonferroni', 'simes', or 'dunnett'");
  }

  size_t s = stg1_inthyp_nr.size();
  std::vector<size_t> inthyp_idx(s);
  IntMatrix inthyp2(s, m);
  std::vector<double> pinter(s);
  for (size_t index = 0; index < s; ++index) {
    size_t i = stg1_inthyp_nr[index];
    inthyp_idx[index] = i;
    for (size_t j = 0; j < m; ++j) {
      inthyp2(index, j) = wgtmat.inthyp(i, j);
    }

    // identify the corresponding index in the adjusted p-values for stage 2
    size_t stg2_index = 0;
    for (size_t j2 = 0; j2 < m2; ++j2) {
      if (wgtmat.inthyp(i, stg2_elemhyp[j2])) {
        stg2_index |= (1 << (m2 - 1 - j2));
      }
    }
    stg2_index = ntests2 - stg2_index;

    // assign the adjusted p-value for the intersection hypothesis
    pinter[index] = adj_p.pinter(0, stg2_index);
  }

  return LocalPValues{std::move(inthyp_idx), std::move(inthyp2),
                      std::move(pinter)};
}

LocalPValues fPCStagewiseCpp(const std::vector<double> &stg2_p,
                             const BoolMatrix &family, const FlatMatrix &corr,
                             const std::vector<size_t> &stg1_inthyp_nr,
                             const std::vector<size_t> &stg2_elemhyp,
                             const std::string &test) {
  return fPCStagewiseCpp(stg2_p, fDefaultWgtmatcpp(family.ncol), family, corr,
                         stg1_inthyp_nr, stg2_elemhyp,
                         fDefaultWgtmatcpp(stg2_p.size()), test);
}

LocalPValues fPCStagewiseCpp(const std::vector<double> &stg2_p,
                             const WeightMatrix &wgtmat,
                             const std::vector<size_t> &stg1_inthyp_nr,
                             const std::vector<size_t> &stg2_elemhyp,
                             const WeightMatrix &stg2_wgtmat,
                             const std::string &test) {
  const BoolMatrix family = fDefaultFamilycpp(wgtmat.inthyp.ncol);
  return fPCStagewiseCpp(stg2_p, wgtmat, family, fDefaultCorrcpp(family),
                         stg1_inthyp_nr, stg2_elemhyp, stg2_wgtmat, test);
}

LocalPValues fPCStagewiseCpp(const std::vector<double> &stg2_p,
                             const WeightMatrix &wgtmat,
                             const BoolMatrix &family,
                             const std::vector<size_t> &stg1_inthyp_nr,
                             const std::vector<size_t> &stg2_elemhyp,
                             const WeightMatrix &stg2_wgtmat,
                             const std::string &test) {
  return fPCStagewiseCpp(stg2_p, wgtmat, family, fDefaultCorrcpp(family),
                         stg1_inthyp_nr, stg2_elemhyp, stg2_wgtmat, test);
}

// [[Rcpp::export]]
Rcpp::List fPCStagewiseRcpp(const std::vector<double> &stg2_p,
                            const Rcpp::Nullable<Rcpp::List> &wgtmat,
                            const Rcpp::Nullable<Rcpp::LogicalMatrix> &family,
                            const Rcpp::Nullable<Rcpp::NumericMatrix> &corr,
                            const std::vector<int> &stg1_inthyp_nr,
                            const std::vector<int> &stg2_elemhyp,
                            const Rcpp::Nullable<Rcpp::List> &stg2_wgtmat,
                            const std::string &test = "dunnett") {
  size_t m;
  if (!wgtmat.isNull()) {
    Rcpp::List wgtmat_list(wgtmat.get());
    ListPtr wgtmat_ptr = listcpp_from_rlist(wgtmat_list);
    m = wgtmat_ptr->get<IntMatrix>("inthyp").ncol;
  } else if (!family.isNull()) {
    Rcpp::LogicalMatrix family_matrix(family.get());
    m = static_cast<size_t>(family_matrix.ncol());
  } else if (!corr.isNull()) {
    Rcpp::NumericMatrix corr_matrix(corr.get());
    m = static_cast<size_t>(corr_matrix.ncol());
  } else {
    throw std::invalid_argument(
        "At least one of wgtmat, family, or corr must be provided");
  }
  BoolMatrix family1 = fDefaultFamilycpp(m);
  if (!family.isNull()) {
    Rcpp::LogicalMatrix family_matrix(family.get());
    family1 = boolmatrix_from_Rmatrix(family_matrix);
  }
  FlatMatrix corr1 = fDefaultCorrcpp(family1);
  if (!corr.isNull()) {
    Rcpp::NumericMatrix corr_matrix(corr.get());
    corr1 = flatmatrix_from_Rmatrix(corr_matrix);
  }
  auto Jplus = validateConvertIdx(stg1_inthyp_nr, (size_t{1} << m) - 1,
                                  "stg1_inthyp_nr");
  auto I2 = validateConvertIdx(stg2_elemhyp, m, "stg2_elemhyp");
  WeightMatrix wgt_pair = fDefaultWgtmatcpp(m);
  if (!wgtmat.isNull()) {
    Rcpp::List wgtmat_list(wgtmat.get());
    ListPtr wgtmat_ptr = listcpp_from_rlist(wgtmat_list);
    wgt_pair.inthyp = wgtmat_ptr->get<IntMatrix>("inthyp");
    wgt_pair.wgtmat = wgtmat_ptr->get<FlatMatrix>("wgtmat");
  }
  WeightMatrix stg2_wgt_pair = fDefaultWgtmatcpp(stg2_p.size());
  if (!stg2_wgtmat.isNull()) {
    Rcpp::List stg2_wgtmat_list(stg2_wgtmat.get());
    ListPtr stg2_wgtmat_ptr = listcpp_from_rlist(stg2_wgtmat_list);
    stg2_wgt_pair.inthyp = stg2_wgtmat_ptr->get<IntMatrix>("inthyp");
    stg2_wgt_pair.wgtmat = stg2_wgtmat_ptr->get<FlatMatrix>("wgtmat");
  }
  auto out = fPCStagewiseCpp(stg2_p, wgt_pair, family1, corr1, Jplus, I2,
                             stg2_wgt_pair, test);
  for (size_t &i : out.inthyp_idx)
    ++i;
  ListCpp result;
  result.push_back(std::move(out.inthyp_idx), "inthyp_idx");
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.pinter), "pinter");
  return Rcpp::wrap(result);
}

PCStage1Result fPCStage1cpp(const LocalPValues &stg1_loc_p,
                            const double alpha1) {

  size_t m = stg1_loc_p.inthyp.ncol;
  size_t ntests = (1 << m) - 1;

  // identify rejected intersection hypotheses
  std::vector<unsigned char> rejint(ntests, 0);
  std::vector<size_t> Jplus;
  for (size_t i = 0; i < ntests; ++i) {
    if (stg1_loc_p.pinter[i] <= alpha1)
      rejint[i] = 1;
    if (!rejint[i])
      Jplus.push_back(i);
  }

  std::vector<unsigned char> rejind(m, 1);
  for (size_t i = 0; i < ntests; ++i) {
    if (!rejint[i]) {
      for (size_t j = 0; j < m; ++j) {
        if (stg1_loc_p.inthyp(i, j)) {
          rejind[j] = 0;
        }
      }
    }
  }

  std::vector<size_t> I1r = which(rejind);

  size_t s = Jplus.size();
  IntMatrix inthyp2(s, m);

  // Main loop over subsets
  for (size_t index = 0; index < s; ++index) {
    size_t i = Jplus[index];
    for (size_t j = 0; j < m; ++j) {
      inthyp2(index, j) = stg1_loc_p.inthyp(i, j);
    }
  } // end subsets

  return PCStage1Result{std::move(I1r), std::move(Jplus), std::move(inthyp2)};
}

// [[Rcpp::export]]
Rcpp::List fPCStage1Rcpp(const Rcpp::List &stg1_loc_p, const double alpha1) {
  ListPtr stg1_loc_p_ptr = listcpp_from_rlist(stg1_loc_p);
  LocalPValues stg1_loc_p_pair;
  stg1_loc_p_pair.inthyp = stg1_loc_p_ptr->get<IntMatrix>("inthyp");
  stg1_loc_p_pair.pinter = stg1_loc_p_ptr->get<std::vector<double>>("pinter");
  auto out = fPCStage1cpp(stg1_loc_p_pair, alpha1);
  for (size_t &i : out.stg1_elemhyp_r_idx)
    ++i;
  for (size_t &i : out.stg1_inthyp_nr_idx)
    ++i;
  ListCpp result;
  result.push_back(std::move(out.stg1_elemhyp_r_idx), "stg1_elemhyp_r_idx");
  result.push_back(std::move(out.stg1_inthyp_nr_idx), "stg1_inthyp_nr_idx");
  result.push_back(std::move(out.inthyp), "inthyp");
  return Rcpp::wrap(result);
}

PCStage2Result fPCRejCpp(const LocalPValues &stg1_loc_p,
                         const LocalPValues &stg2_loc_p,
                         const std::vector<size_t> &stg1_elemhyp_r_idx,
                         const std::vector<size_t> &stg2_elemhyp_idx,
                         const double alpha, const double info_frac) {

  size_t m = stg1_loc_p.inthyp.ncol;
  size_t m2 = stg2_elemhyp_idx.size();
  size_t s = stg2_loc_p.inthyp.nrow;

  double sqrtt = std::sqrt(info_frac);
  double sqrt1minust = std::sqrt(1.0 - info_frac);

  std::vector<double> stg1_pinter(s);
  std::vector<double> stg2_pinter(s);
  std::vector<double> comb_pinter(s);
  std::vector<int> rej_inter(s, 0);
  for (size_t index = 0; index < s; ++index) {
    size_t i = stg2_loc_p.inthyp_idx[index];
    stg1_pinter[index] = stg1_loc_p.pinter[i];
    stg2_pinter[index] = stg2_loc_p.pinter[index];
    comb_pinter[index] =
        1.0 - boost_pnorm(sqrtt * boost_qnorm(1.0 - stg1_pinter[index]) +
                          sqrt1minust * boost_qnorm(1.0 - stg2_pinter[index]));
    if (comb_pinter[index] <= alpha) {
      rej_inter[index] = 1;
    }
  }

  std::vector<int> rej_elem(m);
  for (int i : stg1_elemhyp_r_idx)
    rej_elem[i] = 1;
  for (int i : stg2_elemhyp_idx)
    rej_elem[i] = 1;
  for (size_t index = 0; index < s; ++index) {
    if (!rej_inter[index]) {
      for (size_t j2 = 0; j2 < m2; ++j2) {
        if (stg2_loc_p.inthyp(index, stg2_elemhyp_idx[j2])) {
          rej_elem[stg2_elemhyp_idx[j2]] = 0;
        }
      }
    }
  }

  PCStage2Result result;
  result.stg1_inthyp_nr_idx = stg2_loc_p.inthyp_idx;
  result.stg2_elemhyp_idx = stg2_elemhyp_idx;
  result.inthyp = stg2_loc_p.inthyp;
  result.stg1_pinter = std::move(stg1_pinter);
  result.stg2_pinter = std::move(stg2_pinter);
  result.comb_pinter = std::move(comb_pinter);
  result.rej_elem = std::move(rej_elem);
  return result;
}

// [[Rcpp::export]]
Rcpp::List fPCRejRcpp(const Rcpp::List &stg1_loc_p,
                      const Rcpp::List &stg2_loc_p,
                      const std::vector<int> &stg1_elemhyp_r_idx,
                      const std::vector<int> &stg2_elemhyp_idx,
                      const double alpha, const double info_frac) {

  ListPtr stg1_loc_p_ptr = listcpp_from_rlist(stg1_loc_p);
  LocalPValues stg1_loc_p_tuple;
  stg1_loc_p_tuple.inthyp = stg1_loc_p_ptr->get<IntMatrix>("inthyp");
  stg1_loc_p_tuple.pinter = stg1_loc_p_ptr->get<std::vector<double>>("pinter");
  int m = stg1_loc_p_tuple.inthyp.ncol;
  int ntests = stg1_loc_p_tuple.pinter.size();
  auto stg1_inthyp_idx = stg1_loc_p_ptr->get<std::vector<int>>("inthyp_idx");
  stg1_loc_p_tuple.inthyp_idx =
      validateConvertIdx(stg1_inthyp_idx, ntests, "stg1_inthyp_idx");

  ListPtr stg2_loc_p_ptr = listcpp_from_rlist(stg2_loc_p);
  LocalPValues stg2_loc_p_tuple;
  stg2_loc_p_tuple.inthyp = stg2_loc_p_ptr->get<IntMatrix>("inthyp");
  stg2_loc_p_tuple.pinter = stg2_loc_p_ptr->get<std::vector<double>>("pinter");
  auto stg2_inthyp_idx = stg2_loc_p_ptr->get<std::vector<int>>("inthyp_idx");
  stg2_loc_p_tuple.inthyp_idx =
      validateConvertIdx(stg2_inthyp_idx, ntests, "stg2_inthyp_idx");

  auto I1r = validateConvertIdx(stg1_elemhyp_r_idx, m, "stg1_elemhyp_r_idx");
  auto I2 = validateConvertIdx(stg2_elemhyp_idx, m, "stg2_elemhyp_idx");
  auto out =
      fPCRejCpp(stg1_loc_p_tuple, stg2_loc_p_tuple, I1r, I2, alpha, info_frac);
  for (size_t &i : out.stg1_inthyp_nr_idx)
    ++i;
  for (size_t &i : out.stg2_elemhyp_idx)
    ++i;
  ListCpp result;
  result.push_back(std::move(out.stg1_inthyp_nr_idx), "stg1_inthyp_nr_idx");
  result.push_back(std::move(out.stg2_elemhyp_idx), "stg2_elemhyp_idx");
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.stg1_pinter), "stg1_pinter");
  result.push_back(std::move(out.stg2_pinter), "stg2_pinter");
  result.push_back(std::move(out.comb_pinter), "comb_pinter");
  result.push_back(std::move(out.rej_elem), "rej_elem");
  return Rcpp::wrap(result);
}

// Helper to compute adjusted p-values for Dunnett-based graphical approaches
StageBoundaries fCERStageBoundCpp(const WeightMatrix &wgtmat,
                                  const BoolMatrix &family,
                                  const FlatMatrix &corr, const double alpha,
                                  const double alpha1, const double info_frac) {

  size_t m = wgtmat.inthyp.ncol;
  size_t ntests = (1 << m) - 1;
  size_t nfams = family.nrow;

  if (alpha <= 0.0 || alpha >= 1.0) {
    throw std::invalid_argument("alpha must be in (0, 1)");
  }
  if (alpha1 < 0.0 || alpha1 >= alpha) {
    throw std::invalid_argument("alpha1 must be in [0, alpha)");
  }
  if (info_frac <= 0.0 || info_frac >= 1.0) {
    throw std::invalid_argument("info_frac must be in (0, 1)");
  }

  std::vector<size_t> J, J_h;
  std::vector<std::vector<size_t>> J_h_list(nfams);

  std::vector<double> lower1;
  std::vector<double> upper1;
  std::vector<double> mean1;
  FlatMatrix sigma1;
  std::vector<std::vector<double>> lower1_list(nfams);
  std::vector<std::vector<double>> mean1_list(nfams);
  std::vector<FlatMatrix> sigma1_list(nfams);

  std::vector<double> lower2;
  std::vector<double> upper2;
  std::vector<double> mean2;
  FlatMatrix sigma2;
  std::vector<std::vector<double>> lower2_list(nfams);
  std::vector<std::vector<double>> mean2_list(nfams);
  std::vector<FlatMatrix> sigma2_list(nfams);

  double sqrtt = std::sqrt(info_frac);
  std::vector<double> cJ1(ntests), cJ2(ntests);
  FlatMatrix wcJ1(ntests, m), wcJ2(ntests, m);

  // Main loop over subsets
  for (size_t i = 0; i < ntests; ++i) {
    J.clear();
    for (size_t j = 0; j < m; ++j) {
      if (wgtmat.wgtmat(i, j) > 0.0) {
        J.push_back(j);
      }
    }

    for (size_t h = 0; h < nfams; ++h) {
      J_h.clear();
      for (size_t j : J) {
        if (family(h, j))
          J_h.push_back(j);
      }
      J_h_list[h] = J_h;

      size_t k = J_h.size();
      if (k > 0) {
        lower1.resize(k);
        lower1.assign(k, -POS_INF);
        lower1_list[h] = lower1;

        mean1.resize(k);
        mean1.assign(k, 0.0);
        mean1_list[h] = mean1;

        sigma1.resize(k, k);
        for (size_t t1 = 0; t1 < k; ++t1) {
          for (size_t t2 = 0; t2 < k; ++t2) {
            size_t j1 = J_h[t1];
            size_t j2 = J_h[t2];
            sigma1(t1, t2) = corr(j1, j2);
          }
        }
        sigma1_list[h] = sigma1;

        lower2.resize(2 * k);
        lower2.assign(2 * k, -POS_INF);
        lower2_list[h] = lower2;

        mean2.resize(2 * k);
        mean2.assign(2 * k, 0.0);
        mean2_list[h] = mean2;

        sigma2.resize(2 * k, 2 * k);
        for (size_t t1 = 0; t1 < k; ++t1) {
          for (size_t t2 = 0; t2 < k; ++t2) {
            size_t j1 = J_h[t1];
            size_t j2 = J_h[t2];
            double rho = corr(j1, j2);
            sigma2(t1, t2) = rho;
            sigma2(t1 + k, t2 + k) = rho;
            sigma2(t1, t2 + k) = sqrtt * rho;
            sigma2(t1 + k, t2) = sqrtt * rho;
          }
        }
        sigma2_list[h] = sigma2;
      }
    }

    // solve for c1
    auto f1 = [&](double c) {
      double sum1 = 0.0;
      for (size_t h = 0; h < nfams; ++h) {
        J_h = J_h_list[h];
        size_t k = J_h.size();
        if (k > 0) {
          const auto &lower1ref = lower1_list[h];
          const auto &mean1ref = mean1_list[h];
          const auto &sigma1ref = sigma1_list[h];

          upper1.resize(k);
          for (size_t t = 0; t < k; ++t) {
            size_t j = J_h[t];
            upper1[t] = boost_qnorm(1.0 - wgtmat.wgtmat(i, j) * c);
          }

          double v1 = pmvnormcpp(lower1ref, upper1, mean1ref, sigma1ref).prob;
          sum1 += (1.0 - v1);
        }
      }
      return sum1 - alpha1;
    };

    // alpha1 = 0 means no stage 1 efficacy stopping, so the boundary is 0
    cJ1[i] = (alpha1 > 0.0) ? brent(f1, 0.00001, 0.99999, 1e-5) : 0.0;

    // cJ1[i] == 0 (no stage 1 efficacy stopping) leaves the stage 1
    // dimensions unconstrained (upper bound +Inf); marginalize them out of
    // the integral instead of relying on boost_qnorm's finite clamp for
    // p == 1, which would leave a spurious residual dependence on
    // info_frac through the stage1/stage2 cross-correlation term
    const bool stg1Unconstrained = (cJ1[i] == 0.0);

    // solve for c2
    auto f2 = [&](double c) {
      double sum2 = 0.0;
      for (size_t h = 0; h < nfams; ++h) {
        J_h = J_h_list[h];
        size_t k = J_h.size();
        if (k > 0) {
          double v2;
          if (stg1Unconstrained) {
            const auto &lower1ref = lower1_list[h];
            const auto &mean1ref = mean1_list[h];
            const auto &sigma1ref = sigma1_list[h];

            upper1.resize(k);
            for (size_t t = 0; t < k; ++t) {
              size_t j = J_h[t];
              upper1[t] = boost_qnorm(1.0 - wgtmat.wgtmat(i, j) * c);
            }

            v2 = pmvnormcpp(lower1ref, upper1, mean1ref, sigma1ref).prob;
          } else {
            const auto &lower2ref = lower2_list[h];
            const auto &mean2ref = mean2_list[h];
            const auto &sigma2ref = sigma2_list[h];

            upper2.resize(2 * k);
            for (size_t t = 0; t < k; ++t) {
              size_t j = J_h[t];
              upper2[t] = boost_qnorm(1.0 - wgtmat.wgtmat(i, j) * cJ1[i]);
              upper2[t + k] = boost_qnorm(1.0 - wgtmat.wgtmat(i, j) * c);
            }

            v2 = pmvnormcpp(lower2ref, upper2, mean2ref, sigma2ref).prob;
          }
          sum2 += (1.0 - v2);
        }
      }
      return sum2 - alpha;
    };

    cJ2[i] = brent(f2, 0.00001, 0.99999, 1e-5);

    for (size_t h = 0; h < nfams; ++h) {
      J_h = J_h_list[h];
      size_t k = J_h.size();
      if (k > 0) {
        for (size_t t = 0; t < k; ++t) {
          size_t j = J_h[t];
          wcJ1(i, j) = wgtmat.wgtmat(i, j) * cJ1[i];
          wcJ2(i, j) = wgtmat.wgtmat(i, j) * cJ2[i];
        }
      }
    }
  } // end subsets

  return StageBoundaries{std::move(wgtmat.inthyp), std::move(cJ1),
                         std::move(cJ2), std::move(wcJ1), std::move(wcJ2)};
}

// [[Rcpp::export]]
Rcpp::List fCERStageBoundRcpp(const Rcpp::List &wgtmat,
                              const Rcpp::LogicalMatrix &family,
                              const Rcpp::NumericMatrix &corr,
                              const double alpha, const double alpha1,
                              const double info_frac) {
  ListPtr wgtmat_ptr = listcpp_from_rlist(wgtmat);
  WeightMatrix wgt_pair;
  wgt_pair.inthyp = wgtmat_ptr->get<IntMatrix>("inthyp");
  wgt_pair.wgtmat = wgtmat_ptr->get<FlatMatrix>("wgtmat");
  auto family1 = boolmatrix_from_Rmatrix(family);
  auto corr1 = flatmatrix_from_Rmatrix(corr);
  auto out =
      fCERStageBoundCpp(wgt_pair, family1, corr1, alpha, alpha1, info_frac);
  ListCpp result;
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.stg1_coef), "stg1_coef");
  result.push_back(std::move(out.stg2_coef), "stg2_coef");
  result.push_back(std::move(out.stg1_bnd), "stg1_bnd");
  result.push_back(std::move(out.stg2_bnd), "stg2_bnd");
  return Rcpp::wrap(result);
}

// Helper to compute adjusted p-values for Dunnett-based graphical approaches
CER fCERCerCpp(const std::vector<double> &stg1_p, const WeightMatrix &wgtmat,
               const BoolMatrix &family, const FlatMatrix &corr,
               const double info_frac, const FlatMatrix &stg1_bnd,
               const FlatMatrix &stg2_bnd) {

  size_t m = wgtmat.inthyp.ncol;
  size_t ntests = (1 << m) - 1;
  size_t nfams = family.nrow;

  if (info_frac <= 0.0 || info_frac >= 1.0) {
    throw std::invalid_argument("info_frac must be in (0, 1)");
  }

  // identify rejected intersection hypotheses
  std::vector<unsigned char> rejint(ntests, 0);
  std::vector<size_t> Jplus;
  for (size_t i = 0; i < ntests; ++i) {
    for (size_t j = 0; j < m; ++j) {
      if (wgtmat.inthyp(i, j) && stg1_p[j] <= stg1_bnd(i, j)) {
        rejint[i] = 1;
        break;
      }
    }
    if (!rejint[i])
      Jplus.push_back(i);
  }

  std::vector<unsigned char> rejind(m, 1);
  for (size_t i = 0; i < ntests; ++i) {
    if (!rejint[i]) {
      for (size_t j = 0; j < m; ++j) {
        if (wgtmat.inthyp(i, j)) {
          rejind[j] = 0;
        }
      }
    }
  }

  std::vector<size_t> I1r = which(rejind);

  std::vector<size_t> J, J_h;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> mean;
  FlatMatrix sigma;

  double sqrtt = std::sqrt(info_frac);

  size_t s = Jplus.size();
  std::vector<double> bJ(s);
  IntMatrix inthyp2(s, m);

  // Main loop over subsets
  for (size_t index = 0; index < s; ++index) {
    size_t i = Jplus[index];

    J.clear();
    for (size_t j = 0; j < m; ++j) {
      if (wgtmat.wgtmat(i, j) > 0.0) {
        J.push_back(j);
      }
    }

    double sum2 = 0.0;
    for (size_t h = 0; h < nfams; ++h) {
      J_h.clear();
      for (size_t j : J) {
        if (family(h, j))
          J_h.push_back(j);
      }

      size_t k = J_h.size();
      if (k > 0) {
        lower.resize(k);
        lower.assign(k, -POS_INF);

        sigma.resize(k, k);
        for (size_t t1 = 0; t1 < k; ++t1) {
          for (size_t t2 = 0; t2 < k; ++t2) {
            size_t j1 = J_h[t1];
            size_t j2 = J_h[t2];
            sigma(t1, t2) = (1.0 - info_frac) * corr(j1, j2);
          }
        }

        mean.resize(k);
        upper.resize(k);
        for (size_t t = 0; t < k; ++t) {
          size_t j = J_h[t];
          mean[t] = sqrtt * boost_qnorm(1.0 - stg1_p[j]);
          upper[t] = boost_qnorm(1.0 - stg2_bnd(i, j));
        }

        double v = pmvnormcpp(lower, upper, mean, sigma).prob;
        sum2 += (1.0 - v);
      }
    }

    bJ[index] = sum2;

    for (size_t j = 0; j < m; ++j) {
      inthyp2(index, j) = wgtmat.inthyp(i, j);
    }
  } // end subsets

  return CER{std::move(I1r), std::move(Jplus), std::move(inthyp2),
             std::move(bJ)};
}

// [[Rcpp::export]]
Rcpp::List fCERCerRcpp(const std::vector<double> &stg1_p,
                       const Rcpp::List &wgtmat,
                       const Rcpp::LogicalMatrix &family,
                       const Rcpp::NumericMatrix &corr, const double info_frac,
                       const Rcpp::NumericMatrix &stg1_bnd,
                       const Rcpp::NumericMatrix &stg2_bnd) {
  ListPtr wgtmat_ptr = listcpp_from_rlist(wgtmat);
  WeightMatrix wgt_pair;
  wgt_pair.inthyp = wgtmat_ptr->get<IntMatrix>("inthyp");
  wgt_pair.wgtmat = wgtmat_ptr->get<FlatMatrix>("wgtmat");
  auto family1 = boolmatrix_from_Rmatrix(family);
  auto corr1 = flatmatrix_from_Rmatrix(corr);
  auto stg1_bnd1 = flatmatrix_from_Rmatrix(stg1_bnd);
  auto stg2_bnd1 = flatmatrix_from_Rmatrix(stg2_bnd);
  auto out = fCERCerCpp(stg1_p, wgt_pair, family1, corr1, info_frac, stg1_bnd1,
                        stg2_bnd1);
  // convert to 1-based indices
  for (size_t &i : out.stg1_elemhyp_r_idx)
    ++i;
  for (size_t &i : out.stg1_inthyp_nr_idx)
    ++i;
  ListCpp result;
  result.push_back(std::move(out.stg1_elemhyp_r_idx), "stg1_elemhyp_r_idx");
  result.push_back(std::move(out.stg1_inthyp_nr_idx), "stg1_inthyp_nr_idx");
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.CER), "CER");
  return Rcpp::wrap(result);
}

// Helper to compute adjusted p-values for Dunnett-based graphical approaches
AdjustedBoundaries
fCERNewBoundCpp(const std::vector<double> &stg1_p, const WeightMatrix &wgtmat,
                const BoolMatrix &family, const FlatMatrix &corr,
                const std::vector<size_t> &stg1_inthyp_nr_idx,
                const std::vector<double> &CER,
                const std::vector<size_t> &stg2_elemhyp_idx,
                const WeightMatrix &stg2_wgtmat, const double info_frac_new) {

  size_t m = wgtmat.inthyp.ncol;
  size_t nfams = family.nrow;
  size_t m2 = stg2_elemhyp_idx.size();
  size_t ntests2 = (1 << m2) - 1;
  size_t s = stg1_inthyp_nr_idx.size();

  if (info_frac_new <= 0.0 || info_frac_new >= 1.0) {
    throw std::invalid_argument("info_frac_new must be in (0, 1)");
  }

  std::vector<size_t> J2, J2_h;
  std::vector<std::vector<size_t>> J2_h_list(nfams);

  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> mean;
  FlatMatrix sigma;
  std::vector<std::vector<double>> lower_list(nfams);
  std::vector<FlatMatrix> sigma_list(nfams);

  double sqrttnew = std::sqrt(info_frac_new);

  IntMatrix inthyp2(s, m);
  std::vector<double> cJ2(s);
  FlatMatrix wcJ2(s, m);

  std::vector<size_t> cc(m);
  std::vector<size_t> cc2;
  for (size_t index = 0; index < s; ++index) {
    size_t i = stg1_inthyp_nr_idx[index];
    for (size_t j = 0; j < m; ++j) {
      cc[j] = wgtmat.inthyp(i, j);
      inthyp2(index, j) = cc[j];
    }

    bool withHypInI2 = false;
    for (size_t j : stg2_elemhyp_idx) {
      if (cc[j]) {
        withHypInI2 = true;
        break;
      }
    }

    if (withHypInI2) {
      cc2 = subset(cc, stg2_elemhyp_idx);

      size_t i2 = 0;
      for (size_t j2 = 0; j2 < m2; ++j2) {
        if (cc2[j2])
          i2 |= (1 << (m2 - 1 - j2));
      }
      i2 = ntests2 - i2;

      J2.clear();
      for (size_t j2 = 0; j2 < m2; ++j2) {
        if (stg2_wgtmat.wgtmat(i2, j2) > 0.0) {
          J2.push_back(j2);
        }
      }

      // Mehta et al (2025), equation (31): a conditional error rate above 1
      // is attained by any boundary, so no root finding is needed
      bool unbounded = CER[index] > 1.0;

      for (size_t h = 0; h < nfams; ++h) {
        J2_h.clear();
        for (size_t j2 : J2) {
          if (family(h, stg2_elemhyp_idx[j2]))
            J2_h.push_back(j2);
        }
        J2_h_list[h] = J2_h;

        size_t k = J2_h.size();
        if (k > 0 && !unbounded) {
          lower.resize(k);
          lower.assign(k, -POS_INF);
          lower_list[h] = lower;

          sigma.resize(k, k);
          for (size_t t1 = 0; t1 < k; ++t1) {
            for (size_t t2 = 0; t2 < k; ++t2) {
              size_t j1 = J2_h[t1];
              size_t j2 = J2_h[t2];
              sigma(t1, t2) = (1.0 - info_frac_new) *
                              corr(stg2_elemhyp_idx[j1], stg2_elemhyp_idx[j2]);
            }
          }
          sigma_list[h] = sigma;
        }
      }

      auto f2 = [&](double c) {
        double sum2 = 0.0;
        for (size_t h = 0; h < nfams; ++h) {
          J2_h = J2_h_list[h];
          size_t k = J2_h.size();
          if (k > 0) {
            const auto &lowerref = lower_list[h];
            const auto &sigmaref = sigma_list[h];

            mean.resize(k);
            upper.resize(k);
            for (size_t t = 0; t < k; ++t) {
              size_t j2 = J2_h[t];
              mean[t] =
                  sqrttnew * boost_qnorm(1.0 - stg1_p[stg2_elemhyp_idx[j2]]);
              upper[t] = boost_qnorm(1.0 - stg2_wgtmat.wgtmat(i2, j2) * c);
            }

            double v = pmvnormcpp(lowerref, upper, mean, sigmaref).prob;
            sum2 += (1.0 - v);
          }
        }
        return sum2 - CER[index];
      };

      cJ2[index] = unbounded ? POS_INF : brent(f2, 0.00001, 0.99999, 1e-5);

      for (size_t h = 0; h < nfams; ++h) {
        J2_h = J2_h_list[h];
        size_t k = J2_h.size();
        if (k > 0) {
          for (size_t t = 0; t < k; ++t) {
            size_t j2 = J2_h[t];
            wcJ2(index, stg2_elemhyp_idx[j2]) =
                stg2_wgtmat.wgtmat(i2, j2) * cJ2[index];
          }
        }
      }
    } else {
      cJ2[index] = 0.0;
    }
  }

  return AdjustedBoundaries{std::move(inthyp2), std::move(cJ2),
                            std::move(wcJ2)};
}

// [[Rcpp::export]]
Rcpp::List fCERNewBoundRcpp(
    const std::vector<double> &stg1_p, const Rcpp::List &wgtmat,
    const Rcpp::LogicalMatrix &family, const Rcpp::NumericMatrix &corr,
    const std::vector<int> &stg1_inthyp_nr_idx, const std::vector<double> &CER,
    const std::vector<int> &stg2_elemhyp_idx, const Rcpp::List &stg2_wgtmat,
    const double info_frac_new) {

  ListPtr wgtmat_ptr = listcpp_from_rlist(wgtmat);
  WeightMatrix wgt_pair;
  wgt_pair.inthyp = wgtmat_ptr->get<IntMatrix>("inthyp");
  wgt_pair.wgtmat = wgtmat_ptr->get<FlatMatrix>("wgtmat");

  int m = wgtmat_ptr->get<IntMatrix>("inthyp").ncol;
  int ntests = (1 << m) - 1;

  auto Jplus =
      validateConvertIdx(stg1_inthyp_nr_idx, ntests, "stg1_inthyp_nr_idx");
  auto I2 = validateConvertIdx(stg2_elemhyp_idx, m, "stg2_elemhyp_idx");

  ListPtr wgtmat2_ptr = listcpp_from_rlist(stg2_wgtmat);
  WeightMatrix wgt2_pair;
  wgt2_pair.inthyp = wgtmat2_ptr->get<IntMatrix>("inthyp");
  wgt2_pair.wgtmat = wgtmat2_ptr->get<FlatMatrix>("wgtmat");
  auto family1 = boolmatrix_from_Rmatrix(family);
  auto corr1 = flatmatrix_from_Rmatrix(corr);
  auto out = fCERNewBoundCpp(stg1_p, wgt_pair, family1, corr1, Jplus, CER, I2,
                             wgt2_pair, info_frac_new);
  ListCpp result;
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.stg2_coef_new), "stg2_coef_new");
  result.push_back(std::move(out.stg2_bnd_new), "stg2_bnd_new");
  return Rcpp::wrap(result);
}

std::vector<int> fCERRejCpp(const std::vector<double> &cum_p,
                            const std::vector<size_t> &stg1_elemhyp_r_idx,
                            const std::vector<size_t> &stg2_elemhyp_idx,
                            const IntMatrix &stg2_inthyp,
                            const FlatMatrix &stg2_bnd_new) {

  size_t s = stg2_inthyp.nrow;
  size_t m = stg2_inthyp.ncol;
  size_t m2 = stg2_elemhyp_idx.size();

  std::vector<int> rej_inter(s, 0);
  for (size_t index = 0; index < s; ++index) {
    for (size_t j2 = 0; j2 < m2; ++j2) {
      if (cum_p[j2] <= stg2_bnd_new(index, stg2_elemhyp_idx[j2])) {
        rej_inter[index] = 1;
        break;
      }
    }
  }

  std::vector<int> rej_elem(m);
  for (int i : stg1_elemhyp_r_idx)
    rej_elem[i] = 1;
  for (int i : stg2_elemhyp_idx)
    rej_elem[i] = 1;
  for (size_t index = 0; index < s; ++index) {
    if (!rej_inter[index]) {
      for (size_t j2 = 0; j2 < m2; ++j2) {
        if (stg2_inthyp(index, stg2_elemhyp_idx[j2])) {
          rej_elem[stg2_elemhyp_idx[j2]] = 0;
        }
      }
    }
  }

  return rej_elem;
}

// [[Rcpp::export]]
Rcpp::LogicalVector fCERRejRcpp(const std::vector<double> &cum_p,
                                const std::vector<int> &stg1_elemhyp_r_idx,
                                const std::vector<int> &stg2_elemhyp_idx,
                                const Rcpp::IntegerMatrix &stg2_inthyp,
                                const Rcpp::NumericMatrix &stg2_bnd_new) {

  int m = stg2_inthyp.ncol();
  auto I1r = validateConvertIdx(stg1_elemhyp_r_idx, m, "stg1_elemhyp_r_idx");
  auto I2 = validateConvertIdx(stg2_elemhyp_idx, m, "stg2_elemhyp_idx");
  auto inthyp2 = intmatrix_from_Rmatrix(stg2_inthyp);
  auto wcJ2 = flatmatrix_from_Rmatrix(stg2_bnd_new);

  auto out = fCERRejCpp(cum_p, I1r, I2, inthyp2, wcJ2);
  return Rcpp::wrap(out);
}
