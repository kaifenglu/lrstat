#include "generic_design.h"
#include "mams_design.h"
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
#include <utility>       // pair, make_pair

#include <Rcpp.h>

using std::size_t;


// Compute the p-value given theta, look L, observed z at look L (zL),
// number of active arms M, allocation ratio to common control r,
// whether the correlation is known,
// critical values vector b (length L), and information vector I (length L).
double f_pvalue_mams(const double theta,
                     const size_t M,
                     const double r,
                     const bool corr_known,
                     const size_t L,
                     const std::vector<double>& zL,
                     const FlatMatrix& b,
                     const std::vector<double>& I) {
  // Build the components required by exitprobcpp:
  // upper flatmatrix: first L-1 columns from b, last column = zL
  // theta vector: all = theta scalar
  FlatMatrix upper(M, L);
  if (L > 1) {
    const std::size_t bytes = (L - 1) * M * sizeof(double);
    std::memcpy(upper.data_ptr(), b.data_ptr(), bytes);
  }
  double* last_col = upper.data_ptr() + FlatMatrix::idx_col(0, L - 1, M);
  std::memcpy(last_col, zL.data(), M * sizeof(double));

  std::vector<double> mu(M, theta);
  auto probs = exitprob_mams_cpp(M, r, mu, corr_known, L, upper, I);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  double sum_up = std::accumulate(v.begin(), v.end(), 0.0);
  return sum_up;
}


// Helper to compute the confidence interval at the end of a group sequential trial
DataFrameCpp getCI_mams_cpp(
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t L,
    const std::vector<double>& zL,
    const double IMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const FlatMatrix& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& spendingTime) {

  // Basic argument checks
  if (M < 1) throw std::invalid_argument("M should be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r should be positive");
  if (L <= 0) throw std::invalid_argument("L must be a positive integer");
  if (!none_na(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (!none_na(informationRates))
    throw std::invalid_argument("informationRates must be provided");

  size_t kMax = informationRates.size();
  if (kMax < L)
    throw std::invalid_argument("Insufficient length for informationRates");
  if (informationRates[0] <= 0.0)
    throw std::invalid_argument("informationRates must be positive");
  if (any_nonincreasing(informationRates))
    throw std::invalid_argument("informationRates must be increasing");
  if (informationRates.back() > 1.0)
    throw std::invalid_argument("informationRates must not exceed 1");

  // efficacyStopping: if provided, validate, otherwise default to all ones
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() < kMax)
      throw std::invalid_argument("Insufficient length for efficacyStopping");
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping;
  } else {
    effStopping.assign(kMax, 1);
  }

  // spendingTime: if provided validate, otherwise use informationRates
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() < kMax)
      throw std::invalid_argument("Insufficient length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime.back() > 1.0)
      throw std::invalid_argument("spendingTime must not exceed 1");
    spendTime = spendingTime;
  } else {
    spendTime = informationRates;
  }

  // alpha checks
  if (std::isnan(alpha)) throw std::invalid_argument("alpha must be provided");
  if (alpha < 0.00001 || alpha >= 0.5)
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");

  // typeAlphaSpending to lower-case
  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (!none_na(criticalValues.data) && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument("parameterAlphaSpending must be positive for sfKD");
  }

  if (asf == "of" || asf == "p" || asf == "wt" || asf == "none") {
    if (informationRates.back() != 1.0) {
      throw std::invalid_argument(
          "informationRates must end with 1 for OF, P, WT, or NONE");
    }
    if (spendTime.back() != 1.0) {
      throw std::invalid_argument(
          "spendingTime must end with 1 for OF, P, WT, or NONE");
    }
  }


  // critical values: if not provided, compute using getBound_mams_cpp
  FlatMatrix bMat(L, M); // level M, M-1, ...., 1 critical values
  if (none_na(criticalValues.data)) {
    if (criticalValues.nrow < L)
      throw std::invalid_argument("Insufficient rows for criticalValues");
    if (criticalValues.ncol < M)
      throw std::invalid_argument("Insufficient columns for criticalValues");
    bMat = criticalValues;
  } else {
    for (size_t i = 0; i < M; ++i) {
      size_t level = M - i;
      auto v = getBound_mams_cpp(level, r, corr_known, L, informationRates,
                                 alpha, asf, parameterAlphaSpending,
                                 std::vector<double>{}, spendTime, effStopping);
      flatmatrix_set_column(bMat, i, v);
    }
  }


  // Build full information vector I = IMax * informationRates
  std::vector<double> I(L);
  for (size_t i = 0; i < L; ++i) I[i] = IMax * informationRates[i];

  // sort zL in descending order
  std::vector<size_t> order = seqcpp(0, M-1);
  std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
    return zL[i] > zL[j];
  });


  std::vector<size_t> level(M);
  std::vector<size_t> index(M);
  std::vector<double> pvalue(M);
  std::vector<double> thetahat(M);
  std::vector<double> lower(M);
  std::vector<double> upper(M);

  double sqrtIL = std::sqrt(I[L-1]);
  double cilevel = 1.0 - 2.0 * alpha;
  double target_lower = (1.0 - cilevel) / 2.0;
  double target_upper = (1.0 + cilevel) / 2.0;
  double tol = 1.0e-6;

  // level-M test boundary
  const double* b1col_ptr = bMat.data_ptr(); // column 0
  FlatMatrix b1_matrix(1, L);
  for (size_t j = 0; j < L; ++j) {
    b1_matrix(0, j) = b1col_ptr[j];
  }

  std::vector<double> zL_vec;           // reused for different levels
  std::vector<double> zL1_vec(1);       // reuse single-element vector
  for (size_t h = 0; h < M; ++h) {
    level[h] = M - h;
    index[h] = order[h] + 1; // 1-based index for R
    double zLmax = zL[order[h]];

    const double* bcol_ptr = bMat.data_ptr() + h * bMat.nrow;
    FlatMatrix b_matrix(level[h], L);
    for (size_t j = 0; j < L; ++j) {
      double* colptr = b_matrix.data_ptr() + j * b_matrix.nrow;
      std::fill_n(colptr, level[h], bcol_ptr[j]);
    }

    // p-value at theta = 0
    zL_vec.assign(level[h], zLmax);
    pvalue[h] = f_pvalue_mams(0.0, level[h], r, corr_known, L, zL_vec,
                              b_matrix, I);

    double left = (zLmax - 8.0) / sqrtIL;
    double right = (zLmax + 8.0) / sqrtIL;

    // median estimate thetahat: solve f_pvalue(theta) - 0.5 = 0
    auto f_med = [&](double theta)->double {
      return f_pvalue_mams(theta, level[h], r, corr_known, L, zL_vec,
                           b_matrix, I) - 0.5;
    };
    thetahat[h] = brent(f_med, left, right, tol);

    // lower bound: solve f_pvalue(theta) - (1 - cilevel)/2 = 0
    auto f_lower = [&](double theta)->double {
      return f_pvalue_mams(theta, level[h], r, corr_known, L, zL_vec,
                           b_matrix, I) - target_lower;
    };
    lower[h] = brent(f_lower, left, thetahat[h], tol);

    // upper bound: solve f_pvalue(theta) - (1 + cilevel)/2 = 0
    zL1_vec[0] = zLmax;
    auto f_upper = [&](double theta)->double {
      return f_pvalue_mams(theta, 1, r, corr_known, L, zL1_vec,
                           b1_matrix, I) - target_upper;
    };
    upper[h] = brent(f_upper, thetahat[h], right, tol);
  }

  // Build DataFrameCpp result
  DataFrameCpp df;
  df.push_back(level, "level");
  df.push_back(index, "index");
  df.push_back(pvalue, "pvalue");
  df.push_back(thetahat, "thetahat");
  df.push_back(cilevel, "cilevel");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");

  return df;
}



//' @title Confidence Interval After Trial Termination for a Multi-Arm
//' Multi-Stage Design
//' @description Obtains the p-value, conservative point estimate, and
//' confidence interval after the end of a multi-arm multi-stage trial.
//'
//' @param M Number of active treatment arms.
//' @param r Randomization ratio of each active arm to the common control.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param L The termination look.
//' @param zL The vector of z-test statistics at the termination look.
//' @param IMax Maximum information for any active arm versus the common
//'   control.
//' @param informationRates The information rates up to look \code{L}.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage up to look \code{L}.
//'   Defaults to \code{TRUE} if left unspecified.
//' @param criticalValues The matrix of by-level upper boundaries on the
//'   max z-test statistic scale for efficacy stopping up to look \code{L}.
//'   The first column is for level \code{M}, the second column is for
//'   level \code{M - 1}, and so on, with the last column for level 1.
//'   If left unspecified, the critical values will be computed based
//'   on the specified alpha spending function.
//' @inheritParams param_alpha
//' @param typeAlphaSpending The type of alpha spending for the trial.
//'   One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @inheritParams param_parameterAlphaSpending
//' @param spendingTime The error spending time up to look \code{L}.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//'
//' @details
//' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
//' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
//' and \code{spendingTime} must be of full length \code{kMax}, and
//' \code{informationRates} and \code{spendingTime} must end with 1.
//'
//' @return A data frame with the following components:
//'
//' * \code{level}: Number of individual hypotheses considered for multiplicity.
//'
//' * \code{index}: The treatment arm with max Z among the active arms.
//'
//' * \code{pvalue}: p-value for rejecting the null hypothesis.
//'
//' * \code{thetahat}: Point estimate of the parameter.
//'
//' * \code{cilevel}: Confidence interval level.
//'
//' * \code{lower}: Lower bound of confidence interval.
//'
//' * \code{upper}: Upper bound of confidence interval.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @examples
//' getCI_mams(
//'   L = 2, zL = c(2.075, 2.264),
//'   M = 2, r = 1, corr_known = FALSE,
//'   IMax = 300 / 4, informationRates = c(1/2, 1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF")
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame getCI_mams(
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int L = NA_INTEGER,
    const Rcpp::NumericVector& zL = NA_REAL,
    const double IMax = NA_REAL,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::Nullable<Rcpp::NumericMatrix> criticalValues = R_NilValue,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {

  std::vector<double> zLVec(zL.begin(), zL.end());
  std::vector<double> infoRates(informationRates.begin(), informationRates.end());
  auto effStopping = convertLogicalVector(efficacyStopping);
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());

  // Handle optional matrix safely
  FlatMatrix critValues;
  if (criticalValues.isNotNull()) {
    Rcpp::NumericMatrix cm(criticalValues); // unwrap
    critValues = flatmatrix_from_Rmatrix(cm);
  } else {
    critValues = FlatMatrix(1, 1);
    critValues(0, 0) = std::numeric_limits<double>::quiet_NaN(); // placeholder
  }

  auto result = getCI_mams_cpp(
    static_cast<size_t>(M), r, corr_known, static_cast<size_t>(L), zLVec,
    IMax, infoRates, effStopping, critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, spendTime);
  return Rcpp::wrap(result);
}


// Compute the backward image (J, zJ)
std::pair<size_t, double> f_bwimage_mams(const double theta,
                                         const size_t M,
                                         const double r,
                                         const bool corr_known,
                                         const size_t kMax,
                                         const size_t L,
                                         const std::vector<double>& zL,
                                         const std::vector<double>& b,
                                         const std::vector<double>& I,
                                         const size_t M2,
                                         const double r2,
                                         const size_t L2,
                                         const std::vector<double>& zL2,
                                         const FlatMatrix& b2,
                                         const std::vector<double>& I2) {

  // compute astar for the adapted secondary trial
  double astar = f_pvalue_mams(theta, M2, r2, corr_known, L2, zL2, b2, I2);

  // prepare b1, mu, I1 for the original secondary trial
  size_t k1 = kMax - L;
  std::vector<double> I1(k1);
  for (size_t l = 0; l < k1; ++l) {
    I1[l] = I[l + L] - I[L - 1];
  }

  FlatMatrix b1(M, k1);
  for (size_t i = 0; i < k1; ++i) {
    double r1 = I[L - 1] / I[L + i];
    double sqrt_r1 = std::sqrt(r1);
    double denom = std::sqrt(1.0 - r1);
    double cut = b[L + i];
    double* colptr = b1.data_ptr() + i * M; // contiguous column start
    for (size_t m = 0; m < M; ++m) {
      colptr[m] = (cut - zL[m] * sqrt_r1) / denom;
    }
  }

  // compute exit probabilities for b1
  std::vector<double> mu(M, theta);
  auto probs = exitprob_mams_cpp(M, r, mu, corr_known, k1, b1, I1);
  auto pu = probs.get<std::vector<double>>("exitProbUpper");

  // find interval containing astar
  std::vector<double> cpu(k1);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  size_t j = std::min(findInterval1(astar, cpu) + 1, k1);
  size_t J = L + j; // combined stage index in primary trial numbering

  // find zJ
  double r1 = I[L - 1] / I[L + j - 1];
  std::vector<double> zj(M);
  auto f = [&](double z)->double {
    double sqrt_r1 = std::sqrt(r1);
    double denom = std::sqrt(1.0 - r1);
    for (size_t m = 0; m < M; ++m) {
      zj[m] = (z - zL[m] * sqrt_r1) / denom;
    }
    return f_pvalue_mams(theta, M, r, corr_known, j, zj, b1, I1) - astar;
  };

  double zJ;
  if (j < k1) {
    zJ = brent(f, b[L + j - 1], 8.0, 1e-6);
  } else {
    double zLmin = *std::min_element(zL.begin(), zL.end());
    double lo = -8.0 * std::sqrt(1.0 - r1) + zLmin * std::sqrt(r1);
    zJ = brent(f, lo, 8.0, 1e-6);
  }

  return std::make_pair(J, zJ);
}


// compute backward p-value for adapted trial
double f_bwpvalue_mams(const double theta,
                       const size_t M,
                       const double r,
                       const bool corr_known,
                       const size_t kMax,
                       const size_t L,
                       const std::vector<double>& zL,
                       const std::vector<double>& b,
                       const std::vector<double>& I,
                       const size_t M2,
                       const double r2,
                       const size_t L2,
                       const std::vector<double>& zL2,
                       const FlatMatrix& b2,
                       const std::vector<double>& I2) {
  auto bw = f_bwimage_mams(theta, M, r, corr_known, kMax, L, zL, b, I,
                           M2, r2, L2, zL2, b2, I2);

  size_t J = bw.first;
  std::vector<double> zJ(M, bw.second);
  FlatMatrix bMat(M, kMax);
  for (size_t i = 0; i < kMax; ++i) {
    double* colptr = bMat.data_ptr() + i * M;
    std::fill_n(colptr, M, b[i]);
  }

  return f_pvalue_mams(theta, M, r, corr_known, J, zJ, bMat, I);
}


// Helper to compute confidence interval after the end of an adaptive trial
DataFrameCpp getADCI_mams_cpp(
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t L,
    const std::vector<double>& zL,
    const double IMax,
    const size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const FlatMatrix& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
    const size_t MNew,
    const std::vector<int>& selected,
    const double rNew,
    const size_t Lc,
    const std::vector<double>& zLc,
    const double INew,
    const std::vector<double>& informationRatesNew,
    const std::vector<unsigned char>& efficacyStoppingNew,
    const std::string& typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double>& spendingTimeNew) {

  // Input validation and defaults
  if (M < 1) throw std::invalid_argument("M must be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r must be positive");
  if (L <= 0) throw std::invalid_argument("L must be provided and positive");
  if (!none_na(zL)) throw std::invalid_argument("zL must be provided");
  if (zL.size() != M) throw std::invalid_argument("Invalid length for zL");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

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
    infoRates = informationRates; // copy if provided
  } else {
    infoRates.resize(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      infoRates[i] = static_cast<double>(i + 1) / static_cast<double>(kMax);
    }
  }

  // efficacyStopping: default to all ones if not provided
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy if provided
  } else {
    effStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues.data);
  if (!missingCriticalValues && (criticalValues.ncol != M ||
      criticalValues.nrow != kMax)) {
    throw std::invalid_argument("Invalid dimension for criticalValues");
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 0.5)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 0.5)");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
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


  // Now handle new trial inputs
  std::vector<size_t> selectedNew;
  std::vector<double> infoRatesNew;
  std::vector<unsigned char> effStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (MNew < 1) throw std::invalid_argument("MNew must be at least 1");
  for (auto i : selected) {
    if (i < 1 || i > static_cast<int>(M)) {
      throw std::invalid_argument("Invalid value in selected");
    }
    selectedNew.push_back(static_cast<size_t>(i - 1));
  }
  selectedNew = unique_sorted(selectedNew);
  if (selectedNew.size() != MNew)
    throw std::invalid_argument("Length of selected does not match MNew");

  if (rNew <= 0.0) throw std::invalid_argument("rNew must be positive");

  if (Lc <= L) throw std::invalid_argument("Lc must be greater than L");
  if (!none_na(zLc)) throw std::invalid_argument("zLc must be provided");
  if (std::isnan(INew)) throw std::invalid_argument("INew must be provided");
  if (INew <= 0.0) throw std::invalid_argument("INew must be positive");

  size_t L2 = Lc - L;

  size_t kNew = L2;
  if (MullerSchafer) {
    if (none_na(informationRatesNew)) {
      kNew = informationRatesNew.size();
      if (kNew < L2)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew.back() > 1.0)
        throw std::invalid_argument("informationRatesNew must not exceed 1");
    } else {
      throw std::invalid_argument(
          "informationRatesNew must be provided for MullerSchafer");
    }

    if (none_na(efficacyStoppingNew)) {
      if (efficacyStoppingNew.size() < kNew)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew.back() != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew;
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    if (!(asfNew == "of" || asfNew == "sfof" || asfNew == "sfp" ||
        asfNew == "sfkd" || asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }

    if ((asfNew == "sfkd" || asfNew == "sfhsd") &&
        std::isnan(parameterAlphaSpendingNew))
      throw std::invalid_argument("Missing value for parameterAlphaSpendingNew");

    if (asfNew == "sfkd" && parameterAlphaSpendingNew <= 0.0)
      throw std::invalid_argument(
          "parameterAlphaSpendingNew must be positive for sfKD");

    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() < kNew)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew.back() > 1.0)
        throw std::invalid_argument("spendingTimeNew must not exceed 1");
      spendTimeNew = spendingTimeNew;
    } else {
      spendTimeNew = informationRatesNew;
    }

    if (asfNew == "of" || asfNew == "none") {
      if (informationRatesNew.back() != 1.0) {
        throw std::invalid_argument(
            "informationRatesNew must end with 1 for OF or NONE");
      }
      if (spendTimeNew.back() != 1.0) {
        throw std::invalid_argument(
            "spendingTimeNew must end with 1 for OF or NONE");
      }
    }
  }


  // obtain by-level critical values for the primary trial
  FlatMatrix efficacyBounds1(kMax, M);
  if (!none_na(criticalValues.data)) {
    for (size_t M1 = M; M1 > 0; --M1) {
      auto cut = getBound_mams_cpp(
        M1, r, corr_known, kMax, infoRates, alpha, asf, parameterAlphaSpending,
        std::vector<double>{}, spendTime, effStopping);
      flatmatrix_set_column(efficacyBounds1, M - M1, cut);
    }

    for (size_t m = 0; m < M; ++m) {
      for (size_t i = 0; i < kMax; ++i) {
        if (!effStopping[i]) efficacyBounds1(i, m) = 8.0;
      }
    }
  } else {
    efficacyBounds1 = criticalValues; // copy if provided
  }

  // Primary information vector
  std::vector<double> I(kMax);
  for (size_t i = 0; i < kMax; ++i) I[i] = IMax * informationRates[i];


  // sort zLc in descending order
  std::vector<size_t> order = seqcpp(0, MNew - 1);
  std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
    return zLc[i] > zLc[j];
  });


  size_t k1 = kMax - L;
  std::vector<double> s1(k1);
  for (size_t i = 0; i < k1; ++i) {
    s1[i] = (infoRates[L + i] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
  }

  // secondary trial information for original design
  double INew1 = IMax * (1.0 - infoRates[L - 1]);
  std::vector<double> I1(k1);
  for (size_t i = 0; i < k1; ++i) {
    I1[i] = INew1 * s1[i];
  }

  double IL = IMax * infoRates[L - 1];
  double sqrtIL = std::sqrt(IL);


  size_t k2 = MullerSchafer ? kNew : k1;
  std::vector<double> s2(k2);
  std::string asf2;
  std::vector<unsigned char> effStopping2(k2);

  if (!MullerSchafer) {
    // number of secondary trial looks is the same as the original trial
    // spend the conditional type I error as in the original trial
    s2 = s1;
    asf2 = "user";
    std::memcpy(effStopping2.data(), effStopping.data() + L, k2);
  } else {
    s2 = infoRatesNew;
    asf2 = asfNew;
    effStopping2 = effStoppingNew;
  }

  std::vector<double> I2(k2); // information levels for secondary trial
  std::vector<double> Ic(k2); // information levels for integrated trial
  std::vector<double> sqrtI2(k2), sqrtIc(k2);
  for (size_t i = 0; i < k2; ++i) {
    I2[i] = INew * s2[i];
    Ic[i] = I2[i] + IL;
    sqrtI2[i] = std::sqrt(I2[i]);
    sqrtIc[i] = std::sqrt(Ic[i]);
  }

  std::vector<size_t> level(MNew);
  std::vector<size_t> index(MNew);
  std::vector<double> pvalue(MNew);
  std::vector<double> thetahat(MNew);
  std::vector<double> lower(MNew);
  std::vector<double> upper(MNew);

  double zLmin = *std::min_element(zL.begin(), zL.end());
  double zLmax = *std::max_element(zL.begin(), zL.end());
  double left = (zLmin - efficacyBounds1(L - 1, 0)) / sqrtIL;
  double right = (zLmax + efficacyBounds1(L - 1, 0)) / sqrtIL;
  double tol = 1.0e-6;

  double cilevel = 1.0 - 2.0 * alpha;
  double target_lower = (1.0 - cilevel) / 2.0;
  double target_upper = (1.0 + cilevel) / 2.0;

  // build map from arm index -> position in selectedNew (or -1)
  std::vector<int> idx_in_selected(M, -1);
  for (size_t t = 0; t < selectedNew.size(); ++t)
    idx_in_selected[selectedNew[t]] = static_cast<int>(t);

  // reusable scratch buffers
  std::vector<double> zL1; zL1.reserve(M);
  std::vector<double> zL2; zL2.reserve(MNew);
  std::vector<double> zL_1(1), zL2_1(1);
  FlatMatrix c2_1(1, k2);
  for (size_t h = 0; h < MNew; ++h) {
    level[h] = MNew - h;
    index[h] = selectedNew[order[h]] + 1; // 1-based index in the original trial

    // obtain the boundaries after excluding theta values that ranked higher
    std::vector<size_t> primary;
    std::vector<size_t> selectedNew2;
    primary.reserve(M);
    selectedNew2.reserve(MNew);
    for (size_t j = 0; j < M; ++j) {
      int pos = idx_in_selected[j];
      if (pos >= 0) {
        // pos refers to index in selectedNew; check whether its order is >= h
        for (size_t i = h; i < MNew; ++i) {
          if (pos == static_cast<int>(order[i])) {
            selectedNew2.push_back(j);
            primary.push_back(j);
            break;
          }
        }
      } else {
        primary.push_back(j);
      }
    }

    size_t M1 = primary.size();
    size_t M2 = selectedNew2.size();

    auto critValues = flatmatrix_get_column(efficacyBounds1, M - M1);

    // compute conditional alpha for MAMS with the given M1 hypotheses
    FlatMatrix c1(M1, k1);
    for (size_t i = 0; i < k1; ++i) {
      double* colptr = c1.data_ptr() + i * M1;
      if (effStopping[L + i]) {
        double crit = critValues[L + i];
        double r1 = infoRates[L - 1] / infoRates[L + i];
        double sqrt_r1 = std::sqrt(r1);
        double denom = std::sqrt(1.0 - r1);
        for (size_t j = 0; j < M1; ++j) {
          size_t m = primary[j];
          colptr[j] = (crit - zL[m] * sqrt_r1) / denom;
        }
      } else {
        std::fill_n(colptr, M1, 8.0);
      }
    }

    // conditional type I error
    std::vector<double> zero1(M1, 0.0);
    auto probs = exitprob_mams_cpp(M1, r, zero1, corr_known, k1, c1, I1);
    auto v0 = probs.get<std::vector<double>>("exitProbUpper");
    double c_alpha = std::accumulate(v0.begin(), v0.end(), 0.0);

    std::vector<double> cpu0(k2);
    if (!MullerSchafer) {
      std::partial_sum(v0.begin(), v0.end(), cpu0.begin());
    } else {
      if (asf2 != "none" && asf2 != "of") {
        for (size_t i = 0; i < k2; ++i) {
          cpu0[i] = errorSpentcpp(spendTimeNew[i], c_alpha,
                                  asfNew, parameterAlphaSpendingNew);
        }
      }
    }

    std::vector<double> zero2(M2, 0.0);
    FlatMatrix c2(M2, k2);
    c2.fill(8.0);

    std::vector<double> zscaled(M2);
    for (size_t j = 0; j < M2; ++j) zscaled[j] = zL[selectedNew2[j]] * sqrtIL;

    if (asf2 == "of") {
      auto g = [&c2, &I2, &sqrtI2, &sqrtIc, &zscaled, &zero2, &effStopping2,
                k2, c_alpha, M2, rNew, corr_known]
      (double x)->double {
        double col_const = x * sqrtIc[k2 - 1];
        for (size_t i = 0; i < k2; ++i) {
          double* colptr = c2.data_ptr() + i * M2;
          if (effStopping2[i]) {
            double denom = sqrtI2[i];
            for (size_t j = 0; j < M2; ++j) {
              colptr[j] = (col_const - zscaled[j] ) / denom;
            }
          } else {
            std::fill_n(colptr, M2, 8.0);
          }
        }

        auto probs = exitprob_mams_cpp(M2, rNew, zero2, corr_known, k2, c2, I2);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - c_alpha;
      };

      double cof = brent(g, 0.0, 8.0, 1e-6);
      double col_const = cof * sqrtIc[k2 - 1];
      for (size_t i = 0; i < k2; ++i) {
        double* colptr = c2.data_ptr() + i * M2;
        if (effStopping2[i]) {
          double denom = sqrtI2[i];
          for (size_t j = 0; j < M2; ++j) {
            colptr[j] = (col_const - zscaled[j] ) / denom;
          }
        } else {
          std::fill_n(colptr, M2, 8.0);
        }
      }
    } else if (asf2 == "none") {
      double denom = sqrtI2[k2 - 1];
      auto g = [&c2, &I2, &sqrtIc, &zscaled, &zero2,
                denom, k2, c_alpha, M2, rNew, corr_known]
      (double x)->double {
        double col_const = x * sqrtIc[k2 - 1];
        double* colptr = c2.data_ptr() + (k2 - 1) * M2;
        for (size_t j = 0; j < M2; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }

        auto probs = exitprob_mams_cpp(M2, rNew, zero2, corr_known, k2, c2, I2);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - c_alpha;
      };

      double cof = brent(g, 0.0, 8.0, 1e-6);
      double col_const = cof * sqrtIc[k2 - 1];
      double* colptr = c2.data_ptr() + (k2 - 1) * M2;
      for (size_t j = 0; j < M2; ++j) {
        colptr[j] = (col_const - zscaled[j]) / denom;
      }
    } else {
      for (size_t i = 0; i < k2; ++i) {
        if (!effStopping2[i]) continue;
        double denom = sqrtI2[i];

        auto g = [&c2, &I2, &sqrtIc, &zscaled, &cpu0, &zero2,
                  denom, i, M2, rNew, corr_known]
        (double x)->double {
          double col_const = x * sqrtIc[i];
          double* colptr = c2.data_ptr() + i * M2;
          // update critical values of the secondary trial at current look
          for (size_t j = 0; j < M2; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }

          auto probs = exitprob_mams_cpp(M2, rNew, zero2, corr_known, i + 1, c2, I2);
          auto v = probs.get<std::vector<double>>("exitProbUpper");
          double p0 = std::accumulate(v.begin(), v.end(), 0.0);
          return p0 - cpu0[i];
        };

        double cof = brent(g, 0.0, 8.0, 1e-6);
        double col_const = cof * sqrtIc[i];
        double* colptr = c2.data_ptr() + i * M2;
        for (size_t j = 0; j < M2; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }
      }
    }


    zL1.resize(M1);
    for (size_t i = 0; i < M1; ++i) {
      zL1[i] = zL[primary[i]];
    }

    // construct the Z-statistic vector for the secondary trial
    zL2.resize(M2);
    double col_const = zLc[order[h]] * sqrtIc[L2 - 1];
    double denom = sqrtI2[L2 - 1];
    for (size_t i = 0; i < M2; ++i) {
      zL2[i] = (col_const - zL[selectedNew2[i]] * sqrtIL) / denom;
    }

    pvalue[h] = f_bwpvalue_mams(0.0, M1, r, corr_known, kMax, L, zL1, critValues,
                                I, M2, rNew, L2, zL2, c2, I2);

    auto f_med = [&](double theta)->double {
      return f_bwpvalue_mams(theta, M1, r, corr_known, kMax, L, zL1, critValues,
                              I, M2, rNew, L2, zL2, c2, I2) - 0.5;
    };
    thetahat[h] = brent(f_med, left, right, tol);

    auto f_low = [&](double theta)->double {
      return f_bwpvalue_mams(theta, M1, r, corr_known, kMax, L, zL1, critValues,
                              I, M2, rNew, L2, zL2, c2, I2) - target_lower;
    };
    lower[h] = brent(f_low, left, thetahat[h], tol);

    // find the index in {0, ..., M2-1} corresponding to the current hypothesis
    size_t j0 = 0;
    for (size_t j = 0; j < M2; ++j) {
      if (selectedNew2[j] == index[h] - 1) {
        j0 = j;
        break;
      }
    }

    zL_1[0] = zL[index[h] - 1];
    zL2_1[0] = zL2[j0];
    for (size_t i = 0; i < k2; ++i) {
      c2_1(0, i) = c2(j0, i);
    }

    auto f_high = [&](double theta)->double {
      return f_bwpvalue_mams(theta, 1, r, corr_known, kMax, L, zL_1, critValues,
                              I, 1, rNew, L2, zL2_1, c2_1, I2) - target_upper;
    };
    upper[h] = brent(f_high, thetahat[h], right, tol);
  }

  // Build DataFrameCpp result
  DataFrameCpp df;
  df.push_back(level, "level");
  df.push_back(index, "index");
  df.push_back(pvalue, "pvalue");
  df.push_back(thetahat, "thetahat");
  df.push_back(cilevel, "cilevel");
  df.push_back(lower, "lower");
  df.push_back(upper, "upper");

  return df;
}


//' @title Confidence Interval After Adaptation for a Multi-Arm Multi-Stage
//' Design
//' @description Obtains the p-value, conservative point estimate, and
//' confidence interval after the end of an adaptive multi-arm multi-stage trial.
//'
//' @param M Number of active treatment arms in the primary trial.
//' @param r Randomization ratio of each active arm to the common control
//'   in the primary trial.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is assumed.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistics at the interim adaptation look of
//'   the primary trial.
//' @param IMax Maximum information for any active arm versus the common
//'   control for the primary trial. Must be provided.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param criticalValues The matrix of by-level upper boundaries on the
//'   max z-test statistic scale for efficacy stopping up to look \code{L}
//'   for the primary trial.
//'   The first column is for level \code{M}, the second column is for
//'   level \code{M - 1}, and so on, with the last column for level 1.
//'   If left unspecified, the critical values will be computed based
//'   on the specified alpha spending function.
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
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param MNew The number of active treatment arms in the secondary trial.
//' @param selected The indices of the selected treatment arms for the
//'   secondary trial among the \code{M} active arms in the primary trial.
//' @param rNew The randomization ratio of each active arm to the common control
//'  in the secondary trial.
//' @param Lc The termination look of the integrated trial.
//' @param zLc The z-test statistics at the termination look of the
//'   integrated trial.
//' @param INew The maximum information for any active arm versus the common
//'   control in the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial.
//'   Defaults to \code{TRUE} if left unspecified.
//' @param typeAlphaSpendingNew The type of alpha spending for the secondary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case, it is
//'   the same as \code{informationRatesNew}.
//'
//' @details
//' If typeAlphaSpendingNew is \code{"OF"} or \code{"none"}, then
//' \code{informationRatesNew}, \code{efficacyStoppingNew}, and
//' \code{spendingTimeNew} must be of full length \code{kNew}, and
//' \code{informationRatesNew} and \code{spendingTimeNew} must end with 1.
//'
//' @return A data frame with the following variables:
//'
//' * \code{level}: Number of individual hypotheses considered for multiplicity.
//'
//' * \code{index}: The treatment arm with max Z among the active arms.
//'
//' * \code{pvalue}: p-value for rejecting the null hypothesis.
//'
//' * \code{thetahat}: Point estimate of the parameter.
//'
//' * \code{cilevel}: Confidence interval level.
//'
//' * \code{lower}: Lower bound of confidence interval.
//'
//' * \code{upper}: Upper bound of confidence interval.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @examples
//' getADCI_mams(
//'   M = 2, r = 1, corr_known = FALSE, L = 1, zL = c(2.075, 2.264),
//'   IMax = 300 / 4, kMax = 2, informationRates = c(0.5, 1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   MNew = 1, selected = 2, rNew = 1,
//'   Lc = 2, zLc = 1.667, INew = 374 / 4)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame getADCI_mams(
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int L = NA_INTEGER,
    const Rcpp::NumericVector& zL = NA_REAL,
    const double IMax = NA_REAL,
    const int kMax = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::Nullable<Rcpp::NumericMatrix> criticalValues = R_NilValue,
    const double alpha = 0.25,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int MNew = NA_INTEGER,
    const Rcpp::IntegerVector& selected = NA_INTEGER,
    const double rNew = 1,
    const int Lc = NA_INTEGER,
    const Rcpp::NumericVector& zLc = NA_REAL,
    const double INew = NA_REAL,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL) {

  std::vector<double> zLVec(zL.begin(), zL.end());
  std::vector<double> infoRates(informationRates.begin(), informationRates.end());
  auto effStopping = convertLogicalVector(efficacyStopping);

  // Handle optional matrix safely
  FlatMatrix critValues;
  if (criticalValues.isNotNull()) {
    Rcpp::NumericMatrix cm(criticalValues); // unwrap
    critValues = flatmatrix_from_Rmatrix(cm);
  } else {
    critValues = FlatMatrix(1, 1);
    critValues(0, 0) = std::numeric_limits<double>::quiet_NaN(); // placeholder
  }

  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());
  std::vector<int> selectedNew(selected.begin(), selected.end());
  std::vector<double> zLcVec(zLc.begin(), zLc.end());

  std::vector<double> infoRatesNew(informationRatesNew.begin(),
                                   informationRatesNew.end());
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  std::vector<double> spendTimeNew(spendingTimeNew.begin(), spendingTimeNew.end());

  auto result = getADCI_mams_cpp(
    static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(L), zLVec, IMax, static_cast<size_t>(kMax), infoRates,
    effStopping, critValues, alpha, typeAlphaSpending, parameterAlphaSpending,
    spendTime, MullerSchafer, static_cast<size_t>(MNew), selectedNew, rNew,
    static_cast<size_t>(Lc), zLcVec, INew, infoRatesNew, effStoppingNew,
    typeAlphaSpendingNew, parameterAlphaSpendingNew, spendTimeNew);
  return Rcpp::wrap(result);
}
