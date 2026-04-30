#include "dataframe_list.h"
#include "thread_utils.h"
#include "utilities.h"
#include "generic_design.h"
#include "mvnormr.h"

#include <Rcpp.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using std::size_t;


std::vector<double> exitprob_mams_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t kMax,
    const FlatMatrix& b,
    const std::vector<double>& I) {

  if (M < 1) throw std::invalid_argument("M should be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r should be positive");
  if (theta.size() != M) throw std::invalid_argument("theta should have length M");
  if (kMax < 1) throw std::invalid_argument("kMax should be at least 1");
  if (b.nrow != M) throw std::invalid_argument("b should have M rows");
  if (b.ncol < kMax) throw std::invalid_argument("Insufficient columns for b");

  std::vector<double> Ivec; Ivec.reserve(kMax);
  if (none_na(I)) {
    if (I.size() < kMax) throw std::invalid_argument("Insufficient length for I");
    Ivec.assign(I.begin(), I.begin() + kMax);
    if (Ivec[0] <= 0.0) throw std::invalid_argument("I must be positive");
    if (any_nonincreasing(Ivec)) throw std::invalid_argument("I must be increasing");
  } else {
    Ivec.resize(kMax);
    std::iota(Ivec.begin(), Ivec.end(), 1.0);
  }

  double rho = corr_known ? r / (r + 1.0) : 0;

  std::vector<double> sqrtI(kMax);
  for (size_t k = 0; k < kMax; ++k) {
    sqrtI[k] = std::sqrt(Ivec[k]);
  }

  std::vector<double> pcumreject(kMax, 0.0);

  if (rho == 0) {
    std::vector<double> upper(kMax);
    std::vector<double> lower(kMax, -6.0);
    std::vector<double> zero(kMax, 0.0);
    std::vector<double> p(kMax, 1.0);
    std::vector<double> cumexit(kMax, 0.0);
    for (size_t m = 0; m < M; ++m) {
      for (size_t k = 0; k < kMax; ++k) {
        upper[k] = b(m, k) - theta[m] * sqrtI[k];
        if (upper[k] < -6.0) upper[k] = -6.0;
      }
      auto probs = exitprobcpp(upper, lower, zero, Ivec);
      auto v = probs.get<std::vector<double>>("exitProbUpper");
      std::partial_sum(v.begin(), v.end(), cumexit.begin());
      for (size_t k = 0; k < kMax; ++k) {
        p[k] *= (1.0 - cumexit[k]);
      }
    }
    for (size_t k = 0; k < kMax; ++k) {
      pcumreject[k] = 1.0 - p[k];
    }
  } else if (M == 1) {
    std::vector<double> upper(kMax);
    std::vector<double> lower(kMax, -6.0);
    std::vector<double> zero(kMax, 0.0);
    for (size_t k = 0; k < kMax; ++k) {
      upper[k] = b(0, k) - theta[0] * sqrtI[k];
      if (upper[k] < -6.0) upper[k] = -6.0;
    }
    auto probs = exitprobcpp(upper, lower, zero, Ivec);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    std::partial_sum(v.begin(), v.end(), pcumreject.begin());
  } else if (kMax == 1) {
    double sqrtrho = std::sqrt(rho);
    double sqrt1mr = std::sqrt(1.0 - rho);
    double sqrtI0 = sqrtI[0];
    auto f = [&](double z0) {
      double value = 1.0;
      for (size_t m = 0; m < M; ++m) {
        double upper = (b(m, 0) - theta[m] * sqrtI0 - sqrtrho * z0) / sqrt1mr;
        double mass = boost_pnorm(upper);
        if (mass <= 0.0) return 0.0;
        value *= mass;
      }
      return value * boost_dnorm(z0);
    };

    std::vector<double> breaks = { -6.0, 0.0, 6.0 };
    double p = integrate3(f, breaks, 1e-6);
    pcumreject[0] = 1.0 - p;
  } else {
    double T = Ivec.back();
    double sqrtT = sqrtI.back();
    std::vector<double> s(kMax);
    for (size_t k = 0; k < kMax; ++k) {
      s[k] = Ivec[k] / T;
    }

    std::size_t fullN = M * kMax;
    std::vector<double> upper(fullN);
    for (size_t k = 0; k < kMax; ++k) {
      for (size_t m = 0; m < M; ++m) {
        size_t idx = k * M + m;
        upper[idx] = b(m, k) * sqrt(s[k]) - theta[m] * s[k] * sqrtT;
      }
    }

    std::vector<double> lower(fullN, -6.0);
    std::vector<double> zero(fullN, 0.0);

    // covariance matrix for scaled Brownian motion
    FlatMatrix sigma(fullN, fullN);

    for (size_t k2 = 0; k2 < kMax; ++k2) {
      for (size_t k1 = 0; k1 < kMax; ++k1) {
        size_t k = std::min(k1, k2);
        double val_off = rho * s[k];    // value for off-diagonal (m1 != m2)
        double val_diag = s[k];         // diagonal (m1 == m2)
        // block row base and block col base:
        std::size_t row_base = k1 * M;
        std::size_t col_base = k2 * M;
        // for each column (m2) inside the MxM block, fill contiguous column segment:
        for (size_t m2 = 0; m2 < M; ++m2) {
          // pointer to start of the column (col_index) at row row_base
          double* dst_col = sigma.data_ptr() + (col_base + m2) * fullN + row_base;
          // fill M entries in that column with off-diagonal value
          std::fill_n(dst_col, M, val_off);
          // overwrite the diagonal element (at offset m2 inside this block)
          dst_col[m2] = val_diag;
        }
      }
    }

    std::vector<double> upper1, lower1, zero1;
    FlatMatrix sigma1;
    for (size_t k = 0; k < kMax; ++k) {
      std::size_t nk = M * (k + 1);

      upper1.resize(nk); lower1.resize(nk); zero1.resize(nk);
      std::memcpy(upper1.data(), upper.data(), nk * sizeof(double));
      std::memcpy(lower1.data(), lower.data(), nk * sizeof(double));
      std::memcpy(zero1.data(), zero.data(), nk * sizeof(double));

      sigma1.resize(nk, nk);               // reuses capacity if already large enough
      // copy top-left into sigma1 (column-major): one memcpy per column
      for (std::size_t c = 0; c < nk; ++c) {
        const double* src = sigma.data_ptr() + c * fullN; // start of col c in sigma
        double* dst = sigma1.data_ptr() + c * nk;         // start of col c in sigma1
        std::memcpy(dst, src, nk * sizeof(double));
      }

      auto a = pmvnormcpp(lower1, upper1, zero1, sigma1,
                          1024, 16384, 8, 1e-4, 0.0, 314159, true);
      pcumreject[k] = 1.0 - a.prob;
    }
  }

  std::vector<double> preject(kMax);
  preject[0] = pcumreject[0];
  for (size_t k = 1; k < kMax; ++k) {
    preject[k] = pcumreject[k] - pcumreject[k-1];
  }

  return preject;
}


//' @title Exit Probabilities for a Multi-Arm Multi-Stage Design
//' @description Computes the exit (rejection) probabilities for a multi-arm
//' multi-stage design.
//'
//' @param M Number of active treatment arms.
//' @param r Randomization ratio of each active arm to the common control.
//' @param theta A vector of length \eqn{M} representing the true treatment
//'   effects for each active arm versus the common control.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param kMax Number of sequential looks.
//' @param b A vector of critical values (length \code{kMax}).
//' @param I A vector of information levels (length \code{kMax}) for any active
//'   arm versus the common control.
//'
//' @details
//' The function assumes a multivariate normal distribution for the Wald
//' statistics and all active arms share the same information level.
//'
//' @return A vector \code{exitProb} of length \code{kMax} containing the
//' probability of rejection at each look.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @examples
//'
//' # Setup: 2 active arms vs control and 3 sequential looks.
//'
//' # Information levels: equal spacing over 3 looks based on a maximum of
//' # 95 patients per arm, SD = 1.0
//' I <- 95 / (2 * 1.0^2) * seq(1, 3)/3
//'
//' # O'Brien-Fleming critical values
//' b <- c(3.886562, 2.748214, 2.243907)
//'
//' # Type I error under the global null hypothesis
//' p0 <- exitprob_mams(M = 2, theta = c(0, 0), kMax = 3, b = b, I = I)
//' cumsum(p0)
//'
//' # Power under alternative: Treatment effects of 0.3 and 0.5
//' p1 <- exitprob_mams(M = 2, theta = c(0.3, 0.5), kMax = 3, b = b, I = I)
//' cumsum(p1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector exitprob_mams(
    const int M = NA_INTEGER,
    const double r = 1,
    const Rcpp::NumericVector& theta = NA_REAL,
    const bool corr_known = true,
    const int kMax = NA_INTEGER,
    const Rcpp::NumericVector& b = NA_REAL,
    const Rcpp::NumericVector& I = NA_REAL) {

  size_t M_ = static_cast<size_t>(M);
  size_t kMax_ = static_cast<size_t>(kMax);
  std::vector<double> thetaVec(theta.begin(), theta.end());
  FlatMatrix bMat(M_, kMax_);
  for (size_t k = 0; k < kMax_; ++k) {
    for (size_t m = 0; m < M_; ++m) {
      bMat(m, k) = b[k];
    }
  }
  std::vector<double> IVec(I.begin(), I.end());
  auto v = exitprob_mams_cpp(M_, r, thetaVec, corr_known, kMax_, bMat, IVec);
  return Rcpp::wrap(v);
}


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
    const std::vector<unsigned char>& efficacyStopping) {

  if (M < 1) throw std::invalid_argument("M should be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r should be positive");
  if (k < 1) throw std::invalid_argument("k should be at least 1");

  // infoRates: if missing create 1/k, 2/k, ..., k/k
  size_t kMax = k;
  std::vector<double> infoRates;
  if (none_na(informationRates)) {
    kMax = informationRates.size();
    if (kMax < k)
      throw std::invalid_argument("Insufficient length for informationRates");

    infoRates = informationRates;
    if (infoRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(infoRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (infoRates.back() > 1.0)
      throw std::invalid_argument("informationRates must not exceed 1");
  } else {
    infoRates.resize(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(kMax);
    }
  }

  // spendTime: default to infoRates if missing
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() < kMax)
      throw std::invalid_argument("Insufficient length for spendingTime");

    spendTime = spendingTime;
    if (spendTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendTime.back() > 1.0)
      throw std::invalid_argument("spendingTime must not exceed 1");
  } else {
    spendTime = infoRates;
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() < kMax)
      throw std::invalid_argument("Insufficient length for efficacyStopping");

    effStopping = efficacyStopping;
    if (effStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
  } else {
    effStopping.assign(kMax, 1);
  }

  // asf (alpha spending function) to lower-case
  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  // userAlphaSpending checks when asf == "user"
  if (asf == "user") {
    if (!none_na(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be specified");
    if (userAlphaSpending.size() < kMax)
      throw std::invalid_argument("Insufficient length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending.back() > alpha)
      throw std::invalid_argument("userAlphaSpending must not exceed alpha");
  }

  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }

  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  if (asf == "of" || asf == "p" || asf == "wt" || asf == "none") {
    if (infoRates.back() != 1.0) {
      throw std::invalid_argument(
          "informationRates must end with 1 for OF, P, WT, or NONE");
    }
    if (spendTime.back() != 1.0) {
      throw std::invalid_argument(
          "spendingTime must end with 1 for OF, P, WT, or NONE");
    }
  }


  std::vector<double> criticalValues(kMax);
  FlatMatrix b(M, kMax);
  b.fill(6.0);

  if (asf == "none") {
    std::vector<double> zero(M, 0.0);
    double* colptr = b.data_ptr() + (kMax - 1) * M;
    auto f = [&](double x)->double {
      std::fill_n(colptr, M, x);
      auto v = exitprob_mams_cpp(M, r, zero, corr_known, kMax, b, infoRates);
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    for (size_t i = 0; i < kMax-1; ++i) criticalValues[i] = 6.0;
    criticalValues[kMax-1] = brent(f, 0.0, 6.0, 1e-6);
    return subset(criticalValues, 0, k);
  }

  if (asf == "of" || asf == "p" || asf == "wt") {
    double Delta;
    if (asf == "of") Delta = 0.0;
    else if (asf == "p") Delta = 0.5;
    else Delta = parameterAlphaSpending; // parameterAlphaSpending holds delta for WT

    // for a given multiplier, compute cumulative upper exit probability - alpha
    std::vector<double> u0(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      u0[i] = std::pow(infoRates[i], Delta - 0.5);
    }

    std::vector<double> zero(M, 0.0);

    auto f = [&](double x)->double {
      for (size_t i = 0; i < kMax; ++i) {
        if (!effStopping[i]) continue;
        double val = x * u0[i];

        double* colptr = b.data_ptr() + i * M;
        std::fill_n(colptr, M, val); // contiguous write of M doubles
      }

      auto v = exitprob_mams_cpp(M, r, zero, corr_known, kMax, b, infoRates);
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    double cwt = brent(f, 0.0, 10.0, 1e-6);
    for (size_t i = 0; i < kMax; ++i) {
      double val = cwt * u0[i];
      if (!effStopping[i]) val = 6.0;
      criticalValues[i] = val;
    }
    return subset(criticalValues, 0, k);
  }

  if (asf == "sfof" || asf == "sfp" || asf == "sfkd" ||
      asf == "sfhsd" || asf == "user") {
    // stage 1
    double cumAlpha;
    if (asf == "user") cumAlpha = userAlphaSpending[0];
    else cumAlpha = errorSpentcpp(spendTime[0], alpha, asf, parameterAlphaSpending);

    if (!effStopping[0]) criticalValues[0] = 6.0;
    else {
      std::vector<double> mean(M, 0.0);
      FlatMatrix sigma(M, M);
      for (size_t j = 0; j < M; ++j) {
        for (size_t i = 0; i < M; ++i) {
          sigma(i, j) = (i == j) ? 1.0 : r / (r + 1.0);
        }
      }
      criticalValues[0] = qmvnormcpp(1.0 - cumAlpha, mean, sigma,
                                     1024, 16384, 8, 1e-4, 0.0, 314159, true);
    }

    // Preallocate reusable buffers used by the root-finding lambda
    std::vector<double> theta_vec(M, 0.0);

    // subsequent stages
    for (size_t k1 = 1; k1 < kMax; ++k1) {
      // determine cumulative alpha at this stage
      if (asf == "user") cumAlpha = userAlphaSpending[k1];
      else cumAlpha = errorSpentcpp(spendTime[k1], alpha, asf,
                                    parameterAlphaSpending);

      if (!effStopping[k1]) {
        criticalValues[k1] = 6.0;
        continue;
      }

      // Fill columns k1-1 with their known critical values
      double val = criticalValues[k1 - 1];
      double* col_ptr = b.data_ptr() + (k1 - 1) * M;   // column-major contiguous
      std::fill_n(col_ptr, M, val);              // fast contiguous write

      // Define lambda that only sets the last column of b
      double* last_col = b.data_ptr() + k1 * M;
      auto f = [&](double x)->double {
        // set the last column to the current candidate critical value
        std::fill_n(last_col, M, x);  // fill last column fast
        auto v = exitprob_mams_cpp(M, r, theta_vec, corr_known,
                                   k1 + 1, b, infoRates);
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - cumAlpha;
      };

      double f_6 = f(6.0);
      if (f_6 > 0.0) { // no alpha spent at current visit
        criticalValues[k1] = 6.0;
      } else {
        auto f_for_brent = [&](double x)->double {
          if (x == 6.0) return f_6; // avoid recomputation at 6.0
          return f(x);
        };

        criticalValues[k1] = brent(f_for_brent, 0.0, 6.0, 1e-6);
      }
    }

    return subset(criticalValues, 0, k);
  }

  throw std::invalid_argument("Invalid value for typeAlphaSpending");
}


//' @title Efficacy Boundaries for a Multi-Arm Multi-Stage Design
//' @description Calculates the efficacy stopping boundaries for a multi-arm
//' multi-stage design.
//'
//' @param M Number of active treatment arms.
//' @param r Randomization ratio of each active arm to the common control.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is assumed.
//' @param k The index of the current look.
//' @param informationRates A numeric vector of information rates up to the
//'   current look. Values must be strictly increasing and \eqn{\le 1}.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A numeric vector of length \eqn{k} specifying the
//'   error spending time at each analysis. Values must be strictly increasing
//'   and \eqn{\le 1}. If omitted, defaults to \code{informationRates}.
//' @inheritParams param_efficacyStopping
//'
//' @details
//' The function determines critical values by solving for the boundary that
//' satisfies the alpha-spending requirement.
//'
//' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
//' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
//' and \code{spendingTime} must be of full length \code{kMax}, and
//' \code{informationRates} and \code{spendingTime} must end with 1.
//'
//' @return A numeric vector of length \eqn{k} containing the critical
//' values (on the standard normal Z-scale) for each analysis up to the
//' current look.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @examples
//'
//' # Determine O'Brien-Fleming boundaries for a TSSSD with
//' # 2 active arms and 3 looks.
//' getBound_mams(M = 2, k = 3, informationRates = seq(1, 3)/3,
//'               alpha = 0.025, typeAlphaSpending = "OF")
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector getBound_mams(
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int k = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL) {

  std::vector<double> infoRates(informationRates.begin(), informationRates.end());
  std::vector<double> userAlpha(userAlphaSpending.begin(), userAlphaSpending.end());
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());
  auto effStopping = convertLogicalVector(efficacyStopping);

  auto result = getBound_mams_cpp(
    static_cast<size_t>(M), r, corr_known, static_cast<size_t>(k),
    infoRates, alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlpha, spendTime, effStopping
  );

  return Rcpp::wrap(result);
}


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
    const FlatMatrix& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime) {

  // ----------- Input Validation ----------- //
  if (std::isnan(beta) && std::isnan(IMax))
    throw std::invalid_argument("beta and IMax cannot be missing simultaneously");

  if (!std::isnan(beta) && !std::isnan(IMax))
    throw std::invalid_argument("Only one of beta and IMax should be provided");

  if (!std::isnan(IMax) && IMax <= 0)
    throw std::invalid_argument("IMax must be positive");

  if (M < 1) throw std::invalid_argument("M must be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r must be positive");
  if (!none_na(theta)) throw std::invalid_argument("theta must be provided");
  if (theta.size() != M) throw std::invalid_argument("theta should have length M");
  if (kMax < 1) throw std::invalid_argument("kMax must be at least 1");

  // Alpha and Beta must be within valid ranges
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(beta) && (beta >= 1 - alpha || beta < 0.0001)) {
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");
  }

  std::string unknown = std::isnan(beta) ? "beta" : "IMax";

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (informationRates.size() != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[kMax-1] != 1.0)
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
    if (efficacyStopping[kMax-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues.data);

  if (!missingCriticalValues && criticalValues.nrow != kMax) {
    throw std::invalid_argument("Invalid number of rows for criticalValues");
  }
  if (!missingCriticalValues && criticalValues.ncol != M) {
    throw std::invalid_argument("Invalid number of columns for criticalValues");
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
    if (userAlphaSpending[kMax-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[kMax-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // ----------- End of Input Validation ----------- //

  // set up efficacy bounds
  // construct by-level bounds
  std::vector<int> level(M * kMax);
  std::vector<int> stage(M * kMax);
  std::vector<double> efficacyBounds(M * kMax);
  for (size_t M1 = M; M1 > 0; --M1) {
    for (size_t k = 0; k < kMax; ++k) {
      size_t row = (M - M1) * kMax + k;
      level[row] = M1;
      stage[row] = k + 1;
    }
  }

  bool correctDim = (criticalValues.nrow == kMax && criticalValues.ncol == M);
  if (missingCriticalValues) {
    std::vector<double> cut(kMax);
    for (size_t M1 = M ; M1 > 0; --M1) {
      std::vector<double> zero1(M1, 0.0);
      FlatMatrix b1(M1, kMax);

      bool haybittle = false;
      if (correctDim) {
        cut = flatmatrix_get_column(criticalValues, M - M1);

        if (kMax > 1 && cut.size() == kMax) {
          bool hasNaN = false;
          for (size_t i = 0; i < kMax - 1; ++i) {
            if (std::isnan(cut[i])) { hasNaN = true; break; }
          }
          if (!hasNaN && std::isnan(cut[kMax-1])) haybittle = true;
        }
      }

      if (haybittle) { // Haybittle & Peto
        for (size_t i = 0; i < kMax - 1; ++i) {
          if (!effStopping[i]) cut[i] = 6.0;
          double* colptr = b1.data_ptr() + i * M1;
          std::fill_n(colptr, M1, cut[i]);
        }

        double* last_col = b1.data_ptr() + (kMax - 1) * M1;
        auto f = [&](double x)->double {
          std::fill_n(last_col, M1, x);
          auto v = exitprob_mams_cpp(M1, r, zero1, corr_known, kMax, b1, infoRates);
          double cpu = std::accumulate(v.begin(), v.end(), 0.0);
          return cpu - alpha;
        };

        cut[kMax - 1] = brent(f, 0.0, 6.0, 1e-6);
      } else {
        cut = getBound_mams_cpp(
          M1, r, corr_known, kMax, infoRates, alpha, asf, parameterAlphaSpending,
          userAlphaSpending, spendTime, effStopping);
      }

      std::copy_n(cut.data(), kMax, efficacyBounds.data() + (M - M1) * kMax);
    }
  } else {
    efficacyBounds = criticalValues.data; // copy from input matrix
  }

  std::vector<double> critValues(kMax);
  std::memcpy(critValues.data(), efficacyBounds.data(), kMax * sizeof(double));

  std::vector<double> zero(M, 0.0);
  FlatMatrix b(M, kMax);
  for (size_t i = 0; i < kMax; ++i) {
    double* colptr = b.data_ptr() + i * M;
    std::fill_n(colptr, M, critValues[i]);
  }
  auto v = exitprob_mams_cpp(M, r, zero, corr_known, kMax, b, infoRates);
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent[kMax-1];

  double IMax1 = IMax;
  if (unknown == "IMax") {
    double maxtheta = *std::max_element(theta.begin(), theta.end());
    std::vector<double> theta_norm(M);
    for (size_t i = 0; i < M; ++i) theta_norm[i] = theta[i] / maxtheta;
    std::vector<double> theta1(M);
    auto f = [&](double x)->double {
      for (size_t i = 0; i < M; ++i) theta1[i] = theta_norm[i] * x;
      auto v = exitprob_mams_cpp(M, r, theta1, true, kMax, b, infoRates);
      double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
      return overallReject - (1.0 - beta);
    };

    double drift = brent(f, 0.0, 6.0, 1e-6);
    IMax1 = sq(drift / maxtheta);
    std::vector<double> info(kMax);
    for (size_t i = 0; i < kMax; ++i) info[i] = infoRates[i] * IMax1;
    v = exitprob_mams_cpp(M, r, theta, true, kMax, b, info);
  } else {
    std::vector<double> info = infoRates;
    for (double &x : info) x *= IMax1;
    v = exitprob_mams_cpp(M, r, theta, true, kMax, b, info);
  }

  // output the results
  std::vector<double> information(kMax);
  std::vector<double> efficacyTheta(kMax);
  std::vector<double> efficacyP(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    information[i] = IMax1 * infoRates[i];
    efficacyTheta[i] = critValues[i] / std::sqrt(information[i]);
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
  }

  // stagewise exit probabilities under H1
  std::vector<double> cpu(kMax);
  std::partial_sum(v.begin(), v.end(), cpu.begin());
  double overallReject = cpu[kMax-1];

  for (size_t i = 0; i < kMax; ++i) {
    if (critValues[i] == 6.0) effStopping[i] = 0;
  }


  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha1, "alpha");
  overallResults.push_back(M, "M");
  overallResults.push_back(r, "r");
  overallResults.push_back(corr_known, "corr_known");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(IMax1, "information");

  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(v), "rejectPerStage");
  byStageResults.push_back(std::move(cpu), "cumulativeRejection");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(efficacyTheta), "efficacyTheta");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(information), "information");
  byStageResults.push_back(std::move(effStopping), "efficacyStopping");

  DataFrameCpp byLevelBounds;
  byLevelBounds.push_back(std::move(level), "level");
  byLevelBounds.push_back(std::move(stage), "stage");
  byLevelBounds.push_back(std::move(efficacyBounds), "efficacyBounds");

  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(spendingTime, "spendingTime");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(byLevelBounds), "byLevelBounds");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");
  return result;
}


//' @title Power and Sample Size for a Multi-Arm Multi-Stage Design
//' @description Computes either the maximum information and stopping
//' boundaries for a multi-arm multi-stage design, or
//' the achieved power when the maximum information and stopping boundaries
//' are provided.
//'
//' @param beta Type II error rate. Provide either \code{beta} or \code{IMax};
//'   the other should be missing.
//' @param IMax Maximum information for any active arm versus the common
//'   control. Provide either \code{IMax} or \code{beta}; the other should
//'   be missing.
//' @param theta A vector of length \eqn{M} representing the true treatment
//'   effects for each active arm versus the common control. The global null
//'   is \eqn{\theta_i = 0} for all \eqn{i}, and alternatives are one-sided:
//'   \eqn{\theta_i > 0} for at least one \eqn{i = 1, \ldots, M}.
//' @param M Number of active treatment arms.
//' @param r Randomization ratio of each active arm to the common control.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param kMax Number of sequential looks.
//' @param informationRates A numeric vector of information rates fixed
//'   before the trial. If unspecified, defaults to \eqn{(1:kMax) / kMax}.
//' @inheritParams param_efficacyStopping
//' @param criticalValues The matrix of by-level upper boundaries on the
//'   max z-test statistic scale for efficacy stopping.
//'   The first column is for level \code{M}, the second column is for
//'   level \code{M - 1}, and so on, with the last column for level 1.
//'   If left unspecified, the critical values will be computed based
//'   on the specified alpha spending function.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A numeric vector of length \code{kMax} specifying the
//'   error spending time at each analysis. Values must be strictly increasing
//'   and ends at 1. If omitted, defaults to \code{informationRates}.
//'
//' @return An S3 object of class \code{mams} with the following components:
//'
//' * \code{overallResults}: A data frame containing:
//'     - \code{overallReject}: Overall probability of rejecting the global
//'       null hypothesis.
//'     - \code{alpha}: Overall significance level.
//'     - \code{M}: Number of active arms.
//'     - \code{r}: Randomization ratio per active arm versus control.
//'     - \code{corr_known}: Whether the correlation among Wald statistics
//'       was assumed known.
//'     - \code{kMax}: Number of stages.
//'     - \code{information}: Maximum information for any active arm versus
//'       control.
//'
//' * \code{byStageResults}: A data frame containing:
//'     - \code{informationRates}: Information rates at each analysis.
//'     - \code{efficacyBounds}: Efficacy boundaries on the max Z-scale.
//'     - \code{rejectPerStage}: Probability of efficacy stopping at each stage.
//'     - \code{cumulativeRejection}: Cumulative probability of efficacy stopping.
//'     - \code{cumulativeAlphaSpent}: Cumulative alpha spent.
//'     - \code{efficacyTheta}: Efficacy boundaries on the parameter scale.
//'     - \code{efficacyP}: Efficacy boundaries on the p-value scale.
//'     - \code{information}: Cumulative information for any active arm versus
//'       control at each analysis.
//'     - \code{efficacyStopping}: Indicator of whether efficacy stopping
//'       is permitted at each stage.
//'
//' * \code{settings}: A list of input settings:
//'     - \code{typeAlphaSpending}: Type of alpha spending function.
//'     - \code{parameterAlphaSpending}: Parameter value for the chosen
//'       alpha spending function.
//'     - \code{userAlphaSpending}: User-specified alpha spending values.
//'     - \code{spendingTime}: Error-spending times at each analysis.
//'
//' * \code{byLevelBounds}: A data frame containing the efficacy boundaries
//'   for each level of testing (i.e., number of active arms remaining) and
//'   each stage. Columns include:
//'     - \code{level}: Number of active arms remaining (1 to \eqn{M}).
//'     - \code{stage}: Stage index (1 to \code{kMax}).
//'     - \code{efficacyBounds}: Efficacy boundaries on the max Z-scale
//'       for the given level and stage.
//'
//' @details If \code{corr_known} is \code{FALSE}, critical boundaries are
//' computed assuming independence among the Wald statistics in each stage
//' (a conservative assumption). Power calculations, however, use the
//' correlation implied by the randomization ratio \eqn{r}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @examples
//'
//' # Example 1: obtain the maximum information given power
//' (design1 <- getDesign_mams(
//'   beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   kMax = 3, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesign_mams(
//'   IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   kMax = 3, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getDesign_mams(
    const double beta = NA_REAL,
    const double IMax = NA_REAL,
    const Rcpp::NumericVector& theta = NA_REAL,
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::Nullable<Rcpp::NumericMatrix> criticalValues = R_NilValue,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {

  std::vector<double> thetaVec(theta.begin(), theta.end());
  std::vector<double> infoRates(informationRates.begin(), informationRates.end());
  auto effStopping = convertLogicalVector(efficacyStopping);
  std::vector<double> userAlpha(userAlphaSpending.begin(), userAlphaSpending.end());
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

  auto cpp_result = getDesign_mams_cpp(
    beta, IMax, thetaVec, static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(kMax), infoRates, effStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, spendTime
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "mams";
  return result;
}


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
    const FlatMatrix& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
    const size_t MNew,
    const std::vector<int>& selected,
    const double rNew,
    const size_t kNew,
    const std::vector<double>& informationRatesNew,
    const std::vector<unsigned char>& efficacyStoppingNew,
    const std::string& typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double>& spendingTimeNew) {

  // ----------- Start of Input Validation ----------- //
  if (std::isnan(betaNew) && std::isnan(INew))
    throw std::invalid_argument("betaNew and INew cannot be missing simultaneously");
  if (!std::isnan(betaNew) && !std::isnan(INew))
    throw std::invalid_argument("Only one of betaNew and INew should be provided");
  if (!std::isnan(INew) && INew <= 0.0)
    throw std::invalid_argument("INew must be positive");

  if (M < 1) throw std::invalid_argument("M must be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r must be positive");
  if (L < 1) throw std::invalid_argument("L must be at least 1");
  if (!none_na(zL)) throw std::invalid_argument("zL must be provided");
  if (zL.size() != M) throw std::invalid_argument("Invalid length for zL");
  if (!none_na(theta)) throw std::invalid_argument("theta must be provided");
  if (theta.size() != M) throw std::invalid_argument("Invalid length for theta");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(betaNew) && (betaNew < 0.0001 || betaNew >= 1)) {
    throw std::invalid_argument("betaNew must lie in [0.0001, 1)");
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
    if (informationRates[kMax-1] != 1.0)
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
    if (efficacyStopping[kMax-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues.data);

  if (!missingCriticalValues && criticalValues.nrow != kMax) {
    throw std::invalid_argument("Invalid number of rows for criticalValues");
  }
  if (!missingCriticalValues && criticalValues.ncol != M) {
    throw std::invalid_argument("Invalid number of columns for criticalValues");
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
    if (userAlphaSpending[kMax-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[kMax-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // ----------- New Design Input Validation ----------- //
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

  if (MullerSchafer) {
    if (kNew < 1) throw std::invalid_argument("kNew must be at least 1");

    // informationRatesNew: default to (1:kNew)/kNew if missing
    infoRatesNew.resize(kNew);
    if (none_na(informationRatesNew)) {
      if (informationRatesNew.size() != kNew)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew[kNew-1] != 1.0)
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
      if (efficacyStoppingNew[kNew-1] != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew; // copy
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    if (!(asfNew == "of" || asfNew == "sfof" || asfNew == "sfp" ||
        asfNew == "sfkd" || asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }
    if ((asfNew == "sfkd" || asfNew == "sfhsd") &&
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
      if (spendingTimeNew[kNew-1] != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = infoRatesNew;
    }
  }
  // ----------- End of Input Validation ----------- //

  // obtain critical values for the primary trial
  // by level critical values for the integrated trial
  // construct by-level bounds for the primary trial
  std::vector<int> level1(M * kMax);
  std::vector<int> stage1(M * kMax);
  std::vector<double> efficacyBounds1(M * kMax);
  for (size_t M1 = M; M1 > 0; --M1) {
    for (size_t k = 0; k < kMax; ++k) {
      size_t row = (M - M1) * kMax + k;
      level1[row] = M1;
      stage1[row] = k + 1;
    }
  }

  bool correctDim = (criticalValues.nrow == kMax && criticalValues.ncol == M);
  if (missingCriticalValues) {
    std::vector<double> cut(kMax);
    for (size_t M1 = M ; M1 > 0; --M1) {
      std::vector<double> zero1(M1, 0.0);
      FlatMatrix b1(M1, kMax);

      bool haybittle = false;
      if (correctDim) {
        cut = flatmatrix_get_column(criticalValues, M - M1);

        if (kMax > 1 && cut.size() == kMax) {
          bool hasNaN = false;
          for (size_t i = 0; i < kMax - 1; ++i) {
            if (std::isnan(cut[i])) { hasNaN = true; break; }
          }
          if (!hasNaN && std::isnan(cut[kMax-1])) haybittle = true;
        }
      }

      if (haybittle) { // Haybittle & Peto
        for (size_t i = 0; i < kMax - 1; ++i) {
          if (!effStopping[i]) cut[i] = 6.0;
          double* colptr = b1.data_ptr() + i * M1;
          std::fill_n(colptr, M1, cut[i]);
        }

        double* last_col = b1.data_ptr() + (kMax - 1) * M1;
        auto f = [&](double x)->double {
          std::fill_n(last_col, M1, x);
          auto v = exitprob_mams_cpp(M1, r, zero1, corr_known, kMax, b1, infoRates);
          double cpu = std::accumulate(v.begin(), v.end(), 0.0);
          return cpu - alpha;
        };

        cut[kMax - 1] = brent(f, 0.0, 6.0, 1e-6);
      } else {
        cut = getBound_mams_cpp(
          M1, r, corr_known, kMax, infoRates, alpha, asf, parameterAlphaSpending,
          userAlphaSpending, spendTime, effStopping);
      }

      std::copy_n(cut.data(), kMax, efficacyBounds1.data() + (M - M1) * kMax);
    }
  } else {
    efficacyBounds1 = criticalValues.data; // copy from input matrix
  }

  std::vector<double> critValues(kMax);
  std::memcpy(critValues.data(), efficacyBounds1.data(), kMax * sizeof(double));

  std::vector<double> zero(M, 0.0);
  FlatMatrix b(M, kMax);
  for (size_t i = 0; i < kMax; ++i) {
    double* colptr = b.data_ptr() + i * M;
    std::fill_n(colptr, M, critValues[i]);
  }
  auto v = exitprob_mams_cpp(M, r, zero, corr_known, kMax, b, infoRates);
  double p0 = std::accumulate(v.begin(), v.end(), 0.0);
  double alpha1 = missingCriticalValues ? alpha : p0;

  // compute conditional alpha, conditional power, and predictive power
  size_t k1 = kMax - L;
  std::vector<double> s1(k1);
  FlatMatrix c1(M, k1);
  for (size_t i = 0; i < k1; ++i) {
    s1[i] = (infoRates[L + i] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
    double* colptr = c1.data_ptr() + i * M; // start of column i
    if (effStopping[L + i]) {
      double cut = critValues[L + i];
      double r1 = infoRates[L - 1] / infoRates[L + i];
      double sqrt_r1 = std::sqrt(r1);
      double denom = std::sqrt(1.0 - r1);
      // write contiguous column
      for (size_t m = 0; m < M; ++m) {
        colptr[m] = (cut - zL[m] * sqrt_r1) / denom;
      }
    } else {
      std::fill_n(colptr, M, 6.0);
    }
  }

  // secondary trial information for original design
  double INew1 = IMax * (1.0 - infoRates[L - 1]);
  std::vector<double> I1(k1);
  for (size_t i = 0; i < k1; ++i) {
    I1[i] = INew1 * s1[i];
  }

  // conditional type I error
  auto v0 = exitprob_mams_cpp(M, r, zero, corr_known, k1, c1, I1);
  double c_alpha = std::accumulate(v0.begin(), v0.end(), 0.0);

  // conditional power
  auto v1 = exitprob_mams_cpp(M, r, theta, true, k1, c1, I1);
  double c_power = std::accumulate(v1.begin(), v1.end(), 0.0);

  std::vector<double> information1(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    information1[i] = IMax * infoRates[i];
  }

  // secondary trial design and power calculations
  double IL = IMax * infoRates[L - 1];
  double sqrtIL = std::sqrt(IL);
  std::vector<double> zero2(MNew, 0.0);
  std::vector<double> theta2 = subset(theta, selectedNew);

  size_t k2 = MullerSchafer ? kNew : k1;
  std::vector<double> s2(k2);
  std::string asf2;
  double asfpar2;
  std::vector<double> spendTime2;
  std::vector<unsigned char> effStopping2(k2);
  std::vector<double> cpu0(k2);
  if (!MullerSchafer) {
    // number of secondary trial looks is the same as the original trial
    // spend the conditional type I error as in the original trial
    s2 = s1;
    asf2 = "user";
    asfpar2 = NaN;
    spendTime2 = {NaN};
    std::memcpy(effStopping2.data(), effStopping.data() + L, k2);
    std::partial_sum(v0.begin(), v0.end(), cpu0.begin());
  } else {
    s2 = infoRatesNew;
    asf2 = asfNew;
    asfpar2 = parameterAlphaSpendingNew;
    spendTime2 = spendTimeNew;
    effStopping2 = effStoppingNew;
    if (asf2 != "none" && asf2 != "of") {
      for (size_t i = 0; i < k2; ++i) {
        cpu0[i] = errorSpentcpp(spendTimeNew[i], c_alpha, asfNew, asfpar2);
      }
    }
  }

  std::vector<double> I2(k2); // information levels for secondary trial
  std::vector<double> Ic(k2); // information levels for integrated trial
  std::vector<double> sqrtI2(k2), sqrtIc(k2);
  std::vector<double> critValues2(k2, 6.0);
  FlatMatrix c2(MNew, k2);
  c2.fill(6.0);

  std::vector<double> zscaled(MNew);
  for (size_t j = 0; j < MNew; ++j) zscaled[j] = zL[selectedNew[j]] * sqrtIL;

  if (std::isnan(betaNew)) {
    for (size_t i = 0; i < k2; ++i) {
      I2[i] = INew * s2[i];
      Ic[i] = I2[i] + IL;
      sqrtI2[i] = std::sqrt(I2[i]);
      sqrtIc[i] = std::sqrt(Ic[i]);
    }

    if (asf2 == "of") {
      auto g = [&c2, &I2, &sqrtI2, &sqrtIc, &zscaled, &zero2, &effStopping2,
                k2, c_alpha, MNew, rNew, corr_known]
      (double x)->double {
        double col_const = x * sqrtIc[k2 - 1];
        for (size_t i = 0; i < k2; ++i) {
          double* colptr = c2.data_ptr() + i * MNew;
          if (effStopping2[i]) {
            double denom = sqrtI2[i];
            for (size_t j = 0; j < MNew; ++j) {
              colptr[j] = (col_const - zscaled[j]) / denom;
            }
          } else {
            std::fill_n(colptr, MNew, 6.0);
          }
        }

        auto v = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, c2, I2);
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - c_alpha;
      };

      double cof = brent(g, 0.0, 6.0, 1e-6);
      double col_const = cof * sqrtIc[k2 - 1];
      for (size_t i = 0; i < k2; ++i) {
        double* colptr = c2.data_ptr() + i * MNew;
        if (effStopping2[i]) {
          critValues2[i] = col_const / sqrtIc[i];
          double denom = sqrtI2[i];
          for (size_t j = 0; j < MNew; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }
        } else {
          critValues2[i] = 6.0;
          std::fill_n(colptr, MNew, 6.0);
        }
      }
    } else if (asf2 == "none") {
      for (size_t i = 0; i < k2 - 1; ++i) critValues2[i] = 6.0;
      double denom = sqrtI2[k2 - 1];

      auto g = [&c2, &I2, &sqrtIc, &zscaled, &zero2,
                denom, k2, c_alpha, MNew, rNew, corr_known]
      (double x)->double {
        double* colptr = c2.data_ptr() + (k2 - 1) * MNew;
        double col_const = x * sqrtIc[k2 - 1];
        for (size_t j = 0; j < MNew; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }

        auto v = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, c2, I2);
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - c_alpha;
      };

      double cof = brent(g, 0.0, 6.0, 1e-6);
      critValues2[k2 - 1] = cof;
      double* colptr = c2.data_ptr() + (k2 - 1) * MNew;
      double col_const = cof * sqrtIc[k2 - 1];
      for (size_t j = 0; j < MNew; ++j) {
        colptr[j] = (col_const - zscaled[j]) / denom;
      }
    } else {
      for (size_t i = 0; i < k2; ++i) {
        if (!effStopping2[i]) continue;
        double denom = sqrtI2[i];

        auto g = [&c2, &I2, &sqrtIc, &zscaled, &cpu0, &zero2,
                  denom, i, MNew, rNew, corr_known]
        (double x)->double {
          double col_const = x * sqrtIc[i];
          double* colptr = c2.data_ptr() + i * MNew;
          // update critical values of the secondary trial at current look
          for (size_t j = 0; j < MNew; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }

          auto v = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, i + 1, c2, I2);
          double p0 = std::accumulate(v.begin(), v.end(), 0.0);
          return p0 - cpu0[i];
        };

        double cof = brent(g, 0.0, 6.0, 1e-6);
        double col_const = cof * sqrtIc[i];
        critValues2[i] = cof;
        double* colptr = c2.data_ptr() + i * MNew;
        for (size_t j = 0; j < MNew; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }
      }
    }
  } else {
    // obtain required max information for the secondary trial given target power
    auto f = [&critValues2, &c2, &I2, &Ic, &sqrtI2, &sqrtIc,
              &zscaled, &cpu0, &zero2, &theta2, &s2, &effStopping2,
              betaNew, k2, asf2, c_alpha, IL, MNew, rNew, corr_known]
    (double Inew)->double {
      for (size_t i = 0; i < k2; ++i) {
        I2[i] = Inew * s2[i];
        Ic[i] = I2[i] + IL;
        sqrtI2[i] = std::sqrt(I2[i]);
        sqrtIc[i] = std::sqrt(Ic[i]);
      }

      if (asf2 == "of") {
        auto g = [&c2, &I2, &sqrtI2, &sqrtIc, &zscaled, &zero2, &effStopping2,
                  k2, c_alpha, MNew, rNew, corr_known]
        (double x)->double {
          double col_const = x * sqrtIc[k2 - 1];
          for (size_t i = 0; i < k2; ++i) {
            double* colptr = c2.data_ptr() + i * MNew;
            if (effStopping2[i]) {
              double denom = sqrtI2[i];
              for (size_t j = 0; j < MNew; ++j) {
                colptr[j] = (col_const - zscaled[j]) / denom;
              }
            } else {
              std::fill_n(colptr, MNew, 6.0);
            }
          }

          auto v = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, c2, I2);
          double p0 = std::accumulate(v.begin(), v.end(), 0.0);
          return p0 - c_alpha;
        };

        double cof = brent(g, 0.0, 6.0, 1e-6);
        double col_const = cof * sqrtIc[k2 - 1];
        for (size_t i = 0; i < k2; ++i) {
          double* colptr = c2.data_ptr() + i * MNew;
          if (effStopping2[i]) {
            critValues2[i] = col_const / sqrtIc[i];
            double denom = sqrtI2[i];
            for (size_t j = 0; j < MNew; ++j) {
              colptr[j] = (col_const - zscaled[j]) / denom;
            }
          } else {
            critValues2[i] = 6.0;
            std::fill_n(colptr, MNew, 6.0);
          }
        }
      } else if (asf2 == "none") {
        for (size_t i = 0; i < k2 - 1; ++i) critValues2[i] = 6.0;
        double denom = sqrtI2[k2 - 1];

        auto g = [&c2, &I2, &sqrtIc, &zscaled, &zero2,
                  denom, k2, c_alpha, MNew, rNew, corr_known]
        (double x)->double {
          double col_const = x * sqrtIc[k2 - 1];
          double* colptr = c2.data_ptr() + (k2 - 1) * MNew;
          for (size_t j = 0; j < MNew; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }

          auto v = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, c2, I2);
          double p0 = std::accumulate(v.begin(), v.end(), 0.0);
          return p0 - c_alpha;
        };

        double cof = brent(g, 0.0, 6.0, 1e-6);
        critValues2[k2 - 1] = cof;
        double col_const = cof * sqrtIc[k2 - 1];
        double* colptr = c2.data_ptr() + (k2 - 1) * MNew;
        for (size_t j = 0; j < MNew; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }
      } else {
        for (size_t i = 0; i < k2; ++i) {
          if (!effStopping2[i]) continue;
          double denom = sqrtI2[i];

          auto g = [&c2, &I2, &sqrtIc, &zscaled, &cpu0, &zero2,
                    denom, i, MNew, rNew, corr_known]
          (double x)->double {
            double col_const = x * sqrtIc[i];
            double* colptr = c2.data_ptr() + i * MNew;
            // update critical values of the secondary trial at current look
            for (size_t j = 0; j < MNew; ++j) {
              colptr[j] = (col_const - zscaled[j]) / denom;
            }

            auto v = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, i + 1, c2, I2);
            double p0 = std::accumulate(v.begin(), v.end(), 0.0);
            return p0 - cpu0[i];
          };

          double cof = brent(g, 0.0, 6.0, 1e-6);
          double col_const = cof * sqrtIc[i];
          critValues2[i] = cof;
          double* colptr = c2.data_ptr() + i * MNew;
          for (size_t j = 0; j < MNew; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }
        }
      }

      // compute conditional power of the secondary trial
      auto v1 = exitprob_mams_cpp(MNew, rNew, theta2, true, k2, c2, I2);
      double p1 = std::accumulate(v1.begin(), v1.end(), 0.0);
      return p1 - (1.0 - betaNew);
    };

    double lower = 0.01 * IMax, upper = IMax;
    double fl_val = f(lower), fu_val = f(upper);
    int expand_iter = 0;
    while (fl_val * fu_val > 0.0 && expand_iter < 60) {
      lower = upper; fl_val = fu_val;
      upper *= 2.0;  fu_val = f(upper);
      ++expand_iter;
    }
    if (fl_val * fu_val > 0.0) throw std::runtime_error(
        "Unable to bracket required new information; check inputs");

    auto f_for_brent = [&f, lower, upper, fl_val, fu_val](double x)->double {
      if (x == lower) return fl_val;
      if (x == upper) return fu_val;
      return f(x);
    };

    INew = brent(f_for_brent, lower, upper, 1e-6);
  }

  // compute conditional power of the secondary trial
  auto v2 = exitprob_mams_cpp(MNew, rNew, theta2, true, k2, c2, I2);
  double p2 = std::accumulate(v2.begin(), v2.end(), 0.0);
  std::vector<double> cpu1(k2);
  std::partial_sum(v2.begin(), v2.end(), cpu1.begin());

  std::vector<int> hypothesis2(MNew * k2);
  std::vector<int> stage2(MNew * k2);
  std::vector<double> bounds2(MNew * k2);
  for (size_t i = 0; i < k2; ++i) {
    for (size_t j = 0; j < MNew; ++j) {
      size_t row = j * k2 + i;
      hypothesis2[row] = selectedNew[j] + 1;
      stage2[row] = i + 1;
      bounds2[row] = c2(j, i);
    }
  }

  // integrated trial information and critical values
  size_t kc = L + k2;
  std::vector<double> Ic_full(kc);
  std::copy_n(information1.data(), L, Ic_full.data());
  std::copy_n(Ic.data(), k2, Ic_full.data() + L);

  double IMaxc = Ic_full.back();
  std::vector<double> infoRates_full(kc);
  for (size_t i = 0; i < kc; ++i) infoRates_full[i] = Ic_full[i] / IMaxc;

  std::vector<double> critValues_full(kc);
  std::copy_n(critValues.data(), L, critValues_full.data());
  std::copy_n(critValues2.data(), k2, critValues_full.data() + L);


  // critical values for each intersection hypothesis
  size_t ntests = (1 << MNew) - 1; // 2^MNew - 1
  std::vector<std::string> intersectHyp(ntests * kc);
  std::vector<int> stage(ntests * kc);
  std::vector<double> efficacyBounds(ntests * kc);

  std::vector<int> idx_in_selected(M, -1);
  for (size_t t = 0; t < selectedNew.size(); ++t)
    idx_in_selected[selectedNew[t]] = static_cast<int>(t);

  for (size_t h = 0; h < ntests; ++h) {
    size_t number = ntests - h;

    std::vector<size_t> primary;
    std::vector<size_t> selectedNew2;
    primary.reserve(M); selectedNew2.reserve(MNew);

    // Count M2 by popcount
    size_t M2 = 0;
    for (size_t b = 0; b < MNew; ++b) if ((number >> (MNew - 1 - b)) & 1u) ++M2;
    size_t M1 = M2 + (M - MNew);

    for (size_t j = 0; j < M; ++j) {
      int pos = idx_in_selected[j];
      if (pos >= 0) {
        // arm j is among selectedNew; check corresponding bit in 'number'
        if ((number >> (MNew - 1 - pos)) & 1u) {
          primary.push_back(j);
          selectedNew2.push_back(j);
        }
      } else {
        // non-selected arm always included
        primary.push_back(j);
      }
    }

    std::vector<double> critValues(kMax);
    std::copy_n(efficacyBounds1.data() + (M - M1) * kMax, kMax, critValues.data());

    // compute conditional alpha for MAMS with the given M1 hypotheses
    FlatMatrix c1(M1, k1);
    for (size_t i = 0; i < k1; ++i) {
      double col_const = critValues[L + i];
      double r1 = infoRates[L - 1] / infoRates[L + i];
      double sqrt_r1 = std::sqrt(r1);
      double denom = std::sqrt(1.0 - r1);
      double* colptr = c1.data_ptr() + i * M1; // start of column i
      if (effStopping[L + i]) {
        for (size_t j = 0; j < M1; ++j) {
          size_t m = primary[j];
          colptr[j] = (col_const - zL[m] * sqrt_r1) / denom;
        }
      } else {
        std::fill_n(colptr, M1, 6.0);
      }
    }

    // conditional type I error
    std::vector<double> zero1(M1, 0.0);
    auto v0 = exitprob_mams_cpp(M1, r, zero1, corr_known, k1, c1, I1);
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
    std::vector<double> critValues2(k2, 6.0);
    FlatMatrix c2(M2, k2);
    c2.fill(6.0);

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
              colptr[j] = (col_const - zscaled[j]) / denom;
            }
          } else {
            std::fill_n(colptr, M2, 6.0);
          }
        }

        auto v = exitprob_mams_cpp(M2, rNew, zero2, corr_known, k2, c2, I2);
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - c_alpha;
      };

      double cof = brent(g, 0.0, 6.0, 1e-6);
      double col_const = cof * sqrtIc[k2 - 1];
      for (size_t i = 0; i < k2; ++i) {
        if (effStopping2[i]) {
          critValues2[i] = col_const / sqrtIc[i];
        } else {
          critValues2[i] = 6.0;
        }
      }
    } else if (asf2 == "none") {
      for (size_t i = 0; i < k2 - 1; ++i) critValues2[i] = 6.0;
      double denom = sqrtI2[k2 - 1];

      auto g = [&c2, &I2, &sqrtIc, &zscaled, &zero2,
                denom, k2, c_alpha, M2, rNew, corr_known]
      (double x)->double {
        double col_const = x * sqrtIc[k2 - 1];
        double* colptr = c2.data_ptr() + (k2 - 1) * M2;
        for (size_t j = 0; j < M2; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }

        auto v = exitprob_mams_cpp(M2, rNew, zero2, corr_known, k2, c2, I2);
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - c_alpha;
      };

      double cof = brent(g, 0.0, 6.0, 1e-6);
      critValues2[k2 - 1] = cof;
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

          auto v = exitprob_mams_cpp(M2, rNew, zero2, corr_known, i + 1, c2, I2);
          double p0 = std::accumulate(v.begin(), v.end(), 0.0);
          return p0 - cpu0[i];
        };

        critValues2[i] = brent(g, 0.0, 6.0, 1e-6);
      }
    }

    std::vector<double> critValues_full(kc);
    std::copy_n(critValues.data(), L, critValues_full.data());
    std::copy_n(critValues2.data(), k2, critValues_full.data() + L);
    for (size_t i = 0; i < kc; ++i) {
      size_t row = h * kc + i;

      std::string s;
      s.reserve(3 + selectedNew2.size() * 4); // rough reservation
      s.push_back('{');
      for (size_t t = 0; t < selectedNew2.size(); ++t) {
        if (t) { s.append(", "); }
        s.append(std::to_string(selectedNew2[t] + 1));
      }
      s.push_back('}');
      intersectHyp[row] = std::move(s);

      stage[row] = i + 1;
      efficacyBounds[row] = critValues_full[i];
    }
  }

  ListCpp des1;
  des1.push_back(M, "M");
  des1.push_back(r, "r");
  des1.push_back(corr_known, "corr_known");
  des1.push_back(L, "L");
  des1.push_back(zL, "zL");
  des1.push_back(theta, "theta");
  des1.push_back(IMax, "maxInformation");
  des1.push_back(kMax, "kMax");
  des1.push_back(std::move(infoRates), "informationRates");
  des1.push_back(std::move(critValues), "efficacyBounds");
  des1.push_back(std::move(information1), "information");
  des1.push_back(alpha1, "alpha");
  des1.push_back(c_alpha, "conditionalAlpha");
  des1.push_back(c_power, "conditionalPower");
  des1.push_back(MullerSchafer, "MullerSchafer");

  DataFrameCpp byLevelBounds1;
  byLevelBounds1.push_back(std::move(level1), "level");
  byLevelBounds1.push_back(std::move(stage1), "stage");
  byLevelBounds1.push_back(std::move(efficacyBounds1), "efficacyBounds");
  des1.push_back(std::move(byLevelBounds1), "byLevelBounds");

  ListCpp des2;
  des2.push_back(p2, "overallReject");
  des2.push_back(c_alpha, "alpha");
  des2.push_back(MNew, "M");
  des2.push_back(rNew, "r");
  des2.push_back(selected, "selected");
  des2.push_back(corr_known, "corr_known");
  des2.push_back(k2, "kMax");
  des2.push_back(INew, "maxInformation");
  des2.push_back(std::move(s2), "informationRates");
  des2.push_back(std::move(cpu1), "cumulativeRejection");
  des2.push_back(std::move(cpu0), "cumulativeAlphaSpent");
  des2.push_back(std::move(I2), "information");
  des2.push_back(asf2, "typeAlphaSpending");
  des2.push_back(asfpar2, "parameterAlphaSpending");
  des2.push_back(spendTime2, "spendingTime");

  DataFrameCpp byHypBounds;
  byHypBounds.push_back(std::move(hypothesis2), "hypothesis");
  byHypBounds.push_back(std::move(stage2), "stage");
  byHypBounds.push_back(std::move(bounds2), "efficacyBounds");
  des2.push_back(std::move(byHypBounds), "byHypothesisBounds");

  ListCpp des3;
  des3.push_back(M, "M");
  des3.push_back(r, "r");
  des3.push_back(corr_known, "corr_known");
  des3.push_back(MNew, "MNew");
  des3.push_back(rNew, "rNew");
  des3.push_back(selected, "selected");
  des3.push_back(L, "L");
  des3.push_back(zL, "zL");
  des3.push_back(theta, "theta");
  des3.push_back(IMaxc, "maxInformation");
  des3.push_back(kc, "kMax");
  des3.push_back(std::move(infoRates_full), "informationRates");
  des3.push_back(std::move(critValues_full), "efficacyBounds");
  des3.push_back(std::move(Ic_full), "information");

  DataFrameCpp byInterBounds;
  byInterBounds.push_back(std::move(intersectHyp), "intersectionHypothesis");
  byInterBounds.push_back(std::move(stage), "stage");
  byInterBounds.push_back(std::move(efficacyBounds), "efficacyBounds");
  des3.push_back(std::move(byInterBounds), "byIntersectionBounds");

  ListCpp result;
  result.push_back(std::move(des1), "primaryTrial");
  result.push_back(std::move(des2), "secondaryTrial");
  result.push_back(std::move(des3), "integratedTrial");
  return result;
}


//' @title Adaptive Multi-Arm Multi-Stage Design
//' @description
//' Calculates the conditional power for specified incremental
//' information, given the interim results, parameter value,
//' data-dependent changes in treatment selection,
//' the error spending function, and
//' the number and spacing of interim looks. Conversely,
//' calculates the incremental information required to attain
//' a specified conditional power, given the interim results,
//' parameter value, data-dependent changes in treatment selection,
//' the error spending function, and the number and spacing of interim looks.
//'
//' @param betaNew The type II error for the secondary trial.
//' @param INew The maximum information for any active arm versus the common
//'   control in the secondary trial. Either
//'   \code{betaNew} or \code{INew} should be provided, while the other
//'   must be missing.
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
//' @param theta A vector of length \eqn{M} representing the assumed treatment
//'   effects for each active arm versus the common control. The global null
//'   is \eqn{\theta_i = 0} for all \eqn{i}, and alternatives are one-sided:
//'   \eqn{\theta_i > 0} for at least one \eqn{i = 1, \ldots, M}.
//' @param IMax Maximum information for any active arm versus the common
//'   control for the primary trial. Must be provided.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param criticalValues The matrix of by-level upper boundaries on the
//'   max z-test statistic scale for efficacy stopping for the primary trial.
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
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function,
//'   \code{"user"} for user defined spending, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param userAlphaSpending The user-defined alpha spending for the
//'   primary trial. Represents the cumulative alpha spent up to each stage.
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param MNew Number of active treatment arms in the secondary trial.
//' @param selected The indices of the selected active treatment arms for
//'   the secondary trial.
//' @param rNew Randomization ratio of each active arm to the common control
//'   in the secondary trial.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
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
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRatesNew}.
//'
//' @return An \code{adaptDesign_mams} object with three list components:
//'
//' * \code{primaryTrial}: A list of selected information for the primary
//'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{L},
//'   \code{zL}, \code{theta}, \code{maxInformation}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, \code{information},
//'   \code{alpha}, \code{conditionalAlpha}, \code{conditionalPower},
//'   \code{MullerSchafer}, and \code{byLevelBounds}.
//'
//' * \code{secondaryTrial}: A list of selected information for the secondary
//'   trial, including \code{overallReject}, \code{alpha}, \code{M}, \code{r},
//'   \code{selected}, \code{corr_known}, \code{kMax}, \code{maxInformation},
//'   \code{informationRates}, \code{cumulativeRejection},
//'   \code{cumulativeAlphaSpent}, \code{information},
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{spendingTime}, and \code{byHypothesisBounds}.
//'
//' * \code{integratedTrial}: A list of selected information for the integrated
//'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{MNew},
//'   \code{rNew}, \code{selected}, \code{L}, \code{zL}, \code{theta},
//'   \code{maxInformation}, \code{kMax}, \code{informationRates},
//'   \code{efficacyBounds}, \code{information}, and
//'   \code{byIntersectionBounds}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @seealso \code{\link{getDesign_mams}}
//'
//' @examples
//'
//' # Two active treatment arms are compared with a common control in a
//' # two-look time-to-event design using O'Brien–Fleming–type alpha spending.
//' # Suppose each active arm has a true hazard ratio of 0.75 versus control,
//' # and the total number of events across all three arms at the final analysis
//' # is 486. This corresponds to approximately 324 events for each active arm
//' # versus the common control. Under these assumptions, the trial has about
//' # 80% power to detect the treatment effect in at least one active arm.
//'
//' (des1 <- getDesign_mams(
//'   IMax = 324 / 4, theta = c(-log(0.75), -log(0.75)),
//'   M = 2, r = 1, kMax = 2, informationRates = c(1/2, 1),
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' # Now assume that, at the interim analysis, the observed hazard ratios for
//' # the two active arms versus control are 0.91 and 0.78, respectively. Using
//' # the rule “drop any arm with an observed hazard ratio > 0.9”, arm 1 is
//' # dropped. We then aim to achieve 80% conditional power to detect a hazard
//' # ratio of 0.78 for the remaining arm at the final look. The analysis below
//' # indicates that the required total number of events for arm 2 versus control
//' # at the final analysis should be increased from 324 to 535.
//'
//' (des2 <- adaptDesign_mams(
//'   betaNew = 0.2, M = 2, r = 1, corr_known = FALSE,
//'   L = 1, zL = c(-log(0.91), -log(0.78)) * sqrt(324 / 4 / 2),
//'   theta = c(-log(0.91), -log(0.78)),
//'   IMax = 324 / 4, kMax = 2, informationRates = c(1/2, 1),
//'   alpha = 0.025, typeAlphaSpending = "OF",
//'   MNew = 1, selected = 2, rNew = 1))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List adaptDesign_mams(
    double betaNew = NA_REAL,
    double INew = NA_REAL,
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int L = NA_INTEGER,
    const Rcpp::NumericVector& zL = NA_REAL,
    const Rcpp::NumericVector& theta = NA_REAL,
    const double IMax = NA_REAL,
    const int kMax = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::Nullable<Rcpp::NumericMatrix> criticalValues = R_NilValue,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int MNew = NA_INTEGER,
    const Rcpp::IntegerVector& selected = NA_INTEGER,
    const double rNew = 1,
    const int kNew = NA_INTEGER,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL) {

  auto zLVec = Rcpp::as<std::vector<double>>(zL);
  auto thetaVec = Rcpp::as<std::vector<double>>(theta);
  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto selectedNew = Rcpp::as<std::vector<int>>(selected);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  // Handle optional matrix safely
  FlatMatrix critValues;
  if (criticalValues.isNotNull()) {
    Rcpp::NumericMatrix cm(criticalValues); // unwrap
    critValues = flatmatrix_from_Rmatrix(cm);
  } else {
    critValues = FlatMatrix(1, 1);
    critValues(0, 0) = std::numeric_limits<double>::quiet_NaN(); // placeholder
  }

  auto cpp_result = adaptDesign_mams_cpp(
    betaNew, INew,  static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(L), zLVec, thetaVec, IMax,
    static_cast<size_t>(kMax), infoRates, effStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, spendTime,
    MullerSchafer, static_cast<size_t>(MNew), selectedNew, rNew,
    static_cast<size_t>(kNew), infoRatesNew, effStoppingNew,
    typeAlphaSpendingNew, parameterAlphaSpendingNew, spendTimeNew
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "adaptDesign_mams";
  return result;
}

