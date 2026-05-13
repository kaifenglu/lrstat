#include "seamless_design.h"
#include "dataframe_list.h"
#include "thread_utils.h"
#include "utilities.h"
#include "generic_design.h"
#include "mvnormr.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using std::size_t;


ExitProbSeamless exitprob_seamless_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t K,
    const std::vector<double>& b,
    const std::vector<double>& a,
    const std::vector<double>& I) {

  if (M < 1) {
    throw std::invalid_argument("M should be at least 1");
  }
  if (r <= 0.0) {
    throw std::invalid_argument("r should be positive");
  }
  if (theta.size() != M) {
    throw std::invalid_argument("theta should have length M");
  }
  if (K < 1) {
    throw std::invalid_argument("K should be at least 1");
  }
  if (b.size() < K + 1) {
    throw std::invalid_argument("Insufficient length for b");
  }
  if (a.size() < K + 1) {
    throw std::invalid_argument("Insufficient length for a");
  }
  for (size_t k = 0; k < K + 1; ++k) {
    if (b[k] < a[k]) {
      throw std::invalid_argument("Each b[k] must be >= a[k]");
    }
  }

  std::vector<double> Ivec; Ivec.reserve(K+1);
  if (none_na(I)) {
    if (I.size() < K + 1) throw std::invalid_argument("Insufficient length for I");

    Ivec.assign(I.begin(), I.begin() + K + 1);
    if (Ivec[0] <= 0.0) throw std::invalid_argument("I must be positive");
    if (any_nonincreasing(Ivec)) throw std::invalid_argument("I must be increasing");
  } else {
    Ivec.resize(K+1);
    std::iota(Ivec.begin(), Ivec.end(), 1.0);
  }

  double rho = corr_known ? r / (r + 1.0) : 0;

  std::vector<double> exitProbUpper(K + 1, 0.0);
  std::vector<double> exitProbLower(K + 1, 0.0);
  FlatMatrix exitProbByArmUpper(K + 1, M);
  FlatMatrix exitProbByArmLower(K + 1, M);
  std::vector<double> selectAsBest(M, 0.0);

  // information quantities for each pairwise comparison with control
  std::vector<double> sqrtI(K + 1);
  for (size_t k = 0; k < K + 1; ++k) {
    sqrtI[k] = std::sqrt(Ivec[k]);
  }
  double sqrtI0 = sqrtI[0];

  // mean vector of the Wald statistics for the M active arms in phase 2
  std::vector<double> mean0(M);
  for (size_t m = 0; m < M; ++m) {
    mean0[m] = theta[m] * sqrtI0;
  }

  // covariance matrix of the Wald statistics for the M active arms in phase 2
  FlatMatrix sigma0(M, M);
  sigma0.fill(rho);
  for (size_t i = 0; i < M; ++i) {
    sigma0(i, i) = 1.0;
  }

  // -------- phase 2 upper exit: reject if max Z >= b[0] --------
  {
    std::vector<double> lower0(M, -8.0), upper0(M, b[0]);
    PMVNResult out0 = pmvnormcpp(lower0, upper0, mean0, sigma0,
                                 1024, 16384, 8, 1e-4, 0.0, 314159, true);
    exitProbUpper[0] = 1.0 - out0.prob;
  }

  // -------- phase 2 lower exit: stop for futility if max Z <= a[0] --------
  {
    std::vector<double> lower0(M, -8.0), upper0(M, a[0]);
    PMVNResult out0 = pmvnormcpp(lower0, upper0, mean0, sigma0,
                                 1024, 16384, 8, 1e-4, 0.0, 314159, true);
    exitProbLower[0] = out0.prob;
  }

  // compute the exit probabilities at the K looks in phase 3

  // information increments for the selected arm in phase 3
  std::vector<double> I1(K), sqrtI1(K);
  for (size_t k = 0; k < K; ++k) {
    I1[k] = Ivec[k+1] - Ivec[0];
    sqrtI1[k] = std::sqrt(I1[k]);
  }


  // conditional covariance among non-selected arms in phase 2 given selected arm
  FlatMatrix sigma_m1;
  if (M > 1) {
    sigma_m1.resize(M - 1, M - 1);
    double diag = 1.0 - rho * rho;
    double off = rho * (1.0 - rho);
    for (size_t j = 0; j < M - 1; ++j) {
      double* colptr = sigma_m1.data_ptr() + j * sigma_m1.nrow;
      // fill column with off
      std::fill_n(colptr, M - 1, off);
      colptr[j] = diag;
    }
  }

  std::vector<double> b1(K);
  std::vector<double> a1(K);
  std::vector<double> theta1(K);
  std::vector<double> p1u(K), p1l(K);

  size_t m1 = (M > 1) ? M - 1 : 1;
  std::vector<double> mean(m1);
  std::vector<double> lower(m1, -8.0);
  std::vector<double> upper(m1);

  // integration intervals for z0
  std::vector<double> breaks_select = {-8.0, 8.0};
  std::vector<double> breaks_reject = {b[0], 8.0};
  std::vector<double> breaks_futility = {-8.0, a[0]};
  std::vector<double> breaks_continue = {a[0], b[0]};

  for (size_t m = 0; m < M; ++m) { // loop over the selected arm in phase 2
    double mu = theta[m] * sqrtI0;

    // compute the conditional distribution of the Wald statistics for the
    // M - 1 non-selected arms in phase 2 given that arm m is selected, and
    // then compute the probability of selecting arm m in phase 2 (i.e.,
    // the probability that the Wald statistic for arm m is larger than
    // the Wald statistics for the M - 1 non-selected arms in phase 2)
    // reuse temporary buffers for mean/lower/upper

    // density contribution for arm m being selected as best at phase 2 value z0
    auto f0 = [&theta, &sigma_m1, &mean, &lower, &upper, M, m, mu, sqrtI0, rho]
    (double z0)->double {
      if (M > 1) {
        double delta0 = z0 - mu;
        size_t j = 0;
        for (size_t i = 0; i < M; ++i) {
          if (i == m) continue;
          mean[j++] = theta[i] * sqrtI0 + rho * delta0;
        }
        std::fill_n(upper.data(), M - 1, z0);
        PMVNResult out = pmvnormcpp(lower, upper, mean, sigma_m1,
                                    1024, 16384, 8, 1e-4, 0.0, 314159, true);
        return out.prob * boost_dnorm(z0, mu, 1.0);
      } else {
        return boost_dnorm(z0, mu, 1.0);
      }
    };

    // overall selection probability
    selectAsBest[m] = integrate3(f0, breaks_select, 1e-4);

    // phase 2 efficacy and futility by arm selected
    exitProbByArmUpper(0, m) = integrate3(f0, breaks_reject, 1e-4);
    exitProbByArmLower(0, m) = integrate3(f0, breaks_futility, 1e-4);

    // compute the conditional distribution of the Wald statistics in phase 3
    // given that arm m is selected in phase 2, and then compute the probability
    // of exiting at look k in phase 3 for k = 1,...,K. Note that the conditional
    // distribution the Wald statistics for arm m in phase 3 can be written as
    // a secondary group sequential trial with K looks, where the information
    // level at look k is I1[k], the critical value at look k is b1[k], and
    // the treatment effect is theta[m]
    std::fill_n(theta1.begin(), K, theta[m]);

    auto f1 = [&b1, &a1, &b, &a, &theta1, &sqrtI, &I1, &sqrtI1, K, sqrtI0]
    (double z0, std::vector<double>& outU, std::vector<double>& outL)->void {
      for (size_t k = 0; k < K; ++k) {
        double ub = (b[k + 1] * sqrtI[k + 1] - z0 * sqrtI0) / sqrtI1[k];
        double lb = (a[k + 1] * sqrtI[k + 1] - z0 * sqrtI0) / sqrtI1[k];
        if (ub > 8.0) ub = 8.0;
        if (lb < -8.0) lb = -8.0;
        if (lb > ub) lb = ub;
        b1[k] = ub;
        a1[k] = lb;
      }

      auto probs = exitprobcpp(b1, a1, theta1, I1);
      outU = probs.exitProbUpper;
      outL = probs.exitProbLower;
    };

    // prepare memoization structures; reserve to avoid rehashing
    // memoization keyed by z0
    std::unordered_map<double, size_t> eval_idx;
    eval_idx.reserve(1024);
    std::vector<double> p0_list;
    std::vector<double> p1u_flat;
    std::vector<double> p1l_flat;
    p0_list.reserve(1024);
    p1u_flat.reserve(1024 * K);
    p1l_flat.reserve(1024 * K);

    auto eval_and_get_index = [&](double z0)->size_t {
      auto it = eval_idx.find(z0);
      if (it != eval_idx.end()) return it->second;

      double p0 = f0(z0);
      f1(z0, p1u, p1l);

      size_t idx = p0_list.size();
      eval_idx.emplace(z0, idx);
      p0_list.push_back(p0);

      p1u_flat.resize(p1u_flat.size() + K);
      p1l_flat.resize(p1l_flat.size() + K);
      for (size_t kk = 0; kk < K; ++kk) {
        p1u_flat[idx * K + kk] = p1u[kk];
        p1l_flat[idx * K + kk] = p1l[kk];
      }

      return idx;
    };

    // compute the joint distribution of the Wald statistic for arm m and the
    // Wald statistics for the M - 1 non-selected arms in phase 2, and then
    // compute the probability of selecting arm m in phase 2 and exiting at look k
    // in phase 3 for k = 1,...,K
    // now for each look k we integrate scalar function that reuses cache
    for (size_t k = 0; k < K; ++k) {
      auto hu = [&](double z0)->double {
        size_t idx = eval_and_get_index(z0);
        return p0_list[idx] * p1u_flat[idx * K + k];
      };
      auto hl = [&](double z0)->double {
        size_t idx = eval_and_get_index(z0);
        return p0_list[idx] * p1l_flat[idx * K + k];
      };

      exitProbByArmUpper(k + 1, m) = integrate3(hu, breaks_continue, 1e-4);
      exitProbByArmLower(k + 1, m) = integrate3(hl, breaks_continue, 1e-4);
      exitProbUpper[k + 1] += exitProbByArmUpper(k + 1, m);
      exitProbLower[k + 1] += exitProbByArmLower(k + 1, m);
    }
  }

  return ExitProbSeamless{
    std::move(exitProbUpper),
    std::move(exitProbLower),
    std::move(exitProbByArmUpper),
    std::move(exitProbByArmLower),
    std::move(selectAsBest)
  };
}

ExitProbSeamless exitprob_seamless_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t K,
    const std::vector<double>& b,
    const std::vector<double>& I) {

  if (!none_na(b)) throw std::invalid_argument("b must be provided");
  double amin = std::min(-8.0, *std::min_element(b.begin(), b.end()));
  std::vector<double> a(K + 1, amin);
  return exitprob_seamless_cpp(M, r, theta, corr_known, K, b, a, I);
}


//' @title Exit Probabilities for Phase 2/3 Seamless Design
//' @description Computes the upper and lower exit probabilities for a phase
//' 2/3 seamless design. In Phase 2, multiple active arms are compared
//' against a common control arm. If the phase-2 max-Z statistic crosses the
//' efficacy boundary, the trial stops early for efficacy; if it falls below
//' the futility boundary, the trial stops early for futility. Otherwise, the
//' best-performing arm is selected to proceed to Phase 3, where it is tested
//' against the common control over multiple looks with upper and optional
//' lower stopping boundaries.
//'
//' @param M Number of active treatment arms in Phase 2.
//' @param r Randomization ratio of each active arm to the common control
//'   in Phase 2.
//' @param theta A vector of length \eqn{M} representing the true treatment
//'   effects for each active arm versus the common control.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param K Number of sequential looks in Phase 3.
//' @param b A vector of efficacy boundaries (length \eqn{K+1}). The first
//'   element is the efficacy boundary for the phase-2 max-Z statistic;
//'   the remaining \eqn{K} elements are efficacy boundaries for the selected
//'   arm in Phase 3.
//' @param a An optional vector of futility boundaries (length \eqn{K+1}).
//'   The first element is the futility boundary for the phase-2 max-Z
//'   statistic; the remaining \eqn{K} elements are futility boundaries for
//'   the selected arm in Phase 3. If omitted, no futility stopping is
//'   applied.
//' @param I A vector of information levels (length \eqn{K+1}) for any active
//'   arm versus the common control. The first element is for Phase 2;
//'   the remaining \eqn{K} elements are for the looks in Phase 3.
//'
//' @details
//' The function assumes a multivariate normal distribution for the Wald
//' statistics. The "best" arm is defined as the active arm with the largest
//' Z-statistic at the end of Phase 2 among designs that continue beyond the
//' phase-2 analysis.
//'
//' \strong{Decision Rules:}
//'
//' * \strong{Phase 2 efficacy stop}: reject if the phase-2 max-Z statistic
//'   satisfies \eqn{\max_m Z_m(I_0) \ge b_0}.
//'
//' * \strong{Phase 2 futility stop}: stop for futility if the phase-2 max-Z
//'   statistic satisfies \eqn{\max_m Z_m(I_0) \le a_0}.
//'
//' * \strong{Continue to Phase 3}: if \eqn{a_0 < \max_m Z_m(I_0) < b_0},
//'   select the arm with the largest phase-2 Z-statistic and continue with
//'   that arm only.
//'
//' * \strong{Phase 3 efficacy stop}: at look \eqn{k}, reject if the selected
//'   arm's Z-statistic exceeds the efficacy boundary and no earlier stop has
//'   occurred.
//'
//' * \strong{Phase 3 futility stop}: at look \eqn{k}, stop for futility if
//'   the selected arm's Z-statistic is below the futility boundary and no
//'   earlier stop has occurred.
//'
//' \strong{Design Assumptions:}
//'
//' * All active arms share the same information level in Phase 2.
//'
//' * Exactly one active arm is selected at the end of Phase 2 based on the
//'   largest observed Z-statistic when the trial continues to Phase 3.
//'
//' @return A list containing the following components:
//'
//' * \code{exitProbUpper}: A vector of length \eqn{K + 1}. The first element
//'   is the probability of stopping for efficacy in Phase 2; the remaining
//'   elements are the probabilities of stopping for efficacy at each look in
//'   Phase 3.
//'
//' * \code{exitProbLower}: A vector of length \eqn{K + 1}. The first element
//'   is the probability of stopping for futility in Phase 2; the remaining
//'   elements are the probabilities of stopping for futility at each look in
//'   Phase 3.
//'
//' * \code{exitProbByArmUpper}: A \eqn{(K+1) \times M} matrix. The
//'   \eqn{(k, m)}-th entry gives the probability of stopping for efficacy at
//'   look \eqn{k} given that arm \eqn{m} is selected as best.
//'
//' * \code{exitProbByArmLower}: A \eqn{(K+1) \times M} matrix. The
//'   \eqn{(k, m)}-th entry gives the probability of stopping for futility at
//'   look \eqn{k} given that arm \eqn{m} is selected as best.
//'
//' * \code{selectAsBest}: A vector of length \eqn{M} containing the
//'   probability that each active arm is selected to move on to Phase 3.
//'
//' For backward compatibility, the list also contains:
//'
//' * \code{exitProb}: identical to \code{exitProbUpper}.
//'
//' * \code{exitProbByArm}: identical to \code{exitProbByArmUpper}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive two-stage seamless sequential design for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
//'
//' @examples
//'
//' # Setup: 2 active arms vs control in Phase 2; 1 selected arm vs control
//' # in Phase 3. Phase 3 has 2 sequential looks.
//'
//' # Information levels: equal spacing over 3 looks based on a maximum of
//' # 110 patients per arm, SD = 1.0
//' I <- c(110 / (2 * 1.0^2) * seq(1, 3)/3)
//'
//' # O'Brien-Fleming efficacy boundaries
//' b <- c(3.776605, 2.670463, 2.180424)
//'
//' # No futility stopping
//' p0 <- exitprob_seamless(M = 2, theta = c(0, 0), K = 2, b = b, I = I)
//' cumsum(p0$exitProbUpper)
//'
//' # Add futility stopping
//' a <- c(0, 0.5, b[3])
//' p1 <- exitprob_seamless(M = 2, theta = c(0.3, 0.5), K = 2, b = b, a = a, I = I)
//' cbind(
//'   cumulativeEfficacy = cumsum(p1$exitProbUpper),
//'   cumulativeFutility = cumsum(p1$exitProbLower)
//' )
//'
//' @export
// [[Rcpp::export]]
Rcpp::List exitprob_seamless(
    const int M = NA_INTEGER,
    const double r = 1,
    const Rcpp::NumericVector& theta = NA_REAL,
    const bool corr_known = true,
    const int K = NA_INTEGER,
    SEXP b = R_NilValue,
    SEXP a = R_NilValue,
    SEXP I = R_NilValue) {

  if (M == NA_INTEGER || M < 1) {
    throw std::invalid_argument("M must be a positive integer");
  }
  if (K == NA_INTEGER || K < 1) {
    throw std::invalid_argument("K must be a positive integer");
  }

  size_t M_ = static_cast<size_t>(M);
  size_t K_ = static_cast<size_t>(K);

  std::vector<double> thetaVec(theta.begin(), theta.end());

  std::vector<double> bVec;
  if (b == R_NilValue) {
    throw std::invalid_argument("b must be provided");
  } else {
    Rcpp::NumericVector bv(b);
    bVec.assign(bv.begin(), bv.end());
  }

  std::vector<double> aVec;
  if (a == R_NilValue) {
    double amin = -8.0;
    if (!bVec.empty()) {
      amin = std::min(-8.0, *std::min_element(bVec.begin(), bVec.end()));
    }
    aVec.assign(K_ + 1, amin);
  } else {
    Rcpp::NumericVector av(a);
    aVec.assign(av.begin(), av.end());
  }

  std::vector<double> IVec;
  if (I == R_NilValue) {
    IVec = std::vector<double>(1, NA_REAL);
  } else {
    Rcpp::NumericVector Iv(I);
    IVec.assign(Iv.begin(), Iv.end());
  }

  auto result = exitprob_seamless_cpp(
    M_, r, thetaVec, corr_known, K_, bVec, aVec, IVec);

  ListCpp out;
  out.push_back(result.exitProbUpper, "exitProbUpper");
  out.push_back(result.exitProbLower, "exitProbLower");
  out.push_back(result.exitProbByArmUpper, "exitProbByArmUpper");
  out.push_back(result.exitProbByArmLower, "exitProbByArmLower");
  out.push_back(result.selectAsBest, "selectAsBest");
  return Rcpp::wrap(out);
}


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
    const std::vector<unsigned char>& efficacyStopping) {

  if (M < 1) {
    throw std::invalid_argument("M should be at least 1");
  }
  if (r <= 0.0) {
    throw std::invalid_argument("r should be positive");
  }
  if (k < 1) {
    throw std::invalid_argument("k should be at least 1");
  }

  size_t K = k + 1; // add the look at end of phase 2

  // infoRates: if missing create 1/K, 2/K, ..., k/K
  size_t kMax = K;
  std::vector<double> infoRates(K);
  if (none_na(informationRates)) {
    kMax = informationRates.size();
    if (kMax < K)
      throw std::invalid_argument("Insufficient length for informationRates");

    infoRates = informationRates;
    if (infoRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(infoRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (infoRates.back() > 1.0)
      throw std::invalid_argument("informationRates must not exceed 1");
  } else {
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


  ExitProbSeamless probs;
  std::vector<double> criticalValues(kMax);
  std::vector<double> zero(M, 0.0);

  if (asf == "none") {
    for (size_t i = 0; i < kMax-1; ++i) criticalValues[i] = 8.0;

    auto f = [&](double aval)->double {
      criticalValues[kMax-1] = aval;
      probs = exitprob_seamless_cpp(M, r, zero, corr_known, kMax - 1,
                                    criticalValues, infoRates);
      double cpu = std::accumulate(probs.exitProbUpper.begin(),
                                   probs.exitProbUpper.end(), 0.0);
      return cpu - alpha;
    };

    criticalValues[kMax-1] = brent(f, 0.0, 8.0, 1e-6);
    return subset(criticalValues, 0, K);
  }

  if (asf == "of" || asf == "p" || asf == "wt") {
    double Delta;
    if (asf == "of") Delta = 0.0;
    else if (asf == "p") Delta = 0.5;
    else Delta = parameterAlphaSpending; // parameterAlphaSpending holds delta for WT

    // for a given multiplier, compute cumulative upper exit probability - alpha
    std::vector<double> u(kMax);
    std::vector<double> u0(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      u0[i] = std::pow(infoRates[i], Delta - 0.5);
    }

    auto f = [&](double aval)->double {
      for (size_t i = 0; i < kMax; ++i) {
        u[i] = aval * u0[i];
        if (!effStopping[i]) u[i] = 8.0;
      }

      probs = exitprob_seamless_cpp(M, r, zero, corr_known, kMax - 1,
                                    u, infoRates);
      double cpu = std::accumulate(probs.exitProbUpper.begin(),
                                   probs.exitProbUpper.end(), 0.0);
      return cpu - alpha;
    };

    double cwt = brent(f, 0.0, 10.0, 1e-6);
    for (size_t i = 0; i < kMax; ++i) {
      double val = cwt * u0[i];
      criticalValues[i] = effStopping[i] ? val : 8.0;
    }
    return subset(criticalValues, 0, K);
  }

  if (asf == "sfof" || asf == "sfp" || asf == "sfkd" ||
      asf == "sfhsd" || asf == "user") {
    // stage 1
    double cumAlpha;
    if (asf == "user") cumAlpha = userAlphaSpending[0];
    else cumAlpha = errorSpentcpp(spendTime[0], alpha, asf, parameterAlphaSpending);

    if (!effStopping[0]) criticalValues[0] = 8.0;
    else {
      FlatMatrix sigma(M, M);
      double rho = corr_known ? r / (r + 1.0) : 0;
      sigma.fill(rho);
      for (size_t i = 0; i < M; ++i) {
        sigma(i, i) = 1.0;
      }
      criticalValues[0] = qmvnormcpp(1.0 - cumAlpha, zero, sigma,
                                     1024, 16384, 8, 1e-4, 0.0, 314159, true);
    }

    // subsequent stages
    for (size_t k1 = 1; k1 < kMax; ++k1) {
      // determine cumulative alpha at this stage
      if (asf == "user") cumAlpha = userAlphaSpending[k1];
      else cumAlpha = errorSpentcpp(spendTime[k1], alpha, asf,
                                    parameterAlphaSpending);

      if (!effStopping[k1]) {
        criticalValues[k1] = 8.0;
        continue;
      }

      // Define lambda that only sets the last element of u_vec
      auto f = [&](double aval)->double {
        // set the last element to the current candidate critical value
        criticalValues[k1] = aval;
        probs = exitprob_seamless_cpp(
          M, r, zero, corr_known, k1, criticalValues, infoRates);
        double cpu = std::accumulate(probs.exitProbUpper.begin(),
                                     probs.exitProbUpper.end(), 0.0);
        return cpu - cumAlpha;
      };

      double f_8 = f(8.0);
      if (f_8 > 0.0) { // no alpha spent at current visit
        criticalValues[k1] = 8.0;
      } else {
        auto f_for_brent = [&](double x)->double {
          if (x == 8.0) return f_8; // avoid recomputation at 8.0
          return f(x);
        };

        criticalValues[k1] = brent(f_for_brent, 0.0, 8.0, 1e-6);
      }
    }

    return subset(criticalValues, 0, K);
  }

  throw std::invalid_argument("Invalid value for typeAlphaSpending");
}


//' @title Efficacy Boundaries for Phase 2/3 Seamless Design
//' @description Calculates the efficacy stopping boundaries for a phase 2/3
//' seamless design, accounting for the selection of the best arm
//' at the end of Phase 2 and sequential testing in Phase 3.
//'
//' @param M Number of active treatment arms in Phase 2.
//' @param r Randomization ratio of each active arm to the common control
//'   in Phase 2.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is assumed.
//' @param k The index of the current look in Phase 3.
//' @param informationRates A numeric vector of information rates up to the
//'   current look. Values must be strictly increasing and \eqn{\le 1}.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A numeric vector of length \eqn{k+1} specifying the
//'   error spending time at each analysis. Values must be strictly increasing
//'   and \eqn{\le 1}. If omitted, defaults to \code{informationRates}.
//' @inheritParams param_efficacyStopping
//'
//' @details
//' The function determines critical values by solving for the boundary that
//' satisfies the alpha-spending requirement, given the selection of the
//' "best" arm at the end of Phase 2.
//'
//' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
//' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
//' and \code{spendingTime} must be of full length \code{kMax}, and
//' \code{informationRates} and \code{spendingTime} must end with 1.
//'
//' @return A numeric vector of length \eqn{k + 1} containing the critical
//' values (on the standard normal Z-scale) for each analysis up to the
//' current look.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive two-stage seamless sequential design for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
//'
//' @examples
//'
//' # Determine O'Brien-Fleming boundaries for a seamless design with
//' # 2 active arms in Phase 2 and 2 looks in Phase 3 (3 looks total).
//' getBound_seamless(M = 2, k = 2, informationRates = seq(1, 3)/3,
//'                   alpha = 0.025, typeAlphaSpending = "OF")
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector getBound_seamless(
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
  auto effStopping = convertLogicalVector(efficacyStopping);
  std::vector<double> userAlpha(userAlphaSpending.begin(), userAlphaSpending.end());
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());

  auto result = getBound_seamless_cpp(
    static_cast<size_t>(M), r, corr_known, static_cast<size_t>(k),
    infoRates, alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlpha, spendTime, effStopping
  );

  return Rcpp::wrap(result);
}


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
    const std::vector<unsigned char>& futStopping) {

  if (M < 1) {
    throw std::invalid_argument("M should be at least 1");
  }
  if (r <= 0.0) {
    throw std::invalid_argument("r should be positive");
  }
  if (theta.size() != M) {
    throw std::invalid_argument("theta should have length M");
  }
  if (alpha < 0.00001 || alpha >= 1) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (K < 1) {
    throw std::invalid_argument("K should be at least 1");
  }
  if (critValues.size() < K + 1) {
    throw std::invalid_argument("Insufficient length for critValues");
  }
  if (I.size() < K + 1) {
    throw std::invalid_argument("Insufficient length for I");
  }
  if (st.size() < K + 1) {
    throw std::invalid_argument("Insufficient length for st");
  }
  if (futStopping.size() < K + 1) {
    throw std::invalid_argument("Insufficient length for futStopping");
  }
  if (I[0] <= 0.0) {
    throw std::invalid_argument("I must be positive");
  }
  if (any_nonincreasing(I)) {
    throw std::invalid_argument("I must be increasing");
  }


  double sqrtI0 = std::sqrt(I[0]);
  std::vector<double> lo(M, -8.0), hi(M, 8.0);
  std::vector<double> mean0(M);
  for (size_t m = 0; m < M; ++m) {
    mean0[m] = theta[m] * sqrtI0;
  }

  FlatMatrix sigma0(M, M);
  sigma0.fill(r / (r + 1.0));
  for (size_t i = 0; i < M; ++i) {
    sigma0(i, i) = 1.0;
  }

  ExitProbSeamless probs;
  std::vector<double> futBounds(K + 1, -8.0);
  auto f = [&](double beta)->double {
    std::fill(futBounds.begin(), futBounds.end(), -8.0);

    double eps = 0.0;
    double cb = 0.0;

    // stage 0: futility on phase-2 max-Z scale
    if (futStopping[0]) {
      cb = errorSpentcpp(st[0], beta, bsf, bsfpar);

      auto g = [&](double aval)->double {
        std::fill(hi.begin(), hi.end(), aval);
        PMVNResult out0 = pmvnormcpp(lo, hi, mean0, sigma0,
                                     1024, 16384, 8, 1e-4, 0.0, 314159, true);
        return out0.prob - cb;
      };

      eps = g(critValues[0]);
      if (eps < 0.0) return -1.0;
      futBounds[0] = brent(g, -8.0, critValues[0], 1e-6);
    }

    // stages 1..K
    for (size_t k = 1; k < K + 1; ++k) {
      if (!futStopping[k]) continue;

      cb = errorSpentcpp(st[k], beta, bsf, bsfpar);

      auto g = [&](double aval)->double {
        futBounds[k] = aval;
        probs = exitprob_seamless_cpp(M, r, theta, true, k, critValues,
                                      futBounds, I);
        double cpl = std::accumulate(probs.exitProbLower.begin(),
                                     probs.exitProbLower.end(), 0.0);
        return cpl - cb;
      };

      double bk = critValues[k];
      eps = g(bk);
      double g_minus8 = g(-8.0);

      if (g_minus8 > 0.0) {
        futBounds[k] = -8.0;
      } else if (eps > 0.0) {
        auto g_for_brent = [&](double aval)->double {
          if (aval == -8.0) return g_minus8;
          if (aval == bk) return eps;
          return g(aval);
        };
        futBounds[k] = brent(g_for_brent, -8.0, bk, 1e-6);
      } else if (k < K) {
        return -1.0;
      }
    }

    return eps;
  };

  double v1 = f(0.0001);
  double v2 = f(1.0 - alpha);
  double beta = 0.0;

  if (v1 == -1.0 || (v1 < 0.0 && futBounds[K] == -8.0)) {
    throw std::invalid_argument("Power must be less than 0.9999");
  } else if (v2 > 0.0) {
    throw std::invalid_argument("Power must be greater than alpha");
  } else {
    auto f_for_brent = [&](double x)->double {
      if (x == 0.0001) return v1;
      if (x == 1.0 - alpha) return v2;
      return f(x);
    };

    beta = brent(f_for_brent, 0.0001, 1.0 - alpha, 1e-6);
    futBounds[K] = critValues[K];
    probs = exitprob_seamless_cpp(M, r, theta, true, K, critValues, futBounds, I);
  }

  return GetPowerSeamless{
    1.0 - beta,
    std::move(futBounds),
    std::move(probs)
  };
}



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
    const std::vector<double>& spendingTime) {

  // ----------- Input Validation ----------- //
  if (std::isnan(beta) && std::isnan(IMax)) {
    throw std::invalid_argument("beta and IMax cannot be missing simultaneously");
  }
  if (!std::isnan(beta) && !std::isnan(IMax)) {
    throw std::invalid_argument("Only one of beta and IMax should be provided");
  }
  if (!std::isnan(IMax) && IMax <= 0.0) {
    throw std::invalid_argument("IMax must be positive");
  }
  if (M < 1) {
    throw std::invalid_argument("M must be at least 1");
  }
  if (r <= 0.0) {
    throw std::invalid_argument("r must be positive");
  }
  if (!none_na(theta)) {
    throw std::invalid_argument("theta must be provided");
  }
  if (theta.size() != M) {
    throw std::invalid_argument("theta should have length M");
  }
  if (K < 1) {
    throw std::invalid_argument("K must be at least 1");
  }

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1.0)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(beta) && (beta >= 1.0 - alpha || beta < 0.0001)) {
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");
  }

  std::string unknown = std::isnan(beta) ? "beta" : "IMax";
  size_t kMax = K + 1;

  // informationRates
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (informationRates.size() != kMax) {
      throw std::invalid_argument("Invalid length for informationRates");
    }
    if (informationRates[0] <= 0.0) {
      throw std::invalid_argument("informationRates must be positive");
    }
    if (any_nonincreasing(informationRates)) {
      throw std::invalid_argument("informationRates must be increasing");
    }
    if (informationRates.back() != 1.0) {
      throw std::invalid_argument("informationRates must end with 1");
    }
    infoRates = informationRates;
  } else {
    for (size_t i = 0; i < kMax; ++i) {
      infoRates[i] = static_cast<double>(i + 1) / static_cast<double>(kMax);
    }
  }

  // efficacyStopping
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != kMax) {
      throw std::invalid_argument("Invalid length for efficacyStopping");
    }
    if (efficacyStopping.back() != 1) {
      throw std::invalid_argument("efficacyStopping must end with 1");
    }
    effStopping = efficacyStopping;
  } else {
    effStopping.assign(kMax, 1);
  }

  // futilityStopping
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != kMax) {
      throw std::invalid_argument("Invalid length for futilityStopping");
    }
    if (futilityStopping.back() != 1) {
      throw std::invalid_argument("futilityStopping must end with 1");
    }
    futStopping = futilityStopping;
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds =
      !none_na(futilityBounds) && !none_na(futilityCP) && !none_na(futilityTheta);

  if (!missingCriticalValues && criticalValues.size() != kMax) {
    throw std::invalid_argument("Invalid length for criticalValues");
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }

  std::string asf = typeAlphaSpending;
  for (char& c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues &&
      !(asf == "of" || asf == "p" || asf == "wt" ||
        asf == "sfof" || asf == "sfp" || asf == "sfkd" ||
        asf == "sfhsd" || asf == "user" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument("parameterAlphaSpending must be positive for sfKD");
  }

  if (missingCriticalValues && asf == "user") {
    if (!none_na(userAlphaSpending)) {
      throw std::invalid_argument("userAlphaSpending must be specified");
    }
    if (userAlphaSpending.size() != kMax) {
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    }
    if (userAlphaSpending[0] < 0.0) {
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    }
    if (any_nonincreasing(userAlphaSpending)) {
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    }
    if (userAlphaSpending.back() != alpha) {
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
    }
  }

  if (!missingFutilityBounds) {
    if (none_na(futilityBounds)) {
      if (futilityBounds.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityBounds");
      }
    } else if (none_na(futilityCP)) {
      if (futilityCP.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityCP");
      }
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (futilityCP[i] < 0.0 || futilityCP[i] > 1.0) {
          throw std::invalid_argument("futilityCP must lie in [0, 1]");
        }
      }
    } else {
      if (futilityTheta.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityTheta");
      }
    }
  }

  std::string bsf = typeBetaSpending;
  for (char& c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (unknown == "IMax") {
    if (missingFutilityBounds &&
        !(bsf == "sfof" || bsf == "sfp" || bsf == "sfkd" ||
          bsf == "sfhsd" || bsf == "user" || bsf == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpending");
    }
  } else {
    if (missingFutilityBounds &&
        !(bsf == "sfof" || bsf == "sfp" || bsf == "sfkd" ||
          bsf == "sfhsd" || bsf == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpending");
    }
  }

  if ((bsf == "sfkd" || bsf == "sfhsd") &&
      std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument("parameterBetaSpending must be positive for sfKD");
  }

  if (unknown == "IMax" && bsf == "user") {
    if (!none_na(userBetaSpending)) {
      throw std::invalid_argument("userBetaSpending must be specified");
    }
    if (userBetaSpending.size() != kMax) {
      throw std::invalid_argument("Invalid length of userBetaSpending");
    }
    if (userBetaSpending[0] < 0.0) {
      throw std::invalid_argument("userBetaSpending must be nonnegative");
    }
    if (any_nonincreasing(userBetaSpending)) {
      throw std::invalid_argument("userBetaSpending must be nondecreasing");
    }
    if (userBetaSpending.back() != beta) {
      throw std::invalid_argument("userBetaSpending must end with specified beta");
    }
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax) {
      throw std::invalid_argument("Invalid length for spendingTime");
    }
    if (spendingTime[0] <= 0.0) {
      throw std::invalid_argument("spendingTime must be positive");
    }
    if (any_nonincreasing(spendingTime)) {
      throw std::invalid_argument("spendingTime must be increasing");
    }
    if (spendingTime.back() != 1.0) {
      throw std::invalid_argument("spendingTime must end with 1");
    }
    spendTime = spendingTime;
  } else {
    spendTime = infoRates;
  }

  // ----------- End of Input Validation ----------- //

  // set up efficacy bounds
  ExitProbSeamless probs;
  std::vector<double> zero(M, 0.0);
  std::vector<double> critValues = criticalValues;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (kMax > 1 && criticalValues.size() == kMax) {
      bool hasNaN = false;
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (std::isnan(criticalValues[i])) {
          hasNaN = true;
          break;
        }
      }
      if (!hasNaN && std::isnan(criticalValues[kMax - 1])) haybittle = true;
    }

    if (haybittle) {
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (!effStopping[i]) critValues[i] = 8.0;
      }

      auto f = [&](double aval)->double {
        critValues[kMax - 1] = aval;
        probs = exitprob_seamless_cpp(
          M, r, zero, corr_known, K, critValues, infoRates);
        double cpu = std::accumulate(probs.exitProbUpper.begin(),
                                     probs.exitProbUpper.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax - 1] = brent(f, 0.0, 8.0, 1e-6);
    } else {
      critValues = getBound_seamless_cpp(
        M, r, corr_known, K, infoRates, alpha, asf, parameterAlphaSpending,
        userAlphaSpending, spendTime, effStopping);
    }
  }

  probs = exitprob_seamless_cpp(M, r, zero, corr_known, K, critValues, infoRates);
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(probs.exitProbUpper.begin(),
                   probs.exitProbUpper.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent.back();

  // set up futility bounds
  std::vector<double> futBounds(kMax, NaN);

  if (kMax > 1) {
    if (missingFutilityBounds && bsf == "none") {
      std::fill_n(futBounds.begin(), kMax - 1, -8.0);
      futBounds[kMax - 1] = critValues[kMax - 1];
    } else if (!missingFutilityBounds) {
      if (none_na(futilityBounds)) {
        for (size_t i = 0; i < kMax - 1; ++i) {
          if (futilityBounds[i] > critValues[i]) {
            throw std::invalid_argument(
                "futilityBounds must lie below criticalValues");
          }
        }
        std::copy_n(futilityBounds.begin(), kMax - 1, futBounds.begin());
        futBounds[kMax - 1] = critValues[kMax - 1];
      } else if (none_na(futilityCP)) {
        double c2 = critValues[kMax - 1];
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i]) *
            (c2 - std::sqrt(1 - infoRates[i]) * boost_qnorm(1 - futilityCP[i]));
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityCP values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax - 1] = critValues[kMax - 1];
      } else if (!std::isnan(IMax)) {
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i] * IMax) * futilityTheta[i];
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityTheta values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax - 1] = critValues[kMax - 1];
      }
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  double IMax1 = IMax;
  std::vector<double> information(kMax);

  if (unknown == "IMax") {
    double maxtheta = *std::max_element(theta.begin(), theta.end());
    std::vector<double> lo(M, -8.0), hi(M, critValues[0]);
    std::vector<double> mu0(M);
    FlatMatrix sigma0(M, M);
    sigma0.fill(r / (r + 1.0));
    for (size_t i = 0; i < M; ++i) {
      sigma0(i, i) = 1.0;
    }

    auto f = [&](double x)->double {
      double maxInformation = sq(x / maxtheta);
      for (size_t i = 0; i < kMax; ++i) {
        information[i] = infoRates[i] * maxInformation;
      }
      double sqrtI0 = std::sqrt(information[0]);

      if (!missingFutilityBounds || bsf == "none" || kMax == 1) {
        if (!none_na(futilityBounds) && !none_na(futilityCP) &&
            none_na(futilityTheta)) {
          for (size_t i = 0; i < kMax - 1; ++i) {
            futBounds[i] = std::sqrt(information[i]) * futilityTheta[i];
            if (futBounds[i] > critValues[i]) return -1.0;
          }
          futBounds[kMax - 1] = critValues[kMax - 1];
        }

        probs = exitprob_seamless_cpp(
          M, r, theta, true, K, critValues, futBounds, information);
        double overallReject = std::accumulate(probs.exitProbUpper.begin(),
                                               probs.exitProbUpper.end(), 0.0);
        return (1.0 - overallReject) - beta;
      } else {
        std::fill(futBounds.begin(), futBounds.end(), -8.0);

        double eps = 0.0, cb = 0.0;

        // stage 0: futility on phase-2 max-Z scale
        if (futStopping[0]) {
          cb = (bsf == "user") ? userBetaSpending[0] :
          errorSpentcpp(spendTime[0], beta, bsf, parameterBetaSpending);

          for (size_t m = 0; m < M; ++m) {
            mu0[m] = theta[m] * sqrtI0;
          }

          auto q = pmvnormcpp(lo, hi, mu0, sigma0,
                              1024, 16384, 8, 1e-4, 0.0, 314159);
          eps = q.prob - cb;
          if (eps < 0.0) return -1.0;
          futBounds[0] = qmvnormcpp(cb, mu0, sigma0,
                                    1024, 16384, 8, 1e-4, 0.0, 314159);
        }

        // stages 1..K
        for (size_t k = 1; k < K + 1; ++k) {
          if (!futStopping[k]) continue;

          cb = (bsf == "user") ? userBetaSpending[k] :
            errorSpentcpp(spendTime[k], beta, bsf, parameterBetaSpending);

          auto g = [&](double aval)->double {
            futBounds[k] = aval;
            probs = exitprob_seamless_cpp(M, r, theta, true, k, critValues,
                                          futBounds, information);
            double cpl = std::accumulate(probs.exitProbLower.begin(),
                                         probs.exitProbLower.end(), 0.0);
            return cpl - cb;
          };

          double bk = critValues[k];
          eps = g(bk);
          double g_minus8 = g(-8.0);

          if (g_minus8 > 0.0) {
            futBounds[k] = -8.0;
          } else if (eps > 0.0) {
            auto g_for_brent = [&](double aval)->double {
              if (aval == -8.0) return g_minus8;
              if (aval == bk) return eps;
              return g(aval);
            };
            futBounds[k] = brent(g_for_brent, -8.0, bk, 1e-6);
          } else if (k < K) {
            return -1.0;
          }
        }

        return eps;
      }
    };

    double drift = brent(f, 0.001, 8.0, 1e-6);
    IMax1 = sq(drift / maxtheta);
    futBounds[kMax-1] = critValues[kMax-1];
    probs = exitprob_seamless_cpp(M, r, theta, true, K, critValues,
                                  futBounds, information);
  } else {
    for (size_t i = 0; i < kMax; ++i) information[i] = infoRates[i] * IMax1;

    if (!missingFutilityBounds || bsf == "none" || kMax == 1) {
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (std::isnan(futBounds[i])) futBounds[i] = -8.0;
      }
      futBounds[kMax - 1] = critValues[kMax - 1];
      probs = exitprob_seamless_cpp(M, r, theta, true, K, critValues,
                                    futBounds, information);
    } else {
      auto out = getPower_seamless(
        M, r, theta, alpha1, K, critValues, information, bsf,
        parameterBetaSpending, spendTime, futStopping);
      futBounds = out.futilityBounds;
      probs = out.probs;
    }
  }

  // output quantities
  std::vector<double> efficacyTheta(kMax);
  std::vector<double> futilityThetaOut(kMax);
  std::vector<double> efficacyP(kMax);
  std::vector<double> futilityP(kMax);

  for (size_t i = 0; i < kMax; ++i) {
    efficacyTheta[i] = critValues[i] / std::sqrt(information[i]);
    futilityThetaOut[i] = futBounds[i] / std::sqrt(information[i]);
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
    futilityP[i] = 1.0 - boost_pnorm(futBounds[i]);
  }

  auto pu = probs.exitProbUpper;
  auto pl = probs.exitProbLower;

  std::vector<double> cpu(kMax), cpl(kMax);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());

  double overallReject = cpu.back();

  std::vector<double> ptotal(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    ptotal[i] = pu[i] + pl[i];
  }

  std::vector<double> informationOverall(kMax);
  double informationDelta = (M - 1) * r / (r + 1.0) * information[0];
  for (size_t i = 0; i < kMax; ++i) {
    informationOverall[i] = information[i] + informationDelta;
  }
  double IMaxOverall = informationOverall.back();

  double expectedInformationH1 = std::inner_product(
      ptotal.begin(), ptotal.end(), information.begin(), 0.0);
  double expectedInformationOverallH1 = std::inner_product(
    ptotal.begin(), ptotal.end(), informationOverall.begin(), 0.0);

  auto probsH0 = exitprob_seamless_cpp(M, r, zero, true, K, critValues,
                                       futBounds, infoRates);
  auto puH0 = probsH0.exitProbUpper;
  auto plH0 = probsH0.exitProbLower;

  std::vector<double> cpuH0(kMax), cplH0(kMax);
  std::partial_sum(puH0.begin(), puH0.end(), cpuH0.begin());
  std::partial_sum(plH0.begin(), plH0.end(), cplH0.begin());

  double overallRejectH0 = cpuH0.back();

  std::vector<double> ptotalH0(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    ptotalH0[i] = puH0[i] + plH0[i];
  }
  double expectedInformationH0 = std::inner_product(
      ptotalH0.begin(), ptotalH0.end(), information.begin(), 0.0);
  double expectedInformationOverallH0 = std::inner_product(
    ptotalH0.begin(), ptotalH0.end(), informationOverall.begin(), 0.0);

  for (size_t i = 0; i < kMax; ++i) {
    if (critValues[i] == 8.0) effStopping[i] = 0;
    if (futBounds[i] == -8.0) futStopping[i] = 0;
  }

  auto selectAsBest = probs.selectAsBest;
  auto exitProbByArmUpper = probs.exitProbByArmUpper;
  std::vector<double> powerByArm(M, 0.0);
  for (size_t j = 0; j < M; ++j) {
    double p = 0.0;
    double* colptr = exitProbByArmUpper.data_ptr() + j * exitProbByArmUpper.nrow;
    for (size_t i = 0; i < kMax; ++i) p += colptr[i];
    powerByArm[j] = p;
  }

  std::vector<double> condPowerByArm(M, 0.0);
  for (size_t j = 0; j < M; ++j) {
    if (selectAsBest[j] > 0.0) {
      condPowerByArm[j] = powerByArm[j] / selectAsBest[j];
    }
  }

  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha1, "alpha");
  overallResults.push_back(overallRejectH0, "attainedAlpha");
  overallResults.push_back(M, "M");
  overallResults.push_back(r, "r");
  overallResults.push_back(corr_known, "corr_known");
  overallResults.push_back(K, "K");
  overallResults.push_back(IMax1, "information");
  overallResults.push_back(expectedInformationH1, "expectedInformationH1");
  overallResults.push_back(expectedInformationH0, "expectedInformationH0");
  overallResults.push_back(IMaxOverall, "informationOverall");
  overallResults.push_back(expectedInformationOverallH1,
                           "expectedInformationOverallH1");
  overallResults.push_back(expectedInformationOverallH0,
                           "expectedInformationOverallH0");

  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(futBounds), "futilityBounds");
  byStageResults.push_back(std::move(pu), "rejectPerStage");
  byStageResults.push_back(std::move(pl), "futilityPerStage");
  byStageResults.push_back(std::move(cpu), "cumulativeRejection");
  byStageResults.push_back(std::move(cpl), "cumulativeFutility");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(efficacyTheta), "efficacyTheta");
  byStageResults.push_back(std::move(futilityThetaOut), "futilityTheta");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(futilityP), "futilityP");
  byStageResults.push_back(std::move(information), "information");
  byStageResults.push_back(std::move(informationOverall), "informationOverall");
  byStageResults.push_back(std::move(effStopping), "efficacyStopping");
  byStageResults.push_back(std::move(futStopping), "futilityStopping");
  byStageResults.push_back(std::move(puH0), "rejectPerStageH0");
  byStageResults.push_back(std::move(plH0), "futilityPerStageH0");
  byStageResults.push_back(std::move(cpuH0), "cumulativeRejectionH0");
  byStageResults.push_back(std::move(cplH0), "cumulativeFutilityH0");

  DataFrameCpp byArmResults;
  byArmResults.push_back(theta, "theta");
  byArmResults.push_back(std::move(selectAsBest), "selectAsBest");
  byArmResults.push_back(std::move(powerByArm), "powerByArm");
  byArmResults.push_back(std::move(condPowerByArm), "condPowerByArm");

  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(typeBetaSpending, "typeBetaSpending");
  settings.push_back(parameterBetaSpending, "parameterBetaSpending");
  settings.push_back(userBetaSpending, "userBetaSpending");
  settings.push_back(spendingTime, "spendingTime");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(byArmResults), "byArmResults");
  result.push_back(std::move(settings), "settings");
  return result;
}


//' @title Power and Sample Size for Phase 2/3 Seamless Design
//' @description Computes either the maximum information and stopping
//' boundaries for a phase 2/3 seamless design, or the achieved power when
//' the maximum information and stopping boundaries are provided. Both
//' efficacy and futility stopping can be incorporated.
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
//' @param M Number of active treatment arms in Phase 2.
//' @param r Randomization ratio of each active arm to the common control
//'   in Phase 2.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param K Number of sequential looks in Phase 3.
//' @param informationRates A numeric vector of information rates fixed
//'   before the trial. If unspecified, defaults to \eqn{(1:(K+1)) / (K+1)}.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @param criticalValues The upper boundaries on the max-Z statistic scale
//'   for Phase 2 and the Z statistics for the selected arm in Phase 3.
//'   If missing, boundaries will be computed based on the specified alpha
//'   spending function.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param futilityBounds A numeric vector of length \eqn{K} specifying
//'   futility boundaries on the max-Z scale at the end of Phase 2 and
//'   on the Z scale for the \eqn{K - 1} analyses in Phase 3. The final analysis
//'   uses the efficacy boundary as the futility boundary.
//' @param futilityCP A numeric vector of length \eqn{K} specifying futility
//'   boundaries on the conditional power scale.
//' @param futilityTheta A numeric vector of length \eqn{K} specifying futility
//'   boundaries on the parameter scale.
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param spendingTime A numeric vector of length \eqn{K+1} specifying the
//'   error spending time at each analysis. Values must be strictly increasing
//'   and end at 1. If omitted, defaults to \code{informationRates}.
//'
//' @return An S3 object of class \code{seamless} with the following components:
//'
//' * \code{overallResults}: A data frame containing:
//'     - \code{overallReject}: Overall probability of rejecting the null
//'       hypothesis.
//'     - \code{alpha}: Overall significance level.
//'     - \code{attainedAlpha}: The attained significance level, which may
//'       differ from \code{alpha} in the presence of futility stopping.
//'     - \code{M}: Number of active arms in Phase 2.
//'     - \code{r}: Randomization ratio per active arm versus control in
//'       Phase 2.
//'     - \code{corr_known}: Whether the phase-2 correlation was assumed known.
//'     - \code{K}: Number of looks in Phase 3.
//'     - \code{information}: Maximum information for any active arm versus
//'       control.
//'     - \code{expectedInformationH1}: Expected information under the
//'       alternative.
//'     - \code{expectedInformationH0}: Expected information under the null.
//'     - \code{informationOverall}: Maximum information for the overall study.
//'     - \code{expectedInformationH1}: Expected information under the
//'       alternative for the overall study.
//'     - \code{expectedInformationH0}: Expected information under the null
//'       for the overall study.
//'
//' * \code{byStageResults}: A data frame containing:
//'     - \code{informationRates}: Information rates at each analysis.
//'     - \code{efficacyBounds}: Efficacy boundaries on the Z scale.
//'     - \code{futilityBounds}: Futility boundaries on the Z scale.
//'     - \code{rejectPerStage}: Probability of efficacy stopping at each stage.
//'     - \code{futilityPerStage}: Probability of futility stopping at each
//'       stage.
//'     - \code{cumulativeRejection}: Cumulative probability of efficacy
//'       stopping.
//'     - \code{cumulativeFutility}: Cumulative probability of futility
//'       stopping.
//'     - \code{cumulativeAlphaSpent}: Cumulative alpha spent.
//'     - \code{efficacyTheta}: Efficacy boundaries on the parameter scale.
//'     - \code{futilityTheta}: Futility boundaries on the parameter scale.
//'     - \code{efficacyP}: Efficacy boundaries on the p-value scale.
//'     - \code{futilityP}: Futility boundaries on the p-value scale.
//'     - \code{information}: Cumulative information at each analysis.
//'     - \code{informationOverall}: Cumulative information for the overall
//'       study at each analysis.
//'     - \code{efficacyStopping}: Indicator of whether efficacy stopping is
//'       permitted.
//'     - \code{futilityStopping}: Indicator of whether futility stopping is
//'       permitted.
//'     - \code{rejectPerStageH0}: Probability of efficacy stopping under the
//'       global null.
//'     - \code{futilityPerStageH0}: Probability of futility stopping under the
//'       global null.
//'     - \code{cumulativeRejectionH0}: Cumulative probability of efficacy
//'       stopping under the global null.
//'     - \code{cumulativeFutilityH0}: Cumulative probability of futility
//'       stopping under the global null.
//'
//' * \code{byArmResults}: A data frame containing:
//'     - \code{theta}: Parameter values for the active arms.
//'     - \code{selectAsBest}: Probability an arm is selected as best at the
//'       end of Phase 2.
//'     - \code{powerByArm}: Probability of rejecting the null for each arm by
//'       trial end.
//'     - \code{condPowerByArm}: Conditional power for each arm given it was
//'       selected as the best at the end of Phase 2.
//'
//' * \code{settings}: A list of input settings:
//'     - \code{typeAlphaSpending}: Type of alpha spending function.
//'     - \code{parameterAlphaSpending}: Parameter value for the chosen alpha
//'       spending function.
//'     - \code{userAlphaSpending}: User-specified alpha spending values.
//'     - \code{typeBetaSpending}: Type of beta spending function.
//'     - \code{parameterBetaSpending}: Parameter value for the chosen beta
//'       spending function.
//'     - \code{userBetaSpending}: User-specified beta spending values.
//'     - \code{spendingTime}: Error-spending times at each analysis.
//'
//' @details If \code{corr_known} is \code{FALSE}, critical boundaries are
//' computed assuming independence among the Phase-2 Wald statistics
//' (a conservative assumption). Power calculations, however, use the
//' correlation implied by the randomization ratio \eqn{r}.
//'
//' Futility boundaries may be supplied directly on the Z scale, derived from
//' conditional power, derived from parameter values, or computed from a beta
//' spending function.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive two-stage seamless sequential design for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
//'
//' @examples
//'
//' # Example 1: obtain the maximum information given power with no futility
//' (design1 <- getDesign_seamless(
//'   beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   K = 2, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' # Example 2: obtain power given the maximum information and a futility rule
//' (design2 <- getDesign_seamless(
//'   IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   K = 2, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF",
//'   futilityBounds = c(0.0, 0.5)))
//'
//' # Example 3: derive futility boundaries using beta spending
//' (design3 <- getDesign_seamless(
//'   beta = 0.1, theta = c(-log(0.5), -log(0.7)),
//'   M = 2, r = 1.0, corr_known = FALSE,
//'   K = 2, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   typeBetaSpending = "sfHSD", parameterBetaSpending = -2))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getDesign_seamless(
    const double beta = NA_REAL,
    const double IMax = NA_REAL,
    const Rcpp::NumericVector& theta = NA_REAL,
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int K = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::Nullable<Rcpp::NumericVector> criticalValues = R_NilValue,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBounds = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityCP = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityTheta = R_NilValue,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const Rcpp::NumericVector& userBetaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {

  std::vector<double> thetaVec(theta.begin(), theta.end());
  std::vector<double> infoRates(informationRates.begin(), informationRates.end());
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  std::vector<double> userAlpha(userAlphaSpending.begin(), userAlphaSpending.end());
  std::vector<double> userBeta(userBetaSpending.begin(), userBetaSpending.end());
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());

  std::vector<double> critValues, futBounds, futCP, futTheta;
  if (criticalValues.isNotNull()) {
    critValues = Rcpp::as<std::vector<double>>(criticalValues);
  } else {
    critValues = std::vector<double>(1, NaN);
  }

  if (futilityBounds.isNotNull()) {
    futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  } else {
    futBounds = std::vector<double>(1, NaN);
  }

  if (futilityCP.isNotNull()) {
    futCP = Rcpp::as<std::vector<double>>(futilityCP);
  } else {
    futCP = std::vector<double>(1, NaN);
  }

  if (futilityTheta.isNotNull()) {
    futTheta = Rcpp::as<std::vector<double>>(futilityTheta);
  } else {
    futTheta = std::vector<double>(1, NaN);
  }

  auto cpp_result = getDesign_seamless_cpp(
    beta, IMax, thetaVec, static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(K), infoRates, effStopping, futStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha,
    futBounds, futCP, futTheta, typeBetaSpending,
    parameterBetaSpending, userBeta, spendTime
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "seamless";
  return result;
}


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
    const std::vector<double>& spendingTimeNew) {


  // ----------- Start of Input Validation ----------- //
  if (std::isnan(betaNew) && std::isnan(INew)) {
    throw std::invalid_argument("betaNew and INew cannot be missing simultaneously");
  }
  if (!std::isnan(betaNew) && !std::isnan(INew)) {
    throw std::invalid_argument("Only one of betaNew and INew should be provided");
  }
  if (!std::isnan(INew) && INew <= 0.0) {
    throw std::invalid_argument("INew must be positive");
  }

  if (M < 1) throw std::invalid_argument("M must be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r must be positive");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(theta)) throw std::invalid_argument("theta must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (K <= L) throw std::invalid_argument("K must be greater than L");

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1.0)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(betaNew) && (betaNew < 0.0001 || betaNew >= 1.0)) {
    throw std::invalid_argument("betaNew must lie in [0.0001, 1)");
  }

  size_t kMax = K + 1;

  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (informationRates.size() != kMax) {
      throw std::invalid_argument("Invalid length for informationRates");
    }
    if (informationRates[0] <= 0.0) {
      throw std::invalid_argument("informationRates must be positive");
    }
    if (any_nonincreasing(informationRates)) {
      throw std::invalid_argument("informationRates must be increasing");
    }
    if (informationRates.back() != 1.0) {
      throw std::invalid_argument("informationRates must end with 1");
    }
    infoRates = informationRates;
  } else {
    for (size_t i = 0; i < kMax; ++i) {
      infoRates[i] = static_cast<double>(i + 1) / static_cast<double>(kMax);
    }
  }

  // efficacyStopping
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != kMax) {
      throw std::invalid_argument("Invalid length for efficacyStopping");
    }
    if (efficacyStopping.back() != 1) {
      throw std::invalid_argument("efficacyStopping must end with 1");
    }
    effStopping = efficacyStopping;
  } else {
    effStopping.assign(kMax, 1);
  }

  // futilityStopping
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != kMax) {
      throw std::invalid_argument("Invalid length for futilityStopping");
    }
    if (futilityStopping.back() != 1) {
      throw std::invalid_argument("futilityStopping must end with 1");
    }
    futStopping = futilityStopping;
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds =
      !none_na(futilityBounds) && !none_na(futilityCP) && !none_na(futilityTheta);

  if (!missingCriticalValues && criticalValues.size() != kMax) {
    throw std::invalid_argument("Invalid length for criticalValues");
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }

  std::string asf = typeAlphaSpending;
  for (char& c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues &&
      !(asf == "of" || asf == "p" || asf == "wt" ||
        asf == "sfof" || asf == "sfp" || asf == "sfkd" ||
        asf == "sfhsd" || asf == "user" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }

  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }

  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument("parameterAlphaSpending must be positive for sfKD");
  }

  if (missingCriticalValues && asf == "user") {
    if (!none_na(userAlphaSpending)) {
      throw std::invalid_argument("userAlphaSpending must be specified");
    }
    if (userAlphaSpending.size() != kMax) {
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    }
    if (userAlphaSpending[0] < 0.0) {
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    }
    if (any_nonincreasing(userAlphaSpending)) {
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    }
    if (userAlphaSpending.back() != alpha) {
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
    }
  }

  if (!missingFutilityBounds) {
    if (none_na(futilityBounds)) {
      if (futilityBounds.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityBounds");
      }
    } else if (none_na(futilityCP)) {
      if (futilityCP.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityCP");
      }
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (futilityCP[i] < 0.0 || futilityCP[i] > 1.0) {
          throw std::invalid_argument("futilityCP must lie in [0, 1]");
        }
      }
    } else {
      if (futilityTheta.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityTheta");
      }
    }
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax) {
      throw std::invalid_argument("Invalid length for spendingTime");
    }
    if (spendingTime[0] <= 0.0) {
      throw std::invalid_argument("spendingTime must be positive");
    }
    if (any_nonincreasing(spendingTime)) {
      throw std::invalid_argument("spendingTime must be increasing");
    }
    if (spendingTime.back() != 1.0) {
      throw std::invalid_argument("spendingTime must end with 1");
    }
    spendTime = spendingTime;
  } else {
    spendTime = infoRates;
  }


  // ----------- New Design Input Validation ----------- //
  std::vector<double> infoRatesNew;
  std::vector<unsigned char> effStoppingNew;
  std::vector<unsigned char> futStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char& c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  std::string bsfNew = typeBetaSpendingNew;
  for (char& c : bsfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  size_t k1 = K - L;
  if (!MullerSchafer) {
    infoRatesNew.resize(k1);
    for (size_t i = 0; i < k1; ++i) {
      infoRatesNew[i] =
        (infoRates[L + 1 + i] - infoRates[L]) / (1.0 - infoRates[L]);
    }

    effStoppingNew.resize(k1);
    std::memcpy(effStoppingNew.data(), effStopping.data() + L + 1,
                k1 * sizeof(unsigned char));

    if (none_na(futilityStoppingNew)) {
      if (futilityStoppingNew.size() != k1)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew.back() != 1)
        throw std::invalid_argument("futilityStoppingNew must end with 1");
      futStoppingNew = futilityStoppingNew; // copy
    } else {
      futStoppingNew.resize(k1);
      std::memcpy(futStoppingNew.data(), futStopping.data() + L + 1,
                  k1 * sizeof(unsigned char));
    }

    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() != k1)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew.back() != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = infoRatesNew;
    }
  } else {
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
      if (informationRatesNew.back() != 1.0)
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
      if (efficacyStoppingNew.back() != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew; // copy
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    // futStoppingNew: default to all 1s if missing
    if (none_na(futilityStoppingNew)) {
      if (futilityStoppingNew.size() != kNew)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew.back() != 1)
        throw std::invalid_argument("futilityStoppingNew must end with 1");
      futStoppingNew = futilityStoppingNew; // copy
    } else {
      futStoppingNew.assign(kNew, 1);
    }

    if (!(asfNew == "of" || asfNew == "p" || asfNew == "wt" ||
          asfNew == "sfof" || asfNew == "sfp" ||
          asfNew == "sfkd" || asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }
    if ((asfNew == "wt" || asfNew == "sfkd" || asfNew == "sfhsd") &&
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
      if (spendingTimeNew.back() != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = infoRatesNew;
    }
  }

  size_t k2 = MullerSchafer ? kNew : k1;

  bool missingFutilityBoundsInt =
    !none_na(futilityBoundsInt) && !none_na(futilityCPInt) &&
    !none_na(futilityThetaInt);

  if (!missingFutilityBoundsInt) {
    if (none_na(futilityBoundsInt)) {
      if (futilityBoundsInt.size() < k2 - 1) {
        throw std::invalid_argument("Insufficient length for futilityBoundsInt");
      }
    } else if (none_na(futilityCPInt)) {
      if (futilityCPInt.size() < k2 - 1) {
        throw std::invalid_argument("Insufficient length for futilityCPInt");
      }
      for (size_t i = 0; i < k2 - 1; ++i) {
        if (futilityCPInt[i] < 0.0 || futilityCPInt[i] > 1.0) {
          throw std::invalid_argument("futilityCPInt must lie in [0, 1]");
        }
      }
    } else {
      if (futilityThetaInt.size() < k2 - 1) {
        throw std::invalid_argument("Insufficient length for futilityThetaInt");
      }
    }
  }

  if (std::isnan(INew)) {
    if (missingFutilityBoundsInt &&
        !(bsfNew == "sfof" || bsfNew == "sfp" || bsfNew == "sfkd" ||
        bsfNew == "sfhsd" || bsfNew == "user" || bsfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpendingNew");
    }
  } else {
    if (missingFutilityBoundsInt &&
        !(bsfNew == "sfof" || bsfNew == "sfp" || bsfNew == "sfkd" ||
        bsfNew == "sfhsd" || bsfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpendingNew");
    }
  }

  if ((bsfNew == "sfkd" || bsfNew == "sfhsd") &&
      std::isnan(parameterBetaSpendingNew)) {
    throw std::invalid_argument("Missing value for parameterBetaSpendingNew");
  }
  if (bsfNew == "sfkd" && parameterBetaSpendingNew <= 0.0) {
    throw std::invalid_argument(
        "parameterBetaSpendingNew must be positive for sfKD");
  }

  if (std::isnan(INew) && bsfNew == "user") {
    if (!none_na(userBetaSpendingNew))
      throw std::invalid_argument("userBetaSpendingNew must be specified");
    if (userBetaSpendingNew.size() != k2)
      throw std::invalid_argument("Invalid length of userBetaSpendingNew");
    if (userBetaSpendingNew[0] < 0.0)
      throw std::invalid_argument("userBetaSpendingNew must be nonnegative");
    if (any_nonincreasing(userBetaSpendingNew))
      throw std::invalid_argument("userBetaSpendingNew must be nondecreasing");
    if (userBetaSpendingNew.back() != betaNew)
      throw std::invalid_argument(
          "userBetaSpendingNew must end with specified betaNew");
  }

  // ----------- End of Input Validation ----------- //

  ExitProbResult probs;
  ExitProbSeamless probss;

  // set up efficacy bounds
  std::vector<double> zero(M, 0.0);
  std::vector<double> critValues = criticalValues;
  double alpha1 = alpha;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (kMax > 1 && criticalValues.size() == kMax) {
      bool hasNaN = false;
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[kMax-1])) haybittle = true;
    }

    if (haybittle) { // Haybittle & Peto
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (!effStopping[i]) critValues[i] = 8.0;
      }

      auto f = [&](double aval)->double {
        critValues[kMax-1] = aval;
        probss = exitprob_seamless_cpp(
          M, r, zero, corr_known, K, critValues, infoRates);
        double cpu = std::accumulate(probss.exitProbUpper.begin(),
                                     probss.exitProbUpper.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, 0.0, 8.0, 1e-6);
    } else {
      critValues = getBound_seamless_cpp(
        M, r, corr_known, K, infoRates, alpha, asf, parameterAlphaSpending,
        userAlphaSpending, spendTime, effStopping);
    }
  } else {
    for (size_t i = 0; i < kMax - 1; ++i) {
      if (!effStopping[i]) critValues[i] = 8.0;
    }
    probss = exitprob_seamless_cpp(M, r, zero, corr_known, K, critValues, infoRates);
    alpha1 = std::accumulate(probss.exitProbUpper.begin(),
                             probss.exitProbUpper.end(), 0.0);
  }

  // obtain futility bounds for the primary trial
  std::vector<double> futBounds(kMax);
  if (kMax > 1) {
    if (missingFutilityBounds) {
      std::fill_n(futBounds.begin(), kMax-1, -8.0);
      futBounds[kMax-1] = critValues[kMax-1];
    } else if (!missingFutilityBounds) {
      if (none_na(futilityBounds)) {
        for (size_t i = 0; i < kMax - 1; ++i) {
          if (futilityBounds[i] > critValues[i]) {
            throw std::invalid_argument(
                "futilityBounds must lie below criticalValues");
          }
        }
        std::copy_n(futilityBounds.begin(), kMax-1, futBounds.begin());
        futBounds[kMax-1] = critValues[kMax-1];
      } else if (none_na(futilityCP)) {
        double c2 = critValues[kMax - 1];
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i]) *
            (c2 - std::sqrt(1 - infoRates[i]) * boost_qnorm(1 - futilityCP[i]));
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityCP values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax-1] = critValues[kMax-1];
      } else {
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i] * IMax) * futilityTheta[i];
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityTheta values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax-1] = critValues[kMax-1];
      }
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // information for the primary trial
  std::vector<double> information1(kMax);
  for (size_t i = 0; i < kMax; ++i) information1[i] = IMax * infoRates[i];

  std::vector<double> r1(k1), b1(k1), a1(k1, -8.0), zero1(k1, 0.0);
  for (size_t i = 0; i < k1; ++i) {
    r1[i] = infoRates[L] / infoRates[i + L + 1];
    b1[i] = (critValues[i + L + 1] - std::sqrt(r1[i]) * zL) / std::sqrt(1.0 - r1[i]);
    if (!effStoppingNew[i]) b1[i] = 8.0;
  }

  // conditional type I error
  probs = exitprobcpp(b1, a1, zero1, infoRatesNew);
  auto v0 = probs.exitProbUpper;
  double alphaNew = std::accumulate(v0.begin(), v0.end(), 0.0);

  // conditional power
  for (size_t i = 0; i < k1; ++i) {
    a1[i] = (futBounds[i + L + 1] - std::sqrt(r1[i]) * zL) / std::sqrt(1.0 - r1[i]);
    if (!futStoppingNew[i]) a1[i] = -8.0;
  }

  std::vector<double> I1(k1);
  for (size_t i = 0; i < k1; ++i) {
    I1[i] = information1[i + L + 1] - information1[L];
  }

  std::vector<double> theta1(k1, theta);
  probs = exitprobcpp(b1, a1, theta1, I1);
  double conditionalPower = std::accumulate(probs.exitProbUpper.begin(),
                                            probs.exitProbUpper.end(), 0.0);

  double IL = information1[L];
  double sqrtIL = std::sqrt(IL);
  double zscaled = zL * sqrtIL;

  // critical values for the secondary trial
  std::vector<double> b2;
  std::string asf2;
  double asfpar2;
  std::vector<double> cpu0(k2);
  if (!MullerSchafer) {
    asf2 = "user";
    asfpar2 = NaN;
    std::partial_sum(v0.begin(), v0.end(), cpu0.begin());
    b2 = b1;
  } else {
    asf2 = asfNew;
    asfpar2 = parameterAlphaSpendingNew;
    if (asf2 != "none" && asf2 != "of" && asf2 != "p" && asf2 != "wt") {
      for (size_t i = 0; i < k2; ++i) {
        cpu0[i] = errorSpentcpp(spendTimeNew[i], alphaNew, asf2, asfpar2);
      }
    }
    b2 = getBoundcpp(k2, infoRatesNew, alphaNew, asfNew,
                     parameterAlphaSpendingNew, {NaN},
                     spendTimeNew, effStoppingNew);
  }


  // futility boundaries for the secondary trial
  std::vector<double> a2(k2, -8.0);

  // critical values and futility boundaries for the integrated trial
  std::vector<double> critValues2(k2);
  std::vector<double> futBounds2(k2, -8.0);
  std::vector<double> theta2(k2, theta);

  std::vector<double> I2(k2), Ic(k2), sqrtI2(k2), sqrtIc(k2);

  if (std::isnan(betaNew)) {
    for (size_t i = 0; i < k2; ++i) {
      I2[i] = INew * infoRatesNew[i];
      Ic[i] = I2[i] + IL;
      sqrtI2[i] = std::sqrt(I2[i]);
      sqrtIc[i] = std::sqrt(Ic[i]);
    }

    for (size_t i = 0; i < k2; ++i) {
      critValues2[i] = (b2[i] * sqrtI2[i] + zscaled) / sqrtIc[i];
    }

    // now compute futility bounds if needed
    if (k2 > 1) {
      if (missingFutilityBoundsInt && bsfNew == "none") {
        std::fill_n(futBounds2.begin(), k2 - 1, -8.0);
        futBounds2[k2-1] = critValues2[k2-1];
      } else if (!missingFutilityBoundsInt) {
        if (none_na(futilityBoundsInt)) {
          for (size_t i = 0; i < k2 - 1; ++i) {
            if (futilityBoundsInt[i] > critValues2[i]) {
              throw std::invalid_argument(
                  "futilityBoundsInt must lie below critical values "
                  "for the integrated trial");
            }
          }
          std::copy_n(futilityBoundsInt.begin(), k2-1, futBounds2.begin());
          futBounds2[k2-1] = critValues2[k2-1];
        } else if (none_na(futilityCPInt)) {
          double c2 = critValues2[k2 - 1];
          for (size_t i = 0; i < k2 - 1; ++i) {
            double q = boost_qnorm(1 - futilityCPInt[i]);
            double sc = Ic[i] / Ic[k2 - 1];
            futBounds2[i] = std::sqrt(sc) * (c2 - std::sqrt(1 - sc) * q);
            if (futBounds2[i] > critValues2[i]) {
              throw std::invalid_argument(
                  "futilityCPInt values are too large to be compatible with "
                  "critical values for the integrated trial");
            }
          }
          futBounds2[k2-1] = critValues2[k2-1];
        } else {
          for (size_t i = 0; i < k2 - 1; ++i) {
            futBounds2[i] = std::sqrt(Ic[i]) * futilityThetaInt[i];
            if (futBounds2[i] > critValues2[i]) {
              throw std::invalid_argument(
                  "futilityThetaInt values are too large to be compatible with "
                  "critical values for the integrated trial");
            }
          }
          futBounds2[k2-1] = critValues2[k2-1];
        }
      }
    } else {
      if (missingFutilityBoundsInt) {
        futBounds2 = critValues2;
      }
    }

    if (missingFutilityBoundsInt && bsfNew != "none" && k2 > 1) { // beta-spending
      std::vector<double> wc(k2, 1.0);
      auto out = getPower(
        alphaNew, k2, critValues2, theta2, Ic, bsfNew,
        parameterBetaSpendingNew, spendTimeNew, futStoppingNew,
        wc, IL, theta, zL);
      futBounds2 = out.futilityBounds;
    }

    for (size_t i = 0; i < k2; ++i) {
      a2[i] = (futBounds2[i] * sqrtIc[i] - zscaled) / sqrtI2[i];
    }
  } else {
    // obtain required max information for the secondary trial given target power
    std::vector<double> u; u.reserve(k2);
    std::vector<double> l; l.reserve(k2);

    auto f = [&](double x)->double {
      double Inew = sq(x / theta);
      for (size_t i = 0; i < k2; ++i) {
        I2[i] = Inew * infoRatesNew[i];
        Ic[i] = I2[i] + IL;
        sqrtI2[i] = std::sqrt(I2[i]);
        sqrtIc[i] = std::sqrt(Ic[i]);
      }

      double mu0 = theta * sqrtI2[0];

      for (size_t i = 0; i < k2; ++i) {
        critValues2[i] = (b2[i] * sqrtI2[i] + zscaled) / sqrtIc[i];
      }

      // now compute futility bounds if needed
      if (k2 > 1) {
        if (missingFutilityBoundsInt && bsfNew == "none") {
          std::fill_n(futBounds2.begin(), k2 - 1, -8.0);
          futBounds2[k2-1] = critValues2[k2-1];
        } else if (!missingFutilityBoundsInt) {
          if (none_na(futilityBoundsInt)) {
            for (size_t i = 0; i < k2 - 1; ++i) {
              if (futilityBoundsInt[i] > critValues2[i]) {
                throw std::invalid_argument(
                    "futilityBoundsInt must lie below critical values "
                    "for the integrated trial");
              }
            }
            std::copy_n(futilityBoundsInt.begin(), k2-1, futBounds2.begin());
            futBounds2[k2-1] = critValues2[k2-1];
          } else if (none_na(futilityCPInt)) {
            double c2 = critValues2[k2 - 1];
            for (size_t i = 0; i < k2 - 1; ++i) {
              double q = boost_qnorm(1 - futilityCPInt[i]);
              double sc = Ic[i] / Ic[k2 - 1];
              futBounds2[i] = std::sqrt(sc) * (c2 - std::sqrt(1 - sc) * q);
              if (futBounds2[i] > critValues2[i]) {
                throw std::invalid_argument(
                    "futilityCPInt values are too large to be compatible with "
                    "critical values for the integrated trial");
              }
            }
            futBounds2[k2-1] = critValues2[k2-1];
          } else {
            for (size_t i = 0; i < k2 - 1; ++i) {
              futBounds2[i] = std::sqrt(Ic[i]) * futilityThetaInt[i];
              if (futBounds2[i] > critValues2[i]) {
                throw std::invalid_argument(
                    "futilityThetaInt values are too large to be compatible with "
                    "critical values for the integrated trial");
              }
            }
            futBounds2[k2-1] = critValues2[k2-1];
          }
        }
      } else {
        if (missingFutilityBoundsInt) {
          futBounds2 = critValues2;
        }
      }


      if (!missingFutilityBoundsInt || bsfNew == "none" || k2 == 1) {
        for (size_t i = 0; i < k2 - 1; ++i) {
          a2[i] = (futBounds2[i] * sqrtIc[i] - zscaled) / sqrtI2[i];
        }
        a2[k2 - 1] = b2[k2 - 1];

        probs = exitprobcpp(b2, a2, theta2, I2);
        double overallReject = std::accumulate(probs.exitProbUpper.begin(),
                                               probs.exitProbUpper.end(), 0.0);
        return (1.0 - overallReject) - betaNew;
      } else {
        // initialize futility bound to be updated
        std::fill(futBounds2.begin(), futBounds2.end(), -8.0);
        std::fill(a2.begin(), a2.end(), -8.0);
        double eps = 0.0, cb = 0.0;

        // first stage
        if (futStoppingNew[0]) {
          cb = (bsfNew == "user") ? userBetaSpendingNew[0] :
          errorSpentcpp(spendTimeNew[0], betaNew, bsfNew,
                        parameterBetaSpendingNew);

          eps = boost_pnorm(b2[0] - mu0) - cb;
          if (eps < 0.0) return -1.0; // to decrease beta
          a2[0] = boost_qnorm(cb) + mu0;
          futBounds2[0] = (a2[0] * sqrtI2[0] + zscaled) / sqrtIc[0];
        }

        // subsequent stages
        for (size_t k = 1; k < k2; ++k) {
          if (futStoppingNew[k]) {
            cb = (bsfNew == "user") ? userBetaSpendingNew[k] :
            errorSpentcpp(spendTimeNew[k], betaNew, bsfNew,
                          parameterBetaSpendingNew);

            a2[k-1] = (futBounds2[k-1] * sqrtIc[k-1] - zscaled) / sqrtI2[k-1];

            u.resize(k + 1);
            l.resize(k + 1);
            std::memcpy(u.data(), b2.data(), k * sizeof(double));
            std::memcpy(l.data(), a2.data(), k * sizeof(double));

            // lambda expression for finding futility bound at stage k
            auto g = [&](double aval)->double {
              a2[k] = (aval * sqrtIc[k] - zscaled) / sqrtI2[k];
              l[k] = a2[k];
              u[k] = l[k];
              probs = exitprobcpp(u, l, theta2, I2);
              double cpl = std::accumulate(probs.exitProbLower.begin(),
                                           probs.exitProbLower.end(), 0.0);
              return cpl - cb;
            };

            double bk = critValues2[k];
            eps = g(bk);
            double g_minus8 = g(-8.0);
            if (g_minus8 > 0.0) { // no beta spent at current visit
              futBounds2[k] = -8.0;
            } else if (eps > 0.0) {
              auto g_for_brent = [&](double aval)->double {
                if (aval == -8.0) return g_minus8;  // avoid recomputation at 8.0
                if (aval == bk) return eps;         // avoid recomputation at b[k]
                return g(aval);
              };
              futBounds2[k] = brent(g_for_brent, -8.0, bk, 1e-6);
            } else if (k < k2-1) {
              return -1.0;
            }
          }
        }

        return eps;
      }
    };

    double drift = brent(f, 0.001, 8.0, 1e-6);
    INew = sq(drift / theta);
    futBounds2[k2-1] = critValues2[k2-1];
    a2[k2-1] = b2[k2-1];
  }


  probs = exitprobcpp(b2, a2, theta2, I2);
  std::vector<double> cpu1(k2);
  std::partial_sum(probs.exitProbUpper.begin(),
                   probs.exitProbUpper.end(), cpu1.begin());
  double p2 = cpu1.back();

  std::vector<double> cpu2(k2);
  std::partial_sum(probs.exitProbLower.begin(),
                   probs.exitProbLower.end(), cpu2.begin());


  // combined design
  size_t kc = L + 1 + k2;
  std::vector<double> Ic_full(kc); // cumulative information for the combined design
  std::copy_n(information1.data(), L + 1, Ic_full.data());
  std::copy_n(Ic.data(), k2, Ic_full.data() + L + 1);

  double IMaxc = Ic_full[kc - 1];
  std::vector<double> infoRates_full(kc);
  for (size_t i = 0; i < kc; ++i) infoRates_full[i] = Ic_full[i] / IMaxc;

  std::vector<double> critValues_full(kc);
  std::copy_n(critValues.data(), L + 1, critValues_full.data());
  std::copy_n(critValues2.data(), k2, critValues_full.data() + L + 1);

  std::vector<double> futBounds_full(kc);
  std::copy_n(futBounds.data(), L + 1, futBounds_full.data());
  std::copy_n(futBounds2.data(), k2, futBounds_full.data() + L + 1);

  ListCpp des1;
  des1.push_back(M, "M");
  des1.push_back(r, "r");
  des1.push_back(corr_known, "corr_known");
  des1.push_back(K, "K");
  des1.push_back(L, "L");
  des1.push_back(zL, "zL");
  des1.push_back(theta, "theta");
  des1.push_back(IMax, "maxInformation");
  des1.push_back(kMax, "kMax");
  des1.push_back(std::move(infoRates), "informationRates");
  des1.push_back(std::move(critValues), "efficacyBounds");
  des1.push_back(std::move(futBounds), "futilityBounds");
  des1.push_back(std::move(information1), "information");
  des1.push_back(alpha1, "alpha");
  des1.push_back(alphaNew, "conditionalAlpha");
  des1.push_back(conditionalPower, "conditionalPower");
  des1.push_back(MullerSchafer, "MullerSchafer");

  ListCpp des2;
  des2.push_back(p2, "overallReject");
  des2.push_back(alphaNew, "alpha");
  des2.push_back(k2, "kMax");
  des2.push_back(INew, "maxInformation");
  des2.push_back(std::move(infoRatesNew), "informationRates");
  des2.push_back(std::move(b2), "efficacyBounds");
  des2.push_back(std::move(a2), "futilityBounds");
  des2.push_back(std::move(cpu1), "cumulativeRejection");
  des2.push_back(std::move(cpu2), "cumulativeFutility");
  des2.push_back(std::move(cpu0), "cumulativeAlphaSpent");
  des2.push_back(std::move(I2), "information");
  des2.push_back(asf2, "typeAlphaSpending");
  des2.push_back(asfpar2, "parameterAlphaSpending");
  des2.push_back(bsfNew, "typeBetaSpending");
  des2.push_back(parameterBetaSpendingNew, "parameterBetaSpending");
  des2.push_back(userBetaSpendingNew, "userBetaSpending");
  des2.push_back(spendTimeNew, "spendingTime");

  ListCpp des3;
  des3.push_back(M, "M");
  des3.push_back(r, "r");
  des3.push_back(corr_known, "corr_known");
  des3.push_back(K, "K");
  des3.push_back(L, "L");
  des3.push_back(zL, "zL");
  des3.push_back(theta, "theta");
  des3.push_back(IMaxc, "maxInformation");
  des3.push_back(kc, "kMax");
  des3.push_back(std::move(infoRates_full), "informationRates");
  des3.push_back(std::move(critValues_full), "efficacyBounds");
  des3.push_back(std::move(futBounds_full), "futilityBounds");
  des3.push_back(std::move(Ic_full), "information");

  ListCpp result;
  result.push_back(std::move(des1), "primaryTrial");
  result.push_back(std::move(des2), "secondaryTrial");
  result.push_back(std::move(des3), "integratedTrial");
  return result;
}


//' @title Adaptive Phase 2/3 Seamless Design
//' @description
//' Calculates the conditional power for specified incremental
//' information, given the interim results, parameter value,
//' data-dependent changes in the error spending function, and
//' the number and spacing of interim looks. Conversely,
//' calculates the incremental information required to attain
//' a specified conditional power, given the interim results,
//' parameter value, data-dependent changes in
//' the error spending function, and the number and spacing of interim looks.
//'
//' @param betaNew The type II error for the secondary trial.
//' @param INew The maximum information for the active arm versus the common
//'   control in the secondary trial. Either
//'   \code{betaNew} or \code{INew} should be provided, while the other
//'   must be missing.
//' @param M Number of active treatment arms in Phase 2.
//' @param r Randomization ratio of each active arm to the common control
//'   in Phase 2.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param L The interim adaptation look in Phase 3.
//' @param zL The z-test statistic at the interim adaptation look of
//'   Phase 3.
//' @param theta The assumed treatment effect for the selected arm versus the
//'   common control.
//' @param IMax Maximum information for the active arm versus the common
//'   control for the original trial. Must be provided.
//' @param K Number of sequential looks in Phase 3.
//' @param informationRates A numeric vector of information rates fixed
//'   before the trial. If unspecified, defaults to \eqn{(1:(K+1)) / (K+1)}.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param futilityStopping Indicators of whether futility stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the max z-test statistic
//'   scale for Phase 2 and the z-test statistics for the selected arm
//'   in Phase 3 for the primary trial. If missing, boundaries
//'   will be computed based on the specified alpha spending function.
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
//' @param futilityBounds The lower boundaries on the max-z statistic scale
//'   at end of phase 2 and the z-test statistic scale in phase 3
//'   for futility stopping for the primary trial. Defaults to
//'   \code{rep(-8, kMax-1)} if left unspecified.
//' @param futilityCP The conditional power-based futility bounds for the
//'   primary trial.
//' @param futilityTheta The parameter value-based futility bounds for the
//'   primary trial.
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param futilityStoppingNew The indicators of whether futility stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param typeAlphaSpendingNew The type of alpha spending for the secondary
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
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param futilityBoundsInt The futility boundaries on the z statistic
//'   scale for new stages of the integrated trial.
//' @param futilityCPInt The conditional power-based futility bounds for
//'   new stages of the integrated trial.
//' @param futilityThetaInt The parameter value-based futility bounds for the
//'   new stages of the integrated trial.
//' @param typeBetaSpendingNew The type of beta spending for the secondary
//'   trial. One of the following:
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function,
//'   \code{"user"} for user defined spending, and
//'   \code{"none"} for no early futility stopping.
//'   Defaults to \code{"none"}.
//' @param parameterBetaSpendingNew The parameter value of beta spending
//'   for the secondary trial. Corresponds to \eqn{\rho} for \code{"sfKD"},
//'   and \eqn{\gamma} for \code{"sfHSD"}.
//' @param userBetaSpendingNew The user-defined cumulative beta spending.
//'   Represents the cumulative beta spent up to each stage of the
//'   secondary trial.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRatesNew}.
//'
//' @return An \code{adaptDesign_seamless} object with three list components:
//'
//' * \code{primaryTrial}: A list of selected information for the primary
//'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{K},
//'   \code{L}, \code{zL}, \code{theta}, \code{maxInformation}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
//'   \code{information}, \code{alpha}, \code{conditionalAlpha},
//'   \code{conditionalPower}, and \code{MullerSchafer}.
//'
//' * \code{secondaryTrial}: A list of selected information for the seconary
//'   trial, including \code{overallReject}, \code{alpha}, \code{kMax},
//'   \code{maxInformation}, \code{informationRates}, \code{efficacyBounds},
//'   \code{futilityBounds}, \code{cumulativeRejection},
//'   \code{cumulativeFutility}, \code{cumulativeAlphaSpent},
//'   \code{information}, \code{typeAlphaSpending},
//'   \code{parameterAlphaSpending}, \code{typeBetaSpending},
//'   \code{parameterBetaSpending}, \code{userBetaSpending}, and
//'   \code{spendingTime}.
//'
//' * \code{integratedTrial}: A list of selected information for the integrated
//'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{K},
//'   \code{L}, \code{zL}, \code{theta}, \code{maxInformation}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
//'   and \code{information}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive two-stage seamless sequential design for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
//'
//' @seealso \code{\link{getDesign_seamless}}
//'
//' @examples
//'
//' (des1 <- adaptDesign_seamless(
//'   betaNew = 0.1, M = 2, r = 1, corr_known = FALSE,
//'   L = 1, zL = -log(0.67) * sqrt(80 / 4), theta = -log(0.691),
//'   IMax = 120 / 4, K = 2, informationRates = c(1/3, 2/3, 1),
//'   alpha = 0.025, typeAlphaSpending = "OF", kNew = 1))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List adaptDesign_seamless(
    double betaNew = NA_REAL,
    double INew = NA_REAL,
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int L = NA_INTEGER,
    const double zL = NA_REAL,
    const double theta = NA_REAL,
    const double IMax = NA_REAL,
    const int K = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::Nullable<Rcpp::NumericVector> criticalValues = R_NilValue,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBounds = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityCP = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityTheta = R_NilValue,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int kNew = NA_INTEGER,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBoundsInt = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityCPInt = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityThetaInt = R_NilValue,
    const std::string& typeBetaSpendingNew = "none",
    const double parameterBetaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& userBetaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto futStoppingNew = convertLogicalVector(futilityStoppingNew);
  auto userBetaNew = Rcpp::as<std::vector<double>>(userBetaSpendingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  std::vector<double> critValues, futBounds, futCP, futTheta;
  if (criticalValues.isNotNull()) {
    critValues = Rcpp::as<std::vector<double>>(criticalValues);
  } else {
    critValues = std::vector<double>(1, NaN);
  }

  if (futilityBounds.isNotNull()) {
    futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  } else {
    futBounds = std::vector<double>(1, NaN);
  }

  if (futilityCP.isNotNull()) {
    futCP = Rcpp::as<std::vector<double>>(futilityCP);
  } else {
    futCP = std::vector<double>(1, NaN);
  }

  if (futilityTheta.isNotNull()) {
    futTheta = Rcpp::as<std::vector<double>>(futilityTheta);
  } else {
    futTheta = std::vector<double>(1, NaN);
  }

  std::vector<double> futBoundsInt, futCPInt, futThetaInt;
  if (futilityBoundsInt.isNotNull()) {
    futBoundsInt = Rcpp::as<std::vector<double>>(futilityBoundsInt);
  } else {
    futBoundsInt = std::vector<double>(1, NaN);
  }

  if (futilityCPInt.isNotNull()) {
    futCPInt = Rcpp::as<std::vector<double>>(futilityCPInt);
  } else {
    futCPInt = std::vector<double>(1, NaN);
  }

  if (futilityThetaInt.isNotNull()) {
    futThetaInt = Rcpp::as<std::vector<double>>(futilityThetaInt);
  } else {
    futThetaInt = std::vector<double>(1, NaN);
  }

  auto cpp_result = adaptDesign_seamless_cpp(
    betaNew, INew, static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(L), zL, theta, IMax,
    static_cast<size_t>(K), infoRates, effStopping,
    futStopping, critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, futCP, futTheta, spendTime,
    MullerSchafer, static_cast<size_t>(kNew), infoRatesNew,
    effStoppingNew, futStoppingNew, typeAlphaSpendingNew,
    parameterAlphaSpendingNew, futBoundsInt, futCPInt, futThetaInt,
    typeBetaSpendingNew, parameterBetaSpendingNew, userBetaNew,
    spendTimeNew
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "adaptDesign_seamless";
  return result;
}
