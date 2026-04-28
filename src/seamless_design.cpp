#include "dataframe_list.h"
#include "thread_utils.h"
#include "utilities.h"
#include "generic_design.h"
#include "mvnormr.h"

#include <Rcpp.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using std::size_t;


ListCpp exitprob_seamless_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t K,
    const std::vector<double>& b,
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
  if (b.size() != K + 1) {
    throw std::invalid_argument("b should have length K + 1");
  }

  std::vector<double> Ivec; Ivec.reserve(K+1);
  if (none_na(I)) {
    if (I.size() < K+1) throw std::invalid_argument("Insufficient length for I");

    Ivec.assign(I.begin(), I.begin() + K + 1);
    if (Ivec[0] <= 0.0) throw std::invalid_argument("I must be positive");
    if (any_nonincreasing(Ivec)) throw std::invalid_argument("I must be increasing");
  } else {
    Ivec.resize(K+1);
    std::iota(Ivec.begin(), Ivec.end(), 1.0);
  }

  double rho = corr_known ? r / (r + 1.0) : 0;


  std::vector<double> preject(K + 1, 0.0);
  FlatMatrix preject_by_arm(K + 1, M);
  std::vector<double> select_as_best(M, 0.0);

  // compute the exit probability in phase 2
  std::vector<double> sqrtI(K + 1);
  for (size_t k = 0; k < K + 1; ++k) {
    sqrtI[k] = std::sqrt(Ivec[k]);
  }
  double sqrtI0 = sqrtI[0];
  std::vector<double> mean0(M), lower0(M, -6.0), upper0(M, b[0]);
  for (size_t m = 0; m < M; ++m) {
    mean0[m] = theta[m] * sqrtI0;
  }

  // covariance matrix of the Wald statistics for the M active arms in phase 2
  FlatMatrix sigma0(M, M);
  for (size_t j = 0; j < M; ++j) {
    for (size_t i = 0; i < M; ++i) {
      sigma0(i, j) = (i == j) ? 1.0 : rho;
    }
  }

  PMVNResult out0 = pmvnormcpp(lower0, upper0, mean0, sigma0,
                               1024, 16384, 8, 1e-4, 0.0, 314159, true);
  preject[0] = 1.0 - out0.prob;

  // compute the exit probabilities at the K looks in phase 3
  std::vector<double> I1(K), sqrtI1(K);
  for (size_t k = 0; k < K; ++k) {
    I1[k] = Ivec[k+1] - Ivec[0]; // information for the selected arm in phase 3
    sqrtI1[k] = std::sqrt(I1[k]);
  }

  for (size_t m = 0; m < M; ++m) { // loop over the selected arm in phase 2
    double mu = theta[m] * sqrtI0;

    // compute the conditional distribution of the Wald statistics for the
    // M - 1 non-selected arms in phase 2 given that arm m is selected, and
    // then compute the probability of selecting arm m in phase 2 (i.e.,
    // the probability that the Wald statistic for arm m is larger than
    // the Wald statistics for the M - 1 non-selected arms in phase 2)
    auto f0 = [M, m, mu, theta, sqrtI0, rho](double z0)->double {
      if (M > 1) {
        std::vector<double> mean(M - 1), lower(M - 1, -6.0), upper(M - 1, z0);
        size_t j = 0;
        double delta0 = z0 - mu;
        for (size_t i = 0; i < M; ++i) {
          if (i != m) {
            mean[j] = theta[i] * sqrtI0 + rho * delta0;
            ++j;
          }
        }

        FlatMatrix sigma(M - 1, M - 1);
        for (size_t j = 0; j < M - 1; ++j) {
          for (size_t i = 0; i < M - 1; ++i) {
            sigma(i, j) = (i == j) ? 1.0 - rho * rho : rho * (1.0 - rho);
          }
        }

        PMVNResult out = pmvnormcpp(lower, upper, mean, sigma,
                                    1024, 16384, 8, 1e-4, 0.0, 314159, true);
        // multiply the density of the Wald statistic for arm m at z0
        return out.prob * boost_dnorm(z0, mu, 1.0);
      } else {
        return boost_dnorm(z0, mu, 1.0);
      }
    };

    std::vector<double> breaks0 = {-6.0, 6.0};
    select_as_best[m] = integrate3(f0, breaks0, 1e-4);

    std::vector<double> breaks1 = {b[0], 6.0};
    preject_by_arm(0, m) = integrate3(f0, breaks1, 1e-4);

    // compute the conditional distribution of the Wald statistics in phase 3
    // given that arm m is selected in phase 2, and then compute the probability
    // of exiting at look k in phase 3 for k = 1,...,K. Note that the conditional
    // distribution the Wald statistics for arm m in phase 3 can be written as
    // a secondary group sequential trial with K looks, where the information
    // level at look k is I1[k], the critical value at look k is b1[k], and
    // the treatment effect is theta[m]
    std::vector<double> a1(K, -6.0);
    std::vector<double> theta1(K, theta[m]);
    auto f1 = [K, b, a1, theta1, sqrtI, sqrtI0, sqrtI1, I1](double z0) {
      std::vector<double> b1(K);
      for (size_t k = 0; k < K; ++k) {
        b1[k] = (b[k+1] * sqrtI[k+1] - z0 * sqrtI0) / sqrtI1[k];
        if (b1[k] < -6.0) b1[k] = -6.0;
      }

      ListCpp probs = exitprobcpp(b1, a1, theta1, I1);
      auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
      return exitUpper;
    };

    // compute the joint distribution of the Wald statistic for arm m and the
    // Wald statistics for the M - 1 non-selected arms in phase 2, and then
    // compute the probability of selecting arm m in phase 2 and exiting at look k
    // in phase 3 for k = 1,...,K
    auto g = [K, f0, f1](double z0)->std::vector<double> {
      std::vector<double> vec(K);
      double p0 = f0(z0);
      std::vector<double> p1 = f1(z0);
      for (size_t k = 0; k < K; ++k) {
        vec[k] = p0 * p1[k];
      }
      return vec;
    };

    std::vector<double> breaks = {-6.0, b[0]};
    for (size_t k = 0; k < K; ++k) {
      auto h = [k, g](double z0)->double {
        return g(z0)[k];
      };

      preject_by_arm(k+1, m) = integrate3(h, breaks, 1e-4);

      preject[k+1] += preject_by_arm(k+1, m);
    }
  }

  ListCpp result;
  result.push_back(std::move(preject), "exitProb");
  result.push_back(std::move(preject_by_arm), "exitProbByArm");
  result.push_back(std::move(select_as_best), "selectAsBest");
  return result;
}


//' @title Exit Probabilities for Two-Stage Seamless Sequential Design
//' @description Computes the exit (rejection) probabilities for a two-stage
//' selection and testing design. In Phase 2, multiple active arms are
//' compared against a common control arm. The best-performing arm is
//' selected to proceed to Phase 3, where it is tested against the common
//' control over multiple looks.
//'
//' @param M Number of active treatment arms in Phase 2 (\eqn{M \ge 2}).
//' @param r Randomization ratio of each active arm to the common control
//'   in Phase 2.
//' @param theta A vector of length \eqn{M} representing the true treatment
//'   effects for each active arm versus the common control.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param K Number of sequential looks in Phase 3.
//' @param b A vector of critical values (length \eqn{K+1}). The first element
//'   is for Phase 2; the remaining \eqn{K} elements are for the looks in
//'   Phase 3.
//' @param I A vector of information levels (length \eqn{K+1}) for any active
//'   arm versus the common control. The first element is for Phase 2;
//'   the remaining \eqn{K} elements are for the looks in Phase 3.
//'
//' @details
//' The function assumes a multivariate normal distribution for the Wald
//' statistics. The "best" arm is defined as the active arm with the largest
//' score statistic at the end of Phase 2.
//'
//' \strong{Decision Rules:}
//' * \strong{Phase 2}: The global null hypothesis is rejected if the Wald
//'   statistic for the best arm, \eqn{Z(I_0)}, satisfies \eqn{Z(I_0) \ge b_0}.
//'
//' * \strong{Phase 3}: If the trial continues, the hypothesis is rejected at
//'   look \eqn{k} if \eqn{Z(I_k) \ge b_k} and all previous
//'   looks (including Phase 2) failed to reject.
//'
//' \strong{Design Assumptions:}
//'
//' * All active arms share the same information level in Phase 2.
//'
//' * Exactly one active arm is selected at the end of Phase 2 based on the
//'   largest observed statistic.
//'
//' @return A list containing:
//'
//' * \code{exitProb}: A vector of length \eqn{K + 1}. The first element is the
//' probability of rejection in Phase 2; the remaining elements are the
//' probabilities of rejection at each look in Phase 3.
//'
//' * \code{exitProbByArm}: A \eqn{(K+1) \times M} matrix. The \eqn{(k, m)}-th
//' entry represents the probability that the global null is rejected at
//' look \eqn{k} given that arm \eqn{m} was selected as the best arm.
//'
//' * \code{selectAsBest}: A vector of length \eqn{M} containing the probability
//' that each active arm is selected to move on to Phase 3.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' # Setup: 2 active arms vs control in phase 2; 1 selected arm vs control
//' # in phase 3. Phase 3 has 2 sequential looks.
//'
//' # Information levels: equal spacing over 3 looks based on max 110 patients
//' # per arm, SD = 1.0
//' I <- c(110 / (2 * 1.0^2) * seq(1, 3)/3)
//'
//' # O'Brien-Fleming critical values
//' b <- c(3.776606, 2.670463, 2.180424)
//'
//' # Type I error under the global null hypothesis
//' p0 <- exitprob_seamless(M = 2, theta = c(0, 0), K = 2, b = b, I = I)
//' cumsum(p0$exitProb)
//'
//' # Power under alternative: Treatment effects of 0.3 and 0.5
//' p1 <- exitprob_seamless(M = 2, theta = c(0.3, 0.5), K = 2, b = b, I = I)
//' cumsum(p1$exitProb)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List exitprob_seamless(
    const int M = NA_INTEGER,
    const double r = 1,
    const Rcpp::NumericVector& theta = NA_REAL,
    const bool corr_known = true,
    const int K = NA_INTEGER,
    const Rcpp::NumericVector& b = NA_REAL,
    const Rcpp::NumericVector& I = NA_REAL) {

  std::vector<double> thetaVec(theta.begin(), theta.end());
  std::vector<double> bVec(b.begin(), b.end());
  std::vector<double> IVec(I.begin(), I.end());
  auto result = exitprob_seamless_cpp(
    static_cast<size_t>(M), r, thetaVec, corr_known,
    static_cast<size_t>(K), bVec, IVec);
  return Rcpp::wrap(result);
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


  std::vector<double> criticalValues(kMax);

  if (asf == "none") {
    for (size_t i = 0; i < kMax-1; ++i) criticalValues[i] = 6.0;
    std::vector<double> theta(M, 0.0);

    auto f = [&](double aval)->double {
      criticalValues[kMax-1] = aval;
      auto probs = exitprob_seamless_cpp(M, r, theta, corr_known, k,
                                      criticalValues, infoRates);
      auto v = probs.get<std::vector<double>>("exitProb");
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    criticalValues[kMax-1] = brent(f, 0.0, 6.0, 1e-6);
    return subset(criticalValues, 0, K);
  }

  if (asf == "of" || asf == "p" || asf == "wt") {
    double Delta;
    if (asf == "of") Delta = 0.0;
    else if (asf == "p") Delta = 0.5;
    else Delta = parameterAlphaSpending; // parameterAlphaSpending holds delta for WT

    // for a given multiplier, compute cumulative upper exit probability - alpha
    std::vector<double> u(kMax);
    std::vector<double> theta(M, 0.0);
    std::vector<double> u0(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      u0[i] = std::pow(infoRates[i], Delta - 0.5);
    }

    auto f = [&](double aval)->double {
      for (size_t i = 0; i < kMax; ++i) {
        u[i] = aval * u0[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto probs = exitprob_seamless_cpp(M, r, theta, corr_known, k, u, infoRates);
      auto v = probs.get<std::vector<double>>("exitProb");
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    double cwt = brent(f, 0.0, 10.0, 1e-6);
    for (size_t i = 0; i < kMax; ++i) {
      criticalValues[i] = cwt * u0[i];
      if (!effStopping[i]) criticalValues[i] = 6.0;
    }
    return subset(criticalValues, 0, K);
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
    std::vector<double> u_vec; u_vec.reserve(kMax);
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

      // Ensure reusable buffers have size k1+1 and capacity >= K
      u_vec.resize(k1 + 1);

      // - copy already computed criticalValues[0..k1-1] into u_vec[0..k1-1]
      // the last entry (u_vec[k1]) will be set by the lambda
      std::memcpy(u_vec.data(), criticalValues.data(), k1 * sizeof(double));

      // Define lambda that only sets the last element of u_vec
      auto f = [&](double aval)->double {
        // set the last element to the current candidate critical value
        u_vec[k1] = aval;
        auto probs = exitprob_seamless_cpp(
          M, r, theta_vec, corr_known, k1, u_vec, infoRates);
        auto v = probs.get<std::vector<double>>("exitProb");
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

    return subset(criticalValues, 0, K);
  }

  throw std::invalid_argument("Invalid value for typeAlphaSpending");
}


//' @title Efficacy Boundaries for Two-Stage Seamless Sequential Design
//' @description Calculates the efficacy stopping boundaries for a two-stage
//' seamless sequential design, accounting for the selection of the best arm
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
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime) {

  // ----------- Input Validation ----------- //
  if (std::isnan(beta) && std::isnan(IMax)) {
    throw std::invalid_argument("beta and IMax cannot be missing simultaneously");
  }
  if (!std::isnan(beta) && !std::isnan(IMax)) {
    throw std::invalid_argument("Only one of beta and IMax should be provided");
  }
  if (!std::isnan(IMax) && IMax <= 0) {
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

  // Alpha and Beta must be within valid ranges
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(beta) && (beta >= 1 - alpha || beta < 0.0001)) {
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");
  }

  std::string unknown = std::isnan(beta) ? "beta" : "IMax";

  size_t kMax = K + 1;

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
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);

  if (!missingCriticalValues && criticalValues.size() != kMax) {
    throw std::invalid_argument("Invalid length for criticalValues");
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
    if (userAlphaSpending.back() != alpha)
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
    if (spendingTime.back() != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // ----------- End of Input Validation ----------- //

  // set up efficacy bounds
  std::vector<double> zero(M, 0.0);
  std::vector<double> critValues = criticalValues;
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
      std::vector<double> u(kMax);
      for (size_t i = 0; i < kMax - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        ListCpp probs = exitprob_seamless_cpp(
          M, r, zero, corr_known, K, u, infoRates);
        auto v = probs.get<std::vector<double>>("exitProb");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, 0.0, 6.0, 1e-6);
    } else {
      critValues = getBound_seamless_cpp(
        M, r, corr_known, K, infoRates, alpha, asf, parameterAlphaSpending,
        userAlphaSpending, spendTime, effStopping);
    }
  }

  ListCpp probs = exitprob_seamless_cpp(
    M, r, zero, corr_known, K, critValues, infoRates);
  auto v = probs.get<std::vector<double>>("exitProb");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent[kMax-1];

  double IMax1 = IMax;
  if (unknown == "IMax") {
    double maxtheta = *std::max_element(theta.begin(), theta.end());
    auto f = [&](double aval)->double {
      std::vector<double> theta1(M);
      for (size_t i = 0; i < M; ++i) theta1[i] = theta[i] * aval / maxtheta;
      ListCpp probs = exitprob_seamless_cpp(
        M, r, theta1, true, K, critValues, infoRates);
      auto v = probs.get<std::vector<double>>("exitProb");
      double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
      return overallReject - (1.0 - beta);
    };

    double drift = brent(f, 0.0, 6.0, 1e-6);
    IMax1 = sq(drift / maxtheta);
    std::vector<double> info(kMax);
    for (size_t i = 0; i < kMax; ++i) info[i] = infoRates[i] * IMax1;
    probs = exitprob_seamless_cpp(M, r, theta, true, K, critValues, info);
  } else {
    std::vector<double> info = infoRates;
    for (double &x : info) x *= IMax1;
    probs = exitprob_seamless_cpp(M, r, theta, true, K, critValues, info);
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
  auto pu = probs.get<std::vector<double>>("exitProb");
  std::vector<double> cpu(kMax);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  double overallReject = cpu[kMax-1];

  for (size_t i = 0; i < kMax; ++i) {
    if (critValues[i] == 6.0) effStopping[i] = 0;
  }

  auto selectAsBest = probs.get<std::vector<double>>("selectAsBest");
  auto exitProbByArm = probs.get<FlatMatrix>("exitProbByArm");
  std::vector<double> powerByArm(M, 0.0);
  for (size_t j = 0; j < M; ++j) {
    double p = 0.0;
    for (size_t i = 0; i < kMax; ++i) {
      p += exitProbByArm(i, j);
    }
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
  overallResults.push_back(M, "M");
  overallResults.push_back(r, "r");
  overallResults.push_back(corr_known, "corr_known");
  overallResults.push_back(K, "K");
  overallResults.push_back(IMax1, "information");

  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(pu), "rejectPerStage");
  byStageResults.push_back(std::move(cpu), "cumulativeRejection");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(efficacyTheta), "efficacyTheta");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(information), "information");
  byStageResults.push_back(std::move(effStopping), "efficacyStopping");

  DataFrameCpp byArmResults;
  byArmResults.push_back(theta, "theta");
  byArmResults.push_back(std::move(selectAsBest), "selectAsBest");
  byArmResults.push_back(std::move(powerByArm), "powerByArm");
  byArmResults.push_back(std::move(condPowerByArm), "condPowerByArm");

  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(spendingTime, "spendingTime");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(byArmResults), "byArmResults");
  result.push_back(std::move(settings), "settings");
  return result;
}


//' @title Power and Sample Size for Two-Stage Seamless Sequential Design
//' @description Computes either the maximum information and stopping
//' boundaries for a generic two-stage seamless sequential design, or
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
//' @param M Number of active treatment arms in Phase 2.
//' @param r Randomization ratio of each active arm to the common control
//'   in Phase 2.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics in Phase 2 is derived from the randomization ratio \code{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param K Number of sequential looks in Phase 3.
//' @param informationRates A numeric vector of information rates fixed
//'   before the trial. If unspecified, defaults to \eqn{(1:(K+1)) / (K+1)}.
//' @inheritParams param_efficacyStopping
//' @param criticalValues The upper boundaries on the max z-test statistic
//'   scale for Phase 2 and the z-test statistics for the selected arm
//'   in Phase 3 for the primary trial. If missing, boundaries
//'   will be computed based on the specified alpha spending function.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A numeric vector of length \eqn{K+1} specifying the
//'   error spending time at each analysis. Values must be strictly increasing
//'   and ends at 1. If omitted, defaults to \code{informationRates}.
//'
//' @return An S3 object of class \code{seamless} with these components:
//'
//' * \code{overallResults}: A data frame containing:
//'     - \code{overallReject}: Overall probability of rejecting the null
//'       hypothesis.
//'     - \code{alpha}: Overall significance level.
//'     - \code{M}: Number of active arms in phase 2.
//'     - \code{r}: Randomization ratio per active arm versus control in phase 2.
//'     - \code{corr_known}: Whether the phase-2 correlation was assumed known.
//'     - \code{K}: Number of stages in phase 3.
//'     - \code{information}: Maximum information for any active arm versus control.
//'
//' * \code{byStageResults}: A data frame containing:
//'     - \code{informationRates}: Information rates at each analysis.
//'     - \code{efficacyBounds}: Efficacy boundaries on the Z-scale.
//'     - \code{rejectPerStage}: Probability of efficacy stopping at each stage.
//'     - \code{cumulativeRejection}: Cumulative probability of efficacy stopping.
//'     - \code{cumulativeAlphaSpent}: Cumulative alpha spent.
//'     - \code{efficacyTheta}: Efficacy boundaries on the parameter scale.
//'     - \code{efficacyP}: Efficacy boundaries on the p-value scale.
//'     - \code{information}: Cumulative information for any active arm versus
//'       control at each analysis.
//'     - \code{efficacyStopping}: Indicator of whether efficacy stopping
//'       is permitted.
//'
//' * \code{byArmResults}: A data frame containing:
//'     - \code{theta}: Parameter values for the active arms.
//'     - \code{selectAsBest}: Probability an arm is selected as best in at
//'       the end of phase 2.
//'     - \code{powerByArm}: Probability of rejecting the null for each arm
//'       by trial end.
//'     - \code{condPowerByArm}: Conditional power for each arm given it was
//'       selected as best at the end of phase 2.
//'
//' * \code{settings}: A list of input settings:
//'     - \code{typeAlphaSpending}: Type of alpha spending function.
//'     - \code{parameterAlphaSpending}: Parameter value for the chosen
//'       alpha spending function.
//'     - \code{userAlphaSpending}: User-specified alpha spending values.
//'     - \code{spendingTime}: Error-spending times at each analysis.
//'
//' @details If \code{corr_known} is \code{FALSE}, critical boundaries are
//' computed assuming independence among the phase-2 Wald statistics
//' (a conservative assumption). Power calculations, however, use the
//' correlation implied by the randomization ratio \eqn{r}.
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
//' # Example 1: obtain the maximum information given power
//' (design1 <- getDesign_seamless(
//'   beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   K = 2, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesign_seamless(
//'   IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   K = 2, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF"))
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
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {

  std::vector<double> thetaVec(theta.begin(), theta.end());
  std::vector<double> infoRates(informationRates.begin(), informationRates.end());
  auto effStopping = convertLogicalVector(efficacyStopping);
  std::vector<double> critValues(criticalValues.begin(), criticalValues.end());
  std::vector<double> userAlpha(userAlphaSpending.begin(), userAlphaSpending.end());
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());

  auto cpp_result = getDesign_seamless_cpp(
    beta, IMax, thetaVec, static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(K), infoRates, effStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, spendTime
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
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
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
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(theta)) throw std::invalid_argument("theta must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (K <= L) throw std::invalid_argument("K must be greater than L");

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(betaNew) && (betaNew < 0.0001 || betaNew >= 1)) {
    throw std::invalid_argument("betaNew must lie in [0.0001, 1)");
  }

  size_t kMax = K + 1;

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
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);

  if (!missingCriticalValues && criticalValues.size() != kMax) {
    throw std::invalid_argument("Invalid length for criticalValues");
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
    if (userAlphaSpending.back() != alpha)
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
    if (spendingTime.back() != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // ----------- New Design Input Validation ----------- //
  std::vector<double> infoRatesNew;
  std::vector<unsigned char> effStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

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
      if (informationRatesNew.back() != 1.0)
        throw std::invalid_argument("informationRatesNew must end with 1");
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
  // ----------- End of Input Validation ----------- //

  // set up efficacy bounds
  std::vector<double> zero(M, 0.0);
  std::vector<double> critValues = criticalValues;
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
      std::vector<double> u(kMax);
      for (size_t i = 0; i < kMax - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        ListCpp probs = exitprob_seamless_cpp(
          M, r, zero, corr_known, K, u, infoRates);
        auto v = probs.get<std::vector<double>>("exitProb");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, 0.0, 6.0, 1e-6);
    } else {
      critValues = getBound_seamless_cpp(
        M, r, corr_known, K, infoRates, alpha, asf, parameterAlphaSpending,
        userAlphaSpending, spendTime, effStopping);
    }
  }

  ListCpp probs = exitprob_seamless_cpp(
    M, r, zero, corr_known, K, critValues, infoRates);
  auto v = probs.get<std::vector<double>>("exitProb");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent[kMax-1];

  size_t k1 = K - L;
  std::vector<double> t1(k1), b1(k1), a1(k1, -6.0), zero1(k1, 0.0);
  for (size_t l = 0; l < k1; ++l) {
    t1[l] = (infoRates[l + L + 1] - infoRates[L]) / (1.0 - infoRates[L]);
    double r1 = infoRates[L] / infoRates[l + L + 1];
    b1[l] = (critValues[l + L + 1] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
    if (!effStopping[l + L + 1]) b1[l] = 6.0;
  }

  a1[k1-1] = b1[k1-1];

  // conditional type I error
  ListCpp probs0 = exitprobcpp(b1, a1, zero1, t1);
  auto v0 = probs0.get<std::vector<double>>("exitProbUpper");
  double alphaNew = std::accumulate(v0.begin(), v0.end(), 0.0);


  std::vector<double> I1(k1);
  for (size_t l = 0; l < k1; ++l) {
    I1[l] = IMax * (infoRates[l + L + 1] - infoRates[L]);
  }

  std::vector<double> theta1(k1, theta);
  ListCpp probs1 = exitprobcpp(b1, a1, theta1, I1);
  auto v1 = probs1.get<std::vector<double>>("exitProbUpper");
  double conditionalPower = std::accumulate(v1.begin(), v1.end(), 0.0);

  // information for the primary trial
  std::vector<double> information1(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    information1[i] = IMax * infoRates[i];
  }

  ListCpp des2;
  if (!MullerSchafer) {
    effStoppingNew = subset(effStopping, L + 1, kMax);
    std::vector<unsigned char> futStoppingNew(k1, 1);
    std::string typeBetaSpendingNew = "none";
    double parameterBetaSpendingNew = NaN;
    std::vector<double> userBetaSpendingNew(1, NaN);
    double varianceRatio = 1.0;

    des2 = getDesigncpp(betaNew, INew, theta, k1, t1, effStoppingNew,
                        futStoppingNew, b1, NaN, typeAlphaSpendingNew,
                        parameterAlphaSpendingNew, {NaN}, a1,
                        typeBetaSpendingNew, parameterBetaSpendingNew,
                        userBetaSpendingNew, spendTimeNew, varianceRatio);
  } else {
    if (!std::isnan(betaNew) && betaNew >= 1.0 - alphaNew) {
      throw std::invalid_argument(
          "betaNew must be less than 1 minus the conditional type I error");
    }

    std::vector<unsigned char> futStoppingNew(kNew, 1);
    std::string typeBetaSpendingNew = "none";
    double parameterBetaSpendingNew = NaN;
    std::vector<double> userBetaSpendingNew(1, NaN);
    double varianceRatio = 1.0;

    std::vector<double> b1New(kNew, NaN), a1New(kNew, NaN);
    des2 = getDesigncpp(betaNew, INew, theta, kNew, infoRatesNew, effStoppingNew,
                        futStoppingNew, b1New, alphaNew, typeAlphaSpendingNew,
                        parameterAlphaSpendingNew, {NaN}, a1New,
                        typeBetaSpendingNew, parameterBetaSpendingNew,
                        userBetaSpendingNew, spendTimeNew, varianceRatio);
  }


  // combined design
  DataFrameCpp byStage2 = des2.get<DataFrameCpp>("byStageResults");
  std::vector<double> I2 = byStage2.get<double>("information");
  std::vector<double> b2 = byStage2.get<double>("efficacyBounds");
  std::vector<double> a2 = byStage2.get<double>("futilityBounds");

  size_t k2 = I2.size();
  size_t kc = L + 1 + k2;

  double IL = information1[L];

  std::vector<double> Ic(kc); // cumulative information for the combined design
  for (size_t i = 0; i <= L; ++i) {
    Ic[i] = information1[i];
  }
  for (size_t i = L + 1; i < kc; ++i) {
    Ic[i] = IL + I2[i - L - 1];
  }
  double IMaxc = Ic[kc - 1];

  std::vector<double> tc(kc); // information rates for the combined design
  for (size_t i = 0; i < kc; ++i) {
    tc[i] = Ic[i] / IMaxc;
  }

  std::vector<double> bc(kc); // efficacy bounds for the combined design
  for (size_t i = 0; i <= L; ++i) {
    bc[i] = critValues[i];
  }
  for (size_t i = L + 1; i < kc; ++i) {
    if (b2[i - L - 1] != 6.0) {
      double r1 = IL / Ic[i];
      bc[i] = zL * std::sqrt(r1) + b2[i - L - 1] * std::sqrt(1.0 - r1);
    } else {
      bc[i] = 6.0;
    }
  }

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
  des1.push_back(std::move(information1), "information");
  des1.push_back(alpha1, "alpha");
  des1.push_back(alphaNew, "conditionalAlpha");
  des1.push_back(conditionalPower, "conditionalPower");
  des1.push_back(MullerSchafer, "MullerSchafer");

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
  des3.push_back(std::move(tc), "informationRates");
  des3.push_back(std::move(bc), "efficacyBounds");
  des3.push_back(std::move(Ic), "information");

  ListCpp result;
  result.push_back(std::move(des1), "primaryTrial");
  result.push_back(std::move(des2), "secondaryTrial");
  result.push_back(std::move(des3), "integratedTrial");
  return result;
}


//' @title Adaptive Two-Stage Seamless Sequential Design
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
//' @param INew The maximum information for any active arm versus the common
//'   control in the secondary trial. Either
//'   \code{betaNew} or \code{INew} should be provided, while the other
//'   must be missing.
//' @param M Number of active treatment arms in Phase 2.
//' @param r Randomization ratio of each active arm to the common control
//'   in Phase 2.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics in Phase 2 is derived from the randomization ratio \code{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param L The interim adaptation look in Phase 3.
//' @param zL The z-test statistic at the interim adaptation look of
//'   Phase 3.
//' @param theta The treatment effect for the selected arm versus the
//'   common control.
//' @param IMax Maximum information for any active arm versus the common
//'   control for the original trial. Must be provided.
//' @param K Number of sequential looks in Phase 3.
//' @param informationRates A numeric vector of information rates fixed
//'   before the trial. If unspecified, defaults to \eqn{(1:(K+1)) / (K+1)}.
//' @param efficacyStopping Indicators of whether efficacy stopping is
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
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRatesNew}.
//'
//' @return An \code{adaptDesign_seamless} object with three list components:
//'
//' * \code{primaryTrial}: A list of selected information for the primary
//'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{K},
//'   \code{L}, \code{zL}, \code{theta}, \code{maxInformation}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, \code{information},
//'   \code{alpha}, \code{conditionalAlpha}, \code{conditionalPower},
//'   and \code{MullerSchafer}.
//'
//' * \code{secondaryTrial}: A \code{design} object for the secondary trial.
//'
//' * \code{integratedTrial}: A list of selected information for the integrated
//'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{K},
//'   \code{L}, \code{zL}, \code{theta}, \code{maxInformation}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, and \code{information}.
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
    const double r = NA_REAL,
    const bool corr_known = true,
    const int L = NA_INTEGER,
    const double zL = NA_REAL,
    const double theta = NA_REAL,
    const double IMax = NA_REAL,
    const int K = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int kNew = NA_INTEGER,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  auto cpp_result = adaptDesign_seamless_cpp(
    betaNew, INew, static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(L), zL, theta, IMax,
    static_cast<size_t>(K), infoRates, effStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, spendTime,
    MullerSchafer, static_cast<size_t>(kNew), infoRatesNew, effStoppingNew,
    typeAlphaSpendingNew, parameterAlphaSpendingNew, spendTimeNew
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "adaptDesign_seamless";
  return result;
}

