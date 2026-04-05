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


ListCpp exitprob_tsssd_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t K,
    const std::vector<double>& b,
    const std::vector<double>& I) {

  if (M < 2) {
    throw std::invalid_argument("M should be at least 2");
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

  PMVNResult out0 = pmvnormcpp(lower0, upper0, mean0, sigma0);
  preject[0] = 1.0 - out0.prob;

  // compute the exit probabilities at the K looks in phase 3
  std::vector<double> I1(K), sqrtI1(K);
  for (size_t k = 0; k < K; ++k) {
    I1[k] = Ivec[k+1] - Ivec[0]; // information for the selected arm in phase 3 at look k
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

      PMVNResult out = pmvnormcpp(lower, upper, mean, sigma);
      // multiply the density of the Wald statistic for arm m at z0
      return out.prob * boost_dnorm(z0, mu, 1.0);
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


//' @title Exit probabilities for two-stage seamless sequential design (TSSSD)
//' @description Computes the exit probabilities for a two-stage selection
//' and testing design with a common control arm in phase 2 and the selected
//' active arm and the common control arm in phase 3.
//'
//' @param M number of active arms in phase 2.
//' @param r randomization ratio of each active arm to the common control in phase 2.
//' @param theta vector of true treatment effects for the M active arms in phase 2.
//' @param corr_known whether the correlation between the Wald statistics in
//'   phase 2 is known. If \code{TRUE}, the correlation is determined by
//'   the randomization ratio \code{r} as \eqn{r / (r + 1)};
//'   if \code{FALSE}, the correlation is conservatively set to 0.
//' @param K number of looks in phase 3.
//' @param b vector of critical values for phase 2 and the K looks in phase 3.
//' @param I vector of information levels for phase 2 and the K looks in
//'         phase 3 for any arm versus the common control.
//'
//' @return A list with three components:
//'
//' * exitProb: A vector of length \eqn{K + 1}, where the first element is the
//' exit probability in phase 2, and the \eqn{k}-th element for
//' \eqn{k = 1,\ldots,K} is the exit probability at look \eqn{k} in phase 3.
//'
//' * exitProbByArm: A matrix of dimension \eqn{(K + 1)} by \eqn{M}, where
//' the \eqn{(k, m)}-th  element for \eqn{k = 0,\ldots,K} and
//' \eqn{m = 1,\ldots,M} is the probability of rejecting the
//' global null hypothesis in phase 2 if arm \eqn{m} is selected in phase 2
//' when \eqn{k = 0}, and the probability of rejecting the global null hypothesis
//' at look \eqn{k} in phase 3 if arm \eqn{m} is selected in phase 2 when
//' \eqn{k = 1,...,K}.
//'
//' * selectAsBest: A vector of length \eqn{M}, where the \eqn{m}-th element
//' is the probability of selecting arm \eqn{m} in phase 2 as the best arm
//' among the \eqn{M} active arms in phase 2.
//'
//' @details The exit probability in phase 2 is the probability of rejecting
//' the global null hypothesis at the end of phase 2, and the exit probability
//' at look \eqn{k} in phase 3 is the probability of rejecting the global null
//' hypothesis at look \eqn{k} in phase 3 given that the global null hypothesis
//' is not rejected in phase 2. The exit probabilities are computed under
//' the multivariate normal distribution of the Wald statistics for the \eqn{M}
//' active arms in phase 2 and the selected active arm in phase 3, where
//' the covariance matrix is determined by the randomization ratio \eqn{r}.
//'
//' We reject the global null hypothesis in phase 2 if \deqn{Z(I_0) >= b_0,}
//' where \eqn{Z(I_0)} is the Wald statistic for best arm in phase 2, and
//' we reject the global null hypothesis at look \eqn{k} in phase 3 if
//' \deqn{Z(I_0) < b_0, Z(I_j) < b_j \text{ for } j = 1,\ldots,k-1,
//' \text{ and } Z(I_k) >= b_k,} where \eqn{Z(I_k)} is the Wald statistic
//' for the selected arm at look \eqn{k} in phase 3.
//'
//' The following assumptions are made for the TSSSD design:
//'
//' * Information is the same across active arms in phase 2
//'
//' * At most one active arm can be selected at the end of phase 2
//'
//' * The selected arm is the one with the largest score statistic in phase 2
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' # Two active arms vs. a common control in phase 2, and the selected
//' # active arm vs. the common control in phase 3 with 2 looks in phase 3
//'
//' # equal spacing of information levels across phase 2/3 for any arm vs.
//' # the common control
//'
//' # maximum of 110 patients per arm planned
//' # standard deviation of the endpoint is 1
//' I <- c(110 / (2 * 1.0^2) * seq(1, 3)/3)
//'
//' # O'Brien-Fleming critical values for phase 2/3
//' b <- c(3.776606, 2.670463, 2.180424)
//'
//' # Type I error under the global null hypothesis
//' p0 <- exitprob_tsssd(M = 2, theta = c(0, 0), K = 2, b = b, I = I)
//' cumsum(p0$exitProb)
//'
//' # Power under the alternative hypothesis where the treatment effects for
//' # the two active arms are 0.3 and 0.5, respectively
//' p1 <- exitprob_tsssd(M = 2, theta = c(0.3, 0.5), K = 2, b = b, I = I)
//' cumsum(p1$exitProb)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List exitprob_tsssd(
    const int M = NA_INTEGER,
    const double r = 1.0,
    const Rcpp::NumericVector& theta = NA_REAL,
    const bool corr_known = true,
    const int K = NA_INTEGER,
    const Rcpp::NumericVector& b = NA_REAL,
    const Rcpp::NumericVector& I = NA_REAL) {

  std::vector<double> thetaVec(theta.begin(), theta.end());
  std::vector<double> bVec(b.begin(), b.end());
  std::vector<double> IVec(I.begin(), I.end());
  auto result = exitprob_tsssd_cpp(static_cast<size_t>(M), r, thetaVec, corr_known,
                                   static_cast<size_t>(K), bVec, IVec);
  return Rcpp::wrap(result);
}


std::vector<double> getBound_tsssd_cpp(
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

  if (M < 2) {
    throw std::invalid_argument("M should be at least 2");
  }
  if (r <= 0.0) {
    throw std::invalid_argument("r should be positive");
  }
  if (k < 1) {
    throw std::invalid_argument("k should be at least 1");
  }

  size_t K = k + 1; // add the look at end of phase 2

  // infoRates: if missing create 1/K, 2/K, ..., k/K
  std::vector<double> infoRates(K);
  if (none_na(informationRates)) {
    if (informationRates.size() < K)
      throw std::invalid_argument("Insufficient length for informationRates");

    infoRates.assign(informationRates.begin(), informationRates.begin() + K);
    if (infoRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(infoRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (infoRates[K-1] > 1.0)
      throw std::invalid_argument("informationRates must not exceed 1");
  } else {
    for (size_t i = 0; i < K; ++i) {
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(K);
    }
  }

  // spendTime: default to infoRates if missing
  std::vector<double> spendTime; spendTime.reserve(K);
  if (none_na(spendingTime)) {
    if (spendingTime.size() < K)
      throw std::invalid_argument("Insufficient length for spendingTime");

    spendTime.assign(spendingTime.begin(), spendingTime.begin() + K);
    if (spendTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendTime[K-1] > 1.0)
      throw std::invalid_argument("spendingTime must not exceed 1");
  } else {
    spendTime = infoRates;
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping; effStopping.reserve(K);
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() < K)
      throw std::invalid_argument("Insufficient length for efficacyStopping");

    effStopping.assign(efficacyStopping.begin(), efficacyStopping.begin() + K);
    if (effStopping[K-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
  } else {
    effStopping.assign(K, 1);
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
    if (userAlphaSpending.size() < K)
      throw std::invalid_argument("Insufficient length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[K-1] > alpha)
      throw std::invalid_argument("userAlphaSpending must not exceed alpha");
  }

  if (asf == "of" || asf == "p" || asf == "wt") {
    bool equal_spacing = true;
    for (size_t i = 0; i < K; ++i) {
      double expected = static_cast<double>(i+1) / static_cast<double>(K);
      if (std::fabs(infoRates[i] - expected) > 1e-6) {
        equal_spacing = false; break;
      }
    }
    if (!equal_spacing) {
      thread_utils::push_thread_warning(
        "Equal spacing is used for OF, P, and WT boundaries");
    }
  }

  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }

  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  std::vector<double> criticalValues(K);

  if (asf == "none") {
    for (size_t i = 0; i < K-1; ++i) criticalValues[i] = 6.0;
    criticalValues[K-1] = boost_qnorm(1.0 - alpha);
    return criticalValues;
  }

  if (asf == "of" || asf == "p" || asf == "wt") {
    double Delta;
    if (asf == "of") Delta = 0.0;
    else if (asf == "p") Delta = 0.5;
    else Delta = parameterAlphaSpending; // parameterAlphaSpending holds delta for WT

    // for a given multiplier, compute cumulative upper exit probability - alpha
    std::vector<double> u(K);
    std::vector<double> theta(M, 0.0);
    std::vector<double> I(K), u0(K);
    for (size_t i = 0; i < K; ++i) {
      I[i] = static_cast<double>(i+1) / static_cast<double>(K);
      u0[i] = std::pow(I[i], Delta - 0.5);
    }

    auto f = [&](double aval)->double {
      for (size_t i = 0; i < K; ++i) {
        u[i] = aval * u0[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto probs = exitprob_tsssd_cpp(M, r, theta, corr_known, k, u, I);
      auto v = probs.get<std::vector<double>>("exitProb");
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    double cwt = brent(f, 0.0, 10.0, 1e-4);
    for (size_t i = 0; i < K; ++i) {
      criticalValues[i] = cwt * u0[i];
      if (!effStopping[i]) criticalValues[i] = 6.0;
    }
    return criticalValues;
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
      criticalValues[0] = qmvnormcpp(1.0 - cumAlpha, mean, sigma);
    }

    // Preallocate reusable buffers used by the root-finding lambda
    std::vector<double> u_vec; u_vec.reserve(K);
    std::vector<double> theta_vec(M, 0.0);

    // subsequent stages
    for (size_t k1 = 1; k1 < K; ++k1) {
      // determine cumulative alpha at this stage
      if (asf == "user") cumAlpha = userAlphaSpending[k1];
      else cumAlpha = errorSpentcpp(spendTime[k1], alpha, asf, parameterAlphaSpending);

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
        auto probs = exitprob_tsssd_cpp(M, r, theta_vec, corr_known,
                                        k1, u_vec, infoRates);
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

        criticalValues[k1] = brent(f_for_brent, 0.0, 6.0, 1e-4);
      }
    }

    return criticalValues;
  }

  throw std::invalid_argument("Invalid value for typeAlphaSpending");
}


//' @title Efficacy Boundaries for Two-Stage Seamless Sequential Design
//' @description Obtains the efficacy stopping boundaries for a two-stage
//' seamless sequential design (TSSSD).
//'
//' @param M Number of active arms in phase 2. Must be at least 2.
//' @param r Randomization ratio of each active arm to the common control in
//'  phase 2. Must be positive.
//' @param corr_known whether the correlation between the Wald statistics in
//'   phase 2 is known. If \code{TRUE}, the correlation is determined by
//'   the randomization ratio \code{r} as \eqn{r / (r + 1)};
//'   if \code{FALSE}, the correlation is conservatively set to 0.
//' @param k Look number for the current analysis.
//' @param informationRates Information rates up to the current look. Must be
//'   increasing and less than or equal to 1.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A vector of length \code{k} for the error spending
//'   time at each analysis. Must be increasing and less than or equal to 1.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @inheritParams param_efficacyStopping
//'
//' @details
//' If \code{typeAlphaSpending} is "OF", "P", or "WT", then the boundaries
//' will be based on equally spaced looks.
//'
//' @return A numeric vector of critical values up to the current look.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' getBound_tsssd(M = 2, k = 2, informationRates = seq(1,3)/3,
//'               alpha = 0.025, typeAlphaSpending = "OF")
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector getBound_tsssd(
    const int M = NA_INTEGER,
    const double r = 1.0,
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

  auto result = getBound_tsssd_cpp(
    static_cast<size_t>(M), r, corr_known, static_cast<size_t>(k),
    infoRates, alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlpha, spendTime, effStopping
  );

  return Rcpp::wrap(result);
}


ListCpp getDesign_tsssd_cpp(
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
    const std::vector<double>& spendingTime,
    const double varianceRatio) {

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
  if (M < 2) {
    throw std::invalid_argument("M must be at least 2");
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

  if (varianceRatio <= 0.0) {
    throw std::invalid_argument("varianceRatio must be positive");
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
        ListCpp probs = exitprob_tsssd_cpp(M, r, zero, corr_known,
                                           K, u, infoRates);
        auto v = probs.get<std::vector<double>>("exitProb");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, 0.0, 6.0, 1e-4);
    } else {
      critValues = getBound_tsssd_cpp(
        M, r, corr_known, K, infoRates, alpha, asf, parameterAlphaSpending,
        userAlphaSpending, spendTime, effStopping);
    }
  }

  ListCpp probs = exitprob_tsssd_cpp(M, r, zero, corr_known,
                                     K, critValues, infoRates);
  auto v = probs.get<std::vector<double>>("exitProb");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent[kMax-1];

  // multiplier for the boundaries under the alternative hypothesis
  std::vector<double> w(kMax, std::sqrt(varianceRatio));
  std::vector<double> u(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    u[i] = critValues[i] * w[i];
  }

  double IMax1 = IMax;
  if (unknown == "IMax") {
    double maxtheta = *std::max_element(theta.begin(), theta.end());
    auto f = [&](double aval)->double {
      std::vector<double> theta1(M);
      for (size_t i = 0; i < M; ++i) theta1[i] = theta[i] * aval / maxtheta;
      ListCpp probs = exitprob_tsssd_cpp(M, r, theta1, true,
                                         K, critValues, infoRates);
      auto v = probs.get<std::vector<double>>("exitProb");
      double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
      return overallReject - (1.0 - beta);
    };

    double drift = brent(f, 0.0, 6.0, 1e-4);
    IMax1 = sq(drift / maxtheta);
    std::vector<double> info(kMax);
    for (size_t i = 0; i < kMax; ++i) info[i] = infoRates[i] * IMax1;
    probs = exitprob_tsssd_cpp(M, r, theta, true, K, critValues, info);
  } else {
    std::vector<double> info = infoRates;
    for (double &x : info) x *= IMax1;
    probs = exitprob_tsssd_cpp(M, r, theta, true, K, critValues, info);
  }

  // output the results
  std::vector<double> information(kMax);
  std::vector<double> efficacyTheta(kMax);
  std::vector<double> efficacyP(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    information[i] = IMax1 * infoRates[i];
    efficacyTheta[i] = u[i] / std::sqrt(information[i]);
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
  }

  // stagewise exit probabilities under H1
  auto pu = probs.get<std::vector<double>>("exitProb");
  std::vector<double> cpu(kMax);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  double overallReject = cpu[kMax-1];

  for (size_t i = 0; i < kMax; ++i) {
    if (critValues[i] == 6) effStopping[i] = 0;
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
  overallResults.push_back(static_cast<int>(M), "M");
  overallResults.push_back(r, "r");
  overallResults.push_back(corr_known, "corr_known");
  overallResults.push_back(static_cast<int>(K), "K");
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
  settings.push_back(varianceRatio, "varianceRatio");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(byArmResults), "byArmResults");
  result.push_back(std::move(settings), "settings");
  return result;
}


//' @title Power and Sample Size for a Generic Two-Stage Seamless Sequential Design
//' @description Obtains the maximum information and stopping boundaries
//' for a generic two-stage seamless sequential design design, or obtains
//' the power given the maximum information and stopping boundaries.
//'
//' @param beta The type II error.
//' @param IMax The maximum information for any active arm versus the common
//'   control. Either \code{beta} or \code{IMax} should be provided while
//'   the other one should be missing.
//' @param theta The parameter vector. The global null hypothesis is at
//'   \eqn{\theta_i = 0} for all \eqn{i}, and the alternative hypothesis
//'   is one-sided for \eqn{theta_i > 0} for at least one \eqn{i = 1,\ldots,M}.
//' @param M The number of active arms in phase 2. Must be at least 2.
//' @param r The randomization ratio of each active arm to the common control
//'  in phase 2. Must be positive.
//' @param corr_known whether the correlation between the Wald statistics in
//'   phase 2 is known. If \code{TRUE}, the correlation is determined by
//'   the randomization ratio \code{r} as \eqn{r / (r + 1)};
//'   if \code{FALSE}, the correlation is conservatively set to 0.
//' @param K The maximum number of stages in phase 3. Must be at least 1.
//' @param informationRates The information rates. Fixed prior to the trial.
//'   Defaults to \code{(1:(K+1)) / (K+1)} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A vector of length \code{K+1} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param varianceRatio The ratio of the variance under H0 to the
//'   variance under H1.
//'
//' @return An S3 class \code{design} object with three components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{M}: The number of active arms in phase 2.
//'
//'     - \code{r}: The randomization ratio of each active arm to the
//'       common control in phase 2.
//'
//'     - \code{corr_known}: Whether the correlation between the Wald statistics
//'       in phase 2 is known.
//'
//'     - \code{K}: The number of stages in phase 3.
//'
//'     - \code{information}: The maximum information for any given active arm
//'       versus the common control.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'
//'     - \code{efficacyTheta}: The efficacy boundaries on the parameter
//'       scale.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'
//' * \code{byArmResults}: A data frame containing the following variables:
//'
//'     - \code{theta}: The parameter values for the active arms.
//'
//'     - \code{selectAsBest}: The probability of being selected as the
//'       best arm at the end of phase 2.
//'
//'     - \code{powerByArm}: The probability of rejecting the null hypothesis
//'       for each arm at the end of the trial.
//'
//'     - \code{condPowerByArm}: The conditional power for each arm given that
//'       the arm is selected as the best at the end of phase 2.
//'
//' * \code{settings}: A list containing the following input parameters:
//'
//'     - \code{typeAlphaSpending}: The type of alpha spending.
//'
//'     - \code{parameterAlphaSpending}: The parameter value for alpha
//'       spending.
//'
//'     - \code{userAlphaSpending}: The user defined alpha spending.
//'
//'     - \code{spendingTime}: The error spending time at each analysis.
//'
//'     - \code{varianceRatio}: The ratio of the variance under H0
//'       to the variance under H1.
//'
//' @details When \code{corr_known} is \code{FALSE}, the critical boundaries
//' are computed under the assumption of independence among the Wald
//' statistics in phase 2, whereas the power is evaluated using the
//' correlation implied by the randomization ratio.
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
//' (design1 <- getDesign_tsssd(
//'   beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   K = 2, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesign_tsssd(
//'   IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   K = 2, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getDesign_tsssd(
    const double beta = NA_REAL,
    const double IMax = NA_REAL,
    const Rcpp::NumericVector& theta = NA_REAL,
    const int M = NA_INTEGER,
    const double r = 1.0,
    const bool corr_known = true,
    const int K = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const double varianceRatio = 1) {

  std::vector<double> thetaVec(theta.begin(), theta.end());
  std::vector<double> infoRates(informationRates.begin(), informationRates.end());
  auto effStopping = convertLogicalVector(efficacyStopping);
  std::vector<double> critValues(criticalValues.begin(), criticalValues.end());
  std::vector<double> userAlpha(userAlphaSpending.begin(), userAlphaSpending.end());
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());

  auto cpp_result = getDesign_tsssd_cpp(
    beta, IMax, thetaVec, static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(K), infoRates, effStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha,
    spendTime, varianceRatio
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "tsssd";
  return result;
}
