#include "generic_design.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <algorithm>     // fill, lower_bound, max, min, upper_bound
#include <cctype>        // tolower
#include <cmath>         // exp, fabs, isnan, llround, log, sqrt
#include <cstddef>       // size_t
#include <cstdint>       // int64_t
#include <cstring>       // memcpy
#include <iterator>      // begin, end
#include <mutex>         // lock_guard, mutex
#include <numeric>       // accumulate, iota
#include <stdexcept>     // invalid_argument
#include <string>        // string
#include <utility>       // move
#include <vector>        // vector

#include <Rcpp.h>

using std::size_t;


double errorSpentcpp(const double t,
                     const double error,
                     const std::string& sf,
                     const double sfpar) {
  if (error <= 0 || error >= 1) {
    throw std::invalid_argument("error must be a number between 0 and 1");
  }
  if (t <= 0 || t > 1) {
    throw std::invalid_argument("t must be a number between 0 and 1");
  }

  std::string asf = sf;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  double aval;
  if (asf == "sfp") {
    aval = error * std::log(1.0 + (std::exp(1.0) - 1.0) * t);
  } else if (asf == "sfof") {
    aval = boost_qnorm(1.0 - error / 2.0);
    aval = 2.0 * (1.0 - boost_pnorm(aval / std::sqrt(t)));
  } else if (asf == "sfkd") {
    if (std::isnan(sfpar)) {
      throw std::invalid_argument("sfpar is missing for sfKD");
    } else if (sfpar <= 0) {
      throw std::invalid_argument("sfpar must be positive for sfKD");
    } else {
      aval = error * std::pow(t, sfpar);
    }
  } else if (asf == "sfhsd") {
    if (std::isnan(sfpar)) {
      throw std::invalid_argument("sfpar is missing for sfHSD");
    } else if (sfpar == 0) {
      aval = error * t;
    } else {
      aval = error * (1.0 - std::exp(-sfpar * t)) / (1.0 - std::exp(-sfpar));
    }
  } else {
    throw std::invalid_argument("Invalid spending function");
  }
  return aval;
}


//' @title Error Spending
//' @description Obtains the error spent at given spending times
//' for the specified error spending function.
//'
//' @param t A vector of spending times, typically equal to information
//'   fractions.
//' @param error The total error to spend.
//' @param sf The spending function. One of the following:
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function, and
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function.
//'   Defaults to \code{"sfOF"}.
//' @param sfpar The parameter for the spending function. Corresponds to
//'   \eqn{\rho} for \code{"sfKD"} and \eqn{\gamma} for \code{"sfHSD"}.
//'
//' @details
//' This function implements a variety of error spending functions commonly
//' used in group sequential designs, assuming one-sided hypothesis testing.
//'
//' **O'Brien-Fleming-Type Spending Function**
//'
//' This spending function allocates very little alpha early on and more alpha
//' later in the trial. It is defined as:
//' \deqn{
//' \alpha(t) = 2 - 2\Phi\left(\frac{z_{\alpha/2}}{\sqrt{t}}\right),
//' }
//' where \eqn{\Phi} is the standard normal cumulative distribution function,
//' \eqn{z_{\alpha/2}} is the critical value from the standard normal
//' distribution, and \eqn{t \in [0, 1]} denotes the information fraction.
//'
//' **Pocock-Type Spending Function**
//'
//' This function spends alpha more evenly throughout the study:
//' \deqn{
//' \alpha(t) = \alpha \log(1 + (e - 1)t),
//' }
//' where \eqn{e} is Euler's number (approximately 2.718).
//'
//' **Kim and DeMets Power-Type Spending Function**
//'
//' This family of spending functions is defined as:
//' \deqn{
//' \alpha(t) = \alpha t^{\rho}, \quad \rho > 0.
//' }
//' - When \eqn{\rho = 1}, the function mimics Pocock-type boundaries.
//' - When \eqn{\rho = 3}, it approximates O’Brien-Fleming-type boundaries.
//'
//' **Hwang, Shih, and DeCani Spending Function**
//'
//' This flexible family of functions is given by:
//' \deqn{
//' \alpha(t) =
//' \begin{cases}
//' \alpha \frac{1 - e^{-\gamma t}}{1 - e^{-\gamma}}, & \text{if }
//' \gamma \ne 0 \\ \alpha t, & \text{if } \gamma = 0.
//' \end{cases}
//' }
//' - When \eqn{\gamma = -4}, the spending function resembles
//'   O’Brien-Fleming boundaries.
//' - When \eqn{\gamma = 1}, it resembles Pocock boundaries.
//'
//' @return A vector of errors spent up to the interim look.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' errorSpent(t = 0.5, error = 0.025, sf = "sfOF")
//'
//' errorSpent(t = c(0.5, 0.75, 1), error = 0.025, sf = "sfHSD", sfpar = -4)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> errorSpent(
    const std::vector<double>& t,
    const double error = 0.025,
    const std::string& sf = "sfOF",
    const double sfpar = NA_REAL) {

  size_t n = t.size();
  std::vector<double> result(n);
  for (size_t i = 0; i < n; ++i) {
    result[i] = errorSpentcpp(t[i], error, sf, sfpar);
  }
  return result;
}


ListCpp exitprobcpp(const std::vector<double>& b,
                    const std::vector<double>& a,
                    const std::vector<double>& theta,
                    const std::vector<double>& I) {

  // K is the total number of stages
  if (!none_na(b)) {
    throw std::invalid_argument("b must be provided");
  }
  size_t K = b.size();

  // Integer value controlling grid for numerical integration as in
  // Jennison and Turnbull (2000)
  size_t r = 18;
  size_t r1 = 6 * r - 1;   // size of x1 and shift
  size_t r2 = 12 * r - 3;  // size of z, h vectors when fully expanded

  // Prepare a1, theta1, I1 only if defaults / expansion necessary.
  std::vector<double> a1; a1.reserve(K);
  if (none_na(a)) {
    if (a.size() < K) throw std::invalid_argument("Insufficient length for a");
    a1.assign(a.begin(), a.begin() + K);
  } else {
    a1.assign(K, -6.0);
    a1[K - 1] = b[K - 1]; // set last element to b[K-1]
  }

  // check lower < upper
  for (size_t i = 0; i < K; ++i) {
    if (a1[i] > b[i]) throw std::invalid_argument(
        "Lower bounds (a) must be less than or equal to upper bounds (b)");
  }

  // theta expansion
  std::vector<double> theta1; theta1.reserve(K);
  if (none_na(theta)) {
    if (theta.size() == 1) {
      theta1.assign(K, theta[0]);
    } else if (theta.size() < K) {
      throw std::invalid_argument("Insufficient length for theta");
    } else {
      theta1.assign(theta.begin(), theta.begin() + K);
    }
  } else {
    theta1.assign(K, 0.0);
  }

  // information times expansion / validation
  std::vector<double> I1; I1.reserve(K);
  if (none_na(I)) {
    if (I.size() < K) throw std::invalid_argument("Insufficient length for I");

    I1.assign(I.begin(), I.begin() + K);
    if (I1[0] <= 0.0) throw std::invalid_argument("I must be positive");
    if (any_nonincreasing(I1)) throw std::invalid_argument("I must be increasing");
  } else {
    I1.resize(K);
    std::iota(I1.begin(), I1.end(), 1.0);
  }

  // Precompute shifts (constant across stages)
  std::vector<double> shift(r1);
  for (size_t i = 0; i < r1; ++i) {
    if (i < r - 1) {
      shift[i] = -3.0 - 4.0 * std::log(static_cast<double>(r) / (i + 1.0));
    } else if (i < 5 * r) {
      shift[i] = -3.0 + 3.0 * ( (i + 1.0 - r) / (2.0 * r) );
    } else {
      shift[i] = 3.0 + 4.0 * std::log(static_cast<double>(r) / (6.0 * r - i - 1.0));
    }
  }


  // Precompute sqrt and theta*sqrt/I combos
  std::vector<double> sqrtI(K), thetaSqrtI(K), thetaI(K);
  for (size_t j = 0; j < K; ++j) {
    sqrtI[j] = std::sqrt(I1[j]);
    thetaSqrtI[j] = theta1[j] * sqrtI[j];
    thetaI[j] = theta1[j] * I1[j];
  }


  // dI and dThetaI
  std::vector<double> dI(K), dThetaI(K);
  dI[0] = I1[0];
  dThetaI[0] = thetaI[0];
  for (size_t j = 1; j < K; ++j) {
    dI[j] = I1[j] - I1[j-1];
    dThetaI[j] = thetaI[j] - thetaI[j-1];
  }

  // pre-allocate buffers once, reuse across stages
  std::vector<double> exitProbUpper(K), exitProbLower(K);

  // allocate working arrays at max sizes once
  std::vector<double> x1(r1), x(r1);  // x1 is untrimmed, x is trimmed
  std::vector<double> w(r2), z(r2), h(r2); // z, h for the current stage
  std::vector<double> z0, h0; // z0, h0 for the previous stage
  z0.reserve(r2); h0.reserve(r2);

  size_t m0 = 0; // size of previous z0/h0
  // z0 and h0 are represented by the first m0 entries of z and h vectors.

  for (size_t j = 0; j < K; ++j) {
    double thetaSqrtIj = thetaSqrtI[j];
    double sqrtIj = sqrtI[j];
    double a1j = a1[j];
    double bj = b[j];
    double dThetaIj = dThetaI[j];
    double sqrtdIj = std::sqrt(dI[j]);
    double sqrtI1dIj = std::sqrt(I1[j] / dI[j]);

    // initialize x1 = thetaSqrtI[j] + shift
    for (size_t i = 0; i < r1; ++i) x1[i] = thetaSqrtIj + shift[i];

    // trim off x values outside (a1[j], b[j])
    // find trimming indices using binary search (x1 is sorted b/c shift is monotone)
    auto it1 = std::upper_bound(x1.begin(), x1.end(), a1[j]);
    size_t i1 = static_cast<size_t>(it1 - x1.begin()); // x1[i1-1] <= a1[j] < x1[i1]

    auto it2 = std::lower_bound(x1.begin(), x1.end(), b[j]);
    size_t i2 = static_cast<size_t>(it2 - x1.begin()); // x1[i2-1] < b[j] <= x1[i2]

    // m1 is number of retained x nodes after trimming
    size_t m1 = i2 - i1 + 2;

    // set first and last x values to the bounds
    x[0] = a1[j];
    x[m1 - 1] = b[j];

    // copy interior trimmed x1 values
    for (size_t i = 1; i < m1 - 1; ++i) x[i] = x1[i + i1 - 1];

    // derive z grid (odd + even interleaving)
    size_t m = 2 * m1 - 1;
    // odd points
    for (size_t i = 0; i < m1; ++i) z[2*i] = x[i];
    // even points (midpoints)
    for (size_t i = 0; i < m1 - 1; ++i) z[2*i + 1] = 0.5 * (z[2*i] + z[2*i + 2]);


    // weights w as Simpson-like composite rule (same formulas)
    // first weight
    w[0] = (z[2] - z[0]) / 6.0;

    // interior even indices (i = 2, 4, ..., 2*(m1-2))
    for (size_t i0 = 1; i0 < m1 - 1; ++i0) {
      size_t i = 2 * i0;
      w[i] = (z[i + 2] - z[i - 2]) / 6.0;
    }
    // interior odd indices (i = 1,3,...)
    for (size_t i0 = 1; i0 < m1; ++i0) {
      size_t i = 2 * i0 - 1;
      w[i] = 4.0 * (z[i + 1] - z[i - 1]) / 6.0;
    }
    // last weight
    w[m - 1] = (z[m - 1] - z[m - 3]) / 6.0;


    // first stage is easy
    if (j == 0) {
      // exit probabilities
      exitProbUpper[j] = boost_pnorm(-bj + thetaSqrtIj);
      exitProbLower[j] = boost_pnorm(a1j - thetaSqrtIj);

      // prepare m0, h0, z0 for the next stage
      if (K > 1) {
        m0 = m;
        h0.resize(m0);
        z0.resize(m0);
        for (size_t i = 0; i < m0; ++i) {
          h0[i] = w[i] * boost_dnorm(z[i] - thetaSqrtIj);
        }
        std::memcpy(z0.data(), z.data(), m0 * sizeof(double));
      }
    } else {
      double sqrtIjm1 = sqrtI[j-1];

      // calculate exit probabilities using h0 from the previous stage, (19.8)
      double sumUpper = 0.0, sumLower = 0.0;
      for (size_t i0 = 0; i0 < m0; ++i0) {
        double tupper = (z0[i0] * sqrtIjm1 - bj * sqrtIj + dThetaIj) / sqrtdIj;
        double tlower = (-z0[i0] * sqrtIjm1 + a1j * sqrtIj - dThetaIj) / sqrtdIj;
        sumUpper += h0[i0] * boost_pnorm(tupper);
        sumLower += h0[i0] * boost_pnorm(tlower);
      }
      exitProbUpper[j] = sumUpper;
      exitProbLower[j] = sumLower;

      // prepare m0, h0, z0 for the next stage, equation following (19.7)
      if (j < K-1) {
        for (size_t i = 0; i < m; ++i) {
          double sum = 0.0;
          for (size_t i0 = 0; i0 < m0; ++i0) {
            double t = (z[i] * sqrtIj - z0[i0] * sqrtIjm1 - dThetaIj) / sqrtdIj;
            sum += h0[i0] * boost_dnorm(t);
          }
          h[i] = sum * w[i] * sqrtI1dIj; // factors invariant to i0
        }

        m0 = m;
        h0.resize(m0);
        z0.resize(m0);
        std::memcpy(h0.data(), h.data(), m0 * sizeof(double));
        std::memcpy(z0.data(), z.data(), m0 * sizeof(double));
      }
    }
  }

  // return a list of stagewise exit probabilities
  ListCpp exitProb;
  exitProb.push_back(std::move(exitProbUpper), "exitProbUpper");
  exitProb.push_back(std::move(exitProbLower), "exitProbLower");
  return exitProb;
}


//' @title Stagewise Exit Probabilities
//' @description Obtains the stagewise exit probabilities for both efficacy
//' and futility stopping.
//'
//' @param b Upper boundaries on the z-test statistic scale.
//' @param a Lower boundaries on the z-test statistic scale. Defaults to
//'   \code{c(rep(-6.0, kMax-1), b[kMax])} if left unspecified, where
//'   \code{kMax = length(b)}.
//' @param theta Stagewise parameter of interest, e.g., \code{-U/V} for
//'   weighted log-rank test, where \code{U} is the mean and \code{V} is
//'   the variance of the weighted log-rank test score statistic at each
//'   stage. For proportional hazards and conventional log-rank test, use the
//'   scalar input, \code{theta = -log(HR)}. Defaults to 0 corresponding to
//'   the null hypothesis.
//' @param I Stagewise cumulative information, e.g., \code{V}, the variance
//'   of the weighted log-rank test score statistic at each stage. For
//'   conventional log-rank test, information can be approximated by
//'   \code{phi*(1-phi)*D}, where \code{phi} is the probability of being
//'   allocated to the active arm, and \code{D} is the total number of events
//'   at each stage. Defaults to \code{seq(1, kMax)} if left unspecified.
//'
//' @return A list of stagewise exit probabilities:
//'
//' * \code{exitProbUpper}: The vector of efficacy stopping probabilities
//'
//' * \code{exitProbLower}: The vector of futility stopping probabilities.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' exitprob(b = c(3.471, 2.454, 2.004), theta = -log(0.6),
//'          I = c(50, 100, 150)/4)
//'
//' exitprob(b = c(2.963, 2.359, 2.014),
//'          a = c(-0.264, 0.599, 2.014),
//'          theta = c(0.141, 0.204, 0.289),
//'          I = c(81, 121, 160))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List exitprob(
    const Rcpp::NumericVector& b,
    const Rcpp::NumericVector& a = NA_REAL,
    const Rcpp::NumericVector& theta = 0,
    const Rcpp::NumericVector& I = NA_REAL) {

  auto bv = Rcpp::as<std::vector<double>>(b);
  auto av = Rcpp::as<std::vector<double>>(a);
  auto thetav = Rcpp::as<std::vector<double>>(theta);
  auto Iv = Rcpp::as<std::vector<double>>(I);

  auto out = exitprobcpp(bv, av, thetav, Iv);
  return Rcpp::wrap(out);
}


std::vector<double> getBoundcpp(
    const size_t k,
    const std::vector<double>& informationRates,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const std::vector<unsigned char>& efficacyStopping) {

  if (k <= 0) throw std::invalid_argument("k must be provided and positive");

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

  if (asf == "none") {
    for (size_t i = 0; i < kMax-1; ++i) criticalValues[i] = 6.0;
    criticalValues[kMax-1] = boost_qnorm(1.0 - alpha);
    return subset(criticalValues, 0, k);
  }

  if (asf == "of" || asf == "p" || asf == "wt") {
    double Delta;
    if (asf == "of") Delta = 0.0;
    else if (asf == "p") Delta = 0.5;
    else Delta = parameterAlphaSpending; // parameterAlphaSpending holds delta for WT

    // for a given multiplier, compute cumulative upper exit probability - alpha
    std::vector<double> u(kMax);
    std::vector<double> l(kMax, -6.0);
    std::vector<double> theta(kMax, 0.0);
    std::vector<double> u0(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      u0[i] = std::pow(infoRates[i], Delta - 0.5);
    }

    auto f = [&](double aval)->double {
      for (size_t i = 0; i < kMax; ++i) {
        u[i] = aval * u0[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      ListCpp probs = exitprobcpp(u, l, theta, infoRates);
      auto v = probs.get<std::vector<double>>("exitProbUpper");
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    double cwt = brent(f, 0.0, 10.0, 1e-6);
    for (size_t i = 0; i < kMax; ++i) {
      criticalValues[i] = cwt * u0[i];
      if (!effStopping[i]) criticalValues[i] = 6.0;
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
    else criticalValues[0] = boost_qnorm(1.0 - cumAlpha);

    // Preallocate reusable buffers used by the root-finding lambda
    std::vector<double> u_vec; u_vec.reserve(kMax);
    std::vector<double> l_vec(kMax, -6.0);
    std::vector<double> theta_vec(kMax, 0.0);

    // subsequent stages
    for (size_t k1 = 1; k1 < kMax; ++k1) {
      // determine cumulative alpha at this stage
      if (asf == "user") cumAlpha = userAlphaSpending[k1];
      else cumAlpha = errorSpentcpp(spendTime[k1], alpha,
                                    asf,parameterAlphaSpending);

      if (!effStopping[k1]) {
        criticalValues[k1] = 6.0;
        continue;
      }

      // Ensure reusable buffers have size k1+1 and capacity >= kMax
      u_vec.resize(k1 + 1);

      // - copy already computed criticalValues[0..k1-1] into u_vec[0..k1-1]
      // the last entry (u_vec[k1]) will be set by the lambda
      std::memcpy(u_vec.data(), criticalValues.data(), k1 * sizeof(double));

      // Define lambda that only sets the last element of u_vec
      auto f = [&](double aval)->double {
        // set the last element to the current candidate critical value
        u_vec[k1] = aval;

        // exitprobcpp expects exact-sized vectors
        ListCpp probs = exitprobcpp(u_vec, l_vec, theta_vec, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
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
        criticalValues[k1] = brent(f_for_brent, -5.0, 6.0, 1e-6);
      }
    }

    return subset(criticalValues, 0, k);
  }

  throw std::invalid_argument("Invalid value for typeAlphaSpending");
}


//' @title Efficacy Boundaries for Group Sequential Design
//' @description Obtains the efficacy stopping boundaries for a group
//' sequential design.
//'
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
//' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
//' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
//' and \code{spendingTime} must be of full length \code{kMax}, and
//' \code{informationRates} and \code{spendingTime} must end with 1.
//'
//' @return A numeric vector of critical values up to the current look.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' getBound(k = 2, informationRates = c(0.5,1),
//'          alpha = 0.025, typeAlphaSpending = "sfOF")
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector getBound(
    const int k = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  auto result = getBoundcpp(
    static_cast<size_t>(k), infoRates, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha, spendTime, effStopping
  );

  return Rcpp::wrap(result);
}


BoundCacheAlpha::BoundCacheAlpha(
  size_t k,
  const std::vector<double>& infoRates,
  const std::string& asf,
  double asfpar,
  const std::vector<double>& userAlphaSpending,
  const std::vector<double>& spendTime,
  const std::vector<unsigned char>& effStopping,
  size_t maxEntries,
  int alphaPrecision) :
  k_(k),
  infoRates_(infoRates),
  asf_(asf),
  asfpar_(asfpar),
  userAlphaSpending_(userAlphaSpending),
  spendTime_(spendTime),
  effStopping_(effStopping),
  maxEntries_(maxEntries),
  alphaPrecision_(alphaPrecision) {
}

int64_t BoundCacheAlpha::discretize(double alpha) const {
  double scale = std::pow(10.0, alphaPrecision_);
  return static_cast<int64_t>(std::llround(alpha * scale));
}

std::vector<double> BoundCacheAlpha::get(double alpha) {
  int64_t key = discretize(alpha);
  { std::lock_guard<std::mutex> lg(mu_);
    auto it = map_.find(key);
    if (it != map_.end()) {
      usage_.splice(usage_.begin(), usage_, it->second.lruIt);
      return it->second.value;
    }
  }

  // Compute without holding the lock
  std::vector<double> computed =
    getBoundcpp(k_, infoRates_, alpha, asf_, asfpar_, userAlphaSpending_,
                spendTime_, effStopping_);

  std::lock_guard<std::mutex> lg(mu_);
  auto it2 = map_.find(key);
  if (it2 != map_.end()) {
    usage_.splice(usage_.begin(), usage_, it2->second.lruIt);
    return it2->second.value;
  }

  if (map_.size() >= maxEntries_) {
    int64_t lastkey = usage_.back();
    usage_.pop_back();
    map_.erase(lastkey);
  }
  usage_.push_front(key);
  CacheEntry e;
  e.value = std::move(computed);
  e.lruIt = usage_.begin();
  map_.emplace(key, std::move(e));
  return map_[key].value;
}



ListCpp getPower(
    const double alpha,
    const size_t kMax,
    const std::vector<double>& b,
    const std::vector<double>& theta,
    const std::vector<double>& I,
    const std::string& bsf,
    const double bsfpar,
    const std::vector<double>& st,
    const std::vector<unsigned char>& futilityStopping,
    const std::vector<double>& w) { // w is the sqrt of variance ratio

  double thetaSqrtI0 = theta[0] * std::sqrt(I[0]);
  std::vector<double> a(kMax, 0.0);
  std::vector<double> u(kMax), l(kMax, 0.0);
  for (size_t i = 0; i < kMax; ++i) u[i] = b[i] * w[i];

  // reusable buffers for prefixes
  std::vector<double> u1; u1.reserve(kMax);
  std::vector<double> l1; l1.reserve(kMax);

  auto f = [&](double x) -> double {
    // reset futility bounds
    std::fill(a.begin(), a.end(), -6.0);
    double eps = 0.0;

    // first stage
    double cb = errorSpentcpp(st[0], x, bsf, bsfpar); // cumulative beta spent

    if (futilityStopping[0]) {
      eps = boost_pnorm(u[0] - thetaSqrtI0) - cb;
      if (eps < 0.0) return -1.0; // to decrease beta
      a[0] = (boost_qnorm(cb) + thetaSqrtI0) / w[0];
    }

    // subsequent stages
    for (size_t k = 1; k < kMax; ++k) {
      l[k-1] = a[k-1] * w[k-1];
      if (futilityStopping[k]) {
        cb = errorSpentcpp(st[k], x, bsf, bsfpar);

        u1.resize(k + 1);
        l1.resize(k + 1);
        std::memcpy(u1.data(), u.data(), k * sizeof(double));
        std::memcpy(l1.data(), l.data(), k * sizeof(double));
        u1[k] = 6.0;

        // lambda expression for finding futility bound at stage k
        // it is an increasing function in aval, and we want to find
        // the root where it crosses 0
        auto g = [&](double aval) -> double {
          l1[k] = aval * w[k];
          ListCpp probs = exitprobcpp(u1, l1, theta, I);
          auto v = probs.get<std::vector<double>>("exitProbLower");
          double cpl = std::accumulate(v.begin(), v.end(), 0.0);
          return cpl - cb;
        };

        double bk = b[k];
        eps = g(bk);
        double g_minus6 = g(-6.0);
        if (g_minus6 > 0.0) { // no beta spent at current visit
          a[k] = -6.0;
        } else if (eps > 0.0) {
          auto g_for_brent = [&](double x)->double {
            if (x == -6.0) return g_minus6;  // avoid recomputation at 6.0
            if (x == bk) return eps;         // avoid recomputation at b[k]
            return g(x);
          };
          a[k] = brent(g_for_brent, -6.0, bk, 1e-6);
        } else if (k < kMax-1) {
          return -1.0; // to decrease beta
        } // else it is the final look, a[k] = b[k], so we need eps = g(bk) = 0
          // since eps < 0, it can be used to decrease beta more accurately
      }
    }

    return eps;
  };

  double v1 = f(0.0001), v2 = f(1.0 - alpha);
  double beta = 0.0;
  ListCpp probs;
  if (v1 == -1.0 || (v1 < 0.0 && a[kMax-1] == 0.0)) {
    throw std::invalid_argument("Power must be less than 0.9999");
  } else if (v2 > 0.0) {
    throw std::invalid_argument("Power must be greater than alpha");
  } else {
    beta = brent(f, 0.0001, 1.0 - alpha, 1e-6);
    a[kMax-1] = b[kMax-1];
    l[kMax-1] = a[kMax-1] * w[kMax-1];
    probs = exitprobcpp(u, l, theta, I);
  }

  ListCpp result;
  result.push_back(1.0 - beta, "power");
  result.push_back(std::move(a), "futilityBounds");
  result.push_back(std::move(probs), "probs");
  return result;
}


ListCpp getDesigncpp(const double beta,
                     const double IMax,
                     const double theta,
                     const size_t kMax,
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
  if (std::isnan(theta)) {
    throw std::invalid_argument("theta must be provided");
  }
  if (kMax < 1) {
    throw std::invalid_argument("kMax must be a positive integer");
  }

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

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping[kMax-1] != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds);

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

  if (!missingFutilityBounds) {
    if (!(futilityBounds.size() == kMax - 1 || futilityBounds.size() == kMax)) {
      throw std::invalid_argument("Invalid length for futilityBounds");
    }
  }
  if (!missingCriticalValues && !missingFutilityBounds) {
    for (size_t i = 0; i < kMax - 1; ++i) {
      if (futilityBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      throw std::invalid_argument(
          "futilityBounds must meet criticalValues at the final look");
    }
  }

  std::string bsf = typeBetaSpending;
  for (char &c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (unknown == "IMax") {
    if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
        bsf == "sfkd" || bsf == "sfhsd" || bsf == "user" || bsf == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpending");
    }
  } else {
    if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
        bsf == "sfkd" || bsf == "sfhsd" || bsf == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpending");
    }
  }

  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument ("parameterBetaSpending must be positive for sfKD");
  }

  if (unknown == "IMax" && bsf == "user") {
    if (!none_na(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be specified");
    if (userBetaSpending.size() != kMax)
      throw std::invalid_argument("Invalid length of userBetaSpending");
    if (userBetaSpending[0] < 0.0)
      throw std::invalid_argument("userBetaSpending must be nonnegative");
    if (any_nonincreasing(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be nondecreasing");
    if (userBetaSpending[kMax-1] != beta)
      throw std::invalid_argument("userBetaSpending must end with specified beta");
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
  std::vector<double> l(kMax, -6.0), zero(kMax, 0.0);
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
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent[kMax-1];

  // set up futility bounds
  std::vector<double> futBounds = futilityBounds;
  if (kMax > 1) {
    if (missingFutilityBounds && bsf == "none") {
      futBounds = std::vector<double>(kMax, -6.0);
      futBounds[kMax-1] = critValues[kMax-1];
    } else if (!missingFutilityBounds && futBounds.size() == kMax-1) {
      futBounds.push_back(critValues[kMax-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // multiplier for the boundaries under the alternative hypothesis
  std::vector<double> w(kMax, std::sqrt(varianceRatio));
  std::vector<double> u(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    u[i] = critValues[i] * w[i];
  }

  double beta1 = beta;
  double IMax1 = IMax;
  double drift;
  if (unknown == "IMax") {
    std::vector<double> u1; u1.reserve(kMax);
    std::vector<double> l1; l1.reserve(kMax);
    double sqrtt0 = std::sqrt(infoRates[0]);

    auto f = [&](double aval)->double {
      std::vector<double> delta = std::vector<double>(kMax, aval);

      // compute stagewise exit probabilities
      if (!missingFutilityBounds || bsf == "none" || kMax == 1) {
        for (size_t i = 0; i < kMax; ++i) {
          l[i] = futBounds[i] * w[i];
        }
        ListCpp probs = exitprobcpp(u, l, delta, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
        return overallReject - (1.0 - beta);
      } else {
        // initialize futility bound to be updated
        futBounds = std::vector<double>(kMax);
        double eps = 0.0;

        // first stage
        double cb = (bsf == "user") ? userBetaSpending[0] :
          errorSpentcpp(spendTime[0], beta, bsf, parameterBetaSpending);

        if (!futStopping[0]) {
          futBounds[0] = -6.0;
        } else {
          double dt0 = delta[0] * sqrtt0;
          eps = boost_pnorm(u[0] - dt0) - cb;
          if (eps < 0.0) return -1.0; // to decrease drift
          futBounds[0] = (boost_qnorm(cb) + dt0) / w[0];
        }

        // subsequent stages
        for (size_t k = 1; k < kMax; ++k) {
          l[k-1] = futBounds[k-1] * w[k-1];
          cb = (bsf == "user") ? userBetaSpending[k] :
            errorSpentcpp(spendTime[k], beta, bsf, parameterBetaSpending);

          if (!futStopping[k]) {
            futBounds[k] = -6.0;
          } else {
            u1.resize(k + 1);
            l1.resize(k + 1);

            std::memcpy(u1.data(), u.data(), k * sizeof(double));
            u1[k] = 6.0;
            std::memcpy(l1.data(), l.data(), k * sizeof(double));

            // lambda expression for finding futility bound at stage k
            auto g = [&](double aval)->double {
              l1[k] = aval * w[k];
              ListCpp probs = exitprobcpp(u1, l1, delta, infoRates);
              auto v = probs.get<std::vector<double>>("exitProbLower");
              double cpl = std::accumulate(v.begin(), v.end(), 0.0);
              return cpl - cb;
            };

            double bk = critValues[k];
            eps = g(bk);
            double g_minus6 = g(-6.0);

            if (g_minus6 > 0.0) { // no beta spent at current visit
              futBounds[k] = -6.0;
            } else if (eps > 0.0) {
              auto g_for_brent = [&](double x)->double {
                if (x == -6.0) return g_minus6;  // avoid recomputation at 6.0
                if (x == bk) return eps;         // avoid recomputation at b[k]
                return g(x);
              };

              futBounds[k] = brent(g_for_brent, -6.0, bk, 1e-6);
            } else if (k < kMax-1) {
              return -1.0;
            }
          }
        }

        return eps;
      }
    };

    drift = brent(f, 0.0, 6.0, 1e-6);
    IMax1 = sq(drift / theta);
    futBounds[kMax-1] = critValues[kMax-1];
    l[kMax-1] = futBounds[kMax-1] * w[kMax-1];
    std::vector<double> delta = std::vector<double>(kMax, drift);
    probs = exitprobcpp(u, l, delta, infoRates);
  } else {
    drift = theta * std::sqrt(IMax1);
    std::vector<double> delta = std::vector<double>(kMax, drift);
    if (!missingFutilityBounds || bsf == "none" || kMax == 1) {
      for (size_t i = 0; i < kMax; ++i) {
        l[i] = futBounds[i] * w[i];
      }
      probs = exitprobcpp(u, l, delta, infoRates);
      auto v = probs.get<std::vector<double>>("exitProbUpper");
      double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
      beta1 = 1.0 - overallReject;
    } else {
      ListCpp out = getPower(alpha1, kMax, critValues, delta, infoRates, bsf,
                             parameterBetaSpending, spendTime, futStopping, w);
      double overallReject = out.get<double>("power");
      beta1 = 1.0 - overallReject;
      futBounds = out.get<std::vector<double>>("futilityBounds");
      for (size_t i = 0; i < kMax; ++i) {
        l[i] = futBounds[i] * w[i];
      }
      probs = out.get_list("probs");
    }
  }

  double driftf = boost_qnorm(1.0 - alpha1) * w[0] + boost_qnorm(1.0 - beta1);
  double inflationFactor = sq(drift / driftf);

  // output the results
  std::vector<double> information(kMax);
  std::vector<double> efficacyTheta(kMax);
  std::vector<double> futilityTheta(kMax);
  std::vector<double> efficacyP(kMax);
  std::vector<double> futilityP(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    information[i] = IMax1 * infoRates[i];
    efficacyTheta[i] = u[i] / std::sqrt(information[i]);
    futilityTheta[i] = l[i] / std::sqrt(information[i]);
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
    futilityP[i] = 1.0 - boost_pnorm(futBounds[i]);
  }

  // stagewise exit probabilities under H1
  auto pu = probs.get<std::vector<double>>("exitProbUpper");
  auto pl = probs.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpu(kMax), cpl(kMax);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());
  double overallReject = cpu[kMax-1];
  std::vector<double> ptotal(kMax);
  for (size_t i = 0; i < kMax; ++i) ptotal[i] = pu[i] + pl[i];
  double expectedInformationH1 = std::inner_product(
    ptotal.begin(), ptotal.end(), information.begin(), 0.0);

  // stagewise exit probabilities under H0 with binding futility
  ListCpp probsH0 = exitprobcpp(critValues, futBounds, zero, infoRates);
  auto puH0 = probsH0.get<std::vector<double>>("exitProbUpper");
  auto plH0 = probsH0.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpuH0(kMax), cplH0(kMax);
  std::partial_sum(puH0.begin(), puH0.end(), cpuH0.begin());
  std::partial_sum(plH0.begin(), plH0.end(), cplH0.begin());
  double overallRejectH0 = cpuH0[kMax-1];
  std::vector<double> ptotalH0(kMax);
  for (size_t i = 0; i < kMax; ++i) ptotalH0[i] = puH0[i] + plH0[i];
  double expectedInformationH0 = std::inner_product(
    ptotalH0.begin(), ptotalH0.end(), information.begin(), 0.0);

  for (size_t i = 0; i < kMax; ++i) {
    if (critValues[i] == 6) effStopping[i] = 0;
    if (futBounds[i] == -6) futStopping[i] = 0;
  }

  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha1, "alpha");
  overallResults.push_back(overallRejectH0, "attainedAlpha");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(theta, "theta");
  overallResults.push_back(IMax1, "information");
  overallResults.push_back(expectedInformationH1, "expectedInformationH1");
  overallResults.push_back(expectedInformationH0, "expectedInformationH0");
  overallResults.push_back(drift, "drift");
  overallResults.push_back(inflationFactor, "inflationFactor");

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
  byStageResults.push_back(std::move(futilityTheta), "futilityTheta");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(futilityP), "futilityP");
  byStageResults.push_back(std::move(information), "information");
  byStageResults.push_back(std::move(effStopping), "efficacyStopping");
  byStageResults.push_back(std::move(futStopping), "futilityStopping");
  byStageResults.push_back(std::move(puH0), "rejectPerStageH0");
  byStageResults.push_back(std::move(plH0), "futilityPerStageH0");
  byStageResults.push_back(std::move(cpuH0), "cumulativeRejectionH0");
  byStageResults.push_back(std::move(cplH0), "cumulativeFutilityH0");

  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(typeBetaSpending, "typeBetaSpending");
  settings.push_back(parameterBetaSpending, "parameterBetaSpending");
  settings.push_back(userBetaSpending, "userBetaSpending");
  settings.push_back(spendingTime, "spendingTime");
  settings.push_back(varianceRatio, "varianceRatio");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");
  return result;
}


//' @title Power and Sample Size for a Generic Group Sequential Design
//' @description Obtains the maximum information and stopping boundaries
//' for a generic group sequential design assuming a constant treatment
//' effect, or obtains the power given the maximum information and
//' stopping boundaries.
//'
//' @param beta The type II error.
//' @param IMax The maximum information. Either \code{beta} or \code{IMax}
//'   should be provided while the other one should be missing.
//' @param theta The parameter value. Null hypothesis is at \code{theta = 0},
//'   and the alternative hypothesis is one-sided for \code{theta > 0}.
//' @inheritParams param_kMax
//' @param informationRates The information rates. Fixed prior to the trial.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param spendingTime A vector of length \code{kMax} for the error spending
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
//'     - \code{attainedAlpha}: The attained significance level, which is
//'       different from the overall significance level in the presence of
//'       futility stopping.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{theta}: The parameter value.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedInformationH1}: The expected information under H1.
//'
//'     - \code{expectedInformationH0}: The expected information under H0.
//'
//'     - \code{drift}: The drift parameter, equal to
//'       \code{theta*sqrt(information)}.
//'
//'     - \code{inflationFactor}: The inflation factor (relative to the
//'       fixed design).
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'
//'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{futilityPerStage}: The probability for futility stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeFutility}: The cumulative probability for futility
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'
//'     - \code{efficacyTheta}: The efficacy boundaries on the parameter
//'       scale.
//'
//'     - \code{futilityTheta}: The futility boundaries on the parameter
//'       scale.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{futilityP}: The futility boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'
//'     - \code{futilityStopping}: Whether to allow futility stopping.
//'
//'     - \code{rejectPerStageH0}: The probability for efficacy stopping
//'       under H0.
//'
//'     - \code{futilityPerStageH0}: The probability for futility stopping
//'       under H0.
//'
//'     - \code{cumulativeRejectionH0}: The cumulative probability for
//'       efficacy stopping under H0.
//'
//'     - \code{cumulativeFutilityH0}: The cumulative probability for
//'       futility stopping under H0.
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
//'     - \code{typeBetaSpending}: The type of beta spending.
//'
//'     - \code{parameterBetaSpending}: The parameter value for beta
//'       spending.
//'
//'     - \code{userBetaSpending}: The user defined beta spending.
//'
//'     - \code{spendingTime}: The error spending time at each analysis.
//'
//'     - \code{varianceRatio}: The ratio of the variance under H0
//'       to the variance under H1.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Christopher Jennison, Bruce W. Turnbull.
//' Group Sequential Methods with Applications to Clinical Trials.
//' Chapman & Hall/CRC: Boca Raton, 2000, ISBN:0849303168
//'
//' @examples
//'
//' # Example 1: obtain the maximum information given power
//' (design1 <- getDesign(
//'   beta = 0.2, theta = -log(0.7),
//'   kMax = 2, informationRates = c(0.5,1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   typeBetaSpending = "sfP"))
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesign(
//'   IMax = 72.5, theta = -log(0.7),
//'   kMax = 3, informationRates = c(0.5, 0.75, 1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   typeBetaSpending = "sfP"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getDesign(
    const double beta = NA_REAL,
    const double IMax = NA_REAL,
    const double theta = NA_REAL,
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& futilityBounds = NA_REAL,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const Rcpp::NumericVector& userBetaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const double varianceRatio = 1) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto userBeta = Rcpp::as<std::vector<double>>(userBetaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  auto cpp_result = getDesigncpp(
    beta, IMax, theta, static_cast<size_t>(kMax), infoRates,
    effStopping, futStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, futBounds,
    typeBetaSpending, parameterBetaSpending, userBeta, spendTime, varianceRatio
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "design";
  return result;
}


ListCpp getDesignEquivcpp(const double beta,
                          const double IMax,
                          const double thetaLower,
                          const double thetaUpper,
                          const double theta,
                          const size_t kMax,
                          const std::vector<double>& informationRates,
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
  if (std::isnan(theta)) {
    throw std::invalid_argument("theta must be provided");
  }
  if (std::isnan(thetaLower)) {
    throw std::invalid_argument("thetaLower must be provided");
  }
  if (std::isnan(thetaUpper)) {
    throw std::invalid_argument("thetaUpper must be provided");
  }
  if (thetaLower >= theta) {
    throw std::invalid_argument("thetaLower must be less than theta");
  }
  if (thetaUpper <= theta) {
    throw std::invalid_argument("thetaUpper must be greater than theta");
  }
  if (kMax < 1) {
    throw std::invalid_argument("kMax must be a positive integer");
  }

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

  // ----------- End of Input Validation ----------- //

  std::vector<double> u(kMax), l(kMax, -6.0), zero(kMax, 0.0);
  std::vector<double> critValues = criticalValues;

  // obtain criticalValues
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
        u[i] = criticalValues[i];
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      std::vector<unsigned char> effStopping(kMax, 1);
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  std::vector<double> li(kMax, -6.0), ui(kMax, 6.0);
  ListCpp probs = exitprobcpp(critValues, li, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());

  std::vector<double> efficacyP(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
  }

  // we center the margins at theta so that exitprobcpp can be called
  // with zero drift and with max upper bound for Z statistic of 6.0
  // and min lower bound of -6.0, which is numerically more stable
  double deltaLower = thetaLower - theta;
  double deltaUpper = thetaUpper - theta;

  // obtain IMax if needed
  double IMax1 = IMax;
  std::vector<double> I(kMax);
  std::vector<double> sqrtI(kMax);
  std::vector<double> b(kMax);
  std::vector<double> a(kMax);

  if (unknown == "IMax") {
    auto f = [&](double aval)->double {
      for (size_t i = 0; i < kMax; ++i) {
        I[i] = infoRates[i] * aval;
        sqrtI[i] = std::sqrt(I[i]);
        l[i] = critValues[i] + deltaLower * sqrtI[i];  // reject H10 if Z0 >= l
        u[i] = -critValues[i] + deltaUpper * sqrtI[i]; // reject H20 if Z0 <= u
        b[i] = std::max(l[i], li[i]); // ensure b >= li
        a[i] = std::min(u[i], ui[i]); // ensure ui >= a
      }

      ListCpp probs1 = exitprobcpp(b, li, zero, I);
      ListCpp probs2 = exitprobcpp(ui, a, zero, I);
      auto v1 = probs1.get<std::vector<double>>("exitProbUpper");
      auto v2 = probs2.get<std::vector<double>>("exitProbLower");
      double p1 = std::accumulate(v1.begin(), v1.end(), 0.0);
      double p2 = std::accumulate(v2.begin(), v2.end(), 0.0);

      bool cross = false;
      for (size_t i = 0; i < kMax; ++i) {
        if (l[i] <= u[i]) { cross = true; break; }
      }

      double power;
      if (cross) {
        power = p1 + p2 - 1.0;
      } else {
        ListCpp probs = exitprobcpp(l, u, zero, I);
        auto v1x = probs.get<std::vector<double>>("exitProbUpper");
        auto v2x = probs.get<std::vector<double>>("exitProbLower");
        double p1x = std::accumulate(v1x.begin(), v1x.end(), 0.0);
        double p2x = std::accumulate(v2x.begin(), v2x.end(), 0.0);
        power = p1 + p2 - p1x - p2x;
      }

      return power - (1.0 - beta);
    };

    double z0 = boost_qnorm(1.0 - alpha);
    double z1 = boost_qnorm(1.0 - beta);
    double IMax0 = sq((z0 + z1) / deltaLower);
    double IMaxLower = 0.5 * IMax0;
    double IMaxUpper = 1.5 * IMax0;
    IMax1 = brent(f, IMaxLower, IMaxUpper, 1e-6);
  }

  // cumulative rejection probability under H1 (at theta)
  for (size_t i = 0; i < kMax; ++i) {
    I[i] = infoRates[i] * IMax1;
    sqrtI[i] = std::sqrt(I[i]);
    l[i] = critValues[i] + deltaLower * sqrtI[i];
    u[i] = -critValues[i] + deltaUpper * sqrtI[i];
    b[i] = std::max(l[i], li[i]);
    a[i] = std::min(u[i], ui[i]);
  }

  std::vector<double> cpl(kMax), cpu(kMax);
  ListCpp probs1 = exitprobcpp(b, li, zero, I);
  ListCpp probs2 = exitprobcpp(ui, a, zero, I);
  auto v1 = probs1.get<std::vector<double>>("exitProbUpper");
  auto v2 = probs2.get<std::vector<double>>("exitProbLower");
  std::partial_sum(v1.begin(), v1.end(), cpl.begin());
  std::partial_sum(v2.begin(), v2.end(), cpu.begin());

  size_t k = kMax; // index for the first crossing look (0-based)
  for (size_t i = 0; i < kMax; ++i) {
    if (l[i] <= u[i]) { k = i; break; }
  }

  std::vector<double> cp(kMax);
  if (k == 0) { // crossing at the first look
    for (size_t i = 0; i < kMax; ++i) {
      cp[i] = cpl[i] + cpu[i] - 1.0;
    }
  } else {
    std::vector<double> cplx(k), cpux(k);
    std::vector l1 = subset(l, 0, k);
    std::vector u1 = subset(u, 0, k);
    std::vector d1 = subset(zero, 0, k);
    std::vector I1 = subset(I, 0, k);
    ListCpp probs = exitprobcpp(l1, u1, d1, I1);
    auto v1x = probs.get<std::vector<double>>("exitProbUpper");
    auto v2x = probs.get<std::vector<double>>("exitProbLower");
    std::partial_sum(v1x.begin(), v1x.end(), cplx.begin());
    std::partial_sum(v2x.begin(), v2x.end(), cpux.begin());
    for (size_t i = 0; i < k; ++i) {
      cp[i] = cpl[i] + cpu[i] - cplx[i] - cpux[i];
    }
    for (size_t i = k; i < kMax; ++i) {
      cp[i] = cpl[i] + cpu[i] - 1.0;
    }
  }

  // incremental exit probabilities under H1 (at theta)
  std::vector<double> rejectPerStage(kMax);
  rejectPerStage[0] = cp[0];
  for (size_t i = 1; i < kMax; ++i) {
    rejectPerStage[i] = cp[i] - cp[i-1];
  }

  std::vector<double> q(kMax);
  std::memcpy(q.data(), rejectPerStage.data(), kMax * sizeof(double));
  if (kMax > 1) q[kMax-1] = 1.0 - cp[kMax-2];

  double overallReject = cp[kMax-1];
  double expectedInformationH1 = std::inner_product(
    q.begin(), q.end(), I.begin(), 0.0);

  std::vector<double> efficacyThetaLower(kMax), efficacyThetaUpper(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    double thetaBound = critValues[i] / std::sqrt(I[i]);
    efficacyThetaLower[i] = thetaBound + thetaLower;
    efficacyThetaUpper[i] = -thetaBound + thetaUpper;
  }

  // cumulative attained alpha under H10 (at thetaLower)
  for (size_t i = 0; i < kMax; ++i) {
    l[i] = critValues[i];
    u[i] = -critValues[i] + (thetaUpper - thetaLower) * sqrtI[i];
    a[i] = std::min(u[i], ui[i]);
  }
  ListCpp probsH10 = exitprobcpp(ui, a, zero, I);
  auto vH10 = probsH10.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpuH10(kMax);
  std::partial_sum(vH10.begin(), vH10.end(), cpuH10.begin());
  std::vector<double> cplH10 = cumAlphaSpent;

  std::vector<double> cpH10(kMax);
  if (k == 0) {
    for (size_t i = 0; i < kMax; ++i) {
      cpH10[i] = cplH10[i] + cpuH10[i] - 1.0;
    }
  } else {
    std::vector<double> cplH10x(k), cpuH10x(k);
    std::vector l1 = subset(l, 0, k);
    std::vector u1 = subset(u, 0, k);
    std::vector d1 = subset(zero, 0, k);
    std::vector I1 = subset(I, 0, k);
    ListCpp probs = exitprobcpp(l1, u1, d1, I1);
    auto v1x = probs.get<std::vector<double>>("exitProbUpper");
    auto v2x = probs.get<std::vector<double>>("exitProbLower");
    std::partial_sum(v1x.begin(), v1x.end(), cplH10x.begin());
    std::partial_sum(v2x.begin(), v2x.end(), cpuH10x.begin());
    for (size_t i = 0; i < k; ++i) {
      cpH10[i] = cplH10[i] + cpuH10[i] - cplH10x[i] - cpuH10x[i];
    }
    for (size_t i = k; i < kMax; ++i) {
      cpH10[i] = cplH10[i] + cpuH10[i] - 1.0;
    }
  }

  // incremental exit probabilities under H10
  std::vector<double> qH10(kMax);
  qH10[0] = cpH10[0];
  for (size_t i = 1; i < kMax - 1; ++i) qH10[i] = cpH10[i] - cpH10[i-1];
  if (kMax > 1) qH10[kMax-1] = 1.0 - cpH10[kMax-2];

  double attainedAlphaH10 = cpH10[kMax-1];
  double expectedInformationH10 = std::inner_product(
    qH10.begin(), qH10.end(), I.begin(), 0.0);

  // cumulative attained alpha under H20 (at thetaUpper)
  for (size_t i = 0; i < kMax; ++i) {
    l[i] = critValues[i] + (thetaLower - thetaUpper) * sqrtI[i];
    u[i] = -critValues[i];
    b[i] = std::max(l[i], li[i]);
  }

  ListCpp probsH20 = exitprobcpp(b, li, zero, I);
  auto vH20 = probsH20.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cplH20(kMax);
  std::partial_sum(vH20.begin(), vH20.end(), cplH20.begin());
  std::vector<double> cpuH20 = cumAlphaSpent;

  std::vector<double> cpH20(kMax);
  if (k == 0) {
    for (size_t i = 0; i < kMax; ++i) {
      cpH20[i] = cplH20[i] + cpuH20[i] - 1.0;
    }
  } else {
    std::vector<double> cplH20x(k), cpuH20x(k);
    std::vector l1 = subset(l, 0, k);
    std::vector u1 = subset(u, 0, k);
    std::vector d1 = subset(zero, 0, k);
    std::vector I1 = subset(I, 0, k);
    ListCpp probs = exitprobcpp(l1, u1, d1, I1);
    auto v1x = probs.get<std::vector<double>>("exitProbUpper");
    auto v2x = probs.get<std::vector<double>>("exitProbLower");
    std::partial_sum(v1x.begin(), v1x.end(), cplH20x.begin());
    std::partial_sum(v2x.begin(), v2x.end(), cpuH20x.begin());
    for (size_t i = 0; i < k; ++i) {
      cpH20[i] = cplH20[i] + cpuH20[i] - cplH20x[i] - cpuH20x[i];
    }
    for (size_t i = k; i < kMax; ++i) {
      cpH20[i] = cplH20[i] + cpuH20[i] - 1.0;
    }
  }

  // incremental exit probabilities under H20
  std::vector<double> qH20(kMax);
  qH20[0] = cpH20[0];
  for (size_t i = 1; i < kMax - 1; ++i) qH20[i] = cpH20[i] - cpH20[i-1];
  if (kMax > 1) qH20[kMax-1] = 1.0 - cpH20[kMax-2];

  double attainedAlphaH20 = cpH20[kMax-1];
  double expectedInformationH20 = std::inner_product(
    qH20.begin(), qH20.end(), I.begin(), 0.0);

  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha, "alpha");
  overallResults.push_back(attainedAlphaH10, "attainedAlphaH10");
  overallResults.push_back(attainedAlphaH20, "attainedAlphaH20");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(thetaLower, "thetaLower");
  overallResults.push_back(thetaUpper, "thetaUpper");
  overallResults.push_back(theta, "theta");
  overallResults.push_back(IMax1, "information");
  overallResults.push_back(expectedInformationH1, "expectedInformationH1");
  overallResults.push_back(expectedInformationH10, "expectedInformationH10");
  overallResults.push_back(expectedInformationH20, "expectedInformationH20");

  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(rejectPerStage), "rejectPerStage");
  byStageResults.push_back(std::move(cp), "cumulativeRejection");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(cpH10), "cumulativeAttainedAlphaH10");
  byStageResults.push_back(std::move(cpH20), "cumulativeAttainedAlphaH20");
  byStageResults.push_back(std::move(efficacyThetaLower), "efficacyThetaLower");
  byStageResults.push_back(std::move(efficacyThetaUpper), "efficacyThetaUpper");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(I), "information");

  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(spendingTime, "spendingTime");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");
  return result;
}


//' @title Power and Sample Size for a Generic Group Sequential Equivalence
//' Design
//'
//' @description Obtains the maximum information and stopping boundaries
//' for a generic group sequential equivalence design assuming a constant
//' treatment effect, or obtains the power given the maximum information
//' and stopping boundaries.
//'
//' @param beta The type II error.
//' @param IMax The maximum information. Either \code{beta} or \code{IMax}
//'   should be provided while the other one should be missing.
//' @param thetaLower The parameter value at the lower equivalence limit.
//' @param thetaUpper The parameter value at the upper equivalence limit.
//' @param theta The parameter value under the alternative hypothesis.
//' @inheritParams param_kMax
//' @param informationRates The information rates. Fixed prior to the trial.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_criticalValues
//' @param alpha The significance level for each of the two one-sided
//'   tests, e.g., 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//'
//' @details
//' Consider the equivalence design with two one-sided hypotheses:
//' \deqn{H_{10}: \theta \leq \theta_{10},}
//' \deqn{H_{20}: \theta \geq \theta_{20}.}
//' We reject \eqn{H_{10}} at or before look \eqn{k} if
//' \deqn{Z_{1j} = (\hat{\theta}_j - \theta_{10})\sqrt{I_j}
//' \geq b_j}
//' for some \eqn{j=1,\ldots,k}, where \eqn{\{b_j:j=1,\ldots,K\}} are the
//' critical values associated with the specified alpha-spending function,
//' and \eqn{I_j} is the information for \eqn{\theta} (inverse variance of
//' \eqn{\hat{\theta}}) at the
//' \eqn{j}th look. For example,
//' for estimating the risk difference \eqn{\theta = \pi_1 - \pi_2},
//' \deqn{I_j = \left\{\frac{\pi_1 (1-\pi_1)}{n_{1j}} +
//' \frac{\pi_2(1-\pi_2)}{n_{2j}}\right\}^{-1}.}
//' It follows that
//' \deqn{(Z_{1j} \geq b_j) = (Z_j \geq b_j +
//' \theta_{10}\sqrt{I_j}),}
//' where \eqn{Z_j = \hat{\theta}_j \sqrt{I_j}}.
//'
//' Similarly, we reject \eqn{H_{20}} at or before look \eqn{k} if
//' \deqn{Z_{2j} = (\hat{\theta}_j - \theta_{20})\sqrt{I_j}
//' \leq -b_j} for some \eqn{j=1,\ldots,k}. We have
//' \deqn{(Z_{2j} \leq -b_j) = (Z_j \leq - b_j +
//' \theta_{20}\sqrt{I_j}).}
//'
//' Let \eqn{l_j = b_j + \theta_{10}\sqrt{I_j}},
//' and \eqn{u_j = -b_j + \theta_{20}\sqrt{I_j}}.
//' The cumulative probability to reject \eqn{H_0 = H_{10} \cup H_{20}} at
//' or before look \eqn{k} under the alternative hypothesis \eqn{H_1} is
//' given by
//' \deqn{P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j) \cap
//' \cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right) = p_1 + p_2 - p_{12},}
//' where
//' \deqn{p_1 = P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \geq l_j)\right),}
//' \deqn{p_2 = P_\theta\left(\cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \leq u_j)\right),}
//' and
//' \deqn{p_{12} = P_\theta\left(\cup_{j=1}^{k} (Z_j \geq l_j) \cup
//' (Z_j \leq u_j)\right).}
//' Of note, both \eqn{p_1} and \eqn{p_2} can be evaluated using
//' one-sided exit probabilities for group sequential designs.
//' If there exists \eqn{j\leq k} such that \eqn{l_j \leq u_j}, then
//' \eqn{p_{12} = 1}. Otherwise, \eqn{p_{12}} can be evaluated using
//' two-sided exit probabilities for group sequential designs.
//'
//' Since the equivalent hypothesis is tested using two one-sided tests,
//' the type I error is controlled. To evaluate the attained type I error
//' of the equivalence trial under \eqn{H_{10}} (or \eqn{H_{20}}),
//' we simply fix the control group parameters, update the active
//' treatment group parameters according to the null hypothesis, and
//' use the parameters in the power calculation outlined above.
//'
//' @return An S3 class \code{designEquiv} object with three components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{attainedAlphaH10}: The attained significance level under H10.
//'
//'     - \code{attainedAlphaH20}: The attained significance level under H20.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{thetaLower}: The parameter value at the lower equivalence
//'       limit.
//'
//'     - \code{thetaUpper}: The parameter value at the upper equivalence
//'       limit.
//'
//'     - \code{theta}: The parameter value under the alternative hypothesis.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedInformationH1}: The expected information under H1.
//'
//'     - \code{expectedInformationH10}: The expected information under H10.
//'
//'     - \code{expectedInformationH20}: The expected information under H20.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
//'       each of the two one-sided tests.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
//'       the two one-sided tests.
//'
//'     - \code{cumulativeAttainedAlphaH10}: The cumulative probability for
//'       efficacy stopping under H10.
//'
//'     - \code{cumulativeAttainedAlphaH20}: The cumulative probability for
//'       efficacy stopping under H20.
//'
//'     - \code{efficacyThetaLower}: The efficacy boundaries on the
//'       parameter scale for the one-sided null hypothesis at the
//'       lower equivalence limit.
//'
//'     - \code{efficacyThetaUpper}: The efficacy boundaries on the
//'       parameter scale for the one-sided null hypothesis at the
//'       upper equivalence limit.
//'
//'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
//'       each of the two one-sided tests.
//'
//'     - \code{information}: The cumulative information.
//'
//' * \code{settings}: A list containing the following components:
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
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' # Example 1: obtain the maximum information given power
//' (design1 <- getDesignEquiv(
//'   beta = 0.2, thetaLower = log(0.8), thetaUpper = log(1.25),
//'   kMax = 2, informationRates = c(0.5, 1),
//'   alpha = 0.05, typeAlphaSpending = "sfOF"))
//'
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesignEquiv(
//'   IMax = 72.5, thetaLower = log(0.7), thetaUpper = -log(0.7),
//'   kMax = 3, informationRates = c(0.5, 0.75, 1),
//'   alpha = 0.05, typeAlphaSpending = "sfOF"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getDesignEquiv(
    const double beta = NA_REAL,
    const double IMax = NA_REAL,
    const double thetaLower = NA_REAL,
    const double thetaUpper = NA_REAL,
    const double theta = 0,
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.05,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  auto cpp_result = getDesignEquivcpp(
    beta, IMax, thetaLower, thetaUpper, theta,
    static_cast<size_t>(kMax), infoRates, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, spendTime
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "designEquiv";
  return result;
}


ListCpp adaptDesigncpp(double betaNew,
                       double INew,
                       const size_t L,
                       const double zL,
                       const double theta,
                       const double IMax,
                       const size_t kMax,
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
                       const size_t kNew,
                       const std::vector<double>& informationRatesNew,
                       const std::vector<unsigned char>& efficacyStoppingNew,
                       const std::vector<unsigned char>& futilityStoppingNew,
                       const std::string& typeAlphaSpendingNew,
                       const double parameterAlphaSpendingNew,
                       const std::string& typeBetaSpendingNew,
                       const double parameterBetaSpendingNew,
                       const std::vector<double>& userBetaSpendingNew,
                       const std::vector<double>& spendingTimeNew,
                       const double varianceRatio) {

  // ----------- Input Validation ----------- //
  if (std::isnan(betaNew) && std::isnan(INew)) {
    throw std::invalid_argument("betaNew and INew cannot be missing simultaneously");
  }
  if (!std::isnan(betaNew) && !std::isnan(INew)) {
    throw std::invalid_argument("Only one of betaNew and INew should be provided");
  }
  if (!std::isnan(INew) && INew <= 0.0) {
    throw std::invalid_argument("INew must be positive");
  }
  if (std::isnan(theta)) throw std::invalid_argument("theta must be provided");
  if (L < 1) throw std::invalid_argument("L must be a positive integer");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  // Alpha and Beta must be within valid ranges
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(betaNew) && (betaNew < 0.0001 || betaNew >= 1)) {
    throw std::invalid_argument("betaNew must lie in [0.0001, 1)");
  }

  if (std::isnan(IMax)) {
    throw std::invalid_argument("IMax must be provided");
  }
  if (IMax <= 0.0) {
    throw std::invalid_argument("IMax must be positive");
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

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping[kMax-1] != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds);

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

  if (!missingFutilityBounds) {
    if (!(futilityBounds.size() == kMax - 1 || futilityBounds.size() == kMax)) {
      throw std::invalid_argument("Invalid length for futilityBounds");
    }
  }
  if (!missingCriticalValues && !missingFutilityBounds) {
    for (size_t i = 0; i < kMax - 1; ++i) {
      if (futilityBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      throw std::invalid_argument(
          "futilityBounds must meet criticalValues at the final look");
    }
  }

  std::string bsf = typeBetaSpending;
  for (char &c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
      bsf == "sfkd" || bsf == "sfhsd" || bsf == "none")) {
    throw std::invalid_argument("Invalid value for typeBetaSpending");
  }

  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument ("parameterBetaSpending must be positive for sfKD");
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
  std::vector<double> infoRatesNew;
  std::vector<unsigned char> effStoppingNew;
  std::vector<unsigned char> futStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  std::string bsfNew = typeBetaSpendingNew;
  for (char &c : bsfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (MullerSchafer) {
    if (kNew < 1) {
      throw std::invalid_argument("kNew must be a positive integer");
    }

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

    // futStoppingNew: default to all 1s if missing
    if (none_na(futilityStoppingNew)) {
      if (futilityStoppingNew.size() != kNew)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew[kNew-1] != 1)
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

    if (std::isnan(INew)) {
      if (!(bsfNew == "sfof" || bsfNew == "sfp" ||
          bsfNew == "sfkd" || bsfNew == "sfhsd" || bsfNew == "user" ||
          bsfNew == "none")) {
        throw std::invalid_argument("Invalid value for typeBetaSpendingNew");
      }
    } else {
      if (!(bsfNew == "sfof" || bsfNew == "sfp" ||
          bsfNew == "sfkd" || bsfNew == "sfhsd" || bsfNew == "none")) {
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
      if (userBetaSpendingNew.size() != kNew)
        throw std::invalid_argument("Invalid length of userBetaSpendingNew");
      if (userBetaSpendingNew[0] < 0.0)
        throw std::invalid_argument("userBetaSpendingNew must be nonnegative");
      if (any_nonincreasing(userBetaSpendingNew))
        throw std::invalid_argument("userBetaSpendingNew must be nondecreasing");
      if (userBetaSpendingNew[kNew] != betaNew)
        throw std::invalid_argument(
            "userBetaSpendingNew must end with specified betaNew");
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

  if (varianceRatio <= 0.0) {
    throw std::invalid_argument("varianceRatio must be positive");
  }
  // ----------- End of Input Validation ----------- //

  // obtain critical values for the primary trial
  std::vector<double> l(kMax, -6.0), zero(kMax, 0.0);
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
      std::vector<double> u(kMax);
      for (size_t i = 0; i < kMax - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  } else {
    for (size_t i = 0; i < kMax; ++i) {
      if (!effStopping[i]) critValues[i] = 6.0;
    }
    ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    alpha1 = std::accumulate(v.begin(), v.end(), 0.0);
  }

  // obtain futility bounds for the primary trial
  std::vector<double> futBounds = futilityBounds;
  if (kMax > 1) {
    if (missingFutilityBounds && bsf == "none") {
      futBounds = std::vector<double>(kMax, -6.0);
      futBounds[kMax-1] = critValues[kMax-1];
    } else if (!missingFutilityBounds && futBounds.size() == kMax-1) {
      futBounds.push_back(critValues[kMax-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // multiplier for the boundaries under the alternative hypothesis
  std::vector<double> w(kMax, std::sqrt(varianceRatio));
  if (!none_na(futBounds)) {
    std::vector<double> delta(kMax, theta);
    std::vector<double> information(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      information[i] = IMax * infoRates[i];
    }

    ListCpp out = getPower(alpha1, kMax, critValues, delta, information, bsf,
                           parameterBetaSpending, spendTime, futStopping, w);
    futBounds = out.get<std::vector<double>>("futilityBounds");
  }

  // compute conditional alpha, conditional power, and predictive power
  size_t k1 = kMax - L;

  std::vector<double> t1(k1), r1(k1), b1(k1), a1(k1, -6.0), zero1(k1, 0.0);
  for (size_t l = 0; l < k1; ++l) {
    t1[l] = (infoRates[l + L] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
    r1[l] = infoRates[L - 1] / infoRates[l + L];
    b1[l] = (critValues[l + L] - std::sqrt(r1[l]) * zL) / std::sqrt(1.0 - r1[l]);
    if (!effStopping[l + L]) b1[l] = 6.0;
  }

  // conditional type I error
  ListCpp probs = exitprobcpp(b1, a1, zero1, t1);
  auto v0 = probs.get<std::vector<double>>("exitProbUpper");
  double alphaNew = std::accumulate(v0.begin(), v0.end(), 0.0);

  // conditional power and predictive power
  for (size_t l = 0; l < k1; ++l) {
    a1[l] = (futBounds[l + L] - std::sqrt(r1[l]) * zL) / std::sqrt(1.0 - r1[l]);
    if (!futStopping[l + L]) a1[l] = -6.0;
  }

  std::vector<double> I1(k1);
  for (size_t l = 0; l < k1; ++l) {
    I1[l] = IMax * (infoRates[l + L] - infoRates[L - 1]);
  }

  std::vector<double> theta1(k1, theta);
  ListCpp probs1 = exitprobcpp(b1, a1, theta1, I1);
  auto v1 = probs1.get<std::vector<double>>("exitProbUpper");
  double conditionalPower = std::accumulate(v1.begin(), v1.end(), 0.0);

  // predictive power
  auto f = [&](double theta)->double {
    std::vector<double> theta1(k1, theta);
    ListCpp probs = exitprobcpp(b1, a1, theta1, I1);
    auto v = probs.get<std::vector<double>>("exitProbUpper");
    return std::accumulate(v.begin(), v.end(), 0.0);
  };

  double sigma = 1.0 / std::sqrt(IMax * infoRates[L - 1]);
  double mu = zL * sigma;
  double lower = mu - 6.0 * sigma, upper = mu + 6.0 * sigma;
  double predictivePower = intnorm(f, mu, sigma, lower, upper);

  // information for the primary trial
  std::vector<double> information1(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    information1[i] = IMax * infoRates[i];
  }

  ListCpp des2;
  if (!MullerSchafer) {
    effStoppingNew = subset(effStopping, L, kMax);
    futStoppingNew = subset(futStopping, L, kMax);
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
  size_t kc = L + k2;

  double IL = information1[L - 1];
  double sqrtIL = std::sqrt(IL);

  std::vector<double> Ic(kc); // cumulative information for the combined design
  for (size_t i = 0; i < L; ++i) {
    Ic[i] = information1[i];
  }
  for (size_t i = L; i < kc; ++i) {
    Ic[i] = IL + I2[i - L];
  }
  double IMaxc = Ic[kc - 1];

  std::vector<double> tc(kc); // information rates for the combined design
  for (size_t i = 0; i < kc; ++i) {
    tc[i] = Ic[i] / IMaxc;
  }

  std::vector<double> bc(kc); // efficacy bounds for the combined design
  for (size_t i = 0; i < L; ++i) {
    bc[i] = critValues[i];
  }
  for (size_t i = L; i < kc; ++i) {
    if (b2[i - L] != 6.0) {
      bc[i] = (zL * sqrtIL + b2[i - L] * std::sqrt(I2[i - L])) / std::sqrt(Ic[i]);
    } else {
      bc[i] = 6.0;
    }
  }
  std::vector<double> ac(kc); // futility bounds for the combined design
  for (size_t i = 0; i < L; ++i) {
    ac[i] = futBounds[i];
  }
  for (size_t i = L; i < kc; ++i) {
    if (a2[i - L] != -6.0) {
      ac[i] = (zL * sqrtIL + a2[i - L] * std::sqrt(I2[i - L])) / std::sqrt(Ic[i]);
    } else {
      ac[i] = -6.0;
    }
  }

  ListCpp des1;
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
  des1.push_back(predictivePower, "predictivePower");
  des1.push_back(MullerSchafer, "MullerSchafer");

  ListCpp des3;
  des3.push_back(L, "L");
  des3.push_back(zL, "zL");
  des3.push_back(theta, "theta");
  des3.push_back(IMaxc, "maxInformation");
  des3.push_back(kc, "kMax");
  des3.push_back(std::move(tc), "informationRates");
  des3.push_back(std::move(bc), "efficacyBounds");
  des3.push_back(std::move(ac), "futilityBounds");
  des3.push_back(std::move(Ic), "information");

  ListCpp result;
  result.push_back(std::move(des1), "primaryTrial");
  result.push_back(std::move(des2), "secondaryTrial");
  result.push_back(std::move(des3), "integratedTrial");
  return result;
}


//' @title Adaptive Design at an Interim Look
//' @description
//' Calculates the conditional power for specified incremental
//' information, given the interim results, parameter value,
//' data-dependent changes in the error spending function, and
//' the number and spacing of interim looks. Conversely,
//' calculates the incremental information required to attain
//' a specified conditional power, given the interim results,
//' parameter value, data-dependent changes in the error
//' spending function, and the number and spacing of interim looks.
//'
//' @param betaNew The type II error for the secondary trial.
//' @param INew The maximum information of the secondary trial. Either
//'   \code{betaNew} or \code{INew} should be provided, while the other
//'   must be missing.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of
//'   the primary trial.
//' @param theta The parameter value.
//' @param IMax The maximum information of the primary trial. Must be provided.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param futilityStopping Indicators of whether futility stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping for the primary trial. If missing, boundaries
//'   will be computed based on the specified alpha spending function.
//' @param alpha The significance level of the primary trial.
//'   Defaults to 0.025.
//' @param typeAlphaSpending The type of alpha spending for the primary
//'   trial. One of the following:
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function,
//'   "user" for user defined spending, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to \eqn{\Delta} for "WT",
//'   \eqn{\rho} for "sfKD", and \eqn{\gamma} for "sfHSD".
//' @param userAlphaSpending The user-defined alpha spending for the
//'   primary trial. Represents the cumulative alpha spent up to each stage.
//' @param futilityBounds The lower boundaries on the z-test statistic scale
//'   for futility stopping for the primary trial. Defaults to
//'   \code{rep(-6, kMax-1)} if left unspecified.
//' @param typeBetaSpending The type of beta spending for the primary trial.
//'   One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early futility stopping.
//'   Defaults to "none".
//' @param parameterBetaSpending The parameter value of beta spending
//'   for the primary trial. Corresponds to \eqn{\rho} for "sfKD",
//'   and \eqn{\gamma} for "sfHSD".
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
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to \eqn{\Delta} for "WT",
//'   \eqn{\rho} for "sfKD", and \eqn{\gamma} for "sfHSD".
//' @param typeBetaSpendingNew The type of beta spending for the secondary
//'   trial. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function,
//'   "user" for user defined spending, and
//'   "none" for no early futility stopping.
//'   Defaults to "none".
//' @param parameterBetaSpendingNew The parameter value of beta spending
//'   for the secondary trial. Corresponds to \eqn{\rho} for "sfKD",
//'   and \eqn{\gamma} for "sfHSD".
//' @param userBetaSpendingNew The user-defined cumulative beta spending.
//'   Represents the cumulative beta spent up to each stage of the
//'   secondary trial.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRatesNew}.
//' @param varianceRatio The ratio of the variance under H0 to the
//'   variance under H1.
//'
//' @return An \code{adaptDesign} object with three list components:
//'
//' * \code{primaryTrial}: A list of selected information for the primary
//'   trial, including \code{L}, \code{zL}, \code{theta},
//'   \code{maxInformation}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
//'   \code{information}, \code{alpha}, \code{conditionalAlpha},
//'   \code{conditionalPower}, \code{predictivePower}, and
//'   and \code{MullerSchafer}.
//'
//' * \code{secondaryTrial}: A \code{design} object for the secondary trial.
//'
//' * \code{integratedTrial}: A list of selected information for the integrated
//'   trial, including \code{L}, \code{zL}, \code{theta}, \code{maxInformation},
//'   \code{kMax}, \code{informationRates}, \code{efficacyBounds},
//'   \code{futilityBounds}, and \code{information}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Lu Chi, H. M. James Hung, and Sue-Jane Wang.
//' Modification of sample size in group sequential clinical trials.
//' Biometrics 1999;55:853-857.
//'
//' Hans-Helge Muller and Helmut Schafer.
//' Adaptive group sequential designs for clinical trials:
//' Combining the advantages of adaptive and of
//' classical group sequential approaches.
//' Biometrics 2001;57:886-891.
//'
//' @seealso \code{\link{getDesign}}
//'
//' @examples
//'
//' # two-arm randomized clinical trial with a normally distributed endpoint
//' # 90% power to detect mean difference of 15 with a standard deviation of 50
//' # Design the Stage I Trial with 3 looks and Lan-DeMets O'Brien-Fleming type
//' # spending function
//' delta <- 15
//' sigma <- 50
//'
//' (des1 <- getDesignMeanDiff(
//'   beta = 0.1, meanDiff = delta, stDev = sigma,
//'   kMax = 3, alpha = 0.025, typeAlphaSpending = "sfOF"
//' ))
//'
//' s1 <- des1$byStageResults$informationRates
//' b1 <- des1$byStageResults$efficacyBounds
//' n <- des1$overallResults$numberOfSubjects
//'
//' # Monitoring the Stage I Trial
//' L <- 1
//' nL <- des1$byStageResults$numberOfSubjects[L]
//' deltahat <- 8
//' sigmahat <- 55
//' sedeltahat <- sigmahat * sqrt( 4 / nL)
//' zL <- deltahat / sedeltahat
//'
//' # Making an Adaptive Change: Stage I to Stage II
//' # revised clinically meaningful difference downward to 10 power the study
//' # retain the standard deviation at the design stage
//' # Muller & Schafer (2001) method to design the secondary trial
//' # with 2 looks and Lan-DeMets Pocock type spending function
//' # re-estimate sample size to reach 90% conditional power
//' deltaNew <- 10
//'
//' (des2 <- adaptDesign(
//'   betaNew = 0.1, L = L, zL = zL, theta = deltaNew,
//'   IMax = n / (4 * sigma^2), kMax = 3, informationRates = s1,
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   MullerSchafer = TRUE, kNew = 2, typeAlphaSpendingNew = "sfP"
//' ))
//'
//' INew <- des2$secondaryTrial$overallResults$information
//' (nNew <- ceiling(INew * 4 * sigma^2))
//' (nTotal <- nL + nNew)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List adaptDesign(
    double betaNew = NA_REAL,
    double INew = NA_REAL,
    const int L = NA_INTEGER,
    const double zL = NA_REAL,
    const double theta = NA_REAL,
    const double IMax = NA_REAL,
    const int kMax = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& futilityBounds = NA_REAL,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int kNew = NA_INTEGER,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const std::string& typeBetaSpendingNew = "none",
    const double parameterBetaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& userBetaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL,
    const double varianceRatio = 1.0) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto futStoppingNew = convertLogicalVector(futilityStoppingNew);
  auto userBetaNew = Rcpp::as<std::vector<double>>(userBetaSpendingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  auto cpp_result = adaptDesigncpp(
    betaNew, INew, static_cast<size_t>(L), zL, theta, IMax,
    static_cast<size_t>(kMax), infoRates, effStopping,
    futStopping, critValues, alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlpha, futBounds, typeBetaSpending, parameterBetaSpending, spendTime,
    MullerSchafer, static_cast<size_t>(kNew), infoRatesNew,
    effStoppingNew, futStoppingNew,
    typeAlphaSpendingNew, parameterAlphaSpendingNew, typeBetaSpendingNew,
    parameterBetaSpendingNew, userBetaNew, spendTimeNew, varianceRatio
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "adaptDesign";
  return result;
}

