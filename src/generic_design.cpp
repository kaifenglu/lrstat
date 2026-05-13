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


double errorSpentcpp(
    const double t,
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


namespace {

struct ExitProbCppResult {
  std::vector<double> exitProbUpper;
  std::vector<double> exitProbLower;
};

ExitProbCppResult exitprobcpp_impl(
    const std::vector<double>& b,
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
    a1.assign(K, -8.0);
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

  return ExitProbCppResult{
      std::move(exitProbUpper),
      std::move(exitProbLower)
  };
}

} // namespace


ListCpp exitprobcpp(
    const std::vector<double>& b,
    const std::vector<double>& a,
    const std::vector<double>& theta,
    const std::vector<double>& I) {

  auto probs = exitprobcpp_impl(b, a, theta, I);
  ListCpp exitProb;
  exitProb.push_back(std::move(probs.exitProbUpper), "exitProbUpper");
  exitProb.push_back(std::move(probs.exitProbLower), "exitProbLower");
  return exitProb;
}


double exitprobcpp_cum_upper(
    const std::vector<double>& b,
    const std::vector<double>& a,
    const std::vector<double>& theta,
    const std::vector<double>& I) {
  auto probs = exitprobcpp_impl(b, a, theta, I);
  return std::accumulate(probs.exitProbUpper.begin(), probs.exitProbUpper.end(), 0.0);
}


double exitprobcpp_survival_upper(
    const std::vector<double>& b,
    const std::vector<double>& a,
    const std::vector<double>& theta,
    const std::vector<double>& I) {
  return 1.0 - exitprobcpp_cum_upper(b, a, theta, I);
}


//' @title Stagewise Exit Probabilities
//' @description Obtains the stagewise exit probabilities for both efficacy
//' and futility stopping.
//'
//' @param b Upper boundaries on the z-test statistic scale.
//' @param a Lower boundaries on the z-test statistic scale. Defaults to
//'   \code{c(rep(-8.0, kMax-1), b[kMax])} if left unspecified, where
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


// Compute the p-value given theta, look L, observed z at look L (zL),
// critical values vector b (length L) and information vector I (length L).
double f_pvalue(const double theta,
                const size_t L,
                const double zL,
                const std::vector<double>& b,
                const std::vector<double>& I) {
  // Build the vectors required by exitprobcpp:
  // upper: first L-1 from b, last = zL
  // lower: all -8.0
  // mu: all = theta
  std::vector<double> upper(L);
  if (L > 1) {
    std::memcpy(upper.data(), b.data(), (L - 1) * sizeof(double));
  }
  upper[L - 1] = zL;

  std::vector<double> lower(L, -8.0);
  std::vector<double> mu(L, theta);

  return exitprobcpp_cum_upper(upper, lower, mu, I);
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
    for (size_t i = 0; i < kMax-1; ++i) criticalValues[i] = 8.0;
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
    std::vector<double> l(kMax, -8.0);
    std::vector<double> theta(kMax, 0.0);
    std::vector<double> u0(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      u0[i] = std::pow(infoRates[i], Delta - 0.5);
    }

    auto f = [&](double aval)->double {
      for (size_t i = 0; i < kMax; ++i) {
        u[i] = aval * u0[i];
        if (!effStopping[i]) u[i] = 8.0;
      }

      ListCpp probs = exitprobcpp(u, l, theta, infoRates);
      auto v = probs.get<std::vector<double>>("exitProbUpper");
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    double cwt = brent(f, 0.0, 10.0, 1e-6);
    for (size_t i = 0; i < kMax; ++i) {
      criticalValues[i] = cwt * u0[i];
      if (!effStopping[i]) criticalValues[i] = 8.0;
    }
    return subset(criticalValues, 0, k);
  }

  if (asf == "sfof" || asf == "sfp" || asf == "sfkd" ||
      asf == "sfhsd" || asf == "user") {
    // stage 1
    double cumAlpha;
    if (asf == "user") cumAlpha = userAlphaSpending[0];
    else cumAlpha = errorSpentcpp(spendTime[0], alpha, asf, parameterAlphaSpending);

    if (!effStopping[0]) criticalValues[0] = 8.0;
    else criticalValues[0] = boost_qnorm(1.0 - cumAlpha);

    // Preallocate reusable buffers used by the root-finding lambda
    std::vector<double> u_vec; u_vec.reserve(kMax);
    std::vector<double> l_vec(kMax, -8.0);
    std::vector<double> theta_vec(kMax, 0.0);

    // subsequent stages
    for (size_t k1 = 1; k1 < kMax; ++k1) {
      // determine cumulative alpha at this stage
      if (asf == "user") cumAlpha = userAlphaSpending[k1];
      else cumAlpha = errorSpentcpp(spendTime[k1], alpha,
                                    asf,parameterAlphaSpending);

      if (!effStopping[k1]) {
        criticalValues[k1] = 8.0;
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


      double f_8 = f(8.0);
      if (f_8 > 0.0) { // no alpha spent at current visit
        criticalValues[k1] = 8.0;
      } else {
        auto f_for_brent = [&](double x)->double {
          if (x == 8.0) return f_8; // avoid recomputation at 8.0
          return f(x);
        };
        criticalValues[k1] = brent(f_for_brent, -5.0, 8.0, 1e-6);
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


// The following function computes the conditional power and futility
// bounds of a secondary trial given the interim result of a primary
// trial and the beta-spending function of the secondary trial, where
// IL, thetaL, and zL are the information and Wald statistic at the
// interim look of the primary trial, alphaNew, k2, bsfNew,
// parameterBetaSpendingNew, spendTimeNew,
// futStoppingNew are parameters of the secondary trial, wc, thetac,
// critValues2, futBounds2 and Ic are boundaries and information
// of integrated trial.
ListCpp getPower(
    const double alphaNew,
    const size_t k2,
    const std::vector<double>& critValues2,
    const std::vector<double>& thetac,
    const std::vector<double>& Ic,
    const std::string& bsfNew,
    const double parameterBetaSpendingNew,
    const std::vector<double>& spendTimeNew,
    const std::vector<unsigned char>& futStoppingNew,
    const std::vector<double>& wc,
    const double IL,
    const double thetaL,
    const double zL) { // w is the sqrt of variance ratio

  std::vector<double> I2(k2), sqrtI2(k2), sqrtIc(k2), wsqrtIc(k2);
  for (size_t k = 0; k < k2; ++k) {
    I2[k] = Ic[k] - IL;
    sqrtI2[k] = std::sqrt(I2[k]);
    sqrtIc[k] = std::sqrt(Ic[k]);
    wsqrtIc[k] = wc[k] * sqrtIc[k];
  }

  double sqrtIL = std::sqrt(IL);
  double zscaled = zL * sqrtIL;

  std::vector<double> a2(k2, -8.0);
  std::vector<double> b2(k2);
  for (size_t k = 0; k < k2; ++k) {
    b2[k] = (critValues2[k] * wsqrtIc[k] - zscaled) / sqrtI2[k];
  }

  std::vector<double> theta2(k2);
  for (size_t k = 0; k < k2; ++k) {
    theta2[k] = (thetac[k] * Ic[k] - thetaL * IL) / I2[k];
  }

  double mu0 = theta2[0] * sqrtI2[0];

  // reusable buffers for prefixes
  std::vector<double> u; u.reserve(k2);
  std::vector<double> l; l.reserve(k2);

  ListCpp probs;
  std::vector<double> v;
  std::vector<double> futBounds2(k2);
  auto f = [&](double x) -> double {
    // reset futility bounds
    std::fill(futBounds2.begin(), futBounds2.end(), -8.0);
    std::fill(a2.begin(), a2.end(), -8.0);
    double eps = 0.0, cb = 0.0;

    // first stage
    if (futStoppingNew[0]) {
      cb = errorSpentcpp(spendTimeNew[0], x, bsfNew, parameterBetaSpendingNew);
      eps = boost_pnorm(b2[0] - mu0) - cb;
      if (eps < 0.0) return -1.0; // to decrease beta
      a2[0] = boost_qnorm(cb) + mu0;
      futBounds2[0] = (a2[0] * sqrtI2[0] + zscaled) / wsqrtIc[0];
    }

    // subsequent stages
    for (size_t k = 1; k < k2; ++k) {
      if (futStoppingNew[k]) {
        cb = errorSpentcpp(spendTimeNew[k], x, bsfNew, parameterBetaSpendingNew);

        a2[k-1] = (futBounds2[k-1] * wsqrtIc[k-1] - zscaled) / sqrtI2[k-1];

        u.resize(k + 1);
        l.resize(k + 1);
        std::memcpy(u.data(), b2.data(), k * sizeof(double));
        std::memcpy(l.data(), a2.data(), k * sizeof(double));

        // lambda expression for finding futility bound at stage k
        // it is an increasing function in aval, and we want to find
        // the root where it crosses 0
        auto g = [&](double aval) -> double {
          a2[k] = (aval * wsqrtIc[k] - zscaled) / sqrtI2[k];
          l[k] = a2[k];
          u[k] = l[k];
          probs = exitprobcpp(u, l, theta2, I2);
          v = probs.get<std::vector<double>>("exitProbLower");
          double cpl = std::accumulate(v.begin(), v.end(), 0.0);
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
          return -1.0; // to decrease beta
        } // else it is the final look, a[k] = b[k], so we need eps = g(bk) = 0
          // since eps < 0, it can be used to decrease beta more accurately
      }
    }

    return eps;
  };

  double v1 = f(0.0001), v2 = f(1.0 - alphaNew);
  double beta = 0.0;
  if (v1 == -1.0 || (v1 < 0.0 && futBounds2[k2-1] == -8.0)) {
    throw std::invalid_argument("Power must be less than 0.9999");
  } else if (v2 > 0.0) {
    throw std::invalid_argument("Power must be greater than alpha");
  } else {
    auto f_for_brent = [&](double x)->double {
      if (x == 0.0001) return v1;  // avoid recomputation at 0.0001
      if (x == 1.0 - alphaNew) return v2;  // avoid recomputation at 1.0 - alpha
      return f(x);
    };

    beta = brent(f_for_brent, 0.0001, 1.0 - alphaNew, 1e-6);
    futBounds2[k2-1] = critValues2[k2-1];
    a2[k2-1] = (futBounds2[k2-1] * wsqrtIc[k2-1] - zscaled) / sqrtI2[k2-1];
    probs = exitprobcpp(b2, a2, theta2, I2);
  }

  ListCpp result;
  result.push_back(1.0 - beta, "power");
  result.push_back(std::move(futBounds2), "futilityBounds");
  result.push_back(std::move(probs), "probs");
  return result;
}


ListCpp getDesigncpp(
    const double beta,
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
    const std::vector<double>& futilityCP,
    const std::vector<double>& futilityTheta,
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

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping.back() != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds)
    && !none_na(futilityCP) && !none_na(futilityTheta);

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
    if (userBetaSpending.back() != beta)
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
    if (spendingTime.back() != 1.0)
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
  std::vector<double> l(kMax, -8.0), zero(kMax, 0.0);
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
        if (!effStopping[i]) u[i] = 8.0;
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, -5.0, 8.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  } else {
    for (size_t i = 0; i < kMax; ++i) {
      if (!effStopping[i]) critValues[i] = 8.0;
    }
  }

  ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent[kMax-1];

  std::vector<double> w(kMax, std::sqrt(varianceRatio));

  // set up futility bounds
  std::vector<double> futBounds(kMax, NaN);
  if (kMax > 1) {
    if (missingFutilityBounds && bsf == "none") {
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
        double c2 = critValues[kMax - 1] * w[kMax - 1];
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i]) *
            (c2 - std::sqrt(1 - infoRates[i]) * boost_qnorm(1 - futilityCP[i]));
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityCP values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax-1] = critValues[kMax-1];
      } else if (!std::isnan(IMax)) {
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i] * IMax) * futilityTheta[i] / w[i];
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

  // multiplier for the boundaries under the alternative hypothesis
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

    auto f = [&](double x)->double {
      std::vector<double> delta = std::vector<double>(kMax, x);

      // compute stagewise exit probabilities
      if (!missingFutilityBounds || bsf == "none" || kMax == 1) {
        if (!none_na(futilityBounds) && !none_na(futilityCP) &&
            none_na(futilityTheta)) {
          double maxInformation = sq(x / theta);
          for (size_t i = 0; i < kMax - 1; ++i) {
            futBounds[i] = std::sqrt(infoRates[i] * maxInformation) *
              futilityTheta[i] / w[i];
            if (futBounds[i] > critValues[i]) return -1.0; // to decrease drift
            l[i] = futBounds[i] * w[i];
          }
          futBounds[kMax-1] = critValues[kMax-1];
          l[kMax-1] = futBounds[kMax-1] * w[kMax-1];
        } else {
          for (size_t i = 0; i < kMax; ++i) {
            l[i] = futBounds[i] * w[i];
          }
        }

        probs = exitprobcpp(u, l, delta, infoRates);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
        return (1.0 - overallReject) - beta;
      } else {
        // initialize futility bound to be updated
        std::fill(futBounds.begin(), futBounds.end(), -8.0);
        std::fill(l.begin(), l.end(), -8.0);
        double eps = 0.0, cb = 0.0;

        // first stage
        if (futStopping[0]) {
          cb = (bsf == "user") ? userBetaSpending[0] :
          errorSpentcpp(spendTime[0], beta, bsf, parameterBetaSpending);

          double dt0 = delta[0] * sqrtt0;
          eps = boost_pnorm(u[0] - dt0) - cb;
          if (eps < 0.0) return -1.0; // to decrease drift
          futBounds[0] = (boost_qnorm(cb) + dt0) / w[0];
        }

        // subsequent stages
        for (size_t k = 1; k < kMax; ++k) {
          l[k-1] = futBounds[k-1] * w[k-1];

          if (futStopping[k]) {
            u1.resize(k + 1);
            l1.resize(k + 1);
            std::memcpy(u1.data(), u.data(), k * sizeof(double));
            std::memcpy(l1.data(), l.data(), k * sizeof(double));

            cb = (bsf == "user") ? userBetaSpending[k] :
              errorSpentcpp(spendTime[k], beta, bsf, parameterBetaSpending);

            // lambda expression for finding futility bound at stage k
            auto g = [&](double aval)->double {
              l1[k] = aval * w[k];
              u1[k] = l1[k];
              probs = exitprobcpp(u1, l1, delta, infoRates);
              v = probs.get<std::vector<double>>("exitProbLower");
              double cpl = std::accumulate(v.begin(), v.end(), 0.0);
              return cpl - cb;
            };

            double bk = critValues[k];
            eps = g(bk);
            double g_minus8 = g(-8.0);

            if (g_minus8 > 0.0) { // no beta spent at current visit
              futBounds[k] = -8.0;
            } else if (eps > 0.0) {
              auto g_for_brent = [&](double aval)->double {
                if (aval == -8.0) return g_minus8;  // avoid recomputation at 8.0
                if (aval == bk) return eps;         // avoid recomputation at b[k]
                return g(aval);
              };

              futBounds[k] = brent(g_for_brent, -8.0, bk, 1e-6);
            } else if (k < kMax-1) {
              return -1.0;
            }
          }
        }

        return eps;
      }
    };

    drift = brent(f, 0.0, 8.0, 1e-6);
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
  std::vector<double> futTheta(kMax);
  std::vector<double> efficacyP(kMax);
  std::vector<double> futilityP(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    information[i] = IMax1 * infoRates[i];
    efficacyTheta[i] = u[i] / std::sqrt(information[i]);
    futTheta[i] = l[i] / std::sqrt(information[i]);
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
    futilityP[i] = 1.0 - boost_pnorm(futBounds[i]);
  }

  // stagewise exit probabilities under H1
  auto pu = probs.get<std::vector<double>>("exitProbUpper");
  auto pl = probs.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpu(kMax), cpl(kMax);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());
  double overallReject = cpu.back();
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
    if (critValues[i] == 8) effStopping[i] = 0;
    if (futBounds[i] == -8) futStopping[i] = 0;
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
  byStageResults.push_back(std::move(futTheta), "futilityTheta");
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
//' @param futilityCP The futility bounds on the conditional power scale.
//' @param futilityTheta The futility bounds on the parameter scale.
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param varianceRatio The ratio of the variance under H0 to the
//'   variance under H1.
//'
//' @details
//' The function determines efficacy and futility bounds based on the inputs
//' provided, following a clear priority order.
//'
//' \strong{Efficacy bounds:}
//' If \code{criticalValues} are supplied, they take precedence and all
//' alpha-spending parameters are ignored. Otherwise, efficacy bounds are
//' derived from the specified alpha-spending function.
//'
//' \strong{Futility bounds:}
//' Futility inputs are evaluated in the following order of priority:
//' \enumerate{
//'   \item If \code{futilityBounds} are provided, they override all other
//'   futility-related inputs (\code{futilityCP}, \code{futilityTheta},
//'   and beta-spending parameters).
//'
//'   \item If \code{futilityBounds} are not provided but \code{futilityCP}
//'   is specified, then \code{futilityTheta} and beta-spending parameters
//'   are ignored.
//'
//'   \item If only \code{futilityTheta} is provided, beta-spending parameters
//'   are ignored.
//'
//'   \item If none of \code{futilityBounds}, \code{futilityCP},
//'   or \code{futilityTheta} are specified, futility bounds are computed
//'   using the beta-spending approach.
//' }
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
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const double varianceRatio = 1) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto userBeta = Rcpp::as<std::vector<double>>(userBetaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

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

  auto cpp_result = getDesigncpp(
    beta, IMax, theta, static_cast<size_t>(kMax), infoRates,
    effStopping, futStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha,
    futBounds, futCP, futTheta,
    typeBetaSpending, parameterBetaSpending, userBeta, spendTime, varianceRatio
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "design";
  return result;
}


ListCpp getDesignEquivcpp(
    const double beta,
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
    if (informationRates.back() != 1.0)
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

  std::vector<double> u(kMax), l(kMax, -8.0), zero(kMax, 0.0);
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

      critValues[kMax-1] = brent(f, -5.0, 8.0, 1e-6);
    } else {
      std::vector<unsigned char> effStopping(kMax, 1);
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  std::vector<double> li(kMax, -8.0), ui(kMax, 8.0);
  ListCpp probs = exitprobcpp(critValues, li, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());

  std::vector<double> efficacyP(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
  }

  // we center the margins at theta so that exitprobcpp can be called
  // with zero drift and with max upper bound for Z statistic of 8.0
  // and min lower bound of -8.0, which is numerically more stable
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
    const Rcpp::Nullable<Rcpp::NumericVector> criticalValues = R_NilValue,
    const double alpha = 0.05,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  std::vector<double> critValues;
  if (criticalValues.isNotNull()) {
    critValues = Rcpp::as<std::vector<double>>(criticalValues);
  } else {
    critValues = std::vector<double>(1, NaN);
  }

  auto cpp_result = getDesignEquivcpp(
    beta, IMax, thetaLower, thetaUpper, theta,
    static_cast<size_t>(kMax), infoRates, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, spendTime
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "designEquiv";
  return result;
}


ListCpp adaptDesigncpp(
    double betaNew,
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
  if (L < 1) throw std::invalid_argument("L must be a positive integer");
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(theta)) throw std::invalid_argument("theta must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  // Alpha and Beta must be within valid ranges
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

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping.back() != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds)
    && !none_na(futilityCP) && !none_na(futilityTheta);

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

  size_t k1 = kMax - L;
  if (!MullerSchafer) {
    infoRatesNew.resize(k1);
    for (size_t i = 0; i < k1; ++i) {
      infoRatesNew[i] =
        (infoRates[L + i] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
    }

    effStoppingNew.resize(k1);
    std::memcpy(effStoppingNew.data(), effStopping.data() + L,
                k1 * sizeof(unsigned char));

    if (none_na(futilityStoppingNew)) {
      if (futilityStoppingNew.size() != k1)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew.back() != 1)
        throw std::invalid_argument("futilityStoppingNew must end with 1");
      futStoppingNew = futilityStoppingNew; // copy
    } else {
      futStoppingNew.resize(k1);
      std::memcpy(futStoppingNew.data(), futStopping.data() + L,
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

  bool missingFutilityBoundsInt = !none_na(futilityBoundsInt)
    && !none_na(futilityCPInt) && !none_na(futilityThetaInt);

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
    if (missingFutilityBoundsInt && !(bsfNew == "sfof" || bsfNew == "sfp" ||
        bsfNew == "sfkd" || bsfNew == "sfhsd" || bsfNew == "user" ||
        bsfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpendingNew");
    }
  } else {
    if (missingFutilityBoundsInt && !(bsfNew == "sfof" || bsfNew == "sfp" ||
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


  if (varianceRatio <= 0.0) {
    throw std::invalid_argument("varianceRatio must be positive");
  }
  // ----------- End of Input Validation ----------- //

  ListCpp probs;
  std::vector<double> v;

  // obtain critical values for the primary trial
  std::vector<double> l(kMax, -8.0), zero(kMax, 0.0);
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
        if (!effStopping[i]) u[i] = 8.0;
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        probs = exitprobcpp(u, l, zero, infoRates);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, -5.0, 8.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  } else {
    for (size_t i = 0; i < kMax; ++i) {
      if (!effStopping[i]) critValues[i] = 8.0;
    }
    probs = exitprobcpp(critValues, l, zero, infoRates);
    v = probs.get<std::vector<double>>("exitProbUpper");
    alpha1 = std::accumulate(v.begin(), v.end(), 0.0);
  }

  std::vector<double> w(kMax, std::sqrt(varianceRatio));

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
        double c2 = critValues[kMax - 1] * w[kMax - 1];
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
          futBounds[i] = std::sqrt(infoRates[i] * IMax) * futilityTheta[i] / w[i];
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

  // compute conditional alpha, conditional power, and predictive power
  std::vector<double> r1(k1), b1(k1), a1(k1, -8.0), zero1(k1, 0.0);
  for (size_t i = 0; i < k1; ++i) {
    r1[i] = infoRates[L - 1] / infoRates[i + L];
    b1[i] = (critValues[i + L] - std::sqrt(r1[i]) * zL) / std::sqrt(1.0 - r1[i]);
    if (!effStoppingNew[i]) b1[i] = 8.0;
  }

  // conditional type I error
  probs = exitprobcpp(b1, a1, zero1, infoRatesNew);
  auto v0 = probs.get<std::vector<double>>("exitProbUpper");
  double alphaNew = std::accumulate(v0.begin(), v0.end(), 0.0);

  // conditional power and predictive power
  for (size_t i = 0; i < k1; ++i) {
    a1[i] = (futBounds[i + L] - std::sqrt(r1[i]) * zL) / std::sqrt(1.0 - r1[i]);
    if (!futStoppingNew[i]) a1[i] = -8.0;
  }

  std::vector<double> I1(k1);
  for (size_t i = 0; i < k1; ++i) {
    I1[i] = information1[i + L] - information1[L - 1];
  }

  std::vector<double> theta1(k1, theta);
  probs = exitprobcpp(b1, a1, theta1, I1);
  auto v1 = probs.get<std::vector<double>>("exitProbUpper");
  double conditionalPower = std::accumulate(v1.begin(), v1.end(), 0.0);

  // predictive power
  auto f = [&](double theta)->double {
    std::vector<double> theta1(k1, theta);
    probs = exitprobcpp(b1, a1, theta1, I1);
    v = probs.get<std::vector<double>>("exitProbUpper");
    return std::accumulate(v.begin(), v.end(), 0.0);
  };

  double sigma = 1.0 / std::sqrt(IMax * infoRates[L - 1]);
  double mu = zL * sigma;
  double lower = mu - 8.0 * sigma, upper = mu + 8.0 * sigma;
  double predictivePower = intnorm(f, mu, sigma, lower, upper);

  double IL = information1[L - 1];
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
  std::vector<double> wc(k2, std::sqrt(varianceRatio));
  std::vector<double> theta2(k2, theta);

  std::vector<double> I2(k2), Ic(k2), sqrtI2(k2), sqrtIc(k2), wsqrtIc(k2);

  if (std::isnan(betaNew)) {
    for (size_t i = 0; i < k2; ++i) {
      I2[i] = INew * infoRatesNew[i];
      Ic[i] = I2[i] + IL;
      sqrtI2[i] = std::sqrt(I2[i]);
      sqrtIc[i] = std::sqrt(Ic[i]);
      wsqrtIc[i] = wc[i] * sqrtIc[i];
    }

    for (size_t i = 0; i < k2; ++i) {
      critValues2[i] = (b2[i] * sqrtI2[i] + zscaled) / wsqrtIc[i];
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
          double c2 = critValues2[k2 - 1] * wc[k2 - 1];
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
      ListCpp out = getPower(
        alphaNew, k2, critValues2, theta2, Ic, bsfNew,
        parameterBetaSpendingNew, spendTimeNew, futStoppingNew,
        wc, IL, theta, zL);
      futBounds2 = out.get<std::vector<double>>("futilityBounds");
    }

    for (size_t i = 0; i < k2; ++i) {
      a2[i] = (futBounds2[i] * wsqrtIc[i] - zscaled) / sqrtI2[i];
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
        wsqrtIc[i] = wc[i] * sqrtIc[i];
      }

      double mu0 = theta * sqrtI2[0];

      for (size_t i = 0; i < k2; ++i) {
        critValues2[i] = (b2[i] * sqrtI2[i] + zscaled) / wsqrtIc[i];
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
            double c2 = critValues2[k2 - 1] * wc[k2 - 1];
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
          a2[i] = (futBounds2[i] * wsqrtIc[i] - zscaled) / sqrtI2[i];
        }
        a2[k2 - 1] = b2[k2 - 1];

        probs = exitprobcpp(b2, a2, theta2, I2);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
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
          futBounds2[0] = (a2[0] * sqrtI2[0] + zscaled) / wsqrtIc[0];
        }

        // subsequent stages
        for (size_t k = 1; k < k2; ++k) {
          if (futStoppingNew[k]) {
            cb = (bsfNew == "user") ? userBetaSpendingNew[k] :
            errorSpentcpp(spendTimeNew[k], betaNew, bsfNew,
                          parameterBetaSpendingNew);

            a2[k-1] = (futBounds2[k-1] * wsqrtIc[k-1] - zscaled) / sqrtI2[k-1];

            u.resize(k + 1);
            l.resize(k + 1);
            std::memcpy(u.data(), b2.data(), k * sizeof(double));
            std::memcpy(l.data(), a2.data(), k * sizeof(double));

            // lambda expression for finding futility bound at stage k
            auto g = [&](double aval)->double {
              a2[k] = (aval * wsqrtIc[k] - zscaled) / sqrtI2[k];
              l[k] = a2[k];
              u[k] = l[k];
              probs = exitprobcpp(u, l, theta2, I2);
              v = probs.get<std::vector<double>>("exitProbLower");
              double cpl = std::accumulate(v.begin(), v.end(), 0.0);
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
  v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cpu1(k2);
  std::partial_sum(v.begin(), v.end(), cpu1.begin());
  double p2 = cpu1.back();

  auto v2 = probs.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpu2(k2);
  std::partial_sum(v2.begin(), v2.end(), cpu2.begin());


  // combined design
  size_t kc = L + k2;
  std::vector<double> Ic_full(kc); // cumulative information for the combined design
  std::copy_n(information1.data(), L, Ic_full.data());
  std::copy_n(Ic.data(), k2, Ic_full.data() + L);

  double IMaxc = Ic_full[kc - 1];
  std::vector<double> infoRates_full(kc);
  for (size_t i = 0; i < kc; ++i) infoRates_full[i] = Ic_full[i] / IMaxc;

  std::vector<double> critValues_full(kc);
  std::copy_n(critValues.data(), L, critValues_full.data());
  std::copy_n(critValues2.data(), k2, critValues_full.data() + L);

  std::vector<double> futBounds_full(kc);
  std::copy_n(futBounds.data(), L, futBounds_full.data());
  std::copy_n(futBounds2.data(), k2, futBounds_full.data() + L);

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
//' @param theta The assumed parameter value.
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
//' @param futilityBounds The lower boundaries on the z-test statistic scale
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
//' * \code{secondaryTrial}: A list of selected information for the secondary
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
//' INew <- des2$maxInformation
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
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL,
    const double varianceRatio = 1.0) {

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

  auto cpp_result = adaptDesigncpp(
    betaNew, INew, static_cast<size_t>(L), zL, theta, IMax,
    static_cast<size_t>(kMax), infoRates, effStopping,
    futStopping, critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, futCP, futTheta, spendTime,
    MullerSchafer, static_cast<size_t>(kNew), infoRatesNew,
    effStoppingNew, futStoppingNew, typeAlphaSpendingNew,
    parameterAlphaSpendingNew, futBoundsInt, futCPInt, futThetaInt,
    typeBetaSpendingNew, parameterBetaSpendingNew, userBetaNew,
    spendTimeNew, varianceRatio
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "adaptDesign";
  return result;
}
