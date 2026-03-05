#include "enrollment_event.h"
#include "utilities.h"
#include "dataframe_list.h"

#include <algorithm>    // max, min
#include <cmath>        // exp
#include <cstddef>      // size_t
#include <limits>       // numeric_limits
#include <stdexcept>    // invalid_argument
#include <utility>      // make_pair, pair
#include <vector>       // vector

#include <Rcpp.h>

using std::size_t;


double accrual1(
    const double time,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const double accrualDuration) {

  // up to end of enrollment
  double t = std::max(std::min(time, accrualDuration), 0.0);

  // identify the time interval containing t
  size_t m = findInterval1(t, accrualTime);

  // sum up patients enrolled in each interval up to t
  double n = 0;
  for (size_t j = 0; j < m; ++j) {
    if (j < m - 1) {
      n += accrualIntensity[j] * (accrualTime[j + 1] - accrualTime[j]);
    } else {
      n += accrualIntensity[j] * (t - accrualTime[j]);
    }
  }

  return n;
}


//' @title Number of Enrolled Subjects
//' @description Obtains the number of subjects enrolled by given calendar
//' times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of enrolled subjects.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_accrualDuration
//'
//' @return A vector of total number of subjects enrolled by the
//' specified calendar times.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Example 1: Uniform enrollment with 20 patients per month for 12 months.
//'
//' accrual(time = 3, accrualTime = 0, accrualIntensity = 20,
//'         accrualDuration = 12)
//'
//'
//' # Example 2: Piecewise accrual, 10 patients per month for the first
//' # 3 months, and 20 patients per month thereafter. Patient recruitment
//' # ends at 12 months for the study.
//'
//' accrual(time = c(2, 9), accrualTime = c(0, 3),
//'         accrualIntensity = c(10, 20), accrualDuration = 12)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> accrual(
    const std::vector<double>& time,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const double accrualDuration) {

  size_t k = time.size();
  std::vector<double> n(k);
  for (size_t i = 0; i < k; ++i) {
    n[i] = accrual1(time[i], accrualTime, accrualIntensity, accrualDuration);
  }

  return n;
}


double getAccrualDurationFromN1(
    const double nsubjects,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity) {

  size_t J = accrualTime.size();
  std::vector<double> p(J);
  p[0] = 0;
  for (size_t j = 0; j < J - 1; ++j) {
    p[j + 1] = p[j] + accrualIntensity[j] * (accrualTime[j + 1] - accrualTime[j]);
  }

  size_t m = findInterval1(nsubjects, p) - 1;
  double t = accrualTime[m] + (nsubjects - p[m]) / accrualIntensity[m];
  return t;
}


//' @title Accrual Duration to Enroll Target Number of Subjects
//' @description Obtains the accrual duration to enroll the target number
//' of subjects.
//'
//' @param nsubjects The vector of target number of subjects.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//'
//' @return A vector of accrual durations.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' getAccrualDurationFromN(nsubjects = c(20, 150), accrualTime = c(0, 3),
//'                         accrualIntensity = c(10, 20))
//'
//' @export
// [[Rcpp::export]]
std::vector<double> getAccrualDurationFromN(
    const std::vector<double>& nsubjects,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity) {

  size_t I = nsubjects.size();
  std::vector<double> t(I);
  for (size_t i = 0; i < I; ++i) {
    t[i] = getAccrualDurationFromN1(nsubjects[i], accrualTime, accrualIntensity);
  }

  return t;
}


//' @title Probability of Being at Risk
//' @description Obtains the probability of being at risk at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//'   probability of being at risk.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of being at risk at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' patrisk(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'         lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> patrisk(
    const std::vector<double>& time,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma) {

  size_t J = piecewiseSurvivalTime.size();
  auto lambdax = expand1(lambda, J, "lambda");
  auto gammax = expand1(gamma, J, "gamma");

  // Compute at-risk probability for each time point
  size_t k = time.size();
  std::vector<double> p(k);
  for (size_t i = 0; i < k; ++i) {
    p[i] = patrisk1(time[i], piecewiseSurvivalTime, lambdax, gammax);
  }

  return p;
}


//' @title Probability of Having an Event
//' @description Obtains the probability of having an event at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//'   probability of having an event.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of having an event at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' pevent(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'        lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> pevent(
    const std::vector<double>& time,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma) {

  size_t J = piecewiseSurvivalTime.size();
  auto lambdax = expand1(lambda, J, "lambda");
  auto gammax = expand1(gamma, J, "gamma");

  // Compute cumulative hazard contributions for each time point
  size_t k = time.size();
  std::vector<double> a(k);
  for (size_t i = 0; i < k; ++i) {
    a[i] = pevent1(time[i], piecewiseSurvivalTime, lambdax, gammax);
  }

  return a;
}


//' @title Number of Subjects at Risk
//' @description Obtains the number of subjects at risk at given analysis
//' times for each treatment group.
//'
//' @param t A vector of analysis times at which to calculate the number
//'   of patients at risk.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_maxFollowupTime
//' @param time Calendar time for the analysis.
//'
//' @return A matrix of the number of patients at risk at the specified
//' analysis times (row) for each treatment group (column).
//'
//' @details For a given treatment group \eqn{g} and calendar time \eqn{\tau},
//' the number of patients at risk at analysis time \eqn{t} is calculated as
//' \deqn{\phi_g A(\tau - t) S_g(t) G_g(t),} where \eqn{\phi_g} is the
//' probability of randomization to treatment group \eqn{g},
//' \eqn{A(\tau - t)} is the number of patients enrolled by calendar time
//' \eqn{\tau - t}, \eqn{S_g(t)G_g(t)} is the probability of being at risk at
//' analysis time \eqn{t} for a patient in treatment group \eqn{g}
//' after enrollment. Obviously, \eqn{t < \min(\tau, T_{\rm{fmax}})}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' natrisk(t = c(9, 24), allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, maxFollowupTime = 30, time = 30)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix natrisk(
    const Rcpp::NumericVector& t = NA_REAL,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double maxFollowupTime = NA_REAL,
    const double time = NA_REAL) {

  auto t1 = Rcpp::as<std::vector<double>>(t);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  size_t J = pwSurvT.size();
  auto lambda1x = expand1(lam1, J, "lambda1");
  auto lambda2x = expand1(lam2, J, "lambda2");
  auto gamma1x = expand1(gam1, J, "gamma1");
  auto gamma2x = expand1(gam2, J, "gamma2");

  size_t k = t1.size();
  FlatMatrix n(k, 2);
  for (size_t i = 0; i < k; ++i) {
    auto result = natrisk1(
      t1[i], allocationRatioPlanned, accrualT, accrualInt,
      pwSurvT, lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, maxFollowupTime, time
    );

    n(i, 0) = result.first;  // Patients at risk in active treatment
    n(i, 1) = result.second; // Patients at risk in control
  }

  return Rcpp::wrap(n);
}



//' @title Number of Subjects Having an Event by Calendar Time
//' @description Obtains the number of subjects having an event by given
//' calendar times for each treatment group.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of patients having an event.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients having an event at the
//' specified calendar times (row) for each treatment group (column).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @details For a given treatment group \eqn{g} and calendar time \eqn{\tau},
//' the number of patients having an event by calendar time \eqn{\tau} is
//' calculated as \eqn{I_1 + I_2}, where
//' \deqn{I_1 = \phi_g A(\tau - T_{\rm{fmax}}) P_g(T_{\rm{fmax}}),} and
//' \deqn{I_2 = \phi_g \int_{\tau - T_{\rm{fmax}}}^{\tau} a(u) P_g(\tau - u) du,}
//' where \eqn{\phi_g} is the probability of randomization to treatment group \eqn{g},
//' \eqn{A(\tau - T_{\rm{fmax}})} is the number of patients enrolled by
//' calendar time \eqn{\tau - T_{\rm{fmax}}},
//' \eqn{P_g(T_{\rm{fmax}})} is the probability of having an event by
//' the maximum follow-up time \eqn{T_{\rm{fmax}}} for a patient in
//' treatment group \eqn{g} after enrollment,
//' \eqn{a(u)} is the accrual intensity at calendar time \eqn{u},
//' and \eqn{P_g(\tau - u)} is the probability of having an event by
//' calendar time \eqn{\tau} for a patient in treatment group \eqn{g}
//' enrolled at calendar time \eqn{u}.
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//' nevent(time = c(9, 24), allocationRatioPlanned = 1,
//'        accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'        piecewiseSurvivalTime = c(0, 6),
//'        lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'        gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'        accrualDuration = 12, maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix nevent(
    const Rcpp::NumericVector& time = NA_REAL,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double maxFollowupTime = NA_REAL) {

  auto time1 = Rcpp::as<std::vector<double>>(time);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  size_t J = pwSurvT.size();
  auto lambda1x = expand1(lam1, J, "lambda1");
  auto lambda2x = expand1(lam2, J, "lambda2");
  auto gamma1x = expand1(gam1, J, "gamma1");
  auto gamma2x = expand1(gam2, J, "gamma2");

  size_t k = time1.size();
  FlatMatrix d(k, 2);
  for (size_t i = 0; i < k; ++i) {
    auto result = nevent1(
      time1[i], allocationRatioPlanned, accrualT, accrualInt,
      pwSurvT, lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, maxFollowupTime
    );

    d(i, 0) = result.first;  // Active treatment group
    d(i, 1) = result.second; // Control group
  }

  return Rcpp::wrap(d);
}
