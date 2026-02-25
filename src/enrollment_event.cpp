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
  size_t m = std::max(findInterval1(t, accrualTime), 1);

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


double patrisk1(
    const double time,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma) {

  size_t J = piecewiseSurvivalTime.size();

  // Cumulative hazards for lambda + gamma
  std::vector<double> lamgam(J);
  for (size_t j = 0; j < J; ++j) {
    lamgam[j] = lambda[j] + gamma[j];
  }

  // Find interval containing specified analysis time
  size_t m = std::max(findInterval1(time, piecewiseSurvivalTime), 1);

  // Compute cumulative hazard for the time point
  double a = 0.0;
  for (size_t j = 0; j < m; ++j) {
    if (j < m - 1) {
      // Contribution from intervals fully covered by time[i]
      a += lamgam[j] * (piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j]);
    } else {
      // Contribution from the remaining portion of the last interval
      a += lamgam[j] * (time - piecewiseSurvivalTime[j]);
    }
  }

  return std::exp(-a);
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


double pevent1(
    const double time,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma) {

  size_t J = piecewiseSurvivalTime.size();

  // Compute lambda + gamma
  std::vector<double> lamgam(J);
  for (size_t j = 0; j < J; ++j) {
    lamgam[j] = lambda[j] + gamma[j];
  }

  // Get risk of patients up to each time interval
  auto n = patrisk(piecewiseSurvivalTime, piecewiseSurvivalTime, lambda, gamma);

  size_t m = std::max(findInterval1(time, piecewiseSurvivalTime), 1);

  double a = 0;
  for (size_t j = 0; j < m; ++j) {
    double p;
    if (j < m - 1) {
      // Full interval is covered
      p = lambda[j] / lamgam[j] * (1.0 - std::exp(-lamgam[j] *
        (piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j])));
    } else {
      // Partial interval is covered
      p = lambda[j] / lamgam[j] * (1.0 - std::exp(-lamgam[j] *
        (time - piecewiseSurvivalTime[j])));
    }

    a += n[j] * p;  // Add risk-weighted probability
  }

  return a;
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


//' @title Integrated Event Probability Over an Interval With Constant Hazard
//' @description Obtains the integrated probability of having an event
//' during an interval with constant hazard.
//'
//' @param j The analysis time interval with constant hazard.
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A value for the integrated probability of having an event
//' during an interval with constant hazard for a treatment
//' group with specified piecewise exponential survival and dropout
//' distributions.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' hd(j = 1, t1 = 1, t2 = 3, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double hd(const size_t j,
          const double t1,
          const double t2,
          const std::vector<double>& piecewiseSurvivalTime,
          const std::vector<double>& lambda,
          const std::vector<double>& gamma) {

  size_t j1 = j - 1;

  // lower bound of time interval j for piecewise exponential distribution
  double t0 = piecewiseSurvivalTime[j1];
  std::vector<double> t0_vec = {t0};

  // probability of being at risk at the start of interval j
  double n0 = patrisk1(t0, piecewiseSurvivalTime, lambda, gamma);

  // Compute probability of having an event at the start of interval j
  double d0 = pevent1(t0, piecewiseSurvivalTime, lambda, gamma);

  // Compute total hazard (lambda + gamma)
  size_t J = piecewiseSurvivalTime.size();
  std::vector<double> lamgam(J);
  for (size_t i = 0; i < J; ++i) {
    lamgam[i] = lambda[i] + gamma[i];
  }

  // Integration for conditional probability over (t1, t2)
  double lamgam_j1 = lamgam[j1];
  double exp1 = std::exp(-lamgam_j1 * (t1 - t0));
  double exp2 = std::exp(-lamgam_j1 * (t2 - t0));
  double q1 = (exp1 - exp2) / lamgam_j1;
  double q = lambda[j1] / lamgam_j1 * (t2 - t1 - q1);

  // Sum up the integration for already failed and to-be-failed
  return d0 * (t2 - t1) + n0 * q;
}


//' @title Integrated Event Probability Over an Interval
//' @description Obtains the integration of the probability of having an
//' event during an interval. The specified analysis time interval can span
//' more than one analysis time interval with constant hazard.
//'
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A value for the integrated probability of having an event
//' during an interval for a treatment group with specified
//' piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' pd(t1 = 1, t2 = 8, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double pd(const double t1,
          const double t2,
          const std::vector<double>& piecewiseSurvivalTime,
          const std::vector<double>& lambda,
          const std::vector<double>& gamma) {

  // Identify analysis time intervals containing t1 and t2
  size_t j1 = std::max(findInterval1(t1, piecewiseSurvivalTime) - 1, 0);
  size_t j2 = std::max(findInterval1(t2, piecewiseSurvivalTime) - 1, 0);

  double a = 0.0;

  // Sum up the integrated event probabilities across analysis time intervals
  for (size_t j = j1; j <= j2; ++j) {
    double x = 0.0;
    if (j1 == j2) {
      // Both t1 and t2 are in the same interval
      x = hd(j + 1, t1, t2, piecewiseSurvivalTime, lambda, gamma);
    } else if (j == j1) {
      // First interval
      x = hd(j + 1, t1, piecewiseSurvivalTime[j + 1],
             piecewiseSurvivalTime, lambda, gamma);
    } else if (j == j2) {
      // Last interval
      x = hd(j + 1, piecewiseSurvivalTime[j], t2,
             piecewiseSurvivalTime, lambda, gamma);
    } else {
      // Intermediate intervals
      x = hd(j + 1, piecewiseSurvivalTime[j], piecewiseSurvivalTime[j + 1],
             piecewiseSurvivalTime, lambda, gamma);
    }
    a += x;
  }

  return a;
}


//' @title Number of Patients Enrolled During an Interval and Having an Event
//' by Specified Calendar Times
//' @description Obtains the number of patients who are enrolled during a
//' specified enrollment time interval and have an event by the specified
//' calendar times.
//'
//' @param time Calendar time at which to calculate the number
//'   of patients having an event.
//' @param u1 Lower bound of the accrual time interval.
//' @param u2 Upper bound of the accrual time interval.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return Number of patients who are enrolled during a
//' specified enrollment time interval and have an event by the specified
//' calendar times for a given treatment group had the enrollment being
//' restricted to the treatment group. By definition, we must have
//' \code{time >= u2}.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, 10 patients per month for the first 3 months, and
//' # 20 patients per month thereafter. Piecewise exponential survival with
//' # hazard 0.0533 in the first 6 months, and hazard 0.0309 thereafter,
//' # and 5% dropout by the end of 1 year.
//'
//' ad(time = 9, u1 = 1, u2 = 8, accrualTime = c(0, 3),
//'    accrualIntensity = c(10, 20), piecewiseSurvivalTime=c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double ad(const double time,
          const double u1,
          const double u2,
          const std::vector<double>& accrualTime,
          const std::vector<double>& accrualIntensity,
          const std::vector<double>& piecewiseSurvivalTime,
          const std::vector<double>& lambda,
          const std::vector<double>& gamma) {

  // Identify accrual time intervals containing u1 and u2
  size_t j1 = std::max(findInterval1(u1, accrualTime) - 1, 0);
  size_t j2 = std::max(findInterval1(u2, accrualTime) - 1, 0);

  double a = 0.0;  // Initialize the result with zero

  // Sum up the number of patients with an event across accrual time intervals
  for (size_t j = j1; j <= j2; ++j) {
    double x = 0.0;
    // Check intervals
    if (j1 == j2) {
      // Both u1 and u2 are in the same interval
      x = pd(time - u2, time - u1, piecewiseSurvivalTime, lambda, gamma);
    } else if (j == j1) {
      // First interval
      x = pd(time - accrualTime[j + 1], time - u1, piecewiseSurvivalTime,
             lambda, gamma);
    } else if (j == j2) {
      // Last interval
      x = pd(time - u2, time - accrualTime[j], piecewiseSurvivalTime,
             lambda, gamma);
    } else {
      // Intermediate intervals
      x = pd(time - accrualTime[j + 1], time - accrualTime[j],
             piecewiseSurvivalTime, lambda, gamma);
    }
    // Add the contribution from this interval
    a += accrualIntensity[j] * x;
  }

  return a;
}


std::pair<double, double> natrisk1cpp(
    const double time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime) {

  // truncate the analysis time by the maximum follow-up
  double t = std::min(time, maxFollowupTime);
  double u = std::min(accrualDuration + minFollowupTime - t, accrualDuration);

  // Number of patients enrolled
  double a = accrual1(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Compute probabilities for the active and control groups
  double p1 = patrisk1(t, piecewiseSurvivalTime, lambda1, gamma1);
  double p2 = patrisk1(t, piecewiseSurvivalTime, lambda2, gamma2);

  // Compute values for FlatMatrix directly
  double n1 = phi * a * p1;  // Patients at risk in active treatment
  double n2 = (1.0 - phi) * a * p2;  // Patients at risk in control

  return std::make_pair(n1, n2);
}


FlatMatrix natriskcpp(
    const std::vector<double>& time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime) {

  size_t J = piecewiseSurvivalTime.size();
  auto lambda1x = expand1(lambda1, J, "lambda1");
  auto lambda2x = expand1(lambda2, J, "lambda2");
  auto gamma1x = expand1(gamma1, J, "gamma1");
  auto gamma2x = expand1(gamma2, J, "gamma2");

  size_t k = time.size();
  FlatMatrix n(k, 2);
  for (size_t i = 0; i < k; ++i) {
    auto result = natrisk1cpp(
      time[i], allocationRatioPlanned, accrualTime, accrualIntensity,
      piecewiseSurvivalTime, lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, minFollowupTime, maxFollowupTime
    );

    n(i, 0) = result.first;  // Patients at risk in active treatment
    n(i, 1) = result.second; // Patients at risk in control
  }

  return n;
}


//' @title Number of Subjects at Risk
//' @description Obtains the number of subjects at risk at given analysis
//' times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
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
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients at risk at the specified
//' analysis times (row) for each treatment group (column).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' natrisk(time = c(9, 24), allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, minFollowupTime = 18,
//'         maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix natrisk(
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
    const double minFollowupTime = NA_REAL,
    const double maxFollowupTime = NA_REAL) {

  auto time1 = Rcpp::as<std::vector<double>>(time);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  auto result = natriskcpp(
    time1, allocationRatioPlanned, accrualT, accrualInt,
    pwSurvT, lam1, lam2, gam1, gam2,
    accrualDuration, minFollowupTime, maxFollowupTime
  );

  return Rcpp::wrap(result);
}


std::pair<double, double> nevent1cpp(
    const double time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime) {

  // truncate the analysis time by the maximum follow-up
  double t = std::min(time, maxFollowupTime);
  double u = std::min(accrualDuration + minFollowupTime - t, accrualDuration);

  // Number of patients enrolled
  double a = accrual1(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Compute the probabilities of having events in both treatment groups
  double p1 = pevent1(t, piecewiseSurvivalTime, lambda1, gamma1);
  double p2 = pevent1(t, piecewiseSurvivalTime, lambda2, gamma2);

  // Compute the number of patients having an event in each group
  double u1 = accrualDuration + minFollowupTime;
  double c1 = ad(u1, u, accrualDuration, accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, lambda1, gamma1);
  double c2 = ad(u1, u, accrualDuration, accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, lambda2, gamma2);

  double b1 = a * p1;
  double b2 = a * p2;

  double d1 = phi * (b1 + c1);  // Active treatment group
  double d2 = (1.0 - phi) * (b2 + c2);  // Control group

  return std::make_pair(d1, d2);
}


FlatMatrix neventcpp(
    const std::vector<double>& time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime) {

  size_t J = piecewiseSurvivalTime.size();
  auto lambda1x = expand1(lambda1, J, "lambda1");
  auto lambda2x = expand1(lambda2, J, "lambda2");
  auto gamma1x = expand1(gamma1, J, "gamma1");
  auto gamma2x = expand1(gamma2, J, "gamma2");

  size_t k = time.size();
  FlatMatrix d(k, 2);
  for (size_t i = 0; i < k; ++i) {
    auto result = nevent1cpp(
      time[i], allocationRatioPlanned, accrualTime, accrualIntensity,
      piecewiseSurvivalTime, lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, minFollowupTime, maxFollowupTime
    );

    d(i, 0) = result.first;  // Number of events in active treatment
    d(i, 1) = result.second; // Number of events in control
  }

  return d;
}


//' @title Number of Subjects Having an Event
//' @description Obtains the number of subjects having an event by given
//' analysis times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
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
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients having an event at the
//' specified analysis times (row) for each treatment group (column).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' nevent(time = c(9, 24), allocationRatioPlanned = 1,
//'        accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'        piecewiseSurvivalTime = c(0, 6),
//'        lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'        gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'        accrualDuration = 12, minFollowupTime = 18,
//'        maxFollowupTime = 30)
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
    const double minFollowupTime = NA_REAL,
    const double maxFollowupTime = NA_REAL) {

  auto time1 = Rcpp::as<std::vector<double>>(time);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  auto result = neventcpp(
    time1, allocationRatioPlanned, accrualT, accrualInt,
    pwSurvT, lam1, lam2, gam1, gam2,
    accrualDuration, minFollowupTime, maxFollowupTime
  );

  return Rcpp::wrap(result);
}


std::pair<double, double> nevent21cpp(
    const double time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime) {

  // truncate the analysis time by the maximum follow-up
  double t = std::min(time, accrualDuration + minFollowupTime);
  double u = std::min(std::max(t - maxFollowupTime, 0.0), accrualDuration);
  double w = std::min(t, accrualDuration);

  // Number of patients enrolled
  double a = accrual1(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  double p1 = pevent1(maxFollowupTime, piecewiseSurvivalTime, lambda1, gamma1);
  double p2 = pevent1(maxFollowupTime, piecewiseSurvivalTime, lambda2, gamma2);

  double b1 = a * p1;
  double b2 = a * p2;

  // Compute the number of patients experiencing events in each group
  double c1 = ad(t, u, w, accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, lambda1, gamma1);
  double c2 = ad(t, u, w, accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, lambda2, gamma2);

  double d1 = phi * (b1 + c1);  // Active treatment group
  double d2 = (1.0 - phi) * (b2 + c2);  // Control group

  return std::make_pair(d1, d2);
}


FlatMatrix nevent2cpp(
    const std::vector<double>& time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime) {

  size_t J = piecewiseSurvivalTime.size();
  auto lambda1x = expand1(lambda1, J, "lambda1");
  auto lambda2x = expand1(lambda2, J, "lambda2");
  auto gamma1x = expand1(gamma1, J, "gamma1");
  auto gamma2x = expand1(gamma2, J, "gamma2");

  size_t k = time.size();
  FlatMatrix d(k, 2);
  for (size_t i = 0; i < k; ++i) {
    auto result = nevent21cpp(
      time[i], allocationRatioPlanned, accrualTime, accrualIntensity,
      piecewiseSurvivalTime, lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, minFollowupTime, maxFollowupTime
    );

    d(i, 0) = result.first;  // Active treatment group
    d(i, 1) = result.second; // Control group
  }

  return d;
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
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients having an event at the
//' specified calendar times (row) for each treatment group (column).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//' nevent2(time = c(9, 24), allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, minFollowupTime = 18,
//'         maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix nevent2(
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
    const double minFollowupTime = NA_REAL,
    const double maxFollowupTime = NA_REAL) {

  auto time1 = Rcpp::as<std::vector<double>>(time);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  auto result = nevent2cpp(
    time1, allocationRatioPlanned, accrualT, accrualInt,
    pwSurvT, lam1, lam2, gam1, gam2,
    accrualDuration, minFollowupTime, maxFollowupTime
  );

  return Rcpp::wrap(result);
}
