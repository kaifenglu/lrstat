#include "enrollment_event.h"
#include "utilities.h"
#include "dataframe_list.h"

#include <algorithm>    // max, min
#include <cmath>        // exp
#include <cstddef>      // size_t
#include <limits>       // numeric_limits
#include <stdexcept>    // invalid_argument
#include <vector>       // vector

#include <Rcpp.h>


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
std::vector<double> accrual(const std::vector<double>& time,
                            const std::vector<double>& accrualTime,
                            const std::vector<double>& accrualIntensity,
                            const double accrualDuration) {

  int k = static_cast<int>(time.size());
  std::vector<double> n(k);

  // up to end of enrollment
  std::vector<double> t(k);
  for (int i = 0; i < k; ++i) {
    t[i] = std::max(std::min(time[i], accrualDuration), 0.0);
  }

  // identify the time interval containing t
  std::vector<int> m = findInterval3(t, accrualTime);

  // sum up patients enrolled in each interval up to t
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < m[i]; ++j) {
      if (j < m[i] - 1) {
        n[i] += accrualIntensity[j] * (accrualTime[j + 1] - accrualTime[j]);
      } else {
        n[i] += accrualIntensity[j] * (t[i] - accrualTime[j]);
      }
    }
  }

  return n;
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
  int I = static_cast<int>(nsubjects.size());
  int J = static_cast<int>(accrualTime.size());
  std::vector<double> t(I), p(J);

  p[0] = 0;
  for (int j = 0; j < J - 1; ++j) {
    p[j+1] = p[j] + accrualIntensity[j] * (accrualTime[j+1] - accrualTime[j]);
  }

  std::vector<int> m = findInterval3(nsubjects, p);

  for (int i = 0; i < I; ++i) {
    int j = m[i] - 1;
    t[i] = accrualTime[j] + (nsubjects[i] - p[j]) / accrualIntensity[j];
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
std::vector<double> patrisk(const std::vector<double>& time,
                            const std::vector<double>& piecewiseSurvivalTime,
                            const std::vector<double>& lambda,
                            const std::vector<double>& gamma) {
  std::size_t k = time.size();
  std::size_t J = piecewiseSurvivalTime.size();

  // Validate and replicate lambda and gamma
  std::vector<double> lambdax(J), gammax(J);

  if (lambda.size() == 1) {
    lambdax = std::vector<double>(J, lambda[0]);
  } else if (lambda.size() == J) {
    lambdax = lambda;
  } else {
    throw std::invalid_argument("Invalid length for lambda");
  }

  if (gamma.size() == 1) {
    gammax = std::vector<double>(J, gamma[0]);
  } else if (gamma.size() == J) {
    gammax = gamma;
  } else {
    throw std::invalid_argument("Invalid length for gamma");
  }

  // Cumulative hazards for lambda + gamma
  std::vector<double> lamgam(J);
  for (std::size_t j = 0; j < J; ++j) {
    lamgam[j] = lambdax[j] + gammax[j];
  }

  std::vector<double> cumulativeRisk(k, 0.0);

  // Find intervals containing specified analysis time
  std::vector<int> m = findInterval3(time, piecewiseSurvivalTime);

  // Compute cumulative hazard for each time point
  for (std::size_t i = 0; i < k; ++i) {
    double a = 0.0;  // Hazard accumulator for this time point
    for (int j = 0; j < m[i]; ++j) {
      if (j < m[i] - 1) {
        // Contribution from intervals fully covered by time[i]
        a += lamgam[j] * (piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j]);
      } else {
        // Contribution from the remaining portion of the last interval
        a += lamgam[j] * (time[i] - piecewiseSurvivalTime[j]);
      }
    }
    cumulativeRisk[i] = std::exp(-a);  // Apply exponential decay
  }

  return cumulativeRisk;
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
std::vector<double> pevent(const std::vector<double>& time,
                           const std::vector<double>& piecewiseSurvivalTime,
                           const std::vector<double>& lambda,
                           const std::vector<double>& gamma) {
  std::size_t k = time.size();
  std::size_t J = piecewiseSurvivalTime.size();

  // Validate and replicate lambda and gamma
  std::vector<double> lambdax(J), gammax(J);

  if (lambda.size() == 1) {
    lambdax = std::vector<double>(J, lambda[0]);
  } else if (lambda.size() == J) {
    lambdax = lambda;
  } else {
    throw std::invalid_argument("Invalid length for lambda");
  }

  if (gamma.size() == 1) {
    gammax = std::vector<double>(J, gamma[0]);
  } else if (gamma.size() == J) {
    gammax = gamma;
  } else {
    throw std::invalid_argument("Invalid length for gamma");
  }

  // Compute lambda + gamma
  std::vector<double> lamgam(J);
  for (std::size_t j = 0; j < J; ++j) {
    lamgam[j] = lambdax[j] + gammax[j];
  }

  // Get risk of patients up to each time interval
  std::vector<double> n = patrisk(piecewiseSurvivalTime, piecewiseSurvivalTime,
                                  lambda, gamma);

  std::vector<int> m = findInterval3(time, piecewiseSurvivalTime);
  std::vector<double> a(k, 0.0);

  // Compute cumulative hazard contributions for each time point
  for (std::size_t i = 0; i < k; ++i) {
    double ai = 0.0;  // Accumulator for this time point
    for (int j = 0; j < m[i]; ++j) {
      double p;
      if (j < m[i] - 1) {
        // Full interval is covered
        p = lambdax[j] / lamgam[j] * (1.0 - std::exp(-lamgam[j] *
          (piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j])));
      } else {
        // Partial interval is covered
        p = lambdax[j] / lamgam[j] * (1.0 - std::exp(-lamgam[j] *
          (time[i] - piecewiseSurvivalTime[j])));
      }
      ai += n[j] * p;  // Add risk-weighted probability
    }
    a[i] = ai;
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
double hd(const int j,
          const double t1,
          const double t2,
          const std::vector<double>& piecewiseSurvivalTime,
          const std::vector<double>& lambda,
          const std::vector<double>& gamma) {

  int j1 = j-1;

  // lower bound of time interval j for piecewise exponential distribution
  double t0 = piecewiseSurvivalTime[j1];
  std::vector<double> t0_vec = {t0};

  // probability of being at risk at the start of interval j
  double n0 = patrisk(t0_vec, piecewiseSurvivalTime, lambda, gamma)[0];

  // Compute probability of having an event at the start of interval j
  double d0 = pevent(t0_vec, piecewiseSurvivalTime, lambda, gamma)[0];

  // Prepare lambda and gamma vectors
  std::size_t J = piecewiseSurvivalTime.size();
  std::vector<double> lambdax(J), gammax(J);

  if (lambda.size() == 1) {
    lambdax = std::vector<double>(J, lambda[0]);
  } else if (lambda.size() == J) {
    lambdax = lambda;
  } else {
    throw std::invalid_argument("Invalid length for lambda");
  }

  if (gamma.size() == 1) {
    gammax = std::vector<double>(J, gamma[0]);
  } else if (gamma.size() == J) {
    gammax = gamma;
  } else {
    throw std::invalid_argument("Invalid length for gamma");
  }

  // Compute total hazard (lambda + gamma)
  std::vector<double> lamgam(J);
  for (std::size_t i = 0; i < J; ++i) {
    lamgam[i] = lambdax[i] + gammax[i];
  }

  // Integration for conditional probability over (t1, t2)
  double lamgam_j1 = lamgam[j1];
  double exp1 = std::exp(-lamgam_j1 * (t1 - t0));
  double exp2 = std::exp(-lamgam_j1 * (t2 - t0));
  double q1 = (exp1 - exp2) / lamgam_j1;
  double q = lambdax[j1] / lamgam_j1 * (t2 - t1 - q1);

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
  std::vector<double> t12 = {t1, t2};
  std::vector<int> j12 = findInterval3(t12, piecewiseSurvivalTime);

  int j1 = std::max(j12[0] - 1, 0);  // Ensure index is not less than 0
  int j2 = std::max(j12[1] - 1, 0);

  double a = 0.0;

  // Sum up the integrated event probabilities across analysis time intervals
  for (int j = j1; j <= j2; ++j) {
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
//' @param time A vector of calendar times at which to calculate the number
//'   of patients having an event.
//' @param u1 Lower bound of the accrual time interval.
//' @param u2 Upper bound of the accrual time interval.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of number of patients who are enrolled during a
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
//' ad(time = c(9, 15), u1 = 1, u2 = 8, accrualTime = c(0, 3),
//'    accrualIntensity = c(10, 20), piecewiseSurvivalTime=c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> ad(const std::vector<double>& time,
                       const double u1,
                       const double u2,
                       const std::vector<double>& accrualTime,
                       const std::vector<double>& accrualIntensity,
                       const std::vector<double>& piecewiseSurvivalTime,
                       const std::vector<double>& lambda,
                       const std::vector<double>& gamma) {

  // Identify accrual time intervals containing u1 and u2
  std::vector<double> u12 = {u1, u2};
  std::vector<int> j12 = findInterval3(u12, accrualTime);
  int j1 = std::max(j12[0] - 1, 0);  // 0-based index for j1
  int j2 = std::max(j12[1] - 1, 0);  // 0-based index for j2

  std::size_t k = time.size();
  std::vector<double> a(k, 0.0);  // Initialize the result vector with zeroes

  // Sum up the number of patients with an event across accrual time intervals
  for (std::size_t i = 0; i < k; ++i) {
    double t = time[i];  // Current time
    for (int j = j1; j <= j2; ++j) {
      double x = 0.0;
      // Check intervals
      if (j1 == j2) {
        // Both u1 and u2 are in the same interval
        x = pd(t - u2, t - u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j == j1) {
        // First interval
        x = pd(t - accrualTime[j + 1], t - u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j == j2) {
        // Last interval
        x = pd(t - u2, t - accrualTime[j], piecewiseSurvivalTime, lambda, gamma);
      } else {
        // Intermediate intervals
        x = pd(t - accrualTime[j + 1], t - accrualTime[j],
               piecewiseSurvivalTime, lambda, gamma);
      }
      // Add the contribution from this interval
      a[i] += accrualIntensity[j] * x;
    }
  }

  return a;
}


FlatMatrix natriskcpp(const std::vector<double>& time,
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

  std::size_t k = time.size();

  // truncate the analysis time by the maximum follow-up
  std::vector<double> t(k), u(k);
  for (std::size_t i = 0; i < k; ++i) {
    t[i] = std::min(time[i], maxFollowupTime);
    u[i] = std::min(accrualDuration + minFollowupTime - t[i], accrualDuration);
  }

  // Number of patients enrolled
  std::vector<double> a = accrual(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  FlatMatrix n(k, 2);  // FlatMatrix with column-major storage

  // Compute probabilities for the active and control groups
  std::vector<double> patrisk1 = patrisk(t, piecewiseSurvivalTime, lambda1, gamma1);
  std::vector<double> patrisk2 = patrisk(t, piecewiseSurvivalTime, lambda2, gamma2);

  // Compute values for FlatMatrix directly
  for (std::size_t i = 0; i < k; ++i) {
    n(i, 0) = phi * a[i] * patrisk1[i];  // Patients at risk in active treatment
    n(i, 1) = (1.0 - phi) * a[i] * patrisk2[i];  // Patients at risk in control
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

  auto timecpp = Rcpp::as<std::vector<double>>(time);
  auto accrualTimecpp = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualIntensitycpp = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTimecpp = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lambda1cpp = Rcpp::as<std::vector<double>>(lambda1);
  auto lambda2cpp = Rcpp::as<std::vector<double>>(lambda2);
  auto gamma1cpp = Rcpp::as<std::vector<double>>(gamma1);
  auto gamma2cpp = Rcpp::as<std::vector<double>>(gamma2);

  auto cpp_result = natriskcpp(
    timecpp, allocationRatioPlanned, accrualTimecpp, accrualIntensitycpp,
    piecewiseTimecpp, lambda1cpp, lambda2cpp, gamma1cpp, gamma2cpp,
    accrualDuration, minFollowupTime, maxFollowupTime
  );

  return Rcpp::wrap(cpp_result);
}


FlatMatrix neventcpp(const std::vector<double>& time,
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

  std::size_t k = time.size();

  // truncate the analysis time by the maximum follow-up
  std::vector<double> t(k), u(k);
  for (std::size_t i = 0; i < k; ++i) {
    t[i] = std::min(time[i], maxFollowupTime);
    u[i] = std::min(accrualDuration + minFollowupTime - t[i], accrualDuration);
  }

  // Number of patients enrolled
  std::vector<double> a = accrual(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Prepare FlatMatrix for results (k rows, 2 columns)
  FlatMatrix d(k, 2);  // FlatMatrix with column-major storage

  // Compute the probabilities of having events in both treatment groups
  std::vector<double> pevent1 = pevent(t, piecewiseSurvivalTime, lambda1, gamma1);
  std::vector<double> pevent2 = pevent(t, piecewiseSurvivalTime, lambda2, gamma2);

  // Constant for ad() calculations
  std::vector<double> u1(1, accrualDuration + minFollowupTime);

  // Compute the number of patients having an event in each group
  for (std::size_t i = 0; i < k; ++i) {
    double ad1 = ad(u1, u[i], accrualDuration, accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, lambda1, gamma1)[0];
    double ad2 = ad(u1, u[i], accrualDuration, accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, lambda2, gamma2)[0];

    double d1 = a[i] * pevent1[i];
    double d2 = a[i] * pevent2[i];

    d(i, 0) = phi * (d1 + ad1);  // Active treatment group
    d(i, 1) = (1.0 - phi) * (d2 + ad2);  // Control group
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

  auto timecpp = Rcpp::as<std::vector<double>>(time);
  auto accrualTimecpp = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualIntensitycpp = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTimecpp = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lambda1cpp = Rcpp::as<std::vector<double>>(lambda1);
  auto lambda2cpp = Rcpp::as<std::vector<double>>(lambda2);
  auto gamma1cpp = Rcpp::as<std::vector<double>>(gamma1);
  auto gamma2cpp = Rcpp::as<std::vector<double>>(gamma2);

  auto cpp_result = neventcpp(
    timecpp, allocationRatioPlanned, accrualTimecpp, accrualIntensitycpp,
    piecewiseTimecpp, lambda1cpp, lambda2cpp, gamma1cpp, gamma2cpp,
    accrualDuration, minFollowupTime, maxFollowupTime
  );

  return Rcpp::wrap(cpp_result);
}


FlatMatrix nevent2cpp(const std::vector<double>& time,
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

  std::size_t k = time.size();

  // truncate the analysis time by the maximum follow-up
  std::vector<double> t(k), u(k), w(k);
  for (std::size_t i = 0; i < k; ++i) {
    t[i] = std::min(time[i], accrualDuration + minFollowupTime);
    u[i] = std::min(std::max(t[i] - maxFollowupTime, 0.0), accrualDuration);
    w[i] = std::min(t[i], accrualDuration);
  }

  // Number of patients enrolled
  std::vector<double> a = accrual(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Prepare FlatMatrix for results (k rows, 2 columns)
  FlatMatrix d(k, 2);  // FlatMatrix with column-major storage

  // Precompute probabilities using pevent
  std::vector<double> s(1, maxFollowupTime); // s contains maxFollowupTime
  std::vector<double> d1 = a; // Copy of a
  std::vector<double> d2 = a; // Copy of a

  double pevent1 = pevent(s, piecewiseSurvivalTime, lambda1, gamma1)[0];
  double pevent2 = pevent(s, piecewiseSurvivalTime, lambda2, gamma2)[0];

  for (std::size_t i = 0; i < k; ++i) {
    d1[i] *= pevent1;
    d2[i] *= pevent2;
  }

  // Compute the number of patients experiencing events in each group
  for (std::size_t i = 0; i < k; ++i) {
    std::vector<double> v(1, t[i]);
    double ad1 = ad(v, u[i], w[i], accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, lambda1, gamma1)[0];
    double ad2 = ad(v, u[i], w[i], accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, lambda2, gamma2)[0];

    d(i, 0) = phi * (d1[i] + ad1);  // Active treatment group
    d(i, 1) = (1.0 - phi) * (d2[i] + ad2);  // Control group
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

  auto timecpp = Rcpp::as<std::vector<double>>(time);
  auto accrualTimecpp = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualIntensitycpp = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTimecpp = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lambda1cpp = Rcpp::as<std::vector<double>>(lambda1);
  auto lambda2cpp = Rcpp::as<std::vector<double>>(lambda2);
  auto gamma1cpp = Rcpp::as<std::vector<double>>(gamma1);
  auto gamma2cpp = Rcpp::as<std::vector<double>>(gamma2);

  auto cpp_result = nevent2cpp(
    timecpp, allocationRatioPlanned, accrualTimecpp, accrualIntensitycpp,
    piecewiseTimecpp, lambda1cpp, lambda2cpp, gamma1cpp, gamma2cpp,
    accrualDuration, minFollowupTime, maxFollowupTime
  );

  return Rcpp::wrap(cpp_result);
}


