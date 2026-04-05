#pragma once

#include <vector>

#include "utilities.h"


double accrual1(
    const double time,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const double accrualDuration);

double getAccrualDurationFromN1(
    const double nsubjects,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity);


template <class VLam, class VGam>
inline double patrisk1(
    const double time,
    const std::vector<double>& piecewiseSurvivalTime,
    const VLam& lambda,
    const VGam& gamma) {

  const auto& t = piecewiseSurvivalTime;

  // Find interval containing specified analysis time
  std::size_t m = findInterval1(time, t);

  // Compute cumulative hazard for the time point
  double a = 0.0;
  // Contribution from intervals fully covered
  for (std::size_t j = 0; j < m-1; ++j) {
    a += (lambda[j] + gamma[j]) * (t[j+1] - t[j]);
  }
  // Contribution from the remaining portion of the last interval
  a += (lambda[m-1] + gamma[m-1]) * (time - t[m-1]);

  return std::exp(-a);
}


template <class VLam, class VGam>
inline double pevent1(
    const double time,
    const std::vector<double>& piecewiseSurvivalTime,
    const VLam& lambda,
    const VGam& gamma) {

  const auto& t = piecewiseSurvivalTime;

  std::size_t m = findInterval1(time, t);

  double a = 0;
  // Full interval is covered
  for (std::size_t j = 0; j < m-1; ++j) {
    double n = patrisk1(t[j], t, lambda, gamma);
    double theta = lambda[j] + gamma[j];
    double p = lambda[j] / theta * (1.0 - std::exp(-theta * (t[j + 1] - t[j])));
    a += n * p;  // Add risk-weighted probability
  }

  // Partial interval is covered
  double n = patrisk1(t[m-1], t, lambda, gamma);
  double theta = lambda[m-1] + gamma[m-1];
  double p = lambda[m-1] / theta * (1.0 - std::exp(-theta * (time - t[m-1])));
  a += n * p;

  return a;
}


// Integrate[D(t), {t, t1, t2}], for t1 and t2 in the j-th interval of
// piecewise exponential distribution
template <class VLam, class VGam>
inline double hd(
    const std::size_t j,
    const double t1,
    const double t2,
    const std::vector<double>& piecewiseSurvivalTime,
    const VLam& lambda,
    const VGam& gamma) {

  // lower bound of time interval j for piecewise exponential distribution
  double t0 = piecewiseSurvivalTime[j];

  // probability of being at risk at the start of interval j
  double n0 = patrisk1(t0, piecewiseSurvivalTime, lambda, gamma);

  // Compute probability of having an event at the start of interval j
  double d0 = pevent1(t0, piecewiseSurvivalTime, lambda, gamma);

  // Integration for conditional probability over (t1, t2)
  double theta = lambda[j] + gamma[j];
  double q1 = (std::exp(-theta * (t1 - t0)) - std::exp(-theta * (t2 - t0))) / theta;
  double q = lambda[j] / theta * (t2 - t1 - q1);

  // Sum up the integration for already failed and to-be-failed
  return d0 * (t2 - t1) + n0 * q;
}


// Integration of the probability of having an event during an interval.
// The specified analysis time interval can span more than one analysis
// time interval with constant hazard.
template <class VLam, class VGam>
inline double pd(
    const double t1,
    const double t2,
    const std::vector<double>& piecewiseSurvivalTime,
    const VLam& lambda,
    const VGam& gamma) {

  const std::vector<double>& t = piecewiseSurvivalTime;

  // Identify analysis time intervals containing t1 and t2
  std::size_t j1 = findInterval1(t1, t) - 1;
  std::size_t j2 = findInterval1(t2, t) - 1;

  double a = 0.0;

  // Sum up the integrated event probabilities across analysis time intervals
  if (j1 == j2) {
    // Both t1 and t2 are in the same interval
    a = hd(j1, t1, t2, t, lambda, gamma);
  } else {
    // First interval
    a = hd(j1, t1, t[j1 + 1], t, lambda, gamma);
    for (std::size_t j = j1 + 1; j < j2; ++j) {
      // Intermediate intervals
      a += hd(j, t[j], t[j + 1], t, lambda, gamma);
    }
    // Last interval
    a += hd(j2, t[j2], t2, t, lambda, gamma);
  }

  return a;
}


// Number of Patients Enrolled During an Interval and Having an Event
// by Specified Calendar Times
template <class VLam, class VGam>
inline double ad(
    const double time,
    const double u1,
    const double u2,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const VLam& lambda,
    const VGam& gamma) {

  const std::vector<double>& u = accrualTime;
  const std::vector<double>& t = piecewiseSurvivalTime;

  // Identify accrual time intervals containing u1 and u2
  std::size_t i1 = findInterval1(u1, u) - 1;
  std::size_t i2 = findInterval1(u2, u) - 1;

  double a = 0.0;  // Initialize the result with zero

  // Sum up the number of patients with an event across accrual time intervals
  if (i1 == i2) {
    // Both u1 and u2 are in the same interval
    a = accrualIntensity[i1] * pd(time - u2, time - u1, t, lambda, gamma);
  } else {
    // First interval
    a = accrualIntensity[i1] * pd(time - u[i1 + 1], time - u1, t, lambda, gamma);
    // intermediate intervals
    for (std::size_t i = i1 + 1; i < i2; ++i) {
      a += accrualIntensity[i] * pd(time - u[i + 1], time - u[i], t, lambda, gamma);
    }
    // Last interval
    a += accrualIntensity[i2] * pd(time - u2, time - u[i2], t, lambda, gamma);
  }

  return a;
}


template <class VLam1, class VLam2, class VGam1, class VGam2>
inline std::pair<double, double> natrisk1(
    const double t,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const VLam1& lambda1,
    const VLam2& lambda2,
    const VGam1& gamma1,
    const VGam2& gamma2,
    const double accrualDuration,
    const double maxFollowupTime,
    const double time) {

  // truncate the analysis time by the maximum follow-up
  double t1 = std::min(std::min(t, maxFollowupTime), time);
  double u = time - t1;

  // Number of patients enrolled
  double a = accrual1(u, accrualTime, accrualIntensity, accrualDuration);

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Compute probabilities for the active and control groups
  double p1 = patrisk1(t1, piecewiseSurvivalTime, lambda1, gamma1);
  double p2 = patrisk1(t1, piecewiseSurvivalTime, lambda2, gamma2);

  // Compute values for FlatMatrix directly
  double n1 = phi * a * p1;  // Patients at risk in active treatment
  double n2 = (1.0 - phi) * a * p2;  // Patients at risk in control

  return std::make_pair(n1, n2);
}


template <class VLam1, class VLam2, class VGam1, class VGam2>
inline std::pair<double, double> nevent1(
    const double time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const VLam1& lambda1,
    const VLam2& lambda2,
    const VGam1& gamma1,
    const VGam2& gamma2,
    const double accrualDuration,
    const double maxFollowupTime) {

  // Probability of randomization to the active treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // Number of patients enrolled by calendar time (time - maxFollowupTime)
  double a = accrual1(time - maxFollowupTime, accrualTime, accrualIntensity,
                      accrualDuration);

  double p1 = pevent1(maxFollowupTime, piecewiseSurvivalTime, lambda1, gamma1);
  double p2 = pevent1(maxFollowupTime, piecewiseSurvivalTime, lambda2, gamma2);

  // Calculate the number of patients enrolled during the interval
  // (time - maxFollowupTime, time) and having an event by calendar time time
  double u1 = std::min(std::max(time - maxFollowupTime, 0.0), accrualDuration);
  double u2 = std::min(time, accrualDuration);
  double c1 = ad(time, u1, u2, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, lambda1, gamma1);
  double c2 = ad(time, u1, u2, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, lambda2, gamma2);

  double d1 = phi * (a * p1 + c1);  // Active treatment group
  double d2 = (1.0 - phi) * (a * p2 + c2);  // Control group
  return std::make_pair(d1, d2);
}
