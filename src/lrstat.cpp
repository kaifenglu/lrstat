#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;

//' @title Number of enrolled subjects
//' @description Obtains the number of subjects enrolled by given calendar
//' times.
//'
//' @param time A vector of calendar times at which to calculate the number
//' of enrolled subjects.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_accrualDuration
//'
//' @return A vector of total number of subjects enrolled by the
//' specified calendar times.
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
NumericVector accrual(const NumericVector& time = NA_REAL,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const double accrualDuration = NA_REAL) {
  
  int i, j, k = time.size();
  NumericVector n(k);
  
  // up to end of enrollment
  NumericVector t = pmax(pmin(time, accrualDuration), 0.0);
  
  // identify the time interval containing t
  IntegerVector m = pmax(findInterval2(t, accrualTime), 1);
  
  // sum up patients enrolled in each interval up to t
  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        n[i] += accrualIntensity[j]*(accrualTime[j+1] - accrualTime[j]);
      } else {
        n[i] += accrualIntensity[j]*(t[i] - accrualTime[j]);
      }
    }
  }
  
  return n;
}

//' @title Probability of being at risk
//' @description Obtains the probability of being at risk at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//' probability of being at risk.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of being at risk at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6 months,
//' # and hazard 0.0309 thereafter, and 5% dropout by the end of 1 year.
//'
//' patrisk(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'         lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
NumericVector patrisk(const NumericVector& time = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda = NA_REAL,
                      const NumericVector& gamma = 0) {
  
  // identify the time interval containing the specified analysis time
  IntegerVector m = pmax(findInterval2(time, piecewiseSurvivalTime), 1);
  
  int i, j, k = time.size(), J = lambda.size();
  
  // hazard for failure or dropout
  NumericVector lg(J);
  if (gamma.size() == 1) {
    lg = lambda + gamma[0];
  } else {
    lg = lambda + gamma;
  }
  
  NumericVector t = piecewiseSurvivalTime;
  
  // sum up cumulative hazard up to time
  NumericVector a(k);
  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        a[i] += lg[j]*(t[j+1] - t[j]);
      } else {
        a[i] += lg[j]*(time[i] - t[j]);
      }
    }
  }
  
  return exp(-a);
}


//' @title Probability of having an event
//' @description Obtains the probability of having an event at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//' probability of having an event.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of having an event at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6 months,
//' # and hazard 0.0309 thereafter, and 5% dropout by the end of 1 year.
//'
//' pevent(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'        lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
NumericVector pevent(const NumericVector& time = NA_REAL,
                     const NumericVector& piecewiseSurvivalTime = 0,
                     const NumericVector& lambda = NA_REAL,
                     const NumericVector& gamma = 0) {
  
  // identify the time interval containing the specified analysis time
  IntegerVector m = pmax(findInterval2(time, piecewiseSurvivalTime), 1);
  
  int i, j, k = time.size(), J = lambda.size();
  
  // hazard for failure or dropout
  NumericVector lg(J);
  if (gamma.size() == 1) {
    lg = lambda + gamma[0];
  } else {
    lg = lambda + gamma;
  }
  
  // sum up cumulative hazard up to time
  NumericVector t = piecewiseSurvivalTime;
  NumericVector n = patrisk(t, t, lambda, gamma);
  NumericVector a(k);
  double p;
  
  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        p = lambda[j]/lg[j]*(1 - exp(-lg[j]*(t[j+1] - t[j])));
      } else {
        p = lambda[j]/lg[j]*(1 - exp(-lg[j]*(time[i] - t[j])));
      }
      a[i] += n[j]*p;
    }
  }
  
  return a;
}


//' @title Integrated event probability over an interval with constant hazard
//' @description Obtains the integration probability of having an event
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
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6 months,
//' # and hazard 0.0309 thereafter, and 5% dropout by the end of 1 year.
//'
//' hd(j = 1, t1 = 1, t2 = 3, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double hd(const int j = NA_INTEGER,
          const double t1 = NA_REAL,
          const double t2 = NA_REAL,
          const NumericVector& piecewiseSurvivalTime = 0,
          const NumericVector& lambda = NA_REAL,
          const NumericVector& gamma = 0) {
  
  int j1 = j-1;
  
  // lower bound of time interval j for the piecewise exponential distribution
  NumericVector t0 = NumericVector::create(piecewiseSurvivalTime[j1]);
  
  // probability of being at risk at the start of interval j
  NumericVector n0 = patrisk(t0, piecewiseSurvivalTime, lambda, gamma);
  
  // probability of having an event at the start of interval j
  NumericVector d0 = pevent(t0, piecewiseSurvivalTime, lambda, gamma);
  
  
  int J = lambda.size();
  
  // hazard for failure or dropout
  NumericVector lg(J);
  if (gamma.size() == 1) {
    lg = lambda + gamma[0];
  } else {
    lg = lambda + gamma;
  }
  
  // integration of conditional probability of having an event over (t1,t2)
  // given survival at the start of interval j
  double q1 = (exp(-lg[j1]*(t1-t0[0])) - exp(-lg[j1]*(t2-t0[0])))/lg[j1];
  double q = lambda[j1]/lg[j1] * (t2-t1 - q1);
  
  // sum up the integration for the already failed and to-be-failed
  return d0[0]*(t2-t1) + n0[0]*q;
}


//' @title Integrated event probability over an interval
//' @description Obtains the integration of the probability of having an event
//' during an interval. The specified analysis time interval can span more
//' than one analysis time interval with constant hazard.
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
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6 months,
//' # and hazard 0.0309 thereafter, and 5% dropout by the end of 1 year.
//'
//' pd(t1 = 1, t2 = 8, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double pd(const double t1 = NA_REAL,
          const double t2 = NA_REAL,
          const NumericVector& piecewiseSurvivalTime = 0,
          const NumericVector& lambda = NA_REAL,
          const NumericVector& gamma = 0) {
  
  // identify the analysis time intervals containing t1 and t2
  NumericVector t12 = NumericVector::create(t1, t2);
  IntegerVector j12 = pmax(findInterval2(t12, piecewiseSurvivalTime), 1) - 1;
  
  NumericVector t = piecewiseSurvivalTime;
  
  int j, j1=j12[0], j2=j12[1];
  
  // sum up the integrated event probabilities across analysis time intervals
  double a=0, x;
  for (j=j1; j<=j2; j++) {
    if (j1==j2) {
      x = hd(j+1, t1, t2, t, lambda, gamma);
    } else if (j==j1) {
      x = hd(j+1, t1, t[j+1], t, lambda, gamma);
    } else if (j==j2) {
      x = hd(j+1, t[j], t2, t, lambda, gamma);
    } else {
      x = hd(j+1, t[j], t[j+1], t, lambda, gamma);
    }
    a += x;
  }
  
  return a;
}


//' @title Number of patients enrolled during an interval and having an event
//' by specified calendar times
//' @description Obtains the number of patients who are enrolled during a
//' specified enrollment time interval and have an event by the specified
//' calendar times.
//'
//' @param time A vector of calendar times at which to calculate the number
//' of patients having an event.
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
NumericVector ad(const NumericVector& time = NA_REAL,
                 const double u1 = NA_REAL,
                 const double u2 = NA_REAL,
                 const NumericVector& accrualTime = 0,
                 const NumericVector& accrualIntensity = NA_REAL,
                 const NumericVector& piecewiseSurvivalTime = 0,
                 const NumericVector& lambda = NA_REAL,
                 const NumericVector& gamma = 0) {
  
  // identify the accrual time intervals containing u1 and u2
  NumericVector u12 = NumericVector::create(u1, u2);
  IntegerVector j12 = pmax(findInterval2(u12, accrualTime), 1) - 1;
  
  NumericVector u = accrualTime;
  
  int i, j, j1=j12[0], j2=j12[1], k=time.size();
  
  NumericVector a(k);
  
  // sum up the number of patients with event across accrual time intervals
  double t, x;
  for (i=0; i<k; i++) {
    t = time[i];
    for (j=j1; j<=j2; j++) {
      if (j1==j2) {
        x = pd(t-u2, t-u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j==j1) {
        x = pd(t-u[j+1], t-u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j==j2) {
        x = pd(t-u2, t-u[j], piecewiseSurvivalTime, lambda, gamma);
      } else {
        x = pd(t-u[j+1], t-u[j], piecewiseSurvivalTime, lambda, gamma);
      }
      a[i] += accrualIntensity[j]*x;
    }
  }
  
  return a;
}


//' @title Number of subjects at risk
//' @description Obtains the number of subjects at risk at given analysis
//' times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
//' of patients at risk.
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
//' analysis times for each treatment group.
//'
//' @keywords internal
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
NumericMatrix natrisk(const NumericVector& time = NA_REAL,
                      const double allocationRatioPlanned = 1,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda1 = NA_REAL,
                      const NumericVector& lambda2 = NA_REAL,
                      const NumericVector& gamma1 = 0,
                      const NumericVector& gamma2 = 0,
                      const double accrualDuration = NA_REAL,
                      const double minFollowupTime = NA_REAL,
                      const double maxFollowupTime = NA_REAL) {
  
  // truncate the analysis time by the maximum follow-up
  NumericVector t = pmin(time, maxFollowupTime);
  
  // enrollment time
  NumericVector u = pmin(accrualDuration+minFollowupTime-t, accrualDuration);
  
  // number of patients enrolled
  NumericVector a = accrual(u, accrualTime, accrualIntensity, accrualDuration);
  
  // probability of being randomized to the active treatment group
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);
  
  // number of patients at risk in each treatment group
  int k = time.size();
  NumericMatrix n(k, 2);
  n(_, 0) = phi*a*patrisk(t, piecewiseSurvivalTime, lambda1, gamma1);
  n(_, 1) = (1-phi)*a*patrisk(t, piecewiseSurvivalTime, lambda2, gamma2);
  
  return n;
}


//' @title Number of subjects having an event
//' @description Obtains the number of subjects having an event by given
//' analysis times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
//' of patients having an event.
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
//' @return A matrix of the number of patients having an event at the specified
//' analysis times for each treatment group.
//'
//' @keywords internal
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
NumericMatrix nevent(const NumericVector& time = NA_REAL,
                     const double allocationRatioPlanned = 1,
                     const NumericVector& accrualTime = 0,
                     const NumericVector& accrualIntensity = NA_REAL,
                     const NumericVector& piecewiseSurvivalTime = 0,
                     const NumericVector& lambda1 = NA_REAL,
                     const NumericVector& lambda2 = NA_REAL,
                     const NumericVector& gamma1 = 0,
                     const NumericVector& gamma2 = 0,
                     const double accrualDuration = NA_REAL,
                     const double minFollowupTime = NA_REAL,
                     const double maxFollowupTime = NA_REAL) {
  
  // truncate the analysis time by the maximum follow-up
  NumericVector t = pmin(time, maxFollowupTime);
  
  // enrollment time
  NumericVector u = pmin(accrualDuration+minFollowupTime-t, accrualDuration);
  
  // number of patients enrolled
  NumericVector a = accrual(u, accrualTime, accrualIntensity, accrualDuration);
  
  // probability of being randomized to the active treatment group
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);
  
  // number of patients having an event in each treatment group
  NumericVector u1(1);
  u1[0] = accrualDuration + minFollowupTime;
  
  int i, k = time.size();
  NumericMatrix d(k, 2);
  
  NumericVector d1(k), d2(k);
  d1 = a*pevent(t, piecewiseSurvivalTime, lambda1, gamma1);
  d2 = a*pevent(t, piecewiseSurvivalTime, lambda2, gamma2);
  
  for (i=0; i<k; i++) {
    d(i,0) = phi*(d1[i] + ad(u1, u[i], accrualDuration, accrualTime,
                  accrualIntensity, piecewiseSurvivalTime,
                  lambda1, gamma1)[0]);
    d(i,1) = (1-phi)*(d2[i] + ad(u1, u[i], accrualDuration, accrualTime,
              accrualIntensity, piecewiseSurvivalTime, lambda2, gamma2)[0]);
  }
  
  return d;
}


//' @title Number of subjects having an event by calendar time
//' @description Obtains the number of subjects having an event by given
//' calendar times for each treatment group.
//'
//' @param time A vector of calendar times at which to calculate the number
//' of patients having an event.
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
//' @return A matrix of the number of patients having an event at the specified
//' calendar times for each treatment group.
//'
//' @keywords internal
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
NumericMatrix nevent2(const NumericVector& time = NA_REAL,
                      const double allocationRatioPlanned = 1,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda1 = NA_REAL,
                      const NumericVector& lambda2 = NA_REAL,
                      const NumericVector& gamma1 = 0,
                      const NumericVector& gamma2 = 0,
                      const double accrualDuration = NA_REAL,
                      const double minFollowupTime = NA_REAL,
                      const double maxFollowupTime = NA_REAL) {
  
  // truncate the calendar time by study end
  NumericVector t = pmin(time, accrualDuration + minFollowupTime);
  
  // enrollment time
  NumericVector u = pmin(pmax(t - maxFollowupTime, 0.0), accrualDuration);
  NumericVector w = pmin(t, accrualDuration);
  
  // number of patients enrolled
  NumericVector a = accrual(u, accrualTime, accrualIntensity, accrualDuration);
  
  // probability of being randomized to the active treatment group
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);
  
  // number of patients having an event in each treatment group
  NumericVector s(1), v(1);
  s[0] = maxFollowupTime;
  
  int i, k = time.size();
  NumericMatrix d(k, 2);
  
  NumericVector d1(k), d2(k);
  d1 = a*pevent(s, piecewiseSurvivalTime, lambda1, gamma1)[0];
  d2 = a*pevent(s, piecewiseSurvivalTime, lambda2, gamma2)[0];
  
  for (i=0; i<k; i++) {
    v[0] = t[i];
    d(i,0) = phi*(d1[i] + ad(v, u[i], w[i], accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, lambda1, gamma1)[0]);
    d(i,1) = (1-phi)*(d2[i] + ad(v, u[i], w[i], accrualTime, accrualIntensity,
              piecewiseSurvivalTime, lambda2, gamma2)[0]);
  }
  
  return d;
}


//' @title Number of subjects having an event and log-rank statistic
//' for a hypothesized hazard ratio at a given calendar time
//'
//' @description Obtains the number of subjects having an event in each
//' treatment group by stratum, the mean and variance of weighted log-rank
//' score statistic for a hypothesized hazard ratio at a given calendar time.
//'
//' @param time The calendar time at which to calculate the number
//'  of events and the mean and variance of log-rank test score statistic.
//' @inheritParams param_hazardRatioH0
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @inheritParams param_numSubintervals
//' @param predictEventOnly Whether to predict the number of events only.
//'  Defaults to 0 for obtaining log-rank test score statistic mean
//'  and variance.
//'
//' @return A data frame of the number of subjects enrolled, the number of
//' subjects having an event overall and in each group, the number of
//' subjects who drop out overall and in each group, the mean and
//' variance of weighted log-rank score statistic at the specified
//' calendar time by stratum.
//'
//' @keywords internal
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' lrstat1(time = 22, hazardRatioH0 = 1,
//'         allocationRatioPlanned = 1,
//'         accrualTime = seq(0, 9),
//'         accrualIntensity = c(26/9*seq(1, 9), 26),
//'         piecewiseSurvivalTime = c(0, 6),
//'         stratumFraction = c(0.2, 0.8),
//'         lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'         lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'         gamma1 = -log(1-0.05)/12,
//'         gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 22,
//'         followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
DataFrame lrstat1(const double time = NA_REAL,
                  const double hazardRatioH0 = 1,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = NA_REAL,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const NumericVector& stratumFraction = 1,
                  const NumericVector& lambda1 = NA_REAL,
                  const NumericVector& lambda2 = NA_REAL,
                  const NumericVector& gamma1 = 0,
                  const NumericVector& gamma2 = 0,
                  const double accrualDuration = NA_REAL,
                  const double followupTime = NA_REAL,
                  const bool fixedFollowup = 0,
                  const double rho1 = 0,
                  const double rho2 = 0,
                  const int numSubintervals = 300,
                  const int predictEventOnly = 0) {
  
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi), gamma1x(nsi), gamma2x(nsi);
  
  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nsi);
  } else if (lambda1.size() == nintervals) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }
  
  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nintervals) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }
  
  if (gamma1.size() == 1) {
    gamma1x = rep(gamma1, nsi);
  } else if (gamma1.size() == nintervals) {
    gamma1x = rep(gamma1, nstrata);
  } else if (gamma1.size() == nsi) {
    gamma1x = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }
  
  if (gamma2.size() == 1) {
    gamma2x = rep(gamma2, nsi);
  } else if (gamma2.size() == nintervals) {
    gamma2x = rep(gamma2, nstrata);
  } else if (gamma2.size() == nsi) {
    gamma2x = gamma2;
  } else {
    stop("Invalid length for gamma2");
  }
  
  
  double minFollowupTime = followupTime;
  double maxFollowupTime;
  
  // obtain the follow-up time for the first enrolled subject
  if (fixedFollowup) {
    maxFollowupTime = minFollowupTime;
  } else {
    maxFollowupTime = accrualDuration + minFollowupTime;
  }
  
  IntegerVector l1 = Range(0, nintervals-1);
  IntegerVector q = Range(0, numSubintervals);
  NumericVector q1 = as<NumericVector>(q);
  Range q2 = Range(0, numSubintervals-1), c0 = Range(0,0), c1 = Range(1,1);
  
  double s = std::min(time, accrualDuration + minFollowupTime);
  NumericVector s1 = NumericVector::create(s);
  double a = accrual(s1, accrualTime, accrualIntensity, accrualDuration)[0];
  
  int h, i;
  double frac, accrualDuration0, minFollowupTime0, maxFollowupTime0, inc;
  IntegerVector l(nintervals);
  NumericVector lam1(nintervals), lam2(nintervals);
  NumericVector gam1(nintervals), gam2(nintervals);
  NumericMatrix x(1,2), y(1,2);
  NumericVector nsubjects(nstrata);
  NumericMatrix nevents(nstrata, 2), ndropouts(nstrata, 2);
  NumericVector t(numSubintervals+1);
  NumericMatrix xatrisk(numSubintervals+1, 2);
  NumericMatrix xevent(numSubintervals+1, 2);
  NumericVector atrisk1(numSubintervals), atrisk1x(numSubintervals),
  atrisk2(numSubintervals), atriskt(numSubintervals),
  atrisktx(numSubintervals), event1(numSubintervals), event2(numSubintervals),
  eventt(numSubintervals), km(numSubintervals), w(numSubintervals);
  NumericVector uscore(nstrata), vscore(nstrata), iscore(nstrata);
  NumericVector nevents1(nstrata), nevents2(nstrata), neventst(nstrata);
  NumericVector ndropouts1(nstrata), ndropouts2(nstrata), ndropoutst(nstrata);
  IntegerVector stratum(nstrata);
  NumericVector times(nstrata);
  DataFrame df;
  
  
  for (h=0; h<nstrata; h++) {
    frac = stratumFraction[h];
    l = h*nintervals + l1;
    lam1 = lambda1x[l];
    lam2 = lambda2x[l];
    gam1 = gamma1x[l];
    gam2 = gamma2x[l];
    
    
    // number of events in the stratum at the specified calendar time
    x = nevent2(s1, allocationRatioPlanned, accrualTime,
                frac*accrualIntensity,
                piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
                accrualDuration, minFollowupTime, maxFollowupTime);
    
    y = nevent2(s1, allocationRatioPlanned, accrualTime,
                frac*accrualIntensity,
                piecewiseSurvivalTime, gam1, gam2, lam1, lam2,
                accrualDuration, minFollowupTime, maxFollowupTime);
    
    // obtain number of enrolled subjects and subjects having an event
    nsubjects[h] = frac*a;
    nevents(h, _) = x.row(0);
    ndropouts(h, _) = y.row(0);
    
    // approximate the mean and variance of weighted log-rank test
    // score statistic
    if (predictEventOnly != 1) {
      
      // modify the study design at the calendar time of interest
      accrualDuration0 = std::min(s, accrualDuration);
      minFollowupTime0 = std::max(s - accrualDuration, 0.0);
      maxFollowupTime0 = std::min(s, maxFollowupTime);
      
      // partition the follow-up period into small sub-intervals
      inc = maxFollowupTime0/numSubintervals;
      t = q1*inc;
      
      // obtain number of subjects at risk and the number of subjects having
      // an event at each analysis time point
      xatrisk = natrisk(t, allocationRatioPlanned,
                        accrualTime, frac*accrualIntensity,
                        piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
                        accrualDuration0, minFollowupTime0, maxFollowupTime0);
      
      xevent = nevent(t, allocationRatioPlanned,
                      accrualTime, frac*accrualIntensity,
                      piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
                      accrualDuration0, minFollowupTime0, maxFollowupTime0);
      
      // number of subjects at risk at start of each analysis time interval
      atrisk1 = xatrisk(q2, c0);
      atrisk2 = xatrisk(q2, c1);
      atrisk1x = hazardRatioH0*atrisk1; // adjust with hazard ratio under H0
      atriskt = atrisk1 + atrisk2;
      atrisktx = atrisk1x + atrisk2;
      
      // number of subjects having an event in each analysis time interval
      event1 = diff(xevent(_, 0));
      event2 = diff(xevent(_, 1));
      eventt = event1 + event2;
      
      // Kaplan-Meier estimates of survival probabilities at the start of
      // each analysis time interval
      km[0] = 1;
      for (i=1; i<numSubintervals; i++) {
        km[i] = km[i-1]*(1 - eventt[i-1]/atriskt[i-1]);
      }
      
      // vector of Fleming-Harrington weights
      w = pow(km,rho1)*pow(1-km,rho2);
      
      // mean of the weighted log-rank test score statistic
      uscore[h] = sum(w * (event1 - eventt*atrisk1x/atrisktx));
      
      // variance of the weighted log-rank test score statistic
      vscore[h] = sum(w*w * eventt*atrisk1x*atrisk2/pow(atrisktx,2));
      
      // information of the weighted log-rank test score statistic
      iscore[h] = sum(w * eventt*atrisk1x*atrisk2/pow(atrisktx,2));
    }
  }
  
  // number of subjects having an event in each treatment group and overall
  nevents1 = nevents(_, 0);
  nevents2 = nevents(_, 1);
  neventst = nevents1 + nevents2;
  
  ndropouts1 = ndropouts(_, 0);
  ndropouts2 = ndropouts(_, 1);
  ndropoutst = ndropouts1 + ndropouts2;
  
  // stratum and time
  for (h=0; h<nstrata; h++) {
    stratum[h] = h+1;
    times[h] = s;
  }
  
  
  // output the requested information
  if (predictEventOnly == 1) {
    df = DataFrame::create(_["stratum"] = stratum,
                           _["time"] = times,
                           _["subjects"] = nsubjects,
                           _["nevents"] = neventst,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2,
                           _["ndropouts"] = ndropoutst,
                           _["ndropouts1"] = ndropouts1,
                           _["ndropouts2"] = ndropouts2);
  } else {
    df = DataFrame::create(_["stratum"] = stratum,
                           _["time"] = times,
                           _["subjects"] = nsubjects,
                           _["nevents"] = neventst,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2,
                           _["ndropouts"] = ndropoutst,
                           _["ndropouts1"] = ndropouts1,
                           _["ndropouts2"] = ndropouts2,
                           _["uscore"] = uscore,
                           _["vscore"] = vscore,
                           _["iscore"] = iscore);
  }
  
  return df;
}


//' @title Number of subjects having an event and log-rank statistics
//' @description Obtains the number of subjects accrued, number of events and
//' number of dropouts in each group, mean and variance of weighted log-rank
//' score statistic, estimated hazard ratio from weighted Cox regression
//' and variance of log hazard ratio estimate at given calendar times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'  of events and the mean and variance of log-rank test score statistic.
//' @inheritParams param_hazardRatioH0
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @inheritParams param_numSubintervals
//' @param predictEventOnly Whether to predict the number of events only.
//' Defaults to 0 for obtaining log-rank score statistic mean and variance.
//' Set \code{predictEventOnly = 1} for predicting the number of events only.
//' Set \code{predictEventOnly = 2} for predicting the number of events,
//' calculating the mean and variance of log-rank score statistic, and
//' calculating the estimated hazard ratio and variance of log hazard ratio.
//'
//' @return A data frame of the number of subjects enrolled, the number of
//' subjects having an event overall and in each group, the number of subjects
//' who drop out overall and in each group, the mean and variance of weighted
//' log-rank score statistic, the estimated hazard ratio from weighted Cox
//' regression, and variance of the log hazard ratio estimate at the
//' specified calendar times.
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' lrstat(time = c(22, 40), allocationRatioPlanned = 1,
//'        accrualTime = seq(0, 9),
//'        accrualIntensity = c(26/9*seq(1, 9), 26),
//'        piecewiseSurvivalTime = c(0, 6),
//'        stratumFraction = c(0.2, 0.8),
//'        lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'        lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'        gamma1 = -log(1-0.05)/12,
//'        gamma2 = -log(1-0.05)/12,
//'        accrualDuration = 22,
//'        followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
DataFrame lrstat(const NumericVector& time = NA_REAL,
                 const double hazardRatioH0 = 1,
                 const double allocationRatioPlanned = 1,
                 const NumericVector& accrualTime = 0,
                 const NumericVector& accrualIntensity = NA_REAL,
                 const NumericVector& piecewiseSurvivalTime = 0,
                 const NumericVector& stratumFraction = 1,
                 const NumericVector& lambda1 = NA_REAL,
                 const NumericVector& lambda2 = NA_REAL,
                 const NumericVector& gamma1 = 0,
                 const NumericVector& gamma2 = 0,
                 const double accrualDuration = NA_REAL,
                 const double followupTime = NA_REAL,
                 const bool fixedFollowup = 0,
                 const double rho1 = 0,
                 const double rho2 = 0,
                 const int numSubintervals = 300,
                 const int predictEventOnly = 0) {
  
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi), gamma1x(nsi), gamma2x(nsi);
  
  if (is_true(any(time < 0))) {
    stop("time must be non-negative");
  }
  
  if (hazardRatioH0 <= 0) {
    stop("hazardRatioH0 must be positive");
  }
  
  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }
  
  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }
  
  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }
  
  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }
  
  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }
  
  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  
  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }
  
  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }
  
  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }
  
  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }
  
  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }
  
  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }
  
  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }
  
  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nsi);
  } else if (lambda1.size() == nintervals) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }
  
  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nintervals) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }
  
  if (gamma1.size() == 1) {
    gamma1x = rep(gamma1, nsi);
  } else if (gamma1.size() == nintervals) {
    gamma1x = rep(gamma1, nstrata);
  } else if (gamma1.size() == nsi) {
    gamma1x = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }
  
  if (gamma2.size() == 1) {
    gamma2x = rep(gamma2, nsi);
  } else if (gamma2.size() == nintervals) {
    gamma2x = rep(gamma2, nstrata);
  } else if (gamma2.size() == nsi) {
    gamma2x = gamma2;
  } else {
    stop("Invalid length for gamma2");
  }
  
  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  }
  
  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }
  
  if (R_isnancpp(followupTime)) {
    stop("followupTime must be provided");
  }
  
  if (fixedFollowup && followupTime <= 0) {
    stop("followupTime must be positive for fixed follow-up");
  }
  
  if (!fixedFollowup && followupTime < 0) {
    stop("followupTime must be non-negative for variable follow-up");
  }
  
  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }
  
  if (rho2 < 0) {
    stop("rho2 must be non-negative");
  }
  
  if (numSubintervals <= 0) {
    stop("numSubintervals must be positive");
  }
  
  
  int k = time.size();
  
  DataFrame df;
  
  NumericVector subjects(k), nevents(k), nevents1(k), nevents2(k);
  NumericVector ndropouts(k), ndropouts1(k), ndropouts2(k);
  NumericVector uscore(k), vscore(k), logRankZ(k);
  NumericVector logHR(k), HR(k), vlogHR(k), zlogHR(k);
  
  for (int j=0; j<k; j++) {
    df = lrstat1(time[j], hazardRatioH0, allocationRatioPlanned,
                 accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1x, lambda2x, gamma1x, gamma2x,
                 accrualDuration, followupTime, fixedFollowup,
                 rho1, rho2, numSubintervals, predictEventOnly);
    
    subjects[j] = sum(NumericVector(df[2]));
    nevents[j] = sum(NumericVector(df[3]));
    nevents1[j] = sum(NumericVector(df[4]));
    nevents2[j] = sum(NumericVector(df[5]));
    ndropouts[j] = sum(NumericVector(df[6]));
    ndropouts1[j] = sum(NumericVector(df[7]));
    ndropouts2[j] = sum(NumericVector(df[8]));
    
    if (predictEventOnly != 1) {
      uscore[j] = sum(NumericVector(df[9]));
      vscore[j] = sum(NumericVector(df[10]));
      logRankZ[j] = uscore[j]/sqrt(vscore[j]);
    }
  }
  
  
  // solve for weighted Cox regression estimator
  if (predictEventOnly == 2) {
    double time1 = 0;
    
    auto g = [&time1, allocationRatioPlanned, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1x, lambda2x, gamma1x, gamma2x,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 
              predictEventOnly](double beta)->double {
                double hazardRatio = exp(beta);
                DataFrame df = lrstat1(time1, hazardRatio, 
                                       allocationRatioPlanned,
                                       accrualTime, accrualIntensity,
                                       piecewiseSurvivalTime, stratumFraction,
                                       lambda1x, lambda2x, gamma1x, gamma2x,
                                       accrualDuration, followupTime, 
                                       fixedFollowup,
                                       rho1, rho2, numSubintervals, 
                                       predictEventOnly);
                
                return sum(NumericVector(df[9]));
              };
    
    
    for (int j=0; j<k; j++) {
      time1 = time[j];
      logHR[j] = brent(g, -4.6, 4.6, 0.00001);
      HR[j] = exp(logHR[j]);
      
      DataFrame df = lrstat1(time1, HR[j], allocationRatioPlanned,
                             accrualTime, accrualIntensity,
                             piecewiseSurvivalTime, stratumFraction,
                             lambda1x, lambda2x, gamma1x, gamma2x,
                             accrualDuration, followupTime, fixedFollowup,
                             rho1, rho2, numSubintervals, predictEventOnly);
      
      double vscore1 = sum(NumericVector(df[10]));
      double iscore1 = sum(NumericVector(df[11]));
      
      vlogHR[j] = vscore1/(iscore1*iscore1);
      zlogHR[j] = (logHR[j] - log(hazardRatioH0))/sqrt(vlogHR[j]);
    }
    
    df = DataFrame::create(_["time"] = time,
                           _["subjects"] = subjects,
                           _["nevents"] = nevents,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2,
                           _["ndropouts"] = ndropouts,
                           _["ndropouts1"] = ndropouts1,
                           _["ndropouts2"] = ndropouts2,
                           _["uscore"] = uscore,
                           _["vscore"] = vscore,
                           _["logRankZ"] = logRankZ,
                           _["hazardRatioH0"] = hazardRatioH0,
                           _["HR"] = HR,
                           _["vlogHR"] = vlogHR,
                           _["zlogHR"] = zlogHR);
  } else if (predictEventOnly == 1) {
    df = DataFrame::create(_["time"] = time,
                           _["subjects"] = subjects,
                           _["nevents"] = nevents,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2,
                           _["ndropouts"] = ndropouts,
                           _["ndropouts1"] = ndropouts1,
                           _["ndropouts2"] = ndropouts2);
  } else {
    df = DataFrame::create(_["time"] = time,
                           _["subjects"] = subjects,
                           _["nevents"] = nevents,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2,
                           _["ndropouts"] = ndropouts,
                           _["ndropouts1"] = ndropouts1,
                           _["ndropouts2"] = ndropouts2,
                           _["uscore"] = uscore,
                           _["vscore"] = vscore,
                           _["logRankZ"] = logRankZ,
                           _["hazardRatioH0"] = hazardRatioH0);
  }
  
  
  return df;
  
}


//' @title Calendar times for target number of events
//' @description Obtains the calendar times to reach the target number of
//' subjects having an event.
//'
//' @param nevents A vector of target number of events.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_numSubintervals
//'
//' @return A vector of calendar times expected to yield the target
//' number of events.
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' caltime(nevents = c(24, 80), allocationRatioPlanned = 1,
//'         accrualTime = seq(0, 9),
//'         accrualIntensity = c(26/9*seq(1, 9), 26),
//'         piecewiseSurvivalTime = c(0, 6),
//'         stratumFraction = c(0.2, 0.8),
//'         lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'         lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'         gamma1 = -log(1-0.05)/12,
//'         gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 22,
//'         followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
NumericVector caltime(const NumericVector& nevents = NA_REAL,
                      const double allocationRatioPlanned = 1,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& stratumFraction = 1,
                      const NumericVector& lambda1 = NA_REAL,
                      const NumericVector& lambda2 = NA_REAL,
                      const NumericVector& gamma1 = 0,
                      const NumericVector& gamma2 = 0,
                      const double accrualDuration = NA_REAL,
                      const double followupTime = NA_REAL,
                      const bool fixedFollowup = 0,
                      const int numSubintervals = 300) {
  
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  
  if (is_true(any(nevents <= 0))) {
    stop("nevents must be positive");
  }
  
  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }
  
  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }
  
  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }
  
  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }
  
  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }
  
  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  
  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }
  
  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }
  
  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }
  
  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }
  
  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }
  
  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }
  
  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }
  
  if (lambda1.size() != 1 && lambda1.size() != nintervals &&
      lambda1.size() != nsi) {
    stop("Invalid length for lambda1");
  }
  
  if (lambda2.size() != 1 && lambda2.size() != nintervals &&
      lambda2.size() != nsi) {
    stop("Invalid length for lambda2");
  }
  
  if (gamma1.size() != 1 && gamma1.size() != nintervals &&
      gamma1.size() != nsi) {
    stop("Invalid length for gamma1");
  }
  
  if (gamma2.size() != 1 && gamma2.size() != nintervals &&
      gamma2.size() != nsi) {
    stop("Invalid length for gamma2");
  }
  
  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  }
  
  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }
  
  if (R_isnancpp(followupTime)) {
    stop("followupTime must be provided");
  }
  
  if (fixedFollowup && followupTime <= 0) {
    stop("followupTime must be positive for fixed follow-up");
  }
  
  if (!fixedFollowup && followupTime < 0) {
    stop("followupTime must be non-negative for variable follow-up");
  }
  
  if (numSubintervals <= 0) {
    stop("numSubintervals must be positive");
  }
  
  
  double event;
  
  // Lambda function
  auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup, numSubintervals,
            &event](double t)->double {
              NumericVector t0 = NumericVector::create(t);
              DataFrame lr = lrstat(
                t0, 1, allocationRatioPlanned, accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                0, 0, numSubintervals, 1);
              return sum(NumericVector(lr[2])) - event;
            };
  
  int i, k = nevents.size();
  double studyTime = accrualDuration + followupTime;
  NumericVector time(k);
  
  event = max(nevents);
  if (f(studyTime) < 0.0) {
    stop("followupTime is too short to reach the target number of events");
  }
  
  for (i=0; i<k; i++) {
    // match the predicted number of events to the target
    event = std::max(nevents[i], 0.0);
    time[i] = brent(f, 0.0001, studyTime, 0.0001);
  }
  
  return time;
}



//' @title Error spending functions
//' @description Obtains the error spent at the given information fraction
//' for the specified error spending function.
//'
//' @param t Information fraction for the interim look.
//' @param error Total error to spend.
//' @param sf Spending function. One of the following: "sfOF" for
//'  O'Brien-Fleming type spending function, "sfP" for Pocock type spending
//'  function, "sfKD" for Kim & DeMets spending function, and "sfHSD" for
//'  Hwang, Shi & DeCani spending function. Defaults to "sfOF".
//' @param sfpar Parameter for the spending function. Corresponds to rho for
//'  "sfKD" and gamma for "sfHSD"
//'
//' @return A number for the error spent up to the interim look.
//'
//' @keywords internal
//'
//' @examples
//' errorSpent(t = 0.5, error = 0.025, sf = "sfOF")
//' errorSpent(t = 0.5, error = 0.025, sf = "sfHSD", sfpar = -4)
//'
//' @export
// [[Rcpp::export]]
double errorSpent(const double t, const double error,
                  const String sf = "sfOF",
                  const double sfpar = NA_REAL) {
  if (error <= 0 || error >= 1) {
    stop("error must be a number between 0 and 1");
  }
  if (t <= 0 || t > 1) {
    stop("t must be a number between 0 and 1");
  }
  
  double aval;
  if (sf == "sfP") {
    aval = error*log(1 + (exp(1) - 1)*t);
  } else if (sf == "sfOF") {
    aval = R::qnorm(1-error/2, 0, 1, 1, 0);
    aval = 2*(1 - R::pnorm(aval/sqrt(t), 0, 1, 1, 0));
  } else if (sf == "sfKD") {
    if (R_isnancpp(sfpar)) {
      stop("Parameter sfpar is missing for sfKD");
    } else if (sfpar <= 0) {
      stop ("sfpar must be positive for sfKD");
    } else {
      aval = error*pow(t, sfpar);
    }
  } else if (sf == "sfHSD") {
    if (R_isnancpp(sfpar)) {
      stop("Parameter sfpar is missing for sfHSD");
    } else if (sfpar == 0) {
      aval = error*t;
    } else {
      aval = error*(1 - exp(-sfpar*t))/(1 - exp(-sfpar));
    }
  } else {
    stop("Invalid spending function");
  }
  return aval;
}


// [[Rcpp::export]]
NumericVector getCriticalValues(
    const int kMax = NA_INTEGER,
    const NumericVector& informationRates = NA_REAL,
    const LogicalVector& efficacyStopping = NA_LOGICAL,
    const double alpha = 0.025,
    const String typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const NumericVector& userAlphaSpending = NA_REAL,
    const double hazardRatioH0 = 1,
    const double allocationRatioPlanned = 1,
    const NumericVector& accrualTime = 0,
    const NumericVector& accrualIntensity = 20,
    const NumericVector& piecewiseSurvivalTime = 0,
    const NumericVector& stratumFraction = 1,
    const NumericVector& lambda2 = 0.0533,
    const NumericVector& gamma1 = 0,
    const NumericVector& gamma2 = 0,
    const double accrualDuration = 11.6,
    const double followupTime = 18,
    const bool fixedFollowup = 0,
    const double rho1 = 0,
    const double rho2 = 0,
    const int numSubintervals = 300) {
  
  NumericVector criticalValues(kMax);
  
  // treatment hazard under H0
  NumericVector lambda1 = hazardRatioH0*lambda2;
  
  NumericVector u0(1);
  DataFrame lr;
  NumericVector e0(kMax), time(kMax);
  
  // obtain the total number of events at study end
  u0[0] = accrualDuration + followupTime;
  lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 1);
  
  // obtain the timing of interim analysis
  e0 = sum(NumericVector(lr[2]))*informationRates;
  time = caltime(e0, allocationRatioPlanned,
                 accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 accrualDuration, followupTime, fixedFollowup,
                 numSubintervals);
  
  // obtain the mean and variance of log-rank test score statistic at
  // each stage
  lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 0);
  
  NumericVector vscore = NumericVector(lr[9]);
  
  NumericVector theta(kMax); // critical values under H0, initialized to zero
  
  // information time
  NumericVector t = vscore / (vscore[kMax - 1]);
  
  String asf = typeAlphaSpending;
  double asfpar = parameterAlphaSpending;
  
  if (asf == "none") {
    for (int i=0; i<kMax-1; i++) {
      criticalValues[i] = 6.0;
    }
    criticalValues[kMax-1] = R::qnorm(1-alpha, 0, 1, 1, 0);
  } else if (asf == "OF" || asf == "P" || asf == "WT") {
    double Delta;
    if (asf == "OF") {
      Delta = 0;
    } else if (asf == "P") {
      Delta = 0.5;
    } else {
      Delta = asfpar;
    }
    
    auto f = [kMax, alpha, Delta, theta, vscore, t,
              efficacyStopping] (double aval)->double {
                NumericVector u(kMax), l(kMax);
                for (int i=0; i<kMax; i++) {
                  u[i] = aval*pow(t[i], Delta-0.5);
                  if (!efficacyStopping[i]) u[i] = 6.0;
                  l[i] = -6.0;
                }
                
                List probs = exitprob(u, l, theta, vscore);
                NumericVector pu = NumericVector(probs[0]);
                return sum(pu) - alpha;
              };
    
    double cwt = brent(f, 0, 10, 1e-6);
    for (int i=0; i<kMax; i++) {
      criticalValues[i] = cwt*pow(t[i], Delta-0.5);
      if (!efficacyStopping[i]) criticalValues[i] = 6.0;
    }
  } else if (asf == "sfOF" || asf == "sfP" || asf == "sfKD" ||
    asf == "sfHSD" || asf == "user") {
    
    // stage 1
    double cumAlphaSpent;
    if (asf == "user") {
      cumAlphaSpent = userAlphaSpending[0];
    } else {
      cumAlphaSpent = errorSpent(t[0], alpha, asf, asfpar);
    }
    
    if (!efficacyStopping[0]) {
      criticalValues[0] = 6.0;
    } else {
      criticalValues[0] = R::qnorm(1 - cumAlphaSpent, 0, 1, 1, 0);
    }
    
    
    // lambda expression for finding the critical value at stage k
    int k=0;
    auto f = [&k, &cumAlphaSpent, &criticalValues,
              theta, vscore](double aval)->double {
                NumericVector u(k+1), l(k+1);
                for (int i=0; i<k; i++) {
                  u[i] = criticalValues[i];
                  l[i] = -6.0;
                }
                u[k] = aval;
                l[k] = -6.0;
                
                IntegerVector idx = Range(0,k);
                List probs = exitprob(u, l, theta[idx], vscore[idx]);
                NumericVector cpu = cumsum(NumericVector(probs[0]));
                return cpu[k] - cumAlphaSpent;
              };
    
    // subsequent stages
    for (k=1; k<kMax; k++) {
      if (asf == "user") {
        cumAlphaSpent = userAlphaSpending[k];
      } else {
        cumAlphaSpent = errorSpent(t[k], alpha, asf, asfpar);
      }
      
      if (!efficacyStopping[k]) {
        criticalValues[k] = 6.0;
      } else {
        if (f(6) > 0) { // no alpha spent at current visit
          criticalValues[k] = 6.0;
        } else {
          criticalValues[k] = brent(f, 0, 6, 1e-6);
        }
      }
    }
  } else {
    stop("Invalid type of alpha spending");
  }
  
  return criticalValues;
  
}


// [[Rcpp::export]]
NumericVector getCumAlphaSpent(
    const int kMax = NA_INTEGER,
    const NumericVector& informationRates = NA_REAL,
    const NumericVector& criticalValues = NA_REAL,
    const double hazardRatioH0 = 1,
    const double allocationRatioPlanned = 1,
    const NumericVector& accrualTime = 0,
    const NumericVector& accrualIntensity = 20,
    const NumericVector& piecewiseSurvivalTime = 0,
    const NumericVector& stratumFraction = 1,
    const NumericVector& lambda2 = 0.0533,
    const NumericVector& gamma1 = 0,
    const NumericVector& gamma2 = 0,
    const double accrualDuration = 11.6,
    const double followupTime = 18,
    const bool fixedFollowup = 0,
    const double rho1 = 0,
    const double rho2 = 0,
    const int numSubintervals = 300) {
  
  // treatment hazard under H0
  NumericVector lambda1 = hazardRatioH0*lambda2;
  
  NumericVector u0(1);
  DataFrame lr;
  NumericVector e0(kMax), time(kMax);
  
  // obtain the total number of events at study end
  u0[0] = accrualDuration + followupTime;
  lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 1);
  
  // obtain the timing of interim analysis
  e0 = sum(NumericVector(lr[2]))*informationRates;
  time = caltime(e0, allocationRatioPlanned,
                 accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 accrualDuration, followupTime, fixedFollowup,
                 numSubintervals);
  
  // obtain the mean and variance of log-rank test score statistic at
  // each stage
  lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 0);
  
  NumericVector vscore = NumericVector(lr[9]);
  
  NumericVector theta(kMax); // critical values under H0, initialized to zero
  NumericVector l = rep(-6.0, kMax);
  
  List probs = exitprob(criticalValues, l, theta, vscore);
  NumericVector pu = NumericVector(probs[0]);
  return cumsum(pu);
}



//' @title Log-rank test power
//' @description Estimates the power, stopping probabilities, and expected
//' sample size in a two-group survival design.
//'
//' @inheritParams param_kMax
//' @inheritParams param_informationRates
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @param typeBetaSpending The type of beta spending. One of the following:
//'  "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock type
//'  spending function, "sfKD" for Kim & DeMets spending function,
//'  "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
//'  early futility stopping. Defaults to "none".
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_hazardRatioH0
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @inheritParams param_numSubintervals
//' @inheritParams param_estimateHazardRatio
//' @inheritParams param_typeOfComputation
//' 
//' @return A list of S3 class \code{lrpower} with 3 components:
//'
//' * \code{overallResults} containing the overall rejection probability,
//' overall significance level, maximum and expected number of events,
//' maximum and expected number of dropouts, total and expected number
//' of subjects, maximum and expected study duration, along with
//' input parameters including accrual duration, followup duration,
//' whether a fixed follow-up is used, parameters for the FH weights,
//' allocation ratio, number of stages, and hazard ratio under H0.
//'
//' * \code{byStageResults} containing information rates, efficacy
//' and futility boundaries on the Z-scale, probability for efficacy
//' and futility stopping at the stage, cumulative probability for
//' efficacy and futility stopping by the stage, cumulative alpha spent,
//' expected number of events, number of dropouts, number of subjects,
//' and expected study time, efficacy and futility boundaries on
//' the HR scale and on the p-value scale, information for weighted
//' log-rank test, hazard ratio from weighted Cox regression, and
//' where efficacy and futility stopping are allowed by stage.
//'
//' * \code{settings} containing input parameters such as
//' alpha and beta spending function and parameter values,
//' accrual time, accrual intensity, piecewise survival time, stratum
//' fraction, and hazard rates for survival and dropout by group.
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' lrpower(kMax = 2, informationRates = c(0.8, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         allocationRatioPlanned = 1, accrualTime = seq(0, 9),
//'         accrualIntensity = c(26/9*seq(1, 9), 26),
//'         piecewiseSurvivalTime = c(0, 6),
//'         stratumFraction = c(0.2, 0.8),
//'         lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'         lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'         gamma1 = -log(1-0.05)/12,
//'         gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'         followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
List lrpower(const int kMax = NA_INTEGER,
             const NumericVector& informationRates = NA_REAL,
             const LogicalVector& efficacyStopping = NA_LOGICAL,
             const LogicalVector& futilityStopping = NA_LOGICAL,
             const NumericVector& criticalValues = NA_REAL,
             const double alpha = 0.025,
             const String typeAlphaSpending = "sfOF",
             const double parameterAlphaSpending = NA_REAL,
             const NumericVector& userAlphaSpending = NA_REAL,
             const NumericVector& futilityBounds = NA_REAL,
             const String typeBetaSpending = "none",
             const double parameterBetaSpending = NA_REAL,
             const double hazardRatioH0 = 1,
             const double allocationRatioPlanned = 1,
             const NumericVector& accrualTime = 0,
             const NumericVector& accrualIntensity = 20,
             const NumericVector& piecewiseSurvivalTime = 0,
             const NumericVector& stratumFraction = 1,
             const NumericVector& lambda1 = 0.0309,
             const NumericVector& lambda2 = 0.0533,
             const NumericVector& gamma1 = 0,
             const NumericVector& gamma2 = 0,
             const double accrualDuration = 11.6,
             const double followupTime = 18,
             const bool fixedFollowup = 0,
             const double rho1 = 0,
             const double rho2 = 0,
             const int numSubintervals = 300,
             const bool estimateHazardRatio = 1,
             const String typeOfComputation = "direct") {
  
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  
  String asf = typeAlphaSpending;
  double asfpar = parameterAlphaSpending;
  
  String bsf = typeBetaSpending;
  double bsfpar = parameterBetaSpending;
  
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  
  
  if (R_isnancpp(kMax)) {
    stop("kMax must be provided");
  }
  
  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }
  
  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = as<NumericVector>(tem)/(kMax+0.0);
  }
  
  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }
  
  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    }
  } else {
    futilityStopping1 = rep(1, kMax);
  }
  
  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }
  
  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
  }
  
  if (is_true(any(is_na(criticalValues))) && !(asf=="OF" || asf=="P" ||
      asf=="WT" || asf=="sfOF" || asf=="sfP" ||
      asf=="sfKD" || asf=="sfHSD" || asf=="user" || asf=="none")) {
    stop("Invalid type for alpha spending");
  }
  
  if ((asf=="WT" || asf=="sfKD" || asf=="sfHSD") && R_isnancpp(asfpar)) {
    stop("Missing parameter for the alpha spending function");
  }
  
  if (asf=="sfKD" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }
  
  if (asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegnative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }
  
  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  }
  
  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }
    
    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      stop("futilityBounds and criticalValues must meet at final analysis");
    }
  }
  
  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfOF" || bsf=="sfP" ||
      bsf=="sfKD" || bsf=="sfHSD" || bsf=="none")) {
    stop("Invalid type for beta spending");
  }
  
  if ((bsf=="sfKD" || bsf=="sfHSD") && R_isnancpp(bsfpar)) {
    stop("Missing parameter for the beta spending function");
  }
  
  if (bsf=="sfKD" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }
  
  if (hazardRatioH0 <= 0) {
    stop("hazardRatioH0 must be positive");
  }
  
  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }
  
  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }
  
  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }
  
  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }
  
  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }
  
  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  
  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }
  
  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }
  
  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }
  
  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }
  
  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }
  
  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }
  
  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }
  
  if (lambda1.size() != 1 && lambda1.size() != nintervals &&
      lambda1.size() != nsi) {
    stop("Invalid length for lambda1");
  }
  
  if (lambda2.size() != 1 && lambda2.size() != nintervals &&
      lambda2.size() != nsi) {
    stop("Invalid length for lambda2");
  }
  
  if (gamma1.size() != 1 && gamma1.size() != nintervals &&
      gamma1.size() != nsi) {
    stop("Invalid length for gamma1");
  }
  
  if (gamma2.size() != 1 && gamma2.size() != nintervals &&
      gamma2.size() != nsi) {
    stop("Invalid length for gamma2");
  }
  
  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  }
  
  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }
  
  if (R_isnancpp(followupTime)) {
    stop("followupTime must be provided");
  }
  
  if (fixedFollowup && followupTime <= 0) {
    stop("followupTime must be positive for fixed follow-up");
  }
  
  if (!fixedFollowup && followupTime < 0) {
    stop("followupTime must be non-negative for variable follow-up");
  }
  
  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }
  
  if (rho2 < 0) {
    stop("rho2 must be non-negative");
  }
  
  if (numSubintervals <= 0) {
    stop("numSubintervals must be positive");
  }
  
  
  if (typeOfComputation != "direct" && typeOfComputation != "Schoenfeld") {
    stop("typeOfComputation must be direct or Schoenfeld");
  }
  
  if (typeOfComputation == "Schoenfeld" && (rho1 != 0 || rho2 != 0)) {
    stop("Schoenfeld method can only be used for conventional log-rank test");
  }
  
  double hazardRatio;
  if (typeOfComputation == "Schoenfeld") {
    NumericVector lambda1x = rep(lambda1, nsi/lambda1.size());
    NumericVector lambda2x = rep(lambda2, nsi/lambda2.size());
    NumericVector hrx = lambda1x / lambda2x;
    
    bool proportionalHazards = 1;
    for (int i=1; i<nsi; i++) {
      if (fabs(hrx[i] - hrx[0]) > 1e-8) {
        proportionalHazards = 0;
        break;
      }
    }
    
    if (!proportionalHazards) {
      stop("Schoenfeld method can only be used for proportional hazards");
    } else {
      hazardRatio = hrx[0];
    }
  }
  
  
  NumericVector cumAlphaSpent(kMax);
  if (typeOfComputation == "direct") {
    if (is_true(any(is_na(criticalValues)))) {
      criticalValues1 = getCriticalValues(
        kMax, informationRates1, efficacyStopping1,
        alpha, asf, asfpar, userAlphaSpending, hazardRatioH0,
        allocationRatioPlanned, accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda2, gamma1, gamma2, accrualDuration,
        followupTime, fixedFollowup,
        rho1, rho2, numSubintervals);
    }
    
    cumAlphaSpent = getCumAlphaSpent(
      kMax, informationRates1, criticalValues1, hazardRatioH0,
      allocationRatioPlanned, accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda2, gamma1, gamma2, accrualDuration,
      followupTime, fixedFollowup,
      rho1, rho2, numSubintervals);
    
  } else if (typeOfComputation == "Schoenfeld") {
    if (is_true(any(is_na(criticalValues)))) {
      criticalValues1 = getCriticalValues(
        kMax, informationRates1, efficacyStopping1,
        alpha, asf, asfpar, userAlphaSpending, 1,
        allocationRatioPlanned, accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda2, gamma1, gamma2, accrualDuration,
        followupTime, fixedFollowup,
        rho1, rho2, numSubintervals);
    }
    
    cumAlphaSpent = getCumAlphaSpent(
      kMax, informationRates1, criticalValues1, 1,
      allocationRatioPlanned, accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda2, gamma1, gamma2, accrualDuration,
      followupTime, fixedFollowup,
      rho1, rho2, numSubintervals);
  }
  
  
  bool missingFutilityBounds = is_true(any(is_na(futilityBounds)));
  
  if (kMax > 1) {
    if (missingFutilityBounds && bsf=="none") {
      futilityBounds1 = rep(-6.0, kMax);
      futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    } else if (!missingFutilityBounds && futilityBounds.size() == kMax-1) {
      futilityBounds1.push_back(criticalValues1[kMax-1]);
    } else if (!missingFutilityBounds && futilityBounds.size() < kMax-1) {
      stop("Insufficient length of futilityBounds");
    }
  } else {
    if (missingFutilityBounds) {
      futilityBounds1 = criticalValues1[kMax-1];
    }
  }
  
  
  NumericVector u0(1);
  DataFrame lr;
  NumericVector e0(kMax), time(kMax);
  NumericVector HR(kMax), vlogHR(kMax), hru(kMax), hrl(kMax);
  
  // obtain the total number of events at study end
  u0[0] = accrualDuration + followupTime;
  lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 1);
  
  
  // obtain the timing of interim analysis
  e0 = sum(NumericVector(lr[2]))*informationRates1;
  time = caltime(e0, allocationRatioPlanned,
                 accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 accrualDuration, followupTime, fixedFollowup,
                 numSubintervals);

  // obtain the mean and variance of log-rank test score statistic at
  // each stage
  NumericVector theta(kMax), vscore(kMax);
  
  if (typeOfComputation == "Schoenfeld") {
    theta = rep(-log(hazardRatio/hazardRatioH0), kMax);
    double r1 = allocationRatioPlanned/(allocationRatioPlanned+1);
    vscore = r1*(1-r1)*e0;
    
    if (estimateHazardRatio) {
      HR = rep(hazardRatio, kMax);
      vlogHR = 1/vscore;
    }
    
    lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                rho1, rho2, numSubintervals, 1);
  } else{
    if (estimateHazardRatio) {
      lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  rho1, rho2, numSubintervals, 2);
    } else {
      lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  rho1, rho2, numSubintervals, 0);
    }
    
    if (estimateHazardRatio) {
      HR = NumericVector(lr[12]);
      vlogHR = NumericVector(lr[13]);
    }
    
    NumericVector uscore = NumericVector(lr[8]);
    vscore = NumericVector(lr[9]);
    theta = -uscore/vscore;
  }
  
  NumericVector nsubjects = NumericVector(lr[1]);
  NumericVector nevents = NumericVector(lr[2]);
  NumericVector ndropouts = NumericVector(lr[5]);

  
  // information time
  NumericVector t = vscore / (vscore[kMax - 1]);
  
  // compute the stagewise exit probabilities for efficacy and futility
  
  List probs;
  if (!missingFutilityBounds || bsf=="none" || kMax==1) {
    probs = exitprob(criticalValues1, futilityBounds1, theta, vscore);
  } else {
    auto f = [kMax, criticalValues1, futilityStopping1, &futilityBounds1,
              bsf, bsfpar, theta, vscore, t](double beta)->double {
                // initialize futilityBounds to be updated
                futilityBounds1 = NumericVector(kMax);
                double epsilon;
                
                // first stage
                int k = 0;
                double cumBetaSpent = errorSpent(t[0], beta, bsf, bsfpar);
                if (!futilityStopping1[0]) {
                  futilityBounds1[0] = -6.0;
                } else {
                  epsilon = R::pnorm(criticalValues1[0] -
                    theta[0]*sqrt(vscore[0]), 0, 1, 1, 0) - cumBetaSpent;
                  if (epsilon < 0) return -1.0; // to decrease beta
                  futilityBounds1[0] = R::qnorm(cumBetaSpent, 0, 1, 1, 0) +
                    theta[0]*sqrt(vscore[0]);
                }
                
                // lambda expression for finding the futility bound at stage k
                auto g = [&k, &cumBetaSpent, criticalValues1, &futilityBounds1,
                          theta, vscore](double aval)->double {
                            NumericVector u(k+1);
                            NumericVector l(k+1);
                            for (int i=0; i<k; i++) {
                              u[i] = criticalValues1[i];
                              l[i] = futilityBounds1[i];
                            }
                            u[k] = 6.0;
                            l[k] = aval;
                            
                            IntegerVector idx = Range(0,k);
                            List probs = exitprob(u, l, theta[idx], 
                                                  vscore[idx]);
                            NumericVector cpl = cumsum(NumericVector(probs[1]));
                            return cpl[k] - cumBetaSpent;
                          };
                
                for (k=1; k<kMax; k++) {
                  cumBetaSpent = errorSpent(t[k], beta, bsf, bsfpar);
                  
                  if (!futilityStopping1[k]) {
                    futilityBounds1[k] = -6.0;
                  } else {
                    epsilon = g(criticalValues1[k]);
                    
                    if (g(-6.0) > 0) { // no beta spent at current visit
                      futilityBounds1[k] = -6.0;
                    } else if (epsilon > 0) {
                      futilityBounds1[k] = brent(g, -6.0, criticalValues1[k], 
                                                 1e-6);
                    } else if (k < kMax-1) {
                      return -1.0;
                    }
                  }
                }
                
                return epsilon;
              };
    
    double v1 = f(0.0001), v2 = f(1-alpha);
    
    if (v1 == -1.0 || (v1 < 0 && futilityBounds1[kMax-1] == 0)) {
      stop("Power must be less than 0.9999 to use beta spending");
    } else if (v2 > 0) {
      stop("Power must be greater than alpha to use beta spending");
    } else {
      brent(f, 0.0001, 1-alpha, 1e-6);
      futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    }
    
    probs = exitprob(criticalValues1, futilityBounds1, theta, vscore);
  }
  
  
  
  NumericVector efficacyP(kMax);
  NumericVector futilityP(kMax);
  for (int i=0; i<kMax; i++) {
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
    futilityP[i] = 1 - R::pnorm(futilityBounds1[i], 0, 1, 1, 0);
  }
  
  // stagewise total exit probabilities
  NumericVector pu(kMax), pl(kMax), ptotal(kMax);
  pu = NumericVector(probs[0]);
  pl = NumericVector(probs[1]);
  ptotal = pu + pl;
  
  double overallReject = sum(pu);
  double expectedNumberOfEvents = sum(ptotal*nevents);
  double expectedNumberOfDropouts = sum(ptotal*ndropouts);
  double expectedNumberOfSubjects = sum(ptotal*nsubjects);
  double expectedStudyDuration = sum(ptotal*time);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);
  
  if (estimateHazardRatio) {
    hru = hazardRatioH0*exp(-criticalValues1*sqrt(vlogHR));
    hrl = hazardRatioH0*exp(-futilityBounds1*sqrt(vlogHR));
  }
  
  IntegerVector stageNumber = seq_len(kMax);
  
  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      efficacyStopping1[i] = 0;
    }
    
    if (futilityBounds1[i] == -6) {
      futilityStopping1[i] = 0;
    }
  }
  
  DataFrame byStageResults;
  
  if (estimateHazardRatio) {
    byStageResults = DataFrame::create(
      _["informationRates"] = informationRates1,
      _["efficacyBounds"] = criticalValues1,
      _["futilityBounds"] = futilityBounds1,
      _["rejectPerStage"] = pu,
      _["futilityPerStage"] = pl,
      _["cumulativeRejection"] = cpu,
      _["cumulativeFutility"] = cpl,
      _["cumulativeAlphaSpent"] = cumAlphaSpent,
      _["numberOfEvents"] = nevents,
      _["numberOfDropouts"] = ndropouts,
      _["numberOfSubjects"] = nsubjects,
      _["analysisTime"] = time,
      _["efficacyHR"] = hru,
      _["futilityHR"] = hrl,
      _["efficacyP"] = efficacyP,
      _["futilityP"] = futilityP,
      _["information"] = vscore,
      _["HR"] = HR,
      _["efficacyStopping"] = efficacyStopping1,
      _["futilityStopping"] = futilityStopping1);
  } else {
    byStageResults = DataFrame::create(
      _["informationRates"] = informationRates1,
      _["efficacyBounds"] = criticalValues1,
      _["futilityBounds"] = futilityBounds1,
      _["rejectPerStage"] = pu,
      _["futilityPerStage"] = pl,
      _["cumulativeRejection"] = cpu,
      _["cumulativeFutility"] = cpl,
      _["cumulativeAlphaSpent"] = cumAlphaSpent,
      _["numberOfEvents"] = nevents,
      _["numberOfDropouts"] = ndropouts,
      _["numberOfSubjects"] = nsubjects,
      _["analysisTime"] = time,
      _["efficacyP"] = efficacyP,
      _["futilityP"] = futilityP,
      _["information"] = vscore,
      _["efficacyStopping"] = efficacyStopping1,
      _["futilityStopping"] = futilityStopping1);
  }
  
  
  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfDropouts"] = (ndropouts[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfDropouts"] = expectedNumberOfDropouts,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["rho1"] = rho1,
    _["rho2"] = rho2,
    _["allocationRatioPlanned"] = allocationRatioPlanned,
    _["kMax"] = kMax,
    _["hazardRatioH0"] = hazardRatioH0,
    _["estimateHazardRatio"] = estimateHazardRatio,
    _["typeOfComputation"] = typeOfComputation);
  
  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["typeBetaSpending"] = typeBetaSpending,
    _["parameterBetaSpending"] = parameterBetaSpending,
    _["accrualTime"] = accrualTime,
    _["accrualIntensity"] = accrualIntensity,
    _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
    _["stratumFraction"] = stratumFraction,
    _["lambda1"] = lambda1,
    _["lambda2"] = lambda2,
    _["gamma1"] = gamma1,
    _["gamma2"] = gamma2);
  
  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);
  
  result.attr("class") = "lrpower";
  
  return result;
}


//' @title Log-rank test sample size
//' @description Obtains the needed accrual duration given power and
//' follow-up time, the needed follow-up time given power and
//' accrual duration, or the needed absolute accrual rates given
//' power, accrual duration, follow-up duration, and relative accrual
//' rates in a two-group survival design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @inheritParams param_informationRates
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
//' @inheritParams param_hazardRatioH0
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @inheritParams param_numSubintervals
//' @inheritParams param_estimateHazardRatio
//' @inheritParams param_typeOfComputation
//' @param interval The interval to search for the solution of
//' accrualDuration or followupDuration. Defaults to \code{c(0.001, 240)}.
//' Adjustment may be needed for non-monotone relationship with study power.
//'
//' @return A list of S3 class \code{lrpower}.
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' # Example 1: Obtains accrual duration given power and follow-up duration
//'
//' lrsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              accrualTime = seq(0, 9),
//'              accrualIntensity = c(26/9*seq(1, 9), 26),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12,
//'              accrualDuration = NA,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//'
//' # Example 2: Obtains follow-up duration given power and accrual duration
//'
//' lrsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              accrualTime = seq(0, 9),
//'              accrualIntensity = c(26/9*seq(1, 9), 26),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12,
//'              accrualDuration = 22,
//'              followupTime = NA, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains absolute accrual intensity given power,
//' # accrual duration, follow-up duration, and relative accrual intensity
//'
//' lrsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              accrualTime = seq(0, 9),
//'              accrualIntensity = c(26/9*seq(1, 9), 26),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12,
//'              accrualDuration = 22,
//'              followupTime = 18, fixedFollowup = FALSE)
//'              
//'
//' # Example 4: Non-inferiority trial with fixed follow-up and 
//' # superiority alternative
//' 
//' lr <- lrsamplesize(beta = 0.1, 
//'                    kMax = 3, 
//'                    alpha = 0.025,
//'                    typeAlphaSpending = "sfOF",
//'                    hazardRatioH0 = 1.1,
//'                    accrualTime = c(0, 6),
//'                    accrualIntensity = c(1000, 1500),
//'                    lambda1 = log(2)/48*0.95,
//'                    lambda2 = log(2)/48,
//'                    gamma1 = -log(1-0.08)/12,
//'                    gamma2 = -log(1-0.08)/12,
//'                    accrualDuration = NA,
//'                    followupTime = 18,
//'                    fixedFollowup = 1)
//'                    
//' lr
//'
//' @export
// [[Rcpp::export]]
List lrsamplesize(const double beta = 0.2,
                  const int kMax = NA_INTEGER,
                  const NumericVector& informationRates = NA_REAL,
                  const LogicalVector& efficacyStopping = NA_LOGICAL,
                  const LogicalVector& futilityStopping = NA_LOGICAL,
                  const NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.025,
                  const String typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const NumericVector& futilityBounds = NA_REAL,
                  const String typeBetaSpending = "none",
                  const double parameterBetaSpending = NA_REAL,
                  const NumericVector& userBetaSpending = NA_REAL,
                  const double hazardRatioH0 = 1,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = 20,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const NumericVector& stratumFraction = 1,
                  const NumericVector& lambda1 = 0.0309,
                  const NumericVector& lambda2 = 0.0533,
                  const NumericVector& gamma1 = 0,
                  const NumericVector& gamma2 = 0,
                  double accrualDuration = NA_REAL,
                  double followupTime = 18,
                  const bool fixedFollowup = 0,
                  const double rho1 = 0,
                  const double rho2 = 0,
                  const int numSubintervals = 300,
                  const bool estimateHazardRatio = 1,
                  const String typeOfComputation = "direct",
                  const NumericVector& interval =
                    NumericVector::create(0.001, 240)) {
  
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  
  String asf = typeAlphaSpending;
  double asfpar = parameterAlphaSpending;
  
  String bsf = typeBetaSpending;
  double bsfpar = parameterBetaSpending;
  
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  
  
  if (R_isnancpp(beta)) {
    stop("beta must be provided");
  }
  
  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }
  
  if (R_isnancpp(kMax)) {
    stop("kMax must be provided");
  }
  
  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }
  
  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = as<NumericVector>(tem)/(kMax+0.0);
  }
  
  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }
  
  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    }
  } else {
    futilityStopping1 = rep(1, kMax);
  }
  
  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }
  
  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
  }
  
  if (is_true(any(is_na(criticalValues))) && !(asf=="OF" || asf=="P" ||
      asf=="WT" || asf=="sfOF" || asf=="sfP" ||
      asf=="sfKD" || asf=="sfHSD" || asf=="user" || asf=="none")) {
    stop("Invalid type for alpha spending");
  }
  
  if ((asf=="WT" || asf=="sfKD" || asf=="sfHSD") && R_isnancpp(asfpar)) {
    stop("Missing parameter for the alpha spending function");
  }
  
  if (asf=="sfKD" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }
  
  if (asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegnative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }
  
  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  }
  
  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }
    
    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      stop("futilityBounds and criticalValues must meet at final analysis");
    }
  }
  
  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfOF" || bsf=="sfP" ||
      bsf=="sfKD" || bsf=="sfHSD" || bsf=="user" || bsf=="none")) {
    stop("Invalid type for beta spending");
  }
  
  if ((bsf=="sfKD" || bsf=="sfHSD") && R_isnancpp(bsfpar)) {
    stop("Missing parameter for the beta spending function");
  }
  
  if (bsf=="sfKD" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }
  
  if (bsf=="user") {
    if (is_true(any(is_na(userBetaSpending)))) {
      stop("userBetaSpending must be specified");
    } else if (userBetaSpending.size() < kMax) {
      stop("Insufficient length of userBetaSpending");
    } else if (userBetaSpending[0] < 0) {
      stop("Elements of userBetaSpending must be nonnegnative");
    } else if (kMax > 1 && is_true(any(diff(userBetaSpending) < 0))) {
      stop("Elements of userBetaSpending must be nondecreasing");
    } else if (userBetaSpending[kMax-1] != beta) {
      stop("userBetaSpending must end with specified beta");
    }
  }
  
  if (hazardRatioH0 <= 0) {
    stop("hazardRatioH0 must be positive");
  }
  
  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }
  
  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }
  
  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }
  
  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }
  
  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }
  
  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  
  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }
  
  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }
  
  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }
  
  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }
  
  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }
  
  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }
  
  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }
  
  if (lambda1.size() != 1 && lambda1.size() != nintervals &&
      lambda1.size() != nsi) {
    stop("Invalid length for lambda1");
  }
  
  if (lambda2.size() != 1 && lambda2.size() != nintervals &&
      lambda2.size() != nsi) {
    stop("Invalid length for lambda2");
  }
  
  if (gamma1.size() != 1 && gamma1.size() != nintervals &&
      gamma1.size() != nsi) {
    stop("Invalid length for gamma1");
  }
  
  if (gamma2.size() != 1 && gamma2.size() != nintervals &&
      gamma2.size() != nsi) {
    stop("Invalid length for gamma2");
  }
  
  if (!R_isnancpp(accrualDuration)) {
    if (accrualDuration <= 0) {
      stop("accrualDuration must be positive");
    }
  }
  
  if (!R_isnancpp(followupTime)) {
    if (fixedFollowup && followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }
    
    if (!fixedFollowup && followupTime < 0) {
      stop("followupTime must be non-negative for variable follow-up");
    }
  }
  
  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }
  
  if (rho2 < 0) {
    stop("rho2 must be non-negative");
  }
  
  if (numSubintervals <= 0) {
    stop("numSubintervals must be positive");
  }
  
  
  if (typeOfComputation != "direct" && typeOfComputation != "Schoenfeld") {
    stop("typeOfComputation must be direct or Schoenfeld");
  }
  
  if (typeOfComputation == "Schoenfeld" && (rho1 != 0 || rho2 != 0)) {
    stop("Schoenfeld method can only be used for conventional log-rank test");
  }
  
  double hazardRatio = 1;
  if (typeOfComputation == "Schoenfeld") {
    NumericVector lambda1x = rep(lambda1, nsi/lambda1.size());
    NumericVector lambda2x = rep(lambda2, nsi/lambda2.size());
    NumericVector hrx = lambda1x / lambda2x;
    
    bool proportionalHazards = 1;
    for (int i=1; i<nsi; i++) {
      if (fabs(hrx[i] - hrx[0]) > 1e-8) {
        proportionalHazards = 0;
        break;
      }
    }
    
    if (!proportionalHazards) {
      stop("Schoenfeld method can only be used for proportional hazards");
    } else {
      hazardRatio = hrx[0];
    }
  }
  
  
  if (interval.size() != 2) {
    stop("interval must have 2 elements");
  }
  
  if (interval[0] < 0) {
    stop("lower limit of interval must be positive");
  }
  
  if (interval[0] >= interval[1]) {
    stop("upper limit must be greater than lower limit for interval");
  }
  
  
  bool missingCriticalValues = is_true(any(is_na(criticalValues)));
  bool missingFutilityBounds = is_true(any(is_na(futilityBounds)));
  
  String unknown;
  
  // search for the solution according to the input
  if (R_isnancpp(accrualDuration) && !R_isnancpp(followupTime)) {
    unknown = "accrualDuration";
  } else if (!R_isnancpp(accrualDuration) && R_isnancpp(followupTime)) {
    unknown = "followupTime";
  } else if (!R_isnancpp(accrualDuration) && !R_isnancpp(followupTime)) {
    unknown = "accrualIntensity";
  } else {
    stop("accrualDuration and followupTime cannot be both missing");
  }
  
  
  auto f = [beta, kMax, informationRates1,
            efficacyStopping1, futilityStopping1,
            &criticalValues1, alpha, asf, asfpar, userAlphaSpending,
            &futilityBounds1, bsf, bsfpar, userBetaSpending, hazardRatioH0,
            allocationRatioPlanned, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction, lambda1, lambda2,
            gamma1, gamma2, accrualDuration, followupTime,
            fixedFollowup, rho1, rho2, numSubintervals, 
            typeOfComputation, hazardRatio, 
            nstrata, unknown, missingCriticalValues,
            missingFutilityBounds](double aval)->double {
              
              double dur1=0, dur2=0;
              NumericVector accrualIntensity1 = clone(accrualIntensity);
              
              if (unknown == "accrualDuration") {
                dur1 = aval;
                dur2 = followupTime;
              } else if (unknown == "followupTime") {
                dur1 = accrualDuration;
                dur2 = aval;
              } else if (unknown == "accrualIntensity") {
                dur1 = accrualDuration;
                dur2 = followupTime;
                accrualIntensity1 = aval*accrualIntensity;
              }
              
              if (missingCriticalValues) {
                if (typeOfComputation == "direct") {
                  criticalValues1 = getCriticalValues(
                    kMax, informationRates1, efficacyStopping1,
                    alpha, asf, asfpar, userAlphaSpending, hazardRatioH0,
                    allocationRatioPlanned, accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction, 
                    lambda2, gamma1, gamma2,
                    dur1, dur2, fixedFollowup, rho1, rho2, numSubintervals);
                } else if (typeOfComputation == "Schoenfeld") {
                  criticalValues1 = getCriticalValues(
                    kMax, informationRates1, efficacyStopping1,
                    alpha, asf, asfpar, userAlphaSpending, 1,
                    allocationRatioPlanned, accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction, 
                    lambda2, gamma1, gamma2,
                    dur1, dur2, fixedFollowup, rho1, rho2, numSubintervals);
                }
              }
              
              
              if (kMax > 1) {
                if (missingFutilityBounds && bsf=="none") {
                  futilityBounds1 = rep(-6.0, kMax);
                  futilityBounds1[kMax-1] = criticalValues1[kMax-1];
                } else if (!missingFutilityBounds && futilityBounds1.size() == 
                  kMax-1) {
                  futilityBounds1.push_back(criticalValues1[kMax-1]);
                } else if (!missingFutilityBounds && futilityBounds1.size() < 
                  kMax-1) {
                  stop("Insufficient length of futilityBounds");
                }
              } else {
                if (missingFutilityBounds) {
                  futilityBounds1 = criticalValues1[kMax-1];
                }
              }
              
              NumericVector u0(1);
              DataFrame lr;
              NumericVector e0(kMax), time(kMax);
              
              // obtain the total number of events at study end
              u0[0] = dur1 + dur2;
              lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
                          accrualTime, accrualIntensity1,
                          piecewiseSurvivalTime, stratumFraction,
                          lambda1, lambda2, gamma1, gamma2,
                          dur1, dur2, fixedFollowup,
                          rho1, rho2, numSubintervals, 1);
              
              
              // obtain the timing of interim analysis
              e0 = sum(NumericVector(lr[2]))*informationRates1;
              time = caltime(e0, allocationRatioPlanned, 
                             accrualTime, accrualIntensity1,
                             piecewiseSurvivalTime, stratumFraction,
                             lambda1, lambda2, gamma1, gamma2,
                             dur1, dur2, fixedFollowup);
              

              // obtain the mean and variance of log-rank test score statistic 
              // at each stage
              NumericVector theta(kMax), vscore(kMax);
              
              if (typeOfComputation == "Schoenfeld") {
                theta = rep(-log(hazardRatio/hazardRatioH0), kMax);
                double r1 = allocationRatioPlanned/(allocationRatioPlanned+1);
                vscore = r1*(1-r1)*e0;
              } else {
                lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
                            accrualTime, accrualIntensity1,
                            piecewiseSurvivalTime, stratumFraction,
                            lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup,
                            rho1, rho2, numSubintervals, 0);
                
                NumericVector uscore = NumericVector(lr[8]);
                vscore = NumericVector(lr[9]);
                theta = -uscore/vscore;
              }
              

              // information time
              NumericVector t = vscore / (vscore[kMax - 1]);
              
              // compute the stagewise exit probabilities for efficacy and 
              // futility
              
              if (!missingFutilityBounds || bsf=="none" || kMax==1) {
                List probs = exitprob(criticalValues1, futilityBounds1, theta, 
                                      vscore);
                NumericVector pu = NumericVector(probs[0]);
                double overallReject = sum(pu);
                return overallReject - (1-beta);
              } else {
                // initialize futilityBounds to be updated
                futilityBounds1 = NumericVector(kMax);
                double epsilon;
                
                // first stage
                int k = 0;
                double cumBetaSpent;
                if (bsf=="user") {
                  cumBetaSpent = userBetaSpending[0];
                } else {
                  cumBetaSpent = errorSpent(t[0], beta, bsf, bsfpar);
                }
                
                if (!futilityStopping1[0]) {
                  futilityBounds1[0] = -6.0;
                } else {
                  epsilon = R::pnorm(criticalValues1[0] -
                    theta[0]*sqrt(vscore[0]), 0, 1, 1, 0) - cumBetaSpent;
                  if (epsilon < 0) return -1.0;
                  futilityBounds1[0] = R::qnorm(cumBetaSpent, 0, 1, 1, 0) +
                    theta[0]*sqrt(vscore[0]);
                }
                
                
                // lambda expression for finding the futility bound at stage k
                auto g = [&k, &cumBetaSpent, criticalValues1, &futilityBounds1,
                          theta, vscore](double aval)->double {
                            NumericVector u(k+1);
                            NumericVector l(k+1);
                            for (int i=0; i<k; i++) {
                              u[i] = criticalValues1[i];
                              l[i] = futilityBounds1[i];
                            }
                            u[k] = 6.0;
                            l[k] = aval;
                            
                            IntegerVector idx = Range(0,k);
                            List probs = exitprob(u, l, theta[idx], 
                                                  vscore[idx]);
                            NumericVector cpl = cumsum(NumericVector(probs[1]));
                            return cpl[k] - cumBetaSpent;
                          };
                
                
                for (k=1; k<kMax; k++) {
                  if (bsf == "user") {
                    cumBetaSpent = userBetaSpending[k];
                  } else {
                    cumBetaSpent = errorSpent(t[k], beta, bsf, bsfpar);
                  }
                  
                  if (!futilityStopping1[k]) {
                    futilityBounds1[k] = -6.0;
                  } else {
                    epsilon = g(criticalValues1[k]);
                    
                    if (g(-6.0) > 0) { // no beta spent at the current visit
                      futilityBounds1[k] = -6.0;
                    } else if (epsilon > 0) {
                      futilityBounds1[k] = brent(g, -6.0, criticalValues1[k], 
                                                 1e-6);
                    } else if (k < kMax-1) {
                      return -1.0;
                    }
                    
                  }
                }
                
                return epsilon;
                
              }
            };
  
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  
  if (unknown == "accrualDuration") {
    accrualDuration = brent(f, interval[0], interval[1], 0.0001);
  } else if (unknown == "followupTime") {
    followupTime = brent(f, interval[0], interval[1], 0.0001);
  } else if (unknown == "accrualIntensity") {
    double aval = brent(f, interval[0], interval[1], 0.0001);
    accrualIntensity1 = aval*accrualIntensity;
  }
  
  futilityBounds1[kMax-1] = criticalValues1[kMax-1];
  
  // output the results
  List result = lrpower(kMax, informationRates1,
                        efficacyStopping1, futilityStopping1, criticalValues1,
                        alpha, typeAlphaSpending, parameterAlphaSpending,
                        userAlphaSpending, futilityBounds1,
                        typeBetaSpending, parameterBetaSpending, hazardRatioH0,
                        allocationRatioPlanned, accrualTime, accrualIntensity1,
                        piecewiseSurvivalTime, stratumFraction,
                        lambda1, lambda2, gamma1, gamma2,
                        accrualDuration, followupTime, fixedFollowup,
                        rho1, rho2, numSubintervals, estimateHazardRatio,
                        typeOfComputation);
  return result;
}
