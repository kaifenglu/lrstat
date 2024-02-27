#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;

//' @title Number of enrolled subjects
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


//' @title Accrual duration to enroll target number of subjects
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
NumericVector getAccrualDurationFromN(
    const NumericVector& nsubjects = NA_REAL,
    const NumericVector& accrualTime = 0,
    const NumericVector& accrualIntensity = NA_REAL) {

  int i, j, I = nsubjects.size(), J = accrualTime.size();
  NumericVector t(I), p(J);

  p[0] = 0;
  for (j=0; j<J-1; j++) {
    p[j+1] = p[j] + accrualIntensity[j]*(accrualTime[j+1] - accrualTime[j]);
  }

  IntegerVector m = findInterval2(nsubjects, p);

  for (i=0; i<I; i++) {
    j = m[i] - 1;
    t[i] = accrualTime[j] + (nsubjects[i] - p[j])/accrualIntensity[j];
  }

  return t;
}


//' @title Probability of being at risk
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
//' @keywords internal
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
//'   probability of having an event.
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
double hd(const int j = NA_INTEGER,
          const double t1 = NA_REAL,
          const double t2 = NA_REAL,
          const NumericVector& piecewiseSurvivalTime = 0,
          const NumericVector& lambda = NA_REAL,
          const NumericVector& gamma = 0) {

  int j1 = j-1;

  // lower bound of time interval j for piecewise exponential distribution
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
//' @keywords internal
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
  NumericVector a = accrual(u, accrualTime, accrualIntensity,
                            accrualDuration);

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
//' @keywords internal
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
  NumericVector a = accrual(u, accrualTime, accrualIntensity,
                            accrualDuration);

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
//' @keywords internal
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
  NumericVector a = accrual(u, accrualTime, accrualIntensity,
                            accrualDuration);

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
    d(i,1) = (1-phi)*(d2[i] + ad(v, u[i], w[i], accrualTime,
              accrualIntensity, piecewiseSurvivalTime, lambda2, gamma2)[0]);
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
//'   of events and the mean and variance of log-rank test score statistic.
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
//'   Defaults to 0 for obtaining log-rank test score statistic mean
//'   and variance.
//'
//' @return A data frame of the following variables if
//' \code{predictEventOnly = 1}:
//'
//' * \code{stratum}: The stratum number.
//'
//' * \code{time}: The analysis time since trial start.
//'
//' * \code{subjects}: The number of enrolled subjects.
//'
//' * \code{nevents}: The total number of events.
//'
//' * \code{nevents1}: The number of events in the active treatment group.
//'
//' * \code{nevents2}: The number of events in the control group.
//'
//' * \code{ndropouts}: The total number of dropouts.
//'
//' * \code{ndropouts1}: The number of dropouts in the active treatment
//'   group.
//'
//' * \code{ndropouts2}: The number of dropouts in the control group.
//'
//' * \code{nfmax}: The total number of subjects reaching maximum follow-up.
//'
//' * \code{nfmax1}: The number of subjects reaching maximum follow-up in
//'   the active treatment group.
//'
//' * \code{nfmax2}: The number of subjects reaching maximum follow-up in
//'   the control group.
//'
//' If \code{predictEventOnly = 0}, the following variables will also
//' be included:
//'
//' * \code{uscore}: The numerator of the weighted log-rank test statistic.
//'
//' * \code{vscore}: The variance of the weighted log-rank score statistic
//'   with weight squared.
//'
//' * \code{iscore}: The Fisher information of the weighted log-rank score
//'   statistic.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' lrstat1(time = 22, hazardRatioH0 = 1,
//'         allocationRatioPlanned = 1,
//'         accrualTime = seq(0, 8),
//'         accrualIntensity = 26/9*seq(1, 9),
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
                  const bool predictEventOnly = 0) {

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
  atrisktx(numSubintervals), event1(numSubintervals),
  event2(numSubintervals), eventt(numSubintervals), km(numSubintervals),
  w(numSubintervals);
  NumericVector uscore(nstrata), vscore(nstrata), iscore(nstrata);
  NumericVector nevents1(nstrata), nevents2(nstrata), neventst(nstrata);
  NumericVector ndropouts1(nstrata), ndropouts2(nstrata),
  ndropoutst(nstrata);
  IntegerVector stratum(nstrata);
  NumericVector times(nstrata);
  DataFrame df;

  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);
  NumericVector maxFU = NumericVector::create(maxFollowupTime);
  NumericVector s2 = NumericVector::create(s - maxFollowupTime);
  double a2 = accrual(s2, accrualTime, accrualIntensity, accrualDuration)[0];
  NumericVector nfmax1(nstrata), nfmax2(nstrata), nfmax(nstrata);

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

    // obtain number of subjects censored due to reaching the max follow-up
    double ncom = frac*a2;
    double p1 = patrisk(maxFU, piecewiseSurvivalTime, lam1, gam1)[0];
    double p2 = patrisk(maxFU, piecewiseSurvivalTime, lam2, gam2)[0];
    nfmax1[h] = phi*ncom*p1;
    nfmax2[h] = (1-phi)*ncom*p2;
    nfmax[h] = nfmax1[h] + nfmax2[h];

    // approximate the mean and variance of weighted log-rank test
    // score statistic
    if (!predictEventOnly) {

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
                        accrualDuration0, minFollowupTime0,
                        maxFollowupTime0);

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
  if (predictEventOnly) {
    df = DataFrame::create(_["stratum"] = stratum,
                           _["time"] = times,
                           _["subjects"] = nsubjects,
                           _["nevents"] = neventst,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2,
                           _["ndropouts"] = ndropoutst,
                           _["ndropouts1"] = ndropouts1,
                           _["ndropouts2"] = ndropouts2,
                           _["nfmax"] = nfmax,
                           _["nfmax1"] = nfmax1,
                           _["nfmax2"] = nfmax2);
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
                           _["nfmax"] = nfmax,
                           _["nfmax1"] = nfmax1,
                           _["nfmax2"] = nfmax2,
                           _["uscore"] = uscore,
                           _["vscore"] = vscore,
                           _["iscore"] = iscore);
  }

  return df;
}


//' @title Number of subjects having an event and log-rank statistics
//' @description Obtains the number of subjects accrued, number of events,
//' number of dropouts, and number of subjects reaching the maximum
//' follow-up in each group, mean and variance of weighted log-rank
//' score statistic, estimated hazard ratio from weighted Cox regression
//' and variance of log hazard ratio estimate at given calendar times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of events and the mean and variance of log-rank test score statistic.
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
//' @param predictTarget The target of prediction.
//'   Set \code{predictTarget = 1} to predict the number of events only.
//'   Set \code{predictTarget = 2} (default) to predict the number of events
//'   and log-rank score statistic mean and variance.
//'   Set \code{predictTarget = 3} to predict the number of events,
//'   log-rank score statistic mean and variance, and
//'   hazard ratio and variance of log hazard ratio.
//'
//' @return A data frame containing the following variables if
//' \code{predictTarget = 1}:
//'
//' * \code{time}: The analysis time since trial start.
//'
//' * \code{subjects}: The number of enrolled subjects.
//'
//' * \code{nevents}: The total number of events.
//'
//' * \code{nevents1}: The number of events in the active treatment group.
//'
//' * \code{nevents2}: The number of events in the control group.
//'
//' * \code{ndropouts}: The total number of dropouts.
//'
//' * \code{ndropouts1}: The number of dropouts in the active treatment
//'   group.
//'
//' * \code{ndropouts2}: The number of dropouts in the control group.
//'
//' * \code{nfmax}: The total number of subjects reaching maximum follow-up.
//'
//' * \code{nfmax1}: The number of subjects reaching maximum follow-up in
//'   the active treatment group.
//'
//' * \code{nfmax2}: The number of subjects reaching maximum follow-up in
//'   the control group.
//'
//' If \code{predictTarget = 2}, the following variables will also
//' be included:
//'
//' * \code{uscore}: The numerator of the log-rank test statistic.
//'
//' * \code{vscore}: The variance of the log-rank score test statistic.
//'
//' * \code{logRankZ}: The log-rank test statistic on the Z-scale.
//'
//' * \code{hazardRatioH0}: The hazard ratio under the null hypothesis.
//'
//' Furthermore, if \code{predictTarget = 3}, the following additional
//' variables will also be included:
//'
//' * \code{HR}: The average hazard ratio from weighted Cox regression.
//'
//' * \code{vlogHR}: The variance of log hazard ratio.
//'
//' * \code{zlogHR}: The Z-statistic for log hazard ratio.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' lrstat(time = c(22, 40), allocationRatioPlanned = 1,
//'        accrualTime = seq(0, 8),
//'        accrualIntensity = 26/9*seq(1, 9),
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
                 const int predictTarget = 2) {

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
  NumericVector nfmax(k), nfmax1(k), nfmax2(k);
  NumericVector uscore(k), vscore(k), logRankZ(k);
  NumericVector logHR(k), HR(k), vlogHR(k), zlogHR(k);

  if (predictTarget != 1 && predictTarget != 2 && predictTarget != 3) {
    stop("predictTarget must be equal to 1, 2, or 3");
  }

  bool predictEventOnly = predictTarget == 1;

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
    nfmax[j] = sum(NumericVector(df[9]));
    nfmax1[j] = sum(NumericVector(df[10]));
    nfmax2[j] = sum(NumericVector(df[11]));

    if (predictTarget > 1) {
      uscore[j] = sum(NumericVector(df[12]));
      vscore[j] = sum(NumericVector(df[13]));
      logRankZ[j] = uscore[j]/sqrt(vscore[j]);
    }
  }

  // solve for weighted Cox regression estimator
  if (predictTarget == 3) {
    double time1 = 0;

    auto g = [&time1, allocationRatioPlanned, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1x, lambda2x, gamma1x, gamma2x,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals,
              predictEventOnly](double beta)->double {
                double hazardRatio = exp(beta);
                DataFrame df = lrstat1(
                  time1, hazardRatio, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1x, lambda2x, gamma1x, gamma2x,
                  accrualDuration, followupTime, fixedFollowup,
                  rho1, rho2, numSubintervals, predictEventOnly);

                return sum(NumericVector(df[12]));
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

      double vscore1 = sum(NumericVector(df[13]));
      double iscore1 = sum(NumericVector(df[14]));

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
                           _["nfmax"] = nfmax,
                           _["nfmax1"] = nfmax1,
                           _["nfmax2"] = nfmax2,
                           _["uscore"] = uscore,
                           _["vscore"] = vscore,
                           _["logRankZ"] = logRankZ,
                           _["hazardRatioH0"] = hazardRatioH0,
                           _["HR"] = HR,
                           _["vlogHR"] = vlogHR,
                           _["zlogHR"] = zlogHR);
  } else if (predictTarget == 1) {
    df = DataFrame::create(_["time"] = time,
                           _["subjects"] = subjects,
                           _["nevents"] = nevents,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2,
                           _["ndropouts"] = ndropouts,
                           _["ndropouts1"] = ndropouts1,
                           _["ndropouts2"] = ndropouts2,
                           _["nfmax"] = nfmax,
                           _["nfmax1"] = nfmax1,
                           _["nfmax2"] = nfmax2);
  } else {
    df = DataFrame::create(_["time"] = time,
                           _["subjects"] = subjects,
                           _["nevents"] = nevents,
                           _["nevents1"] = nevents1,
                           _["nevents2"] = nevents2,
                           _["ndropouts"] = ndropouts,
                           _["ndropouts1"] = ndropouts1,
                           _["ndropouts2"] = ndropouts2,
                           _["nfmax"] = nfmax,
                           _["nfmax1"] = nfmax1,
                           _["nfmax2"] = nfmax2,
                           _["uscore"] = uscore,
                           _["vscore"] = vscore,
                           _["logRankZ"] = logRankZ,
                           _["hazardRatioH0"] = hazardRatioH0);
  }

  return df;
}


//' @title Kaplan-Meier estimate of milestone survival
//'
//' @description Obtains the Kaplan-Meier estimate of milestone survival
//' probability and associated variance estimate using the Greenwood formula
//' by treatment group and by stratum at given analysis time.
//'
//' @param time The calendar time for data cut.
//' @param milestone The milestone time at which to calculate the
//'   Kaplan-Meier estimate of survival probability.
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
//' @return A data frame containing the following variables:
//'
//' * \code{stratum}: The stratum.
//'
//' * \code{time}: The calendar time since trial start.
//'
//' * \code{subjects}: The enrolled number of subjects.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{surv1}: The milestone survival probability for the treatment
//'   group.
//'
//' * \code{surv2}: The milestone survival probability for the control group.
//'
//' * \code{vsurv1}: The variance for \code{surv1}.
//'
//' * \code{vsurv2}: The variance for \code{surv2}.
//'
//' * \code{survdiff}: The difference in milestone survival probabilities,
//'   i.e., \code{surv1 - surv2}.
//'
//' * \code{vsurvdiff}: The variance for \code{survdiff}.
//'
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' kmest1(time = 40,
//'        milestone = 18,
//'        allocationRatioPlanned = 1,
//'        accrualTime = seq(0, 8),
//'        accrualIntensity = 26/9*seq(1, 9),
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
DataFrame kmest1(const double time = NA_REAL,
                 const double milestone = NA_REAL,
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

  // modify the study design at the calendar time of interest
  double accrualDuration0 = std::min(s, accrualDuration);
  double minFollowupTime0 = std::max(s - accrualDuration, 0.0);
  double maxFollowupTime0 = std::min(s, maxFollowupTime);

  // partition the follow-up period into small sub-intervals
  double maxTime = std::min(milestone, maxFollowupTime0);
  double inc = maxTime/numSubintervals;
  NumericVector t = q1*inc;

  int h, i;
  double frac, km1, km2, vm1, vm2;
  IntegerVector l(nintervals);
  NumericVector lam1(nintervals), lam2(nintervals);
  NumericVector gam1(nintervals), gam2(nintervals);
  NumericVector nsubjects(nstrata);
  NumericMatrix xatrisk(numSubintervals+1, 2);
  NumericMatrix xevent(numSubintervals+1, 2);
  NumericVector atrisk1(numSubintervals), atrisk2(numSubintervals),
                event1(numSubintervals), event2(numSubintervals);
  NumericVector surv1(nstrata), surv2(nstrata), survdiff(nstrata);
  NumericVector vsurv1(nstrata), vsurv2(nstrata), vsurvdiff(nstrata);
  IntegerVector stratum(nstrata);
  NumericVector calTime(nstrata), mileTime(nstrata);
  DataFrame df;

  for (h=0; h<nstrata; h++) {
    frac = stratumFraction[h];
    l = h*nintervals + l1;
    lam1 = lambda1x[l];
    lam2 = lambda2x[l];
    gam1 = gamma1x[l];
    gam2 = gamma2x[l];

    // obtain number of enrolled subjects
    nsubjects[h] = frac*a;


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

    // number of subjects having an event in each analysis time interval
    event1 = diff(xevent(_, 0));
    event2 = diff(xevent(_, 1));

    // Kaplan-Meier estimates of survival probabilities
    km1 = 1;
    km2 = 1;
    for (i=0; i<numSubintervals; i++) {
      km1 = km1*(1 - event1[i]/atrisk1[i]);
      km2 = km2*(1 - event2[i]/atrisk2[i]);
    }

    vm1 = 0;
    vm2 = 0;
    for (i=0; i<numSubintervals; i++) {
      vm1 = vm1 + event1[i]/(atrisk1[i]*(atrisk1[i] - event1[i]));
      vm2 = vm2 + event2[i]/(atrisk2[i]*(atrisk2[i] - event2[i]));
    }

    surv1[h] = km1;
    surv2[h] = km2;
    vsurv1[h] = km1*km1*vm1;
    vsurv2[h] = km2*km2*vm2;
    survdiff[h] = surv1[h] - surv2[h];
    vsurvdiff[h] = vsurv1[h] + vsurv2[h];
  }

  // stratum and time
  for (h=0; h<nstrata; h++) {
    stratum[h] = h+1;
    calTime[h] = s;
    mileTime[h] = maxTime;
  }

  // output the requested information
  df = DataFrame::create(_["stratum"] = stratum,
                         _["time"] = calTime,
                         _["subjects"] = nsubjects,
                         _["milestone"] = maxTime,
                         _["surv1"] = surv1,
                         _["surv2"] = surv2,
                         _["vsurv1"] = vsurv1,
                         _["vsurv2"] = vsurv2,
                         _["survdiff"] = survdiff,
                         _["vsurvdiff"] = vsurvdiff);

  return df;
}


//' @title Stratified difference in milestone survival
//' @description Obtains the stratified Kaplan-Meier estimate of
//'   milestone survival probabilities and difference in milestone
//'   survival at given calendar times and milestone time.
//'
//' @param time A vector of calendar times at which to calculate the
//'   milestone survival.
//' @param milestone The milestone time at which to calculate the
//'   Kaplan-Meier estimate of survival probability.
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
//' @return A data frame containing the following variables:
//'
//' * \code{time}: The calendar time at which to calculate the milestone
//'   survival.
//'
//' * \code{subjects}: The enrolled number of subjects.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{surv1}: The milestone survival probability for the treatment
//'   group.
//'
//' * \code{surv2}: The milestone survival probability for the control group.
//'
//' * \code{vsurv1}: The variance for \code{surv1}.
//'
//' * \code{vsurv2}: The variance for \code{surv2}.
//'
//' * \code{survdiff}: The difference in milestone survival probabilities,
//'   i.e., \code{surv1 - surv2}.
//'
//' * \code{vsurvdiff}: The variance for \code{survdiff}.
//'
//' * \code{survdiffZ}: The Z-statistic value, i.e.,
//'   \code{survdiff/sqrt(vsurvdiff)}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' kmest(time = c(22, 40),
//'       milestone = 18,
//'       allocationRatioPlanned = 1,
//'       accrualTime = seq(0, 8),
//'       accrualIntensity = 26/9*seq(1, 9),
//'       piecewiseSurvivalTime = c(0, 6),
//'       stratumFraction = c(0.2, 0.8),
//'       lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'       lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'       gamma1 = -log(1-0.05)/12,
//'       gamma2 = -log(1-0.05)/12,
//'       accrualDuration = 22,
//'       followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
DataFrame kmest(const NumericVector& time = NA_REAL,
                const double milestone = NA_REAL,
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
  NumericVector lambda1x(nsi), lambda2x(nsi), gamma1x(nsi), gamma2x(nsi);

  if (is_true(any(time <= 0))) {
    stop("time must be positive");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
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

  if (numSubintervals <= 0) {
    stop("numSubintervals must be positive");
  }

  int k = time.size();
  DataFrame df;
  NumericVector calTime(k), mileTime(k), subjects(k),
                surv1(k), surv2(k), vsurv1(k), vsurv2(k),
                survdiff(k), vsurvdiff(k), survdiffZ(k);

  for (int j=0; j<k; j++) {
    df = kmest1(time[j], milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1x, lambda2x, gamma1x, gamma2x,
                accrualDuration, followupTime, fixedFollowup,
                numSubintervals);

    calTime[j] = max(NumericVector(df[1]));
    subjects[j] = sum(NumericVector(df[2]));
    mileTime[j] = max(NumericVector(df[3]));
    surv1[j] = sum(stratumFraction*NumericVector(df[4]));
    surv2[j] = sum(stratumFraction*NumericVector(df[5]));
    vsurv1[j] = sum(stratumFraction*stratumFraction*NumericVector(df[6]));
    vsurv2[j] = sum(stratumFraction*stratumFraction*NumericVector(df[7]));
    survdiff[j] = sum(stratumFraction*NumericVector(df[8]));
    vsurvdiff[j] = sum(stratumFraction*stratumFraction*NumericVector(df[9]));
    survdiffZ[j] = survdiff[j]/sqrt(vsurvdiff[j]);
  }

  df = DataFrame::create(_["time"] = calTime,
                         _["subjects"] = subjects,
                         _["milestone"] = mileTime,
                         _["surv1"] = surv1,
                         _["surv2"] = surv2,
                         _["vsurv1"] = vsurv1,
                         _["vsurv2"] = vsurv2,
                         _["survdiff"] = survdiff,
                         _["vsurvdiff"] = vsurvdiff,
                         _["survdiffZ"] = survdiffZ);

  return df;
}


//' @title Calendar times for target number of events
//' @description Obtains the calendar times needed to reach the target
//' number of subjects experiencing an event.
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
//'
//' @return A vector of calendar times expected to yield the target
//' number of events.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' caltime(nevents = c(24, 80), allocationRatioPlanned = 1,
//'         accrualTime = seq(0, 8),
//'         accrualIntensity = 26/9*seq(1, 9),
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
                      const bool fixedFollowup = 0) {

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


  double event;

  // Lambda function
  auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &event](double t)->double {
              NumericVector t0 = NumericVector::create(t);
              DataFrame lr = lrstat(
                t0, 1, allocationRatioPlanned, accrualTime,
                accrualIntensity, piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                0, 0, 1, 1);
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


//' @title Range of accrual duration for target number of events
//' @description Obtains a range of accrual duration to reach the
//' target number of events.
//'
//' @param nevents The target number of events.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param npoints The number of accrual duration time points.
//'   Defaults to 23.
//' @param interval The interval to search for the solution of
//'   accrualDuration. Defaults to \code{c(0.001, 240)}.
//'
//' @return A data frame of the following variables:
//'
//' * \code{nevents}: The target number of events.
//'
//' * \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//' * \code{accrualDuration}: The accrual duration.
//'
//' * \code{subjects}: The total number of subjects.
//'
//' * \code{studyDuration}: The study duration.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' getDurationFromNevents(
//'   nevents = 80, allocationRatioPlanned = 1,
//'   accrualTime = seq(0, 8),
//'   accrualIntensity = 26/9*seq(1, 9),
//'   piecewiseSurvivalTime = c(0, 6),
//'   stratumFraction = c(0.2, 0.8),
//'   lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'   gamma1 = -log(1-0.05)/12,
//'   gamma2 = -log(1-0.05)/12,
//'   fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
DataFrame getDurationFromNevents(
    const double nevents = NA_REAL,
    const double allocationRatioPlanned = 1,
    const NumericVector& accrualTime = 0,
    const NumericVector& accrualIntensity = NA_REAL,
    const NumericVector& piecewiseSurvivalTime = 0,
    const NumericVector& stratumFraction = 1,
    const NumericVector& lambda1 = NA_REAL,
    const NumericVector& lambda2 = NA_REAL,
    const NumericVector& gamma1 = 0,
    const NumericVector& gamma2 = 0,
    const double followupTime = 18,
    const bool fixedFollowup = 0,
    const int npoints = 23,
    const NumericVector& interval =
      NumericVector::create(0.001, 240)) {

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;

  if (R_isnancpp(nevents)) {
    stop("nevents must be provided");
  }

  if (nevents <= 0) {
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


  if (fixedFollowup && !R_isnancpp(followupTime) && followupTime <= 0) {
    stop("followupTime must be positive for fixed follow-up");
  }

  if (npoints < 2) {
    stop("npoints must be greater than or equal to 2");
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


  NumericVector t(2);

  // Lambda function
  if (!fixedFollowup) {
    auto fmin = [allocationRatioPlanned, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 fixedFollowup, nevents](double t)->double {
                   NumericVector t0(1);
                   t0[0] = t + 1000;
                   DataFrame lr = lrstat(
                     t0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, stratumFraction,
                     lambda1, lambda2, gamma1, gamma2,
                     t, 1000, fixedFollowup, 0, 0, 1, 1);
                   return sum(NumericVector(lr[2])) - nevents;
                 };

    t[0] = brent(fmin, interval[0], interval[1], 0.001);

    auto fmax = [allocationRatioPlanned, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 fixedFollowup, nevents](double t)->double {
                   NumericVector t0(1);
                   t0[0] = t;
                   DataFrame lr = lrstat(
                     t0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, stratumFraction,
                     lambda1, lambda2, gamma1, gamma2,
                     t, 0, fixedFollowup, 0, 0, 1, 1);
                   return sum(NumericVector(lr[2])) - nevents;
                 };

    t[1] = brent(fmax, t[0], interval[1], 0.001);
  } else {
    auto fmin = [allocationRatioPlanned, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 followupTime, fixedFollowup, nevents](double t)->double {
                   NumericVector t0(1);
                   t0[0] = t + followupTime;
                   DataFrame lr = lrstat(
                     t0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, stratumFraction,
                     lambda1, lambda2, gamma1, gamma2,
                     t, followupTime, fixedFollowup, 0, 0, 1, 1);
                   return sum(NumericVector(lr[2])) - nevents;
                 };

    t[0] = brent(fmin, interval[0], interval[1], 0.001);

    auto fmax = [allocationRatioPlanned, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 followupTime, fixedFollowup, nevents](double t)->double {
                   NumericVector t0(1);
                   t0[0] = t;
                   DataFrame lr = lrstat(
                     t0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, stratumFraction,
                     lambda1, lambda2, gamma1, gamma2,
                     t, followupTime, fixedFollowup, 0, 0, 1, 1);
                   return sum(NumericVector(lr[2])) - nevents;
                 };

    t[1] = brent(fmax, t[0], interval[1], 0.001);
  }

  NumericVector time(1), bigd(1);
  bigd[0] = nevents;

  NumericVector ta(npoints), n(npoints), ts(npoints);
  double dt = (t[1] - t[0])/(npoints - 1);

  for (int i=0; i<npoints; i++) {
    ta[i] = t[0] + i*dt;
    time = ta[i];

    if (!fixedFollowup) {
      if (i==0) {
        ts[i] = ta[i] + 1000;
      } else if (i == npoints - 1){
        ts[i] = ta[i];
      } else {
        time = caltime(bigd, allocationRatioPlanned,
                       accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, stratumFraction,
                       lambda1, lambda2, gamma1, gamma2,
                       ta[i], 1000, fixedFollowup);
        ts[i] = time[0];
      }
    } else {
      if (i==0) {
        ts[i] = ta[i] + followupTime;
      } else if (i == npoints - 1) {
        ts[i] = ta[i];
      } else {
        time = caltime(bigd, allocationRatioPlanned,
                       accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, stratumFraction,
                       lambda1, lambda2, gamma1, gamma2,
                       ta[i], followupTime, fixedFollowup);
        ts[i] = time[0];
      }
    }
  }

  n = accrual(ta, accrualTime, accrualIntensity, 1000);

  DataFrame df;

  if (!fixedFollowup) {
    df = DataFrame::create(_["nevents"] = rep(nevents, npoints),
                           _["fixedFollowup"] = rep(fixedFollowup, npoints),
                           _["accrualDuration"] = ta,
                           _["subjects"] = n,
                           _["studyDuration"] = ts);
  } else {
    df = DataFrame::create(_["nevents"] = rep(nevents, npoints),
                           _["fixedFollowup"] = rep(fixedFollowup, npoints),
                           _["followupTime"] = rep(followupTime, npoints),
                           _["accrualDuration"] = ta,
                           _["subjects"] = n,
                           _["studyDuration"] = ts);
  }

  return df;
}


//' @title Log-rank test power
//' @description Estimates the power, stopping probabilities, and expected
//' sample size in a two-group survival design.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates in terms of number
//'   of events for the conventional log-rank test and in terms of
//'   the actual information for weighted log-rank tests.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @param typeBetaSpending The type of beta spending. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock
//'   type spending function, "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
//'   early futility stopping. Defaults to "none".
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
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{lrpower} object with 4 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numberOfDropouts}: The total number of dropouts.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfDropouts}: The expected number of dropouts.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//'     - \code{rho1}: The first parameter of the Fleming-Harrington family
//'       of weighted log-rank test.
//'
//'     - \code{rho2}: The second parameter of the Fleming-Harrington family
//'       of weighted log-rank test.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{hazardRatioH0}: The hazard ratio under the null hypothesis.
//'
//'     - \code{typeOfComputation}: The type of computation,
//'       either "direct" for the direct approximation method,
//'       or "schoenfeld" for the Schoenfeld method.
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
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyHR}: The efficacy boundaries on the hazard ratio
//'       scale if \code{estimateHazardRatio}.
//'
//'     - \code{futilityHR}: The futility boundaries on the hazard ratio
//'       scale if \code{estimateHazardRatio}.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{futilityP}: The futility boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{HR}: The average hazard ratio.
//'
//'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'
//'     - \code{futilityStopping}: Whether to allow futility stopping.
//'
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{typeBetaSpending},
//'   \code{parameterBetaSpending}, \code{allocationRatioPlanned},
//'   \code{accrualTime}, \code{accuralIntensity},
//'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   \code{estimateHazardRatio}, and \code{spendingTime}.
//'
//' * \code{byTreatmentCounts}: A list containing the following counts by
//'   treatment group:
//'
//'     - \code{numberOfEvents1}: The number of events by stage for
//'       the treatment group.
//'
//'     - \code{numberOfDropouts1}: The number of dropouts by stage for
//'       the treatment group.
//'
//'     - \code{numberOfSubjects1}: The number of subjects by stage for
//'       the treatment group.
//'
//'     - \code{numberOfEvents2}: The number of events by stage for
//'       the control group.
//'
//'     - \code{numberOfDropouts2}: The number of dropouts by stage for
//'       the control group.
//'
//'     - \code{numberOfSubjects2}: The number of subjects by stage for
//'       the control group.
//'
//'     - \code{expectedNumberOfEvents1}: The expected number of events for
//'       the treatment group.
//'
//'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
//'       for the treatment group.
//'
//'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
//'       for the treatment group.
//'
//'     - \code{expectedNumberOfEvents2}: The expected number of events for
//'       control group.
//'
//'     - \code{expectedNumberOfDropouts2}: The expected number of dropouts
//'       for the control group.
//'
//'     - \code{expectedNumberOfSubjects2}: The expected number of subjects
//'       for the control group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survival, and 5% dropout by
//' # the end of 1 year.
//'
//' lrpower(kMax = 2, informationRates = c(0.8, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'         accrualIntensity = 26/9*seq(1, 9),
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
List lrpower(const int kMax = 1,
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
             const String typeOfComputation = "direct",
             const NumericVector& spendingTime = NA_REAL,
             const double studyDuration = NA_REAL) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double bsfpar = parameterBetaSpending;

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;


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
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
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
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
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

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
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

  if (fixedFollowup && R_isnancpp(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
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

  std::string su = typeOfComputation;
  std::for_each(su.begin(), su.end(), [](char & c) {
    c = std::tolower(c);
  });

  if (su != "direct" && su != "schoenfeld") {
    stop("typeOfComputation must be direct or Schoenfeld");
  }

  if (su == "schoenfeld" && (rho1 != 0 || rho2 != 0)) {
    stop("Schoenfeld method can only be used for ordinary log-rank test");
  }

  double hazardRatio = 1;
  if (su == "schoenfeld") {
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

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));
  alpha1 = cumAlphaSpent[kMax - 1];

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

  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || R_isnancpp(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }
  u0[0] = studyDuration1;

  // obtain the timing of interim analysis
  if (rho1 == 0 && rho2 == 0) { // conventional log-rank test
    lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                rho1, rho2, numSubintervals, 1);

    e0 = sum(NumericVector(lr[2]))*informationRates1;
    time = caltime(e0, allocationRatioPlanned,
                   accrualTime, accrualIntensity,
                   piecewiseSurvivalTime, stratumFraction,
                   lambda1, lambda2, gamma1, gamma2,
                   accrualDuration, followupTime, fixedFollowup);
  } else {
    lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                rho1, rho2, numSubintervals, 2);

    double maxInformation = sum(NumericVector(lr[12]));
    double information1;

    auto f = [hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals,
              &information1](double aval)->double {
                NumericVector u0(1, aval);
                DataFrame lr = lrstat(
                  u0, hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  rho1, rho2, numSubintervals, 2);
                return sum(NumericVector(lr[12])) - information1;
              };

    for (int i=0; i<kMax; i++) {
      information1 = maxInformation*informationRates1[i];
      time[i] = brent(f, 0.0001, studyDuration1, 0.0001);
    };
  }

  // obtain mean and variance of log-rank test score statistic at each stage
  NumericVector theta(kMax), vscore(kMax);

  double r1 = allocationRatioPlanned/(allocationRatioPlanned+1);

  if (su == "schoenfeld") {
    theta = rep(-log(hazardRatio/hazardRatioH0), kMax);

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
  } else {
    if (estimateHazardRatio) {
      lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  rho1, rho2, numSubintervals, 3);
    } else {
      lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  rho1, rho2, numSubintervals, 2);
    }

    if (estimateHazardRatio) {
      HR = NumericVector(lr[15]);
      vlogHR = NumericVector(lr[16]);
    }

    NumericVector uscore = NumericVector(lr[11]);
    vscore = NumericVector(lr[12]);

    theta = -uscore/vscore;
  }

  NumericVector nsubjects = NumericVector(lr[1]);
  NumericVector nsubjects1 = r1*nsubjects;
  NumericVector nsubjects2 = (1-r1)*nsubjects;
  NumericVector nevents = NumericVector(lr[2]);
  NumericVector nevents1 = NumericVector(lr[3]);
  NumericVector nevents2 = NumericVector(lr[4]);
  NumericVector ndropouts = NumericVector(lr[5]);
  NumericVector ndropouts1 = NumericVector(lr[6]);
  NumericVector ndropouts2 = NumericVector(lr[7]);

  // compute the stagewise exit probabilities for efficacy and futility
  if (!missingFutilityBounds || bsf=="none" || kMax==1) {
    probs = exitprobcpp(criticalValues1, futilityBounds1, theta, vscore);
  } else {
    NumericVector w(kMax, 1.0);
    List out = getPower(alpha1, kMax, criticalValues1, theta, vscore,
                        bsf, bsfpar, spendingTime1, futilityStopping1, w);
    futilityBounds1 = out[1];
    probs = out[2];
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
  double expectedNumberOfEvents1 = sum(ptotal*nevents1);
  double expectedNumberOfDropouts1 = sum(ptotal*ndropouts1);
  double expectedNumberOfSubjects1 = sum(ptotal*nsubjects1);
  double expectedNumberOfEvents2 = sum(ptotal*nevents2);
  double expectedNumberOfDropouts2 = sum(ptotal*ndropouts2);
  double expectedNumberOfSubjects2 = sum(ptotal*nsubjects2);
  double expectedStudyDuration = sum(ptotal*time);
  double expectedInformation = sum(ptotal*vscore);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  if (estimateHazardRatio) {
      hru = hazardRatioH0*exp(-criticalValues1*sqrt(vlogHR));
      hrl = hazardRatioH0*exp(-futilityBounds1*sqrt(vlogHR));
  }

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      hru[i] = NA_REAL;
      efficacyStopping1[i] = 0;
    }

    if (futilityBounds1[i] == -6) {
      hrl[i] = NA_REAL;
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
    _["information"] = (vscore[kMax-1]),
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfDropouts"] = expectedNumberOfDropouts,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["rho1"] = rho1,
    _["rho2"] = rho2,
    _["kMax"] = kMax,
    _["hazardRatioH0"] = hazardRatioH0,
    _["typeOfComputation"] = typeOfComputation);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["typeBetaSpending"] = typeBetaSpending,
    _["parameterBetaSpending"] = parameterBetaSpending,
    _["allocationRatioPlanned"] = allocationRatioPlanned,
    _["accrualTime"] = accrualTime,
    _["accrualIntensity"] = accrualIntensity,
    _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
    _["stratumFraction"] = stratumFraction,
    _["lambda1"] = lambda1,
    _["lambda2"] = lambda2,
    _["gamma1"] = gamma1,
    _["gamma2"] = gamma2,
    _["estimateHazardRatio"] = estimateHazardRatio,
    _["spendingTime"] = spendingTime);

  List byTreatmentCounts = List::create(
    _["numberOfEvents1"] = nevents1,
    _["numberOfDropouts1"] = ndropouts1,
    _["numberOfSubjects1"] = nsubjects1,
    _["numberOfEvents2"] = nevents2,
    _["numberOfDropouts2"] = ndropouts2,
    _["numberOfSubjects2"] = nsubjects2,
    _["expectedNumberOfEvents1"] = expectedNumberOfEvents1,
    _["expectedNumberOfDropouts1"] = expectedNumberOfDropouts1,
    _["expectedNumberOfSubjects1"] = expectedNumberOfSubjects1,
    _["expectedNumberOfEvents2"] = expectedNumberOfEvents2,
    _["expectedNumberOfDropouts2"] = expectedNumberOfDropouts2,
    _["expectedNumberOfSubjects2"] = expectedNumberOfSubjects2);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings,
    _["byTreatmentCounts"] = byTreatmentCounts);

  result.attr("class") = "lrpower";

  return result;
}


//' @title Power and sample size for a generic group sequential design
//' @description Obtains the maximum information and stopping boundaries
//' for a generic group sequential design assuming a constant treatment
//' effect, or obtains the power given the maximum information and
//' stopping boundaries.
//'
//' @param beta The type II error.
//' @param IMax The maximum information. Either \code{beta} or \code{IMax}
//'   should be provided while the other one should be missing.
//' @param theta The parameter value.
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
//'   variance under H1. Defaults to 1.
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
//'     - \code{calculationTarget}: The calculation target, \code{beta} or
//'       \code{IMax}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Christopher Jennison, Bruce W. Turnbull. Group Sequential Methods with
//' Applications to Clinical Trials. Chapman & Hall/CRC: Boca Raton, 2000,
//' ISBN:0849303168
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
List getDesign(const double beta = NA_REAL,
               const double IMax = NA_REAL,
               const double theta = NA_REAL,
               const int kMax = 1,
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
               const NumericVector& spendingTime = NA_REAL,
               const double varianceRatio = 1) {

  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector spendingTime1 = clone(spendingTime);

  double alpha1 = alpha;
  double beta1 = beta;
  double IMax1 = IMax;
  double drift, inflationFactor;

  String unknown;

  if (R_isnancpp(beta) && R_isnancpp(IMax)) {
    stop("beta and IMax cannot be both missing");
  }

  if (!R_isnancpp(beta) && !R_isnancpp(IMax)) {
    stop("Only one of beta and IMax should be provided");
  }

  if (!R_isnancpp(IMax)) {
    if (IMax <= 0) {
      stop("IMax must be positive");
    }
    unknown = "beta";
  } else if (!R_isnancpp(beta)) {
    unknown = "IMax";
  }

  if (R_isnancpp(theta)) {
    stop("theta must be provided");
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
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
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
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if ((unknown == "IMax") && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }


  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
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

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double bsfpar = parameterBetaSpending;

  if (unknown == "IMax") {
    if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
        bsf=="sfkd" || bsf=="sfhsd" || bsf=="user" || bsf=="none")) {
      stop("Invalid value for typeBetaSpending");
    }
  } else {
    if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
        bsf=="sfkd" || bsf=="sfhsd" || bsf=="none")) {
      stop("Invalid value for typeBetaSpending");
    }
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (unknown=="IMax" && bsf=="user") {
    if (is_true(any(is_na(userBetaSpending)))) {
      stop("userBetaSpending must be specified");
    } else if (userBetaSpending.size() < kMax) {
      stop("Insufficient length of userBetaSpending");
    } else if (userBetaSpending[0] < 0) {
      stop("Elements of userBetaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userBetaSpending) < 0))) {
      stop("Elements of userBetaSpending must be nondecreasing");
    } else if (userBetaSpending[kMax-1] != beta) {
      stop("userBetaSpending must end with specified beta");
    }
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (varianceRatio <= 0) {
    stop("varianceRatio must be positive");
  }


  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  alpha1 = sum(NumericVector(probs[0]));

  bool missingFutilityBounds = is_true(any(is_na(futilityBounds)));
  if (kMax > 1) {
    if (missingFutilityBounds && bsf=="none") {
      futilityBounds1 = rep(-6.0, kMax);
      futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    } else if (!missingFutilityBounds &&
      futilityBounds1.size() == kMax-1) {
      futilityBounds1.push_back(criticalValues1[kMax-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futilityBounds1 = criticalValues1[kMax-1];
    }
  }

  // multiplier for the boundaries under the alternative hypothesis
  NumericVector w = rep(sqrt(varianceRatio), kMax);

  NumericVector t = informationRates1;
  NumericVector st = spendingTime1;
  NumericVector theta1(kMax);
  if (unknown == "IMax") {
    auto f = [beta, kMax, t, futilityStopping1,
              criticalValues1, &futilityBounds1,
              bsf, bsfpar, userBetaSpending, st, w,
              missingFutilityBounds](double aval)->double {

                NumericVector theta1 = rep(aval, kMax);

                // compute stagewise exit probabilities
                if (!missingFutilityBounds || bsf=="none" || kMax==1) {
                  List probs = exitprobcpp(
                    criticalValues1*w, futilityBounds1*w, theta1, t);
                  NumericVector pu = NumericVector(probs[0]);
                  double overallReject = sum(pu);
                  return overallReject - (1-beta);
                } else {
                  // initialize futility bound to be updated
                  futilityBounds1 = NumericVector(kMax);
                  double epsilon;

                  // first stage
                  int k = 0;
                  double cumBetaSpent;
                  if (bsf=="user") {
                    cumBetaSpent = userBetaSpending[0];
                  } else {
                    cumBetaSpent = errorSpentcpp(st[0], beta, bsf, bsfpar);
                  }

                  if (!futilityStopping1[0]) {
                    futilityBounds1[0] = -6.0;
                  } else {
                    epsilon = R::pnorm(criticalValues1[0]*w[0] -
                      theta1[0]*sqrt(t[0]), 0, 1, 1, 0) - cumBetaSpent;
                    if (epsilon < 0) return -1.0;
                    futilityBounds1[0] = (R::qnorm(cumBetaSpent, 0, 1, 1, 0)
                                            + theta1[0]*sqrt(t[0]))/w[0];
                  }


                  // lambda expression for finding futility bound at stage k
                  auto g = [&k, &cumBetaSpent, criticalValues1,
                            &futilityBounds1, theta1, w,
                            t](double aval)->double {
                              NumericVector u(k+1), l(k+1);
                              for (int i=0; i<k; i++) {
                                u[i] = criticalValues1[i]*w[i];
                                l[i] = futilityBounds1[i]*w[i];
                              }
                              u[k] = 6.0;
                              l[k] = aval*w[k];

                              IntegerVector idx = Range(0,k);
                              List probs = exitprobcpp(
                                u, l, theta1[idx], t[idx]);
                              double cpl = sum(NumericVector(probs[1]));
                              return cpl - cumBetaSpent;
                            };


                  for (k=1; k<kMax; k++) {
                    if (bsf == "user") {
                      cumBetaSpent = userBetaSpending[k];
                    } else {
                      cumBetaSpent = errorSpentcpp(st[k], beta, bsf, bsfpar);
                    }

                    if (!futilityStopping1[k]) {
                      futilityBounds1[k] = -6.0;
                    } else {
                      epsilon = g(criticalValues1[k]);

                      if (g(-6.0) > 0) { // no beta spent at current visit
                        futilityBounds1[k] = -6.0;
                      } else if (epsilon > 0) {
                        futilityBounds1[k] = brent(
                          g, -6.0, criticalValues1[k], 1e-6);
                      } else if (k < kMax-1) {
                        return -1.0;
                      }
                    }
                  }

                  return epsilon;
                }
              };

    drift = brent(f, 0, 6, 0.0001);
    IMax1 = pow(drift/theta, 2);
    futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    theta1 = rep(drift, kMax);
    probs = exitprobcpp(criticalValues1*w, futilityBounds1*w, theta1, t);
  } else {
    drift = theta*sqrt(IMax1);
    theta1 = rep(drift, kMax);

    if (!missingFutilityBounds || bsf=="none" || kMax==1) {
      probs = exitprobcpp(criticalValues1*w, futilityBounds1*w, theta1, t);
      beta1 = 1 - sum(NumericVector(probs[0]));
    } else {
      List out = getPower(alpha1, kMax, criticalValues1, theta1, t,
                          bsf, bsfpar, st, futilityStopping1, w);

      beta1 = out[0];
      futilityBounds1 = out[1];
      probs = out[2];
    }
  }

  double driftf = R::qnorm(1-alpha1, 0, 1, 1, 0)*w[0] +
    R::qnorm(1-beta1, 0, 1, 1, 0);
  inflationFactor = pow(drift/driftf, 2);


  // output the results
  NumericVector information(kMax);
  NumericVector efficacyTheta(kMax);
  NumericVector futilityTheta(kMax);
  NumericVector efficacyP(kMax);
  NumericVector futilityP(kMax);
  for (int i=0; i<kMax; i++) {
    information[i] = IMax1*informationRates1[i];
    efficacyTheta[i] = criticalValues1[i]/sqrt(information[i])*w[i];
    futilityTheta[i] = futilityBounds1[i]/sqrt(information[i])*w[i];
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
    futilityP[i] = 1 - R::pnorm(futilityBounds1[i], 0, 1, 1, 0);
  }

  // stagewise exit probabilities under H1
  NumericVector pu(kMax), pl(kMax), ptotal(kMax);
  pu = NumericVector(probs[0]);
  pl = NumericVector(probs[1]);
  ptotal = pu + pl;

  double expectedInformationH1 = sum(ptotal*information);

  double overallReject = sum(pu);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  // cumulative alpha spent under H0 with non-binding futility
  NumericVector futilityBounds0(kMax, -6.0), theta0(kMax);
  List probs0 = exitprobcpp(criticalValues1, futilityBounds0, theta0, t);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs0[0]));

  // stagewise exit probabilities under H0 with binding futility
  probs0 = exitprobcpp(criticalValues1, futilityBounds1, theta0, t);
  NumericVector pu0(kMax), pl0(kMax), ptotal0(kMax);
  pu0 = NumericVector(probs0[0]);
  pl0 = NumericVector(probs0[1]);
  ptotal0 = pu0 + pl0;

  double expectedInformationH0 = sum(ptotal0*information);

  double overallRejectH0 = sum(pu0);
  NumericVector cpu0 = cumsum(pu0);
  NumericVector cpl0 = cumsum(pl0);

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      efficacyStopping1[i] = 0;
    }
    if (futilityBounds1[i] == -6) {
      futilityStopping1[i] = 0;
    }
  }

  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["futilityBounds"] = futilityBounds1,
    _["rejectPerStage"] = pu,
    _["futilityPerStage"] = pl,
    _["cumulativeRejection"] = cpu,
    _["cumulativeFutility"] = cpl,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["efficacyTheta"] = efficacyTheta,
    _["futilityTheta"] = futilityTheta,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["information"] = information,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1,
    _["rejectPerStageH0"] = pu0,
    _["futilityPerStageH0"] = pl0,
    _["cumulativeRejectionH0"] = cpu0,
    _["cumulativeFutilityH0"] = cpl0);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["attainedAlpha"] = overallRejectH0,
    _["kMax"] = kMax,
    _["theta"] = theta,
    _["information"] = IMax1,
    _["expectedInformationH1"] = expectedInformationH1,
    _["expectedInformationH0"] = expectedInformationH0,
    _["drift"] = drift,
    _["inflationFactor"] = inflationFactor);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["typeBetaSpending"] = typeBetaSpending,
    _["parameterBetaSpending"] = parameterBetaSpending,
    _["userBetaSpending"] = userBetaSpending,
    _["spendingTime"] = spendingTime,
    _["varianceRatio"] = varianceRatio,
    _["calculationTarget"] = unknown);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "design";

  return result;
}


//' @title Power and sample size for a generic group sequential equivalence
//' design
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
//'   tests. Defaults to 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param varianceRatioH10 The ratio of the variance under H10 to
//'   the variance under H1. Defaults to 1.
//' @param varianceRatioH20 The ratio of the variance under H20 to
//'   the variance under H1. Defaults to 1.
//' @param varianceRatioH12 The ratio of the variance under H10 to
//'   the variance under H20. Defaults to 1.
//' @param varianceRatioH21 The ratio of the variance under H20 to
//'   the variance under H10. Defaults to 1.
//'
//' @details
//' Consider the equivalence design with two one-sided hypotheses:
//' \deqn{H_{10}: \theta \leq \theta_{10},}
//' \deqn{H_{20}: \theta \geq \theta_{20}.}
//' We reject \eqn{H_{10}} at or before look \eqn{k} if
//' \deqn{Z_{1j} = (\hat{\theta}_j - \theta_{10})\sqrt{\frac{n_j}{v_{10}}}
//' \geq b_j}
//' for some \eqn{j=1,\ldots,k}, where \eqn{\{b_j:j=1,\ldots,K\}} are the
//' critical values associated with the specified alpha-spending function,
//' and \eqn{v_{10}} is the null variance of
//' \eqn{\hat{\theta}} based on the restricted maximum likelihood (reml)
//' estimate of model parameters subject to the constraint imposed by
//' \eqn{H_{10}} for one sampling unit drawn from \eqn{H_1}. For example,
//' for estimating the risk difference \eqn{\theta = \pi_1 - \pi_2},
//' the asymptotic limits of the
//' reml estimates of \eqn{\pi_1} and \eqn{\pi_2} subject to the constraint
//' imposed by \eqn{H_{10}} are given by
//' \deqn{(\tilde{\pi}_1, \tilde{\pi}_2) = f(\theta_{10}, r, r\pi_1,
//' 1-r, (1-r)\pi_2),}
//' where \eqn{f(\theta_0, n_1, y_1, n_2, y_2)} is the function to obtain
//' the reml of \eqn{\pi_1} and \eqn{\pi_2} subject to the constraint that
//' \eqn{\pi_1-\pi_2 = \theta_0} with observed data
//' \eqn{(n_1, y_1, n_2, y_2)} for the number of subjects and number of
//' responses in the active treatment and control groups,
//' \eqn{r} is the randomization probability for the active treatment
//' group, and \deqn{v_{10} = \frac{\tilde{\pi}_1 (1-\tilde{\pi}_1)}{r} +
//' \frac{\tilde{\pi}_2 (1-\tilde{\pi}_2)}{1-r}.}
//'
//' Let \eqn{I_j = n_j/v_1} denote the information for \eqn{\theta} at the
//' \eqn{j}th look, where
//' \deqn{v_{1} = \frac{\pi_1 (1-\pi_1)}{r} + \frac{\pi_2 (1-\pi_2)}{1-r}}
//' denotes the variance of \eqn{\hat{\theta}} under \eqn{H_1} for one
//' sampling unit. It follows that
//' \deqn{(Z_{1j} \geq b_j) = (Z_j \geq w_{10} b_j +
//' (\theta_{10}-\theta)\sqrt{I_j}),}
//' where \eqn{Z_j = (\hat{\theta}_j - \theta)\sqrt{I_j}}, and
//' \eqn{w_{10} = \sqrt{v_{10}/v_1}}.
//'
//' Similarly, we reject \eqn{H_{20}} at or before look \eqn{k} if
//' \deqn{Z_{2j} = (\hat{\theta}_j - \theta_{20})\sqrt{\frac{n_j}{v_{20}}}
//' \leq -b_j} for some \eqn{j=1,\ldots,k}, where \eqn{v_{20}} is the null
//' variance of \eqn{\hat{\theta}} based on the reml estimate of model
//' parameters subject to the constraint imposed by \eqn{H_{20}} for
//' one sampling unit drawn from \eqn{H_1}. We have
//' \deqn{(Z_{2j} \leq -b_j) = (Z_j \leq -w_{20} b_j +
//' (\theta_{20}-\theta)\sqrt{I_j}),}
//' where \eqn{w_{20} = \sqrt{v_{20}/v_1}}.
//'
//' Let \eqn{l_j = w_{10}b_j + (\theta_{10}-\theta)\sqrt{I_j}},
//' and \eqn{u_j = -w_{20}b_j + (\theta_{20}-\theta)\sqrt{I_j}}.
//' The cumulative probability to reject \eqn{H_0 = H_{10} \cup H_{20}} at
//' or before look \eqn{k} under the alternative hypothesis \eqn{H_1} is
//' given by
//' \deqn{P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j) \cap
//' \cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right) = p_1 + p_2 + p_{12},}
//' where
//' \deqn{p_1 = P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \geq l_j)\right),}
//' \deqn{p_2 = P_\theta\left(\cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \leq u_j)\right),}
//' and
//' \deqn{p_{12} = P_\theta\left(\cup_{j=1}^{k} \{(Z_j \geq l_j) \cup
//' (Z_j \leq u_j)\}\right).}
//' Of note, both \eqn{p_1} and \eqn{p_2} can be evaluated using
//' one-sided exit probabilities for group sequential designs.
//' If there exists \eqn{j\leq k} such that \eqn{l_j \leq u_j}, then
//' \eqn{p_{12} = 1}. Otherwise, \eqn{p_{12}} can be evaluated using
//' two-sided exit probabilities for group sequential designs.
//'
//' To evaluate the type I error of the equivalence trial under
//' \eqn{H_{10}}, we first match the information under \eqn{H_{10}}
//' with the information under \eqn{H_1}. For example, for estimating
//' the risk difference for two independent samples, the sample size
//' \eqn{n_{10}} under \eqn{H_{10}} must satisfy
//' \deqn{\frac{1}{n_{10}}\left(\frac{(\pi_2 + \theta_{10})
//' (1 - \pi_2 - \theta_{10})}{r} + \frac{\pi_2 (1-\pi_2)}{1-r}\right)
//' = \frac{1}{n}\left(\frac{\pi_1(1-\pi_1)}{r} +
//' \frac{\pi_2 (1-\pi_2)}{1-r}\right).}
//' Then we obtain the reml estimates of \eqn{\pi_1} and \eqn{\pi_2}
//' subject to the constraint imposed by \eqn{H_{20}} for one sampling
//' unit drawn from \eqn{H_{10}},
//' \deqn{(\tilde{\pi}_{10}, \tilde{\pi}_{20}) = f(\theta_{20}, r,
//' r(\pi_2 + \theta_{10}), 1-r, (1-r)\pi_2).}
//' Let \eqn{t_j} denote the information fraction at look \eqn{j}.
//' Define \deqn{\tilde{v}_1 = \frac{(\pi_2 + \theta_{10})
//' (1-\pi_2 -\theta_{10})}{r} + \frac{\pi_2 (1-\pi_2)}{1-r},} and
//' \deqn{\tilde{v}_{20} = \frac{\tilde{\pi}_{10}(1-\tilde{\pi}_{10})}{r} +
//' \frac{\tilde{\pi}_{20} (1-\tilde{\pi}_{20})}{1-r}.}
//'
//' The cumulative rejection probability under \eqn{H_{10}} at or before
//' look \eqn{k} is given by
//' \deqn{P_{\theta_{10}}\left(\cup_{j=1}^{k} \{(\hat{\theta}_j - \theta_{10})
//' \sqrt{n_{10} t_j/\tilde{v}_1} \geq b_j\} \cap
//' \cup_{j=1}^{k} \{(\hat{\theta}_j - \theta_{20})
//' \sqrt{n_{10} t_j/\tilde{v}_{20}} \leq -b_j\}\right) =
//' q_1 + q_2 + q_{12},}
//' where
//' \deqn{q_1 = P_{\theta_{10}}\left(\cup_{j=1}^{k}
//' \{(\hat{\theta}_j - \theta_{10})
//' \sqrt{n_{10} t_j/\tilde{v}_1} \geq b_j\}\right) =
//' P_{\theta_{10}}\left(\cup_{j=1}^{k} (Z_j \geq b_j)\right),}
//' \deqn{q_2 = P_{\theta_{10}}\left(\cup_{j=1}^{k}
//' \{(\hat{\theta}_j - \theta_{20})
//' \sqrt{n_{10} t_j/\tilde{v}_{20}} \leq -b_j\}\right) =
//' P_{\theta_{10}}\left(\cup_{j=1}^{k} (Z_j \leq -b_j w_{21} +
//' (\theta_{20} - \theta_{10})\sqrt{I_j})\right),}
//' and
//' \deqn{q_{12} = P_{\theta_{10}}\left(\cup_{j=1}^{k}
//' \{(Z_j \geq b_j) \cup (Z_j \leq -w_{21} b_j +
//' (\theta_{20} - \theta_{10})\sqrt{I_j})\}\right).}
//' Here \eqn{Z_j = (\hat{\theta}_j - \theta_{10}) \sqrt{I_j}}, and
//' \eqn{w_{21} = \sqrt{\tilde{v}_{20}/\tilde{v}_1}}.
//' Of note, \eqn{q_1}, \eqn{q_2}, and \eqn{q_{12}}
//' can be evaluated using group sequential exit probabilities.
//' Similarly, we can define \eqn{\tilde{v}_2}, \eqn{\tilde{v}_{10}},
//' and \eqn{w_{12} = \sqrt{\tilde{v}_{10}/\tilde{v}_2}}, and
//' evaluate the type I error under \eqn{H_{20}}.
//'
//' The variance ratios correspond to
//' \deqn{\text{varianceRatioH10} = v_{10}/v_1,}
//' \deqn{\text{varianceRatioH20} = v_{20}/v_1,}
//' \deqn{\text{varianceRatioH12} = \tilde{v}_{10}/\tilde{v}_2,}
//' \deqn{\text{varianceRatioH21} = \tilde{v}_{20}/\tilde{v}_1.}
//' If the alternative variance is used, then the variance ratios
//' are all equal to 1.
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
//'     - \code{varianceRatioH10}: The ratio of the variance under H10 to
//'       the variance under H1.
//'
//'     - \code{varianceRatioH20}: The ratio of the variance under H20 to
//'       the variance under H1.
//'
//'     - \code{varianceRatioH12}: The ratio of the variance under H10 to
//'       the variance under H20.
//'
//'     - \code{varianceRatioH21}: The ratio of the variance under H20 to
//'       the variance under H10.
//'
//'     - \code{calculationTarget}: The calculation target, \code{beta} or
//'       \code{IMax}.
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
List getDesignEquiv(const double beta = NA_REAL,
                    const double IMax = NA_REAL,
                    const double thetaLower = NA_REAL,
                    const double thetaUpper = NA_REAL,
                    const double theta = 0,
                    const int kMax = 1,
                    const NumericVector& informationRates = NA_REAL,
                    const NumericVector& criticalValues = NA_REAL,
                    const double alpha = 0.05,
                    const String typeAlphaSpending = "sfOF",
                    const double parameterAlphaSpending = NA_REAL,
                    const NumericVector& userAlphaSpending = NA_REAL,
                    const NumericVector& spendingTime = NA_REAL,
                    const double varianceRatioH10 = 1,
                    const double varianceRatioH20 = 1,
                    const double varianceRatioH12 = 1,
                    const double varianceRatioH21 = 1) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector spendingTime1 = clone(spendingTime);

  double IMax1 = IMax;

  String unknown;

  if (R_isnancpp(beta) && R_isnancpp(IMax)) {
    stop("beta and IMax cannot be both missing");
  }

  if (!R_isnancpp(beta) && !R_isnancpp(IMax)) {
    stop("Only one of beta and IMax should be provided");
  }

  if (!R_isnancpp(IMax)) {
    if (IMax <= 0) {
      stop("IMax must be positive");
    }
    unknown = "beta";
  } else if (!R_isnancpp(beta)) {
    unknown = "IMax";
  }

  if (R_isnancpp(thetaLower)) {
    stop("thetaLower must be provided");
  }

  if (R_isnancpp(thetaUpper)) {
    stop("thetaUpper must be provided");
  }

  if (thetaLower >= theta) {
    stop("thetaLower must be less than theta");
  }

  if (thetaUpper <= theta) {
    stop("thetaUpper must be greater than theta");
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

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if ((unknown == "IMax") && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }


  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (varianceRatioH10 <= 0) {
    stop("varianceRatioH10 must be positive");
  }

  if (varianceRatioH20 <= 0) {
    stop("varianceRatioH20 must be positive");
  }

  if (varianceRatioH12 <= 0) {
    stop("varianceRatioH12 must be positive");
  }

  if (varianceRatioH21 <= 0) {
    stop("varianceRatioH21 must be positive");
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      LogicalVector efficacyStopping1(kMax, 1);
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, li, zero, informationRates1);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));

  NumericVector efficacyP(kMax);
  for (int i=0; i<kMax; i++) {
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
  }

  double wH10 = sqrt(varianceRatioH10), wH20 = sqrt(varianceRatioH20);
  double wH12 = sqrt(varianceRatioH12), wH21 = sqrt(varianceRatioH21);

  // calculate cumulative rejection probability under H1
  NumericVector t = informationRates1;
  NumericVector b = criticalValues1;
  double deltaLower = thetaLower - theta;
  double deltaUpper = thetaUpper - theta;

  // obtain IMax if needed
  if (unknown == "IMax") {
    auto f = [beta, kMax, t, b, wH10, wH20, deltaLower, deltaUpper,
              li, ui, zero](double aval)->double {
                NumericVector I = t*aval;
                NumericVector l = b*wH10 + deltaLower*sqrt(I);
                NumericVector u = -b*wH20 + deltaUpper*sqrt(I);

                List probs1 = exitprobcpp(pmax(l, li), li, zero, I);
                List probs2 = exitprobcpp(ui, pmin(u, ui), zero, I);

                double cpl = sum(NumericVector(probs1[0]));
                double cpu = sum(NumericVector(probs2[1]));

                double power;
                if (is_true(any(l <= u))) {
                  power = cpl + cpu - 1;
                } else {
                  List a = exitprobcpp(l, u, zero, I);
                  double ca = sum(NumericVector(a[0]) + NumericVector(a[1]));
                  power = cpl + cpu - ca;
                }

                return power - (1-beta);
              };

    double z0 = R::qnorm(1-alpha, 0, 1, 1, 0);
    double z1 = R::qnorm(1-beta, 0, 1, 1, 0);
    double IMax10 = pow((z0*wH10 + z1)/deltaLower, 2);
    double IMax20 = pow((z0*wH20 + z1)/deltaUpper, 2);
    double IMaxLower = 0.5*std::min(IMax10, IMax20);
    double IMaxUpper = 1.5*std::max(IMax10, IMax20);
    IMax1 = brent(f, IMaxLower, IMaxUpper, 0.0001);
  }

  // obtain cumulative rejection probabilities under H1
  NumericVector I = t*IMax1;
  NumericVector l = b*wH10 + deltaLower*sqrt(I);
  NumericVector u = -b*wH20 + deltaUpper*sqrt(I);

  List probs1 = exitprobcpp(pmax(l, li), li, zero, I);
  List probs2 = exitprobcpp(ui, pmin(u, ui), zero, I);

  NumericVector cpl = cumsum(NumericVector(probs1[0]));
  NumericVector cpu = cumsum(NumericVector(probs2[1]));

  IntegerVector k = which(l >= u);
  NumericVector cp(kMax);
  if (k.size() == 0) {
    cp = cpl + cpu - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    List a = exitprobcpp(l[idx], u[idx], zero[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cp[i] = cpl[i] + cpu[i] - ca[i];
      } else {
        cp[i] = cpl[i] + cpu[i] - 1;
      }
    }
  }

  // incremental exit probabilities under H1
  NumericVector q(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      q[i] = cp[i];
    } else if (i<kMax-1) {
      q[i] = cp[i] - cp[i-1];
    } else {
      q[i] = 1 - cp[i-1];
    }
  }

  NumericVector rejectPerStage(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      rejectPerStage[i] = cp[i];
    } else {
      rejectPerStage[i] = cp[i] - cp[i-1];
    }
  }

  double overallReject = cp[kMax-1];
  double expectedInformationH1 = sum(q*I);

  NumericVector efficacyThetaLower = b/sqrt(I)*wH10 + thetaLower;
  NumericVector efficacyThetaUpper = -b/sqrt(I)*wH20 + thetaUpper;


  // cumulative rejection probability under H10
  NumericVector lH10 = b;
  NumericVector uH10 = -b*wH21 + (thetaUpper - thetaLower)*sqrt(I);
  List probs2H10 = exitprobcpp(ui, pmin(uH10, ui), zero, I);
  NumericVector cpuH10 = cumsum(NumericVector(probs2H10[1]));
  NumericVector cplH10 = cumAlphaSpent;

  IntegerVector kH10 = which(lH10 >= uH10);
  NumericVector cpH10(kMax);
  if (kH10.size() == 0) {
    cpH10 = cplH10 + cpuH10 - 1;
  } else {
    int K = max(kH10);
    IntegerVector idx = Range(0, K);
    List aH10 = exitprobcpp(lH10[idx], uH10[idx], zero[idx], I[idx]);
    NumericVector caH10 = cumsum(NumericVector(aH10[0]) +
      NumericVector(aH10[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH10[i] = cplH10[i] + cpuH10[i] - caH10[i];
      } else {
        cpH10[i] = cplH10[i] + cpuH10[i] - 1;
      }
    }
  }

  // incremental exit probabilities under H10
  NumericVector qH10(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      qH10[i] = cpH10[i];
    } else if (i<kMax-1) {
      qH10[i] = cpH10[i] - cpH10[i-1];
    } else {
      qH10[i] = 1 - cpH10[i-1];
    }
  }

  double attainedAlphaH10 = cpH10[kMax-1];
  double expectedInformationH10 = sum(qH10*I);


  // cumulative rejection probability under H20
  NumericVector lH20 = b*wH12 + (thetaLower - thetaUpper)*sqrt(I);
  NumericVector uH20 = -b;
  List probs1H20 = exitprobcpp(pmax(lH20, li), li, zero, I);
  NumericVector cplH20 = cumsum(NumericVector(probs1H20[0]));
  NumericVector cpuH20 = cumAlphaSpent;

  IntegerVector kH20 = which(lH20 >= uH20);
  NumericVector cpH20(kMax);
  if (kH20.size() == 0) {
    cpH20 = cplH20 + cpuH20 - 1;
  } else {
    int K = max(kH20);
    IntegerVector idx = Range(0, K);
    List aH20 = exitprobcpp(lH20[idx], uH20[idx], zero[idx], I[idx]);
    NumericVector caH20 = cumsum(NumericVector(aH20[0]) +
      NumericVector(aH20[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH20[i] = cplH20[i] + cpuH20[i] - caH20[i];
      } else {
        cpH20[i] = cplH20[i] + cpuH20[i] - 1;
      }
    }
  }

  // incremental exit probabilities under H20
  NumericVector qH20(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      qH20[i] = cpH20[i];
    } else if (i<kMax-1) {
      qH20[i] = cpH20[i] - cpH20[i-1];
    } else {
      qH20[i] = 1 - cpH20[i-1];
    }
  }

  double attainedAlphaH20 = cpH20[kMax-1];
  double expectedInformationH20 = sum(qH20*I);


  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["rejectPerStage"] = rejectPerStage,
    _["cumulativeRejection"] = cp,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["cumulativeAttainedAlphaH10"] = cpH10,
    _["cumulativeAttainedAlphaH20"] = cpH20,
    _["efficacyThetaLower"] = efficacyThetaLower,
    _["efficacyThetaUpper"] = efficacyThetaUpper,
    _["efficacyP"] = efficacyP,
    _["information"] = I);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = alpha,
    _["attainedAlphaH10"] = attainedAlphaH10,
    _["attainedAlphaH20"] = attainedAlphaH20,
    _["kMax"] = kMax,
    _["thetaLower"] = thetaLower,
    _["thetaUpper"] = thetaUpper,
    _["theta"] = theta,
    _["information"] = IMax1,
    _["expectedInformationH1"] = expectedInformationH1,
    _["expectedInformationH10"] = expectedInformationH10,
    _["expectedInformationH20"] = expectedInformationH20);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["spendingTime"] = spendingTime,
    _["varianceRatioH10"] = varianceRatioH10,
    _["varianceRatioH20"] = varianceRatioH20,
    _["varianceRatioH12"] = varianceRatioH12,
    _["varianceRatioH21"] = varianceRatioH21,
    _["calculationTarget"] = unknown);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "designEquiv";

  return result;
}


//' @title Adaptive design at an interim look
//' @description Obtains the conditional power for specified incremental
//' information given the interim results, parameter value, and
//' data-dependent changes in the error spending function, and the number
//' and spacing of interim looks. Conversely, obtains the incremental
//' information needed to attain a specified conditional power given
//' the interim results, parameter value, and data-dependent changes
//' in the error spending function, and the number and spacing of
//' interim looks.
//'
//' @param betaNew The type II error for the secondary trial.
//' @param INew The maximum information of the secondary trial. Either
//'   \code{betaNew} or \code{INew} should be provided while the other one
//'   should be missing.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of
//'   the primary trial.
//' @param theta The parameter value.
//' @param IMax The maximum information of the primary trial. Must be
//'   provided if \code{futilityBounds} is missing and
//'   \code{typeBetaSpending} is not equal to "none", or
//'   if conditional power calculation is desired.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to true
//'   if left unspecified.
//' @param futilityStopping Indicators of whether futility stopping is
//'   allowed at each stage of the primary trial. Defaults to true
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping for the primary trial.
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
//'   for the primary trial. Corresponds to Delta for "WT", rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param userAlphaSpending The user defined alpha spending for the primary
//'   trial. Cumulative alpha spent up to each stage.
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
//'   for the primary trial. Corresponds to rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to true
//'   if left unspecified.
//' @param futilityStoppingNew The indicators of whether futility stopping is
//'   allowed at each look of the secondary trial. Defaults to true
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
//'   for the secondary trial. Corresponds to Delta for "WT", rho for "sfKD",
//'   and gamma for "sfHSD".
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
//'   for the secondary trial. Corresponds to rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param userBetaSpendingNew The user defined cumulative beta spending.
//'   Cumulative beta spent up to each stage of the secondary trial.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRatesNew}.
//' @param varianceRatio The ratio of the variance under H0 to the
//'   variance under H1.
//'
//' @return An \code{adaptDesign} object with two list components:
//'
//' * \code{primaryTrial}: A list of selected information for the primary
//'   trial, including \code{L}, \code{zL}, \code{theta}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
//'   and \code{MullerSchafer}.
//'
//' * \code{secondaryTrial}: A \code{design} object for the secondary trial.
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
//' # original group sequential design with 90% power to detect delta = 6
//' delta = 6
//' sigma = 17
//' n = 282
//' (des1 = getDesign(IMax = n/(4*sigma^2), theta = delta, kMax = 3,
//'                   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'                   parameterAlphaSpending = -4))
//'
//' # interim look results
//' L = 1
//' n1 = n/3
//' delta1 = 4.5
//' sigma1 = 20
//' zL = delta1/sqrt(4/n1*sigma1^2)
//'
//' t = des1$byStageResults$informationRates
//'
//' # conditional power with sample size increase
//' (des2 = adaptDesign(
//'   betaNew = NA, INew = 420/(4*sigma1^2),
//'   L, zL, theta = delta1,
//'   IMax = n/(4*sigma1^2), kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4))
//'
//' # Muller & Schafer (2001) method to design the secondary trial:
//' # 3-look gamma(-2) spending with 84% power at delta = 4.5 and sigma = 20
//' (des2 = adaptDesign(
//'   betaNew = 0.16, INew = NA,
//'   L, zL, theta = delta1,
//'   IMax = n/(4*sigma1^2), kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpendingNew = "sfHSD",
//'   parameterAlphaSpendingNew = -2))
//'
//' # incremental sample size for sigma = 20
//' (nNew = 4*sigma1^2*des2$secondaryTrial$overallResults$information)
//'
//' @export
// [[Rcpp::export]]
List adaptDesign(double betaNew = NA_REAL,
                 double INew = NA_REAL,
                 const int L = NA_INTEGER,
                 const double zL = NA_REAL,
                 const double theta = NA_REAL,
                 const double IMax = NA_REAL,
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
                 const NumericVector& spendingTime = NA_REAL,
                 const bool MullerSchafer = 0,
                 const int kNew = NA_INTEGER,
                 const NumericVector& informationRatesNew = NA_REAL,
                 const LogicalVector& efficacyStoppingNew = NA_LOGICAL,
                 const LogicalVector& futilityStoppingNew = NA_LOGICAL,
                 const String typeAlphaSpendingNew = "sfOF",
                 const double parameterAlphaSpendingNew = NA_REAL,
                 const String typeBetaSpendingNew = "none",
                 const double parameterBetaSpendingNew = NA_REAL,
                 const NumericVector& userBetaSpendingNew = NA_REAL,
                 const NumericVector& spendingTimeNew = NA_REAL,
                 const double varianceRatio = 1) {

  NumericVector t = clone(informationRates);
  LogicalVector es = clone(efficacyStopping);
  LogicalVector fs = clone(futilityStopping);
  NumericVector b = clone(criticalValues);
  NumericVector a = clone(futilityBounds);
  NumericVector st = clone(spendingTime);

  NumericVector tNew = clone(informationRatesNew);
  LogicalVector esNew = clone(efficacyStoppingNew);
  LogicalVector fsNew = clone(futilityStoppingNew);
  NumericVector stNew = clone(spendingTimeNew);

  double alpha1 = alpha;

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double bsfpar = parameterBetaSpending;

  std::string asfNew = typeAlphaSpendingNew;
  std::for_each(asfNew.begin(), asfNew.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfparNew = parameterAlphaSpendingNew;

  std::string bsfNew = typeBetaSpendingNew;
  std::for_each(bsfNew.begin(), bsfNew.end(), [](char & c) {
    c = std::tolower(c);
  });

  double bsfparNew = parameterBetaSpendingNew;

  if (R_isnancpp(betaNew) && R_isnancpp(INew)) {
    stop("betaNew and INew cannot be both missing");
  }

  if (!R_isnancpp(betaNew) && !R_isnancpp(INew)) {
    stop("Only one of betaNew and INew should be provided");
  }

  if (!R_isnancpp(betaNew) && betaNew < 0.0001 && betaNew >= 1) {
    stop("betaNew must be greater than or equal to 0.0001 and less than 1");
  }

  if (!R_isnancpp(INew) && INew <= 0) {
    stop("INew must be positive");
  }

  if (R_isnancpp(L)) {
    stop("L must be provided");
  }

  if (L < 1) {
    stop("L must be a positive integer");
  }

  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }

  if (R_isnancpp(theta)) {
    stop("theta must be provided");
  }

  if (R_isnancpp(kMax)) {
    stop("kMax must be provided");
  }

  if (kMax <= L) {
    stop("kMax must be greater than L");
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
    t = as<NumericVector>(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    es = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
    }
  } else {
    fs = rep(1, kMax);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
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

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    st = clone(t);
  }


  if (MullerSchafer) {
    if (R_isnancpp(kNew)) {
      stop("kNew must be provided");
    }

    if (is_false(any(is_na(informationRatesNew)))) {
      if (informationRatesNew.size() != kNew) {
        stop("Invalid length for informationRatesNew");
      } else if (informationRatesNew[0] <= 0) {
        stop("Elements of informationRatesNew must be positive");
      } else if (kNew > 1 && is_true(any(diff(informationRatesNew) <= 0))) {
        stop("Elements of informationRatesNew must be increasing");
      } else if (informationRatesNew[kNew-1] != 1) {
        stop("informationRatesNew must end with 1");
      }
    } else {
      IntegerVector tem = seq_len(kNew);
      tNew = as<NumericVector>(tem)/(kNew+0.0);
    }

    if (is_false(any(is_na(efficacyStoppingNew)))) {
      if (efficacyStoppingNew.size() != kNew) {
        stop("Invalid length for efficacyStoppingNew");
      } else if (efficacyStoppingNew[kNew-1] != 1) {
        stop("efficacyStoppingNew must end with 1");
      } else if (is_false(all((efficacyStoppingNew == 1) |
        (efficacyStoppingNew == 0)))) {
        stop("Elements of efficacyStoppingNew must be 1 or 0");
      }
    } else {
      esNew = rep(1, kNew);
    }

    if (is_false(any(is_na(futilityStoppingNew)))) {
      if (futilityStoppingNew.size() != kNew) {
        stop("Invalid length for futilityStoppingNew");
      } else if (futilityStoppingNew[kNew-1] != 1) {
        stop("futilityStoppingNew must end with 1");
      } else if (is_false(all((futilityStoppingNew == 1) |
        (futilityStoppingNew == 0)))) {
        stop("Elements of futilityStoppingNew must be 1 or 0");
      }
    } else {
      fsNew = rep(1, kNew);
    }

    if (!(asfNew=="of" || asfNew=="p" || asfNew=="wt" ||
        asfNew=="sfof" || asfNew=="sfp" ||
        asfNew=="sfkd" || asfNew=="sfhsd" || asfNew=="none")) {
      stop("Invalid value for typeAlphaSpendingNew");
    }

    if ((asfNew=="wt" || asfNew=="sfkd" || asfNew=="sfhsd") &&
        R_isnancpp(asfparNew)) {
      stop("Missing value for parameterAlphaSpendingNew");
    }

    if (asfNew=="sfkd" && asfparNew <= 0) {
      stop ("parameterAlphaSpendingNew must be positive for sfKD");
    }

    if (R_isnancpp(INew) && !(bsfNew=="sfof" || bsfNew=="sfp" ||
        bsfNew=="sfkd" || bsfNew=="sfhsd" ||
        bsfNew=="user" || bsfNew=="none")) {
      stop("Invalid value for typeBetaSpendingNew");
    } else if (!(bsfNew=="sfof" || bsfNew=="sfp" || bsfNew=="sfkd" ||
      bsfNew=="sfhsd" || bsfNew=="none")) {
      stop("Invalid value for typeBetaSpendingNew");
    }

    if ((bsfNew=="sfkd" || bsfNew=="sfhsd") && R_isnancpp(bsfparNew)) {
      stop("Missing value for parameterBetaSpendingNew");
    }

    if (bsfNew=="sfkd" && bsfparNew <= 0) {
      stop ("parameterBetaSpendingNew must be positive for sfKD");
    }

    if (R_isnancpp(INew) && bsfNew=="user") {
      if (is_true(any(is_na(userBetaSpendingNew)))) {
        stop("userBetaSpendingNew must be specified");
      } else if (userBetaSpendingNew.size() < kNew) {
        stop("Insufficient length of userBetaSpendingNew");
      } else if (userBetaSpendingNew[0] < 0) {
        stop("Elements of userBetaSpendingNew must be nonnegative");
      } else if (kNew > 1 && is_true(any(diff(userBetaSpendingNew) < 0))) {
        stop("Elements of userBetaSpendingNew must be nondecreasing");
      } else if (userBetaSpendingNew[kNew] != betaNew) {
        stop("userBetaSpendingNew must end with specified betaNew");
      }
    }

    if (is_false(any(is_na(spendingTimeNew)))) {
      if (spendingTimeNew.size() != kNew) {
        stop("Invalid length for spendingTimeNew");
      } else if (spendingTimeNew[0] <= 0) {
        stop("Elements of spendingTimeNew must be positive");
      } else if (kNew > 1 && is_true(any(diff(spendingTimeNew) <= 0))) {
        stop("Elements of spendingTimeNew must be increasing");
      } else if (spendingTimeNew[kNew-1] != 1) {
        stop("spendingTimeNew must end with 1");
      }
    } else {
      stNew = clone(tNew);
    }
  }

  if (varianceRatio <= 0) {
    stop("varianceRatio must be positive");
  }

  NumericVector w = rep(sqrt(varianceRatio), kMax);


  // obtain critical values for the primary trial
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, t, es, criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!es[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, t);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      b[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      b = getBoundcpp(kMax, t, alpha, asf, asfpar, userAlphaSpending,
                      st, es);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(b, l, zero, t);
  alpha1 = sum(NumericVector(probs[0]));

  // obtain futility bounds for the primary trial
  if (kMax > 1) {
    if (is_true(any(is_na(futilityBounds))) && bsf=="none") {
      a = rep(-6.0, kMax);
      a[kMax-1] = b[kMax-1];
    } else if (is_false(any(is_na(futilityBounds))) && a.size() == kMax-1) {
      a.push_back(b[kMax-1]);
    }
  } else {
    if (is_true(any(is_na(futilityBounds)))) {
      a = b[kMax-1];
    }
  }

  if (is_true(any(is_na(a)))) {
    if (R_isnancpp(IMax)) {
      stop("IMax must be provided");
    }

    if (IMax <= 0) {
      stop("IMax must be positive");
    }

    NumericVector theta1(kMax, theta);
    List out = getPower(alpha1, kMax, b, theta1, IMax*t, bsf, bsfpar,
                        st, fs, w);
    a = out[1];
  }

  int k1 = kMax - L;
  double alphaNew, conditionalPower, predictivePower;

  NumericVector t1(k1), r1(k1), b1(k1), a1(k1, -6.0), theta0(k1);
  for (int l=0; l<k1; l++) {
    t1[l] = (t[l+L] - t[L-1])/(1 - t[L-1]);
    r1[l] = t[L-1]/t[l+L];
    b1[l] = (b[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
    if (!es[l+L]) b1[l] = 6.0;
  }

  // conditional type I error
  probs = exitprobcpp(b1, a1, theta0, t1);
  alphaNew = sum(NumericVector(probs[0]));

  // conditional power
  for (int l=0; l<k1; l++) {
    a1[l] = (a[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
    if (!fs[l+L]) a1[l] = -6.0;
  }

  if (!R_isnancpp(IMax)) {
    double sigma = 1/sqrt(IMax*t[L-1]);
    double mu = zL*sigma;
    NumericVector theta1(k1, mu);

    NumericVector I1(k1);
    for (int l=0; l<k1; l++) {
      I1[l] = IMax*(t[l+L] - t[L-1]);
    }

    probs = exitprobcpp(b1, a1, theta1, I1);
    conditionalPower = sum(NumericVector(probs[0]));

    // predictive power
    auto f = [k1, b1, a1, I1](double theta)->double {
      NumericVector theta1(k1, theta);
      List probs = exitprobcpp(b1, a1, theta1, I1);
      return sum(NumericVector(probs[0]));
    };

    double lower = mu - 6*sigma, upper = mu + 6*sigma;
    predictivePower = intnorm(f, mu, sigma, lower, upper);
  } else {
    conditionalPower = NA_REAL;
    predictivePower = NA_REAL;
  }

  List des1 = List::create(
    _["L"] = L,
    _["zL"] = zL,
    _["theta"] = theta,
    _["kMax"] = kMax,
    _["informationRates"] = t,
    _["efficacyBounds"] = b,
    _["futilityBounds"] = a,
    _["conditionalAlpha"] = alphaNew,
    _["conditionalPower"] = conditionalPower,
    _["predictivePower"] = predictivePower,
    _["MullerSchafer"] = MullerSchafer);


  List des2;

  if (!MullerSchafer) {
    IntegerVector idx = Range(L, kMax-1);
    LogicalVector esNew = es[idx];
    LogicalVector fsNew = fs[idx];

    des2 = getDesign(betaNew, INew, theta, k1, t1, esNew, fsNew,
                     b1, NA_REAL, typeAlphaSpendingNew,
                     parameterAlphaSpendingNew, 0,
                     a1, typeBetaSpendingNew, parameterBetaSpendingNew,
                     userBetaSpendingNew, stNew, varianceRatio);
  } else {
    if (!R_isnancpp(betaNew) && betaNew >= 1-alphaNew) {
      stop("betaNew must be less than 1 minus conditional type I error");
    }

    NumericVector b1New(kNew, NA_REAL), a1New(kNew, NA_REAL);

    des2 = getDesign(betaNew, INew, theta, kNew, tNew, esNew, fsNew,
                     b1New, alphaNew, typeAlphaSpendingNew,
                     parameterAlphaSpendingNew, 0,
                     a1New, typeBetaSpendingNew, parameterBetaSpendingNew,
                     userBetaSpendingNew, stNew, varianceRatio);
  }

  List result = List::create(
    _["primaryTrial"] = des1,
    _["secondaryTrial"] = des2);

  result.attr("class") = "adaptDesign";

  return result;
}


//' @title Get the required number of events given hazard ratio
//' @description Obtains the required number of events given the hazard
//' ratios under the null and alternative hypotheses for a group
//' sequential design.
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
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @inheritParams param_hazardRatioH0
//' @param hazardRatio Hazard ratio under the alternative hypothesis
//'   for the active treatment versus control. Defaults to 0.5.
//' @inheritParams param_allocationRatioPlanned
//' @param rounding Whether to round up the number of events.
//'   Defaults to 1 for rounding.
//'
//' @return The required number of events.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' getNeventsFromHazardRatio(
//'   beta = 0.2, kMax = 2,
//'   informationRates = c(0.5,1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   typeBetaSpending = "sfP",
//'   hazardRatio = 0.673)
//'
//' @export
// [[Rcpp::export]]
double getNeventsFromHazardRatio(
    const double beta = 0.2,
    const int kMax = 1,
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
    const NumericVector& spendingTime = NA_REAL,
    const double hazardRatioH0 = 1,
    const double hazardRatio = 0.5,
    const double allocationRatioPlanned = 1,
    const bool rounding = 1) {

  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (hazardRatioH0 <= 0) {
    stop("hazardRatioH0 must be positive");
  }

  if (hazardRatio <= 0) {
    stop("hazardRatio must be positive");
  }

  if (allocationRatioPlanned <= 0) {
    stop("allocationRatioPlanned must be positive");
  }

  double theta = fabs(-log(hazardRatio/hazardRatioH0));
  List design = getDesign(beta, NA_REAL, theta, kMax, informationRates,
                          efficacyStopping, futilityStopping,
                          criticalValues, alpha, typeAlphaSpending,
                          parameterAlphaSpending, userAlphaSpending,
                          futilityBounds, typeBetaSpending,
                          parameterBetaSpending, userBetaSpending,
                          spendingTime, 1);

  DataFrame byStageResults = as<DataFrame>(design["byStageResults"]);
  NumericVector information = byStageResults["information"];
  double maxInformation = information[kMax-1];
  double r1 = allocationRatioPlanned/(1+allocationRatioPlanned);
  double D = maxInformation/(r1*(1-r1));
  if (rounding) D = std::ceil(D);
  return D;
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
//' @param informationRates The information rates in terms of number
//'   of events for the conventional log-rank test and in terms of
//'   the actual information for weighted log-rank tests.
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
//'   accrualDuration, followupTime, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}. Adjustment
//'   may be needed for non-monotone relationship with study power.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size and events.
//'   Defaults to 1 for sample size rounding.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{lrpower} object under the
//' alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{lrpower} object under the
//' null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{lrpower}}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survival, and 5% dropout by
//' # the end of 1 year.
//'
//' # Example 1: Obtains accrual duration given power and follow-up duration
//'
//' lrsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
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
//'              accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
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
//'              accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12,
//'              accrualDuration = 22,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
List lrsamplesize(const double beta = 0.2,
                  const int kMax = 1,
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
                  double followupTime = NA_REAL,
                  const bool fixedFollowup = 0,
                  const double rho1 = 0,
                  const double rho2 = 0,
                  const int numSubintervals = 300,
                  const bool estimateHazardRatio = 1,
                  const String typeOfComputation = "direct",
                  const NumericVector& interval =
                    NumericVector::create(0.001, 240),
                    const NumericVector& spendingTime = NA_REAL,
                    const bool rounding = 1) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });

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
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
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
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
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

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="user" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(futilityBounds))) && bsf=="user") {
    if (is_true(any(is_na(userBetaSpending)))) {
      stop("userBetaSpending must be specified");
    } else if (userBetaSpending.size() < kMax) {
      stop("Insufficient length of userBetaSpending");
    } else if (userBetaSpending[0] < 0) {
      stop("Elements of userBetaSpending must be nonnegative");
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

  if (fixedFollowup && R_isnancpp(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
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

  std::string su = typeOfComputation;
  std::for_each(su.begin(), su.end(), [](char & c) {
    c = std::tolower(c);
  });

  if (su != "direct" && su != "schoenfeld") {
    stop("typeOfComputation must be direct or Schoenfeld");
  }

  if (su == "schoenfeld" && (rho1 != 0 || rho2 != 0)) {
    stop("Schoenfeld method can only be used for ordinary log-rank test");
  }

  double hazardRatio = 1;
  if (su == "schoenfeld") {
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

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  alpha1 = sum(NumericVector(probs[0]));

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

  if (su == "schoenfeld") {
    double delta = -log(hazardRatio/hazardRatioH0);

    List design = getDesign(
      beta, NA_REAL, delta, kMax, informationRates1,
      efficacyStopping1, futilityStopping1, criticalValues1,
      alpha1, asf, asfpar, userAlphaSpending, futilityBounds1,
      bsf, bsfpar, userBetaSpending, spendingTime1, 1);

    DataFrame byStageResults = as<DataFrame>(design["byStageResults"]);
    criticalValues1 = byStageResults["efficacyBounds"];
    futilityBounds1 = byStageResults["futilityBounds"];

    DataFrame overallResults = as<DataFrame>(design["overallResults"]);
    double maxInformation = overallResults["information"];
    double r1 = allocationRatioPlanned/(allocationRatioPlanned+1);
    double D = maxInformation/(r1*(1-r1));

    auto f = [hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              unknown, D](double aval)-> double{
                NumericVector accrualIntensity1 = clone(accrualIntensity);
                double dur1=0, dur2=0;

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

                // obtain the total number of events at study end
                NumericVector u0(1, dur1 + dur2);
                DataFrame lr = lrstat(
                  u0, hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  dur1, dur2, fixedFollowup, 0, 0, 1, 1);

                return sum(NumericVector(lr[2])) - D;
              };

    if (unknown == "accrualDuration") {
      accrualDuration = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 0.0001);
      accrualIntensity1 = aval*accrualIntensity;
    }
  } else {
    auto f = [beta, kMax, informationRates1,
              futilityStopping1, criticalValues1,
              &futilityBounds1, bsf, bsfpar, userBetaSpending,
              hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals,
              spendingTime1, unknown,
              missingFutilityBounds](double aval)->double {
                NumericVector accrualIntensity1 = clone(accrualIntensity);
                double dur1=0, dur2=0;

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

                NumericVector u0(1);
                DataFrame lr;
                NumericVector e0(kMax), time(kMax);

                double studyDuration1 = dur1 + dur2;
                u0[0] = studyDuration1;

                // obtain the timing of interim analysis
                if (rho1 == 0 && rho2 == 0) { // conventional log-rank test
                  lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
                              accrualTime, accrualIntensity,
                              piecewiseSurvivalTime, stratumFraction,
                              lambda1, lambda2, gamma1, gamma2,
                              dur1, dur2, fixedFollowup,
                              rho1, rho2, numSubintervals, 1);

                  e0 = sum(NumericVector(lr[2]))*informationRates1;
                  time = caltime(e0, allocationRatioPlanned,
                                 accrualTime, accrualIntensity,
                                 piecewiseSurvivalTime, stratumFraction,
                                 lambda1, lambda2, gamma1, gamma2,
                                 dur1, dur2, fixedFollowup);
                } else {
                  lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
                              accrualTime, accrualIntensity,
                              piecewiseSurvivalTime, stratumFraction,
                              lambda1, lambda2, gamma1, gamma2,
                              dur1, dur2, fixedFollowup,
                              rho1, rho2, numSubintervals, 2);

                  double maxInformation = sum(NumericVector(lr[12]));
                  double information1;

                  auto g = [hazardRatioH0, allocationRatioPlanned,
                            accrualTime, accrualIntensity,
                            piecewiseSurvivalTime, stratumFraction,
                            lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup,
                            rho1, rho2, numSubintervals,
                            &information1](double aval)->double {
                              NumericVector u0(1, aval);
                              DataFrame lr = lrstat(
                                u0, hazardRatioH0, allocationRatioPlanned,
                                accrualTime, accrualIntensity,
                                piecewiseSurvivalTime, stratumFraction,
                                lambda1, lambda2, gamma1, gamma2,
                                dur1, dur2, fixedFollowup,
                                rho1, rho2, numSubintervals, 2);
                              return sum(NumericVector(lr[12])) -
                                information1;
                            };

                  for (int i=0; i<kMax; i++) {
                    information1 = maxInformation*informationRates1[i];
                    time[i] = brent(g, 0.0001, studyDuration1, 0.0001);
                  };
                }

                // obtain the mean and variance of log-rank test score
                // statistic at each stage
                NumericVector theta(kMax), vscore(kMax);
                lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
                            accrualTime, accrualIntensity1,
                            piecewiseSurvivalTime, stratumFraction,
                            lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup,
                            rho1, rho2, numSubintervals, 2);

                NumericVector uscore = NumericVector(lr[11]);
                vscore = NumericVector(lr[12]);

                theta = -uscore/vscore;

                // information time and spending time
                NumericVector t = vscore / (vscore[kMax - 1]);
                NumericVector st = spendingTime1;

                // compute stagewise exit probabilities
                if (!missingFutilityBounds || bsf=="none" || kMax==1) {
                  List probs = exitprobcpp(criticalValues1, futilityBounds1,
                                           theta, vscore);
                  double overallReject = sum(NumericVector(probs[0]));
                  return overallReject - (1-beta);
                } else {
                  // initialize futility bounds to be updated
                  futilityBounds1 = NumericVector(kMax);
                  double epsilon;

                  // first stage
                  int k = 0;
                  double cumBetaSpent;
                  if (bsf == "user") {
                    cumBetaSpent = userBetaSpending[0];
                  } else {
                    cumBetaSpent = errorSpentcpp(st[0], beta, bsf, bsfpar);
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

                  // lambda expression for finding futility bound at stage k
                  auto g = [&k, &cumBetaSpent, criticalValues1,
                            &futilityBounds1, theta,
                            vscore](double aval)->double {
                              NumericVector u(k+1), l(k+1);
                              for (int i=0; i<k; i++) {
                                u[i] = criticalValues1[i];
                                l[i] = futilityBounds1[i];
                              }
                              u[k] = 6.0;
                              l[k] = aval;

                              IntegerVector idx = Range(0,k);
                              List probs = exitprobcpp(u, l, theta[idx],
                                                       vscore[idx]);
                              double cpl = sum(NumericVector(probs[1]));
                              return cpl - cumBetaSpent;
                            };

                  for (k=1; k<kMax; k++) {
                    if (bsf == "user") {
                      cumBetaSpent = userBetaSpending[k];
                    } else {
                      cumBetaSpent = errorSpentcpp(st[k], beta, bsf, bsfpar);
                    }

                    if (!futilityStopping1[k]) {
                      futilityBounds1[k] = -6.0;
                    } else {
                      epsilon = g(criticalValues1[k]);

                      if (g(-6.0) > 0) { // no beta spent at current visit
                        futilityBounds1[k] = -6.0;
                      } else if (epsilon > 0) {
                        futilityBounds1[k] = brent(
                          g, -6.0, criticalValues1[k], 1e-6);
                      } else if (k < kMax-1) {
                        return -1.0;
                      }
                    }
                  }

                  return epsilon;
                }
              };

    if (unknown == "accrualDuration") {
      accrualDuration = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 0.0001);
      accrualIntensity1 = aval*accrualIntensity;
    }
  }

  futilityBounds1[kMax-1] = criticalValues1[kMax-1];

  // output the results
  List resultH1, resultH0, result;

  if (rounding) {
    NumericVector u(1);
    u[0] = accrualDuration + followupTime;
    DataFrame lr = lrstat(u, hazardRatioH0, allocationRatioPlanned,
                          accrualTime, accrualIntensity1,
                          piecewiseSurvivalTime, stratumFraction,
                          lambda1, lambda2, gamma1, gamma2,
                          accrualDuration, followupTime, fixedFollowup,
                          0, 0, 1, 1);

    // round up the total number of events
    double D0 = sum(NumericVector(lr[2]));
    double D = std::ceil(D0);

    // adjust design parameters to obtain integer number of events
    double n0, n, studyDuration;
    if (!fixedFollowup) {
      n0 = sum(NumericVector(lr[1]));
      n = std::ceil(n0);

      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        double aval = n/n0;
        accrualIntensity1 = aval*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        u = getAccrualDurationFromN(ns, accrualTime, accrualIntensity1);
        accrualDuration = u[0];
      }

      // adjust follow-up time to obtain integer number of events
      auto h = [hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                D](double aval)->double {
                  NumericVector u(1);
                  u[0] = accrualDuration + aval*followupTime;

                  DataFrame lr = lrstat(
                    u, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, aval*followupTime, fixedFollowup,
                    0, 0, 1, 1);
                  return sum(NumericVector(lr[2])) - D;
                };

      double aval = brent(h, 0, 1.1, 1e-6);
      followupTime = aval*followupTime;
      studyDuration = accrualDuration + followupTime;
    } else {
      // adjust accrual intensity or duration to obtain int number of events
      if (unknown == "accrualIntensity") {
        double aval = D/D0;
        accrualIntensity1 = aval*accrualIntensity1;
      } else {
        auto h = [hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  D](double aval)->double {
                    NumericVector u(1);
                    u[0] = aval*accrualDuration + followupTime;

                    DataFrame lr = lrstat(
                      u, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      aval*accrualDuration, followupTime, fixedFollowup,
                      0, 0, 1, 1);
                    return sum(NumericVector(lr[2])) - D;
                  };

        double aval = brent(h, 1, 1.1, 1e-6);
        accrualDuration = aval*accrualDuration;
      }

      NumericVector u(1, accrualDuration);
      n0 = accrual(u, accrualTime, accrualIntensity1, accrualDuration)[0];

      // round up the sample size
      n = std::ceil(n0);

      if (unknown == "accrualIntensity") {
        double aval = n/n0;
        accrualIntensity1 = aval*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        u = getAccrualDurationFromN(ns, accrualTime, accrualIntensity1);
        accrualDuration = u[0];
      }

      // adjust study duration to obtain integer number of events
      auto h = [hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                D](double aval)->double {
                  NumericVector u(1);
                  u[0] = accrualDuration + aval*followupTime;

                  DataFrame lr = lrstat(
                    u, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup,
                    0, 0, 1, 1);
                  return sum(NumericVector(lr[2])) - D;
                };

      double aval = brent(h, 0, 1, 1e-6);
      studyDuration = accrualDuration + aval*followupTime;
    }

    // update information rates to calculate new boundaries
    NumericVector nevents(kMax), information(kMax), time(kMax);

    if (rho1 == 0 && rho2 == 0) {
      nevents = floor(D*informationRates1 + 0.5);
      informationRates1 = nevents/nevents[kMax-1];
    } else {
      u[0] = studyDuration;
      lr = lrstat(u, hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  rho1, rho2, numSubintervals, 2);

      double maxInformation = sum(NumericVector(lr[12]));
      double information1;

      auto f = [hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                rho1, rho2, numSubintervals,
                &information1](double aval)->double {
                  NumericVector u0(1, aval);
                  DataFrame lr = lrstat(
                    u0, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup,
                    rho1, rho2, numSubintervals, 2);
                  return sum(NumericVector(lr[12])) - information1;
                };

      for (int i=0; i<kMax; i++) {
        information1 = maxInformation*informationRates1[i];
        time[i] = brent(f, 0.0001, studyDuration, 0.0001);
      };

      lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  rho1, rho2, numSubintervals, 1);

      nevents = floor(NumericVector(lr[2]) + 0.5); // round up events
      time = caltime(nevents, allocationRatioPlanned,
                     accrualTime, accrualIntensity1,
                     piecewiseSurvivalTime, stratumFraction,
                     lambda1, lambda2, gamma1, gamma2,
                     accrualDuration, followupTime, fixedFollowup);

      lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  rho1, rho2, numSubintervals, 2);

      information = NumericVector(lr[12]);
      informationRates1 = information/maxInformation;
    }

    // recalculate boundaries
    if (bsf != "user") {
      resultH1 = lrpower(
        kMax, informationRates1,
        efficacyStopping1, futilityStopping1, criticalValues,
        alpha1, typeAlphaSpending, parameterAlphaSpending,
        userAlphaSpending, futilityBounds,
        typeBetaSpending, parameterBetaSpending, hazardRatioH0,
        allocationRatioPlanned, accrualTime, accrualIntensity1,
        piecewiseSurvivalTime, stratumFraction,
        lambda1, lambda2, gamma1, gamma2,
        accrualDuration, followupTime, fixedFollowup,
        rho1, rho2, numSubintervals, estimateHazardRatio,
        typeOfComputation, spendingTime, studyDuration);
    } else {
      resultH1 = lrpower(
        kMax, informationRates1,
        efficacyStopping1, futilityStopping1, criticalValues,
        alpha1, typeAlphaSpending, parameterAlphaSpending,
        userAlphaSpending, futilityBounds1,
        typeBetaSpending, parameterBetaSpending, hazardRatioH0,
        allocationRatioPlanned, accrualTime, accrualIntensity1,
        piecewiseSurvivalTime, stratumFraction,
        lambda1, lambda2, gamma1, gamma2,
        accrualDuration, followupTime, fixedFollowup,
        rho1, rho2, numSubintervals, estimateHazardRatio,
        typeOfComputation, spendingTime, studyDuration);
    }
  } else {
    double studyDuration = accrualDuration + followupTime;

    resultH1 = lrpower(
      kMax, informationRates1,
      efficacyStopping1, futilityStopping1, criticalValues1,
      alpha1, typeAlphaSpending, parameterAlphaSpending,
      userAlphaSpending, futilityBounds1,
      typeBetaSpending, parameterBetaSpending, hazardRatioH0,
      allocationRatioPlanned, accrualTime, accrualIntensity1,
      piecewiseSurvivalTime, stratumFraction,
      lambda1, lambda2, gamma1, gamma2,
      accrualDuration, followupTime, fixedFollowup,
      rho1, rho2, numSubintervals, estimateHazardRatio,
      typeOfComputation, spendingTime, studyDuration);
  }

  // obtain results under H0 by matching the total number of events
  // for conventional log-rank test and maximum information for
  // weighted log-rank tests
  DataFrame overallResults = as<DataFrame>(resultH1["overallResults"]);
  DataFrame byStageResults = as<DataFrame>(resultH1["byStageResults"]);
  double D = overallResults["numberOfEvents"];
  NumericVector information = byStageResults["information"];
  double maxInformation = information[kMax-1];
  double aval, studyDuration;

  if (!fixedFollowup) {
    auto h = [hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals,
              D, maxInformation](double aval)->double {
                NumericVector u(1);
                u[0] = accrualDuration + aval*followupTime;

                if (rho1 == 0 && rho2 == 0) {
                  DataFrame lr = lrstat(
                    u, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                    accrualDuration, aval*followupTime, fixedFollowup,
                    rho1, rho2, numSubintervals, 1);
                  return sum(NumericVector(lr[2])) - D;
                } else {
                  DataFrame lr = lrstat(
                    u, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                    accrualDuration, aval*followupTime, fixedFollowup,
                    rho1, rho2, numSubintervals, 2);
                  return sum(NumericVector(lr[12])) - maxInformation;
                }
              };

    if (h(0) < 0) { // adjust the follow-up time
      aval = brent(h, 0, 2, 1e-6);
      followupTime = aval*followupTime;
      studyDuration = accrualDuration + followupTime;
    } else { // adjust the accrual duration
      auto g = [hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda2, gamma1, gamma2,
                accrualDuration, fixedFollowup,
                rho1, rho2, numSubintervals,
                D, maxInformation](double aval)->double {
                  NumericVector u(1);
                  u[0] = aval*accrualDuration;

                  if (rho1 == 0 && rho2 == 0) {
                    DataFrame lr = lrstat(
                      u, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval*accrualDuration, 0, fixedFollowup,
                      rho1, rho2, numSubintervals, 1);
                    return sum(NumericVector(lr[2])) - D;
                  } else {
                    DataFrame lr = lrstat(
                      u, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval*accrualDuration, 0, fixedFollowup,
                      rho1, rho2, numSubintervals, 2);
                    return sum(NumericVector(lr[12])) - maxInformation;
                  }
                };

      aval = brent(g, 0.1, 1, 1e-6);
      accrualDuration = aval*accrualDuration;
      followupTime = 0;
      studyDuration = accrualDuration + followupTime;
    }
  } else { // fixed follow-up
    auto h = [hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals,
              D, maxInformation](double aval)->double {
                NumericVector u(1);
                u[0] = accrualDuration + aval*followupTime;

                if (rho1 == 0 && rho2 == 0) {
                  DataFrame lr = lrstat(
                    u, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup,
                    rho1, rho2, numSubintervals, 1);
                  return sum(NumericVector(lr[2])) - D;
                } else {
                  DataFrame lr = lrstat(
                    u, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup,
                    rho1, rho2, numSubintervals, 2);
                  return sum(NumericVector(lr[12])) - maxInformation;
                }
              };

    if (h(0) < 0) { // adjust the study duration
      aval = brent(h, 0, 2, 1e-6);
      studyDuration = accrualDuration + aval*followupTime;
    } else { // adjust the accrual duration
      auto g = [hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                rho1, rho2, numSubintervals,
                D, maxInformation](double aval)->double {
                  NumericVector u(1);
                  u[0] = aval*accrualDuration;

                  if (rho1 == 0 && rho2 == 0) {
                    DataFrame lr = lrstat(
                      u, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval*accrualDuration, followupTime, fixedFollowup,
                      rho1, rho2, numSubintervals, 1);
                    return sum(NumericVector(lr[2])) - D;
                  } else {
                    DataFrame lr = lrstat(
                      u, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval*accrualDuration, followupTime, fixedFollowup,
                      rho1, rho2, numSubintervals, 2);
                    return sum(NumericVector(lr[12])) - maxInformation;
                  }
                };

      aval = brent(g, 0.1, 1, 1e-6);
      accrualDuration = aval*accrualDuration;
      studyDuration = accrualDuration;
    }
  }

  // use the same stopping boundaries as under H1
  criticalValues1 = byStageResults["efficacyBounds"];
  futilityBounds1 = byStageResults["futilityBounds"];

  resultH0 = lrpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending, hazardRatioH0,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    rho1, rho2, numSubintervals, estimateHazardRatio,
    typeOfComputation, spendingTime, studyDuration);

  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}


//' @title Number of events and information for negative binomial rate ratio
//' @description Obtains the number of subjects accrued, number of events,
//' number of dropouts, number of subjects reaching the maximum
//' follow-up, total exposure, and information for log rate in each group,
//' rate ratio, variance, information, and Wald test statistic of
//' log rate ratio at given calendar times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of events and the Wald test statistic.
//' @param rateRatioH0 Rate ratio under the null hypothesis.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @param kappa1 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the active treatment group.
//' @param kappa2 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the control group.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group.
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param nullVariance Whether to calculate the variance for log rate ratio
//'   under the null hypothesis.
//'
//' @details
//' The probability mass function for a negative binomial distribution with
//' dispersion parameter \eqn{\kappa_i} and rate parameter \eqn{\lambda_i}
//' is given by
//' \deqn{P(Y_{ij} = y) = \frac{\Gamma(y+1/\kappa_i)}{\Gamma(1/\kappa_i) y!}
//' \left(\frac{1}{1 + \kappa_i \lambda_i t_{ij}}\right)^{1/\kappa_i}
//' \left(\frac{\kappa_i \lambda_i t_{ij}}
//' {1 + \kappa_i \lambda_i t_{ij}}\right)^{y},}
//' where \eqn{Y_{ij}} is the event count for subject \eqn{j} in
//' treatment group \eqn{i}, and \eqn{t_{ij}} is the exposure time for
//' the subject. If \eqn{\kappa_i=0}, the negative binomial distribution
//' reduces to the Poisson distribution.
//'
//' For treatment group \eqn{i}, let \eqn{\beta_i = \log(\lambda_i)}.
//' The likelihood for \eqn{\{(\kappa_i, \beta_i):i=1,2\}} can be written as
//' \deqn{l = \sum_{i=1}^{2}\sum_{j=1}^{n_{i}}
//' \{\log \Gamma(y_{ij} + 1/\kappa_i) - \log \Gamma(1/\kappa_i) + y_{ij}
//' (\log(\kappa_i) + \beta_i) - (y_{ij} + 1/\kappa_i)
//' \log(1+ \kappa_i \exp(\beta_i) t_{ij})\}.}
//' It follows that
//' \deqn{\frac{\partial l}{\partial \beta_i} = \sum_{j=1}^{n_i}
//' \left\{y_{ij} - (y_{ij} + 1/\kappa_i)
//' \frac{\kappa_i \exp(\beta_i) t_{ij}}
//' {1 + \kappa_i \exp(\beta_i)t_{ij}}\right\},}
//' and
//' \deqn{-\frac{\partial^2 l}{\partial \beta_i^2} =
//' \sum_{j=1}^{n_i} (y_{ij} + 1/\kappa_i) \frac{\kappa_i \lambda_i t_{ij}}
//' {(1 + \kappa_i \lambda_i t_{ij})^2}.}
//' The Fisher information for \eqn{\beta_i} is
//' \deqn{E\left(-\frac{\partial^2 l}{\partial \beta_i^2}\right)
//' = n_i E\left(\frac{\lambda_i t_{ij}}
//' {1 + \kappa_i \lambda_i t_{ij}}\right).}
//' In addition, we can show that
//' \deqn{E\left(-\frac{\partial^2 l}
//' {\partial \beta_i \partial \kappa_i}\right) = 0.}
//' Therefore, the variance of \eqn{\hat{\beta}_i} is
//' \deqn{Var(\hat{\beta}_i) = \frac{1}{n_i} \left\{
//' E\left(\frac{\lambda_i t_{ij}}{1 + \kappa_i \lambda_i t_{ij}}\right)
//' \right\}^{-1}.}
//'
//' To evaluate the integral, we need to obtain the distribution of the
//' exposure time,
//' \deqn{t_{ij} = \min(\tau - W_{ij}, C_{ij}, T_{fmax}),}
//' where \eqn{\tau} denotes the calendar time since trial start,
//' \eqn{W_{ij}} denotes the enrollment time for subject \eqn{j}
//' in treatment group \eqn{i}, \eqn{C_{ij}} denotes the time to dropout
//' after enrollment for subject \eqn{j} in treatment group \eqn{i}, and
//' \eqn{T_{fmax}} denotes the maximum follow-up time for
//' all subjects. Therefore,
//' \deqn{P(t_{ij} \geq t) = P(W_{ij} \leq \tau - t)P(C_{ij} \geq t)
//' I(t\leq T_{fmax}).}
//' Let \eqn{H} denote the distribution function of the enrollment time,
//' and \eqn{G_i} denote the survival function of the dropout time for
//' treatment group \eqn{i}. By the change of variables, we have
//' \deqn{E\left(\frac{\lambda_i t_{ij}}{1 + \kappa_i \lambda_i t_{ij}}
//' \right) = \int_{0}^{\tau \wedge T_{fmax}}
//' \frac{\lambda_i}{(1 + \kappa_i \lambda_i t)^2} H(\tau - t) G_i(t) dt.}
//' A numerical integration algorithm for a univariate function can be
//' used to evaluate the above integral.
//'
//' For the restricted maximum likelihood (reml) estimate of
//' \eqn{(\beta_1,\beta_2)} subject to the
//' constraint that \eqn{\beta_1 - \beta_2 = \Delta}, we express the
//' log-likelihood in terms of \eqn{(\beta_2,\Delta,\kappa_1,\kappa_2)},
//' and takes the derivative of the log-likelihood function with respect
//' to \eqn{\beta_2}. The resulting score equation has asymptotic limit
//' \deqn{E\left(\frac{\partial l}{\partial \beta_2}\right) = s_1 + s_2,}
//' where
//' \deqn{s_1 = n r E\left\{\lambda1_1 t_{1j} - \left(\lambda_1t_{1j}
//' + \frac{1}{\kappa_1}\right) \frac{\kappa_1 e^{\tilde{\beta}_2 +
//' \Delta}t_{1j}}{1 + \kappa_1 e^{\tilde{\beta}_2 +\Delta}t_{1j}}\right\},}
//' and
//' \deqn{s_2 = n (1-r) E\left\{\lambda_2 t_{2j} -
//' \left(\lambda_2 t_{2j} + \frac{1}{\kappa_2}\right)
//' \frac{\kappa_2 e^{\tilde{\beta}_2} t_{2j}}
//' {1 + \kappa_2 e^{\tilde{\beta}_2}t_{2j}}\right\}.}
//' Here \eqn{r} is the randomization probability for the active
//' treatment group. The asymptotic limit of the reml of \eqn{\beta_2}
//' is the solution \eqn{\tilde{\beta}_2} to
//' \eqn{E\left(\frac{\partial l}{\partial \beta_2}\right) = 0.}
//'
//' @return A list with two components:
//'
//' * \code{resultsUnderH1}: A data frame containing the following variables:
//'
//'     - \code{time}: The analysis time since trial start.
//'
//'     - \code{subjects}: The number of enrolled subjects.
//'
//'     - \code{nevents}: The total number of events.
//'
//'     - \code{nevents1}: The number of events in the active treatment
//'       group.
//'
//'     - \code{nevents2}: The number of events in the control group.
//'
//'     - \code{ndropouts}: The total number of dropouts.
//'
//'     - \code{ndropouts1}: The number of dropouts in the active treatment
//'       group.
//'
//'     - \code{ndropouts2}: The number of dropouts in the control group.
//'
//'     - \code{nfmax}: The total number of subjects reaching maximum
//'       follow-up.
//'
//'     - \code{nfmax1}: The number of subjects reaching maximum follow-up
//'       in the active treatment group.
//'
//'     - \code{nfmax2}: The number of subjects reaching maximum follow-up
//'       in the control group.
//'
//'     - \code{exposure}: The total exposure time.
//'
//'     - \code{exposure1}: The exposure time for the active treatment group.
//'
//'     - \code{exposure2}: The exposure time for the control group.
//'
//'     - \code{information1}: The Fisher information for the log rate
//'       parameter for the active treatment group.
//'
//'     - \code{information2}: The Fisher information for the log rate
//'       parameter for the control group.
//'
//'     - \code{rateRatio}: The rate ratio of the active treatment group
//'       versus the control group.
//'
//'     - \code{vlogRR}: The variance of log rate ratio.
//'
//'     - \code{information}: The information of log rate ratio.
//'
//'     - \code{zlogRR}: The Z-statistic for log rate ratio.
//'
//' * \code{resultsUnderH0} when \code{nullVariance = TRUE}: A data frame
//'   with the following variables:
//'
//'     - \code{time}: The analysis time since trial start.
//'
//'     - \code{lambda1H0}: The restricted maximum likelihood estimate
//'       of the event rate for the active treatment group.
//'
//'     - \code{lambda2H0}: The restricted maximum likelihood estimate
//'       of the event rate for the control group.
//'
//'     - \code{information1H0}: The Fisher information for the log rate
//'       parameter for the active treatment group under H0.
//'
//'     - \code{information2H0}: The Fisher information for the log rate
//'       parameter for the control group under H0.
//'
//'     - \code{vlogRRH0}: The variance of log rate ratio under H0.
//'
//'     - \code{informationH0}: The information of log rate ratio under H0.
//'
//'     - \code{zlogRRH0}: The Z-statistic for log rate ratio under H0.
//'
//'     - \code{varianceRatio}: The ratio of the variance under H0 versus
//'       the variance under H1.
//'
//'     - \code{rateRatioH0}: The rate ratio under H0.
//'
//'     - \code{lambda1}: The true event rate for the active treatment group.
//'
//'     - \code{lambda2}: The true event rate for the control group.
//'
//'     - \code{rateRatio}: The true rate ratio.
//'
//' * \code{resultsUnderH0} when \code{nullVariance = FALSE}: A data frame
//'   with the following variables:
//'
//'     - \code{time}: The analysis time since trial start.
//'
//'     - \code{varianceRatio}: Equal to 1.
//'
//'     - \code{rateRatioH0}: The rate ratio under H0.
//'
//'     - \code{lambda1}: The true event rate for the active treatment group.
//'
//'     - \code{lambda2}: The true event rate for the control group.
//'
//'     - \code{rateRatio}: The true rate ratio.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Example 1: Variable follow-up design
//'
//' nbstat(time = c(1, 1.25, 2, 3, 4),
//'        accrualIntensity = 1956/1.25,
//'        kappa1 = 5,
//'        kappa2 = 5,
//'        lambda1 = 0.7*0.125,
//'        lambda2 = 0.125,
//'        gamma1 = 0,
//'        gamma2 = 0,
//'        accrualDuration = 1.25,
//'        followupTime = 2.75)
//'
//' # Example 2: Fixed follow-up design
//'
//' nbstat(time = c(0.5, 1, 1.5, 2),
//'        accrualIntensity = 220/1.5,
//'        kappa1 = 3,
//'        kappa2 = 3,
//'        lambda1 = 0.5*8.4,
//'        lambda2 = 8.4,
//'        gamma1 = 0,
//'        gamma2 = 0,
//'        accrualDuration = 1.5,
//'        followupTime = 0.5,
//'        fixedFollowup = 1,
//'        nullVariance = 1)
//'
//' @export
// [[Rcpp::export]]
List nbstat(const NumericVector& time = NA_REAL,
            const double rateRatioH0 = 1,
            const double allocationRatioPlanned = 1,
            const NumericVector& accrualTime = 0,
            const NumericVector& accrualIntensity = NA_REAL,
            const NumericVector& piecewiseSurvivalTime = 0,
            const double kappa1 = NA_REAL,
            const double kappa2 = NA_REAL,
            const double lambda1 = NA_REAL,
            const double lambda2 = NA_REAL,
            const NumericVector& gamma1 = 0,
            const NumericVector& gamma2 = 0,
            const double accrualDuration = NA_REAL,
            const double followupTime = NA_REAL,
            const bool fixedFollowup = 0,
            const bool nullVariance = 0) {

  int nintervals = piecewiseSurvivalTime.size();
  NumericVector lam1(nintervals), lam2(nintervals);
  NumericVector gam1(nintervals), gam2(nintervals);

  if (is_true(any(time < 0))) {
    stop("time must be non-negative");
  }

  if (rateRatioH0 <= 0) {
    stop("rateRatioH0 must be positive");
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

  if (kappa1 < 0) {
    stop("kappa1 must be nonnegative");
  }

  if (kappa2 < 0) {
    stop("kappa2 must be nonnegative");
  }

  if (lambda1 < 0) {
    stop("lambda1 must be positive");
  }

  if (lambda2 < 0) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (gamma1.size() == 1) {
    gam1 = rep(gamma1, nintervals);
  } else if (gamma1.size() == nintervals) {
    gam1 = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() == 1) {
    gam2 = rep(gamma2, nintervals);
  } else if (gamma2.size() == nintervals) {
    gam2 = gamma2;
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

  double minFollowupTime = followupTime;
  double maxFollowupTime;

  // obtain the follow-up time for the first enrolled subject
  if (fixedFollowup) {
    maxFollowupTime = minFollowupTime;
  } else {
    maxFollowupTime = accrualDuration + minFollowupTime;
  }

  int k = time.size();
  NumericVector nsubjects(k), nevents(k), nevents1(k), nevents2(k);
  NumericVector ndropouts(k), ndropouts1(k), ndropouts2(k);
  NumericVector nfmax(k), nfmax1(k), nfmax2(k);
  NumericVector exposure(k), exposure1(k), exposure2(k);
  NumericVector information(k), information1(k), information2(k);
  NumericVector rateRatio(k), vlogRR(k), zlogRR(k);

  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);
  NumericVector maxFU = NumericVector::create(maxFollowupTime);

  // number of subjects enrolled
  NumericVector s1 = pmin(time, accrualDuration + minFollowupTime);
  nsubjects = accrual(s1, accrualTime, accrualIntensity, accrualDuration);

  // number of dropouts
  if (max(gam1) + max(gam2) > 0.0) {
    NumericMatrix y = nevent2(s1, allocationRatioPlanned, accrualTime,
                              accrualIntensity, piecewiseSurvivalTime,
                              gam1, gam2, lam1, lam2, accrualDuration,
                              minFollowupTime, maxFollowupTime);
    ndropouts1 = y(_,0);
    ndropouts2 = y(_,1);
    ndropouts = ndropouts1 + ndropouts2;
  }

  // number of subjects reaching maximum follow-up adjusting for dropouts
  double p1 = patrisk(maxFU, piecewiseSurvivalTime, lam1, gam1)[0];
  double p2 = patrisk(maxFU, piecewiseSurvivalTime, lam2, gam2)[0];
  NumericVector s = s1 - maxFollowupTime;
  NumericVector n = accrual(s, accrualTime, accrualIntensity,
                            accrualDuration);
  nfmax1 = phi*n*p1;
  nfmax2 = (1-phi)*n*p2;
  nfmax = nfmax1 + nfmax2;

  double tau = 0;
  double r = phi, kappa = kappa1, lambda = lambda1;
  NumericVector lam = lam1, gam = gam1;

  auto f_ex = [&tau, &r, accrualTime, accrualIntensity,
               piecewiseSurvivalTime, &lam, &gam,
               accrualDuration](double t)->double {
                 NumericVector s = NumericVector::create(tau - t);
                 double a = accrual(s, accrualTime, accrualIntensity,
                                    accrualDuration)[0];
                 s[0] = t;
                 double p = patrisk(s, piecewiseSurvivalTime, lam,
                                    gam)[0];
                 return r*a*p;
               };


  auto f_info = [&tau, &r, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, &kappa, &lambda, &lam, &gam,
                 accrualDuration](double t)->double {
                   double u = lambda/pow(1 + kappa*lambda*t, 2);
                   NumericVector s = NumericVector::create(tau - t);
                   double a = accrual(s, accrualTime, accrualIntensity,
                                      accrualDuration)[0];
                   s[0] = t;
                   double p = patrisk(s, piecewiseSurvivalTime, lam,
                                      gam)[0];
                   return r*u*a*p;
                 };

  for (int i=0; i<k; i++) {
    tau = time[i];
    r = phi; kappa = kappa1; lambda = lambda1;
    lam = lam1; gam = gam1;
    exposure1[i] = qromb(f_ex, 0.0, std::min(tau, maxFollowupTime));
    information1[i] = qromb(f_info, 0.0, std::min(tau, maxFollowupTime));

    r = 1-phi; kappa = kappa2; lambda = lambda2;
    lam = lam2; gam = gam2;
    exposure2[i] = qromb(f_ex, 0.0, std::min(tau, maxFollowupTime));
    information2[i] = qromb(f_info, 0.0, std::min(tau, maxFollowupTime));
  }

  exposure = exposure1 + exposure2;
  nevents1 = lambda1*exposure1;
  nevents2 = lambda2*exposure2;
  nevents = nevents1 + nevents2;

  rateRatio = rep(lambda1/lambda2, k);
  vlogRR = 1/information1 + 1/information2;
  information = 1/vlogRR;
  zlogRR = log(rateRatio/rateRatioH0)/sqrt(vlogRR);

  DataFrame resultsUnderH1, resultsUnderH0;
  if (!nullVariance) {
    resultsUnderH1 = DataFrame::create(
      _["time"] = time,
      _["subjects"] = nsubjects,
      _["nevents"] = nevents,
      _["nevents1"] = nevents1,
      _["nevents2"] = nevents2,
      _["ndropouts"] = ndropouts,
      _["ndropouts1"] = ndropouts1,
      _["ndropouts2"] = ndropouts2,
      _["nfmax"] = nfmax,
      _["nfmax1"] = nfmax1,
      _["nfmax2"] = nfmax2,
      _["exposure"] = exposure,
      _["exposure1"] = exposure1,
      _["exposure2"] = exposure2,
      _["information1"] = information1,
      _["information2"] = information2,
      _["rateRatio"] = rateRatio,
      _["vlogRR"] = vlogRR,
      _["information"] = information,
      _["zlogRR"] = zlogRR);

    resultsUnderH0 = DataFrame::create(
      _["time"] = time,
      _["varianceRatio"] = rep(1.0, k),
      _["rateRatioH0"] = rep(rateRatioH0, k),
      _["lambda1"] = rep(lambda1, k),
      _["lambda2"] = rep(lambda2, k),
      _["rateRatio"] = rateRatio);
  } else {
    NumericVector lambda1H0(k), lambda2H0(k);
    NumericVector information1H0(k), information2H0(k);

    auto f = [&tau, rateRatioH0, phi, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, lam1, lam2, gam1, gam2,
              accrualDuration, maxFollowupTime](double aval)->double {
                double r = phi, kappa = kappa1, lambda = aval*rateRatioH0;
                NumericVector lam = lam1, gam = gam1;

                auto g = [&tau, &r, accrualTime,
                          accrualIntensity, piecewiseSurvivalTime,
                          &kappa, &lambda, &lam, &gam,
                          accrualDuration](double t)->double {
                            double u = 1/pow(1 + kappa*lambda*t, 2);
                            NumericVector s = NumericVector::create(tau-t);
                            double a = accrual(s, accrualTime,
                                               accrualIntensity,
                                               accrualDuration)[0];
                            s[0] = t;
                            double p = patrisk(s, piecewiseSurvivalTime,
                                               lam, gam)[0];
                            return r*u*a*p;
                          };

                double a1 = qromb(g, 0.0, std::min(tau, maxFollowupTime));

                r = 1-phi; kappa = kappa2; lambda = aval;
                lam = lam2; gam = gam2;
                double a2 = qromb(g, 0.0, std::min(tau, maxFollowupTime));

                return phi*(lambda1 - aval*rateRatioH0)*a1 +
                  (1-phi)*(lambda2 - aval)*a2;
              };


    for (int i=0; i<k; i++) {
      tau = time[i];

      double t1 = exposure1[i], t2 = exposure2[i];
      double a = (phi*kappa2 + (1-phi)*kappa1)*rateRatioH0*t1*t2;
      double b = -(phi*t1*(kappa2*lambda1*t2 - rateRatioH0) +
                   (1-phi)*t2*(kappa1*lambda2*rateRatioH0*t1 - 1));
      double c = -(phi*lambda1*t1 + (1-phi)*lambda2*t2);

      double init;
      if (kappa1 == 0 && kappa2 == 0) {
        init = -c/b;
      } else {
        init = (-b + sqrt(b*b - 4*a*c))/(2*a);
      }

      lambda2H0[i] = brent(f, 0.5*init, 1.5*init, 1e-6);
      lambda1H0[i] = rateRatioH0*lambda2H0[i];

      r = phi; kappa = kappa1; lambda = lambda1H0[i];
      lam = lam1; gam = gam1;
      information1H0[i] = qromb(f_info, 0.0, std::min(tau, maxFollowupTime));

      r = 1-phi; kappa = kappa2; lambda = lambda2H0[i];
      lam = lam2; gam = gam2;
      information2H0[i] = qromb(f_info, 0.0, std::min(tau, maxFollowupTime));
    }

    NumericVector vlogRRH0 = 1/information1H0 + 1/information2H0;
    NumericVector informationH0 = 1/vlogRRH0;
    NumericVector zlogRRH0 = log(rateRatio/rateRatioH0)/sqrt(vlogRRH0);
    NumericVector varianceRatio = vlogRRH0/vlogRR;


    resultsUnderH1 = DataFrame::create(
      _["time"] = time,
      _["subjects"] = nsubjects,
      _["nevents"] = nevents,
      _["nevents1"] = nevents1,
      _["nevents2"] = nevents2,
      _["ndropouts"] = ndropouts,
      _["ndropouts1"] = ndropouts1,
      _["ndropouts2"] = ndropouts2,
      _["nfmax"] = nfmax,
      _["nfmax1"] = nfmax1,
      _["nfmax2"] = nfmax2,
      _["exposure"] = exposure,
      _["exposure1"] = exposure1,
      _["exposure2"] = exposure2,
      _["information1"] = information1,
      _["information2"] = information2,
      _["rateRatio"] = rateRatio,
      _["vlogRR"] = vlogRR,
      _["information"] = information,
      _["zlogRR"] = zlogRR);

    resultsUnderH0 = DataFrame::create(
      _["time"] = time,
      _["lambda1H0"] = lambda1H0,
      _["lambda2H0"] = lambda2H0,
      _["information1H0"] = information1H0,
      _["information2H0"] = information2H0,
      _["vlogRRH0"] = vlogRRH0,
      _["informationH0"] = informationH0,
      _["zlogRRH0"] = zlogRRH0,
      _["varianceRatio"] = varianceRatio,
      _["rateRatioH0"] = rep(rateRatioH0, k),
      _["lambda1"] = rep(lambda1, k),
      _["lambda2"] = rep(lambda2, k),
      _["rateRatio"] = rateRatio);
  }

  return List::create(
    _["resultsUnderH1"] = resultsUnderH1,
    _["resultsUnderH0"] = resultsUnderH0);
}


//' @title Negative binomial rate ratio power
//' @description Estimates the power for negative binomial rate ratio test.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @param typeBetaSpending The type of beta spending. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock
//'   type spending function, "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
//'   early futility stopping. Defaults to "none".
//' @inheritParams param_parameterBetaSpending
//' @param rateRatioH0 Rate ratio under the null hypothesis.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @param kappa1 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the active treatment group.
//' @param kappa2 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the control group.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group.
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//' @param nullVariance Whether to calculate the variance for log rate ratio
//'   under the null hypothesis.
//'
//' @return An S3 class \code{nbpower} object with 4 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numberOfDropouts}: The total number of dropouts.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{exposure}: The total exposure.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfDropouts}: The expected number of dropouts.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedExposure}: The expected exposure.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{rateRatioH0}: The rate ratio under the null hypothesis.
//'
//'     - \code{rateRatio}: The rate ratio.
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
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{exposure}: The exposure.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRateRatio}: The efficacy boundaries on the rate
//'       ratio scale.
//'
//'     - \code{futilityRateRatio}: The futility boundaries on the rate
//'       ratio scale.
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
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{typeBetaSpending},
//'   \code{parameterBetaSpending}, \code{allocationRatioPlanned},
//'   \code{accrualTime}, \code{accuralIntensity},
//'   \code{piecewiseSurvivalTime}, \code{kappa1}, \code{kappa2},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   \code{accrualDuration}, \code{followupTime}, \code{fixedFollowup},
//'   \code{spendingTime}, and \code{nullVariance}.
//'
//' * \code{byTreatmentCounts}: A list containing the following counts by
//'   treatment group:
//'
//'     - \code{numberOfEvents1}: The number of events by stage for
//'       the treatment group.
//'
//'     - \code{numberOfDropouts1}: The number of dropouts by stage for
//'       the treatment group.
//'
//'     - \code{numberOfSubjects1}: The number of subjects by stage for
//'       the treatment group.
//'
//'     - \code{exposure1}: The exposure by stage for the treatment group.
//'
//'     - \code{numberOfEvents2}: The number of events by stage for
//'       the control group.
//'
//'     - \code{numberOfDropouts2}: The number of dropouts by stage for
//'       the control group.
//'
//'     - \code{numberOfSubjects2}: The number of subjects by stage for
//'       the control group.
//'
//'     - \code{exposure2}: The exposure by stage for the control group.
//'
//'     - \code{expectedNumberOfEvents1}: The expected number of events for
//'       the treatment group.
//'
//'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
//'       for the treatment group.
//'
//'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
//'       for the treatment group.
//'
//'     - \code{expectedExposure1}: The expected exposure for the treatment
//'       group.
//'
//'     - \code{expectedNumberOfEvents2}: The expected number of events for
//'       control group.
//'
//'     - \code{expectedNumberOfDropouts2}: The expected number of dropouts
//'       for the control group.
//'
//'     - \code{expectedNumberOfSubjects2}: The expected number of subjects
//'       for the control group.
//'
//'     - \code{expectedExposure2}: The expected exposure for the control
//'       group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{nbstat}}
//'
//' @examples
//' # Example 1: Variable follow-up design
//'
//' nbpower(kMax = 2, informationRates = c(0.5, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         accrualIntensity = 1956/1.25,
//'         kappa1 = 5, kappa2 = 5,
//'         lambda1 = 0.0875, lambda2 = 0.125,
//'         gamma1 = 0, gamma2 = 0,
//'         accrualDuration = 1.25,
//'         followupTime = 2.75, fixedFollowup = FALSE,
//'         nullVariance = 1)
//'
//' # Example 2: Fixed follow-up design
//'
//' nbpower(kMax = 2, informationRates = c(0.5, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         accrualIntensity = 220/1.5,
//'         kappa1 = 3, kappa2 = 3,
//'         lambda1 = 0.5*8.4, lambda2 = 8.4,
//'         gamma1 = 0, gamma2 = 0,
//'         accrualDuration = 1.5,
//'         followupTime = 0.5, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List nbpower(const int kMax = 1,
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
             const double rateRatioH0 = 1,
             const double allocationRatioPlanned = 1,
             const NumericVector& accrualTime = 0,
             const NumericVector& accrualIntensity = 1500,
             const NumericVector& piecewiseSurvivalTime = 0,
             const double kappa1 = 5,
             const double kappa2 = 5,
             const double lambda1 = 0.0875,
             const double lambda2 = 0.125,
             const NumericVector& gamma1 = 0,
             const NumericVector& gamma2 = 0,
             const double accrualDuration = 1.25,
             const double followupTime = 2.75,
             const bool fixedFollowup = 0,
             const NumericVector& spendingTime = NA_REAL,
             const double studyDuration = NA_REAL,
             const bool nullVariance = 0) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double bsfpar = parameterBetaSpending;

  int nintervals = piecewiseSurvivalTime.size();

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
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
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
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
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

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (rateRatioH0 <= 0) {
    stop("rateRatioH0 must be positive");
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

  if (kappa1 < 0) {
    stop("kappa1 must be nonnegative");
  }

  if (kappa2 < 0) {
    stop("kappa2 must be nonnegative");
  }

  if (lambda1 <= 0) {
    stop("lambda1 must be positive");
  }

  if (lambda2 <= 0) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals) {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals) {
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

  if (fixedFollowup && R_isnancpp(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
  }


  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));
  alpha1 = cumAlphaSpent[kMax - 1];

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
  List na;
  DataFrame nb, nc;
  NumericVector I(kMax), time(kMax);

  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || R_isnancpp(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }
  u0[0] = studyDuration1;

  // obtain the timing of interim analysis
  na = nbstat(u0, 1, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup, 0);
  nb = as<DataFrame>(na["resultsUnderH1"]);
  double maxInformation = sum(NumericVector(nb[18]));
  I = maxInformation*informationRates1;

  double information = 0;

  // Lambda function
  auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, kappa1, kappa2,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information](double t)->double {
              NumericVector u0(1, t);
              List na = nbstat(
                u0, 1, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, kappa1, kappa2,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup, 0);
              DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
              return sum(NumericVector(nb[18])) - information;
            };

  for (int i=0; i<kMax-1; i++) {
    // match the predicted information to the target
    information = std::max(I[i], 0.0);
    time[i] = brent(f, 0.0001, studyDuration1, 0.0001);
  }
  time[kMax-1] = studyDuration1;

  // obtain the variance ratio
  na = nbstat(time, rateRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              nullVariance);
  nb = as<DataFrame>(na["resultsUnderH1"]);
  nc = as<DataFrame>(na["resultsUnderH0"]);
  NumericVector varianceRatio = as<NumericVector>(nc["varianceRatio"]);
  NumericVector w = sqrt(varianceRatio);

  NumericVector rateRatio = nb[16];
  NumericVector theta(kMax);

  bool directionUpper = lambda1/lambda2 > rateRatioH0;
  if (directionUpper) {
    theta = log(rateRatio/rateRatioH0);
  } else {
    theta = -log(rateRatio/rateRatioH0);
  }

  double r1 = allocationRatioPlanned/(1+allocationRatioPlanned);

  NumericVector nsubjects = NumericVector(nb[1]);
  NumericVector nsubjects1 = r1*nsubjects;
  NumericVector nsubjects2 = (1-r1)*nsubjects;
  NumericVector nevents = NumericVector(nb[2]);
  NumericVector nevents1 = NumericVector(nb[3]);
  NumericVector nevents2 = NumericVector(nb[4]);
  NumericVector ndropouts = NumericVector(nb[5]);
  NumericVector ndropouts1 = NumericVector(nb[6]);
  NumericVector ndropouts2 = NumericVector(nb[7]);
  NumericVector exposure = NumericVector(nb[11]);
  NumericVector exposure1 = NumericVector(nb[12]);
  NumericVector exposure2 = NumericVector(nb[13]);


  // compute the stagewise exit probabilities for efficacy and futility
  if (!missingFutilityBounds || bsf=="none" || kMax==1) {
    probs = exitprobcpp(criticalValues1*w, futilityBounds1*w, theta, I);
  } else {
    List out = getPower(alpha1, kMax, criticalValues1, theta, I,
                        bsf, bsfpar, spendingTime1, futilityStopping1, w);
    futilityBounds1 = out[1];
    probs = out[2];
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
  double expectedExposure = sum(ptotal*exposure);
  double expectedNumberOfEvents1 = sum(ptotal*nevents1);
  double expectedNumberOfDropouts1 = sum(ptotal*ndropouts1);
  double expectedNumberOfSubjects1 = sum(ptotal*nsubjects1);
  double expectedExposure1 = sum(ptotal*exposure1);
  double expectedNumberOfEvents2 = sum(ptotal*nevents2);
  double expectedNumberOfDropouts2 = sum(ptotal*ndropouts2);
  double expectedNumberOfSubjects2 = sum(ptotal*nsubjects2);
  double expectedExposure2 = sum(ptotal*exposure2);
  double expectedStudyDuration = sum(ptotal*time);
  double expectedInformation = sum(ptotal*I);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  NumericVector rru(kMax), rrl(kMax);
  if (directionUpper) {
    rru = rateRatioH0*exp(criticalValues1/sqrt(I)*w);
    rrl = rateRatioH0*exp(futilityBounds1/sqrt(I)*w);
  } else {
    rru = rateRatioH0*exp(-criticalValues1/sqrt(I)*w);
    rrl = rateRatioH0*exp(-futilityBounds1/sqrt(I)*w);
  }

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      rru[i] = NA_REAL;
      efficacyStopping1[i] = 0;
    }

    if (futilityBounds1[i] == -6) {
      rrl[i] = NA_REAL;
      futilityStopping1[i] = 0;
    }
  }


  DataFrame byStageResults = DataFrame::create(
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
    _["exposure"] = exposure,
    _["analysisTime"] = time,
    _["efficacyRateRatio"] = rru,
    _["futilityRateRatio"] = rrl,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["information"] = I,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1);


  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfDropouts"] = (ndropouts[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["exposure"] = (exposure[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfDropouts"] = expectedNumberOfDropouts,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedExposure"] = expectedExposure,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["kMax"] = kMax,
    _["rateRatioH0"] = rateRatioH0,
    _["rateRatio"] = lambda1/lambda2);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["typeBetaSpending"] = typeBetaSpending,
    _["parameterBetaSpending"] = parameterBetaSpending,
    _["allocationRatioPlanned"] = allocationRatioPlanned,
    _["accrualTime"] = accrualTime,
    _["accrualIntensity"] = accrualIntensity,
    _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
    _["kappa1"] = kappa1,
    _["kappa2"] = kappa2,
    _["lambda1"] = lambda1,
    _["lambda2"] = lambda2,
    _["gamma1"] = gamma1,
    _["gamma2"] = gamma2,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["spendingTime"] = spendingTime,
    _["nullVariance"] = nullVariance);

  List byTreatmentCounts = List::create(
    _["numberOfEvents1"] = nevents1,
    _["numberOfDropouts1"] = ndropouts1,
    _["numberOfSubjects1"] = nsubjects1,
    _["exposure1"] = exposure1,
    _["numberOfEvents2"] = nevents2,
    _["numberOfDropouts2"] = ndropouts2,
    _["numberOfSubjects2"] = nsubjects2,
    _["exposure2"] = exposure2,
    _["expectedNumberOfEvents1"] = expectedNumberOfEvents1,
    _["expectedNumberOfDropouts1"] = expectedNumberOfDropouts1,
    _["expectedNumberOfSubjects1"] = expectedNumberOfSubjects1,
    _["expectedExposure1"] = expectedExposure1,
    _["expectedNumberOfEvents2"] = expectedNumberOfEvents2,
    _["expectedNumberOfDropouts2"] = expectedNumberOfDropouts2,
    _["expectedNumberOfSubjects2"] = expectedNumberOfSubjects2,
    _["expectedExposure2"] = expectedExposure2);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings,
    _["byTreatmentCounts"] = byTreatmentCounts);

  result.attr("class") = "nbpower";

  return result;
}


//' @title Negative binomial rate ratio sample size
//' @description Obtains the needed accrual duration given power and
//' follow-up time, the needed follow-up time given power and
//' accrual duration, or the needed absolute accrual rates given
//' power, accrual duration, follow-up duration, and relative accrual
//' rates in a two-group negative binomial design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
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
//' @param rateRatioH0 Rate ratio under the null hypothesis.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @param kappa1 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the active treatment group.
//' @param kappa2 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the control group.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group.
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}. Adjustment
//'   may be needed for non-monotone relationship with study power.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//' @param nullVariance Whether to calculate the variance for log rate ratio
//'   under the null hypothesis.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{nbpower} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{nbpower} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{nbpower}}
//'
//' @examples
//' # Example 1: Obtains follow-up duration given power, accrual intensity,
//' # and accrual duration for variable follow-up
//'
//' nbsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.5, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              accrualIntensity = 1956/1.25,
//'              kappa1 = 5, kappa2 = 5,
//'              lambda1 = 0.0875, lambda2 = 0.125,
//'              gamma1 = 0, gamma2 = 0,
//'              accrualDuration = 1.25,
//'              followupTime = NA, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains accrual intensity given power, accrual duration, and
//' # follow-up duration for variable follow-up
//'
//' nbsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.5, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              kappa1 = 5, kappa2 = 5,
//'              lambda1 = 0.0875, lambda2 = 0.125,
//'              gamma1 = 0, gamma2 = 0,
//'              accrualDuration = 1.25,
//'              followupTime = 2.25, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains accrual duration given power, accrual intensity, and
//' # follow-up duration for fixed follow-up
//'
//' nbsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.5, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              accrualIntensity = 1667,
//'              kappa1 = 5, kappa2 = 5,
//'              lambda1 = 0.0875, lambda2 = 0.125,
//'              gamma1 = 0, gamma2 = 0,
//'              accrualDuration = NA,
//'              followupTime = 0.5, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List nbsamplesize(const double beta = 0.2,
                  const int kMax = 1,
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
                  const double rateRatioH0 = 1,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = 1500,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const double kappa1 = 5,
                  const double kappa2 = 5,
                  const double lambda1 = 0.0875,
                  const double lambda2 = 0.125,
                  const NumericVector& gamma1 = 0,
                  const NumericVector& gamma2 = 0,
                  double accrualDuration = NA_REAL,
                  double followupTime = NA_REAL,
                  const bool fixedFollowup = 0,
                  const NumericVector& interval =
                    NumericVector::create(0.001, 240),
                    const NumericVector& spendingTime = NA_REAL,
                    const bool rounding = 1,
                    const bool nullVariance = 0) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double bsfpar = parameterBetaSpending;

  int nintervals = piecewiseSurvivalTime.size();

  if (R_isnancpp(beta)) {
    stop("beta must be provided");
  }

  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
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
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
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
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
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

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="user" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(futilityBounds))) && bsf=="user") {
    if (is_true(any(is_na(userBetaSpending)))) {
      stop("userBetaSpending must be specified");
    } else if (userBetaSpending.size() < kMax) {
      stop("Insufficient length of userBetaSpending");
    } else if (userBetaSpending[0] < 0) {
      stop("Elements of userBetaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userBetaSpending) < 0))) {
      stop("Elements of userBetaSpending must be nondecreasing");
    } else if (userBetaSpending[kMax-1] != beta) {
      stop("userBetaSpending must end with specified beta");
    }
  }

  if (rateRatioH0 <= 0) {
    stop("rateRatioH0 must be positive");
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

  if (kappa1 < 0) {
    stop("kappa1 must be nonnegative");
  }

  if (kappa2 < 0) {
    stop("kappa2 must be nonnegative");
  }

  if (lambda1 <= 0) {
    stop("lambda1 must be positive");
  }

  if (lambda2 <= 0) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals) {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals) {
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

  if (fixedFollowup && R_isnancpp(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
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

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  alpha1 = sum(NumericVector(probs[0]));

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


  double theta1 = fabs(log(lambda1/lambda2) - log(rateRatioH0));
  NumericVector theta(kMax, theta1);

  double maxInformation;
  if (!nullVariance) {
    List design = getDesign(
      beta, NA_REAL, theta1, kMax, informationRates1,
      efficacyStopping1, futilityStopping1, criticalValues1,
      alpha1, asf, asfpar, userAlphaSpending, futilityBounds1,
      bsf, bsfpar, userBetaSpending, spendingTime1, 1);

    DataFrame byStageResults = as<DataFrame>(design["byStageResults"]);
    criticalValues1 = byStageResults["efficacyBounds"];

    DataFrame overallResults = as<DataFrame>(design["overallResults"]);
    maxInformation = overallResults["information"];

    auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              unknown, maxInformation](double aval)-> double{
                NumericVector accrualIntensity1 = clone(accrualIntensity);
                double dur1=0, dur2=0;

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

                // obtain the maximum information at study end
                NumericVector u0(1, dur1 + dur2);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda1, lambda2, gamma1, gamma2,
                  dur1, dur2, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
              };

    if (unknown == "accrualDuration") {
      accrualDuration = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 0.0001);
      accrualIntensity1 = aval*accrualIntensity;
    }
  } else {
    auto f = [beta, kMax, informationRates1,
              futilityStopping1, criticalValues1, alpha1,
              futilityBounds1, bsf, bsfpar,
              rateRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              spendingTime1, nullVariance, missingFutilityBounds,
              theta, unknown](double aval)-> double{
                NumericVector accrualIntensity1 = clone(accrualIntensity);
                double dur1=0, dur2=0;

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

                double studyDuration = dur1 + dur2;

                // obtain the timing of interim analysis
                NumericVector u0(1, studyDuration), time(kMax);
                List na = nbstat(u0, 1, allocationRatioPlanned,
                                 accrualTime, accrualIntensity1,
                                 piecewiseSurvivalTime, kappa1, kappa2,
                                 lambda1, lambda2, gamma1, gamma2,
                                 dur1, dur2, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                double maxInformation = sum(NumericVector(nb[18]));
                NumericVector I = maxInformation*informationRates1;

                double information = 0;

                // Lambda function
                auto g = [allocationRatioPlanned, accrualTime,
                          accrualIntensity1, piecewiseSurvivalTime,
                          kappa1, kappa2, lambda1, lambda2,
                          gamma1, gamma2, dur1, dur2, fixedFollowup,
                          &information](double t)->double {
                            NumericVector u0(1,t);
                            List na = nbstat(u0, 1, allocationRatioPlanned,
                                             accrualTime, accrualIntensity1,
                                             piecewiseSurvivalTime, kappa1,
                                             kappa2, lambda1, lambda2,
                                             gamma1, gamma2, dur1, dur2,
                                             fixedFollowup, 0);
                            DataFrame nb =
                              as<DataFrame>(na["resultsUnderH1"]);
                            return sum(NumericVector(nb[18])) - information;
                          };

                for (int i=0; i<kMax-1; i++) {
                  // match the predicted information to the target
                  information = std::max(I[i], 0.0);
                  time[i] = brent(g, 0.0001, studyDuration, 0.0001);
                }
                time[kMax-1] = studyDuration;

                // obtain the variance ratio
                na = nbstat(time, rateRatioH0, allocationRatioPlanned,
                            accrualTime, accrualIntensity1,
                            piecewiseSurvivalTime, kappa1, kappa2,
                            lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup, nullVariance);
                DataFrame nc = as<DataFrame>(na["resultsUnderH0"]);
                NumericVector varianceRatio =
                  as<NumericVector>(nc["varianceRatio"]);
                NumericVector w = sqrt(varianceRatio);

                // compute the stagewise exit probabilities
                List probs;
                if (!missingFutilityBounds || bsf=="none" || kMax==1) {
                  probs = exitprobcpp(criticalValues1*w, futilityBounds1*w,
                                      theta, I);
                } else {
                  List out = getPower(alpha1, kMax, criticalValues1, theta,
                                      I, bsf, bsfpar, spendingTime1,
                                      futilityStopping1, w);
                  probs = out[2];
                }

                return sum(NumericVector(probs[0])) - (1-beta);
              };

    if (unknown == "accrualDuration") {
      accrualDuration = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 0.0001);
      accrualIntensity1 = aval*accrualIntensity;
    }

    NumericVector u0(1, accrualDuration + followupTime);
    List na = nbstat(u0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity1,
                     piecewiseSurvivalTime, kappa1, kappa2,
                     lambda1, lambda2, gamma1, gamma2,
                     accrualDuration, followupTime, fixedFollowup, 0);
    DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
    maxInformation = sum(NumericVector(nb[18]));
  }

  double studyDuration = accrualDuration + followupTime;


  // output the results
  List resultH1, resultH0, result;

  if (rounding) {
    NumericVector u0(1, studyDuration);
    double n0 = accrual(u0, accrualTime, accrualIntensity1,
                        accrualDuration)[0];
    double n = std::ceil(n0);

    if (n - n0 > 1e-4) {
      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        double aval = n/n0;
        accrualIntensity1 = aval*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      if (!fixedFollowup) {
        // adjust follow-up time to obtain the target maximum information
        auto h = [allocationRatioPlanned, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval*followupTime);
                    List na = nbstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, kappa1, kappa2,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval*followupTime, fixedFollowup, 0);
                    DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[18])) - maxInformation;
                  };

        double aval = brent(h, 0, 1.1, 1e-6);
        followupTime = aval*followupTime;
        studyDuration = accrualDuration + followupTime;
      } else {
        // adjust study duration to obtain the target maximum information
        auto h = [allocationRatioPlanned, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval*followupTime);
                    List na = nbstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, kappa1, kappa2,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup, 0);
                    DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[18])) - maxInformation;
                  };

        double aval = brent(h, 0, 1, 1e-6);
        studyDuration = accrualDuration + aval*followupTime;
      }
    }
  }

  resultH1 = nbpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending, rateRatioH0,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, kappa1, kappa2,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration, nullVariance);


  // obtain results under H0 by matching the maximum information
  double aval;
  if (!fixedFollowup) {
    auto h = [rateRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval*followupTime);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda2*rateRatioH0, lambda2, gamma1, gamma2,
                  accrualDuration, aval*followupTime, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
              };

    if (h(0) < 0) { // adjust the follow-up time
      aval = brent(h, 0, 2, 1e-6);
      followupTime = aval*followupTime;
      studyDuration = accrualDuration + followupTime;
    } else { // adjust the accrual duration
      auto g = [rateRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, kappa1, kappa2,
                lambda2, gamma1, gamma2,
                accrualDuration, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval*accrualDuration);
                  List na = nbstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, kappa1, kappa2,
                    lambda2*rateRatioH0, lambda2, gamma1, gamma2,
                    aval*accrualDuration, 0, fixedFollowup, 0);
                  DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[18])) - maxInformation;
                };

      aval = brent(g, 0.1, 1, 1e-6);
      accrualDuration = aval*accrualDuration;
      followupTime = 0;
      studyDuration = accrualDuration + followupTime;
    }
  } else { // fixed follow-up
    auto h = [rateRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval*followupTime);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda2*rateRatioH0, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
              };

    if (h(1) < 0) { // increase accrual duration
      auto g = [rateRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, kappa1, kappa2,
                lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval*accrualDuration + followupTime);
                  List na = nbstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, kappa1, kappa2,
                    lambda2*rateRatioH0, lambda2, gamma1, gamma2,
                    aval*accrualDuration, followupTime, fixedFollowup, 0);
                  DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[18])) - maxInformation;
                };
      double aval = brent(g, 1, 2, 1e-6);
      accrualDuration = aval*accrualDuration;
      studyDuration = accrualDuration + followupTime;
    } else if (h(0) < 0) { // shorten study duration
      double aval = brent(h, 0, 1, 1e-6);
      studyDuration = accrualDuration + aval*followupTime;
    } else { // decrease accrual duration
      auto g = [rateRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, kappa1, kappa2,
                lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval*accrualDuration);
                  List na = nbstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, kappa1, kappa2,
                    lambda2*rateRatioH0, lambda2, gamma1, gamma2,
                    aval*accrualDuration, followupTime, fixedFollowup, 0);
                  DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[18])) - maxInformation;
                };

      double aval = brent(g, 0.1, 1, 1e-6);
      accrualDuration = aval*accrualDuration;
      studyDuration = accrualDuration;
    }
  }

  // use the same stopping boundaries as under H1
  resultH0 = nbpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending, rateRatioH0,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, kappa1, kappa2,
    lambda2*rateRatioH0, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration, 0);

  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}


//' @title One-sample negative binomial rate power
//' @description Estimates the power, stopping probabilities, and expected
//' sample size in a one-group negative binomial design.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @param typeBetaSpending The type of beta spending. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock
//'   type spending function, "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
//'   early futility stopping. Defaults to "none".
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @param kappa The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) of the negative binomial
//'   distribution.
//' @param lambdaH0 The rate parameter of the negative binomial distribution
//'   under the null hypothesis.
//' @param lambda The rate parameter of the negative binomial distribution
//'   under the alternative hypothesis.
//' @param gamma The hazard rate for exponential dropout or a vector of
//'   hazard rates for piecewise exponential dropout. Defaults to 0 for
//'   no dropout.
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{nbpower1s} object with 3 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numberOfDropouts}: The total number of dropouts.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{exposure}: The total exposure.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfDropouts}: The expected number of dropouts.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedExposure}: The expected exposure.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{kMax}: The number of stages.
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
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{exposure}: The exposure.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRate}: The efficacy boundaries on the rate scale.
//'
//'     - \code{futilityRate}: The futility boundaries on the rate scale.
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
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{typeBetaSpending},
//'   \code{parameterBetaSpending}, \code{accrualTime},
//'   \code{accuralIntensity}, \code{piecewiseSurvivalTime},
//'   \code{kappa}, \code{lambdaH0}, \code{lambda}, \code{gamma},
//'   \code{accrualDuration}, \code{followupTime}, \code{fixedFollowup},
//'   and \code{spendingTime}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{nbstat}}
//'
//' @examples
//' # Example 1: Variable follow-up design
//'
//' nbpower1s(kMax = 2, informationRates = c(0.5, 1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           accrualIntensity = 500,
//'           kappa = 5, lambdaH0 = 0.125, lambda = 0.0875,
//'           gamma = 0, accrualDuration = 1.25,
//'           followupTime = 2.75, fixedFollowup = FALSE)
//'
//' # Example 2: Fixed follow-up design
//'
//' nbpower1s(kMax = 2, informationRates = c(0.5, 1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           accrualIntensity = 40,
//'           kappa = 3, lambdaH0 = 8.4, lambda = 0.5*8.4,
//'           gamma = 0, accrualDuration = 1.5,
//'           followupTime = 0.5, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List nbpower1s(const int kMax = 1,
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
               const NumericVector& accrualTime = 0,
               const NumericVector& accrualIntensity = 1500,
               const NumericVector& piecewiseSurvivalTime = 0,
               const double kappa = 5,
               const double lambdaH0 = 0.125,
               const double lambda = 0.0875,
               const NumericVector& gamma = 0,
               const double accrualDuration = 1.25,
               const double followupTime = 2.75,
               const bool fixedFollowup = 0,
               const NumericVector& spendingTime = NA_REAL,
               const double studyDuration = NA_REAL) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double bsfpar = parameterBetaSpending;

  int nintervals = piecewiseSurvivalTime.size();

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
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
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
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
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

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
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

  if (kappa < 0) {
    stop("kappa must be nonnegative");
  }

  if (lambdaH0 <= 0) {
    stop("lambdaH0 must be positive");
  }

  if (lambda <= 0) {
    stop("lambda must be positive");
  }

  if (is_true(any(gamma < 0))) {
    stop("gamma must be non-negative");
  }

  if (gamma.size() != 1 && gamma.size() != nintervals) {
    stop("Invalid length for gamma");
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

  if (fixedFollowup && R_isnancpp(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
  }


  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));
  alpha1 = cumAlphaSpent[kMax - 1];

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
  List na;
  DataFrame nb;
  NumericVector I(kMax), time(kMax);

  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || R_isnancpp(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }
  u0[0] = studyDuration1;

  // obtain the timing of interim analysis using the twin treatment group
  na = nbstat(u0, 1, 1, accrualTime, 2*accrualIntensity,
              piecewiseSurvivalTime, kappa, kappa,
              lambda, lambda, gamma, gamma,
              accrualDuration, followupTime, fixedFollowup);
  nb = as<DataFrame>(na["resultsUnderH1"]);
  double maxInformation = sum(NumericVector(nb[14]));
  I = maxInformation*informationRates1;

  double information = 0;

  // Lambda function
  auto f = [accrualTime, accrualIntensity,
            piecewiseSurvivalTime, kappa, lambda, gamma,
            accrualDuration, followupTime, fixedFollowup,
            &information](double t)->double {
              NumericVector u0(1, t);
              List na = nbstat(
                u0, 1, 1, accrualTime, 2*accrualIntensity,
                piecewiseSurvivalTime, kappa, kappa,
                lambda, lambda, gamma, gamma,
                accrualDuration, followupTime, fixedFollowup, 0);
              DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
              return sum(NumericVector(nb[14])) - information;
            };

  for (int i=0; i<kMax; i++) {
    // match the predicted information to the target
    information = std::max(I[i], 0.0);
    time[i] = brent(f, 0.0001, studyDuration1, 0.0001);
  }

  na = nbstat(time, 1, 1, accrualTime, 2*accrualIntensity,
              piecewiseSurvivalTime, kappa, kappa,
              lambda, lambda, gamma, gamma,
              accrualDuration, followupTime, fixedFollowup, 0);
  nb = as<DataFrame>(na["resultsUnderH1"]);

  bool directionUpper = lambda > lambdaH0;
  double theta1 = fabs(log(lambda/lambdaH0));
  NumericVector theta = rep(theta1, kMax);

  NumericVector nsubjects = NumericVector(nb[1])/2.0;
  NumericVector nevents = NumericVector(nb[3]);
  NumericVector ndropouts = NumericVector(nb[6]);
  NumericVector exposure = NumericVector(nb[12]);

  // compute the stagewise exit probabilities for efficacy and futility
  if (!missingFutilityBounds || bsf=="none" || kMax==1) {
    probs = exitprobcpp(criticalValues1, futilityBounds1, theta, I);
  } else {
    NumericVector w(kMax, 1.0);
    List out = getPower(alpha1, kMax, criticalValues1, theta, I,
                        bsf, bsfpar, spendingTime1, futilityStopping1, w);
    futilityBounds1 = out[1];
    probs = out[2];
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
  double expectedExposure = sum(ptotal*exposure);
  double expectedStudyDuration = sum(ptotal*time);
  double expectedInformation = sum(ptotal*I);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  NumericVector rateu(kMax), ratel(kMax);

  if (directionUpper) {
    rateu = lambdaH0*exp(criticalValues1/sqrt(I));
    ratel = lambdaH0*exp(futilityBounds1/sqrt(I));
  } else {
    rateu = lambdaH0*exp(-criticalValues1/sqrt(I));
    ratel = lambdaH0*exp(-futilityBounds1/sqrt(I));
  }

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      rateu[i] = NA_REAL;
      efficacyStopping1[i] = 0;
    }

    if (futilityBounds1[i] == -6) {
      ratel[i] = NA_REAL;
      futilityStopping1[i] = 0;
    }
  }


  DataFrame byStageResults = DataFrame::create(
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
    _["exposure"] = exposure,
    _["analysisTime"] = time,
    _["efficacyRate"] = rateu,
    _["futilityRate"] = ratel,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["information"] = I,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1);


  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfDropouts"] = (ndropouts[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["exposure"] = (exposure[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfDropouts"] = expectedNumberOfDropouts,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedExposure"] = expectedExposure,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["kMax"] = kMax);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["typeBetaSpending"] = typeBetaSpending,
    _["parameterBetaSpending"] = parameterBetaSpending,
    _["accrualTime"] = accrualTime,
    _["accrualIntensity"] = accrualIntensity,
    _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
    _["kappa"] = kappa,
    _["lambdaH0"] = lambdaH0,
    _["lambda"] = lambda,
    _["gamma"] = gamma,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["spendingTime"] = spendingTime);


  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "nbpower1s";

  return result;
}


//' @title One-sample negative binomial rate sample size
//' @description Obtains the needed accrual duration given power and
//' follow-up time, the needed follow-up time given power and
//' accrual duration, or the needed absolute accrual rates given
//' power, accrual duration, follow-up duration, and relative accrual
//' rates in a one-group negative binomial design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
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
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @param kappa The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) of the negative binomial
//'   distribution.
//' @param lambdaH0 The rate parameter of the negative binomial distribution
//'   under the null hypothesis.
//' @param lambda The rate parameter of the negative binomial distribution
//'   under the alternative hypothesis.
//' @param gamma The hazard rate for exponential dropout or a vector of
//'   hazard rates for piecewise exponential dropout. Defaults to 0 for
//'   no dropout.
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}. Adjustment
//'   may be needed for non-monotone relationship with study power.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{nbpower1s} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{nbpower1s} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{nbpower1s}}
//'
//' @examples
//' # Example 1: Obtains follow-up duration given power, accrual intensity,
//' # and accrual duration for variable follow-up
//'
//' nbsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.5, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                accrualIntensity = 500,
//'                kappa = 5, lambdaH0 = 0.125, lambda = 0.0875,
//'                gamma = 0, accrualDuration = 1.25,
//'                followupTime = NA, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains accrual intensity given power, accrual duration, and
//' # follow-up duration for variable follow-up
//'
//' nbsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.5, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                kappa = 5, lambdaH0 = 0.125, lambda = 0.0875,
//'                gamma = 0, accrualDuration = 1.25,
//'                followupTime = 2.25, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains accrual duration given power, accrual intensity, and
//' # follow-up duration for fixed follow-up
//'
//' nbsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.5, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                accrualIntensity = 40,
//'                kappa = 3, lambdaH0 = 8.4, lambda = 4.2,
//'                gamma = 0, accrualDuration = NA,
//'                followupTime = 0.5, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List nbsamplesize1s(const double beta = 0.2,
                    const int kMax = 1,
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
                    const NumericVector& accrualTime = 0,
                    const NumericVector& accrualIntensity = 1500,
                    const NumericVector& piecewiseSurvivalTime = 0,
                    const double kappa = 5,
                    const double lambdaH0 = 0.125,
                    const double lambda = 0.0875,
                    const NumericVector& gamma = 0,
                    double accrualDuration = NA_REAL,
                    double followupTime = NA_REAL,
                    const bool fixedFollowup = 0,
                    const NumericVector& interval =
                    NumericVector::create(0.001, 240),
                    const NumericVector& spendingTime = NA_REAL,
                    const bool rounding = 1) {

  double alpha1 = alpha;
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double bsfpar = parameterBetaSpending;

  int nintervals = piecewiseSurvivalTime.size();

  if (R_isnancpp(beta)) {
    stop("beta must be provided");
  }

  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
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
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
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
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
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

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="user" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(futilityBounds))) && bsf=="user") {
    if (is_true(any(is_na(userBetaSpending)))) {
      stop("userBetaSpending must be specified");
    } else if (userBetaSpending.size() < kMax) {
      stop("Insufficient length of userBetaSpending");
    } else if (userBetaSpending[0] < 0) {
      stop("Elements of userBetaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userBetaSpending) < 0))) {
      stop("Elements of userBetaSpending must be nondecreasing");
    } else if (userBetaSpending[kMax-1] != beta) {
      stop("userBetaSpending must end with specified beta");
    }
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


  if (kappa < 0) {
    stop("kappa must be nonnegative");
  }

  if (lambdaH0 <= 0) {
    stop("lambdaH0 must be positive");
  }

  if (lambda <= 0) {
    stop("lambda must be positive");
  }

  if (is_true(any(gamma < 0))) {
    stop("gamma must be non-negative");
  }

  if (gamma.size() != 1 && gamma.size() != nintervals) {
    stop("Invalid length for gamma");
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

  if (fixedFollowup && R_isnancpp(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
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

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  alpha1 = sum(NumericVector(probs[0]));

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

  double theta1 = fabs(log(lambda/lambdaH0));

  List design = getDesign(
    beta, NA_REAL, theta1, kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, asf, asfpar, userAlphaSpending, futilityBounds1,
    bsf, bsfpar, userBetaSpending, spendingTime1, 1);

  DataFrame byStageResults = as<DataFrame>(design["byStageResults"]);
  criticalValues1 = byStageResults["efficacyBounds"];
  futilityBounds1 = byStageResults["futilityBounds"];

  DataFrame overallResults = as<DataFrame>(design["overallResults"]);
  double maxInformation = overallResults["information"];

  auto f = [accrualTime, accrualIntensity,
            piecewiseSurvivalTime, kappa, lambda, gamma,
            accrualDuration, followupTime, fixedFollowup,
            unknown, maxInformation](double aval)-> double{
              NumericVector accrualIntensity1 = clone(accrualIntensity);
              double dur1=0, dur2=0;

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

              // obtain the maximum information at study end
              NumericVector u0(1, dur1 + dur2);
              List na = nbstat(
                u0, 1, 1, accrualTime, 2*accrualIntensity1,
                piecewiseSurvivalTime, kappa, kappa,
                lambda, lambda, gamma, gamma,
                dur1, dur2, fixedFollowup, 0);
              DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
              return sum(NumericVector(nb[14])) - maxInformation;
            };

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
  List resultH1, resultH0, result;
  double studyDuration = accrualDuration + followupTime;

  if (rounding) {
    NumericVector u0(1, studyDuration);
    double n0 = accrual(u0, accrualTime, accrualIntensity1,
                        accrualDuration)[0];
    double n = std::ceil(n0);

    if (n - n0 > 1e-4) {
      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        double aval = n/n0;
        accrualIntensity1 = aval*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      if (!fixedFollowup) {
        // adjust follow-up time to obtain the target maximum information
        auto h = [accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa, lambda, gamma,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval*followupTime);
                    List na = nbstat(
                      u0, 1, 1, accrualTime, 2*accrualIntensity1,
                      piecewiseSurvivalTime, kappa, kappa,
                      lambda, lambda, gamma, gamma,
                      accrualDuration, aval*followupTime, fixedFollowup, 0);
                    DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[14])) - maxInformation;
                  };

        double aval = brent(h, 0, 2, 1e-6);
        followupTime = aval*followupTime;
        studyDuration = accrualDuration + followupTime;
      } else {
        // adjust study duration to obtain the target maximum information
        auto h = [accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa, lambda, gamma,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval*followupTime);
                    List na = nbstat(
                      u0, 1, 1, accrualTime, 2*accrualIntensity1,
                      piecewiseSurvivalTime, kappa, kappa,
                      lambda, lambda, gamma, gamma,
                      accrualDuration, followupTime, fixedFollowup, 0);
                    DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[14])) - maxInformation;
                  };

        double aval = brent(h, 0.1, 1, 1e-6);
        studyDuration = accrualDuration + aval*followupTime;
      }
    }
  }

  resultH1 = nbpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    accrualTime, accrualIntensity1, piecewiseSurvivalTime,
    kappa, lambdaH0, lambda, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  // obtain results under H0 by matching the maximum information
  overallResults = as<DataFrame>(resultH1["overallResults"]);
  byStageResults = as<DataFrame>(resultH1["byStageResults"]);
  NumericVector information = byStageResults["information"];
  maxInformation = information[kMax-1];
  double aval;

  if (!fixedFollowup) {
    auto h = [lambdaH0, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, kappa, gamma,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval*followupTime);
                List na = nbstat(
                  u0, 1, 1, accrualTime, 2*accrualIntensity1,
                  piecewiseSurvivalTime, kappa, kappa,
                  lambdaH0, lambdaH0, gamma, gamma,
                  accrualDuration, aval*followupTime, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[14])) - maxInformation;
              };

    if (h(0) < 0) { // adjust the follow-up time
      aval = brent(h, 0, 2, 1e-6);
      followupTime = aval*followupTime;
      studyDuration = accrualDuration + followupTime;
    } else { // adjust the accrual duration
      auto g = [lambdaH0, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, kappa, gamma,
                accrualDuration, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval*accrualDuration);
                  List na = nbstat(
                    u0, 1, 1, accrualTime, 2*accrualIntensity1,
                    piecewiseSurvivalTime, kappa, kappa,
                    lambdaH0, lambdaH0, gamma, gamma,
                    aval*accrualDuration, 0, fixedFollowup, 0);
                  DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[14])) - maxInformation;
                };

      aval = brent(g, 0.1, 1, 1e-6);
      accrualDuration = aval*accrualDuration;
      followupTime = 0;
      studyDuration = accrualDuration + followupTime;
    }
  } else { // fixed follow-up
    auto h = [lambdaH0, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, kappa, gamma,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval*followupTime);
                List na = nbstat(
                  u0, 1, 1, accrualTime, 2*accrualIntensity1,
                  piecewiseSurvivalTime, kappa, kappa,
                  lambdaH0, lambdaH0, gamma, gamma,
                  accrualDuration, followupTime, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[14])) - maxInformation;
              };

    if (h(1) < 0) { // increase the accrual duration
      auto g = [lambdaH0, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, kappa, gamma,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval*accrualDuration + followupTime);
                  List na = nbstat(
                    u0, 1, 1, accrualTime, 2*accrualIntensity1,
                    piecewiseSurvivalTime, kappa, kappa,
                    lambdaH0, lambdaH0, gamma, gamma,
                    aval*accrualDuration, followupTime, fixedFollowup, 0);
                  DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[14])) - maxInformation;
                };

      aval = brent(g, 1, 2, 1e-6);
      accrualDuration = aval*accrualDuration;
      studyDuration = accrualDuration + followupTime;
    } else if (h(0) < 0) { // shorten the study duration
      aval = brent(h, 0, 1, 1e-6);
      studyDuration = accrualDuration + aval*followupTime;
    } else { // adjust the accrual duration
      auto g = [lambdaH0, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, kappa, gamma,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval*accrualDuration);
                  List na = nbstat(
                    u0, 1, 1, accrualTime, 2*accrualIntensity1,
                    piecewiseSurvivalTime, kappa, kappa,
                    lambdaH0, lambdaH0, gamma, gamma,
                    aval*accrualDuration, followupTime, fixedFollowup);
                  DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[14])) - maxInformation;
                };

      aval = brent(g, 0.1, 1, 1e-6);
      accrualDuration = aval*accrualDuration;
      studyDuration = accrualDuration;
    }
  }

  // use the same stopping boundaries as under H1
  criticalValues1 = byStageResults["efficacyBounds"];
  futilityBounds1 = byStageResults["futilityBounds"];

  resultH0 = nbpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    accrualTime, accrualIntensity1, piecewiseSurvivalTime,
    kappa, lambdaH0, lambdaH0, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}


//' @title Power for equivalence in negative binomial rate ratio
//' @description Obtains the power for equivalence in negative binomial
//' rate ratio.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_criticalValues
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param rateRatioLower The lower equivalence limit of rate ratio.
//' @param rateRatioUpper The upper equivalence limit of rate ratio.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @param kappa1 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the active treatment group.
//' @param kappa2 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the control group.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group.
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//' @param nullVariance Whether to calculate the variance for log rate ratio
//'   under the null hypothesis.
//'
//' @return An S3 class \code{nbpowerequiv} object with 4 components:
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
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numberOfDropouts}: The total number of dropouts.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{exposure}: The total exposure.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfDropouts}: The expected number of dropouts.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedExposure}: The expected exposure.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{rateRatioLower}: The lower equivalence limit of rate ratio.
//'
//'     - \code{rateRatioUpper}: The upper equivalence limit of rate ratio.
//'
//'     - \code{rateRatio}: The rate ratio.
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
//'     - \code{cumulativeAttainedAlphaH10}: The cumulative alpha attained
//'       under \code{H10}.
//'
//'     - \code{cumulativeAttainedAlphaH20}: The cumulative alpha attained
//'       under \code{H20}.
//'
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{exposure}: The exposure.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRateRatioLower}: The efficacy boundaries on the
//'       rate ratio scale for the one-sided null hypothesis at the
//'       lower equivalence limit.
//'
//'     - \code{efficacyRateRatioUpper}: The efficacy boundaries on the
//'       rate ratio scale for the one-sided null hypothesis at the
//'       upper equivalence limit.
//'
//'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
//'       each of the two one-sided tests.
//'
//'     - \code{information}: The cumulative information.
//'
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{allocationRatioPlanned},
//'   \code{accrualTime}, \code{accuralIntensity},
//'   \code{piecewiseSurvivalTime}, \code{kappa1}, \code{kappa2},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   \code{accrualDuration}, \code{followupTime}, \code{fixedFollowup},
//'   \code{spendingTime}, \code{nullVariance}, and \code{varianceRatios}.
//'   The \code{varianceRatios} is a data frame with the following
//'   variables:
//'
//'     - \code{varianceRatioH10}: The ratio of the variance under
//'       \code{H10} to the variance under \code{H1}.
//'
//'     - \code{varianceRatioH20}: The ratio of the variance under
//'       \code{H20} to the variance under \code{H1}.
//'
//'     - \code{varianceRatioH12}: The ratio of the variance under
//'       \code{H10} to the variance under \code{H20}.
//'
//'     - \code{varianceRatioH21}: The ratio of the variance under
//'       \code{H20} to the variance under \code{H10}.//'
//'
//' * \code{byTreatmentCounts}: A list containing the following counts by
//'   treatment group:
//'
//'     - \code{numberOfEvents1}: The number of events by stage for
//'       the treatment group.
//'
//'     - \code{numberOfDropouts1}: The number of dropouts by stage for
//'       the treatment group.
//'
//'     - \code{numberOfSubjects1}: The number of subjects by stage for
//'       the treatment group.
//'
//'     - \code{exposure1}: The exposure by stage for the treatment group.
//'
//'     - \code{numberOfEvents2}: The number of events by stage for
//'       the control group.
//'
//'     - \code{numberOfDropouts2}: The number of dropouts by stage for
//'       the control group.
//'
//'     - \code{numberOfSubjects2}: The number of subjects by stage for
//'       the control group.
//'
//'     - \code{exposure2}: The exposure by stage for the control group.
//'
//'     - \code{expectedNumberOfEvents1}: The expected number of events for
//'       the treatment group.
//'
//'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
//'       for the treatment group.
//'
//'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
//'       for the treatment group.
//'
//'     - \code{expectedExposure1}: The expected exposure for the treatment
//'       group.
//'
//'     - \code{expectedNumberOfEvents2}: The expected number of events for
//'       control group.
//'
//'     - \code{expectedNumberOfDropouts2}: The expected number of dropouts
//'       for the control group.
//'
//'     - \code{expectedNumberOfSubjects2}: The expected number of subjects
//'       for the control group.
//'
//'     - \code{expectedExposure2}: The expected exposure for the control
//'       group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{nbstat}}
//'
//' @examples
//'
//' # Example 1: Variable follow-up design
//' nbpowerequiv(kMax = 2, informationRates = c(0.5, 1),
//'              alpha = 0.05, typeAlphaSpending = "sfOF",
//'              rateRatioLower = 2/3, rateRatioUpper = 3/2,
//'              accrualIntensity = 1956/1.25,
//'              kappa1 = 5, kappa2 = 5,
//'              lambda1 = 0.125, lambda2 = 0.125,
//'              gamma1 = 0, gamma2 = 0,
//'              accrualDuration = 1.25,
//'              followupTime = 2.75, fixedFollowup = FALSE,
//'              nullVariance = 1)
//'
//' # Example 2: Fixed follow-up design
//' nbpowerequiv(kMax = 2, informationRates = c(0.5, 1),
//'              alpha = 0.05, typeAlphaSpending = "sfOF",
//'              rateRatioLower = 0.5, rateRatioUpper = 2,
//'              accrualIntensity = 220/1.5,
//'              kappa1 = 3, kappa2 = 3,
//'              lambda1 = 8.4, lambda2 = 8.4,
//'              gamma1 = 0, gamma2 = 0,
//'              accrualDuration = 1.5,
//'              followupTime = 0.5, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List nbpowerequiv(const int kMax = 1,
                  const NumericVector& informationRates = NA_REAL,
                  const NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.05,
                  const String typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const double rateRatioLower = NA_REAL,
                  const double rateRatioUpper = NA_REAL,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = 1500,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const double kappa1 = 5,
                  const double kappa2 = 5,
                  const double lambda1 = 0.125,
                  const double lambda2 = 0.125,
                  const NumericVector& gamma1 = 0,
                  const NumericVector& gamma2 = 0,
                  const double accrualDuration = 1.25,
                  const double followupTime = 2.75,
                  const bool fixedFollowup = 0,
                  const NumericVector& spendingTime = NA_REAL,
                  const double studyDuration = NA_REAL,
                  const bool nullVariance = 0) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  int nintervals = piecewiseSurvivalTime.size();

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


  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (R_isnancpp(rateRatioLower)) {
    stop("rateRatioLower must be provided");
  }

  if (R_isnancpp(rateRatioUpper)) {
    stop("rateRatioUpper must be provided");
  }

  if (rateRatioLower <= 0) {
    stop("rateRatioLower must be positive");
  }

  double rateRatio = lambda1/lambda2;

  if (rateRatioLower >= rateRatio) {
    stop("rateRatioLower must be less than lambda1/lambda2");
  }

  if (rateRatioUpper <= rateRatio) {
    stop("rateRatioUpper must be greater than lambda1/lambda2");
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

  if (kappa1 < 0) {
    stop("kappa1 must be nonnegative");
  }

  if (kappa2 < 0) {
    stop("kappa2 must be nonnegative");
  }

  if (lambda1 <= 0) {
    stop("lambda1 must be positive");
  }

  if (lambda2 <= 0) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals) {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals) {
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

  if (fixedFollowup && R_isnancpp(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      LogicalVector efficacyStopping1(kMax, 1);
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector efficacyP(kMax);
  for (int i=0; i<kMax; i++) {
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
  }

  NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, li, zero, informationRates1);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));

  double theta10 = log(rateRatioLower), theta20 = log(rateRatioUpper);
  NumericVector theta = rep(log(rateRatio), kMax);

  NumericVector u0(1);
  List na;
  DataFrame nb, nc;
  NumericVector I(kMax), time(kMax);

  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || R_isnancpp(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }
  u0[0] = studyDuration1;

  // obtain the timing of interim analysis
  na = nbstat(u0, 1, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup, 0);
  nb = as<DataFrame>(na["resultsUnderH1"]);
  double maxInformation = sum(NumericVector(nb[18]));
  I = maxInformation*informationRates1;

  double information = 0;

  // Lambda function
  auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, kappa1, kappa2,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information](double t)->double {
              NumericVector u0(1, t);
              List na = nbstat(
                u0, 1, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, kappa1, kappa2,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup, 0);
              DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
              return sum(NumericVector(nb[18])) - information;
            };

  for (int i=0; i<kMax-1; i++) {
    // match the predicted information to the target
    information = std::max(I[i], 0.0);
    time[i] = brent(f, 0.0001, studyDuration1, 0.0001);
  }
  time[kMax-1] = studyDuration1;

  // obtain the variance ratio for the lower equivalence limit
  na = nbstat(time, rateRatioLower, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              nullVariance);
  nb = as<DataFrame>(na["resultsUnderH1"]);
  nc = as<DataFrame>(na["resultsUnderH0"]);
  NumericVector varianceRatioH10 = as<NumericVector>(nc["varianceRatio"]);

  // obtain the variance ratio for the upper equivalence limit
  na = nbstat(time, rateRatioUpper, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              nullVariance);
  nc = as<DataFrame>(na["resultsUnderH0"]);
  NumericVector varianceRatioH20 = as<NumericVector>(nc["varianceRatio"]);
  NumericVector wH10 = sqrt(varianceRatioH10);
  NumericVector wH20 = sqrt(varianceRatioH20);


  // calculate cumulative rejection probability under H1
  NumericVector b = criticalValues1;
  NumericVector l = b*wH10 + theta10*sqrt(I);
  NumericVector u = -b*wH20 + theta20*sqrt(I);

  List probs1 = exitprobcpp(pmax(l, li), li, theta, I);
  List probs2 = exitprobcpp(ui, pmin(u, ui), theta, I);

  NumericVector cpl = cumsum(NumericVector(probs1[0]));
  NumericVector cpu = cumsum(NumericVector(probs2[1]));


  // identify the last look with l[k] >= u[k] if it exists
  IntegerVector k = which(l >= u);
  NumericVector cp(kMax);
  if (k.size() == 0) {
    cp = cpl + cpu - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    List a = exitprobcpp(l[idx], u[idx], theta[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cp[i] = cpl[i] + cpu[i] - ca[i];
      } else {
        cp[i] = cpl[i] + cpu[i] - 1;
      }
    }
  }

  // incremental exit probabilities under H1
  NumericVector q(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      q[i] = cp[i];
    } else if (i<kMax-1) {
      q[i] = cp[i] - cp[i-1];
    } else {
      q[i] = 1 - cp[i-1];
    }
  }

  NumericVector rejectPerStage(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      rejectPerStage[i] = cp[i];
    } else {
      rejectPerStage[i] = cp[i] - cp[i-1];
    }
  }

  NumericVector efficacyRateRatioLower = exp(theta10 + b/sqrt(I)*wH10);
  NumericVector efficacyRateRatioUpper = exp(theta20 - b/sqrt(I)*wH20);

  // calculate cumulative rejection under H10
  // match the maximum information under H10 with that under H1
  double accrualDurationH10, followupTimeH10, studyDurationH10;
  if (!fixedFollowup) {
    auto h = [rateRatioLower, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval*followupTime);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda2*rateRatioLower, lambda2, gamma1, gamma2,
                  accrualDuration, aval*followupTime, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
              };

    double aval = brent(h, 1, 2, 1e-6); // increase the follow-up time
    followupTimeH10 = aval*followupTime;
    accrualDurationH10 = accrualDuration;
    studyDurationH10 = accrualDurationH10 + followupTimeH10;
  } else { // fixed follow-up
    auto h = [rateRatioLower, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, aval*accrualDuration + followupTime);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda2*rateRatioLower, lambda2, gamma1, gamma2,
                  aval*accrualDuration, followupTime, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
              };
    double aval = brent(h, 1, 2, 1e-6); // increase the accrual duration
    accrualDurationH10 = aval*accrualDuration;
    followupTimeH10 = followupTime;
    studyDurationH10 = accrualDurationH10 + followupTimeH10;
  }

  // obtain the timing of interim analysis under H10
  auto fH10 = [rateRatioLower, allocationRatioPlanned,
               accrualTime, accrualIntensity,
               piecewiseSurvivalTime, kappa1, kappa2,
               lambda2, gamma1, gamma2,
               accrualDurationH10, followupTimeH10, fixedFollowup,
               &information](double t)->double {
                 NumericVector u0(1, t);
                 List na = nbstat(
                   u0, 1, allocationRatioPlanned,
                   accrualTime, accrualIntensity,
                   piecewiseSurvivalTime, kappa1, kappa2,
                   lambda2*rateRatioLower, lambda2, gamma1, gamma2,
                   accrualDurationH10, followupTimeH10, fixedFollowup, 0);
                 DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                 return sum(NumericVector(nb[18])) - information;
               };

  NumericVector timeH10(kMax);
  for (int i=0; i<kMax-1; i++) {
    // match the predicted information to the target
    information = std::max(I[i], 0.0);
    timeH10[i] = brent(fH10, 0.0001, studyDurationH10, 0.0001);
  }
  timeH10[kMax-1] = studyDurationH10;


  // obtain the variance ratio for H20 with respect to H10
  List naH21 = nbstat(timeH10, rateRatioUpper, allocationRatioPlanned,
                      accrualTime, accrualIntensity,
                      piecewiseSurvivalTime, kappa1, kappa2,
                      lambda2*rateRatioLower, lambda2, gamma1, gamma2,
                      accrualDurationH10, followupTimeH10, fixedFollowup,
                      nullVariance);

  DataFrame ncH21 = as<DataFrame>(naH21["resultsUnderH0"]);
  NumericVector varianceRatioH21 =
    as<NumericVector>(ncH21["varianceRatio"]);
  NumericVector wH21 = sqrt(varianceRatioH21);

  NumericVector bH10 = -b*wH21 + (theta20 - theta10)*sqrt(I);
  List probsH10 = exitprobcpp(ui, pmin(bH10, ui), zero, I);

  NumericVector cplH10 = cumsum(NumericVector(probsH10[1]));
  NumericVector cpuH10 = cumAlphaSpent;

  // identify the last look with b[k] > bH10[k] if it exists
  IntegerVector kH10 = which(b > bH10);
  NumericVector cpH10(kMax);
  if (kH10.size() == 0) {
    cpH10 = cplH10 + cpuH10 - 1;
  } else {
    int K = max(kH10);
    IntegerVector idx = Range(0, K);
    List aH10 = exitprobcpp(b[idx], bH10[idx], zero[idx], I[idx]);
    NumericVector caH10 = cumsum(NumericVector(aH10[0]) +
      NumericVector(aH10[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH10[i] = cplH10[i] + cpuH10[i] - caH10[i];
      } else {
        cpH10[i] = cplH10[i] + cpuH10[i] - 1;
      }
    }
  }


  // calculate cumulative rejection under H20
  // match the maximum information under H20 with that under H1
  double accrualDurationH20, followupTimeH20, studyDurationH20;
  if (!fixedFollowup) {
    auto h = [rateRatioUpper, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval*followupTime);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda2*rateRatioUpper, lambda2, gamma1, gamma2,
                  accrualDuration, aval*followupTime, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
              };

    if (h(0) < 0) { // adjust the follow-up time
      double aval = brent(h, 0, 1, 1e-6);
      followupTimeH20 = aval*followupTime;
      accrualDurationH20 = accrualDuration;
      studyDurationH20 = accrualDurationH20 + followupTimeH20;
    } else { // adjust the accrual duration
      auto g = [rateRatioUpper, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, kappa1, kappa2,
                lambda2, gamma1, gamma2,
                accrualDuration, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval*accrualDuration);
                  List na = nbstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, kappa1, kappa2,
                    lambda2*rateRatioUpper, lambda2, gamma1, gamma2,
                    aval*accrualDuration, 0, fixedFollowup, 0);
                  DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[18])) - maxInformation;
                };

      double aval = brent(g, 0.1, 1, 1e-6);
      accrualDurationH20 = aval*accrualDuration;
      followupTimeH20 = 0;
      studyDurationH20 = accrualDurationH20 + followupTimeH20;
    }
  } else { // fixed follow-up
    auto h = [rateRatioUpper, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval*followupTime);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda2*rateRatioUpper, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
              };

    if (h(0) < 0) { // adjust the study duration
      accrualDurationH20 = accrualDuration;
      followupTimeH20 = followupTime;
      double aval = brent(h, 0, 1, 1e-6);
      studyDurationH20 = accrualDurationH20 + aval*followupTimeH20;
    } else { // adjust the accrual duration
      auto g = [rateRatioUpper, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, kappa1, kappa2,
                lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval*accrualDuration);
                  List na = nbstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, kappa1, kappa2,
                    lambda2*rateRatioUpper, lambda2, gamma1, gamma2,
                    aval*accrualDuration, followupTime, fixedFollowup, 0);
                  DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[18])) - maxInformation;
                };

      double aval = brent(g, 0.1, 1, 1e-6);
      accrualDurationH20 = aval*accrualDuration;
      followupTimeH20 = followupTime;
      studyDurationH20 = accrualDurationH20;
    }
  }

  // obtain the timing of interim analysis under H20
  auto fH20 = [rateRatioUpper, allocationRatioPlanned,
               accrualTime, accrualIntensity,
               piecewiseSurvivalTime, kappa1, kappa2,
               lambda2, gamma1, gamma2,
               accrualDurationH20, followupTimeH20, fixedFollowup,
               &information](double t)->double {
                 NumericVector u0(1, t);
                 List na = nbstat(
                   u0, 1, allocationRatioPlanned,
                   accrualTime, accrualIntensity,
                   piecewiseSurvivalTime, kappa1, kappa2,
                   lambda2*rateRatioUpper, lambda2, gamma1, gamma2,
                   accrualDurationH20, followupTimeH20, fixedFollowup, 0);
                 DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                 return sum(NumericVector(nb[18])) - information;
               };

  NumericVector timeH20(kMax);
  for (int i=0; i<kMax-1; i++) {
    // match the predicted information to the target
    information = std::max(I[i], 0.0);
    timeH20[i] = brent(fH20, 0.0001, studyDurationH20, 0.0001);
  }
  timeH20[kMax-1] = studyDurationH20;

  // obtain the variance ratio for H10 with respect to H20
  List naH12 = nbstat(timeH20, rateRatioLower, allocationRatioPlanned,
                      accrualTime, accrualIntensity,
                      piecewiseSurvivalTime, kappa1, kappa2,
                      lambda2*rateRatioUpper, lambda2, gamma1, gamma2,
                      accrualDurationH20, followupTimeH20, fixedFollowup,
                      nullVariance);

  DataFrame ncH12 = as<DataFrame>(naH12["resultsUnderH0"]);
  NumericVector varianceRatioH12 = as<NumericVector>(ncH12["varianceRatio"]);
  NumericVector wH12 = sqrt(varianceRatioH12);

  NumericVector bH20 = b*wH12 + (theta10 - theta20)*sqrt(I);

  List probsH20 = exitprobcpp(pmax(bH20, li), li, zero, I);

  NumericVector cpuH20 = cumsum(NumericVector(probsH20[0]));
  NumericVector cplH20 = cumAlphaSpent;

  // identify the last look with bH20[k] >= -b[k] if it exists
  IntegerVector kH20 = which(bH20 > -b);
  NumericVector cpH20(kMax);
  if (kH20.size() == 0) {
    cpH20 = cplH20 + cpuH20 - 1;
  } else {
    int K = max(kH20);
    IntegerVector idx = Range(0, K);
    NumericVector uH20 = -b;
    List aH20 = exitprobcpp(bH20[idx], uH20[idx], zero[idx], I[idx]);
    NumericVector caH20 = cumsum(NumericVector(aH20[0]) +
      NumericVector(aH20[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH20[i] = cplH20[i] + cpuH20[i] - caH20[i];
      } else {
        cpH20[i] = cplH20[i] + cpuH20[i] - 1;
      }
    }
  }


  double r1 = allocationRatioPlanned/(1+allocationRatioPlanned);

  NumericVector nsubjects = NumericVector(nb[1]);
  NumericVector nsubjects1 = r1*nsubjects;
  NumericVector nsubjects2 = (1-r1)*nsubjects;
  NumericVector nevents = NumericVector(nb[2]);
  NumericVector nevents1 = NumericVector(nb[3]);
  NumericVector nevents2 = NumericVector(nb[4]);
  NumericVector ndropouts = NumericVector(nb[5]);
  NumericVector ndropouts1 = NumericVector(nb[6]);
  NumericVector ndropouts2 = NumericVector(nb[7]);
  NumericVector exposure = NumericVector(nb[11]);
  NumericVector exposure1 = NumericVector(nb[12]);
  NumericVector exposure2 = NumericVector(nb[13]);

  double overallReject = cp[kMax-1];
  double attainedAlphaH10 = cpH10[kMax-1];
  double attainedAlphaH20 = cpH20[kMax-1];

  double expectedNumberOfEvents = sum(q*nevents);
  double expectedNumberOfEvents1 = sum(q*nevents1);
  double expectedNumberOfEvents2 = sum(q*nevents2);
  double expectedNumberOfDropouts = sum(q*ndropouts);
  double expectedNumberOfDropouts1 = sum(q*ndropouts1);
  double expectedNumberOfDropouts2 = sum(q*ndropouts2);
  double expectedNumberOfSubjects = sum(q*nsubjects);
  double expectedNumberOfSubjects1 = sum(q*nsubjects1);
  double expectedNumberOfSubjects2 = sum(q*nsubjects2);
  double expectedExposure = sum(q*exposure);
  double expectedExposure1 = sum(q*exposure1);
  double expectedExposure2 = sum(q*exposure2);
  double expectedStudyDuration = sum(q*time);
  double expectedInformation = sum(q*I);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = alpha,
    _["attainedAlphaH10"] = attainedAlphaH10,
    _["attainedAlphaH20"] = attainedAlphaH20,
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfDropouts"] = (ndropouts[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["exposure"] = (exposure[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfDropouts"] = expectedNumberOfDropouts,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedExposure"] = expectedExposure,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["kMax"] = kMax,
    _["rateRatioLower"] = rateRatioLower,
    _["rateRatioUpper"] = rateRatioUpper,
    _["rateRatio"] = lambda1/lambda2);

  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["rejectPerStage"] = rejectPerStage,
    _["cumulativeRejection"] = cp,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["cumulativeAttainedAlphaH10"] = cpH10,
    _["cumulativeAttainedAlphaH20"] = cpH20,
    _["numberOfEvents"] = nevents,
    _["numberOfDropouts"] = ndropouts,
    _["numberOfSubjects"] = nsubjects,
    _["exposure"] = exposure,
    _["analysisTime"] = time,
    _["efficacyRateRatioLower"] = efficacyRateRatioLower,
    _["efficacyRateRatioUpper"] = efficacyRateRatioUpper,
    _["efficacyP"] = efficacyP,
    _["information"] = I);

  DataFrame varianceRatios = DataFrame::create(
    _["varianceRatioH10"] = varianceRatioH10,
    _["varianceRatioH20"] = varianceRatioH20,
    _["varianceRatioH12"] = varianceRatioH12,
    _["varianceRatioH21"] = varianceRatioH21);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["allocationRatioPlanned"] = allocationRatioPlanned,
    _["accrualTime"] = accrualTime,
    _["accrualIntensity"] = accrualIntensity,
    _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
    _["kappa1"] = kappa1,
    _["kappa2"] = kappa2,
    _["lambda1"] = lambda1,
    _["lambda2"] = lambda2,
    _["gamma1"] = gamma1,
    _["gamma2"] = gamma2,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["spendingTime"] = spendingTime,
    _["nullVariance"] = nullVariance,
    _["varianceRatios"] = varianceRatios);

  List byTreatmentCounts = List::create(
    _["numberOfEvents1"] = nevents1,
    _["numberOfDropouts1"] = ndropouts1,
    _["numberOfSubjects1"] = nsubjects1,
    _["exposure1"] = exposure1,
    _["numberOfEvents2"] = nevents2,
    _["numberOfDropouts2"] = ndropouts2,
    _["numberOfSubjects2"] = nsubjects2,
    _["exposure2"] = exposure2,
    _["expectedNumberOfEvents1"] = expectedNumberOfEvents1,
    _["expectedNumberOfDropouts1"] = expectedNumberOfDropouts1,
    _["expectedNumberOfSubjects1"] = expectedNumberOfSubjects1,
    _["expectedExposure1"] = expectedExposure1,
    _["expectedNumberOfEvents2"] = expectedNumberOfEvents2,
    _["expectedNumberOfDropouts2"] = expectedNumberOfDropouts2,
    _["expectedNumberOfSubjects2"] = expectedNumberOfSubjects2,
    _["expectedExposure2"] = expectedExposure2);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings,
    _["byTreatmentCounts"] = byTreatmentCounts);

  result.attr("class") = "nbpowerequiv";

  return result;
}



//' @title Sample size for equivalence in negative binomial rate ratio
//' @description Obtains the sample size for equivalence in negative binomial
//' rate ratio.
//'
//' @param beta The type II error.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_criticalValues
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param rateRatioLower The lower equivalence limit of rate ratio.
//' @param rateRatioUpper The upper equivalence limit of rate ratio.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @param kappa1 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the active treatment group.
//' @param kappa2 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the control group.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group.
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}. Adjustment
//'   may be needed for non-monotone relationship with study power.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//' @param nullVariance Whether to calculate the variance for log rate ratio
//'   under the null hypothesis.
//'
//' @return An S3 class \code{nbpowerequiv} object
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{nbpowerequiv}}
//'
//' @examples
//'
//' # Example 1: Variable follow-up design and solve for follow-up time
//' nbsamplesizeequiv(beta = 0.1, kMax = 2, informationRates = c(0.5, 1),
//'                   alpha = 0.05, typeAlphaSpending = "sfOF",
//'                   rateRatioLower = 2/3, rateRatioUpper = 3/2,
//'                   accrualIntensity = 1956/1.25,
//'                   kappa1 = 5, kappa2 = 5,
//'                   lambda1 = 0.125, lambda2 = 0.125,
//'                   gamma1 = 0, gamma2 = 0,
//'                   accrualDuration = 1.25,
//'                   followupTime = NA, fixedFollowup = FALSE,
//'                   nullVariance = 1)
//'
//' # Example 2: Fixed follow-up design and solve for accrual duration
//' nbsamplesizeequiv(beta = 0.2, kMax = 2, informationRates = c(0.5, 1),
//'                   alpha = 0.05, typeAlphaSpending = "sfOF",
//'                   rateRatioLower = 0.5, rateRatioUpper = 2,
//'                   accrualIntensity = 220/1.5,
//'                   kappa1 = 3, kappa2 = 3,
//'                   lambda1 = 8.4, lambda2 = 8.4,
//'                   gamma1 = 0, gamma2 = 0,
//'                   accrualDuration = NA,
//'                   followupTime = 0.5, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List nbsamplesizeequiv(const double beta = 0.2,
                       const int kMax = 1,
                       const NumericVector& informationRates = NA_REAL,
                       const NumericVector& criticalValues = NA_REAL,
                       const double alpha = 0.05,
                       const String typeAlphaSpending = "sfOF",
                       const double parameterAlphaSpending = NA_REAL,
                       const NumericVector& userAlphaSpending = NA_REAL,
                       const double rateRatioLower = NA_REAL,
                       const double rateRatioUpper = NA_REAL,
                       const double allocationRatioPlanned = 1,
                       const NumericVector& accrualTime = 0,
                       const NumericVector& accrualIntensity = 1500,
                       const NumericVector& piecewiseSurvivalTime = 0,
                       const double kappa1 = 5,
                       const double kappa2 = 5,
                       const double lambda1 = 0.125,
                       const double lambda2 = 0.125,
                       const NumericVector& gamma1 = 0,
                       const NumericVector& gamma2 = 0,
                       double accrualDuration = NA_REAL,
                       double followupTime = NA_REAL,
                       const bool fixedFollowup = 0,
                       const NumericVector& interval =
                         NumericVector::create(0.001, 240),
                         const NumericVector& spendingTime = NA_REAL,
                         const bool rounding = 1,
                         const bool nullVariance = 0) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  int nintervals = piecewiseSurvivalTime.size();

  if (R_isnancpp(beta)) {
    stop("beta must be provided");
  }

  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
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


  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (R_isnancpp(rateRatioLower)) {
    stop("rateRatioLower must be provided");
  }

  if (R_isnancpp(rateRatioUpper)) {
    stop("rateRatioUpper must be provided");
  }

  if (rateRatioLower <= 0) {
    stop("rateRatioLower must be positive");
  }

  double rateRatio = lambda1/lambda2;

  if (rateRatioLower >= rateRatio) {
    stop("rateRatioLower must be less than lambda1/lambda2");
  }

  if (rateRatioUpper <= rateRatio) {
    stop("rateRatioUpper must be greater than lambda1/lambda2");
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

  if (kappa1 < 0) {
    stop("kappa1 must be nonnegative");
  }

  if (kappa2 < 0) {
    stop("kappa2 must be nonnegative");
  }

  if (lambda1 <= 0) {
    stop("lambda1 must be positive");
  }

  if (lambda2 <= 0) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (gamma1.size() != 1 && gamma1.size() != nintervals) {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() != 1 && gamma2.size() != nintervals) {
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

  if (fixedFollowup && R_isnancpp(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
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

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      LogicalVector efficacyStopping1(kMax, 1);
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }


  std::string unknown;
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


  NumericVector b = criticalValues1;
  NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);
  double theta10 = log(rateRatioLower), theta20 = log(rateRatioUpper);
  NumericVector theta = rep(log(rateRatio), kMax);

  double maxInformation;

  if (!nullVariance) {
    List design = getDesignEquiv(
      beta, NA_REAL, log(rateRatioLower), log(rateRatioUpper),
      log(rateRatio), kMax, informationRates1, criticalValues1,
      alpha, asf, asfpar, userAlphaSpending, spendingTime1,
      1, 1, 1, 1);

    DataFrame overallResults = as<DataFrame>(design["overallResults"]);
    maxInformation = overallResults["information"];

    auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              unknown, maxInformation](double aval)-> double{
                NumericVector accrualIntensity1 = clone(accrualIntensity);
                double dur1=0, dur2=0;

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

                // obtain the maximum information at study end
                NumericVector u0(1, dur1 + dur2);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda1, lambda2, gamma1, gamma2,
                  dur1, dur2, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
              };

    if (unknown == "accrualDuration") {
      accrualDuration = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 0.0001);
      accrualIntensity1 = aval*accrualIntensity;
    }
  } else {
    auto f = [beta, kMax, informationRates1, b, li, ui,
              rateRatioLower, rateRatioUpper, theta10, theta20, theta,
              allocationRatioPlanned, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, kappa1, kappa2,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              nullVariance, unknown](double aval)-> double{
                NumericVector accrualIntensity1 = clone(accrualIntensity);
                double dur1=0, dur2=0;

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

                double studyDuration = dur1 + dur2;

                // obtain the timing of interim analysis
                NumericVector u0(1, studyDuration), time(kMax);
                List na = nbstat(u0, 1, allocationRatioPlanned,
                                 accrualTime, accrualIntensity1,
                                 piecewiseSurvivalTime, kappa1, kappa2,
                                 lambda1, lambda2, gamma1, gamma2,
                                 dur1, dur2, fixedFollowup, 0);
                DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                double maxInformation = sum(NumericVector(nb[18]));
                NumericVector I = maxInformation*informationRates1;

                double information = 0;

                // Lambda function
                auto g = [allocationRatioPlanned, accrualTime,
                          accrualIntensity1, piecewiseSurvivalTime,
                          kappa1, kappa2, lambda1, lambda2,
                          gamma1, gamma2, dur1, dur2, fixedFollowup,
                          &information](double t)->double {
                            NumericVector u0(1,t);
                            List na = nbstat(u0, 1, allocationRatioPlanned,
                                             accrualTime, accrualIntensity1,
                                             piecewiseSurvivalTime, kappa1,
                                             kappa2, lambda1, lambda2,
                                             gamma1, gamma2, dur1, dur2,
                                             fixedFollowup, 0);
                            DataFrame nb =
                              as<DataFrame>(na["resultsUnderH1"]);
                            return sum(NumericVector(nb[18])) - information;
                          };

                for (int i=0; i<kMax-1; i++) {
                  // match the predicted information to the target
                  information = std::max(I[i], 0.0);
                  time[i] = brent(g, 0.0001, studyDuration, 0.0001);
                }
                time[kMax-1] = studyDuration;


                // obtain the variance ratio for the lower equivalence limit
                na = nbstat(time, rateRatioLower, allocationRatioPlanned,
                            accrualTime, accrualIntensity1,
                            piecewiseSurvivalTime, kappa1, kappa2,
                            lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup, nullVariance);
                DataFrame nc = as<DataFrame>(na["resultsUnderH0"]);
                NumericVector varianceRatioH10 =
                  as<NumericVector>(nc["varianceRatio"]);

                // obtain the variance ratio for the upper equivalence limit
                na = nbstat(time, rateRatioUpper, allocationRatioPlanned,
                            accrualTime, accrualIntensity1,
                            piecewiseSurvivalTime, kappa1, kappa2,
                            lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup, nullVariance);
                nc = as<DataFrame>(na["resultsUnderH0"]);
                NumericVector varianceRatioH20 =
                  as<NumericVector>(nc["varianceRatio"]);

                NumericVector wH10 = sqrt(varianceRatioH10);
                NumericVector wH20 = sqrt(varianceRatioH20);

                // calculate cumulative rejection probability under H1
                NumericVector l = b*wH10 + theta10*sqrt(I);
                NumericVector u = -b*wH20 + theta20*sqrt(I);

                List probs1 = exitprobcpp(pmax(l, li), li, theta, I);
                List probs2 = exitprobcpp(ui, pmin(u, ui), theta, I);

                double cpl = sum(NumericVector(probs1[0]));
                double cpu = sum(NumericVector(probs2[1]));

                double power;
                if (is_true(any(l <= u))) {
                  power = cpl + cpu - 1;
                } else {
                  List a = exitprobcpp(l, u, theta, I);
                  double p = sum(NumericVector(a[0]) + NumericVector(a[1]));
                  power = cpl + cpu - p;
                }

                return power - (1-beta);
              };

    if (unknown == "accrualDuration") {
      accrualDuration = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 0.0001);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 0.0001);
      accrualIntensity1 = aval*accrualIntensity;
    }

    NumericVector u0(1, accrualDuration + followupTime);
    List na = nbstat(u0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity1,
                     piecewiseSurvivalTime, kappa1, kappa2,
                     lambda1, lambda2, gamma1, gamma2,
                     accrualDuration, followupTime, fixedFollowup, 0);
    DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
    maxInformation = sum(NumericVector(nb[18]));
  }

  double studyDuration = accrualDuration + followupTime;


  // output the results
  if (rounding) {
    NumericVector u0(1, studyDuration);
    double n0 = accrual(u0, accrualTime, accrualIntensity1,
                        accrualDuration)[0];
    double n = std::ceil(n0);

    if (n - n0 > 1e-4) {
      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        double aval = n/n0;
        accrualIntensity1 = aval*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      if (!fixedFollowup) {
        // adjust follow-up time to obtain the target maximum information
        auto h = [allocationRatioPlanned, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval*followupTime);
                    List na = nbstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, kappa1, kappa2,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval*followupTime, fixedFollowup, 0);
                    DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[18])) - maxInformation;
                  };

        double aval = brent(h, 0, 1.1, 1e-6);
        followupTime = aval*followupTime;
        studyDuration = accrualDuration + followupTime;
      } else {
        // adjust study duration to obtain the target maximum information
        auto h = [allocationRatioPlanned, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, kappa1, kappa2,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval*followupTime);
                    List na = nbstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, kappa1, kappa2,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup, 0);
                    DataFrame nb = as<DataFrame>(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[18])) - maxInformation;
                  };

        double aval = brent(h, 0, 1, 1e-6);
        studyDuration = accrualDuration + aval*followupTime;
      }
    }
  }

  List result = nbpowerequiv(
    kMax, informationRates1, criticalValues1,
    alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, rateRatioLower, rateRatioUpper,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, kappa1, kappa2,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration, nullVariance);

  return result;
}
