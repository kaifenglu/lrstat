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
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
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
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
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
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
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
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
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
//' calendar times (row) for each treatment group (column).
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
//' * \code{ndropouts1}: The number of events for the active treatment group.
//' 
//' * \code{ndropouts2}: The number of events for the control group.
//' 
//' * \code{nfmax}: The total number of subjects reaching maximum follow-up.
//' 
//' * \code{nfmax1}: The number of subjects reaching maximum follow-up in 
//' the active treatment group. 
//' 
//' * \code{nfmax2}: The number of subjects reaching maximum follow-up in 
//' the control group. 
//' 
//' If \code{predictEventOnly = 0}, the following variables will also 
//' be included:
//' 
//' * \code{uscore}: The numerator of the weighted log-rank test statistic.
//' 
//' * \code{vscore}: The variance of the weighted log-rank score statistic 
//' with weight squared.
//' 
//' * \code{iscore}: The Fisher information of the weighted log-rank score 
//' statistic. 
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
  atrisktx(numSubintervals), event1(numSubintervals), event2(numSubintervals),
  eventt(numSubintervals), km(numSubintervals), w(numSubintervals);
  NumericVector uscore(nstrata), vscore(nstrata), iscore(nstrata);
  NumericVector nevents1(nstrata), nevents2(nstrata), neventst(nstrata);
  NumericVector ndropouts1(nstrata), ndropouts2(nstrata), ndropoutst(nstrata);
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
//' @param predictTarget The target of prediction. 
//' Set \code{predictTarget = 1} to predict the number of events only.
//' Set \code{predictTarget = 2} (default) to predict the umber of events 
//' and log-rank score statistic mean and variance.
//' Set \code{predictTarget = 3} to predict the number of events,
//' log-rank score statistic mean and variance, and
//' hazard ratio and variance of log hazard ratio.
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
//' * \code{ndropouts1}: The number of events for the active treatment group.
//' 
//' * \code{ndropouts2}: The number of events for the control group.
//' 
//' * \code{nfmax}: The total number of subjects reaching maximum follow-up.
//' 
//' * \code{nfmax1}: The number of subjects reaching maximum follow-up in 
//' the active treatment group. 
//' 
//' * \code{nfmax2}: The number of subjects reaching maximum follow-up in 
//' the control group. 
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
                DataFrame df = lrstat1(time1, hazardRatio, 
                                       allocationRatioPlanned,
                                       accrualTime, accrualIntensity,
                                       piecewiseSurvivalTime, stratumFraction,
                                       lambda1x, lambda2x, gamma1x, gamma2x,
                                       accrualDuration, followupTime, 
                                       fixedFollowup,
                                       rho1, rho2, numSubintervals, 
                                       predictEventOnly);
                
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
//' * \code{surv1}: The milestone survival probability for the treatment group.
//' 
//' * \code{surv2}: The milestone survival probability for the control group.
//' 
//' * \code{vsurv1}: The variance for \code{surv1}.
//' 
//' * \code{vsurv2}: The variance for \code{surv2}.
//' 
//' * \code{survdiff}: The difference in milestone survival probabilities, 
//' i.e., \code{surv1 - surv2}.
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
//' survival.
//' 
//' * \code{subjects}: The enrolled number of subjects. 
//' 
//' * \code{milestone}: The milestone time relative to randomization. 
//' 
//' * \code{surv1}: The milestone survival probability for the treatment group.
//' 
//' * \code{surv2}: The milestone survival probability for the control group.
//' 
//' * \code{vsurv1}: The variance for \code{surv1}.
//' 
//' * \code{vsurv2}: The variance for \code{surv2}.
//' 
//' * \code{survdiff}: The difference in milestone survival probabilities, 
//' i.e., \code{surv1 - surv2}.
//' 
//' * \code{vsurvdiff}: The variance for \code{survdiff}.
//' 
//' * \code{survdiffZ}: The Z-statistic value, i.e., 
//' \code{survdiff/sqrt(vsurvdiff)}.
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
//'       accrualTime = seq(0, 9),
//'       accrualIntensity = c(26/9*seq(1, 9), 26),
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
                t0, 1, allocationRatioPlanned, accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
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
    const int numSubintervals = 300,
    const NumericVector& spendingTime = NA_REAL) {
  
  NumericVector st = clone(spendingTime);
  
  
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
                 accrualDuration, followupTime, fixedFollowup);
  
  // obtain the mean and variance of log-rank test score statistic at
  // each stage
  lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 2);
  
  NumericVector vscore = NumericVector(lr[12]);
  
  NumericVector theta(kMax); // mean values under H0, initialized to zero
  
  // information time
  NumericVector t = vscore / (vscore[kMax - 1]);
  if (is_true(any(is_na(spendingTime)))) {
    st = clone(t);
  }
  
  
  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  return getBoundcpp(kMax, t, alpha, asf, parameterAlphaSpending, 
                     userAlphaSpending, st, efficacyStopping);
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
                 accrualDuration, followupTime, fixedFollowup);
  
  // obtain the mean and variance of log-rank test score statistic at
  // each stage
  lr = lrstat(time, hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              rho1, rho2, numSubintervals, 2);
  
  NumericVector vscore = NumericVector(lr[12]);
  
  NumericVector theta(kMax); // mean values under H0, initialized to zero
  NumericVector l = rep(-6.0, kMax);
  
  List probs = exitprobcpp(criticalValues, l, theta, vscore);
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
//' @param spendingTime A vector of length \code{kMax} for the error spending 
//'   time at each analysis. Defaults to missing, in which case, it is the 
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design. 
//'   Defaults to missing, which is to be replaced with the sum of 
//'   \code{accrualDuration} and \code{followupTime}. If provided, 
//'   the value is allowed to be less than the sum of \code{accrualDuration} 
//'   and \code{followupTime}.
//' 
//'   
//' @return An S3 class \code{lrpower} object with 4 components:
//'
//'* \code{overallResults}: A data frame containing the following variables: 
//'
//'   - \code{overallReject}: The overall rejection probability.
//'   
//'   - \code{alpha}: The overall significance level.
//'   
//'   - \code{numberOfEvents}: The total number of events. 
//'   
//'   - \code{numberOfDropouts}: The total number of dropouts. 
//'   
//'   - \code{numbeOfSubjects}: The total number of subjects.
//'   
//'   - \code{studyDuration}: The total study duration. 
//'   
//'   - \code{expectedNumberOfEvents}: The expected number of events.
//'   
//'   - \code{expectedNumberOfDropouts}: The expected number of dropouts. 
//'   
//'   - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'   
//'   - \code{expectedStudyDuration}: The expected study duration.
//'   
//'   - \code{accrualDuration}: The accrual duration.
//'   
//'   - \code{followupTime}: The follow-up duration. 
//'   
//'   - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'   
//'   - \code{rho1}: The first parameter of the Fleming-Harrington family 
//'   of weighted log-rank test.
//'   
//'   - \code{rho2}: The second parameter of the Fleming-Harrington family 
//'   of weighted log-rank test.
//'   
//'   - \code{allocationRatioPlanned}: Allocation ratio for the active 
//'   treatment versus control. 
//'   
//'   - \code{kMax}: The number of stages.
//'   
//'   - \code{hazardRatioH0}: The hazard ratio under the null hypothesis.
//'   
//'   - \code{etimateHazardRatio}: Whether to estimate the hazard ratio.
//'   
//'   - \code{typeOfComputation}: The type of computation, 
//'   either "direct" for the direct approximation method, 
//'   or "schoenfeld" for the Schoenfeld method.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//' 
//'   - \code{informationRates}: The information rates.
//'   
//'   - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'   
//'   - \code{futilityBounds}: The futility boundaries on the Z-scale.
//'   
//'   - \code{rejectPerStage}: The probability for efficacy stopping.
//'   
//'   - \code{futilityPerStage}: The probability for futility stopping.
//'   
//'   - \code{cumulativeRejection}: The cumulative probability for efficacy 
//'   stopping.
//'   
//'   - \code{cumulativeFutility}: The cumulative probability for futility 
//'   stopping.
//'   
//'   - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'   
//'   - \code{numberOfEvents}: The number of events.
//'   
//'   - \code{numberOfDropouts}: The number of dropouts.
//'   
//'   - \code{numberOfSubjects}: The number of subjects.
//'   
//'   - \code{analysisTime}: The average time since trial start.
//'   
//'   - \code{efficacyHR}: The efficacy boundaries on the hazard ratio scale.
//'   
//'   - \code{futilityHR}: The futility boundaries on the hazard ratio scale.
//'   
//'   - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'   
//'   - \code{futilityP}: The futility boundaries on the p-value scale.
//'   
//'   - \code{information}: The cumulative information.
//'   
//'   - \code{HR}: The average hazard ratio. 
//'   
//'   - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'   
//'   - \code{futilityStopping}: Whether to allow futility stopping.
//'
//' * \code{settings}: A list containing the following input parameters: 
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{typeBetaSpending},
//'   \code{parameterBetaSpending}, \code{userBetaSpending},
//'   \code{accrualTime}, \code{accuralIntensity},
//'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2}, and
//'   \code{spendingTime}.
//'
//' * \code{byTreatmentCounts}: A list containing the following counts by 
//' treatment group:
//' 
//'   - \code{numberOfEvents1}: The number of events by stage for 
//'   the treatment group.
//'   
//'   - \code{numberOfDropouts1}: The number of dropouts by stage for 
//'   the treatment group.
//'   
//'   - \code{numberOfSubjects1}: The number of subjects by stage for 
//'   the treatment group.
//'   
//'   - \code{numberOfEvents2}: The number of events by stage for 
//'   the control group.
//'   
//'   - \code{numberOfDropouts2}: The number of dropouts by stage for 
//'   the control group.
//'   
//'   - \code{numberOfSubjects2}: The number of subjects by stage for 
//'   the control group.   
//'   
//'   - \code{expectedNumberOfEvents1}: The expected number of events for 
//'   the treatment group.
//'   
//'   - \code{expectedNumberOfDropouts1}: The expected number of dropouts for 
//'   the treatment group.
//'   
//'   - \code{expectedNumberOfSubjects1}: The expected number of subjects for 
//'   the treatment group.
//'   
//'   - \code{expectedNumberOfEvents2}: The expected number of events for 
//'   control group.
//'   
//'   - \code{expectedNumberOfDropouts2}: The expected number of dropouts for 
//'   the control group.
//'   
//'   - \code{expectedNumberOfSubjects2}: The expected number of subjects for 
//'   the control group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//' 
//' @examples
//' # Piecewise accrual, piecewise exponential survival, and 5% dropout by
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
  NumericVector st = clone(spendingTime);
  
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
    stop("Schoenfeld method can only be used for conventional log-rank test");
  }
  
  double hazardRatio;
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
  
  if (fixedFollowup && !R_isnancpp(studyDuration) && 
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }
  
  if (fixedFollowup && !R_isnancpp(studyDuration) && 
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }
  
  
  NumericVector cumAlphaSpent(kMax);
  if (su == "direct") {
    if (is_true(any(is_na(criticalValues)))) {
      criticalValues1 = getCriticalValues(
        kMax, informationRates1, efficacyStopping1,
        alpha1, asf, asfpar, userAlphaSpending, hazardRatioH0,
        allocationRatioPlanned, accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda2, gamma1, gamma2, accrualDuration,
        followupTime, fixedFollowup,
        rho1, rho2, numSubintervals, st);
    }
    
    cumAlphaSpent = getCumAlphaSpent(
      kMax, informationRates1, criticalValues1, hazardRatioH0,
      allocationRatioPlanned, accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda2, gamma1, gamma2, accrualDuration,
      followupTime, fixedFollowup,
      rho1, rho2, numSubintervals);
    
  } else if (su == "schoenfeld") {
    if (is_true(any(is_na(criticalValues)))) {
      criticalValues1 = getCriticalValues(
        kMax, informationRates1, efficacyStopping1,
        alpha1, asf, asfpar, userAlphaSpending, 1,
        allocationRatioPlanned, accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda2, gamma1, gamma2, accrualDuration,
        followupTime, fixedFollowup,
        rho1, rho2, numSubintervals, st);
    }
    
    cumAlphaSpent = getCumAlphaSpent(
      kMax, informationRates1, criticalValues1, 1,
      allocationRatioPlanned, accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda2, gamma1, gamma2, accrualDuration,
      followupTime, fixedFollowup,
      rho1, rho2, numSubintervals);
  }
  
  alpha1 = cumAlphaSpent[kMax - 1];
  
  if (!R_isnancpp(alpha1)) {
    if (alpha1 < 0.00001 || alpha1 >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
  }
  
  if (is_true(any(is_na(criticalValues))) && asf=="user") {
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
  if (!fixedFollowup || R_isnancpp(studyDuration)) {
    u0[0] = accrualDuration + followupTime;
  } else {
    u0[0] = studyDuration;
  }
  
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
                 accrualDuration, followupTime, fixedFollowup);
  
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
  } else{
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
  
  
  // information time
  NumericVector t = vscore / (vscore[kMax - 1]);
  
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
  
  // compute the stagewise exit probabilities for efficacy and futility
  
  List probs;
  if (!missingFutilityBounds || bsf=="none" || kMax==1) {
    probs = exitprobcpp(criticalValues1, futilityBounds1, theta, vscore);
  } else {
    List out = getPower(alpha1, kMax, criticalValues1, theta, vscore, 
                        bsf, bsfpar, st, futilityStopping1);
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
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);
  
  if (estimateHazardRatio) {
    hru = hazardRatioH0*exp(-criticalValues1*sqrt(vlogHR));
    hrl = hazardRatioH0*exp(-futilityBounds1*sqrt(vlogHR));
  }
  
  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      hru[i] = NA_REAL;
    }
    
    if (futilityBounds1[i] == -6) {
      hrl[i] = NA_REAL;
    }
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
    _["gamma2"] = gamma2,
    _["spendingTime"] = st);
  
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


//' @title Get group sequential design
//' @description Obtains the maximum information and stopping boundaries 
//' for a generic group sequential design assuming a constant treatment 
//' effect, or obtains the power given the maximum information and 
//' stopping boundaries.
//'
//' @param beta The type II error.
//' @param IMax The maximum information. Either \code{beta} or \code{IMax} 
//' should be provided while the other one should be missing.
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
//'
//' @return An S3 class \code{design} object with three components:
//' 
//' * \code{overallResults}: A data frame containing the following variables:
//' 
//'   - \code{overallReject}: The overall rejection probability.
//'   
//'   - \code{alpha}: The overall significance level.
//'   
//'   - \code{kMax}: The number of stages.
//'   
//'   - \code{theta}: The parameter value.
//'   
//'   - \code{maxInformation}: The maximum information.
//'   
//'   - \code{expectedInformationH1}: The expected information under H1.
//'   
//'   - \code{expectedInformationH0}: The expected information under H0.
//'   
//'   - \code{drift}: The drift parameter, equal to 
//'   \code{theta*sqrt(maxInformation)}.
//'   
//'   - \code{inflationFactor}: The inflation factor (relative to the 
//'   fixed design).
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//' 
//'   - \code{informationRates}: The information rates.
//'   
//'   - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'   
//'   - \code{futilityBounds}: The futility boundaries on the Z-scale.
//'   
//'   - \code{rejectPerStage}: The probability for efficacy stopping.
//'   
//'   - \code{futilityPerStage}: The probability for futility stopping.
//'   
//'   - \code{cumulativeRejection}: The cumulative probability for efficacy 
//'   stopping.
//'   
//'   - \code{cumulativeFutility}: The cumulative probability for futility 
//'   stopping.
//'   
//'   - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'   
//'   - \code{efficacyTheta}: The efficacy boundaries on the parameter scale.
//'   
//'   - \code{futilityTheta}: The futility boundaries on the parameter scale.
//'   
//'   - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'   
//'   - \code{futilityP}: The futility boundaries on the p-value scale.
//'   
//'   - \code{information}: The cumulative information.
//'   
//'   - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'   
//'   - \code{futilityStopping}: Whether to allow futility stopping.
//'
//' * \code{settings}: A list containing the following input parameters: 
//' 
//'   - \code{typeAlphaSpending}: The type of alpha spending. 
//'   
//'   - \code{parameterAlphaSpending}: The parameter value for alpha spending.
//'   
//'   - \code{userAlphaSpending}: The user defined alpha spending. 
//'   
//'   - \code{typeBetaSpending}: The type of beta spending. 
//'   
//'   - \code{parameterBetaSpending}: The parameter value for beta spending. 
//'   
//'   - \code{userBetaSpending}: The user defined beta spending.
//'   
//'   - \code{spendingTime}: The error spending time at each analysis. 
//'   
//'   - \code{calculationTarget}: The calculation target, \code{beta} or 
//'   \code{IMax}.
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
//' getDesign(beta = 0.2, theta = -log(0.7),
//'           kMax = 2, informationRates = c(0.5,1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           typeBetaSpending = "sfP")
//'           
//' 
//' # Example 2: obtain power given the maximum information
//' getDesign(IMax = 72.5, theta = -log(0.7),
//'           kMax = 3, informationRates = c(0.5, 0.75, 1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           typeBetaSpending = "sfP")
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
               const NumericVector& spendingTime = NA_REAL) {
  
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
    
    NumericVector u(kMax), l(kMax, -6.0), theta0(kMax);
    for (int i=0; i<kMax; i++) {
      u[i] = criticalValues[i];
      if (!efficacyStopping1[i]) u[i] = 6.0;
    }
    
    List probs = exitprobcpp(u, l, theta0, informationRates1);
    alpha1 = sum(NumericVector(probs[0]));
  }
  
  if (!R_isnancpp(alpha1)) {
    if (alpha1 < 0.00001 || alpha1 >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
  } else {
    stop("alpha must be provided for missing criticalValues");
  }
  
  if ((unknown == "IMax") && (beta >= 1-alpha1 || beta < 0.0001)) {
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
      stop("Elements of userBetaSpending must be nonnegnative");
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
  
  if (is_true(any(is_na(criticalValues)))) {
    criticalValues1 = getBoundcpp(kMax, informationRates1, alpha1, 
                                  asf, asfpar, userAlphaSpending, 
                                  spendingTime1, efficacyStopping1);
  }
  
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
  
  
  List probs;
  NumericVector t = informationRates1;
  NumericVector st = spendingTime1;
  NumericVector theta1(kMax);
  if (unknown == "IMax") {
    auto f = [beta, kMax, t, futilityStopping1,
              criticalValues1, &futilityBounds1, 
              bsf, bsfpar, userBetaSpending, st, 
              missingFutilityBounds](double aval)->double {
                
                NumericVector theta1 = rep(aval, kMax);

                // compute stagewise exit probabilities
                if (!missingFutilityBounds || bsf=="none" || kMax==1) {
                  List probs = exitprobcpp(
                    criticalValues1, futilityBounds1, theta1, t);
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
                    epsilon = R::pnorm(criticalValues1[0] -
                      theta1[0]*sqrt(t[0]), 0, 1, 1, 0) - cumBetaSpent;
                    if (epsilon < 0) return -1.0;
                    futilityBounds1[0] = R::qnorm(cumBetaSpent, 0, 1, 1, 0) +
                      theta1[0]*sqrt(t[0]);
                  }
                  
                  
                  // lambda expression for finding futility bound at stage k
                  auto g = [&k, &cumBetaSpent, criticalValues1, 
                            &futilityBounds1, theta1, t](double aval)->double {
                              NumericVector u(k+1), l(k+1);
                              for (int i=0; i<k; i++) {
                                u[i] = criticalValues1[i];
                                l[i] = futilityBounds1[i];
                              }
                              u[k] = 6.0;
                              l[k] = aval;
                              
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
    probs = exitprobcpp(criticalValues1, futilityBounds1, theta1, t);
  } else {
    drift = theta*sqrt(IMax1);
    theta1 = rep(drift, kMax);

    if (!missingFutilityBounds || bsf=="none" || kMax==1) {
      probs = exitprobcpp(criticalValues1, futilityBounds1, theta1, t);
      beta1 = 1 - sum(NumericVector(probs[0]));
    } else {
      List out = getPower(alpha1, kMax, criticalValues1, theta1, t,
                          bsf, bsfpar, st, futilityStopping1);
      
      beta1 = out[0];
      futilityBounds1 = out[1];
      probs = out[2];
    }
  }
  
  double driftf = R::qnorm(1-alpha1, 0, 1, 1, 0) +
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
    efficacyTheta[i] = criticalValues1[i]/sqrt(information[i]);
    futilityTheta[i] = futilityBounds1[i]/sqrt(information[i]);
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
  
  IntegerVector stageNumber = seq_len(kMax);
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
    _["maxInformation"] = IMax1,
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
    _["spendingTime"] = spendingTime1,
    _["calculationTarget"] = unknown);
  
  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);
  
  result.attr("class") = "design";
  
  return result;
}


//' @title Adaptive design at an interim look
//' @description Obtains the conditional power for specified incremental 
//' information given the interim results, parameter value, and data-dependent 
//' changes in the error spending function and the number and spacing of 
//' interim looks. Conversely, obtain the incremental information needed 
//' to attain a specified conditional power given the interim results,
//' parameter value, and data-dependent changes in the error spending 
//' function and the number and spacing of interim looks.
//'
//' @param betaNew The type II error for the secondary trial.
//' @param INew The maximum information of the secondary trial. Either 
//' \code{betaNew} or \code{INew} should be provided while the other one 
//' should be missing.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of 
//'   the primary trial.
//' @param theta The parameter value.
//' @param IMax The maximum information of the primary trial. Must be 
//'   provided if \code{futilityBounds} is missing and 
//'   \code{typeBetaSpending} is not equal to "none".
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
//'
//' @return An \code{adaptDesign} object with two list components: 
//' 
//' * \code{primaryTrial}: A list of selected information for the primary 
//' trial, including \code{L}, \code{zL}, \code{theta}, \code{kMax}, 
//' \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
//' and \code{MullerSchafer}.
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
//' # conditional power for original design at estimated parameter value
//' (des2 = adaptDesign(
//'   betaNew = NA, INew = (n-n1)/(4*sigma1^2), 
//'   L, zL, theta = delta1, 
//'   kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD", 
//'   parameterAlphaSpending = -4))
//' 
//' # conditional power with sample size increase
//' (des2 = adaptDesign(
//'   betaNew = NA, INew = 420/(4*sigma1^2), 
//'   L, zL, theta = delta1, 
//'   kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD", 
//'   parameterAlphaSpending = -4))
//' 
//' # Muller & Schafer (2001) method to design the secondary trial: 
//' # 3-look gamma(-2) spending with 84% power at delta = 4.5 and sigma = 20
//' (des2 = adaptDesign(
//'   betaNew = 0.16, INew = NA, 
//'   L, zL, theta = delta1,
//'   kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD", 
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpendingNew = "sfHSD", 
//'   parameterAlphaSpendingNew = -2))
//'   
//' # incremental sample size for sigma = 20
//' (nNew = 4*sigma1^2*des2$secondaryTrial$overallResults$maxInformation)
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
                 const NumericVector& spendingTimeNew = NA_REAL) {
  
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
    
    NumericVector u(kMax), l(kMax, -6.0), theta0(kMax);
    for (int i=0; i<kMax; i++) {
      u[i] = criticalValues[i];
      if (!es[i]) u[i] = 6.0;
    }
    
    List probs = exitprobcpp(u, l, theta0, t);
    alpha1 = sum(NumericVector(probs[0]));
  }
  
  if (!R_isnancpp(alpha1)) {
    if (alpha1 < 0.00001 || alpha1 >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
  } else {
    stop("alpha must be provided for missing criticalValues");
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
        stop("Elements of userBetaSpendingNew must be nonnegnative");
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
  
  
  // obtain critical values for the primary trial
  if (is_true(any(is_na(criticalValues)))) {
    b = getBoundcpp(kMax, t, alpha1, asf, asfpar, userAlphaSpending, st, es);
  }
  
  // obtain futility bounds for the primary trial
  if (kMax > 1) {
    if (is_true(any(is_na(futilityBounds))) && bsf=="none") {
      a = rep(-6.0, kMax);
      a[kMax-1] = b[kMax-1];
    } else if (is_false(any(is_na(futilityBounds))) && 
      a.size() == kMax-1) {
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
    List out = getPower(alpha1, kMax, b, theta1, IMax*t, bsf, bsfpar, st, fs);
    a = out[1];
  }
  
  
  List des1 = List::create(
    _["L"] = L,
    _["zL"] = zL,
    _["theta"] = theta,
    _["kMax"] = kMax,
    _["informationRates"] = t,
    _["efficacyBounds"] = b,
    _["futilityBounds"] = a,
    _["MullerSchafer"] = MullerSchafer);
  
  
  List des2;
  int k1 = kMax - L;
  NumericVector t1(k1), r1(k1), b1(k1), a1(k1, -6.0), theta0(k1);
  for (int l=0; l<k1; l++) {
    t1[l] = (t[l+L] - t[L-1])/(1 - t[L-1]);
    r1[l] = t[L-1]/t[l+L];
    b1[l] = (b[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
    if (!es[l+L]) b1[l] = 6.0;
  }
  
  if (!MullerSchafer) {
    for (int l=0; l<k1; l++) {
      a1[l] = (a[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
      if (!fs[l+L]) a1[l] = -6.0;
    }
    
    IntegerVector idx = Range(L, kMax-1);
    LogicalVector esNew = es[idx];
    LogicalVector fsNew = fs[idx];
    
    des2 = getDesign(betaNew, INew, theta, k1, t1, esNew, fsNew, 
                     b1, NA_REAL, typeAlphaSpendingNew, 
                     parameterAlphaSpendingNew, 0,
                     a1, typeBetaSpendingNew, parameterBetaSpendingNew, 
                     userBetaSpendingNew, stNew);
  } else {
    List probs = exitprobcpp(b1, a1, theta0, t1);
    double alphaNew = sum(NumericVector(probs[0]));
    
    if (!R_isnancpp(betaNew) && betaNew >= 1-alphaNew) {
      stop("betaNew must be less than 1 minus conditional type I error");
    }
    
    NumericVector b1New(kNew, NA_REAL), a1New(kNew, NA_REAL);
    
    des2 = getDesign(betaNew, INew, theta, kNew, tNew, esNew, fsNew, 
                     b1New, alphaNew, typeAlphaSpendingNew, 
                     parameterAlphaSpendingNew, 0,
                     a1New, typeBetaSpendingNew, parameterBetaSpendingNew, 
                     userBetaSpendingNew, stNew);
  }
  
  List result = List::create(
    _["primaryTrial"] = des1,
    _["secondaryTrial"] = des2);
  
  result.attr("class") = "adaptDesign";
  
  return result;
}


//' @title Get the required number of events from hazard ratios
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
  
  List design = getDesign(beta, NA_REAL, 1, kMax, informationRates, 
                          efficacyStopping, futilityStopping, 
                          criticalValues, alpha, typeAlphaSpending,
                          parameterAlphaSpending, userAlphaSpending, 
                          futilityBounds, typeBetaSpending, 
                          parameterBetaSpending, userBetaSpending, 
                          spendingTime);
  
  DataFrame overallResults = as<DataFrame>(design["overallResults"]);
  double drift = overallResults["drift"];
  
  double delta = -log(hazardRatio) + log(hazardRatioH0);
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);
  
  double D = pow(drift,2)/(phi*(1-phi)*pow(delta,2));
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
//' accrualDuration, followupDuration, or the proportionality constant 
//' of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' Adjustment may be needed for non-monotone relationship with study power.
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
                  double followupTime = 18,
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
  
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  
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
    if (alpha < 0.00001 || alpha >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
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
  
  std::string su = typeOfComputation;
  std::for_each(su.begin(), su.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  if (su != "direct" && su != "schoenfeld") {
    stop("typeOfComputation must be direct or Schoenfeld");
  }
  
  if (su == "schoenfeld" && (rho1 != 0 || rho2 != 0)) {
    stop("Schoenfeld method can only be used for conventional log-rank test");
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
  
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  
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
  }
  
  if (su == "schoenfeld") {
    List design = getDesign(beta, NA_REAL, 1, kMax, informationRates1, 
                            efficacyStopping1, 
                            futilityStopping1, criticalValues1, alpha, 
                            asf, asfpar, userAlphaSpending, futilityBounds1, 
                            bsf, bsfpar, userBetaSpending, spendingTime);
    DataFrame byStageResults = as<DataFrame>(design["byStageResults"]);
    criticalValues1 = byStageResults["efficacyBounds"];
    futilityBounds1 = byStageResults["futilityBounds"];
    
    DataFrame overallResults = as<DataFrame>(design["overallResults"]);
    double drift = overallResults["drift"];
    double delta = -log(hazardRatio/hazardRatioH0);
    double r1 = allocationRatioPlanned/(allocationRatioPlanned+1);
    double D = pow(drift/delta, 2)/(r1*(1-r1));
    
    auto f = [hazardRatioH0, allocationRatioPlanned, accrualTime, 
              accrualIntensity, piecewiseSurvivalTime, stratumFraction, 
              lambda1, lambda2, gamma1, gamma2, accrualDuration, 
              followupTime, fixedFollowup, unknown, D](double aval)-> double{
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
                
                NumericVector u0(1);
                DataFrame lr;
                
                // obtain the total number of events at study end
                u0[0] = dur1 + dur2;
                lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
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
              efficacyStopping1, futilityStopping1,
              &criticalValues1, alpha, asf, asfpar, userAlphaSpending,
              &futilityBounds1, bsf, bsfpar, userBetaSpending, hazardRatioH0,
              allocationRatioPlanned, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction, lambda1, lambda2,
              gamma1, gamma2, accrualDuration, followupTime,
              fixedFollowup, rho1, rho2, numSubintervals, 
              su, spendingTime, unknown, missingCriticalValues,
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
                  criticalValues1 = getCriticalValues(
                    kMax, informationRates1, efficacyStopping1,
                    alpha, asf, asfpar, userAlphaSpending, hazardRatioH0,
                    allocationRatioPlanned, accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction, 
                    lambda2, gamma1, gamma2,
                    dur1, dur2, fixedFollowup, rho1, rho2, 
                    numSubintervals, spendingTime);
                }
                
                if (kMax > 1) {
                  if (missingFutilityBounds && bsf=="none") {
                    futilityBounds1 = rep(-6.0, kMax);
                    futilityBounds1[kMax-1] = criticalValues1[kMax-1];
                  } else if (!missingFutilityBounds && 
                    futilityBounds1.size() == kMax-1) {
                    futilityBounds1.push_back(criticalValues1[kMax-1]);
                  } else if (!missingFutilityBounds && 
                    futilityBounds1.size() < kMax-1) {
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
                
                // information time
                NumericVector t = vscore / (vscore[kMax - 1]);
                
                NumericVector st = clone(spendingTime);
                if (is_true(any(is_na(spendingTime)))) {
                  st = clone(t);
                }
                
                // compute stagewise exit probabilities
                if (!missingFutilityBounds || bsf=="none" || kMax==1) {
                  List probs = exitprobcpp(criticalValues1, futilityBounds1, 
                                           theta, vscore);
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
                      
                      if (g(-6.0) > 0) { // no beta spent at the current visit
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
    NumericVector nevents = floor(D*informationRates1 + 0.5);
    
    // new information rates for boundary calculations
    informationRates1 = nevents/nevents[kMax-1];
    
    // sample size before rounding up
    double n;
    
    if (fixedFollowup) {
      // adjust accrual intensity or duration to obtain int number of events
      if (unknown == "accrualIntensity") {
        double aval = D/D0;
        accrualIntensity1 = aval*accrualIntensity1;
      } else {
        auto h = [hazardRatioH0, allocationRatioPlanned, accrualTime, 
                  accrualIntensity1, piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2, accrualDuration, 
                  followupTime, fixedFollowup, D](double aval)->double {
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
      
      u[0] = accrualDuration;
      NumericVector ns = accrual(u, accrualTime, accrualIntensity1, u[0]);
      n = ns[0];      
    } else {
      n = sum(NumericVector(lr[1]));
    }
  
    // round up the sample size
    NumericVector frac = NumericVector::create(
      1/1, 1/2, 1/3, 1/4, 1/5, 2/1, 2/3, 2/5, 3/1, 3/2, 3/4, 3/5, 
      4/1, 4/3, 4/5, 5/1, 5/2, 5/3, 5/4);
    
    IntegerVector blocksize = IntegerVector::create(
      2, 3, 4, 5, 6, 3, 5, 7, 4, 5, 7, 8,
      5, 7, 9, 6, 7, 8, 9);
    
    int i = which_min(abs(frac - allocationRatioPlanned));
    int b = blocksize[i];
    
    double n0 = n;
    n = std::ceil(n/b)*b;
    
    // adjust accrual intensity or duration to obtain int number of subjects
    if (unknown == "accrualIntensity") {
      double aval = n/n0;
      accrualIntensity1 = aval*accrualIntensity1;
    } else {
      NumericVector ns(1);
      ns[0] = n;
      u = getAccrualDurationFromN(ns, accrualTime, accrualIntensity1);
      accrualDuration = u[0];
    }
    
    // adjust study duration to obtain integer number of events
    auto h = [hazardRatioH0, allocationRatioPlanned, accrualTime, 
              accrualIntensity1, piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2, accrualDuration, 
              followupTime, fixedFollowup, D](double aval)->double {
                NumericVector u(1);
                u[0] = accrualDuration + aval*followupTime;
                double followupTime1 = fixedFollowup ? followupTime : 
                  aval*followupTime;
                
                DataFrame lr = lrstat(
                  u, hazardRatioH0, allocationRatioPlanned, 
                  accrualTime, accrualIntensity1, 
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime1, 
                  fixedFollowup, 0, 0, 1, 1);
                return sum(NumericVector(lr[2])) - D;
              };
    
    double aval;
    if (!fixedFollowup) {
      aval = brent(h, 0.5, 1.1, 1e-6);
    } else {
      aval = brent(h, 0.5, 1, 1e-6);
    }
    
    double studyDuration = accrualDuration + aval*followupTime;
    if (!fixedFollowup) {
      followupTime = aval*followupTime; 
    }
    
    // recalculate boundaries
    resultH1 = lrpower(kMax, informationRates1,
                       efficacyStopping1, futilityStopping1, criticalValues,
                       alpha, typeAlphaSpending, parameterAlphaSpending,
                       userAlphaSpending, futilityBounds,
                       typeBetaSpending, parameterBetaSpending, hazardRatioH0,
                       allocationRatioPlanned, accrualTime, accrualIntensity1,
                       piecewiseSurvivalTime, stratumFraction,
                       lambda1, lambda2, gamma1, gamma2,
                       accrualDuration, followupTime, fixedFollowup,
                       rho1, rho2, numSubintervals, estimateHazardRatio,
                       typeOfComputation, spendingTime, studyDuration);
  } else {
    resultH1 = lrpower(kMax, informationRates1,
                       efficacyStopping1, futilityStopping1, criticalValues1,
                       alpha, typeAlphaSpending, parameterAlphaSpending,
                       userAlphaSpending, futilityBounds1,
                       typeBetaSpending, parameterBetaSpending, hazardRatioH0,
                       allocationRatioPlanned, accrualTime, accrualIntensity1,
                       piecewiseSurvivalTime, stratumFraction,
                       lambda1, lambda2, gamma1, gamma2,
                       accrualDuration, followupTime, fixedFollowup,
                       rho1, rho2, numSubintervals, estimateHazardRatio,
                       typeOfComputation, spendingTime);
  }
  
  // results under H0 by matching the total number of events
  // adjust study duration to obtain integer number of events
  DataFrame overallResults = as<DataFrame>(resultH1["overallResults"]);
  double D = overallResults["numberOfEvents"];
  double aval, studyDuration;
  
  // fix accrualDuration
  auto h = [hazardRatioH0, allocationRatioPlanned, accrualTime, 
            accrualIntensity1, piecewiseSurvivalTime, stratumFraction,
            lambda2, gamma1, gamma2, accrualDuration, 
            followupTime, fixedFollowup, D](double aval)->double {
              NumericVector u(1);
              u[0] = accrualDuration + aval*followupTime;
              double followupTime1 = fixedFollowup ? followupTime : 
                aval*followupTime;
              
              DataFrame lr = lrstat(
                u, hazardRatioH0, allocationRatioPlanned, 
                accrualTime, accrualIntensity1, 
                piecewiseSurvivalTime, stratumFraction,
                lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                accrualDuration, followupTime1, 
                fixedFollowup, 0, 0, 1, 1);
              return sum(NumericVector(lr[2])) - D;
            };
  
  if (!fixedFollowup || h(0) < 0) { // reduce study duration
    aval = brent(h, 0, 1, 1e-6);
    studyDuration = accrualDuration + aval*followupTime;
    if (!fixedFollowup) {
      followupTime = aval*followupTime;
    }
  } else { // reduce accrualDuration for fixed follow-up
    auto g = [hazardRatioH0, allocationRatioPlanned, accrualTime, 
              accrualIntensity1, piecewiseSurvivalTime, stratumFraction,
              lambda2, gamma1, gamma2, accrualDuration, 
              followupTime, fixedFollowup, D](double aval)->double {
                NumericVector u(1);
                u[0] = aval*accrualDuration;
                
                DataFrame lr = lrstat(
                  u, hazardRatioH0, allocationRatioPlanned, 
                  accrualTime, accrualIntensity1, 
                  piecewiseSurvivalTime, stratumFraction,
                  lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                  aval*accrualDuration, followupTime, 
                  fixedFollowup, 0, 0, 1, 1);
                return sum(NumericVector(lr[2])) - D;
              };
    
    aval = brent(g, 1e-6, 1, 1e-6);
    accrualDuration = aval*accrualDuration;
    studyDuration = accrualDuration;
  }
  
  resultH0 = lrpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending, hazardRatioH0,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    rho1, rho2, numSubintervals, 0,
    typeOfComputation, spendingTime, studyDuration);
  
  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}
