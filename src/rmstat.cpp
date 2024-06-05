#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;


//' @title Restricted mean survival time
//' @description Obtains the restricted mean survival time over an interval.
//'
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//'
//' @return The integral of the survival function from \code{t1} to \code{t2}
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' rmst(t1 = 0, t2 = 7, piecewiseSurvivalTime = c(0, 6),
//'      lambda = c(0.0533, 0.0309))
//'
//' @export
// [[Rcpp::export]]
double rmst(const double t1 = 0,
            const double t2 = NA_REAL,
            const NumericVector& piecewiseSurvivalTime = 0,
            const NumericVector& lambda = NA_REAL) {

  if (R_isnancpp(t2)) {
    stop("t2 must be provided");
  }

  if (t1 < 0) {
    stop("t1 must be non-negative");
  }

  if (t2 < t1) {
    stop("t1 must be less than or equal to t2");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (piecewiseSurvivalTime.size() > 1 &&
      is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(is_na(lambda)))) {
    stop("lambda must be provided");
  }

  if (is_true(any(lambda < 0))) {
    stop("lambda must be non-negative");
  }

  NumericVector t = piecewiseSurvivalTime;

  // identify the time interval containing the specified analysis time
  NumericVector time = NumericVector::create(t1, t2);
  IntegerVector m = findInterval3(time, t) - 1;

  double s1, s2, ch = 0.0, aval = 0.0;
  for (int j=0; j<=m[0]-1; j++) {
    ch += lambda[j]*(t[j+1] - t[j]);
  }

  if (m[0] == m[1]) {
    s1 = exp(-lambda[m[0]]*(t1 - t[m[0]]));
    s2 = exp(-lambda[m[0]]*(t2 - t[m[0]]));
    aval = exp(-ch)*(s1 - s2)/lambda[m[0]];
  } else {
    s1 = exp(-lambda[m[0]]*(t1 - t[m[0]]));
    s2 = exp(-lambda[m[0]]*(t[m[0]+1] - t[m[0]]));
    aval = aval + exp(-ch)*(s1 - s2)/lambda[m[0]];

    for (int j=m[0]+1; j<=m[1]-1; j++) {
      ch += lambda[j-1]*(t[j] - t[j-1]);
      s2 = exp(-lambda[j]*(t[j+1] - t[j]));
      aval = aval + exp(-ch)*(1.0 - s2)/lambda[j];
    }

    ch += lambda[m[1]-1]*(t[m[1]] - t[m[1]-1]);
    s2 = exp(-lambda[m[1]]*(t2 - t[m[1]]));
    aval = aval + exp(-ch)*(1.0 - s2)/lambda[m[1]];
  }

  return aval;
}


// define the integrand function for covrmst
typedef struct {
  double time;
  double tau1;
  double tau2;
  double phi;
  NumericVector accrualTime;
  NumericVector accrualIntensity;
  NumericVector piecewiseSurvivalTime;
  NumericVector lambda;
  NumericVector gamma;
  double accrualDuration;
} rmparams;


void f_rm(double *x, int n, void *ex) {
  rmparams *param = (rmparams *) ex;
  NumericVector u0(n), a1(n), a2(n);
  for (int i=0; i<n; i++) {
    u0[i] = x[i];
    a1[i] = rmst(u0[i], param->tau1, param->piecewiseSurvivalTime,
                 param->lambda);
    a2[i] = rmst(u0[i], param->tau2, param->piecewiseSurvivalTime,
                 param->lambda);
  }

  IntegerVector j = findInterval3(u0, param->piecewiseSurvivalTime) - 1;
  NumericVector lambda = param->lambda[j];
  NumericVector p = patrisk(u0, param->piecewiseSurvivalTime, param->lambda,
                            param->gamma);
  u0 = param->time - u0;
  NumericVector N = accrual(u0, param->accrualTime, param->accrualIntensity,
                            param->accrualDuration);
  u0 = a1*a2*lambda/((param->phi)*N*p);
  for (int i=0; i<n; i++) {
    x[i] = u0[i];
  }
}


//' @title Covariance between restricted mean survival times
//' @description Obtains the covariance between restricted mean survival
//' times at two different time points.
//'
//' @param t2 The calendar time for analysis 2.
//' @param tau1 The milestone time for analysis 1.
//' @param tau2 The milestone time for analysis 2.
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
//' @return The covariance between the restricted mean survival times
//' for each treatment group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' covrmst(t2 = 25, tau1 = 16, tau2 = 18, allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
NumericVector covrmst(const double t2 = NA_REAL,
                      const double tau1 = NA_REAL,
                      const double tau2 = NA_REAL,
                      const double allocationRatioPlanned = 1,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda1 = NA_REAL,
                      const NumericVector& lambda2 = NA_REAL,
                      const NumericVector& gamma1 = 0,
                      const NumericVector& gamma2 = 0,
                      const double accrualDuration = NA_REAL,
                      const double maxFollowupTime = NA_REAL) {

  if (R_isnancpp(t2)) {
    stop("t2 must be provided");
  }

  if (R_isnancpp(tau1)) {
    stop("tau1 must be provided");
  }

  if (tau1 <= 0) {
    stop("tau1 must be positive");
  }

  if (R_isnancpp(tau2)) {
    stop("tau2 must be provided");
  }

  if (tau2 < tau1) {
    stop("tau2 must be greater than or equal to tau1");
  }

  if (t2 <= tau2) {
    stop("t2 must be greater than tau2");
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

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
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

  if (piecewiseSurvivalTime.size() > 1 &&
      is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
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

  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (R_isnancpp(maxFollowupTime)) {
    stop("maxFollowupTime must be provided");
  }

  if (tau2 > maxFollowupTime) {
    stop("tau2 must be less than or equal to maxFollowupTime");
  }


  NumericVector t = piecewiseSurvivalTime;
  double phi = allocationRatioPlanned/(1 + allocationRatioPlanned);

  double tol = 1e-6;
  double umax = std::min(tau1, maxFollowupTime);

  rmparams param1 = {t2, tau1, tau2, phi, accrualTime, accrualIntensity,
                     t, lambda1, gamma1, accrualDuration};
  rmparams param2 = {t2, tau1, tau2, 1-phi, accrualTime, accrualIntensity,
                     t, lambda2, gamma2, accrualDuration};

  double q1 = quad(f_rm, &param1, 0.0, umax, tol)[0];
  double q2 = quad(f_rm, &param2, 0.0, umax, tol)[0];

  return NumericVector::create(q1, q2);
}


//' @title Restricted mean survival time by stratum
//'
//' @description Obtains the restricted mean survival time and associated
//' variance by treatment group and by stratum at a given calendar time.
//'
//' @param time The calendar time for data cut.
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
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
//' @return A data frame containing the following variables:
//'
//' * \code{stratum}: The stratum.
//'
//' * \code{time}: The calendar time since trial start.
//'
//' * \code{subjects}: The number of enrolled subjects.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{rmst1}: The restricted mean survival time for the treatment
//'   group.
//'
//' * \code{rmst2}: The restricted mean survival time for the control group.
//'
//' * \code{rmstDiff}: The difference in restricted mean survival times,
//'   i.e., \code{rmst1 - rmst2}.
//'
//' * \code{vrmst1}: The variance for \code{rmst1}.
//'
//' * \code{vrmst2}: The variance for \code{rmst2}.
//'
//' * \code{vrmstDiff}: The variance for \code{rmstDiff}.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' rmstat1(time = 40,
//'         milestone = 18,
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
DataFrame rmstat1(const double time = NA_REAL,
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
                  const bool fixedFollowup = 0) {

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


  NumericVector ss(1, time);
  double a = accrual(ss, accrualTime, accrualIntensity, accrualDuration)[0];

  // obtain the follow-up time for the first enrolled subject
  double maxFollowupTime;
  if (fixedFollowup) {
    maxFollowupTime = followupTime;
  } else {
    maxFollowupTime = accrualDuration + followupTime;
  }


  IntegerVector l1 = Range(0, nintervals-1);
  IntegerVector l(nintervals);
  NumericVector lam1(nintervals), lam2(nintervals);
  NumericVector gam1(nintervals), gam2(nintervals);
  NumericVector nsubjects(nstrata);
  NumericVector rmst1(nstrata), rmst2(nstrata), rmstDiff(nstrata);
  NumericVector vrmst1(nstrata), vrmst2(nstrata), vrmstDiff(nstrata);
  IntegerVector stratum(nstrata);
  NumericVector calTime(nstrata), mileTime(nstrata);
  DataFrame df;

  for (int h=0; h<nstrata; h++) {
    stratum[h] = h+1;
    calTime[h] = time;
    mileTime[h] = milestone;

    double frac = stratumFraction[h];
    l = h*nintervals + l1;
    lam1 = lambda1x[l];
    lam2 = lambda2x[l];
    gam1 = gamma1x[l];
    gam2 = gamma2x[l];

    // obtain number of enrolled subjects
    nsubjects[h] = frac*a;

    rmst1[h] = rmst(0, milestone, piecewiseSurvivalTime, lam1);
    rmst2[h] = rmst(0, milestone, piecewiseSurvivalTime, lam2);
    rmstDiff[h] = rmst1[h] - rmst2[h];

    NumericVector v = covrmst(
      time, milestone, milestone, allocationRatioPlanned,
      accrualTime, frac*accrualIntensity, piecewiseSurvivalTime,
      lam1, lam2, gam1, gam2, accrualDuration, maxFollowupTime);

    vrmst1[h] = v[0];
    vrmst2[h] = v[1];
    vrmstDiff[h] = v[0] + v[1];
  }

  // output the requested information
  df = DataFrame::create(_["stratum"] = stratum,
                         _["time"] = calTime,
                         _["subjects"] = nsubjects,
                         _["milestone"] = mileTime,
                         _["rmst1"] = rmst1,
                         _["rmst2"] = rmst2,
                         _["rmstDiff"] = rmstDiff,
                         _["vrmst1"] = vrmst1,
                         _["vrmst2"] = vrmst2,
                         _["vrmstDiff"] = vrmstDiff);

  return df;
}


//' @title Stratified difference in restricted mean survival times
//' @description Obtains the stratified restricted mean survival times
//' and difference in restricted mean survival times at given calendar
//' times.
//'
//' @param time A vector of calendar times at which to calculate the
//'   restricted mean survival time.
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
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
//' @return A data frame containing the following variables:
//'
//' * \code{time}: The calendar time at which to calculate the restricted
//'   mean survival time.
//'
//' * \code{subjects}: The number of enrolled subjects.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{rmst1}: The restricted mean survival time for the treatment
//'   group.
//'
//' * \code{rmst2}: The restricted mean survival time for the control group.
//'
//' * \code{rmstDiff}: The difference in restricted mean survival times,
//'   i.e., \code{rmst1 - rmst2}.
//'
//' * \code{vrmst1}: The variance for \code{rmst1}.
//'
//' * \code{vrmst2}: The variance for \code{rmst2}.
//'
//' * \code{vrmstDiff}: The variance for \code{rmstDiff}.
//'
//' * \code{information}: The information for \code{rmstDiff}, equal to
//'   \code{1/vrmstDiff}.
//'
//' * \code{rmstDiffZ}: The Z-statistic value, i.e.,
//'   \code{rmstDiff/sqrt(vrmstDiff)}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' rmstat(time = c(22, 40),
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
DataFrame rmstat(const NumericVector& time = NA_REAL,
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
                 const bool fixedFollowup = 0) {

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi), gamma1x(nsi), gamma2x(nsi);

  if (is_true(any(is_na(time)))) {
    stop("time must be provided");
  }

  if (is_true(any(time <= 0))) {
    stop("time must be positive");
  }

  if (R_isnancpp(milestone)) {
    stop("milestone must be provided");
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

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
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

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
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

  if (is_true(any(time > accrualDuration + followupTime))) {
    stop("time cannot exceed accrualDuration + followupTime");
  }

  if (is_true(any(time <= milestone))) {
    stop("time must be greater than milestone");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }


  int k = time.size();
  NumericVector calTime(k), mileTime(k), subjects(k),
  rmst1(k), rmst2(k), vrmst1(k), vrmst2(k),
  rmstDiff(k), vrmstDiff(k), information(k), rmstDiffZ(k);
  DataFrame df;

  for (int j=0; j<k; j++) {
    df = rmstat1(time[j], milestone, allocationRatioPlanned,
                 accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1x, lambda2x, gamma1x, gamma2x,
                 accrualDuration, followupTime, fixedFollowup);

    calTime[j] = max(NumericVector(df[1]));
    subjects[j] = sum(NumericVector(df[2]));
    mileTime[j] = max(NumericVector(df[3]));
    rmst1[j] = sum(stratumFraction*NumericVector(df[4]));
    rmst2[j] = sum(stratumFraction*NumericVector(df[5]));
    rmstDiff[j] = sum(stratumFraction*NumericVector(df[6]));
    vrmst1[j] = sum(stratumFraction*stratumFraction*NumericVector(df[7]));
    vrmst2[j] = sum(stratumFraction*stratumFraction*NumericVector(df[8]));
    vrmstDiff[j] = sum(stratumFraction*stratumFraction*NumericVector(df[9]));
    information[j] = 1.0/vrmstDiff[j];
    rmstDiffZ[j] = rmstDiff[j]/sqrt(vrmstDiff[j]);
  }

  df = DataFrame::create(_["time"] = calTime,
                         _["subjects"] = subjects,
                         _["milestone"] = mileTime,
                         _["rmst1"] = rmst1,
                         _["rmst2"] = rmst2,
                         _["rmstDiff"] = rmstDiff,
                         _["vrmst1"] = vrmst1,
                         _["vrmst2"] = vrmst2,
                         _["vrmstDiff"] = vrmstDiff,
                         _["information"] = information,
                         _["rmstDiffZ"] = rmstDiffZ);

  return df;
}


//' @title Power for difference in restricted mean survival times
//' @description Estimates the power for testing the difference in
//' restricted mean survival times in a two-sample survival design.
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
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffH0 The difference in restricted mean survival times
//'   under the null hypothesis. Defaults to 0 for superiority test.
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
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{rmpower} object with 3 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{drift}: The drift parameter, equal to
//'       \code{(rmstDiff - rmstDiffH0)*sqrt(information)}.
//'
//'     - \code{inflationFactor}: The inflation factor (relative to the
//'       fixed design).
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
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
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{milestone}: The milestone time relative to randomization.
//'
//'     - \code{rmstDiffH0}: The difference in restricted mean survival
//'       times under the null hypothesis.
//'
//'     - \code{rmst1}: The restricted mean survival time for the
//'       treatment group.
//'
//'     - \code{rmst2}: The restricted mean survival time for the
//'       control group.
//'
//'     - \code{rmstDiff}: The difference in restricted mean survival times,
//'       equal to \code{rmst1 - rmst2}.
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
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRmstDiff}: The efficacy boundaries on the restricted
//'       mean survival time difference scale.
//'
//'     - \code{futilityRmstDiff}: The futility boundaries on the restricted
//'       mean survival time difference scale.
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
//'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   and \code{spendingTime}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survival, and 5% dropout by
//' # the end of 1 year.
//'
//' rmpower(kMax = 2, informationRates = c(0.8, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         milestone = 18,
//'         allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'         accrualIntensity = 100/9*seq(1, 9),
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
List rmpower(const int kMax = 1,
             const NumericVector& informationRates = NA_REAL,
             const LogicalVector& efficacyStopping = NA_LOGICAL,
             const LogicalVector& futilityStopping = NA_LOGICAL,
             const NumericVector& criticalValues = NA_REAL,
             const double alpha = 0.025,
             const std::string typeAlphaSpending = "sfOF",
             const double parameterAlphaSpending = NA_REAL,
             const NumericVector& userAlphaSpending = NA_REAL,
             const NumericVector& futilityBounds = NA_REAL,
             const std::string typeBetaSpending = "none",
             const double parameterBetaSpending = NA_REAL,
             const double milestone = NA_REAL,
             const double rmstDiffH0 = 0,
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
    informationRates1 = NumericVector(tem)/(kMax+0.0);
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
    stop("parameterBetaSpending must be positive for sfKD");
  }

  if (R_isnancpp(milestone)) {
    stop("milestone must be provided");
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

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
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

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
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

  if (milestone >= accrualDuration + followupTime) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
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


  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || R_isnancpp(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }

  // obtain the timing of interim analysis
  DataFrame rm;
  NumericVector time(kMax), rdu(kMax), rdl(kMax);
  NumericVector u0(1, studyDuration1);
  rm = rmstat(u0, milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup);

  double maxInformation = sum(NumericVector(rm[9]));
  double rmst1 = sum(NumericVector(rm[3]));
  double rmst2 = sum(NumericVector(rm[4]));
  double rmstDiff = sum(NumericVector(rm[5]));
  NumericVector theta(kMax, rmstDiff);
  NumericVector I = maxInformation*informationRates1;

  double information1;
  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double aval)->double {
              NumericVector u0(1, aval);
              DataFrame rm = rmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup);
              return sum(NumericVector(rm[9])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    information1 = I[i];
    time[i] = brent(f, milestone + 1.0e-6, studyDuration1, 1.0e-6);
  };
  time[kMax-1] = studyDuration1;

  NumericVector nsubjects = accrual(time, accrualTime, accrualIntensity,
                                    accrualDuration);

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
  double expectedNumberOfSubjects = sum(ptotal*nsubjects);
  double expectedStudyDuration = sum(ptotal*time);
  double expectedInformation = sum(ptotal*I);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  rdu = rmstDiffH0 + criticalValues1/sqrt(I);
  rdl = rmstDiffH0 + futilityBounds1/sqrt(I);

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      rdu[i] = NA_REAL;
      efficacyStopping1[i] = 0;
    }

    if (futilityBounds1[i] == -6) {
      rdl[i] = NA_REAL;
      futilityStopping1[i] = 0;
    }
  }

  double drift = (rmstDiff - rmstDiffH0)*sqrt(maxInformation);
  double driftf = R::qnorm(1-alpha1, 0, 1, 1, 0) +
    R::qnorm(overallReject, 0, 1, 1, 0);
  double inflationFactor = pow(drift/driftf, 2);


  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["futilityBounds"] = futilityBounds1,
    _["rejectPerStage"] = pu,
    _["futilityPerStage"] = pl,
    _["cumulativeRejection"] = cpu,
    _["cumulativeFutility"] = cpl,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["numberOfSubjects"] = nsubjects,
    _["analysisTime"] = time,
    _["efficacyRmstDiff"] = rdu,
    _["futilityRmstDiff"] = rdl,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["information"] = I,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["drift"] = drift,
    _["inflationFactor"] = inflationFactor,
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["kMax"] = kMax,
    _["milestone"] = milestone,
    _["rmstDiffH0"] = rmstDiffH0,
    _["rmst1"] = rmst1,
    _["rmst2"] = rmst2,
    _["rmstDiff"] = rmstDiff);

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
    _["spendingTime"] = spendingTime);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "rmpower";

  return result;
}


//' @title Sample size for difference in restricted mean survival times
//' @description Obtains the needed accrual duration given power,
//' accrual intensity, and follow-up time, the needed follow-up time
//' given power, accrual intensity, and accrual duration, or the needed
//' absolute accrual intensity given power, relative accrual intensity,
//' accrual duration, and follow-up time in a two-group survival design.
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
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffH0 The difference in restricted mean survival times
//'   under the null hypothesis. Defaults to 0 for superiority test.
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
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupTime, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{rmpower} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{rmpower} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmpower}}
//'
//' @examples
//' # Example 1: Obtains follow-up time given power, accrual intensity,
//' # and accrual duration for variable follow-up. Of note, the power
//' # reaches the maximum when the follow-up time equals milestone.
//'
//' rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'              followupTime = NA, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains accrual intensity given power, accrual duration, and
//' # follow-up time for variable follow-up
//'
//' rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains accrual duration given power, accrual intensity, and
//' # follow-up time for fixed follow-up
//'
//' rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = NA,
//'              followupTime = 18, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List rmsamplesize(const double beta = 0.2,
                  const int kMax = 1,
                  const NumericVector& informationRates = NA_REAL,
                  const LogicalVector& efficacyStopping = NA_LOGICAL,
                  const LogicalVector& futilityStopping = NA_LOGICAL,
                  const NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.025,
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const NumericVector& futilityBounds = NA_REAL,
                  const std::string typeBetaSpending = "none",
                  const double parameterBetaSpending = NA_REAL,
                  const NumericVector& userBetaSpending = NA_REAL,
                  const double milestone = NA_REAL,
                  const double rmstDiffH0 = 0,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = NA_REAL,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const NumericVector& stratumFraction = 1,
                  const NumericVector& lambda1 = NA_REAL,
                  const NumericVector& lambda2 = NA_REAL,
                  const NumericVector& gamma1 = 0,
                  const NumericVector& gamma2 = 0,
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

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi);


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
    informationRates1 = NumericVector(tem)/(kMax+0.0);
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

  if (R_isnancpp(milestone)) {
    stop("milestone must be provided");
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

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
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

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
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

  if (!R_isnancpp(accrualDuration) && !R_isnancpp(followupTime) &&
      (milestone >= accrualDuration + followupTime)) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
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


  NumericVector rmsts1(nstrata), rmsts2(nstrata);
  IntegerVector l1 = Range(0, nintervals-1);
  for (int h=0; h<nstrata; h++) {
    l = h*nintervals + l1;
    NumericVector lam1 = lambda1x[l];
    NumericVector lam2 = lambda2x[l];
    rmsts1[h] = rmst(0, milestone, piecewiseSurvivalTime, lam1);
    rmsts2[h] = rmst(0, milestone, piecewiseSurvivalTime, lam2);
  }

  double rmst1 = sum(stratumFraction*rmsts1);
  double rmst2 = sum(stratumFraction*rmsts2);
  double theta1 = rmst1 - rmst2 - rmstDiffH0;
  NumericVector theta(kMax, theta1);

  List design = getDesign(
    beta, NA_REAL, theta1, kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, asf, asfpar, userAlphaSpending, futilityBounds1,
    bsf, bsfpar, userBetaSpending, spendingTime1, 1);

  DataFrame byStageResults = DataFrame(design["byStageResults"]);
  futilityBounds1 = byStageResults["futilityBounds"];

  DataFrame overallResults = DataFrame(design["overallResults"]);
  double maxInformation = overallResults["information"];
  double studyDuration;

  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
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
              DataFrame rm = rmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                dur1, dur2, fixedFollowup);
              return sum(NumericVector(rm[9])) - maxInformation;
            };

  if (unknown == "accrualDuration") {
    double lower = std::max(milestone - followupTime, 0.0) + 1.0e-6;
    accrualDuration = brent(f, lower, interval[1], 1.0e-6);
    studyDuration = accrualDuration + followupTime;
  } else if (unknown == "followupTime") {
    if (f(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    aval, followupTime, fixedFollowup);
                  return sum(NumericVector(rm[9])) - maxInformation;
                };

      accrualDuration = brent(g, accrualDuration, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (!fixedFollowup) { // adjust follow-up time
      double lower =  std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(f, lower, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // fixed follow-up
      followupTime = milestone;
      // adjust study duration to obtain the target max information
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, accrualDuration + aval);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup);
                  return sum(NumericVector(rm[9])) - maxInformation;
                };

      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      double aval = brent(g, lower, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    }
  } else if (unknown == "accrualIntensity") {
    double aval = brent(f, interval[0], interval[1], 1.0e-6);
    accrualIntensity1 = aval*accrualIntensity;
    studyDuration = accrualDuration + followupTime;
  }


  // output the results
  List resultH1, resultH0, result;

  if (rounding) {
    NumericVector u0(1, studyDuration);
    double n0 = accrual(u0, accrualTime, accrualIntensity1,
                        accrualDuration)[0];
    double n = std::ceil(n0);

    if (n - n0 > 1e-6) {
      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      if (!fixedFollowup) {
        // adjust follow-up time to obtain the target maximum information
        auto h = [milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval, fixedFollowup);
                    return sum(NumericVector(rm[9])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // the follow-up time cannot be adjusted for fixed follow-up
        // adjust study duration to obtain the target maximum information
        auto h = [milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup);
                    return sum(NumericVector(rm[9])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        double aval = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }


  resultH1 = rmpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, rmstDiffH0,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  overallResults = DataFrame(resultH1["overallResults"]);
  maxInformation = overallResults["information"];


  // obtain results under H0 by matching the maximum information
  // first find the hazard rate for the treatment group that yields
  // the specified difference in restricted mean survival times under H0
  auto frmst = [milestone, piecewiseSurvivalTime, stratumFraction,
                nintervals, nstrata, l1, lambda2x, rmst2,
                rmstDiffH0](double aval)-> double {
                  NumericVector rmsts1(nstrata);
                  for (int h=0; h<nstrata; h++) {
                    IntegerVector l = h*nintervals + l1;
                    NumericVector lam2 = lambda2x[l];
                    rmsts1[h] = rmst(0, milestone, piecewiseSurvivalTime,
                                     aval*lam2);
                  }
                  double rmst1 = sum(stratumFraction*rmsts1);
                  return rmst1 - rmst2 - rmstDiffH0;
                };

  double aval = brent(frmst, 0.1, 1.9, 1.0e-6);
  NumericVector lambda1H0 = aval*lambda2;

  if (!fixedFollowup) {
    auto h = [milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambda1H0, lambda2, gamma1, gamma2,
              accrualDuration, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                DataFrame rm = rmstat(
                  u0, milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1H0, lambda2, gamma1, gamma2,
                  accrualDuration, aval, fixedFollowup);
                return sum(NumericVector(rm[9])) - maxInformation;
              };

    if (h(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1H0, lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1H0, lambda2, gamma1, gamma2,
                    aval, followupTime, fixedFollowup);
                  return sum(NumericVector(rm[9])) - maxInformation;
                };

      double lower = accrualDuration, upper = 2.0*lower;
      while (g(upper) <= 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(g, lower, upper, 1.0e-6);
    } else if ((accrualDuration <= milestone) || (h(0) < 0)) {
      // adjust follow-up time
      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(h, lower, milestone, 1.0e-6);
    } else {
      // adjust accrual duration
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1H0, lambda2, gamma1, gamma2,
                fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1H0, lambda2, gamma1, gamma2,
                    aval, 0, fixedFollowup);
                  return sum(NumericVector(rm[9])) - maxInformation;
                };

      double lower = milestone + 1.0e-6;
      accrualDuration = brent(g, lower, accrualDuration, 1.0e-6);
      followupTime = 0.0;
    }
    studyDuration = accrualDuration + followupTime;
  } else { // fixed follow-up
    // cannot adjust the followupTime
    auto h = [milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambda1H0, lambda2, gamma1, gamma2,
              followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, aval + followupTime);
                DataFrame rm = rmstat(
                  u0, milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1H0, lambda2, gamma1, gamma2,
                  aval, followupTime, fixedFollowup);
                return sum(NumericVector(rm[9])) - maxInformation;
              };

    if (h(accrualDuration) < 0) { // increase accrual duration
      double lower = accrualDuration, upper = 2.0*lower;
      while (h(upper) <= 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(h, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // decrease study duration
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1H0, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, accrualDuration + aval);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1H0, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup);
                  return sum(NumericVector(rm[9])) - maxInformation;
                };

      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      double aval = brent(g, lower, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    }
  }


  // use the same stopping boundaries as under H1
  resultH0 = rmpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, rmstDiffH0,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda1H0, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}



//' @title Power for one-sample restricted mean survival time
//' @description Estimates the power, stopping probabilities, and expected
//' sample size in a one-group survival design.
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
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstH0 The restricted mean survival time under the null
//'   hypothesis.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param lambda A vector of hazard rates for the event in each analysis
//'  time interval by stratum under the alternative hypothesis.
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
//' @return An S3 class \code{rmpower1s} object with 3 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{drift}: The drift parameter, equal to
//'       \code{(rmst - rmstH0)*sqrt(information)}.
//'
//'     - \code{inflationFactor}: The inflation factor (relative to the
//'       fixed design).
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
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
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{milestone}: The milestone time relative to randomization.
//'
//'     - \code{rmstH0}: The restricted mean survival time under the null
//'       hypothesis.
//'
//'     - \code{rmst}: The restricted mean survival time under the
//'       alternative hypothesis.
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
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRmst}: The efficacy boundaries on the restricted
//'       mean survival time.
//'
//'     - \code{futilityRmst}: The futility boundaries on the restricted
//'       mean survival time.
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
//'   \code{stratumFraction}, \code{lambda}, \code{gamma},
//'   and \code{spendingTime}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmstat}}
//'
//' @examples
//'
//' rmpower1s(kMax = 2, informationRates = c(0.8, 1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           milestone = 18, rmstH0 = 10,
//'           accrualTime = seq(0, 8),
//'           accrualIntensity = 26/9*seq(1, 9),
//'           piecewiseSurvivalTime = c(0, 6),
//'           stratumFraction = c(0.2, 0.8),
//'           lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'           gamma = -log(1-0.05)/12, accrualDuration = 22,
//'           followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
//'
// [[Rcpp::export]]
List rmpower1s(const int kMax = 1,
               const NumericVector& informationRates = NA_REAL,
               const LogicalVector& efficacyStopping = NA_LOGICAL,
               const LogicalVector& futilityStopping = NA_LOGICAL,
               const NumericVector& criticalValues = NA_REAL,
               const double alpha = 0.025,
               const std::string typeAlphaSpending = "sfOF",
               const double parameterAlphaSpending = NA_REAL,
               const NumericVector& userAlphaSpending = NA_REAL,
               const NumericVector& futilityBounds = NA_REAL,
               const std::string typeBetaSpending = "none",
               const double parameterBetaSpending = NA_REAL,
               const double milestone = NA_REAL,
               const double rmstH0 = NA_REAL,
               const NumericVector& accrualTime = 0,
               const NumericVector& accrualIntensity = NA_REAL,
               const NumericVector& piecewiseSurvivalTime = 0,
               const NumericVector& stratumFraction = 1,
               const NumericVector& lambda = NA_REAL,
               const NumericVector& gamma = 0,
               const double accrualDuration = NA_REAL,
               const double followupTime = NA_REAL,
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
    informationRates1 = NumericVector(tem)/(kMax+0.0);
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

  if (R_isnancpp(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (R_isnancpp(rmstH0)) {
    stop("rmstH0 must be provided");
  }

  if (rmstH0 <= 0) {
    stop("rmstH0 must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
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

  if (is_true(any(is_na(lambda)))) {
    stop("lambda must be provided");
  }

  if (is_true(any(lambda < 0))) {
    stop("lambda must be non-negative");
  }

  if (is_true(any(gamma < 0))) {
    stop("gamma must be non-negative");
  }

  if (lambda.size() != 1 && lambda.size() != nintervals &&
      lambda.size() != nsi) {
    stop("Invalid length for lambda");
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

  if (milestone >= accrualDuration + followupTime) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      (milestone >= studyDuration)) {
    stop("milestone cannot exceed studyDuration for fixed follow-up");
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


  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || R_isnancpp(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }


  // obtain the timing of interim analysis using the twin treatment group
  DataFrame rm;
  NumericVector time(kMax), rmstu(kMax), rmstl(kMax);
  NumericVector u0(1, studyDuration1);
  rm = rmstat(u0, milestone, 1,
              accrualTime, 2.0*accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda, lambda, gamma, gamma,
              accrualDuration, followupTime, fixedFollowup);

  double maxInformation = 2.0*sum(NumericVector(rm[9]));
  double rmstH1 = sum(NumericVector(rm[3]));
  NumericVector theta = rep(rmstH1 - rmstH0, kMax);
  NumericVector I = maxInformation*informationRates1;

  double information1;
  auto f = [milestone, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda, gamma,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double aval)->double {
              NumericVector u0(1, aval);
              DataFrame rm = rmstat(
                u0, milestone, 1,
                accrualTime, 2.0*accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda, lambda, gamma, gamma,
                accrualDuration, followupTime, fixedFollowup);
              return 2.0*sum(NumericVector(rm[9])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    // match the predicted information to the target
    information1 = std::max(I[i], 0.0);
    time[i] = brent(f, milestone + 1.0e-6, studyDuration1, 1.0e-6);
  }
  time[kMax-1] = studyDuration1;

  NumericVector nsubjects = accrual(time, accrualTime, accrualIntensity,
                                    accrualDuration);

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
  double expectedNumberOfSubjects = sum(ptotal*nsubjects);
  double expectedStudyDuration = sum(ptotal*time);
  double expectedInformation = sum(ptotal*I);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  rmstu = rmstH0 + criticalValues1/sqrt(I);
  rmstl = rmstH0 + futilityBounds1/sqrt(I);

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      rmstu[i] = NA_REAL;
      efficacyStopping1[i] = 0;
    }

    if (futilityBounds1[i] == -6) {
      rmstl[i] = NA_REAL;
      futilityStopping1[i] = 0;
    }
  }

  double drift = (rmstH1 - rmstH0)*sqrt(maxInformation);
  double driftf = R::qnorm(1-alpha1, 0, 1, 1, 0) +
    R::qnorm(overallReject, 0, 1, 1, 0);
  double inflationFactor = pow(drift/driftf, 2);

  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["futilityBounds"] = futilityBounds1,
    _["rejectPerStage"] = pu,
    _["futilityPerStage"] = pl,
    _["cumulativeRejection"] = cpu,
    _["cumulativeFutility"] = cpl,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["numberOfSubjects"] = nsubjects,
    _["analysisTime"] = time,
    _["efficacyRmst"] = rmstu,
    _["futilityRmst"] = rmstl,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["information"] = I,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["drift"] = drift,
    _["inflationFactor"] = inflationFactor,
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["kMax"] = kMax,
    _["milestone"] = milestone,
    _["rmstH0"] = rmstH0,
    _["rmst"] = rmstH1);

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
    _["lambda"] = lambda,
    _["gamma"] = gamma,
    _["spendingTime"] = spendingTime);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "rmpower1s";

  return result;
}


//' @title Sample size for one-sample restricted mean survival time
//' @description Obtains the needed accrual duration given power and
//' follow-up time, the needed follow-up time given power and
//' accrual duration, or the needed absolute accrual rates given
//' power, accrual duration, follow-up duration, and relative accrual
//' rates in a one-group survival design.
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
//' @param milestone The milestone time at which to calculate the
//'   restricted survival time.
//' @param rmstH0 The restricted mean survival time under the null
//'   hypothesis.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param lambda A vector of hazard rates for the event in each analysis
//'  time interval by stratum under the alternative hypothesis.
//' @param gamma The hazard rate for exponential dropout or a vector of
//'   hazard rates for piecewise exponential dropout. Defaults to 0 for
//'   no dropout.
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{rmpower1s} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{rmpower1s} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmpower1s}}
//'
//' @examples
//' # Example 1: Obtains follow-up duration given power, accrual intensity,
//' # and accrual duration for variable follow-up
//'
//' rmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, rmstH0 = 10,
//'                accrualTime = seq(0, 8),
//'                accrualIntensity = 26/9*seq(1, 9),
//'                piecewiseSurvivalTime = c(0, 6),
//'                stratumFraction = c(0.2, 0.8),
//'                lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'                gamma = -log(1-0.05)/12, accrualDuration = 22,
//'                followupTime = NA, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains accrual intensity given power, accrual duration, and
//' # follow-up duration for variable follow-up
//'
//' rmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, rmstH0 = 10,
//'                accrualTime = seq(0, 8),
//'                accrualIntensity = 26/9*seq(1, 9),
//'                piecewiseSurvivalTime = c(0, 6),
//'                stratumFraction = c(0.2, 0.8),
//'                lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'                gamma = -log(1-0.05)/12, accrualDuration = 22,
//'                followupTime = 18, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains accrual duration given power, accrual intensity, and
//' # follow-up duration for fixed follow-up
//'
//' rmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, rmstH0 = 10,
//'                accrualTime = seq(0, 8),
//'                accrualIntensity = 26/9*seq(1, 9),
//'                piecewiseSurvivalTime = c(0, 6),
//'                stratumFraction = c(0.2, 0.8),
//'                lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'                gamma = -log(1-0.05)/12, accrualDuration = NA,
//'                followupTime = 18, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
List rmsamplesize1s(const double beta = 0.2,
                    const int kMax = 1,
                    const NumericVector& informationRates = NA_REAL,
                    const LogicalVector& efficacyStopping = NA_LOGICAL,
                    const LogicalVector& futilityStopping = NA_LOGICAL,
                    const NumericVector& criticalValues = NA_REAL,
                    const double alpha = 0.025,
                    const std::string typeAlphaSpending = "sfOF",
                    const double parameterAlphaSpending = NA_REAL,
                    const NumericVector& userAlphaSpending = NA_REAL,
                    const NumericVector& futilityBounds = NA_REAL,
                    const std::string typeBetaSpending = "none",
                    const double parameterBetaSpending = NA_REAL,
                    const NumericVector& userBetaSpending = NA_REAL,
                    const double milestone = NA_REAL,
                    const double rmstH0 = NA_REAL,
                    const NumericVector& accrualTime = 0,
                    const NumericVector& accrualIntensity = NA_REAL,
                    const NumericVector& piecewiseSurvivalTime = 0,
                    const NumericVector& stratumFraction = 1,
                    const NumericVector& lambda = NA_REAL,
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

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  NumericVector lambdax(nsi);


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
    informationRates1 = NumericVector(tem)/(kMax+0.0);
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

  if (R_isnancpp(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (R_isnancpp(rmstH0)) {
    stop("rmstH0 must be provided");
  }

  if (rmstH0 <= 0) {
    stop("rmstH0 must be positive");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
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

  if (is_true(any(is_na(lambda)))) {
    stop("lambda must be provided");
  }

  if (is_true(any(lambda < 0))) {
    stop("lambda must be non-negative");
  }

  if (is_true(any(gamma < 0))) {
    stop("gamma must be non-negative");
  }

  if (lambda.size() == 1) {
    lambdax = rep(lambda, nsi);
  } else if (lambda.size() == nintervals) {
    lambdax = rep(lambda, nstrata);
  } else if (lambda.size() == nsi) {
    lambdax = lambda;
  } else {
    stop("Invalid length for lambda");
  }

  if (lambda.size() != 1 && lambda.size() != nintervals &&
      lambda.size() != nsi) {
    stop("Invalid length for lambda");
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

  if (!R_isnancpp(accrualDuration) && !R_isnancpp(followupTime) &&
      (milestone >= accrualDuration + followupTime)) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
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


  NumericVector rmsts(nstrata);
  IntegerVector l1 = Range(0, nintervals-1);
  for (int h=0; h<nstrata; h++) {
    l = h*nintervals + l1;
    NumericVector lam = lambdax[l];
    rmsts[h] = rmst(0, milestone, piecewiseSurvivalTime, lam);
  }

  double rmstH1 = sum(stratumFraction*rmsts);
  double theta1 = rmstH1 - rmstH0;
  NumericVector theta(kMax, theta1);

  List design = getDesign(
    beta, NA_REAL, theta1, kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, asf, asfpar, userAlphaSpending, futilityBounds1,
    bsf, bsfpar, userBetaSpending, spendingTime1, 1);

  DataFrame byStageResults = DataFrame(design["byStageResults"]);
  futilityBounds1 = byStageResults["futilityBounds"];

  DataFrame overallResults = DataFrame(design["overallResults"]);
  double maxInformation = overallResults["information"];
  double studyDuration;

  auto f = [milestone, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda, gamma,
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
              DataFrame rm = rmstat(
                u0, milestone, 1,
                accrualTime, 2.0*accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda, lambda, gamma, gamma,
                dur1, dur2, fixedFollowup);
              return 2.0*sum(NumericVector(rm[9])) - maxInformation;
            };

  if (unknown == "accrualDuration") {
    double lower = std::max(milestone - followupTime, 0.0) + 1.0e-6;
    accrualDuration = brent(f, lower, interval[1], 1.0e-6);
    studyDuration = accrualDuration + followupTime;
  } else if (unknown == "followupTime") {
    if (f(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda, gamma,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda, lambda, gamma, gamma,
                    aval, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(rm[9])) - maxInformation;
                };

      accrualDuration = brent(g, accrualDuration, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (!fixedFollowup) { // adjust follow-up time
      double lower =  std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(f, lower, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // fixed follow-up
      // adjust study duration to obtain the target max information
      followupTime = milestone;
      auto g = [milestone, accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda, gamma,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, accrualDuration + aval);
                  DataFrame rm = rmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda, lambda, gamma, gamma,
                    accrualDuration, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(rm[9])) - maxInformation;
                };

      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      double aval = brent(g, lower, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    }
  } else if (unknown == "accrualIntensity") {
    double aval = brent(f, interval[0], interval[1], 1.0e-6);
    accrualIntensity1 = aval*accrualIntensity;
    studyDuration = accrualDuration + followupTime;
  }


  // output the results
  List resultH1, resultH0, result;

  if (rounding) {
    NumericVector u0(1, studyDuration);
    double n0 = accrual(u0, accrualTime, accrualIntensity1,
                        accrualDuration)[0];
    double n = std::ceil(n0);

    if (n - n0 > 1e-6) {
      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      if (!fixedFollowup) {
        // adjust follow-up time to obtain the target maximum information
        auto h = [milestone, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda, gamma,
                  accrualDuration, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, 1,
                      accrualTime, 2.0*accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda, lambda, gamma, gamma,
                      accrualDuration, aval, fixedFollowup);
                    return 2.0*sum(NumericVector(rm[9])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // the follow-up time cannot be adjusted for fixed follow-up
        // adjust study duration to obtain the target maximum information
        auto h = [milestone, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda, gamma,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, 1,
                      accrualTime, 2.0*accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda, lambda, gamma, gamma,
                      accrualDuration, followupTime, fixedFollowup);
                    return 2.0*sum(NumericVector(rm[9])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        double aval = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }


  resultH1 = rmpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, rmstH0,
    accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  overallResults = DataFrame(resultH1["overallResults"]);
  maxInformation = overallResults["information"];


  // obtain results under H0 by matching the maximum information
  // first find the hazard rate that yields
  // the specified milestone survival probability under H0
  auto fsurv = [milestone, piecewiseSurvivalTime, stratumFraction,
                nintervals, nstrata, l1, lambdax,
                rmstH0](double aval)-> double {
                  NumericVector rmsts(nstrata);
                  for (int h=0; h<nstrata; h++) {
                    IntegerVector l = h*nintervals + l1;
                    NumericVector lam = lambdax[l];
                    rmsts[h] = rmst(0, milestone, piecewiseSurvivalTime,
                                    aval*lam);
                  }
                  double rmstH1 = sum(stratumFraction*rmsts);
                  return rmstH1 - rmstH0;
                };

  double aval = brent(fsurv, 0.1, 1.9, 1.0e-6);
  NumericVector lambdaH0 = aval*lambda;

  if (!fixedFollowup) {
    auto h = [milestone, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambdaH0, gamma,
              accrualDuration, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                DataFrame rm = rmstat(
                  u0, milestone, 1,
                  accrualTime, 2.0*accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambdaH0, lambdaH0, gamma, gamma,
                  accrualDuration, aval, fixedFollowup);
                return 2.0*sum(NumericVector(rm[9])) - maxInformation;
              };

    if (h(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambdaH0, gamma,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambdaH0, lambdaH0, gamma, gamma,
                    aval, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(rm[9])) - maxInformation;
                };

      double lower = accrualDuration, upper = 2.0*lower;
      while (g(upper) <= 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(g, lower, upper, 1.0e-6);
    } else if ((accrualDuration <= milestone) || (h(0) < 0)) {
      // adjust follow-up time
      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(h, lower, milestone, 1.0e-6);
    } else {
      // adjust accrual duration
      auto g = [milestone, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambdaH0, gamma,
                fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  DataFrame rm = rmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambdaH0, lambdaH0, gamma, gamma,
                    aval, 0, fixedFollowup);
                  return 2.0*sum(NumericVector(rm[9])) - maxInformation;
                };

      double lower = milestone + 1.0e-6;
      accrualDuration = brent(g, lower, accrualDuration, 1.0e-6);
      followupTime = 0.0;
    }
    studyDuration = accrualDuration + followupTime;
  } else { // fixed follow-up
    // cannot adjust the followupTime
    auto h = [milestone, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambdaH0, gamma,
              followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, aval + followupTime);
                DataFrame rm = rmstat(
                  u0, milestone, 1,
                  accrualTime, 2.0*accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambdaH0, lambdaH0, gamma, gamma,
                  aval, followupTime, fixedFollowup);
                return 2.0*sum(NumericVector(rm[9])) - maxInformation;
              };

    if (h(accrualDuration) < 0) { // increase accrual duration
      double lower = accrualDuration, upper = 2.0*lower;
      while (h(upper) <= 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(h, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // decrease study duration
      auto g = [milestone, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambdaH0, gamma,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, accrualDuration + aval);
                  DataFrame rm = rmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambdaH0, lambdaH0, gamma, gamma,
                    accrualDuration, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(rm[9])) - maxInformation;
                };

      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      double aval = brent(g, lower, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    }
  }


  // use the same stopping boundaries as under H1
  resultH0 = rmpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, rmstH0,
    accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambdaH0, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}



//' @title Power for equivalence in restricted mean survival time difference
//' @description Obtains the power for equivalence in restricted mean
//' survival time difference.
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
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffLower The lower equivalence limit of restricted mean
//'   survival time difference.
//' @param rmstDiffUpper The upper equivalence limit of restricted mean
//'   survival time difference.
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
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{rmpowerequiv} object with 3 components:
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
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{milestone}: The milestone time relative to randomization.
//'
//'     - \code{rmstDiffLower}: The lower equivalence limit of restricted
//'       mean survival time difference.
//'
//'     - \code{rmstDiffUpper}: The upper equivalence limit of restricted
//'       mean survival time difference.
//'
//'     - \code{rmst1}: The restricted mean survival time for the
//'       treatment group.
//'
//'     - \code{rmst2}: The restricted mean survival time for the
//'       control group.
//'
//'     - \code{rmstDiff}: The restricted mean survival time difference.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
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
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRmstDiffLower}: The efficacy boundaries on the
//'       restricted mean survival time difference scale for the one-sided
//'       null hypothesis at the lower equivalence limit.
//'
//'     - \code{efficacyRmstDiffUpper}: The efficacy boundaries on the
//'       restricted mean survival time difference scale for the one-sided
//'       null hypothesis at the upper equivalence limit.
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
//'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   and \code{spendingTime}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmstat}}
//'
//' @examples
//'
//' rmpowerequiv(kMax = 2, informationRates = c(0.5, 1),
//'              alpha = 0.05, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              rmstDiffLower = -2, rmstDiffUpper = 2,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
List rmpowerequiv(const int kMax = 1,
                  const NumericVector& informationRates = NA_REAL,
                  const NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.05,
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const double milestone = NA_REAL,
                  const double rmstDiffLower = NA_REAL,
                  const double rmstDiffUpper = NA_REAL,
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
                  const NumericVector& spendingTime = NA_REAL,
                  const double studyDuration = NA_REAL) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

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
    informationRates1 = NumericVector(tem)/(kMax+0.0);
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

  if (R_isnancpp(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (R_isnancpp(rmstDiffLower)) {
    stop("rmstDiffLower must be provided");
  }

  if (R_isnancpp(rmstDiffUpper)) {
    stop("rmstDiffUpper must be provided");
  }

  if (rmstDiffLower >= rmstDiffUpper) {
    stop("rmstDiffLower must be less than rmstDiffUpper");
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

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
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

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
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

  if (milestone >= accrualDuration + followupTime) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      (milestone >= studyDuration)) {
    stop("milestone cannot exceed studyDuration for fixed follow-up");
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

  // obtain the study duration
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || R_isnancpp(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }

  // obtain the timing of interim analysis
  NumericVector time(kMax);
  NumericVector u0(1, studyDuration1);
  DataFrame rm = rmstat(u0, milestone, allocationRatioPlanned,
                        accrualTime, accrualIntensity,
                        piecewiseSurvivalTime, stratumFraction,
                        lambda1, lambda2, gamma1, gamma2,
                        accrualDuration, followupTime, fixedFollowup);

  double rmst1 = sum(NumericVector(rm[3]));
  double rmst2 = sum(NumericVector(rm[4]));
  double rmstDiff = sum(NumericVector(rm[5]));
  NumericVector theta(kMax, rmstDiff);
  double maxInformation = sum(NumericVector(rm[9]));
  NumericVector I = maxInformation*informationRates1;
  double information1;

  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double aval)->double {
              NumericVector u0(1, aval);
              DataFrame rm = rmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup);
              return sum(NumericVector(rm[9])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    information1 = I[i];
    time[i] = brent(f, milestone + 1.0e-6, studyDuration1, 1.0e-6);
  };
  time[kMax-1] = studyDuration1;

  // calculate cumulative rejection probability under H1
  double theta10 = rmstDiffLower, theta20 = rmstDiffUpper;
  NumericVector b = criticalValues1;
  NumericVector l = b + theta10*sqrt(I);
  NumericVector u = -b + theta20*sqrt(I);

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

  NumericVector efficacyRmstDiffLower = theta10 + b/sqrt(I);
  NumericVector efficacyRmstDiffUpper = theta20 - b/sqrt(I);

  // calculate cumulative rejection under H10
  NumericVector bH10 = -b + (theta20 - theta10)*sqrt(I);
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
  NumericVector bH20 = b + (theta10 - theta20)*sqrt(I);
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

  NumericVector nsubjects = accrual(time, accrualTime, accrualIntensity,
                                    accrualDuration);

  double overallReject = cp[kMax-1];
  double attainedAlphaH10 = cpH10[kMax-1];
  double attainedAlphaH20 = cpH20[kMax-1];
  double expectedNumberOfSubjects = sum(q*nsubjects);
  double expectedStudyDuration = sum(q*time);
  double expectedInformation = sum(q*I);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = alpha,
    _["attainedAlphaH10"] = attainedAlphaH10,
    _["attainedAlphaH20"] = attainedAlphaH20,
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["kMax"] = kMax,
    _["milestone"] = milestone,
    _["rmstDiffLower"] = rmstDiffLower,
    _["rmstDiffUpper"] = rmstDiffUpper,
    _["rmst1"] = rmst1,
    _["rmst2"] = rmst2,
    _["rmstDiff"] = rmstDiff,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup);

  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["rejectPerStage"] = rejectPerStage,
    _["cumulativeRejection"] = cp,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["cumulativeAttainedAlphaH10"] = cpH10,
    _["cumulativeAttainedAlphaH20"] = cpH20,
    _["numberOfSubjects"] = nsubjects,
    _["analysisTime"] = time,
    _["efficacyRmstDiffLower"] = efficacyRmstDiffLower,
    _["efficacyRmstDiffUpper"] = efficacyRmstDiffUpper,
    _["efficacyP"] = efficacyP,
    _["information"] = I);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["allocationRatioPlanned"] = allocationRatioPlanned,
    _["accrualTime"] = accrualTime,
    _["accrualIntensity"] = accrualIntensity,
    _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
    _["stratumFraction"] = stratumFraction,
    _["lambda1"] = lambda1,
    _["lambda2"] = lambda2,
    _["gamma1"] = gamma1,
    _["gamma2"] = gamma2,
    _["spendingTime"] = spendingTime);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "rmpowerequiv";

  return result;
}


//' @title Sample size for equivalence in restricted mean survival time
//' difference
//' @description Obtains the sample size for equivalence in restricted
//' mean survival time difference.
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
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffLower The lower equivalence limit of restricted mean
//'   survival time difference.
//' @param rmstDiffUpper The upper equivalence limit of restricted mean
//'   survival time difference.
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
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return An S3 class \code{rmpowerequiv} object
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmpowerequiv}}
//'
//' @examples
//'
//' rmsamplesizeequiv(beta = 0.1, kMax = 2, informationRates = c(0.5, 1),
//'                   alpha = 0.05, typeAlphaSpending = "sfOF",
//'                   milestone = 18,
//'                   rmstDiffLower = -2, rmstDiffUpper = 2,
//'                   allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'                   accrualIntensity = 26/9*seq(1, 9),
//'                   piecewiseSurvivalTime = c(0, 6),
//'                   stratumFraction = c(0.2, 0.8),
//'                   lambda1 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'                   gamma1 = -log(1-0.05)/12,
//'                   gamma2 = -log(1-0.05)/12, accrualDuration = NA,
//'                   followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
List rmsamplesizeequiv(const double beta = 0.2,
                       const int kMax = 1,
                       const NumericVector& informationRates = NA_REAL,
                       const NumericVector& criticalValues = NA_REAL,
                       const double alpha = 0.05,
                       const std::string typeAlphaSpending = "sfOF",
                       const double parameterAlphaSpending = NA_REAL,
                       const NumericVector& userAlphaSpending = NA_REAL,
                       const double milestone = NA_REAL,
                       const double rmstDiffLower = NA_REAL,
                       const double rmstDiffUpper = NA_REAL,
                       const double allocationRatioPlanned = 1,
                       const NumericVector& accrualTime = 0,
                       const NumericVector& accrualIntensity = NA_REAL,
                       const NumericVector& piecewiseSurvivalTime = 0,
                       const NumericVector& stratumFraction = 1,
                       const NumericVector& lambda1 = NA_REAL,
                       const NumericVector& lambda2 = NA_REAL,
                       const NumericVector& gamma1 = 0,
                       const NumericVector& gamma2 = 0,
                       double accrualDuration = NA_REAL,
                       double followupTime = NA_REAL,
                       const bool fixedFollowup = 0,
                       const NumericVector& interval =
                         NumericVector::create(0.001, 240),
                         const NumericVector& spendingTime = NA_REAL,
                         const bool rounding = 1) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector accrualIntensity1 = clone(accrualIntensity);
  NumericVector spendingTime1 = clone(spendingTime);

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });

  double asfpar = parameterAlphaSpending;

  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi);


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
    informationRates1 = NumericVector(tem)/(kMax+0.0);
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

  if (R_isnancpp(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (R_isnancpp(rmstDiffLower)) {
    stop("rmstDiffLower must be provided");
  }

  if (R_isnancpp(rmstDiffUpper)) {
    stop("rmstDiffUpper must be provided");
  }

  if (rmstDiffLower >= rmstDiffUpper) {
    stop("rmstDiffLower must be less than rmstDiffUpper");

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

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
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

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
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

  NumericVector rmsts1(nstrata), rmsts2(nstrata);
  IntegerVector l1 = Range(0, nintervals-1);
  for (int h=0; h<nstrata; h++) {
    IntegerVector l = h*nintervals + l1;
    NumericVector lam1 = lambda1x[l];
    NumericVector lam2 = lambda2x[l];
    rmsts1[h] = rmst(0, milestone, piecewiseSurvivalTime, lam1);
    rmsts2[h] = rmst(0, milestone, piecewiseSurvivalTime, lam2);
  }

  double rmst1 = sum(stratumFraction*rmsts1);
  double rmst2 = sum(stratumFraction*rmsts2);
  double rmstDiff = rmst1 - rmst2;
  double theta10 = rmstDiffLower, theta20 = rmstDiffUpper;
  NumericVector theta(kMax, rmstDiff);

  List design = getDesignEquiv(
    beta, NA_REAL, theta10, theta20, rmstDiff,
    kMax, informationRates1, criticalValues1,
    alpha, asf, asfpar, userAlphaSpending, spendingTime1,
    1, 1, 1, 1);

  DataFrame overallResults = DataFrame(design["overallResults"]);
  double maxInformation = overallResults["information"];
  double studyDuration;

  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
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
              DataFrame rm = rmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                dur1, dur2, fixedFollowup);
              return sum(NumericVector(rm[9])) - maxInformation;
            };

  if (unknown == "accrualDuration") {
    double lower = std::max(milestone - followupTime, 0.0) + 1.0e-6;
    accrualDuration = brent(f, lower, interval[1], 1.0e-6);
    studyDuration = accrualDuration + followupTime;
  } else if (unknown == "followupTime") {
    if (f(milestone) < 0) {
      std::string str1 = "NOTE: The required information cannot be ";
      std::string str2 = "attained by increasing followupTime alone.";
      std::string str3 = "NOTE: accrualDuration is also increased to ";
      std::string str4 = "attain the required information.";
      Rcout << str1 + str2 << "\n";
      Rcout << str3 + str4 << "\n";

      followupTime = milestone;
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                maxInformation](double aval)-> double{
                  NumericVector u0(1, aval + followupTime);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    aval, followupTime, fixedFollowup);
                  return sum(NumericVector(rm[9])) - maxInformation;
                };

      accrualDuration = brent(g, accrualDuration, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (!fixedFollowup) { // adjust follow-up time
      double lower =  std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      followupTime = brent(f, lower, interval[1], 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // fixed follow-up
      // adjust study duration to obtain the target max information
      followupTime = milestone;
      auto g = [milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, accrualDuration + aval);
                  DataFrame rm = rmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup);
                  return sum(NumericVector(rm[9])) - maxInformation;
                };

      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      double aval = brent(g, lower, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    }
  } else if (unknown == "accrualIntensity") {
    double aval = brent(f, interval[0], interval[1], 1.0e-6);
    accrualIntensity1 = aval*accrualIntensity;
    studyDuration = accrualDuration + followupTime;
  }


  // output the results
  if (rounding) {
    NumericVector u0(1, studyDuration);
    double n0 = accrual(u0, accrualTime, accrualIntensity1,
                        accrualDuration)[0];
    double n = std::ceil(n0);

    if (n - n0 > 1e-6) {
      // adjust accrual intensity or duration to obtain int # of subjects
      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      if (!fixedFollowup) {
        // adjust follow-up time to obtain the target maximum information
        auto h = [milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval, fixedFollowup);
                    return sum(NumericVector(rm[9])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // the follow-up time cannot be adjusted for fixed follow-up
        // adjust study duration to obtain the target maximum information
        auto h = [milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    DataFrame rm = rmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup);
                    return sum(NumericVector(rm[9])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        double aval = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }


  List result = rmpowerequiv(
    kMax, informationRates1, criticalValues1,
    alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, milestone, rmstDiffLower, rmstDiffUpper,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  return result;
}


//' @title Estimate of restricted mean survival time
//' @description Obtains the estimate of restricted means survival time
//' for each stratum.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{rep}: The replication for by-group processing.
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param rep The name of the replication variable in the input data.
//' @param stratum The name of the stratum variable in the input data.
//' @param time The name of the time variable in the input data.
//' @param event The name of the event variable in the input data.
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param confint The level of the two-sided confidence interval for
//'   the survival probabilities. Defaults to 0.95.
//' @param biascorrection Whether to apply bias correction for the
//'   variance estimate. Defaults to no bias correction.
//'
//' @return A data frame with the following variables:
//'
//' * \code{rep}: The replication.
//'
//' * \code{stratum}: The stratum variable.
//'
//' * \code{size}: The number of subjects in the stratum.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{rmst}: The estimate of restricted mean survival time.
//'
//' * \code{stderr}: The standard error of the estimated rmst.
//'
//' * \code{lower}: The lower bound of confidence interval if requested.
//'
//' * \code{upper}: The upper bound of confidence interval if requested.
//'
//' * \code{confint}: The level of confidence interval if requested.
//'
//' * \code{biascorrection}: Whether to apply bias correction for the
//'   variance estimate.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' rmest(data = aml, stratum = "x",
//'       time = "time", event = "status", milestone = 24)
//'
//' @export
// [[Rcpp::export]]
DataFrame rmest(const DataFrame data,
                const std::string rep = "rep",
                const std::string stratum = "stratum",
                const std::string time = "time",
                const std::string event = "event",
                const double milestone = NA_REAL,
                const double confint = 0.95,
                const bool biascorrection = 0) {
  int h, i, j, k, n = data.nrows();

  bool has_rep = hasVariable(data, rep);
  bool has_stratum = hasVariable(data, stratum);
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);

  if (!has_time) {
    stop("data must contain the time variable");
  }

  if (!has_event) {
    stop("data must contain the event variable");
  }

  NumericVector timen = data[time];
  NumericVector eventn = data[event];

  if (is_true(any(timen <= 0))) {
    stop("time must be positive for each subject");
  }

  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }


  // create the numeric rep variable
  IntegerVector repn(n);
  IntegerVector repwi;
  NumericVector repwn;
  CharacterVector repwc;
  if (!has_rep) {
    repn.fill(1);
  } else {
    if (TYPEOF(data[rep]) == INTSXP) {
      IntegerVector repv = data[rep];
      repwi = unique(repv);  // Get unique levels
      std::sort(repwi.begin(), repwi.end());
      repn = match(repv, repwi);  // Map codes to levels
    } else if (TYPEOF(data[rep]) == REALSXP) {
      NumericVector repv = data[rep];
      repwn = unique(repv);
      std::sort(repwn.begin(), repwn.end());
      repn = match(repv, repwn);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      CharacterVector repv = data[rep];
      repwc = unique(repv);
      std::sort(repwc.begin(), repwc.end());
      repn = match(repv, repwc);
    } else {
      stop("incorrect type for the rep variable in the input data");
    }
  }


  // create the numeric stratum variable
  IntegerVector stratumn(n);
  IntegerVector stratumwi;
  NumericVector stratumwn;
  CharacterVector stratumwc;
  if (!has_stratum) {
    stratumn.fill(1);
  } else {
    if (TYPEOF(data[stratum]) == INTSXP) {
      IntegerVector stratumv = data[stratum];
      stratumwi = unique(stratumv);  // Get unique levels
      std::sort(stratumwi.begin(), stratumwi.end());
      stratumn = match(stratumv, stratumwi);  // Map codes to levels
    } else if (TYPEOF(data[stratum]) == REALSXP) {
      NumericVector stratumv = data[stratum];
      stratumwn = unique(stratumv);
      std::sort(stratumwn.begin(), stratumwn.end());
      stratumn = match(stratumv, stratumwn);
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      CharacterVector stratumv = data[stratum];
      stratumwc = unique(stratumv);
      std::sort(stratumwc.begin(), stratumwc.end());
      stratumn = match(stratumv, stratumwc);
    } else {
      stop("incorrect type for the stratum variable in the input data");
    }
  }


  if (R_isnancpp(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (confint <= 0 || confint >= 1) {
    stop("confint must lie between 0 and 1");
  }


  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  repn = repn[order];
  stratumn = stratumn[order];
  timen = timen[order];
  eventn = eventn[order];

  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (i=1; i<n; i++) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }

  int nreps = idx.size();
  idx.push_back(n);

  IntegerVector rep0(n, NA_INTEGER);
  IntegerVector stratum0(n), size0(n);
  NumericVector rmst0(n), stderr0(n);
  NumericVector lower0(n), upper0(n);

  double z = R::qnorm((1.0 + confint)/2.0, 0, 1, 1, 0);

  bool noerr = 1;
  int index = 0;
  for (h=0; h<nreps; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = q1.size();

    IntegerVector stratum1 = stratumn[q1];
    NumericVector time1 = timen[q1];
    NumericVector event1 = eventn[q1];

    // sort by stratum, time, and event with event in descending order
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return (stratum1[i] < stratum1[j]) ||
        ((stratum1[i] == stratum1[j]) && (time1[i] < time1[j])) ||
        ((stratum1[i] == stratum1[j]) && (time1[i] == time1[j]) &&
        (event1[i] > event1[j]));
    });

    stratum1 = stratum1[order1];
    time1 = time1[order1];
    event1 = event1[order1];

    // identify the locations of the unique values of stratum
    IntegerVector idx1(1,0);
    for (i=1; i<n1; i++) {
      if (stratum1[i] != stratum1[i-1]) {
        idx1.push_back(i);
      }
    }

    int nstrata = idx1.size();
    idx1.push_back(n1);

    for (i=0; i<nstrata; i++) {
      IntegerVector q2 = Range(idx1[i], idx1[i+1]-1);
      NumericVector time2 = time1[q2];
      NumericVector event2 = event1[q2];
      int n2 = q2.size();

      if (milestone > max(time2)) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else if (TYPEOF(data[rep]) == INTSXP) {
          reperr = " " + rep + " = " + std::to_string(repwi[repn[idx[h]]-1]);
        } else if (TYPEOF(data[rep]) == REALSXP) {
          reperr = " " + rep + " = " + std::to_string(repwn[repn[idx[h]]-1]);
        } else {
          reperr = " " + rep + " = " + repwc[repn[idx[h]]-1];
        }

        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else if (TYPEOF(data[stratum]) == INTSXP) {
          stratumerr = " " + stratum + " = " +
            std::to_string(stratumwi[stratum1[idx1[i]]-1]);
        } else if (TYPEOF(data[stratum]) == REALSXP) {
          stratumerr = " " + stratum + " = " +
            std::to_string(stratumwn[stratum1[idx1[i]]-1]);
        } else {
          stratumerr = " " + stratum + " = " +
            stratumwc[stratum1[idx1[i]]-1];
        }

        std::string str1 = "The milestone is larger than";
        std::string str2 = "the largest observed time";
        std::string errmsg = str1 + " " + str2;
        if (!reperr.empty() || !stratumerr.empty()) {
          errmsg = errmsg + ":" + reperr + stratumerr;
        }

        if (noerr) {
          Rcout << errmsg << "\n";
          Rcout << "Additional warning messages are suppressed" << "\n";
          noerr = 0;
        }

        continue;
      }

      NumericVector time0(n2, NA_REAL);
      NumericVector nrisk0(n2), nevent0(n2), surv0(n2);
      int index1 = 0;

      double t, nrisk, nevent, surv = 1.0;
      bool cache = 0;
      for (j=0; j<n2; j++) {
        if (((j == 0) && (event2[j] == 1)) ||
            ((j >= 1) && (event2[j] == 1) && (time2[j] > time2[j-1]))) {
          // new event
          // add the info for the previous event
          if (cache) {
            surv = surv*(1.0 - nevent/nrisk);

            time0[index1] = t;
            nrisk0[index1] = nrisk;
            nevent0[index1] = nevent;
            surv0[index1] = surv;

            index1++;
          }

          // update the buffer for the current event time
          t = time2[j];
          nrisk = n2-j;
          nevent = 1;

          cache = 1;
        } else if ((j >= 1) && (event2[j] == 1) && (event2[j-1] == 1) &&
          (time2[j] == time2[j-1])) { // tied event
          nevent = nevent + 1;
        } else if ((j >= 1) && (event2[j] == 0) && (event2[j-1] == 1)) {
          // new censoring
          // add the info for the previous event
          surv = surv*(1.0 - nevent/nrisk);

          time0[index1] = t;
          nrisk0[index1] = nrisk;
          nevent0[index1] = nevent;
          surv0[index1] = surv;

          index1++;

          // empty the cache for the current event time
          cache = 0;
        }
      }

      // add the info for the last event
      if (cache) {
        surv = surv*(1.0 - nevent/nrisk);

        time0[index1] = t;
        nrisk0[index1] = nrisk;
        nevent0[index1] = nevent;
        surv0[index1] = surv;

        index1++;
      }

      // only keep nonmissing records
      int N;
      LogicalVector sub = !is_na(time0);
      if (is_false(any(sub))) { // no event
        N = 0;
        time0 = NumericVector::create(0.0);
        nrisk0 = NumericVector::create(n2);
        nevent0 = NumericVector::create(0.0);
        surv0 = NumericVector::create(1.0);
      } else { // at least 1 event
        time0 = time0[sub];
        nrisk0 = nrisk0[sub];
        nevent0 = nevent0[sub];
        surv0 = surv0[sub];

        // locate the latest event time before milestone
        NumericVector milestone1(1, milestone);
        N = findInterval3(milestone1, time0)[0];

        // prepend time zero information
        time0.push_front(0.0);
        nrisk0.push_front(n2);
        nevent0.push_front(0.0);
        surv0.push_front(1.0);
      }

      // replace the last time of interest with milestone
      if (N == time0.size() - 1) {
        time0.push_back(milestone);
      } else {
        time0[N+1] = milestone;
      }

      // calculate the partial sum of the trapezoid integration
      NumericVector rmstx(N+1);
      rmstx[0] = surv0[0]*(time0[1] - time0[0]);
      for (k=1; k<=N; k++) {
        rmstx[k] = rmstx[k-1] + surv0[k]*(time0[k+1] - time0[k]);
      }

      // calculate rmst and its variance
      double u = rmstx[N];
      double v = 0.0;
      for (k=1; k<=N; k++) {
        // rmst from the kth event time to milestone
        double a = u - rmstx[k-1];
        // do not add variance if the largest observed time is an event time
        if (nrisk0[k] > nevent0[k]) {
          v += nevent0[k]*a*a/(nrisk0[k]*(nrisk0[k] - nevent0[k]));
        }
      }

      // apply bias correction if requested
      if (biascorrection) {
        double m1 = 0;
        for (k=1; k<=N; k++) {
          m1 += nevent0[k];
        }

        if (m1 <= 1.0) {
          std::string reperr;
          if (!has_rep) {
            reperr = "";
          } else if (TYPEOF(data[rep]) == INTSXP) {
            reperr = " " + rep + " = " +
              std::to_string(repwi[repn[idx[h]]-1]);
          } else if (TYPEOF(data[rep]) == REALSXP) {
            reperr = " " + rep + " = " +
              std::to_string(repwn[repn[idx[h]]-1]);
          } else {
            reperr = " " + rep + " = " + repwc[repn[idx[h]]-1];
          }

          std::string stratumerr;
          if (!has_stratum) {
            stratumerr = "";
          } else if (TYPEOF(data[stratum]) == INTSXP) {
            stratumerr = " " + stratum + " = " +
              std::to_string(stratumwi[stratum1[idx1[i]]-1]);
          } else if (TYPEOF(data[stratum]) == REALSXP) {
            stratumerr = " " + stratum + " = " +
              std::to_string(stratumwn[stratum1[idx1[i]]-1]);
          } else {
            stratumerr = " " + stratum + " = " +
              stratumwc[stratum1[idx1[i]]-1];
          }

          std::string str1 = "Bias correction is not done due to no or";
          std::string str2 = "only 1 event before the milestone time:";
          std::string errmsg = str1 + " " + str2;
          if (!reperr.empty() || !stratumerr.empty()) {
            errmsg = errmsg + ":" + reperr + stratumerr;
          }

          if (noerr) {
            Rcout << errmsg << "\n";
            Rcout << "Additional warning messages are suppressed" << "\n";
            noerr = 0;
          }
        } else {
          v = m1/(m1 - 1.0)*v;
        }
      }

      rep0[index] = repn[idx[h]];
      stratum0[index] = stratum1[idx1[i]];
      size0[index] = n2;
      rmst0[index] = u;
      stderr0[index] = sqrt(v);
      lower0[index] = u - z*stderr0[index];
      upper0[index] = u + z*stderr0[index];

      index++;
    }
  }

  // only keep nonmissing records
  LogicalVector sub = !is_na(rep0);
  if (is_false(any(sub))) {
    stop("no replication enables valid inference");
  }

  rep0 = rep0[sub];
  stratum0 = stratum0[sub];
  size0 = size0[sub];
  rmst0 = rmst0[sub];
  stderr0 = stderr0[sub];
  lower0 = lower0[sub];
  upper0 = upper0[sub];


  DataFrame result = DataFrame::create(
      _["size"] = size0,
      _["milestone"] = milestone,
      _["rmst"] = rmst0,
      _["stderr"] = stderr0,
      _["lower"] = lower0,
      _["upper"] = upper0,
      _["confint"] = confint,
      _["biascorrection"] = biascorrection);

  if (has_rep) {
    if (TYPEOF(data[rep]) == INTSXP) {
      result.push_back(repwi[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      result.push_back(repwn[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      result.push_back(repwc[rep0-1], rep);
    }
  }

  if (has_stratum) {
    if (TYPEOF(data[stratum]) == INTSXP) {
      result.push_back(stratumwi[stratum0-1], stratum);
    } else if (TYPEOF(data[stratum]) == REALSXP) {
      result.push_back(stratumwn[stratum0-1], stratum);
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      result.push_back(stratumwc[stratum0-1], stratum);
    }
  }

  return result;
}


//' @title Estimate of restricted mean survival time difference
//' @description Obtains the estimate of restricted mean survival time
//' difference between two treatment groups.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{rep}: The replication for by-group processing.
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{treat}: The treatment.
//'
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param rep The name of the replication variable in the input data.
//' @param stratum The name of the stratum variable in the input data.
//' @param treat The name of the treatment variable in the input data.
//' @param time The name of the time variable in the input data.
//' @param event The name of the event variable in the input data.
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffH0 The difference in restricted mean survival times
//'   under the null hypothesis. Defaults to 0 for superiority test.
//' @param confint The level of the two-sided confidence interval for
//'   the difference in restricted mean survival times. Defaults to 0.95.
//' @param biascorrection Whether to apply bias correction for the
//'   variance estimate of individual restricted mean survival times.
//'   Defaults to no bias correction.
//'
//' @return A data frame with the following variables:
//'
//' * \code{rep}: The replication number.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{rmstDiffH0}: The difference in restricted mean survival times
//'   under the null hypothesis.
//'
//' * \code{rmst1}: The estimated restricted mean survival time for
//'   the treatment group.
//'
//' * \code{rmst2}: The estimated restricted mean survival time for
//'   the control group.
//'
//' * \code{rmstDiff}: The estimated difference in restricted mean
//'   survival times.
//'
//' * \code{vrmst1}: The variance for rmst1.
//'
//' * \code{vrmst2}: The variance for rmst2.
//'
//' * \code{vrmstDiff}: The variance for rmstDiff.
//'
//' * \code{rmstDiffZ}: The Z-statistic value.
//'
//' * \code{rmstDiffPValue}: The one-sided p-value.
//'
//' * \code{lower}: The lower bound of confidence interval.
//'
//' * \code{upper}: The upper bound of confidence interval.
//'
//' * \code{confint}: The level of confidence interval.
//'
//' * \code{biascorrection}: Whether to apply bias correction for the
//'   variance estimate of individual restricted mean survival times.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' df <- rmdiff(data = rawdata, rep = "iterationNumber",
//'              stratum = "stratum", treat = "treatmentGroup",
//'              time = "timeUnderObservation", event = "event",
//'              milestone = 12)
//' head(df)
//'
//' @export
// [[Rcpp::export]]
DataFrame rmdiff(const DataFrame data,
                 const std::string rep = "rep",
                 const std::string stratum = "stratum",
                 const std::string treat = "treat",
                 const std::string time = "time",
                 const std::string event = "event",
                 const double milestone = NA_REAL,
                 const double rmstDiffH0 = 0,
                 const double confint = 0.95,
                 const bool biascorrection = 0) {
  int h, i, n = data.nrows();

  bool has_rep = hasVariable(data, rep);
  bool has_stratum = hasVariable(data, stratum);
  bool has_treat = hasVariable(data, treat);
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);

  if (!has_treat) {
    stop("data must contain the treat variable");
  }

  if (!has_time) {
    stop("data must contain the time variable");
  }

  if (!has_event) {
    stop("data must contain the event variable");
  }


  // create the numeric treat variable
  IntegerVector treatn(n);
  IntegerVector treatwi;
  NumericVector treatwn;
  CharacterVector treatwc;
  if (TYPEOF(data[treat]) == LGLSXP) {
    LogicalVector treatv = data[treat];
    treatn = 2 - treatv;
    treatwi = unique(treatn);
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwi = 2 - treatwi;
  } else if (TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    treatwi = unique(treatv);  // Get unique levels
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    std::sort(treatwi.begin(), treatwi.end());
    treatn = match(treatv, treatwi);  // Map codes to levels
  } else if (TYPEOF(data[treat]) == REALSXP) {
    NumericVector treatv = data[treat];
    treatwn = unique(treatv);
    if (treatwn.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    std::sort(treatwn.begin(), treatwn.end());
    treatn = match(treatv, treatwn);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    CharacterVector treatv = data[treat];
    treatwc = unique(treatv);
    if (treatwc.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    std::sort(treatwc.begin(), treatwc.end());
    treatn = match(treatv, treatwc);
  } else {
    stop("incorrect type for the treat variable in the input data");
  }


  NumericVector timen = data[time];
  NumericVector eventn = data[event];

  if (is_true(any(timen <= 0))) {
    stop("time must be positive for each subject");
  }

  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }


  // create the numeric rep variable
  IntegerVector repn(n);
  IntegerVector repwi;
  NumericVector repwn;
  CharacterVector repwc;
  if (!has_rep) {
    repn.fill(1);
  } else {
    if (TYPEOF(data[rep]) == INTSXP) {
      IntegerVector repv = data[rep];
      repwi = unique(repv);  // Get unique levels
      std::sort(repwi.begin(), repwi.end());
      repn = match(repv, repwi);  // Map codes to levels
    } else if (TYPEOF(data[rep]) == REALSXP) {
      NumericVector repv = data[rep];
      repwn = unique(repv);
      std::sort(repwn.begin(), repwn.end());
      repn = match(repv, repwn);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      CharacterVector repv = data[rep];
      repwc = unique(repv);
      std::sort(repwc.begin(), repwc.end());
      repn = match(repv, repwc);
    } else {
      stop("incorrect type for the rep variable in the input data");
    }
  }


  // create the numeric stratum variable
  IntegerVector stratumn(n);
  IntegerVector stratumwi;
  NumericVector stratumwn;
  CharacterVector stratumwc;
  if (!has_stratum) {
    stratumn.fill(1);
  } else {
    if (TYPEOF(data[stratum]) == INTSXP) {
      IntegerVector stratumv = data[stratum];
      stratumwi = unique(stratumv);  // Get unique levels
      std::sort(stratumwi.begin(), stratumwi.end());
      stratumn = match(stratumv, stratumwi);  // Map codes to levels
    } else if (TYPEOF(data[stratum]) == REALSXP) {
      NumericVector stratumv = data[stratum];
      stratumwn = unique(stratumv);
      std::sort(stratumwn.begin(), stratumwn.end());
      stratumn = match(stratumv, stratumwn);
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      CharacterVector stratumv = data[stratum];
      stratumwc = unique(stratumv);
      std::sort(stratumwc.begin(), stratumwc.end());
      stratumn = match(stratumv, stratumwc);
    } else {
      stop("incorrect type for the stratum variable in the input data");
    }
  }


  if (R_isnancpp(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (confint <= 0 || confint >= 1) {
    stop("confint must lie between 0 and 1");
  }


  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  repn = repn[order];
  stratumn = stratumn[order];
  treatn = treatn[order];
  timen = timen[order];
  eventn = eventn[order];

  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (i=1; i<n; i++) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }

  int nreps = idx.size();
  idx.push_back(n);

  IntegerVector rep0(nreps, NA_INTEGER);
  NumericVector rmst10(nreps), rmst20(nreps), rmstDiff0(nreps);
  NumericVector vrmst10(nreps), vrmst20(nreps), vrmstDiff0(nreps);
  NumericVector rmstDiffZ0(nreps), rmstDiffPValue0(nreps);
  NumericVector lower0(nreps), upper0(nreps);

  double z = R::qnorm((1.0 + confint)/2.0, 0, 1, 1, 0);

  bool noerr = 1;
  int index = 0;
  for (h=0; h<nreps; h++) {
    bool skip = 0;
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);

    IntegerVector stratum1 = stratumn[q1];
    IntegerVector treat1 = treatn[q1];
    NumericVector time1 = timen[q1];
    NumericVector event1 = eventn[q1];

    DataFrame dfin = DataFrame::create(
      _["stratum"] = stratum1,
      _["treat"] = treat1,
      _["time"] = time1,
      _["event"] = event1);

    DataFrame dfout = rmest(dfin, "stratum", "treat", "time", "event",
                            milestone, 0.95, biascorrection);

    IntegerVector stratum2 = dfout["stratum"];
    IntegerVector treat2 = dfout["treat"];
    IntegerVector treatsize = dfout["size"];
    NumericVector rmstime2 = dfout["rmst"];
    NumericVector stderr2 = dfout["stderr"];
    int n2 = stratum2.size();

    // identify the locations of the unique values of stratum
    IntegerVector idx2(1,0);
    for (i=1; i<n2; i++) {
      if (stratum2[i] != stratum2[i-1]) {
        idx2.push_back(i);
      }
    }

    int nstrata = idx2.size();
    idx2.push_back(n2);

    IntegerVector m(nstrata, 0); // number of subjects in each stratum
    for (i=0; i<nstrata; i++) {
      int j1 = idx2[i], j2 = idx2[i+1] - 1;
      if ((treat2[j1] != 1) || (treat2[j2] != 2)) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else if (TYPEOF(data[rep]) == INTSXP) {
          reperr = " " + rep + " = " + std::to_string(repwi[repn[idx[h]]-1]);
        } else if (TYPEOF(data[rep]) == REALSXP) {
          reperr = " " + rep + " = " + std::to_string(repwn[repn[idx[h]]-1]);
        } else {
          reperr = " " + rep + " = " + repwc[repn[idx[h]]-1];
        }

        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else if (TYPEOF(data[stratum]) == INTSXP) {
          stratumerr = " " + stratum + " = " +
            std::to_string(stratumwi[stratum2[j1]-1]);
        } else if (TYPEOF(data[stratum]) == REALSXP) {
          stratumerr = " " + stratum + " = " +
            std::to_string(stratumwn[stratum2[j1]-1]);
        } else {
          stratumerr = " " + stratum + " = " + stratumwc[stratum2[j1]-1];
        }

        int k = treat2[j1] != 1 ? 0 : 1;
        std::string treaterr;
        if ((TYPEOF(data[treat]) == LGLSXP) ||
            (TYPEOF(data[treat]) == INTSXP)) {
          treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
        } else if (TYPEOF(data[treat]) == REALSXP) {
          treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
        } else {
          treaterr = " " + treat + " = " + treatwc[k];
        }

        std::string str1 = "The data set does not contain";
        std::string errmsg = str1 + treaterr;
        if (!reperr.empty() || !stratumerr.empty()) {
          errmsg = errmsg + ":" + reperr + stratumerr;
        }

        if (noerr) {
          Rcout << errmsg << "\n";
          Rcout << "Additional warning messages are suppressed" << "\n";
          noerr = 0;
        }

        skip = 1;
        break;
      }

      m[i] += treatsize[j1] + treatsize[j2];
    }

    // skip the replication if there is a stratum without both treatments
    if (skip) continue;

    double M = sum(m);
    NumericVector p(nstrata);

    double rmst1 = 0.0, rmst2 = 0.0, vrmst1 = 0.0, vrmst2 = 0.0;
    for (i=0; i<nstrata; i++) {
      p[i] = m[i]/M; // fraction of subjects in the stratum
      IntegerVector q = Range(idx2[i], idx2[i+1]-1);
      NumericVector rmst = rmstime2[q];
      NumericVector stderrx = stderr2[q];
      NumericVector vrmst = stderrx*stderrx;

      rmst1 += p[i]*rmst[0];
      rmst2 += p[i]*rmst[1];
      vrmst1 += pow(p[i],2)*vrmst[0];
      vrmst2 += pow(p[i],2)*vrmst[1];
    }

    rep0[index] = repn[idx[h]];
    rmst10[index] = rmst1;
    rmst20[index] = rmst2;
    vrmst10[index] = vrmst1;
    vrmst20[index] = vrmst2;
    rmstDiff0[index] = rmst1 - rmst2;
    vrmstDiff0[index] = vrmst1 + vrmst2;
    double sermstDiff = sqrt(vrmstDiff0[index]);
    rmstDiffZ0[index] = (rmstDiff0[index] - rmstDiffH0)/sermstDiff;
    rmstDiffPValue0[index] = 1.0 - R::pnorm(rmstDiffZ0[index], 0, 1, 1, 0);
    lower0[index] = rmstDiff0[index] - z*sermstDiff;
    upper0[index] = rmstDiff0[index] + z*sermstDiff;

    index++;
  }

  // only keep nonmissing records
  LogicalVector sub = !is_na(rep0);
  if (is_false(any(sub))) {
    stop("no replication enables valid inference");
  }

  rep0 = rep0[sub];
  rmst10 = rmst10[sub];
  rmst20 = rmst20[sub];
  rmstDiff0 = rmstDiff0[sub];
  vrmst10 = vrmst10[sub];
  vrmst20 = vrmst20[sub];
  vrmstDiff0 = vrmstDiff0[sub];
  rmstDiffZ0 = rmstDiffZ0[sub];
  rmstDiffPValue0 = rmstDiffPValue0[sub];
  lower0 = lower0[sub];
  upper0 = upper0[sub];


  DataFrame result = DataFrame::create(
    _["milestone"] = milestone,
    _["rmstDiffH0"] = rmstDiffH0,
    _["rmst1"] = rmst10,
    _["rmst2"] = rmst20,
    _["rmstDiff"] = rmstDiff0,
    _["vrmst1"] = vrmst10,
    _["vrmst2"] = vrmst20,
    _["vrmstDiff"] = vrmstDiff0,
    _["rmstDiffZ"] = rmstDiffZ0,
    _["rmstDiffPValue"] = rmstDiffPValue0,
    _["lower"] = lower0,
    _["upper"] = upper0,
    _["confint"] = confint,
    _["biascorrection"] = biascorrection);

  if (has_rep) {
    if (TYPEOF(data[rep]) == INTSXP) {
      result.push_back(repwi[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      result.push_back(repwn[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      result.push_back(repwc[rep0-1], rep);
    }
  }

  return result;
}

