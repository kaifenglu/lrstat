#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


// define the integrand function for kmstat1
struct kmparams {
  double time;
  double phi;
  NumericVector accrualTime;
  NumericVector accrualIntensity;
  NumericVector piecewiseSurvivalTime;
  NumericVector lambda;
  NumericVector gamma;
  double accrualDuration;
};


void f_km(double *x, int n, void *ex) {
  kmparams *param = (kmparams *) ex;
  NumericVector u0(n);
  for (int i=0; i<n; i++) {
    u0[i] = x[i];
  }
  IntegerVector j = findInterval3(u0, param->piecewiseSurvivalTime) - 1;
  NumericVector lambda = param->lambda[j];
  NumericVector p = patrisk(u0, param->piecewiseSurvivalTime, param->lambda,
                            param->gamma);
  u0 = param->time - u0;
  NumericVector N = accrual(u0, param->accrualTime, param->accrualIntensity,
                            param->accrualDuration);
  u0 = lambda/((param->phi)*N*p);
  for (int i=0; i<n; i++) {
    x[i] = u0[i];
  }
}


//' @title Milestone Survival Probability by Stratum
//'
//' @description Obtains the milestone survival probability and associated
//' variance by treatment group and by stratum at a given calendar time.
//'
//' @param time The calendar time for data cut.
//' @param milestone The milestone time at which to calculate the
//'   survival probability.
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
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{nmilestone}: The total number of subjects reaching milestone.
//'
//' * \code{nmilestone1}: The number of subjects reaching milestone
//'   in the active treatment group.
//'
//' * \code{nmiletone2}: The number of subjects reaching milestone
//'   in the control group.
//'
//' * \code{surv1}: The milestone survival probability for the treatment
//'   group.
//'
//' * \code{surv2}: The milestone survival probability for the control group.
//'
//' * \code{survDiff}: The difference in milestone survival probabilities,
//'   i.e., \code{surv1 - surv2}.
//'
//' * \code{vsurv1}: The variance for \code{surv1}.
//'
//' * \code{vsurv2}: The variance for \code{surv2}.
//'
//' * \code{vsurvDiff}: The variance for \code{survDiff}.
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
//' kmstat1(time = 40,
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
DataFrame kmstat1(const double time = NA_REAL,
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

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
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

  // obtain the follow-up time for the first enrolled subject
  double maxFollowupTime;
  if (fixedFollowup) {
    maxFollowupTime = followupTime;
  } else {
    maxFollowupTime = accrualDuration + followupTime;
  }

  NumericVector ss(1, time);
  double a = accrual(ss, accrualTime, accrualIntensity, accrualDuration)[0];
  double phi = allocationRatioPlanned/(1 + allocationRatioPlanned);

  IntegerVector l1 = Range(0, nintervals-1);
  IntegerVector l(nintervals);
  NumericVector lam1(nintervals), lam2(nintervals);
  NumericVector gam1(nintervals), gam2(nintervals);
  NumericMatrix x(1,2), y(1,2);
  NumericVector nsubjects(nstrata);
  NumericMatrix nevents(nstrata, 2), ndropouts(nstrata, 2);
  NumericVector surv1(nstrata), surv2(nstrata), survDiff(nstrata);
  NumericVector vsurv1(nstrata), vsurv2(nstrata), vsurvDiff(nstrata);
  IntegerVector stratum(nstrata);
  NumericVector calTime(nstrata), mileTime(nstrata);
  DataFrame df;

  NumericVector milestone1(1, milestone);
  NumericVector t = piecewiseSurvivalTime;
  int jstar = findInterval3(milestone1, t)[0] - 1;
  double tol = 1e-6;

  NumericVector tt = NumericVector::create(time - milestone);
  double a2 = accrual(tt, accrualTime, accrualIntensity, accrualDuration)[0];
  NumericVector nmilestone1(nstrata), nmilestone2(nstrata);
  NumericVector nmilestone(nstrata);

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

    // number of events in the stratum at the specified calendar time
    x = nevent2(ss, allocationRatioPlanned, accrualTime,
                frac*accrualIntensity,
                piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
                accrualDuration, followupTime, maxFollowupTime);

    y = nevent2(ss, allocationRatioPlanned, accrualTime,
                frac*accrualIntensity,
                piecewiseSurvivalTime, gam1, gam2, lam1, lam2,
                accrualDuration, followupTime, maxFollowupTime);

    // obtain number of enrolled subjects and subjects having an event
    nsubjects[h] = frac*a;
    nevents(h, _) = x.row(0);
    ndropouts(h, _) = y.row(0);

    // obtain number of subjects censored due to reaching the max follow-up
    double ncom = frac*a2;
    double p1 = patrisk(milestone, piecewiseSurvivalTime, lam1, gam1)[0];
    double p2 = patrisk(milestone, piecewiseSurvivalTime, lam2, gam2)[0];
    nmilestone1[h] = phi*ncom*p1;
    nmilestone2[h] = (1-phi)*ncom*p2;
    nmilestone[h] = nmilestone1[h] + nmilestone2[h];

    double s1 = 0.0, s2 = 0.0;
    for (int j=0; j<=jstar; j++) {
      if (j<jstar) {
        s1 += lam1[j]*(t[j+1] - t[j]);
        s2 += lam2[j]*(t[j+1] - t[j]);
      } else {
        s1 += lam1[j]*(milestone - t[j]);
        s2 += lam2[j]*(milestone - t[j]);
      }
    }
    s1 = exp(-s1);
    s2 = exp(-s2);
    surv1[h] = s1;
    surv2[h] = s2;
    survDiff[h] = surv1[h] - surv2[h];

    kmparams param1 = {time, phi, accrualTime, frac*accrualIntensity,
                       t, lam1, gam1, accrualDuration};
    kmparams param2 = {time, 1-phi, accrualTime, frac*accrualIntensity,
                       t, lam2, gam2, accrualDuration};

    double q1 = quad(f_km, &param1, 0.0, milestone, tol)[0];
    double q2 = quad(f_km, &param2, 0.0, milestone, tol)[0];

    vsurv1[h] = s1*s1*q1;
    vsurv2[h] = s2*s2*q2;
    vsurvDiff[h] = vsurv1[h] + vsurv2[h];
  }

  // number of subjects having an event in each treatment group and overall
  NumericVector nevents1 = nevents(_, 0);
  NumericVector nevents2 = nevents(_, 1);
  NumericVector neventst = nevents1 + nevents2;

  NumericVector ndropouts1 = ndropouts(_, 0);
  NumericVector ndropouts2 = ndropouts(_, 1);
  NumericVector ndropoutst = ndropouts1 + ndropouts2;

  // output the requested information
  df = DataFrame::create(_["stratum"] = stratum,
                         _["time"] = calTime,
                         _["subjects"] = nsubjects,
                         _["nevents"] = neventst,
                         _["nevents1"] = nevents1,
                         _["nevents2"] = nevents2,
                         _["ndropouts"] = ndropoutst,
                         _["ndropouts1"] = ndropouts1,
                         _["ndropouts2"] = ndropouts2,
                         _["milestone"] = mileTime,
                         _["nmilestone"] = nmilestone,
                         _["nmilestone1"] = nmilestone1,
                         _["nmilestone2"] = nmilestone2,
                         _["surv1"] = surv1,
                         _["surv2"] = surv2,
                         _["survDiff"] = survDiff,
                         _["vsurv1"] = vsurv1,
                         _["vsurv2"] = vsurv2,
                         _["vsurvDiff"] = vsurvDiff);

  return df;
}


//' @title Stratified Difference in Milestone Survival Probabilities
//' @description Obtains the stratified milestone survival probabilities
//' and difference in milestone survival probabilities at given
//' calendar times.
//'
//' @param time A vector of calendar times for data cut.
//' @param milestone The milestone time at which to calculate the
//'   survival probability.
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
//' * \code{time}: The calendar time since trial start.
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
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{nmilestone}: The total number of subjects reaching milestone.
//'
//' * \code{nmilestone1}: The number of subjects reaching milestone
//'   in the active treatment group.
//'
//' * \code{nmiletone2}: The number of subjects reaching milestone
//'   in the control group.
//'
//' * \code{surv1}: The milestone survival probability for the treatment
//'   group.
//'
//' * \code{surv2}: The milestone survival probability for the control group.
//'
//' * \code{survDiff}: The difference in milestone survival probabilities,
//'   i.e., \code{surv1 - surv2}.
//'
//' * \code{vsurv1}: The variance for \code{surv1}.
//'
//' * \code{vsurv2}: The variance for \code{surv2}.
//'
//' * \code{vsurvDiff}: The variance for \code{survDiff}.
//'
//' * \code{information}: The information for \code{survDiff}, equal to
//'   \code{1/vsurvDiff}.
//'
//' * \code{survDiffZ}: The Z-statistic value, i.e.,
//'   \code{survDiff/sqrt(vsurvDiff)}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' kmstat(time = c(22, 40),
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
DataFrame kmstat(const NumericVector& time = NA_REAL,
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

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi), gamma1x(nsi), gamma2x(nsi);

  if (is_true(any(is_na(time)))) {
    stop("time must be provided");
  }

  if (is_true(any(time <= 0))) {
    stop("time must be positive");
  }

  if (std::isnan(milestone)) {
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

  if (std::isnan(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (std::isnan(followupTime)) {
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


  int k = static_cast<int>(time.size());
  NumericVector calTime(k), mileTime(k), subjects(k),
  nevents(k), nevents1(k), nevents2(k),
  ndropouts(k), ndropouts1(k), ndropouts2(k),
  nmilestone(k), nmilestone1(k), nmilestone2(k),
  surv1(k), surv2(k), vsurv1(k), vsurv2(k),
  survDiff(k), vsurvDiff(k), information(k), survDiffZ(k);
  DataFrame df;

  for (int j=0; j<k; j++) {
    df = kmstat1(time[j], milestone, allocationRatioPlanned,
                 accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1x, lambda2x, gamma1x, gamma2x,
                 accrualDuration, followupTime, fixedFollowup);

    calTime[j] = max(NumericVector(df[1]));
    subjects[j] = sum(NumericVector(df[2]));
    nevents[j] = sum(NumericVector(df[3]));
    nevents1[j] = sum(NumericVector(df[4]));
    nevents2[j] = sum(NumericVector(df[5]));
    ndropouts[j] = sum(NumericVector(df[6]));
    ndropouts1[j] = sum(NumericVector(df[7]));
    ndropouts2[j] = sum(NumericVector(df[8]));
    mileTime[j] = max(NumericVector(df[9]));
    nmilestone[j] = sum(NumericVector(df[10]));
    nmilestone1[j] = sum(NumericVector(df[11]));
    nmilestone2[j] = sum(NumericVector(df[12]));
    surv1[j] = sum(stratumFraction*NumericVector(df[13]));
    surv2[j] = sum(stratumFraction*NumericVector(df[14]));
    survDiff[j] = sum(stratumFraction*NumericVector(df[15]));
    vsurv1[j] = sum(stratumFraction*stratumFraction*NumericVector(df[16]));
    vsurv2[j] = sum(stratumFraction*stratumFraction*NumericVector(df[17]));
    vsurvDiff[j] = sum(stratumFraction*stratumFraction*NumericVector(df[18]));
    information[j] = 1.0/vsurvDiff[j];
    survDiffZ[j] = survDiff[j]/sqrt(vsurvDiff[j]);
  }

  df = DataFrame::create(_["time"] = calTime,
                         _["subjects"] = subjects,
                         _["nevents"] = nevents,
                         _["nevents1"] = nevents1,
                         _["nevents2"] = nevents2,
                         _["ndropouts"] = ndropouts,
                         _["ndropouts1"] = ndropouts1,
                         _["ndropouts2"] = ndropouts2,
                         _["milestone"] = mileTime,
                         _["nmilestone"] = nmilestone,
                         _["nmilestone1"] = nmilestone1,
                         _["nmilestone2"] = nmilestone2,
                         _["surv1"] = surv1,
                         _["surv2"] = surv2,
                         _["survDiff"] = survDiff,
                         _["vsurv1"] = vsurv1,
                         _["vsurv2"] = vsurv2,
                         _["vsurvDiff"] = vsurvDiff,
                         _["information"] = information,
                         _["survDiffZ"] = survDiffZ);

  return df;
}


//' @title Power for Difference in Milestone Survival Probabilities
//' @description Estimates the power for testing the difference in
//' milestone survival probabilities in a two-sample survival design.
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
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survDiffH0 The difference in milestone survival probabilities
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
//' @return An S3 class \code{kmpower} object with 4 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
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
//'     - \code{survDiffH0}: The difference in milestone survival
//'       probabilities under the null hypothesis.
//'
//'     - \code{surv1}: The milestone survival probability for the
//'       treatment group.
//'
//'     - \code{surv2}: The milestone survival probability for the
//'       control group.
//'
//'     - \code{survDiff}: The difference in milestone survival
//'       probabilities, equal to \code{surv1 - surv2}.
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
//'     - \code{numberOfMilestone}: The number of subjects reaching
//'       milestone.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacySurvDiff}: The efficacy boundaries on the survival
//'       difference scale.
//'
//'     - \code{futilitySurvDiff}: The futility boundaries on the survival
//'       difference scale.
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
//'     - \code{numberOfMilestone1}: The number of subjects reaching
//'       milestone by stage for the active treatment group.
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
//'     - \code{numberOfMilestone2}: The number of subjects reaching
//'       milestone by stage for the control group.
//'
//'     - \code{expectedNumberOfEvents1}: The expected number of events for
//'       the treatment group.
//'
//'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfMilestone1}: The expected number of subjects
//'       reaching milestone for the active treatment group.
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
//'     - \code{expectedNumberOfMilestone2}: The expected number of subjects
//'       reaching milestone for the control group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survival, and 5% dropout by
//' # the end of 1 year.
//'
//' kmpower(kMax = 2, informationRates = c(0.8, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         milestone = 18,
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
List kmpower(const int kMax = 1,
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
             const double survDiffH0 = 0,
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
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
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

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop("parameterBetaSpending must be positive for sfKD");
  }

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (survDiffH0 <= -1 || survDiffH0 >= 1) {
    stop("survDiffH0 must lie between -1 and 1");
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

  if (std::isnan(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (std::isnan(followupTime)) {
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

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  if (milestone >= accrualDuration + followupTime) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      (milestone >= studyDuration)) {
    stop("milestone cannot exceed studyDuration for fixed follow-up");
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

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
  if (!fixedFollowup || std::isnan(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }

  // obtain the timing of interim analysis
  DataFrame km;
  NumericVector time(kMax), sdu(kMax), sdl(kMax);
  NumericVector u0(1, studyDuration1);
  km = kmstat(u0, milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup);

  double maxInformation = sum(NumericVector(km[18]));
  double surv1 = sum(NumericVector(km[12]));
  double surv2 = sum(NumericVector(km[13]));
  double survDiff = sum(NumericVector(km[14]));
  NumericVector theta(kMax, survDiff - survDiffH0);
  NumericVector I = maxInformation*informationRates1;

  double information1;
  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double aval)->double {
              NumericVector u0(1, aval);
              DataFrame km = kmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup);
              return sum(NumericVector(km[18])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    information1 = I[i];
    time[i] = brent(f, milestone + 1.0e-6, studyDuration1, 1.0e-6);
  };
  time[kMax-1] = studyDuration1;

  km = kmstat(time, milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup);

  double phi = allocationRatioPlanned/(allocationRatioPlanned+1);
  NumericVector nsubjects = NumericVector(km[1]);
  NumericVector nsubjects1 = phi*nsubjects;
  NumericVector nsubjects2 = (1-phi)*nsubjects;
  NumericVector nevents = NumericVector(km[2]);
  NumericVector nevents1 = NumericVector(km[3]);
  NumericVector nevents2 = NumericVector(km[4]);
  NumericVector ndropouts = NumericVector(km[5]);
  NumericVector ndropouts1 = NumericVector(km[6]);
  NumericVector ndropouts2 = NumericVector(km[7]);
  NumericVector nmilestone = NumericVector(km[9]);
  NumericVector nmilestone1 = NumericVector(km[10]);
  NumericVector nmilestone2 = NumericVector(km[11]);

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
  double expectedNumberOfSubjects = sum(ptotal*nsubjects);
  double expectedNumberOfEvents1 = sum(ptotal*nevents1);
  double expectedNumberOfDropouts1 = sum(ptotal*ndropouts1);
  double expectedNumberOfSubjects1 = sum(ptotal*nsubjects1);
  double expectedNumberOfMilestone1 = sum(ptotal*nmilestone1);
  double expectedNumberOfEvents2 = sum(ptotal*nevents2);
  double expectedNumberOfDropouts2 = sum(ptotal*ndropouts2);
  double expectedNumberOfSubjects2 = sum(ptotal*nsubjects2);
  double expectedNumberOfMilestone2 = sum(ptotal*nmilestone2);
  double expectedStudyDuration = sum(ptotal*time);
  double expectedInformation = sum(ptotal*I);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  sdu = survDiffH0 + criticalValues1/sqrt(I);
  sdl = survDiffH0 + futilityBounds1/sqrt(I);

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      sdu[i] = NA_REAL;
      efficacyStopping1[i] = 0;
    }

    if (futilityBounds1[i] == -6) {
      sdl[i] = NA_REAL;
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
    _["numberOfMilestone"] = nmilestone,
    _["analysisTime"] = time,
    _["efficacySurvDiff"] = sdu,
    _["futilitySurvDiff"] = sdl,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["information"] = I,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["kMax"] = kMax,
    _["milestone"] = milestone,
    _["survDiffH0"] = survDiffH0,
    _["surv1"] = surv1,
    _["surv2"] = surv2,
    _["survDiff"] = survDiff);

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

  List byTreatmentCounts = List::create(
    _["numberOfEvents1"] = nevents1,
    _["numberOfDropouts1"] = ndropouts1,
    _["numberOfSubjects1"] = nsubjects1,
    _["numberOfMilestone1"] = nmilestone1,
    _["numberOfEvents2"] = nevents2,
    _["numberOfDropouts2"] = ndropouts2,
    _["numberOfSubjects2"] = nsubjects2,
    _["numberOfMilestone2"] = nmilestone2,
    _["expectedNumberOfEvents1"] = expectedNumberOfEvents1,
    _["expectedNumberOfDropouts1"] = expectedNumberOfDropouts1,
    _["expectedNumberOfSubjects1"] = expectedNumberOfSubjects1,
    _["expectedNumberOfMilestone1"] = expectedNumberOfMilestone1,
    _["expectedNumberOfEvents2"] = expectedNumberOfEvents2,
    _["expectedNumberOfDropouts2"] = expectedNumberOfDropouts2,
    _["expectedNumberOfSubjects2"] = expectedNumberOfSubjects2,
    _["expectedNumberOfMilestone2"] = expectedNumberOfMilestone2);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings,
    _["byTreatmentCounts"] = byTreatmentCounts);

  result.attr("class") = "kmpower";

  return result;
}


//' @title Sample Size for Difference in Milestone Survival Probabilities
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
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survDiffH0 The difference in milestone survival probabilities
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
//' * \code{resultsUnderH1}: An S3 class \code{kmpower} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{kmpower} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{kmpower}}
//'
//' @examples
//' # Example 1: Obtains follow-up time given power, accrual intensity,
//' # and accrual duration for variable follow-up. Of note, the power
//' # reaches the maximum when the follow-up time equals milestone.
//'
//' kmsamplesize(beta = 0.25, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
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
//' kmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
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
//' kmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
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
List kmsamplesize(const double beta = 0.2,
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
                  const double survDiffH0 = 0,
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
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
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

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
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

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (survDiffH0 <= -1 || survDiffH0 >= 1) {
    stop("survDiffH0 must lie between -1 and 1");
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

  if (!std::isnan(accrualDuration)) {
    if (accrualDuration <= 0) {
      stop("accrualDuration must be positive");
    }
  }

  if (!std::isnan(followupTime)) {
    if (fixedFollowup && followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }

    if (!fixedFollowup && followupTime < 0) {
      stop("followupTime must be non-negative for variable follow-up");
    }
  }

  if (fixedFollowup && std::isnan(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
  }

  if (!std::isnan(accrualDuration) && !std::isnan(followupTime) &&
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
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

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
  if (std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualDuration";
  } else if (!std::isnan(accrualDuration) && std::isnan(followupTime)) {
    unknown = "followupTime";
  } else if (!std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualIntensity";
  } else {
    stop("accrualDuration and followupTime cannot be both missing");
  }


  NumericVector survs1(nstrata), survs2(nstrata);
  IntegerVector l1 = Range(0, nintervals-1);
  NumericVector zerogam(nintervals);
  NumericVector u0(1, milestone);
  for (int h=0; h<nstrata; h++) {
    l = h*nintervals + l1;
    NumericVector lam1 = lambda1x[l];
    NumericVector lam2 = lambda2x[l];
    survs1[h] = patrisk(u0, piecewiseSurvivalTime, lam1, zerogam)[0];
    survs2[h] = patrisk(u0, piecewiseSurvivalTime, lam2, zerogam)[0];
  }

  double surv1 = sum(stratumFraction*survs1);
  double surv2 = sum(stratumFraction*survs2);
  double theta1 = surv1 - surv2 - survDiffH0;
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
              DataFrame km = kmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                dur1, dur2, fixedFollowup);
              return sum(NumericVector(km[18])) - maxInformation;
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
                  DataFrame km = kmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    aval, followupTime, fixedFollowup);
                  return sum(NumericVector(km[18])) - maxInformation;
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
                  DataFrame km = kmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup);
                  return sum(NumericVector(km[18])) - maxInformation;
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
    double n = std::ceil(n0 - 1.0e-12);

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
                    DataFrame km = kmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval, fixedFollowup);
                    return sum(NumericVector(km[18])) - maxInformation;
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
                    DataFrame km = kmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup);
                    return sum(NumericVector(km[18])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        double aval = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }


  resultH1 = kmpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, survDiffH0,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  overallResults = DataFrame(resultH1["overallResults"]);
  maxInformation = overallResults["information"];


  // obtain results under H0 by matching the maximum information
  // first find the hazard rate for the active treatment group that yields
  // the specified difference in milestone survival probabilities under H0
  auto fsurv = [milestone, piecewiseSurvivalTime, stratumFraction,
                nintervals, nstrata, l1, lambda2x, zerogam, surv2,
                survDiffH0](double aval)-> double {
                  NumericVector survs1(nstrata);
                  NumericVector u0(1, milestone);
                  for (int h=0; h<nstrata; h++) {
                    IntegerVector l = h*nintervals + l1;
                    NumericVector lam2 = lambda2x[l];
                    survs1[h] = patrisk(u0, piecewiseSurvivalTime,
                                        aval*lam2, zerogam)[0];
                  }
                  double surv1 = sum(stratumFraction*survs1);
                  return surv1 - surv2 - survDiffH0;
                };

  double lower = 0.001, upper = 1.999;
  while (fsurv(upper) > 0) {
    lower = upper;
    upper = 2.0*upper;
  }
  double aval = brent(fsurv, lower, upper, 1.0e-6);
  NumericVector lambda1H0 = aval*lambda2;

  if (!fixedFollowup) {
    auto h = [milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambda1H0, lambda2, gamma1, gamma2,
              accrualDuration, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                DataFrame km = kmstat(
                  u0, milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1H0, lambda2, gamma1, gamma2,
                  accrualDuration, aval, fixedFollowup);
                return sum(NumericVector(km[18])) - maxInformation;
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
                  DataFrame km = kmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1H0, lambda2, gamma1, gamma2,
                    aval, followupTime, fixedFollowup);
                  return sum(NumericVector(km[18])) - maxInformation;
                };

      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (g(upper) < 0) {
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
                  DataFrame km = kmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1H0, lambda2, gamma1, gamma2,
                    aval, 0, fixedFollowup);
                  return sum(NumericVector(km[18])) - maxInformation;
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
                DataFrame km = kmstat(
                  u0, milestone, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1H0, lambda2, gamma1, gamma2,
                  aval, followupTime, fixedFollowup);
                return sum(NumericVector(km[18])) - maxInformation;
              };

    if (h(accrualDuration) < 0) { // increase accrual duration
      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (h(upper) < 0) {
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
                  DataFrame km = kmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1H0, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup);
                  return sum(NumericVector(km[18])) - maxInformation;
                };

      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      double aval = brent(g, lower, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    }
  }


  // use the same stopping boundaries as under H1
  resultH0 = kmpower(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, survDiffH0,
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


//' @title Power for One-Sample Milestone Survival Probability
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
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survH0 The milestone survival probability under the null
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
//' @return An S3 class \code{kmpower1s} object with 3 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{drift}: The drift parameter, equal to
//'       \code{(surv - survH0)*sqrt(information)}.
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
//'     - \code{milestone}: The milestone time to calculate the survival
//'       probability.
//'
//'     - \code{survH0}: The milestone survival probability under the null
//'       hypothesis.
//'
//'     - \code{surv}: The milestone survival probability under the
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
//'     - \code{efficacySurv}: The efficacy boundaries on the milestone
//'       survival probability scale.
//'
//'     - \code{futilitySurv}: The futility boundaries on the milestone
//'       survival probability scale.
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
//' @seealso \code{\link{kmstat}}
//'
//' @examples
//'
//' kmpower1s(kMax = 2, informationRates = c(0.8, 1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           milestone = 18, survH0 = 0.30,
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
List kmpower1s(const int kMax = 1,
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
               const double survH0 = NA_REAL,
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
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
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

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (std::isnan(survH0)) {
    stop("survH0 must be provided");
  }

  if (survH0 <= 0 || survH0 >= 1) {
    stop("survH0 must lie between 0 and 1");
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

  if (std::isnan(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (std::isnan(followupTime)) {
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

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  if (milestone >= accrualDuration + followupTime) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      (milestone >= studyDuration)) {
    stop("milestone cannot exceed studyDuration for fixed follow-up");
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

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
  if (!fixedFollowup || std::isnan(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }


  // obtain the timing of interim analysis using the twin treatment group
  DataFrame km;
  NumericVector time(kMax), survu(kMax), survl(kMax);
  NumericVector u0(1, studyDuration1);
  km = kmstat(u0, milestone, 1,
              accrualTime, 2.0*accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda, lambda, gamma, gamma,
              accrualDuration, followupTime, fixedFollowup);

  double maxInformation = 2.0*sum(NumericVector(km[18]));
  double surv = sum(NumericVector(km[12]));
  NumericVector theta = rep(surv - survH0, kMax);
  NumericVector I = maxInformation*informationRates1;

  double information1;
  auto f = [milestone, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda, gamma,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double aval)->double {
              NumericVector u0(1, aval);
              DataFrame km = kmstat(
                u0, milestone, 1,
                accrualTime, 2.0*accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda, lambda, gamma, gamma,
                accrualDuration, followupTime, fixedFollowup);
              return 2.0*sum(NumericVector(km[18])) - information1;
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

  survu = survH0 + criticalValues1/sqrt(I);
  survl = survH0 + futilityBounds1/sqrt(I);

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      survu[i] = NA_REAL;
      efficacyStopping1[i] = 0;
    }

    if (futilityBounds1[i] == -6) {
      survl[i] = NA_REAL;
      futilityStopping1[i] = 0;
    }
  }

  double drift = (surv - survH0)*sqrt(maxInformation);
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
    _["efficacySurv"] = survu,
    _["futilitySurv"] = survl,
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
    _["survH0"] = survH0,
    _["surv"] = surv);

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

  result.attr("class") = "kmpower1s";

  return result;
}


//' @title Sample Size for One-Sample Milestone Survival Probability
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
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survH0 The milestone survival probability under the null
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
//' * \code{resultsUnderH1}: An S3 class \code{kmpower1s} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{kmpower1s} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{kmpower1s}}
//'
//' @examples
//' # Example 1: Obtains follow-up duration given power, accrual intensity,
//' # and accrual duration for variable follow-up
//'
//' kmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, survH0 = 0.30,
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
//' kmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, survH0 = 0.30,
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
//' kmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, survH0 = 0.30,
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
List kmsamplesize1s(const double beta = 0.2,
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
                    const double survH0 = NA_REAL,
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
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
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

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
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

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (std::isnan(survH0)) {
    stop("survH0 must be provided");
  }

  if (survH0 <= 0 || survH0 >= 1) {
    stop("survH0 must lie between 0 and 1");
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

  if (!std::isnan(accrualDuration)) {
    if (accrualDuration <= 0) {
      stop("accrualDuration must be positive");
    }
  }

  if (!std::isnan(followupTime)) {
    if (fixedFollowup && followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }

    if (!fixedFollowup && followupTime < 0) {
      stop("followupTime must be non-negative for variable follow-up");
    }
  }

  if (fixedFollowup && std::isnan(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
  }

  if (!std::isnan(accrualDuration) && !std::isnan(followupTime) &&
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
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

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
  if (std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualDuration";
  } else if (!std::isnan(accrualDuration) && std::isnan(followupTime)) {
    unknown = "followupTime";
  } else if (!std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualIntensity";
  } else {
    stop("accrualDuration and followupTime cannot be both missing");
  }


  NumericVector survs(nstrata);
  IntegerVector l1 = Range(0, nintervals-1);
  NumericVector zerogam(nintervals);
  NumericVector u0(1, milestone);
  for (int h=0; h<nstrata; h++) {
    l = h*nintervals + l1;
    NumericVector lam = lambdax[l];
    survs[h] = patrisk(u0, piecewiseSurvivalTime, lam, zerogam)[0];
  }

  double surv = sum(stratumFraction*survs);
  double theta1 = surv - survH0;
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
              DataFrame km = kmstat(
                u0, milestone, 1,
                accrualTime, 2.0*accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda, lambda, gamma, gamma,
                dur1, dur2, fixedFollowup);
              return 2.0*sum(NumericVector(km[18])) - maxInformation;
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
                  DataFrame km = kmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda, lambda, gamma, gamma,
                    aval, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(km[18])) - maxInformation;
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
                  DataFrame km = kmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda, lambda, gamma, gamma,
                    accrualDuration, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(km[18])) - maxInformation;
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
    double n = std::ceil(n0 - 1.0e-12);

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
                    DataFrame km = kmstat(
                      u0, milestone, 1,
                      accrualTime, 2.0*accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda, lambda, gamma, gamma,
                      accrualDuration, aval, fixedFollowup);
                    return 2.0*sum(NumericVector(km[18])) - maxInformation;
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
                    DataFrame km = kmstat(
                      u0, milestone, 1,
                      accrualTime, 2.0*accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda, lambda, gamma, gamma,
                      accrualDuration, followupTime, fixedFollowup);
                    return 2.0*sum(NumericVector(km[18])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        double aval = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }


  resultH1 = kmpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, survH0,
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
                nintervals, nstrata, l1, lambdax, zerogam,
                survH0](double aval)-> double {
                  NumericVector survs(nstrata);
                  NumericVector u0(1, milestone);
                  for (int h=0; h<nstrata; h++) {
                    IntegerVector l = h*nintervals + l1;
                    NumericVector lam = lambdax[l];
                    survs[h] = patrisk(u0, piecewiseSurvivalTime,
                                       aval*lam, zerogam)[0];
                  }
                  double surv = sum(stratumFraction*survs);
                  return surv - survH0;
                };

  double lower = 0.001, upper = 1.999;
  while (fsurv(upper) > 0) {
    lower = upper;
    upper = 2.0*upper;
  }
  double aval = brent(fsurv, lower, upper, 1.0e-6);
  NumericVector lambdaH0 = aval*lambda;

  if (!fixedFollowup) {
    auto h = [milestone, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambdaH0, gamma,
              accrualDuration, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                DataFrame km = kmstat(
                  u0, milestone, 1,
                  accrualTime, 2.0*accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambdaH0, lambdaH0, gamma, gamma,
                  accrualDuration, aval, fixedFollowup);
                return 2.0*sum(NumericVector(km[18])) - maxInformation;
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
                  DataFrame km = kmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambdaH0, lambdaH0, gamma, gamma,
                    aval, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(km[18])) - maxInformation;
                };

      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (g(upper) < 0) {
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
                lambdaH0, gamma, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  DataFrame km = kmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambdaH0, lambdaH0, gamma, gamma,
                    aval, 0, fixedFollowup);
                  return 2.0*sum(NumericVector(km[18])) - maxInformation;
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
                DataFrame km = kmstat(
                  u0, milestone, 1,
                  accrualTime, 2.0*accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambdaH0, lambdaH0, gamma, gamma,
                  aval, followupTime, fixedFollowup);
                return 2.0*sum(NumericVector(km[18])) - maxInformation;
              };

    if (h(accrualDuration) < 0) { // increase accrual duration
      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (h(upper) < 0) {
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
                  DataFrame km = kmstat(
                    u0, milestone, 1,
                    accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambdaH0, lambdaH0, gamma, gamma,
                    accrualDuration, followupTime, fixedFollowup);
                  return 2.0*sum(NumericVector(km[18])) - maxInformation;
                };

      double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
      double aval = brent(g, lower, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    }
  }


  // use the same stopping boundaries as under H1
  resultH0 = kmpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending,
    milestone, survH0,
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



//' @title Power for Equivalence in Milestone Survival Probability Difference
//' @description Obtains the power for equivalence in milestone survival
//' probability difference.
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
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survDiffLower The lower equivalence limit of milestone survival
//'   probability difference.
//' @param survDiffUpper The upper equivalence limit of milestone survival
//'   probability difference.
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
//' @return An S3 class \code{kmpowerequiv} object with 4 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numbeOfEvents}: The total number of events.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
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
//'     - \code{survDiffLower}: The lower equivalence limit of milestone
//'       survival probability difference.
//'
//'     - \code{survDiffUpper}: The upper equivalence limit of milestone
//'       survival probability difference.
//'
//'     - \code{surv1}: The milestone survival probability for the
//'       treatment group.
//'
//'     - \code{surv2}: The milestone survival probability for the
//'       control group.
//'
//'     - \code{survDiff}: The milestone survival probability difference.
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
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{numberOfMilestone}: The number of subjects reaching
//'       milestone.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacySurvDiffLower}: The efficacy boundaries on the
//'       milestone survival probability difference scale for the one-sided
//'       null hypothesis at the lower equivalence limit.
//'
//'     - \code{efficacySurvDiffUpper}: The efficacy boundaries on the
//'       milestone survival probability difference scale for the one-sided
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
//'     - \code{numberOfMilestone1}: The number of subjects reaching
//'       milestone by stage for the active treatment group.
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
//'     - \code{numberOfMilestone2}: The number of subjects reaching
//'       milestone by stage for the control group.
//'
//'     - \code{expectedNumberOfEvents1}: The expected number of events for
//'       the treatment group.
//'
//'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfMilestone1}: The expected number of subjects
//'       reaching milestone for the active treatment group.
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
//'     - \code{expectedNumberOfMilestone2}: The expected number of subjects
//'       reaching milestone for the control group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{kmstat}}
//'
//' @examples
//'
//' kmpowerequiv(kMax = 2, informationRates = c(0.5, 1),
//'              alpha = 0.05, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              survDiffLower = -0.13, survDiffUpper = 0.13,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
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
List kmpowerequiv(const int kMax = 1,
                  const NumericVector& informationRates = NA_REAL,
                  const NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.05,
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const double milestone = NA_REAL,
                  const double survDiffLower = NA_REAL,
                  const double survDiffUpper = NA_REAL,
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
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
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

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (std::isnan(survDiffLower)) {
    stop("survDiffLower must be provided");
  }

  if (std::isnan(survDiffUpper)) {
    stop("survDiffUpper must be provided");
  }

  if (survDiffLower <= -1) {
    stop("survDiffLower must be greater than -1");
  }

  if (survDiffUpper >= 1) {
    stop("survDiffUpper must be less than 1");
  }

  if (survDiffLower >= survDiffUpper) {
    stop("survDiffLower must be less than survDiffUpper");
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

  if (std::isnan(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
  }

  if (std::isnan(followupTime)) {
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

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  if (milestone >= accrualDuration + followupTime) {
    stop("milestone must be less than accrualDuration + followupTime");
  }

  if (fixedFollowup && (milestone > followupTime)) {
    stop("milestone cannot exceed followupTime for fixed follow-up");
  }

  if (fixedFollowup && !std::isnan(studyDuration) &&
      (milestone >= studyDuration)) {
    stop("milestone cannot exceed studyDuration for fixed follow-up");
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

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
  if (!fixedFollowup || std::isnan(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }

  // obtain the timing of interim analysis
  NumericVector time(kMax);
  NumericVector u0(1, studyDuration1);
  DataFrame km = kmstat(u0, milestone, allocationRatioPlanned,
                        accrualTime, accrualIntensity,
                        piecewiseSurvivalTime, stratumFraction,
                        lambda1, lambda2, gamma1, gamma2,
                        accrualDuration, followupTime, fixedFollowup);

  double surv1 = sum(NumericVector(km[12]));
  double surv2 = sum(NumericVector(km[13]));
  double survDiff = sum(NumericVector(km[14]));
  NumericVector theta(kMax, survDiff);
  double maxInformation = sum(NumericVector(km[18]));
  NumericVector I = maxInformation*informationRates1;
  double information1;

  auto f = [milestone, allocationRatioPlanned,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double aval)->double {
              NumericVector u0(1, aval);
              DataFrame km = kmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup);
              return sum(NumericVector(km[18])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    information1 = I[i];
    time[i] = brent(f, milestone + 1.0e-6, studyDuration1, 1.0e-6);
  };
  time[kMax-1] = studyDuration1;

  km = kmstat(time, milestone, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup);

  double phi = allocationRatioPlanned/(allocationRatioPlanned+1);
  NumericVector nsubjects = NumericVector(km[1]);
  NumericVector nsubjects1 = phi*nsubjects;
  NumericVector nsubjects2 = (1-phi)*nsubjects;
  NumericVector nevents = NumericVector(km[2]);
  NumericVector nevents1 = NumericVector(km[3]);
  NumericVector nevents2 = NumericVector(km[4]);
  NumericVector ndropouts = NumericVector(km[5]);
  NumericVector ndropouts1 = NumericVector(km[6]);
  NumericVector ndropouts2 = NumericVector(km[7]);
  NumericVector nmilestone = NumericVector(km[9]);
  NumericVector nmilestone1 = NumericVector(km[10]);
  NumericVector nmilestone2 = NumericVector(km[11]);

  // calculate cumulative rejection probability under H1
  NumericVector theta10 = rep(survDiffLower, kMax);
  NumericVector theta20 = rep(survDiffUpper, kMax);
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

  NumericVector efficacySurvDiffLower = theta10 + b/sqrt(I);
  NumericVector efficacySurvDiffUpper = theta20 - b/sqrt(I);

  // calculate cumulative rejection under H10
  List probs2H10 = exitprobcpp(ui, pmin(u, ui), theta10, I);

  NumericVector cplH10 = cumAlphaSpent;
  NumericVector cpuH10 = cumsum(NumericVector(probs2H10[1]));

  NumericVector cpH10(kMax);
  if (k.size() == 0) {
    cpH10 = cplH10 + cpuH10 - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    List a = exitprobcpp(l[idx], u[idx], theta10[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH10[i] = cplH10[i] + cpuH10[i] - ca[i];
      } else {
        cpH10[i] = cplH10[i] + cpuH10[i] - 1;
      }
    }
  }

  // calculate cumulative rejection under H20
  List probs1H20 = exitprobcpp(pmax(l, li), li, theta20, I);

  NumericVector cplH20 = cumsum(NumericVector(probs1H20[0]));
  NumericVector cpuH20 = cumAlphaSpent;

  NumericVector cpH20(kMax);
  if (k.size() == 0) {
    cpH20 = cplH20 + cpuH20 - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    NumericVector uH20 = -b;
    List a = exitprobcpp(l[idx], u[idx], theta20[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH20[i] = cplH20[i] + cpuH20[i] - ca[i];
      } else {
        cpH20[i] = cplH20[i] + cpuH20[i] - 1;
      }
    }
  }

  double overallReject = cp[kMax-1];
  double expectedNumberOfEvents = sum(q*nevents);
  double expectedNumberOfSubjects = sum(q*nsubjects);
  double expectedNumberOfEvents1 = sum(q*nevents1);
  double expectedNumberOfDropouts1 = sum(q*ndropouts1);
  double expectedNumberOfSubjects1 = sum(q*nsubjects1);
  double expectedNumberOfMilestone1 = sum(q*nmilestone1);
  double expectedNumberOfEvents2 = sum(q*nevents2);
  double expectedNumberOfDropouts2 = sum(q*ndropouts2);
  double expectedNumberOfSubjects2 = sum(q*nsubjects2);
  double expectedNumberOfMilestone2 = sum(q*nmilestone2);
  double expectedStudyDuration = sum(q*time);
  double expectedInformation = sum(q*I);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = alpha,
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["kMax"] = kMax,
    _["milestone"] = milestone,
    _["survDiffLower"] = survDiffLower,
    _["survDiffUpper"] = survDiffUpper,
    _["surv1"] = surv1,
    _["surv2"] = surv2,
    _["survDiff"] = survDiff,
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
    _["numberOfEvents"] = nevents,
    _["numberOfDropouts"] = ndropouts,
    _["numberOfSubjects"] = nsubjects,
    _["numberOfMilestone"] = nmilestone,
    _["analysisTime"] = time,
    _["efficacySurvDiffLower"] = efficacySurvDiffLower,
    _["efficacySurvDiffUpper"] = efficacySurvDiffUpper,
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

  List byTreatmentCounts = List::create(
    _["numberOfEvents1"] = nevents1,
    _["numberOfDropouts1"] = ndropouts1,
    _["numberOfSubjects1"] = nsubjects1,
    _["numberOfMilestone1"] = nmilestone1,
    _["numberOfEvents2"] = nevents2,
    _["numberOfDropouts2"] = ndropouts2,
    _["numberOfSubjects2"] = nsubjects2,
    _["numberOfMilestone2"] = nmilestone2,
    _["expectedNumberOfEvents1"] = expectedNumberOfEvents1,
    _["expectedNumberOfDropouts1"] = expectedNumberOfDropouts1,
    _["expectedNumberOfSubjects1"] = expectedNumberOfSubjects1,
    _["expectedNumberOfMilestone1"] = expectedNumberOfMilestone1,
    _["expectedNumberOfEvents2"] = expectedNumberOfEvents2,
    _["expectedNumberOfDropouts2"] = expectedNumberOfDropouts2,
    _["expectedNumberOfSubjects2"] = expectedNumberOfSubjects2,
    _["expectedNumberOfMilestone2"] = expectedNumberOfMilestone2);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings,
    _["byTreatmentCounts"] = byTreatmentCounts);

  result.attr("class") = "kmpowerequiv";

  return result;
}


//' @title Sample Size for Equivalence in Milestone Survival Probability
//' Difference
//' @description Obtains the sample size for equivalence in milestone
//' survival probability difference.
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
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survDiffLower The lower equivalence limit of milestone survival
//'   probability difference.
//' @param survDiffUpper The upper equivalence limit of milestone survival
//'   probability difference.
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
//' @return An S3 class \code{kmpowerequiv} object
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{kmpowerequiv}}
//'
//' @examples
//'
//' kmsamplesizeequiv(beta = 0.1, kMax = 2, informationRates = c(0.5, 1),
//'                   alpha = 0.05, typeAlphaSpending = "sfOF",
//'                   milestone = 18,
//'                   survDiffLower = -0.13, survDiffUpper = 0.13,
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
List kmsamplesizeequiv(const double beta = 0.2,
                       const int kMax = 1,
                       const NumericVector& informationRates = NA_REAL,
                       const NumericVector& criticalValues = NA_REAL,
                       const double alpha = 0.05,
                       const std::string typeAlphaSpending = "sfOF",
                       const double parameterAlphaSpending = NA_REAL,
                       const NumericVector& userAlphaSpending = NA_REAL,
                       const double milestone = NA_REAL,
                       const double survDiffLower = NA_REAL,
                       const double survDiffUpper = NA_REAL,
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
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
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

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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

  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }

  if (milestone <= 0) {
    stop("milestone must be positive");
  }

  if (std::isnan(survDiffLower)) {
    stop("survDiffLower must be provided");
  }

  if (std::isnan(survDiffUpper)) {
    stop("survDiffUpper must be provided");
  }

  if (survDiffLower <= -1) {
    stop("survDiffLower must be greater than -1");
  }

  if (survDiffUpper >= 1) {
    stop("survDiffUpper must be less than 1");
  }

  if (survDiffLower >= survDiffUpper) {
    stop("survDiffLower must be less than survDiffUpper");
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

  if (!std::isnan(accrualDuration)) {
    if (accrualDuration <= 0) {
      stop("accrualDuration must be positive");
    }
  }

  if (!std::isnan(followupTime)) {
    if (fixedFollowup && followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }

    if (!fixedFollowup && followupTime < 0) {
      stop("followupTime must be non-negative for variable follow-up");
    }
  }

  if (fixedFollowup && std::isnan(followupTime)) {
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
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

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
  if (std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualDuration";
  } else if (!std::isnan(accrualDuration) && std::isnan(followupTime)) {
    unknown = "followupTime";
  } else if (!std::isnan(accrualDuration) && !std::isnan(followupTime)) {
    unknown = "accrualIntensity";
  } else {
    stop("accrualDuration and followupTime cannot be both missing");
  }


  NumericVector b = criticalValues1;
  NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);

  NumericVector survs1(nstrata), survs2(nstrata);
  IntegerVector l1 = Range(0, nintervals-1);
  NumericVector zerogam(nintervals);
  NumericVector u0(1, milestone);
  for (int h=0; h<nstrata; h++) {
    IntegerVector l = h*nintervals + l1;
    NumericVector lam1 = lambda1x[l];
    NumericVector lam2 = lambda2x[l];
    survs1[h] = patrisk(u0, piecewiseSurvivalTime, lam1, zerogam)[0];
    survs2[h] = patrisk(u0, piecewiseSurvivalTime, lam2, zerogam)[0];
  }

  double surv1 = sum(stratumFraction*survs1);
  double surv2 = sum(stratumFraction*survs2);
  double survDiff = surv1 - surv2;
  double theta10 = survDiffLower, theta20 = survDiffUpper;
  NumericVector theta(kMax, survDiff);

  List design = getDesignEquiv(
    beta, NA_REAL, theta10, theta20, survDiff,
    kMax, informationRates1, criticalValues1,
    alpha, asf, asfpar, userAlphaSpending, spendingTime1);

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
              DataFrame km = kmstat(
                u0, milestone, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                dur1, dur2, fixedFollowup);
              return sum(NumericVector(km[18])) - maxInformation;
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
                  DataFrame km = kmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    aval, followupTime, fixedFollowup);
                  return sum(NumericVector(km[18])) - maxInformation;
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
                  DataFrame km = kmstat(
                    u0, milestone, allocationRatioPlanned,
                    accrualTime, accrualIntensity,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup);
                  return sum(NumericVector(km[18])) - maxInformation;
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
    double n = std::ceil(n0 - 1.0e-12);

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
                    DataFrame km = kmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval, fixedFollowup);
                    return sum(NumericVector(km[18])) - maxInformation;
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
                    DataFrame km = kmstat(
                      u0, milestone, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup);
                    return sum(NumericVector(km[18])) - maxInformation;
                  };

        double lower = std::max(milestone - accrualDuration, 0.0) + 1.0e-6;
        double aval = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }


  List result = kmpowerequiv(
    kMax, informationRates1, criticalValues1,
    alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, milestone, survDiffLower, survDiffUpper,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  return result;
}
