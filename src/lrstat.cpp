#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;


// define the integrand functions for lrstat1
typedef struct {
  double hazardRatioH0;
  double allocationRatioPlanned;
  NumericVector accrualTime;
  NumericVector accrualIntensity;
  NumericVector piecewiseSurvivalTime;
  NumericVector lambda1;
  NumericVector lambda2;
  NumericVector gamma1;
  NumericVector gamma2;
  NumericVector zero;
  double rho1;
  double rho2;
  double phi;
  double accrualDuration0;
  double minFollowupTime0;
  double maxFollowupTime0;
} lrparams;


void f_uscore(double *x, int n, void *ex) {
  lrparams *param = (lrparams *) ex;
  NumericVector u0(n);
  for (int i=0; i<n; i++) {
    u0[i] = x[i];
  }
  NumericMatrix xatrisk = natrisk(
    u0, param->allocationRatioPlanned, param->accrualTime,
    param->accrualIntensity, param->piecewiseSurvivalTime,
    param->lambda1, param->lambda2, param->gamma1, param->gamma2,
    param->accrualDuration0, param->minFollowupTime0,
    param->maxFollowupTime0);
  NumericVector r1 = xatrisk(_, 0), r2 = xatrisk(_, 1);
  IntegerVector j = findInterval3(u0, param->piecewiseSurvivalTime) - 1;
  NumericVector w(n), N(n), lam1(n), lam2(n), d(n);
  if (param->rho1 != 0.0 || param->rho2 != 0.0) {
    NumericVector s1(n), s2(n), s(n);
    s1 = patrisk(u0, param->piecewiseSurvivalTime, param->lambda1,
                 param->zero);
    s2 = patrisk(u0, param->piecewiseSurvivalTime, param->lambda2,
                 param->zero);
    s = (param->phi)*s1 + (1.0 - (param->phi))*s2;
    w = pow(s, param->rho1)*pow(1.0-s, param->rho2);
  } else {
    w.fill(1.0);
  }
  N = (r1*(param->hazardRatioH0))*r2/(r1*(param->hazardRatioH0) + r2);
  lam1 = param->lambda1[j];
  lam2 = param->lambda2[j];
  d = lam1/(param->hazardRatioH0) - lam2;
  u0 = w*N*d;
  for (int i=0; i<n; i++) {
    x[i] = u0[i];
  }
}


void f_vscore(double *x, int n, void *ex) {
  lrparams *param = (lrparams *) ex;
  NumericVector u0(n);
  for (int i=0; i<n; i++) {
    u0[i] = x[i];
  }
  NumericMatrix xatrisk = natrisk(
    u0, param->allocationRatioPlanned, param->accrualTime,
    param->accrualIntensity, param->piecewiseSurvivalTime,
    param->lambda1, param->lambda2, param->gamma1, param->gamma2,
    param->accrualDuration0, param->minFollowupTime0,
    param->maxFollowupTime0);
  NumericVector r1 = xatrisk(_, 0), r2 = xatrisk(_, 1);
  IntegerVector j = findInterval3(u0, param->piecewiseSurvivalTime) - 1;
  NumericVector w(n), N(n), lam1(n), lam2(n), d(n);
  if (param->rho1 != 0.0 || param->rho2 != 0.0) {
    NumericVector s1(n), s2(n), s(n);
    s1 = patrisk(u0, param->piecewiseSurvivalTime, param->lambda1,
                 param->zero);
    s2 = patrisk(u0, param->piecewiseSurvivalTime, param->lambda2,
                 param->zero);
    s = (param->phi)*s1 + (1.0 - (param->phi))*s2;
    w = pow(s, param->rho1)*pow(1.0-s, param->rho2);
  } else {
    w.fill(1.0);
  }
  N = (r1*(param->hazardRatioH0))*r2/pow(r1*(param->hazardRatioH0) + r2, 2);
  lam1 = param->lambda1[j];
  lam2 = param->lambda2[j];
  d = r1*lam1 + r2*lam2;
  u0 = w*w*N*d;
  for (int i=0; i<n; i++) {
    x[i] = u0[i];
  }
}


void f_iscore(double *x, int n, void *ex) {
  lrparams *param = (lrparams *) ex;
  NumericVector u0(n);
  for (int i=0; i<n; i++) {
    u0[i] = x[i];
  }
  NumericMatrix xatrisk = natrisk(
    u0, param->allocationRatioPlanned, param->accrualTime,
    param->accrualIntensity, param->piecewiseSurvivalTime,
    param->lambda1, param->lambda2, param->gamma1, param->gamma2,
    param->accrualDuration0, param->minFollowupTime0,
    param->maxFollowupTime0);
  NumericVector r1 = xatrisk(_, 0), r2 = xatrisk(_, 1);
  IntegerVector j = findInterval3(u0, param->piecewiseSurvivalTime) - 1;
  NumericVector w(n), N(n), lam1(n), lam2(n), d(n);
  if (param->rho1 != 0.0 || param->rho2 != 0.0) {
    NumericVector s1(n), s2(n), s(n);
    s1 = patrisk(u0, param->piecewiseSurvivalTime, param->lambda1,
                 param->zero);
    s2 = patrisk(u0, param->piecewiseSurvivalTime, param->lambda2,
                 param->zero);
    s = (param->phi)*s1 + (1.0 - (param->phi))*s2;
    w = pow(s, param->rho1)*pow(1.0-s, param->rho2);
  } else {
    w.fill(1.0);
  }
  N = (r1*(param->hazardRatioH0))*r2/pow(r1*(param->hazardRatioH0) + r2, 2);
  lam1 = param->lambda1[j];
  lam2 = param->lambda2[j];
  d = r1*lam1 + r2*lam2;
  u0 = w*N*d;
  for (int i=0; i<n; i++) {
    x[i] = u0[i];
  }
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


  // obtain the follow-up time for the first enrolled subject
  double maxFollowupTime;
  if (fixedFollowup) {
    maxFollowupTime = followupTime;
  } else {
    maxFollowupTime = accrualDuration + followupTime;
  }

  IntegerVector q = Range(0, numSubintervals);
  NumericVector q1 = NumericVector(q);
  Range q2 = Range(0, numSubintervals-1);
  Range c0 = Range(0,0), c1 = Range(1,1);

  NumericVector ss(1, time);
  double a = accrual(ss, accrualTime, accrualIntensity, accrualDuration)[0];
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);

  double frac, accrualDuration0, minFollowupTime0, maxFollowupTime0, inc;
  IntegerVector l1 = Range(0, nintervals-1);
  IntegerVector l(nintervals);
  NumericVector lam1(nintervals), lam2(nintervals);
  NumericVector gam1(nintervals), gam2(nintervals);
  NumericMatrix x(1,2), y(1,2);
  NumericVector nsubjects(nstrata);
  NumericMatrix nevents(nstrata, 2), ndropouts(nstrata, 2);
  NumericMatrix xatrisk(numSubintervals+1, 2);
  NumericMatrix xevent(numSubintervals+1, 2);
  NumericVector atrisk1(numSubintervals), atrisk1x(numSubintervals);
  NumericVector atrisk2(numSubintervals);
  NumericVector atriskt(numSubintervals), atrisktx(numSubintervals);
  NumericVector event1(numSubintervals), event2(numSubintervals);
  NumericVector eventt(numSubintervals);
  NumericVector km(numSubintervals), w(numSubintervals);
  NumericVector uscore(nstrata), vscore(nstrata), iscore(nstrata);
  NumericVector nevents1(nstrata), nevents2(nstrata), neventst(nstrata);
  NumericVector ndropouts1(nstrata), ndropouts2(nstrata);
  NumericVector ndropoutst(nstrata);
  IntegerVector stratum(nstrata);
  NumericVector times(nstrata);
  DataFrame df;

  NumericVector maxFU(1, maxFollowupTime);
  NumericVector tt = NumericVector::create(time - maxFollowupTime);
  double a2 = accrual(tt, accrualTime, accrualIntensity, accrualDuration)[0];
  NumericVector nfmax1(nstrata), nfmax2(nstrata), nfmax(nstrata);

  bool dqagx = is_true(all(gamma1x == gamma2x)) || (rho1 == 0 && rho2 == 0);
  double tol = 1e-6;

  for (int h=0; h<nstrata; h++) {
    stratum[h] = h+1;
    times[h] = time;

    frac = stratumFraction[h];
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
    double p1 = patrisk(maxFU, piecewiseSurvivalTime, lam1, gam1)[0];
    double p2 = patrisk(maxFU, piecewiseSurvivalTime, lam2, gam2)[0];
    nfmax1[h] = phi*ncom*p1;
    nfmax2[h] = (1-phi)*ncom*p2;
    nfmax[h] = nfmax1[h] + nfmax2[h];

    // approximate the mean and variance of weighted log-rank test
    // score statistic
    if (!predictEventOnly) {

      // modify the study design at the calendar time of interest
      accrualDuration0 = std::min(time, accrualDuration);
      minFollowupTime0 = std::max(time - accrualDuration, 0.0);
      maxFollowupTime0 = std::min(time, maxFollowupTime);

      if (!dqagx) {
        // partition the follow-up period into small sub-intervals
        inc = maxFollowupTime0/numSubintervals;
        NumericVector t = q1*inc;

        // obtain number of subjects at risk and number of subjects having
        // an event at each analysis time point
        xatrisk = natrisk(t, allocationRatioPlanned,
                          accrualTime, frac*accrualIntensity,
                          piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
                          accrualDuration0, minFollowupTime0,
                          maxFollowupTime0);

        xevent = nevent(t, allocationRatioPlanned,
                        accrualTime, frac*accrualIntensity,
                        piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
                        accrualDuration0, minFollowupTime0,
                        maxFollowupTime0);

        // number of subjects at risk at start of each analysis time interval
        atrisk1 = xatrisk(q2, c0);
        atrisk2 = xatrisk(q2, c1);
        atrisk1x = hazardRatioH0*atrisk1; // adjust with HR under H0
        atriskt = atrisk1 + atrisk2;
        atrisktx = atrisk1x + atrisk2;

        // number of subjects having an event in each analysis time interval
        event1 = diff(xevent(_, 0));
        event2 = diff(xevent(_, 1));
        eventt = event1 + event2;

        // Kaplan-Meier estimates of survival probabilities at the start of
        // each analysis time interval
        if (!(rho1 == 0.0 && rho2 == 0.0)) {
          km[0] = 1;
          for (int i=1; i<numSubintervals; i++) {
            km[i] = km[i-1]*(1 - eventt[i-1]/atriskt[i-1]);
          }

          // vector of Fleming-Harrington weights
          w = pow(km,rho1)*pow(1-km,rho2);
        } else {
          w.fill(1.0);
        }

        // mean of the weighted log-rank test score statistic
        uscore[h] = sum(w * (event1 - eventt*atrisk1x/atrisktx));

        // variance of the weighted log-rank test score statistic
        vscore[h] = sum(w*w * eventt*atrisk1x*atrisk2/pow(atrisktx,2));

        // information of the weighted log-rank test score statistic
        iscore[h] = sum(w * eventt*atrisk1x*atrisk2/pow(atrisktx,2));
      } else {
        NumericVector zero(nintervals);
        lrparams param = {hazardRatioH0, allocationRatioPlanned,
                          accrualTime, frac*accrualIntensity,
                          piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
                          zero, rho1, rho2, phi, accrualDuration0,
                          minFollowupTime0, maxFollowupTime0};

        uscore[h] = quad(f_uscore, &param, 0.0, maxFollowupTime0, tol)[0];
        vscore[h] = quad(f_vscore, &param, 0.0, maxFollowupTime0, tol)[0];
        iscore[h] = quad(f_iscore, &param, 0.0, maxFollowupTime0, tol)[0];
      }
    }
  }

  // number of subjects having an event in each treatment group and overall
  nevents1 = nevents(_, 0);
  nevents2 = nevents(_, 1);
  neventst = nevents1 + nevents2;

  ndropouts1 = ndropouts(_, 0);
  ndropouts2 = ndropouts(_, 1);
  ndropoutst = ndropouts1 + ndropouts2;

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

  if (is_true(any(is_na(time)))) {
    stop("time must be provided");
  }

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
      logHR[j] = brent(g, -4.6, 4.6, 1.0e-6);
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

  if (is_true(any(is_na(nevents)))) {
    stop("nevents must be provided");
  }

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
    time[i] = brent(f, 0.0, studyTime, 1.0e-6);
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
//' @param followupTime Follow-up time for the last enrolled subjects.
//'   Must be provided for fixed follow-up design.
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
//' * \code{followupTime}: The follow-up time for the last enrolled subject.
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
    const double followupTime = NA_REAL,
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

  if (fixedFollowup && R_isnancpp(followupTime)) {
    stop("followupTime must be provided for fixed follow-up");
  }

  if (fixedFollowup && followupTime <= 0) {
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
                   NumericVector u0(1, t + 1000);
                   DataFrame lr = lrstat(
                     u0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, stratumFraction,
                     lambda1, lambda2, gamma1, gamma2,
                     t, 1000, fixedFollowup, 0, 0, 1, 1);
                   return sum(NumericVector(lr[2])) - nevents;
                 };

    t[0] = brent(fmin, interval[0], interval[1], 1.0e-6);

    auto fmax = [allocationRatioPlanned, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 fixedFollowup, nevents](double t)->double {
                   NumericVector u0(1, t);
                   DataFrame lr = lrstat(
                     u0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, stratumFraction,
                     lambda1, lambda2, gamma1, gamma2,
                     t, 0, fixedFollowup, 0, 0, 1, 1);
                   return sum(NumericVector(lr[2])) - nevents;
                 };

    t[1] = brent(fmax, t[0], interval[1], 1.0e-6);
  } else {
    auto fmin = [allocationRatioPlanned, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 followupTime, fixedFollowup, nevents](double t)->double {
                   NumericVector u0(1, t + followupTime);
                   DataFrame lr = lrstat(
                     u0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, stratumFraction,
                     lambda1, lambda2, gamma1, gamma2,
                     t, followupTime, fixedFollowup, 0, 0, 1, 1);
                   return sum(NumericVector(lr[2])) - nevents;
                 };

    t[0] = brent(fmin, interval[0], interval[1], 1.0e-6);

    auto fmax = [allocationRatioPlanned, accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 followupTime, fixedFollowup, nevents](double t)->double {
                   NumericVector u0(1, t);
                   DataFrame lr = lrstat(
                     u0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, stratumFraction,
                     lambda1, lambda2, gamma1, gamma2,
                     t, followupTime, fixedFollowup, 0, 0, 1, 1);
                   return sum(NumericVector(lr[2])) - nevents;
                 };

    t[1] = brent(fmax, t[0], interval[1], 1.0e-6);
  }

  NumericVector bigd(1, nevents);
  NumericVector ta(npoints), n(npoints), ts(npoints), tf(npoints);
  double dt = (t[1] - t[0])/(npoints - 1);

  for (int i=0; i<npoints; i++) {
    ta[i] = t[0] + i*dt;

    if (!fixedFollowup) {
      if (i==0) {
        ts[i] = ta[i] + 1000;
      } else if (i == npoints - 1){
        ts[i] = ta[i];
      } else {
        ts[i] = caltime(bigd, allocationRatioPlanned,
                        accrualTime, accrualIntensity,
                        piecewiseSurvivalTime, stratumFraction,
                        lambda1, lambda2, gamma1, gamma2,
                        ta[i], 1000, fixedFollowup)[0];
      }
      tf[i] = ts[i] - ta[i];
    } else {
      if (i==0) {
        ts[i] = ta[i] + followupTime;
      } else if (i == npoints - 1) {
        ts[i] = ta[i];
      } else {
        ts[i] = caltime(bigd, allocationRatioPlanned,
                        accrualTime, accrualIntensity,
                        piecewiseSurvivalTime, stratumFraction,
                        lambda1, lambda2, gamma1, gamma2,
                        ta[i], followupTime, fixedFollowup)[0];
      }
      tf[i] = followupTime;
    }
  }

  n = accrual(ta, accrualTime, accrualIntensity, 1000);

  DataFrame df = DataFrame::create(
    _["nevents"] = rep(nevents, npoints),
    _["fixedFollowup"] = rep(fixedFollowup, npoints),
    _["accrualDuration"] = ta,
    _["subjects"] = n,
    _["followupTime"] = tf,
    _["studyDuration"] = ts);

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
//'     - \code{followupTime}: The follow-up time.
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
             const std::string typeAlphaSpending = "sfOF",
             const double parameterAlphaSpending = NA_REAL,
             const NumericVector& userAlphaSpending = NA_REAL,
             const NumericVector& futilityBounds = NA_REAL,
             const std::string typeBetaSpending = "none",
             const double parameterBetaSpending = NA_REAL,
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
             const bool estimateHazardRatio = 1,
             const std::string typeOfComputation = "direct",
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

    for (int i=0; i<kMax-1; i++) {
      information1 = maxInformation*informationRates1[i];
      time[i] = brent(f, 1.0e-6, studyDuration1, 1.0e-6);
    };
    time[kMax-1] = studyDuration1;
  }


  // obtain mean and variance of log-rank test score statistic at each stage
  NumericVector theta(kMax), vscore(kMax);

  double phi = allocationRatioPlanned/(allocationRatioPlanned+1);

  if (su == "schoenfeld") {
    theta = rep(-log(hazardRatio/hazardRatioH0), kMax);

    vscore = phi*(1-phi)*e0;

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
  NumericVector nsubjects1 = phi*nsubjects;
  NumericVector nsubjects2 = (1-phi)*nsubjects;
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
    const std::string typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const NumericVector& userAlphaSpending = NA_REAL,
    const NumericVector& futilityBounds = NA_REAL,
    const std::string typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const NumericVector& userBetaSpending = NA_REAL,
    const NumericVector& spendingTime = NA_REAL,
    const double hazardRatioH0 = 1,
    const double hazardRatio = NA_REAL,
    const double allocationRatioPlanned = 1,
    const bool rounding = 1) {

  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (hazardRatioH0 <= 0) {
    stop("hazardRatioH0 must be positive");
  }

  if (R_isnancpp(hazardRatio)) {
    stop("hazardRatio must be provided");
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

  DataFrame byStageResults = DataFrame(design["byStageResults"]);
  NumericVector information = byStageResults["information"];
  double maxInformation = information[kMax-1];
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);
  double D = maxInformation/(phi*(1-phi));
  if (rounding) D = std::ceil(D);
  return D;
}


//' @title Log-rank test sample size
//' @description Obtains the needed accrual duration given power and
//' follow-up time, the needed follow-up time given power and
//' accrual duration, or the needed absolute accrual rates given
//' power, accrual duration, follow-up time, and relative accrual
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
//' # Example 1: Obtains accrual duration given power and follow-up time
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
//' # Example 2: Obtains follow-up time given power and accrual duration
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
//' # accrual duration, follow-up time, and relative accrual intensity
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
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const NumericVector& futilityBounds = NA_REAL,
                  const std::string typeBetaSpending = "none",
                  const double parameterBetaSpending = NA_REAL,
                  const NumericVector& userBetaSpending = NA_REAL,
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
                  double accrualDuration = NA_REAL,
                  double followupTime = NA_REAL,
                  const bool fixedFollowup = 0,
                  const double rho1 = 0,
                  const double rho2 = 0,
                  const int numSubintervals = 300,
                  const bool estimateHazardRatio = 1,
                  const std::string typeOfComputation = "direct",
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

  if (su == "schoenfeld") {
    double delta = -log(hazardRatio/hazardRatioH0);

    List design = getDesign(
      beta, NA_REAL, delta, kMax, informationRates1,
      efficacyStopping1, futilityStopping1, criticalValues1,
      alpha1, asf, asfpar, userAlphaSpending, futilityBounds1,
      bsf, bsfpar, userBetaSpending, spendingTime1, 1);

    DataFrame byStageResults = DataFrame(design["byStageResults"]);
    futilityBounds1 = byStageResults["futilityBounds"];

    DataFrame overallResults = DataFrame(design["overallResults"]);
    double maxInformation = overallResults["information"];
    double phi = allocationRatioPlanned/(allocationRatioPlanned+1);
    double D = maxInformation/(phi*(1-phi));

    auto f = [allocationRatioPlanned,
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
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  dur1, dur2, fixedFollowup, 0, 0, 1, 1);

                return sum(NumericVector(lr[2])) - D;
              };

    if (unknown == "accrualDuration") {
      accrualDuration = brent(f, interval[0], interval[1], 1.0e-6);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 1.0e-6);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 1.0e-6);
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
                  lr = lrstat(u0, 1, allocationRatioPlanned,
                              accrualTime, accrualIntensity1,
                              piecewiseSurvivalTime, stratumFraction,
                              lambda1, lambda2, gamma1, gamma2,
                              dur1, dur2, fixedFollowup,
                              rho1, rho2, numSubintervals, 1);

                  e0 = sum(NumericVector(lr[2]))*informationRates1;
                  time = caltime(e0, allocationRatioPlanned,
                                 accrualTime, accrualIntensity1,
                                 piecewiseSurvivalTime, stratumFraction,
                                 lambda1, lambda2, gamma1, gamma2,
                                 dur1, dur2, fixedFollowup);
                } else {
                  lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
                              accrualTime, accrualIntensity1,
                              piecewiseSurvivalTime, stratumFraction,
                              lambda1, lambda2, gamma1, gamma2,
                              dur1, dur2, fixedFollowup,
                              rho1, rho2, numSubintervals, 2);

                  double maxInformation = sum(NumericVector(lr[12]));
                  double information1;

                  auto g = [hazardRatioH0, allocationRatioPlanned,
                            accrualTime, accrualIntensity1,
                            piecewiseSurvivalTime, stratumFraction,
                            lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup,
                            rho1, rho2, numSubintervals,
                            &information1](double aval)->double {
                              NumericVector u0(1, aval);
                              DataFrame lr = lrstat(
                                u0, hazardRatioH0, allocationRatioPlanned,
                                accrualTime, accrualIntensity1,
                                piecewiseSurvivalTime, stratumFraction,
                                lambda1, lambda2, gamma1, gamma2,
                                dur1, dur2, fixedFollowup,
                                rho1, rho2, numSubintervals, 2);
                              return sum(NumericVector(lr[12])) -
                                information1;
                            };

                  for (int i=0; i<kMax-1; i++) {
                    information1 = maxInformation*informationRates1[i];
                    time[i] = brent(g, 1.0e-6, studyDuration1, 1.0e-6);
                  };
                  time[kMax-1] = studyDuration1;
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
                          g, -6.0, criticalValues1[k], 1.0e-6);
                      } else if (k < kMax-1) {
                        return -1.0;
                      }
                    }
                  }

                  return epsilon;
                }
              };

    if (unknown == "accrualDuration") {
      accrualDuration = brent(f, interval[0], interval[1], 1.0e-6);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 1.0e-6);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 1.0e-6);
      accrualIntensity1 = aval*accrualIntensity;
    }
  }

  futilityBounds1[kMax-1] = criticalValues1[kMax-1];


  // output the results
  List resultH1, resultH0, result;

  if (rounding) {
    NumericVector u0(1, accrualDuration + followupTime);
    DataFrame lr = lrstat(u0, 1, allocationRatioPlanned,
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
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      // adjust follow-up time to obtain integer number of events
      auto h = [allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, fixedFollowup,
                D](double aval)->double {
                  NumericVector u0(1, accrualDuration + aval);
                  DataFrame lr = lrstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, aval, fixedFollowup,
                    0, 0, 1, 1);
                  return sum(NumericVector(lr[2])) - D;
                };

      followupTime = brent(h, 0.0, 1.1*followupTime, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else {
      // adjust accrual intensity or duration to obtain int number of events
      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (D/D0)*accrualIntensity1;
      } else {
        auto h = [allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  followupTime, fixedFollowup,
                  D](double aval)->double {
                    NumericVector u0(1, aval + followupTime);
                    DataFrame lr = lrstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      aval, followupTime, fixedFollowup,
                      0, 0, 1, 1);
                    return sum(NumericVector(lr[2])) - D;
                  };

        accrualDuration = brent(h, accrualDuration, 1.1*accrualDuration,
                                1.0e-6);
      }

      NumericVector u0(1, accrualDuration);
      n0 = accrual(u0, accrualTime, accrualIntensity1, accrualDuration)[0];

      // round up the sample size
      n = std::ceil(n0);

      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      // adjust study duration to obtain integer number of events
      auto h = [allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                D](double aval)->double {
                  NumericVector u0(1, accrualDuration + aval);
                  DataFrame lr = lrstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup,
                    0, 0, 1, 1);
                  return sum(NumericVector(lr[2])) - D;
                };

      double aval = brent(h, 0.0, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    }


    // update information rates to calculate new boundaries
    NumericVector nevents(kMax), information(kMax), time(kMax);

    if (rho1 == 0 && rho2 == 0) {
      nevents = floor(D*informationRates1 + 0.5);
      informationRates1 = nevents/nevents[kMax-1];
    } else {
      u0[0] = studyDuration;
      lr = lrstat(u0, hazardRatioH0, allocationRatioPlanned,
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

      for (int i=0; i<kMax-1; i++) {
        information1 = maxInformation*informationRates1[i];
        time[i] = brent(f, 1.0e-6, studyDuration, 1.0e-6);
      };
      time[kMax-1] = studyDuration;

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
  DataFrame overallResults = DataFrame(resultH1["overallResults"]);
  DataFrame byStageResults = DataFrame(resultH1["byStageResults"]);
  double D = overallResults["numberOfEvents"];
  NumericVector information = byStageResults["information"];
  double maxInformation = information[kMax-1];
  double studyDuration;

  if (!fixedFollowup) {
    auto h = [hazardRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              lambda2, gamma1, gamma2,
              accrualDuration, fixedFollowup,
              rho1, rho2, numSubintervals,
              D, maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                if (rho1 == 0 && rho2 == 0) {
                  DataFrame lr = lrstat(
                    u0, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                    accrualDuration, aval, fixedFollowup,
                    rho1, rho2, numSubintervals, 1);
                  return sum(NumericVector(lr[2])) - D;
                } else {
                  DataFrame lr = lrstat(
                    u0, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                    accrualDuration, aval, fixedFollowup,
                    rho1, rho2, numSubintervals, 2);
                  return sum(NumericVector(lr[12])) - maxInformation;
                }
              };

    if (h(0) < 0) { // adjust the follow-up time
      double lower = 0.0, upper = followupTime;
      while (h(upper) < 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      followupTime = brent(h, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // adjust the accrual duration
      auto g = [hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda2, gamma1, gamma2,
                fixedFollowup,
                rho1, rho2, numSubintervals,
                D, maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  if (rho1 == 0 && rho2 == 0) {
                    DataFrame lr = lrstat(
                      u0, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval, 0, fixedFollowup,
                      rho1, rho2, numSubintervals, 1);
                    return sum(NumericVector(lr[2])) - D;
                  } else {
                    DataFrame lr = lrstat(
                      u0, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval, 0, fixedFollowup,
                      rho1, rho2, numSubintervals, 2);
                    return sum(NumericVector(lr[12])) - maxInformation;
                  }
                };

      accrualDuration = brent(g, 1.0e-6, accrualDuration, 1.0e-6);
      followupTime = 0.0;
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
                NumericVector u0(1, accrualDuration + aval);
                if (rho1 == 0 && rho2 == 0) {
                  DataFrame lr = lrstat(
                    u0, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup,
                    rho1, rho2, numSubintervals, 1);
                  return sum(NumericVector(lr[2])) - D;
                } else {
                  DataFrame lr = lrstat(
                    u0, hazardRatioH0, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup,
                    rho1, rho2, numSubintervals, 2);
                  return sum(NumericVector(lr[12])) - maxInformation;
                }
              };

    if (h(followupTime) < 0) { // increase the accrual duration
      auto g = [hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                rho1, rho2, numSubintervals,
                D, maxInformation](double aval)->double {
                  NumericVector u0(1, aval + followupTime);
                  if (rho1 == 0 && rho2 == 0) {
                    DataFrame lr = lrstat(
                      u0, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval, followupTime, fixedFollowup,
                      rho1, rho2, numSubintervals, 1);
                    return sum(NumericVector(lr[2])) - D;
                  } else {
                    DataFrame lr = lrstat(
                      u0, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval, followupTime, fixedFollowup,
                      rho1, rho2, numSubintervals, 2);
                    return sum(NumericVector(lr[12])) - maxInformation;
                  }
                };

      double lower = accrualDuration, upper = 2.0*lower;
      while (g(upper) < 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(g, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (h(0) < 0) { // decrease the study duration
      double aval = brent(h, 0.0, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    } else { // decrease the accrual duration
      auto g = [hazardRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                rho1, rho2, numSubintervals,
                D, maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  if (rho1 == 0 && rho2 == 0) {
                    DataFrame lr = lrstat(
                      u0, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval, followupTime, fixedFollowup,
                      rho1, rho2, numSubintervals, 1);
                    return sum(NumericVector(lr[2])) - D;
                  } else {
                    DataFrame lr = lrstat(
                      u0, hazardRatioH0, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda2*hazardRatioH0, lambda2, gamma1, gamma2,
                      aval, followupTime, fixedFollowup,
                      rho1, rho2, numSubintervals, 2);
                    return sum(NumericVector(lr[12])) - maxInformation;
                  }
                };

      accrualDuration = brent(g, 1.0e-6, accrualDuration, 1.0e-6);
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


//' @title Power for equivalence in hazard ratio
//' @description Obtains the power for equivalence in hazard ratio.
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
//' @param hazardRatioLower The lower equivalence limit of hazard ratio.
//' @param hazardRatioUpper The upper equivalence limit of hazard ratio.
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
//' @return An S3 class \code{lrpowerequiv} object with 4 components:
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
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{hazardRatioLower}: The lower equivalence limit of hazard
//'       ratio.
//'
//'     - \code{hazardRatioUpper}: The upper equivalence limit of hazard
//'       ratio.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up time.
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
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyHRLower}: The efficacy boundaries on the
//'       hazard ratio scale for the one-sided null hypothesis
//'       at the lower equivalence limit.
//'
//'     - \code{efficacyHRUpper}: The efficacy boundaries on the
//'       hazard ratio scale for the one-sided null hypothesis
//'       at the upper equivalence limit.
//'
//'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
//'       each of the two one-sided tests.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{HR}: The average hazard ratio.
//'
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{allocationRatioPlanned},
//'   \code{accrualTime}, \code{accuralIntensity},
//'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   \code{typeOfComputation}, and \code{spendingTime}.
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
//' @seealso \code{\link{rmstat}}
//'
//' @examples
//'
//' lrpowerequiv(kMax = 2, informationRates = c(0.5, 1),
//'              alpha = 0.05, typeAlphaSpending = "sfOF",
//'              hazardRatioLower = 0.71, hazardRatioUpper = 1.4,
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
List lrpowerequiv(const int kMax = 1,
                  const NumericVector& informationRates = NA_REAL,
                  const NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.05,
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const double hazardRatioLower = NA_REAL,
                  const double hazardRatioUpper = NA_REAL,
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
                  const std::string typeOfComputation = "direct",
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

  if (R_isnancpp(hazardRatioLower)) {
    stop("hazardRatioLower must be provided");
  }

  if (R_isnancpp(hazardRatioUpper)) {
    stop("hazardRatioUpper must be provided");
  }

  if (hazardRatioLower <= 0) {
    stop("hazardRatioLower must be positive");
  }

  if (hazardRatioLower >= hazardRatioUpper) {
    stop("hazardRatioLower must be less than hazardRatioUpper");
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

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration < accrualDuration) {
    stop("studyDuration must be greater than or equal to accrualDuration");
  }

  if (fixedFollowup && !R_isnancpp(studyDuration) &&
      studyDuration > accrualDuration + followupTime) {
    stop("studyDuration cannot exceed accrualDuration + followupTime");
  }

  std::string su = typeOfComputation;
  std::for_each(su.begin(), su.end(), [](char & c) {
    c = std::tolower(c);
  });

  if (su != "direct" && su != "schoenfeld") {
    stop("typeOfComputation must be direct or Schoenfeld");
  }

  double hazardRatio = 1.0;
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
  NumericVector e0(kMax), time(kMax);
  NumericVector u0(1, studyDuration1);
  DataFrame lr = lrstat(u0, 1.0, allocationRatioPlanned,
                        accrualTime, accrualIntensity,
                        piecewiseSurvivalTime, stratumFraction,
                        lambda1, lambda2, gamma1, gamma2,
                        accrualDuration, followupTime, fixedFollowup,
                        0, 0, 1, 1);

  e0 = sum(NumericVector(lr[2]))*informationRates1;
  time = caltime(e0, allocationRatioPlanned,
                 accrualTime, accrualIntensity,
                 piecewiseSurvivalTime, stratumFraction,
                 lambda1, lambda2, gamma1, gamma2,
                 accrualDuration, followupTime, fixedFollowup);

  double phi = allocationRatioPlanned/(1.0 + allocationRatioPlanned);
  NumericVector HR(kMax), theta(kMax), I(kMax);

  if (su == "schoenfeld") {
    HR = rep(hazardRatio, kMax);
    theta = log(HR);
    I = phi*(1-phi)*e0;

    lr = lrstat(time, 1, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                0, 0, 1, 1);
  } else {
    lr = lrstat(time, 1, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                0, 0, 1, 3);

    HR = NumericVector(lr[15]);
    theta = log(HR);
    I = 1.0/NumericVector(lr[16]);
  }

  NumericVector nsubjects = NumericVector(lr[1]);
  NumericVector nsubjects1 = phi*nsubjects;
  NumericVector nsubjects2 = (1-phi)*nsubjects;
  NumericVector nevents = NumericVector(lr[2]);
  NumericVector nevents1 = NumericVector(lr[3]);
  NumericVector nevents2 = NumericVector(lr[4]);
  NumericVector ndropouts = NumericVector(lr[5]);
  NumericVector ndropouts1 = NumericVector(lr[6]);
  NumericVector ndropouts2 = NumericVector(lr[7]);


  // calculate cumulative rejection probability under H1
  double theta10 = log(hazardRatioLower), theta20 = log(hazardRatioUpper);
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

  NumericVector efficacyHRLower = exp(theta10 + b/sqrt(I));
  NumericVector efficacyHRUpper = exp(theta20 - b/sqrt(I));

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

  double overallReject = cp[kMax-1];
  double attainedAlphaH10 = cpH10[kMax-1];
  double attainedAlphaH20 = cpH20[kMax-1];
  double expectedNumberOfEvents = sum(q*nevents);
  double expectedNumberOfDropouts = sum(q*ndropouts);
  double expectedNumberOfSubjects = sum(q*nsubjects);
  double expectedNumberOfEvents1 = sum(q*nevents1);
  double expectedNumberOfDropouts1 = sum(q*ndropouts1);
  double expectedNumberOfSubjects1 = sum(q*nsubjects1);
  double expectedNumberOfEvents2 = sum(q*nevents2);
  double expectedNumberOfDropouts2 = sum(q*ndropouts2);
  double expectedNumberOfSubjects2 = sum(q*nsubjects2);
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
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = (I[kMax-1]),
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfDropouts"] = expectedNumberOfDropouts,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["kMax"] = kMax,
    _["hazardRatioLower"] = hazardRatioLower,
    _["hazardRatioUpper"] = hazardRatioUpper,
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
    _["analysisTime"] = time,
    _["efficacyHRLower"] = efficacyHRLower,
    _["efficacyHRUpper"] = efficacyHRUpper,
    _["efficacyP"] = efficacyP,
    _["information"] = I,
    _["HR"] = hazardRatio);

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
    _["typeOfComputation"] = typeOfComputation,
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

  result.attr("class") = "lrpowerequiv";

  return result;
}


//' @title Sample size for equivalence in hazard ratio
//' @description Obtains the sample size for equivalence in hazard ratio.
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
//' @param hazardRatioLower The lower equivalence limit of hazard ratio.
//' @param hazardRatioUpper The upper equivalence limit of hazard ratio.
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
//' @inheritParams param_typeOfComputation
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return An S3 class \code{lrpowerequiv} object
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{lrpowerequiv}}
//'
//' @examples
//'
//' lrsamplesizeequiv(kMax = 2, informationRates = c(0.5, 1),
//'                   alpha = 0.05, typeAlphaSpending = "sfOF",
//'                   hazardRatioLower = 0.71, hazardRatioUpper = 1.4,
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
List lrsamplesizeequiv(const double beta = 0.2,
                       const int kMax = 1,
                       const NumericVector& informationRates = NA_REAL,
                       const NumericVector& criticalValues = NA_REAL,
                       const double alpha = 0.05,
                       const std::string typeAlphaSpending = "sfOF",
                       const double parameterAlphaSpending = NA_REAL,
                       const NumericVector& userAlphaSpending = NA_REAL,
                       const double hazardRatioLower = NA_REAL,
                       const double hazardRatioUpper = NA_REAL,
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
                       const std::string typeOfComputation = "direct",
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

  if (R_isnancpp(hazardRatioLower)) {
    stop("hazardRatioLower must be provided");
  }

  if (R_isnancpp(hazardRatioUpper)) {
    stop("hazardRatioUpper must be provided");
  }

  if (hazardRatioLower <= 0) {
    stop("hazardRatioLower must be positive");
  }

  if (hazardRatioLower >= hazardRatioUpper) {
    stop("hazardRatioLower must be less than hazardRatioUpper");
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

  std::string su = typeOfComputation;
  std::for_each(su.begin(), su.end(), [](char & c) {
    c = std::tolower(c);
  });

  if (su != "direct" && su != "schoenfeld") {
    stop("typeOfComputation must be direct or Schoenfeld");
  }

  double hazardRatio = 1.0;
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


  double phi = allocationRatioPlanned/(1.0 + allocationRatioPlanned);
  NumericVector b = criticalValues1;
  NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);
  double theta10 = log(hazardRatioLower), theta20 = log(hazardRatioUpper);

  if (su == "schoenfeld") {
    List design = getDesignEquiv(
      beta, NA_REAL, theta10, theta20, log(hazardRatio),
      kMax, informationRates1, criticalValues1,
      alpha, asf, asfpar, userAlphaSpending, spendingTime1,
      1, 1, 1, 1);

    DataFrame overallResults = DataFrame(design["overallResults"]);
    double maxInformation = overallResults["information"];
    double D = maxInformation/(phi*(1-phi));

    auto f = [allocationRatioPlanned,
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
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  dur1, dur2, fixedFollowup, 0, 0, 1, 1);

                return sum(NumericVector(lr[2])) - D;
              };

    if (unknown == "accrualDuration") {
      accrualDuration = brent(f, interval[0], interval[1], 1.0e-6);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 1.0e-6);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 1.0e-6);
      accrualIntensity1 = aval*accrualIntensity;
    }
  } else {
    auto f = [beta, kMax, informationRates1, criticalValues1,
              theta10, theta20, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              unknown, li, ui](double aval)->double {
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

                DataFrame lr;
                NumericVector e0(kMax), time(kMax);

                double studyDuration1 = dur1 + dur2;
                NumericVector u0(1, studyDuration1);

                // obtain the timing of interim analysis
                lr = lrstat(u0, 1, allocationRatioPlanned,
                            accrualTime, accrualIntensity,
                            piecewiseSurvivalTime, stratumFraction,
                            lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup,
                            0, 0, 1, 1);

                e0 = sum(NumericVector(lr[2]))*informationRates1;
                time = caltime(e0, allocationRatioPlanned,
                               accrualTime, accrualIntensity,
                               piecewiseSurvivalTime, stratumFraction,
                               lambda1, lambda2, gamma1, gamma2,
                               dur1, dur2, fixedFollowup);

                // obtain the mean and information at each stage
                lr = lrstat(time, 1, allocationRatioPlanned,
                            accrualTime, accrualIntensity1,
                            piecewiseSurvivalTime, stratumFraction,
                            lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup,
                            0, 0, 1, 3);

                NumericVector HR = NumericVector(lr[15]);
                NumericVector theta = log(HR);
                NumericVector I = 1.0/NumericVector(lr[16]);

                NumericVector b = criticalValues1;
                NumericVector l = b + theta10*sqrt(I);
                NumericVector u = -b + theta20*sqrt(I);

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
      accrualDuration = brent(f, interval[0], interval[1], 1.0e-6);
    } else if (unknown == "followupTime") {
      followupTime = brent(f, interval[0], interval[1], 1.0e-6);
    } else if (unknown == "accrualIntensity") {
      double aval = brent(f, interval[0], interval[1], 1.0e-6);
      accrualIntensity1 = aval*accrualIntensity;
    }
  }


  // output the results
  List result;
  if (rounding) {
    NumericVector u0(1, accrualDuration + followupTime);
    DataFrame lr = lrstat(u0, 1, allocationRatioPlanned,
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
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      // adjust follow-up time to obtain integer number of events
      auto h = [allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, fixedFollowup,
                D](double aval)->double {
                  NumericVector u0(1, accrualDuration + aval);
                  DataFrame lr = lrstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, aval, fixedFollowup,
                    0, 0, 1, 1);
                  return sum(NumericVector(lr[2])) - D;
                };

      followupTime = brent(h, 0.0, 1.1*followupTime, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else {
      // adjust accrual intensity or duration to obtain int number of events
      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (D/D0)*accrualIntensity1;
      } else {
        auto h = [allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  lambda1, lambda2, gamma1, gamma2,
                  followupTime, fixedFollowup,
                  D](double aval)->double {
                    NumericVector u0(1, aval + followupTime);
                    DataFrame lr = lrstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      lambda1, lambda2, gamma1, gamma2,
                      aval, followupTime, fixedFollowup,
                      0, 0, 1, 1);
                    return sum(NumericVector(lr[2])) - D;
                  };

        accrualDuration = brent(h, accrualDuration, 1.1*accrualDuration,
                                1.0e-6);
      }

      NumericVector u0(1, accrualDuration);
      n0 = accrual(u0, accrualTime, accrualIntensity1, accrualDuration)[0];

      // round up the sample size
      n = std::ceil(n0);

      if (unknown == "accrualIntensity") {
        accrualIntensity1 = (n/n0)*accrualIntensity1;
      } else {
        NumericVector ns(1, n);
        accrualDuration = getAccrualDurationFromN(ns, accrualTime,
                                                  accrualIntensity1)[0];
      }

      // adjust study duration to obtain integer number of events
      auto h = [allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup,
                D](double aval)->double {
                  NumericVector u0(1, accrualDuration + aval);
                  DataFrame lr = lrstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    lambda1, lambda2, gamma1, gamma2,
                    accrualDuration, followupTime, fixedFollowup,
                    0, 0, 1, 1);
                  return sum(NumericVector(lr[2])) - D;
                };

      double aval = brent(h, 0.0, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    }

    // update information rates to calculate new boundaries
    NumericVector nevents = floor(D*informationRates1 + 0.5);
    informationRates1 = nevents/nevents[kMax-1];

    // recalculate boundaries
    result = lrpowerequiv(
        kMax, informationRates1, criticalValues,
        alpha, typeAlphaSpending, parameterAlphaSpending,
        userAlphaSpending, hazardRatioLower, hazardRatioUpper,
        allocationRatioPlanned, accrualTime, accrualIntensity1,
        piecewiseSurvivalTime, stratumFraction,
        lambda1, lambda2, gamma1, gamma2,
        accrualDuration, followupTime, fixedFollowup,
        typeOfComputation, spendingTime, studyDuration);
  } else {
    double studyDuration = accrualDuration + followupTime;

    result = lrpowerequiv(
      kMax, informationRates1, criticalValues1,
      alpha, typeAlphaSpending, parameterAlphaSpending,
      userAlphaSpending, hazardRatioLower, hazardRatioUpper,
      allocationRatioPlanned, accrualTime, accrualIntensity1,
      piecewiseSurvivalTime, stratumFraction,
      lambda1, lambda2, gamma1, gamma2,
      accrualDuration, followupTime, fixedFollowup,
      typeOfComputation, spendingTime, studyDuration);
  }

  return result;
}


//' @title Log-rank test of survival curve difference
//' @description Obtains the log-rank test using the Fleming-Harrington
//' family of weights.
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
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//'
//' @return A data frame with the following variables:
//'
//' * \code{rep}: The replication.
//'
//' * \code{uscore}: The numerator of the log-rank test statistic.
//'
//' * \code{vscore}: The variance of the log-rank score test statistic.
//'
//' * \code{logRankZ}: The Z-statistic value.
//'
//' * \code{logRankPValue}: The one-sided p-value.
//'
//' * \code{rho1}: The first parameter of the Fleming-Harrington weights.
//'
//' * \code{rho2}: The second parameter of the Fleming-Harrington weights.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' df <- lrtest(data = rawdata, rep = "iterationNumber",
//'              stratum = "stratum", treat = "treatmentGroup",
//'              time = "timeUnderObservation", event = "event",
//'              rho1 = 0.5, rho2 = 0)
//' head(df)
//'
//' @export
// [[Rcpp::export]]
DataFrame lrtest(const DataFrame data,
                 const std::string rep = "rep",
                 const std::string stratum = "stratum",
                 const std::string treat = "treat",
                 const std::string time = "time",
                 const std::string event = "event",
                 const double rho1 = 0,
                 const double rho2 = 0) {

  int h, i, j, k, n = data.nrows();

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


  // create the numeric treatment variable
  IntegerVector treatn(n);
  IntegerVector treatwn;
  CharacterVector treatwc;
  if (TYPEOF(data[treat]) == LGLSXP) {
    LogicalVector treatv = data[treat];
    treatn = 2 - treatv;
    treatwn = IntegerVector::create(1,0);
  } else if (TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    IntegerVector treatw(n);

    // sort the treatment variable
    IntegerVector order = seq(0, n-1);
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return treatv[i] < treatv[j];
    });

    treatw = treatv[order];

    // identify the locations of the unique values
    IntegerVector idx(1,0);
    for (i=1; i<n; i++) {
      if (treatw[i] != treatw[i-1]) {
        idx.push_back(i);
      }
    }

    treatwn = treatw[idx];

    // check whether there are only two treatment values
    if (idx.size() != 2) {
      stop("treat must have two and only two distinct values");
    }

    // code the treatment variable
    for (i=0; i<n; i++) {
      treatn[i] = treatv[i] == treatw[idx[0]] ? 1 : 2;
    }
  } else if (TYPEOF(data[treat]) == STRSXP) {
    CharacterVector treatv = data[treat];
    CharacterVector treatw(n);

    // sort the treatment variable
    IntegerVector order = seq(0, n-1);
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return treatv[i] < treatv[j];
    });

    treatw = treatv[order];

    // identify the locations of the unique values
    IntegerVector idx(1,0);
    for (i=1; i<n; i++) {
      if (treatw[i] != treatw[i-1]) {
        idx.push_back(i);
      }
    }

    treatwc = treatw[idx];

    // check whether there are only two treatment values
    if (idx.size() != 2) {
      stop("treat must have two and only two distinct values");
    }

    // code the treatment variable
    for (i=0; i<n; i++) {
      treatn[i] = treatv[i] == treatw[idx[0]] ? 1 : 2;
    }
  } else {
    stop("incorrect type for the treatment variable in the input data");
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
  IntegerVector repwn;
  CharacterVector repwc;
  if (!has_rep) {
    repn.fill(1);
  } else {
    if (TYPEOF(data[rep]) == INTSXP) {
      IntegerVector repv = data[rep];
      IntegerVector repw(n);

      // sort the rep variable
      IntegerVector order = seq(0, n-1);
      std::sort(order.begin(), order.end(), [&](int i, int j) {
        return repv[i] < repv[j];
      });

      repw = repv[order];

      // identify the locations of the unique values
      IntegerVector idx(1,0);
      for (i=1; i<n; i++) {
        if (repw[i] != repw[i-1]) {
          idx.push_back(i);
        }
      }

      repwn = repw[idx]; // unique numeric values

      // code the rep variable
      for (i=0; i<n; i++) {
        for (j=0; j<idx.size(); j++) {
          if (repv[i] == repwn[j]) {
            repn[i] = j+1;
            break;
          }
        }
      }
    } else if (TYPEOF(data[rep]) == STRSXP) {
      CharacterVector repv = data[rep];
      CharacterVector repw(n);

      // sort the rep variable
      IntegerVector order = seq(0, n-1);
      std::sort(order.begin(), order.end(), [&](int i, int j) {
        return repv[i] < repv[j];
      });

      repw = repv[order];

      // identify the locations of the unique values
      IntegerVector idx(1,0);
      for (i=1; i<n; i++) {
        if (repw[i] != repw[i-1]) {
          idx.push_back(i);
        }
      }

      repwc = repw[idx]; // unique character values

      // code the rep variable
      for (i=0; i<n; i++) {
        for (j=0; j<idx.size(); j++) {
          if (repv[i] == repwc[j]) {
            repn[i] = j+1;
            break;
          }
        }
      }
    } else {
      stop("incorrect type for the replication variable in the input data");
    }
  }


  // create the numeric stratum variable
  IntegerVector stratumn(n);
  IntegerVector stratumwn;
  CharacterVector stratumwc;
  if (!has_stratum) {
    stratumn.fill(1);
  } else {
    if (TYPEOF(data[stratum]) == INTSXP) {
      IntegerVector stratumv = data[stratum];
      IntegerVector stratumw(n);

      // sort the stratum variable
      IntegerVector order = seq(0, n-1);
      std::sort(order.begin(), order.end(), [&](int i, int j) {
        return stratumv[i] < stratumv[j];
      });

      stratumw = stratumv[order];

      // identify the locations of the unique values
      IntegerVector idx(1,0);
      for (i=1; i<n; i++) {
        if (stratumw[i] != stratumw[i-1]) {
          idx.push_back(i);
        }
      }

      stratumwn = stratumw[idx]; // unique numeric values

      // code the stratum variable
      for (i=0; i<n; i++) {
        for (j=0; j<idx.size(); j++) {
          if (stratumv[i] == stratumwn[j]) {
            stratumn[i] = j+1;
            break;
          }
        }
      }
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      CharacterVector stratumv = data[stratum];
      CharacterVector stratumw(n);

      // sort the stratum variable
      IntegerVector order = seq(0, n-1);
      std::sort(order.begin(), order.end(), [&](int i, int j) {
        return stratumv[i] < stratumv[j];
      });

      stratumw = stratumv[order];

      // identify the locations of the unique values
      IntegerVector idx(1,0);
      for (i=1; i<n; i++) {
        if (stratumw[i] != stratumw[i-1]) {
          idx.push_back(i);
        }
      }

      stratumwc = stratumw[idx]; // unique character values

      // code the stratum variable
      for (i=0; i<n; i++) {
        for (j=0; j<idx.size(); j++) {
          if (stratumv[i] == stratumwc[j]) {
            stratumn[i] = j+1;
            break;
          }
        }
      }
    } else {
      stop("incorrect type for the stratum variable in the input data");
    }
  }


  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }

  if (rho2 < 0) {
    stop("rho2 must be non-negative");
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
  NumericVector uscore0(nreps), vscore0(nreps);
  NumericVector logRankZ0(nreps), logRankPValue0(nreps);

  bool noerr = 1;
  int index = 0;
  for (h=0; h<nreps; h++) {
    bool skip = 0;
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = q1.size();

    IntegerVector stratum1 = stratumn[q1];
    IntegerVector treat1 = treatn[q1];
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
    treat1 = treat1[order1];
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

    double uscore = 0.0, vscore = 0.0;
    for (i=0; i<nstrata; i++) {
      IntegerVector q = Range(idx1[i], idx1[i+1]-1);
      IntegerVector treat2 = treat1[q];
      NumericVector time2 = time1[q];
      NumericVector event2 = event1[q];
      int n2 = q.size();

      // running index of number of subjects at risk
      int nriskx = n2;
      int nrisk1x = sum(treat2 == 1), nrisk2x = sum(treat2 == 2);

      if ((nrisk1x == 0) || (nrisk2x == 0)) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else if (TYPEOF(data[rep]) == INTSXP) {
          reperr = " " + rep + " = " + std::to_string(repwn[repn[idx[h]]-1]);
        } else {
          reperr = " " + rep + " = " + repwc[repn[idx[h]]-1];
        }

        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else if (TYPEOF(data[stratum]) == INTSXP) {
          stratumerr = " " + stratum + " = " +
            std::to_string(stratumwn[stratum1[idx1[i]]-1]);
        } else {
          stratumerr = " " + stratum + " = " +
            stratumwc[stratum1[idx1[i]]-1];
        }


        k = nrisk1x == 0 ? 0 : 1;
        std::string treaterr;
        if ((TYPEOF(data[treat]) == LGLSXP) ||
            (TYPEOF(data[treat]) == INTSXP)) {
          treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
        } else {
          treaterr = " " + treat + " = " + treatwc[k];
        }

        std::string str1 = "Warning: The data set does not contain";
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


      double nrisk = nrisk1x, nrisk1 = nrisk1x, nrisk2 = nrisk2x;
      double nevent, nevent1, nevent2;
      double surv = 1.0;

      NumericVector nrisk0(n2, NA_REAL), nevent0(n2);
      NumericVector nrisk10(n2), nevent10(n2);
      NumericVector nrisk20(n2), nevent20(n2);
      NumericVector surv0(n2);

      int index1 = 0;
      bool cache = 0;
      for (j=0; j<n2; j++) {
        if (((j == 0) && (event2[j] == 1)) ||
            ((j >= 1) && (event2[j] == 1) && (time2[j] > time2[j-1]))) {
          // new event
          // add the info for the previous event
          if (cache) {
            surv = surv*(1.0 - nevent/nrisk);

            nrisk0[index1] = nrisk;
            nevent0[index1] = nevent;
            nrisk10[index1] = nrisk1;
            nevent10[index1] = nevent1;
            nrisk20[index1] = nrisk2;
            nevent20[index1] = nevent2;
            surv0[index1] = surv;

            index1++;
          }

          // update the cache for the current event time
          nrisk = nriskx;
          nrisk1 = nrisk1x;
          nrisk2 = nrisk2x;
          nevent = 1;
          nevent1 = (treat2[j] == 1);
          nevent2 = (treat2[j] == 2);

          cache = 1;
        } else if ((j >= 1) && (event2[j] == 1) && (event2[j-1] == 1) &&
          (time2[j] == time2[j-1])) { // tied event
          nevent = nevent + 1;
          nevent1 = nevent1 + (treat2[j] == 1);
          nevent2 = nevent2 + (treat2[j] == 2);
        } else if ((j >= 1) && (event2[j] == 0) && (event2[j-1] == 1)) {
          // new censoring
          // add the info for the previous event
          surv = surv*(1.0 - nevent/nrisk);

          nrisk0[index1] = nrisk;
          nevent0[index1] = nevent;
          nrisk10[index1] = nrisk1;
          nevent10[index1] = nevent1;
          nrisk20[index1] = nrisk2;
          nevent20[index1] = nevent2;
          surv0[index1] = surv;

          index1++;

          // empty the cache for the current event time
          cache = 0;
        }

        nriskx--;
        if (treat2[j] == 1) nrisk1x--;
        if (treat2[j] == 2) nrisk2x--;
      }

      // add the info for the last event
      if (cache) {
        surv = surv*(1.0 - nevent/nrisk);

        nrisk0[index1] = nrisk;
        nevent0[index1] = nevent;
        nrisk10[index1] = nrisk1;
        nevent10[index1] = nevent1;
        nrisk20[index1] = nrisk2;
        nevent20[index1] = nevent2;
        surv0[index1] = surv;

        index1++;
      }

      // only keep nonmissing records
      double uscore1 = 0.0, vscore1 = 0.0;
      LogicalVector sub = !is_na(nrisk0);
      if (is_true(any(sub))) { // at least 1 event
        nrisk0 = nrisk0[sub];
        nevent0 = nevent0[sub];
        nrisk10 = nrisk10[sub];
        nevent10 = nevent10[sub];
        nrisk20 = nrisk20[sub];
        nevent20 = nevent20[sub];
        surv0 = surv0[sub];

        int K = sum(sub);
        NumericVector w(K);
        surv0.push_front(1.0);
        for (k=0; k<K; k++) {
          double w1, w2;
          if (surv0[k] == 1.0) {
            w1 = 1.0;
            w2 = rho2 > 0.0 ? 0.0 : 1.0;
          } else if (surv0[k] == 0.0) {
            w1 = rho1 > 0.0 ? 0.0 : 1.0;
            w2 = 1.0;
          } else {
            w1 = pow(surv0[k], rho1);
            w2 = pow(1.0 - surv0[k], rho2);
          }
          w[k] = w1*w2;
        }

        for (k=0; k<K; k++) {
          uscore1 += w[k]*(nevent10[k] - nevent0[k]*nrisk10[k]/nrisk0[k]);
          if (nrisk0[k] > 1.0) {
            vscore1 += pow(w[k],2)*nevent0[k]*(nrisk0[k]-nevent0[k])*
              nrisk10[k]*nrisk20[k]/(pow(nrisk0[k],2)*(nrisk0[k]-1.0));
          }
        }
      }

      uscore += uscore1;
      vscore += vscore1;
    }

    // skip the replication if there is a stratum without both treatments
    if (skip) continue;

    rep0[index] = repn[idx[h]];
    uscore0[index] = uscore;
    vscore0[index] = vscore;
    logRankZ0[index] = uscore/sqrt(vscore);
    logRankPValue0[index] = R::pnorm(logRankZ0[index], 0, 1, 1, 0);

    index++;
  }

  // only keep nonmissing records
  LogicalVector sub = !is_na(rep0);
  if (is_false(any(sub))) {
    stop("no replication enables valid inference");
  }

  rep0 = rep0[sub];
  uscore0 = uscore0[sub];
  vscore0 = vscore0[sub];
  logRankZ0 = logRankZ0[sub];
  logRankPValue0 = logRankPValue0[sub];


  DataFrame result;

  if (!has_rep) {
    result = DataFrame::create(
      _[rep] = rep0,
      _["uscore"] = uscore0,
      _["vscore"] = vscore0,
      _["logRankZ"] = logRankZ0,
      _["logRankPValue"] = logRankPValue0,
      _["rho1"] = rho1,
      _["rho2"] = rho2);
  } else if (TYPEOF(data[rep]) == INTSXP) {
    IntegerVector rep0n = repwn[rep0-1];

    result = DataFrame::create(
      _[rep] = rep0n,
      _["uscore"] = uscore0,
      _["vscore"] = vscore0,
      _["logRankZ"] = logRankZ0,
      _["logRankPValue"] = logRankPValue0,
      _["rho1"] = rho1,
      _["rho2"] = rho2);
  } else if (TYPEOF(data[rep]) == STRSXP) {
    CharacterVector rep0c = repwc[rep0-1];

    result = DataFrame::create(
      _[rep] = rep0c,
      _["uscore"] = uscore0,
      _["vscore"] = vscore0,
      _["logRankZ"] = logRankZ0,
      _["logRankPValue"] = logRankPValue0,
      _["rho1"] = rho1,
      _["rho2"] = rho2);
  }

  return result;
}

