#include "utilities.h"

using namespace Rcpp;


// define the integrand functions
struct nbparams {
  double tau;
  double phi;
  NumericVector accrualTime;
  NumericVector accrualIntensity;
  NumericVector piecewiseSurvivalTime;
  double kappa;
  double lambda;
  NumericVector zero;
  NumericVector gam;
  double accrualDuration;
};


void f_ex(double *x, int n, void *ex) {
  nbparams *param = (nbparams *) ex;
  NumericVector u0(n);
  for (int i=0; i<n; i++) {
    u0[i] = x[i];
  }
  NumericVector p = patrisk(u0, param->piecewiseSurvivalTime, param->zero,
                            param->gam);
  u0 = param->tau - u0;
  NumericVector N = accrual(u0, param->accrualTime, param->accrualIntensity,
                            param->accrualDuration);
  u0 = (param->phi)*N*p;
  for (int i=0; i<n; i++) {
    x[i] = u0[i];
  }
}


void f_info(double *x, int n, void *ex) {
  nbparams *param = (nbparams *) ex;
  NumericVector u0(n);
  for (int i=0; i<n; i++) {
    u0[i] = x[i];
  }
  NumericVector p = patrisk(u0, param->piecewiseSurvivalTime, param->zero,
                            param->gam);
  NumericVector u = param->lambda/pow(1.0 +
    (param->kappa)*(param->lambda)*u0, 2);

  u0 = param->tau - u0;
  NumericVector N = accrual(u0, param->accrualTime, param->accrualIntensity,
                            param->accrualDuration);

  u0 = (param->phi)*u*N*p;
  for (int i=0; i<n; i++) {
    x[i] = u0[i];
  }
}



//' @title Negative Binomial Rate Ratio by Stratum
//'
//' @description Obtains the number of subjects accrued, number of events,
//' number of dropouts, number of subjects reaching the maximum
//' follow-up, total exposure, rate and variance for log rate in each group,
//' rate ratio and variance for log rate ratio by stratum at a given
//' calendar time.
//'
//' @param time The calendar time for data cut.
//' @param rateRatioH0 Rate ratio under the null hypothesis.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param kappa1 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the active treatment
//'   group by stratum.
//' @param kappa2 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the control group by
//'   stratum.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group by stratum.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group by stratum.
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param nullVariance Whether to calculate the variance for log rate ratio
//'   under the null hypothesis.
//'
//' @return A list with two components:
//'
//' * \code{resultsUnderH1}: A data frame containing the following variables:
//'
//'     - \code{stratum}: The stratum.
//'
//'     - \code{time}: The calendar time since trial start.
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
//'     - \code{rateRatio}: The rate ratio of the active treatment group
//'       versus the control group.
//'
//'     - \code{vlogRate1}: The variance for the log rate parameter for the
//'       active treatment group.
//'
//'     - \code{vlogRate2}: The variance for the log rate parameter for the
//'       control group.
//'
//'     - \code{vlogRR}: The variance of log rate ratio.
//'
//' * \code{resultsUnderH0} when \code{nullVariance = TRUE}: A data frame
//'   with the following variables:
//'
//'     - \code{stratum}: The stratum.
//'
//'     - \code{time}: The analysis time since trial start.
//'
//'     - \code{lambda1H0}: The restricted maximum likelihood estimate
//'       of the event rate for the active treatment group.
//'
//'     - \code{lambda2H0}: The restricted maximum likelihood estimate
//'       of the event rate for the control group.
//'
//'     - \code{rateRatioH0}: The rate ratio under H0.
//'
//'     - \code{vlogRate1H0}: The variance for the log rate parameter for
//'       the active treatment group under H0.
//'
//'     - \code{vlogRate2H0}: The variance for the log rate parameter for
//'       the control group under H0.
//'
//'     - \code{vlogRRH0}: The variance of log rate ratio under H0.
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
//'     - \code{stratum}: The stratum.
//'
//'     - \code{time}: The analysis time since trial start.
//'
//'     - \code{rateRatioH0}: The rate ratio under H0.
//'
//'     - \code{lambda1}: The true event rate for the active treatment group.
//'
//'     - \code{lambda2}: The true event rate for the control group.
//'
//'     - \code{rateRatio}: The true rate ratio.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Example 1: Variable follow-up design
//'
//' nbstat1(time = 2,
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
//' nbstat1(time = 1.8,
//'        accrualIntensity = 220/1.5,
//'        stratumFraction = c(0.2, 0.8),
//'        kappa1 = 3,
//'        kappa2 = 3,
//'        lambda1 = c(0.5*8.4, 0.7*10.2),
//'        lambda2 = c(8.4, 10.2),
//'        gamma1 = 0.05,
//'        gamma2 = 0.05,
//'        accrualDuration = 1.5,
//'        followupTime = 0.5,
//'        fixedFollowup = 1,
//'        nullVariance = 1)
//'
//' @export
// [[Rcpp::export]]
List nbstat1(const double time = NA_REAL,
             const double rateRatioH0 = 1,
             const double allocationRatioPlanned = 1,
             const NumericVector& accrualTime = 0,
             const NumericVector& accrualIntensity = NA_REAL,
             const NumericVector& piecewiseSurvivalTime = 0,
             const NumericVector& stratumFraction = 1,
             const NumericVector& kappa1 = NA_REAL,
             const NumericVector& kappa2 = NA_REAL,
             const NumericVector& lambda1 = NA_REAL,
             const NumericVector& lambda2 = NA_REAL,
             const NumericVector& gamma1 = 0,
             const NumericVector& gamma2 = 0,
             const double accrualDuration = NA_REAL,
             const double followupTime = NA_REAL,
             const bool fixedFollowup = 0,
             const bool nullVariance = 0) {

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;
  NumericVector kappa1x(nstrata), kappa2x(nstrata);
  NumericVector lambda1x(nstrata), lambda2x(nstrata);
  NumericVector gamma1x(nsi), gamma2x(nsi);

  if (kappa1.size() == 1) {
    kappa1x = rep(kappa1, nstrata);
  } else if (kappa1.size() == nstrata) {
    kappa1x = kappa1;
  } else {
    stop("Invalid length for kappa1");
  }

  if (kappa2.size() == 1) {
    kappa2x = rep(kappa2, nstrata);
  } else if (kappa2.size() == nstrata) {
    kappa2x = kappa2;
  } else {
    stop("Invalid length for kappa2");
  }

  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nstrata) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nstrata) {
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

  // number of subjects enrolled
  NumericVector ss(1, time);
  double a = accrual(ss, accrualTime, accrualIntensity, accrualDuration)[0];
  double phi = allocationRatioPlanned/(1.0 + allocationRatioPlanned);

  double frac, k1, k2, lam1, lam2;
  IntegerVector l1 = Range(0, nintervals-1);
  IntegerVector l(nintervals);
  NumericVector zero(nintervals), gam1(nintervals), gam2(nintervals);
  NumericVector stratum(nstrata), times(nstrata), nsubjects(nstrata);
  NumericVector nevents(nstrata), nevents1(nstrata), nevents2(nstrata);
  NumericVector ndropouts(nstrata), ndropouts1(nstrata), ndropouts2(nstrata);
  NumericVector nfmax(nstrata), nfmax1(nstrata), nfmax2(nstrata);
  NumericVector exposure(nstrata), exposure1(nstrata), exposure2(nstrata);
  NumericVector rate1(nstrata), rate2(nstrata), rateRatio(nstrata);
  NumericVector information1(nstrata), information2(nstrata);
  NumericVector vlogRate1(nstrata), vlogRate2(nstrata), vlogRR(nstrata);
  NumericMatrix y(1,2);

  NumericVector maxFU(1, maxFollowupTime);
  NumericVector tt(1, time - maxFollowupTime);
  double a2 = accrual(tt, accrualTime, accrualIntensity, accrualDuration)[0];

  double tol = 1.0e-6;
  double upper = std::min(time, maxFollowupTime);
  NumericVector t = piecewiseSurvivalTime;

  for (int h=0; h<nstrata; h++) {
    stratum[h] = h+1;
    times[h] = time;

    k1 = kappa1x[h];
    k2 = kappa2x[h];
    lam1 = lambda1x[h];
    lam2 = lambda2x[h];

    frac = stratumFraction[h];
    l = h*nintervals + l1;
    gam1 = gamma1x[l];
    gam2 = gamma2x[l];

    // number of enrolled subjects
    nsubjects[h] = frac*a;

    // number of dropouts
    if (max(gam1) + max(gam2) > 0.0) {
      y = nevent2(ss, allocationRatioPlanned, accrualTime,
                  frac*accrualIntensity, piecewiseSurvivalTime,
                  gam1, gam2, zero, zero, accrualDuration,
                  followupTime, maxFollowupTime);
      ndropouts1[h] = y(0,0);
      ndropouts2[h] = y(0,1);
      ndropouts[h] = ndropouts1[h] + ndropouts2[h];
    }

    // number of subjects reaching maximum follow-up
    double ncom = frac*a2;
    double p1 = patrisk(maxFU, piecewiseSurvivalTime, zero, gam1)[0];
    double p2 = patrisk(maxFU, piecewiseSurvivalTime, zero, gam2)[0];
    nfmax1[h] = phi*ncom*p1;
    nfmax2[h] = (1-phi)*ncom*p2;
    nfmax[h] = nfmax1[h] + nfmax2[h];

    // number of events
    nbparams param1 = {time, phi, accrualTime, frac*accrualIntensity,
                       t, k1, lam1, zero, gam1, accrualDuration};
    nbparams param2 = {time, 1-phi, accrualTime, frac*accrualIntensity,
                       t, k2, lam2, zero, gam2, accrualDuration};

    exposure1[h] = quad(f_ex, &param1, 0.0, upper, tol)[0];
    exposure2[h] = quad(f_ex, &param2, 0.0, upper, tol)[0];
    exposure[h] = exposure1[h] + exposure2[h];

    information1[h] = quad(f_info, &param1, 0.0, upper, tol)[0];
    information2[h] = quad(f_info, &param2, 0.0, upper, tol)[0];

    nevents1[h] = lam1*exposure1[h];
    nevents2[h] = lam2*exposure2[h];
    nevents[h] = nevents1[h] + nevents2[h];

    rateRatio[h] = lam1/lam2;
    vlogRate1[h] = 1.0/information1[h];
    vlogRate2[h] = 1.0/information2[h];
    vlogRR[h] = vlogRate1[h] + vlogRate2[h];
  }

  DataFrame resultsUnderH1 = DataFrame::create(
    _["stratum"] = stratum,
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
    _["rateRatio"] = rateRatio,
    _["vlogRate1"] = vlogRate1,
    _["vlogRate2"] = vlogRate2,
    _["vlogRR"] = vlogRR);

  DataFrame resultsUnderH0;
  if (!nullVariance) {
    resultsUnderH0 = DataFrame::create(
      _["stratum"] = stratum,
      _["time"] = time,
      _["rateRatioH0"] = rateRatioH0,
      _["lambda1"] = lambda1x,
      _["lambda2"] = lambda2x,
      _["rateRatio"] = rateRatio);
  } else {
    NumericVector lambda1H0(nstrata), lambda2H0(nstrata);
    NumericVector information1H0(nstrata), information2H0(nstrata);

    auto f = [time, rateRatioH0, phi, accrualTime, &frac, accrualIntensity,
              t, &k1, &k2, &lam1, &lam2, zero, &gam1, &gam2,
              accrualDuration, maxFollowupTime, tol](double aval)->double {

                nbparams param1 = {time, phi, accrualTime,
                                   frac*accrualIntensity,
                                   t, k1, aval*rateRatioH0,
                                   zero, gam1, accrualDuration};

                nbparams param2 = {time, 1-phi, accrualTime,
                                   frac*accrualIntensity,
                                   t, k2, aval,
                                   zero, gam2, accrualDuration};

                double upper = std::min(time, maxFollowupTime);

                double a1 = quad(f_info, &param1, 0.0, upper, tol)[0];
                double a2 = quad(f_info, &param2, 0.0, upper, tol)[0];

                return phi*(lam1/(aval*rateRatioH0) - 1.0)*a1 +
                  (1-phi)*(lam2/aval - 1.0)*a2;
              };

    for (int h=0; h<nstrata; h++) {
      k1 = kappa1x[h];
      k2 = kappa2x[h];
      lam1 = lambda1x[h];
      lam2 = lambda2x[h];

      frac = stratumFraction[h];
      l = h*nintervals + l1;
      gam1 = gamma1x[l];
      gam2 = gamma2x[l];

      double t1 = exposure1[h]/(frac*a*phi);
      double t2 = exposure2[h]/(frac*a*(1-phi));
      double a = (phi*k2 + (1-phi)*k1)*rateRatioH0*t1*t2;
      double b = -(phi*t1*(k2*lam1*t2 - rateRatioH0) +
                   (1-phi)*t2*(k1*lam2*rateRatioH0*t1 - 1));
      double c = -(phi*lam1*t1 + (1-phi)*lam2*t2);

      double init;
      if (k1 == 0 && k2 == 0) {
        init = -c/b;
      } else {
        init = (-b + sqrt(b*b - 4*a*c))/(2*a);
      }

      lambda2H0[h] = brent(f, 0.5*init, 1.5*init, 1.0e-6);
      lambda1H0[h] = rateRatioH0*lambda2H0[h];

      nbparams param1 = {time, phi, accrualTime, frac*accrualIntensity,
                         t, k1, lambda1H0[h], zero, gam1,
                         accrualDuration};
      nbparams param2 = {time, 1-phi, accrualTime, frac*accrualIntensity,
                         t, k2, lambda2H0[h], zero, gam2,
                         accrualDuration};

      information1H0[h] = quad(f_info, &param1, 0.0, upper, tol)[0];
      information2H0[h] = quad(f_info, &param2, 0.0, upper, tol)[0];
    }

    NumericVector vlogRate1H0 = 1.0/information1H0;
    NumericVector vlogRate2H0 = 1.0/information2H0;
    NumericVector vlogRRH0 =  vlogRate1H0 + vlogRate2H0;

    resultsUnderH0 = DataFrame::create(
      _["stratum"] = stratum,
      _["time"] = time,
      _["lambda1H0"] = lambda1H0,
      _["lambda2H0"] = lambda2H0,
      _["rateRatioH0"] = rateRatioH0,
      _["vlogRate1H0"] = vlogRate1H0,
      _["vlogRate2H0"] = vlogRate2H0,
      _["vlogRRH0"] = vlogRRH0,
      _["lambda1"] = lambda1x,
      _["lambda2"] = lambda2x,
      _["rateRatio"] = rateRatio);
  }

  return List::create(
    _["resultsUnderH1"] = resultsUnderH1,
    _["resultsUnderH0"] = resultsUnderH0);
}


//' @title Negative Binomial Rate Ratio
//' @description Obtains the number of subjects accrued, number of events,
//' number of dropouts, number of subjects reaching the maximum
//' follow-up, total exposure, and variance for log rate in each group,
//' rate ratio, variance, and Wald test statistic of
//' log rate ratio at given calendar times.
//'
//' @param time A vector of calendar times for data cut.
//' @param rateRatioH0 Rate ratio under the null hypothesis.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param kappa1 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the active treatment
//'   group by stratum.
//' @param kappa2 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the control group by
//'   stratum.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group by stratum.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group by stratum.
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
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
//' The log-likelihood for \eqn{\{(\kappa_i, \beta_i):i=1,2\}}
//' can be written as
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
//' \deqn{s_1 = n r E\left\{\lambda_1 t_{1j} - \left(\lambda_1t_{1j}
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
//'     - \code{rateRatio}: The rate ratio of the active treatment group
//'       versus the control group.
//'
//'     - \code{vlogRate1}: The variance for the log rate
//'       parameter for the active treatment group.
//'
//'     - \code{vlogRate2}: The variance for the log rate
//'       parameter for the control group.
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
//'     - \code{rateRatioH0}: The rate ratio under H0.
//'
//'     - \code{vlogRate1H0}: The variance for the log rate
//'       parameter for the active treatment group under H0.
//'
//'     - \code{vlogRate2H0}: The variance for the log rate
//'       parameter for the control group under H0.
//'
//'     - \code{vlogRRH0}: The variance of log rate ratio under H0.
//'
//'     - \code{informationH0}: The information of log rate ratio under H0.
//'
//'     - \code{zlogRRH0}: The Z-statistic for log rate ratio with variance
//'       evaluated under H0.
//'
//'     - \code{varianceRatio}: The ratio of the variance under H0 versus
//'       the variance under H1.
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
//'     - \code{rateRatioH0}: The rate ratio under H0.
//'
//'     - \code{varianceRatio}: Equal to 1.
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
//'        stratumFraction = c(0.2, 0.8),
//'        kappa1 = 3,
//'        kappa2 = 3,
//'        lambda1 = c(0.5*8.4, 0.6*10.5),
//'        lambda2 = c(8.4, 10.5),
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
            const NumericVector& stratumFraction = 1,
            const NumericVector& kappa1 = NA_REAL,
            const NumericVector& kappa2 = NA_REAL,
            const NumericVector& lambda1 = NA_REAL,
            const NumericVector& lambda2 = NA_REAL,
            const NumericVector& gamma1 = 0,
            const NumericVector& gamma2 = 0,
            const double accrualDuration = NA_REAL,
            const double followupTime = NA_REAL,
            const bool fixedFollowup = 0,
            const bool nullVariance = 0) {

  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;
  NumericVector kappa1x(nstrata), kappa2x(nstrata);
  NumericVector lambda1x(nstrata), lambda2x(nstrata);
  NumericVector gamma1x(nsi), gamma2x(nsi);

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

  if (is_true(any(is_na(kappa1)))) {
    stop("kappa1 must be provided");
  }

  if (is_true(any(is_na(kappa2)))) {
    stop("kappa2 must be provided");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(kappa1 < 0))) {
    stop("kappa1 must be non-negative");
  }

  if (is_true(any(kappa2 < 0))) {
    stop("kappa2 must be non-negative");
  }

  if (is_true(any(lambda1 <= 0))) {
    stop("lambda1 must be positive");
  }

  if (is_true(any(lambda2 <= 0))) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (kappa1.size() == 1) {
    kappa1x = rep(kappa1, nstrata);
  } else if (kappa1.size() == nstrata) {
    kappa1x = kappa1;
  } else {
    stop("Invalid length for kappa1");
  }

  if (kappa2.size() == 1) {
    kappa2x = rep(kappa2, nstrata);
  } else if (kappa2.size() == nstrata) {
    kappa2x = kappa2;
  } else {
    stop("Invalid length for kappa2");
  }

  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nstrata) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nstrata) {
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


  int k = static_cast<int>(time.size());
  NumericVector subjects(k);
  NumericVector nevents(k), nevents1(k), nevents2(k);
  NumericVector ndropouts(k), ndropouts1(k), ndropouts2(k);
  NumericVector nfmax(k), nfmax1(k), nfmax2(k);
  NumericVector exposure(k), exposure1(k), exposure2(k);
  NumericVector rateRatio(k);
  NumericVector vlogRate1(k), vlogRate2(k), vlogRR(k);
  NumericVector information(k), zlogRR(k);
  NumericVector w = stratumFraction;
  NumericVector lam1(k),lam2(k), varianceRatio(k);
  NumericVector lam1H0(k), lam2H0(k);
  NumericVector vlogRate1H0(k), vlogRate2H0(k), vlogRRH0(k);
  NumericVector informationH0(k), zlogRRH0(k);

  List results;
  DataFrame resultsUnderH1, resultsUnderH0, dfH1, dfH0;
  for (int j=0; j<k; j++) {
    results = nbstat1(
      time[j], rateRatioH0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      kappa1x, kappa2x, lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup, nullVariance);

    dfH1 = DataFrame(results["resultsUnderH1"]);
    dfH0 = DataFrame(results["resultsUnderH0"]);

    subjects[j] = sum(NumericVector(dfH1[2]));
    nevents[j] = sum(NumericVector(dfH1[3]));
    nevents1[j] = sum(NumericVector(dfH1[4]));
    nevents2[j] = sum(NumericVector(dfH1[5]));
    ndropouts[j] = sum(NumericVector(dfH1[6]));
    ndropouts1[j] = sum(NumericVector(dfH1[7]));
    ndropouts2[j] = sum(NumericVector(dfH1[8]));
    nfmax[j] = sum(NumericVector(dfH1[9]));
    nfmax1[j] = sum(NumericVector(dfH1[10]));
    nfmax2[j] = sum(NumericVector(dfH1[11]));
    exposure[j] = sum(NumericVector(dfH1[12]));
    exposure1[j] = sum(NumericVector(dfH1[13]));
    exposure2[j] = sum(NumericVector(dfH1[14]));
    rateRatio[j] = exp(sum(w*log(NumericVector(dfH1[15]))));
    vlogRate1[j] = sum(w*w*NumericVector(dfH1[16]));
    vlogRate2[j] = sum(w*w*NumericVector(dfH1[17]));
    vlogRR[j] = sum(w*w*NumericVector(dfH1[18]));
    information[j] = 1.0/vlogRR[j];
    zlogRR[j] = (log(rateRatio[j]) - log(rateRatioH0))/sqrt(vlogRR[j]);

    if (!nullVariance) {
      lam1[j] = exp(sum(w*log(NumericVector(dfH0[3]))));
      lam2[j] = exp(sum(w*log(NumericVector(dfH0[4]))));
      varianceRatio[j] = 1.0;
    } else {
      lam1H0[j] = exp(sum(w*log(NumericVector(dfH0[2]))));
      lam2H0[j] = exp(sum(w*log(NumericVector(dfH0[3]))));
      vlogRate1H0[j] = sum(w*w*NumericVector(dfH0[5]));
      vlogRate2H0[j] = sum(w*w*NumericVector(dfH0[6]));
      vlogRRH0[j] = sum(w*w*NumericVector(dfH0[7]));
      informationH0[j] = 1.0/vlogRRH0[j];
      zlogRRH0[j] = (log(rateRatio[j]) - log(rateRatioH0))/sqrt(vlogRRH0[j]);
      lam1[j] = exp(sum(w*log(NumericVector(dfH0[8]))));
      lam2[j] = exp(sum(w*log(NumericVector(dfH0[9]))));
      varianceRatio[j] = vlogRRH0[j]/vlogRR[j];
    }
  }

  resultsUnderH1 = DataFrame::create(
    _["time"] = time,
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
    _["exposure"] = exposure,
    _["exposure1"] = exposure1,
    _["exposure2"] = exposure2,
    _["rateRatio"] = rateRatio,
    _["vlogRate1"] = vlogRate1,
    _["vlogRate2"] = vlogRate2,
    _["vlogRR"] = vlogRR,
    _["information"] = information,
    _["zlogRR"] = zlogRR);

  if (!nullVariance) {
    resultsUnderH0 = DataFrame::create(
      _["time"] = time,
      _["rateRatioH0"] = rateRatioH0,
      _["varianceRatio"] = varianceRatio,
      _["lambda1"] = lam1,
      _["lambda2"] = lam2,
      _["rateRatio"] = rateRatio);
  } else {
    resultsUnderH0 = DataFrame::create(
      _["time"] = time,
      _["lambda1H0"] = lam1H0,
      _["lambda2H0"] = lam2H0,
      _["rateRatioH0"] = rateRatioH0,
      _["vlogRate1H0"] = vlogRate1H0,
      _["vlogRate2H0"] = vlogRate2H0,
      _["vlogRRH0"] = vlogRRH0,
      _["informationH0"] = informationH0,
      _["zlogRRH0"] = zlogRRH0,
      _["varianceRatio"] = varianceRatio,
      _["lambda1"] = lam1,
      _["lambda2"] = lam2,
      _["rateRatio"] = rateRatio);
  }

  return List::create(
    _["resultsUnderH1"] = resultsUnderH1,
    _["resultsUnderH0"] = resultsUnderH0);
}


//' @title Power for Negative Binomial Rate Ratio
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
//' @inheritParams param_stratumFraction
//' @param kappa1 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the active treatment
//'   group by stratum.
//' @param kappa2 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the control group by
//'   stratum.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group by stratum.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group by stratum.
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
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
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
//'         stratumFraction = c(0.2, 0.8),
//'         kappa1 = 5, kappa2 = 5,
//'         lambda1 = c(0.7*0.125, 0.75*0.25),
//'         lambda2 = c(0.125, 0.25),
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
             const std::string typeAlphaSpending = "sfOF",
             const double parameterAlphaSpending = NA_REAL,
             const NumericVector& userAlphaSpending = NA_REAL,
             const NumericVector& futilityBounds = NA_REAL,
             const std::string typeBetaSpending = "none",
             const double parameterBetaSpending = NA_REAL,
             const double rateRatioH0 = 1,
             const double allocationRatioPlanned = 1,
             const NumericVector& accrualTime = 0,
             const NumericVector& accrualIntensity = NA_REAL,
             const NumericVector& piecewiseSurvivalTime = 0,
             const NumericVector& stratumFraction = 1,
             const NumericVector& kappa1 = NA_REAL,
             const NumericVector& kappa2 = NA_REAL,
             const NumericVector& lambda1 = NA_REAL,
             const NumericVector& lambda2 = NA_REAL,
             const NumericVector& gamma1 = 0,
             const NumericVector& gamma2 = 0,
             const double accrualDuration = NA_REAL,
             const double followupTime = NA_REAL,
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

  if (is_true(any(is_na(kappa1)))) {
    stop("kappa1 must be provided");
  }

  if (is_true(any(is_na(kappa2)))) {
    stop("kappa2 must be provided");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(kappa1 < 0))) {
    stop("kappa1 must be non-negative");
  }

  if (is_true(any(kappa2 < 0))) {
    stop("kappa2 must be non-negative");
  }

  if (is_true(any(lambda1 <= 0))) {
    stop("lambda1 must be positive");
  }

  if (is_true(any(lambda2 <= 0))) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (kappa1.size() != 1 && kappa1.size() != nstrata) {
    stop("Invalid length for kappa1");
  }

  if (kappa2.size() != 1 && kappa2.size() != nstrata) {
    stop("Invalid length for kappa2");
  }

  if (lambda1.size() != 1 && lambda1.size() != nstrata) {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() != 1 && lambda2.size() != nstrata) {
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
  List na;
  DataFrame nb, nc;
  NumericVector time(kMax), rru(kMax), rrl(kMax);
  NumericVector u0(1, studyDuration1);
  na = nbstat(u0, 1, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup, 0);
  nb = DataFrame(na["resultsUnderH1"]);
  double maxInformation = sum(NumericVector(nb[18]));
  double rateRatio = sum(NumericVector(nb[14]));
  double theta1 = -(log(rateRatio) - log(rateRatioH0));
  NumericVector theta(kMax, theta1);
  NumericVector I = maxInformation*informationRates1;

  double information1;
  auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double t)->double {
              NumericVector u0(1, t);
              List na = nbstat(
                u0, 1, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup, 0);
              DataFrame nb = DataFrame(na["resultsUnderH1"]);
              return sum(NumericVector(nb[18])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    // match the predicted information to the target
    information1 = std::max(I[i], 0.0);
    time[i] = brent(f, 1.0e-6, studyDuration1, 1.0e-6);
  }
  time[kMax-1] = studyDuration1;

  // obtain the variance ratio
  na = nbstat(time, rateRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              nullVariance);
  nb = DataFrame(na["resultsUnderH1"]);
  nc = DataFrame(na["resultsUnderH0"]);
  NumericVector varianceRatio = as<NumericVector>(nc["varianceRatio"]);
  NumericVector w = sqrt(varianceRatio);

  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);
  NumericVector nsubjects = NumericVector(nb[1]);
  NumericVector nsubjects1 = phi*nsubjects;
  NumericVector nsubjects2 = (1-phi)*nsubjects;
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

  rru = rateRatioH0*exp(-criticalValues1/sqrt(I)*w);
  rrl = rateRatioH0*exp(-futilityBounds1/sqrt(I)*w);

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
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["kMax"] = kMax,
    _["rateRatioH0"] = rateRatioH0,
    _["rateRatio"] = rateRatio);

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
    _["kappa1"] = kappa1,
    _["kappa2"] = kappa2,
    _["lambda1"] = lambda1,
    _["lambda2"] = lambda2,
    _["gamma1"] = gamma1,
    _["gamma2"] = gamma2,
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


//' @title Sample Size for Negative Binomial Rate Ratio
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
//' @inheritParams param_stratumFraction
//' @param kappa1 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the active treatment
//'   group by stratum.
//' @param kappa2 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the control group by
//'   stratum.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group by stratum.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group by stratum.
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
//'              accrualIntensity = 100,
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
//'              stratumFraction = c(0.2, 0.8),
//'              kappa1 = 5, kappa2 = 5,
//'              lambda1 = c(0.7*0.125, 0.75*0.25),
//'              lambda2 = c(0.125, 0.25),
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
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const NumericVector& futilityBounds = NA_REAL,
                  const std::string typeBetaSpending = "none",
                  const double parameterBetaSpending = NA_REAL,
                  const NumericVector& userBetaSpending = NA_REAL,
                  const double rateRatioH0 = 1,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = NA_REAL,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const NumericVector& stratumFraction = 1,
                  const NumericVector& kappa1 = NA_REAL,
                  const NumericVector& kappa2 = NA_REAL,
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
  NumericVector lambda1x(nstrata), lambda2x(nstrata);

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

  if (is_true(any(is_na(kappa1)))) {
    stop("kappa1 must be provided");
  }

  if (is_true(any(is_na(kappa2)))) {
    stop("kappa2 must be provided");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(kappa1 < 0))) {
    stop("kappa1 must be non-negative");
  }

  if (is_true(any(kappa2 < 0))) {
    stop("kappa2 must be non-negative");
  }

  if (is_true(any(lambda1 <= 0))) {
    stop("lambda1 must be positive");
  }

  if (is_true(any(lambda2 <= 0))) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (kappa1.size() != 1 && kappa1.size() != nstrata) {
    stop("Invalid length for kappa1");
  }

  if (kappa2.size() != 1 && kappa2.size() != nstrata) {
    stop("Invalid length for kappa2");
  }

  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nstrata) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nstrata) {
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


  double rateRatio = exp(sum(stratumFraction*log(lambda1x/lambda2x)));
  double theta1 = -(log(rateRatio) - log(rateRatioH0));
  NumericVector theta(kMax, theta1);

  double maxInformation;
  if (!nullVariance) {
    List design = getDesign(
      beta, NA_REAL, theta1, kMax, informationRates1,
      efficacyStopping1, futilityStopping1, criticalValues1,
      alpha1, asf, asfpar, userAlphaSpending, futilityBounds1,
      bsf, bsfpar, userBetaSpending, spendingTime1, 1);

    DataFrame byStageResults = DataFrame(design["byStageResults"]);
    futilityBounds1 = byStageResults["futilityBounds"];

    DataFrame overallResults = DataFrame(design["overallResults"]);
    maxInformation = overallResults["information"];

    auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
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
                  piecewiseSurvivalTime, stratumFraction,
                  kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                  dur1, dur2, fixedFollowup, 0);
                DataFrame nb = DataFrame(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
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
              rateRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              spendingTime1, nullVariance, theta, unknown,
              missingFutilityBounds](double aval)-> double{
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

                double studyDuration1 = dur1 + dur2;

                // obtain the timing of interim analysis
                NumericVector u0(1, studyDuration1), time(kMax);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                  dur1, dur2, fixedFollowup, 0);
                DataFrame nb = DataFrame(na["resultsUnderH1"]);
                double maxInformation = sum(NumericVector(nb[18]));
                NumericVector I = maxInformation*informationRates1;
                double information1;

                auto g = [allocationRatioPlanned,
                          accrualTime, accrualIntensity1,
                          piecewiseSurvivalTime, stratumFraction,
                          kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                          dur1, dur2, fixedFollowup,
                          &information1](double aval)->double {
                            NumericVector u0(1, aval);
                            List na = nbstat(
                              u0, 1, allocationRatioPlanned,
                              accrualTime, accrualIntensity1,
                              piecewiseSurvivalTime, stratumFraction,
                              kappa1, kappa2, lambda1, lambda2,
                              gamma1, gamma2,
                              dur1, dur2, fixedFollowup, 0);

                            DataFrame nb = DataFrame(na["resultsUnderH1"]);
                            return sum(NumericVector(nb[18])) - information1;
                          };

                for (int i=0; i<kMax-1; i++) {
                  // match the predicted information to the target
                  information1 = std::max(I[i], 0.0);
                  time[i] = brent(g, 1.0e-6, studyDuration1, 1.0e-6);
                }
                time[kMax-1] = studyDuration1;

                // obtain the variance ratio
                na = nbstat(time, rateRatioH0, allocationRatioPlanned,
                            accrualTime, accrualIntensity1,
                            piecewiseSurvivalTime, stratumFraction,
                            kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                            dur1, dur2, fixedFollowup, nullVariance);
                DataFrame nc = DataFrame(na["resultsUnderH0"]);
                NumericVector varianceRatio =
                  as<NumericVector>(nc["varianceRatio"]);
                NumericVector w = sqrt(varianceRatio);

                NumericVector st = spendingTime1;

                // compute the stagewise exit probabilities
                List probs;
                if (!missingFutilityBounds || bsf=="none" || kMax==1) {
                  probs = exitprobcpp(criticalValues1*w, futilityBounds1*w,
                                      theta, I);
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
                    epsilon = R::pnorm(criticalValues1[0]*w[0] -
                      theta[0]*sqrt(I[0]), 0, 1, 1, 0) - cumBetaSpent;
                    if (epsilon < 0) return -1.0;
                    futilityBounds1[0] = (R::qnorm(cumBetaSpent, 0, 1, 1, 0)
                                            + theta[0]*sqrt(I[0]))/w[0];
                  }

                  // lambda expression for finding futility bound at stage k
                  auto g = [&k, &cumBetaSpent, criticalValues1,
                            &futilityBounds1, theta, I,
                            w](double aval)->double {
                              NumericVector u(k+1), l(k+1);
                              for (int i=0; i<k; i++) {
                                u[i] = criticalValues1[i]*w[i];
                                l[i] = futilityBounds1[i]*w[i];
                              }
                              u[k] = 6.0;
                              l[k] = aval*w[k];

                              IntegerVector idx = Range(0,k);
                              List probs = exitprobcpp(u, l, theta[idx],
                                                       I[idx]);
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

    NumericVector u0(1, accrualDuration + followupTime);
    List na = nbstat(u0, 1, allocationRatioPlanned,
                     accrualTime, accrualIntensity1,
                     piecewiseSurvivalTime, stratumFraction,
                     kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                     accrualDuration, followupTime, fixedFollowup, 0);
    DataFrame nb = DataFrame(na["resultsUnderH1"]);
    maxInformation = sum(NumericVector(nb[18]));
  }

  futilityBounds1[kMax-1] = criticalValues1[kMax-1];


  // output the results
  List resultH1, resultH0, result;
  double studyDuration = accrualDuration + followupTime;

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
        auto h = [allocationRatioPlanned, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    List na = nbstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval, fixedFollowup, 0);
                    DataFrame nb = DataFrame(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[18])) - maxInformation;
                  };

        double lower = 0.0, upper = 1.1*followupTime;
        while (h(upper) < 0) {
          lower = upper;
          upper = 2.0*upper;
        }
        followupTime = brent(h, lower, upper, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // adjust study duration to obtain the target maximum information
        auto h = [allocationRatioPlanned, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    List na = nbstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup, 0);
                    DataFrame nb = DataFrame(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[18])) - maxInformation;
                  };

        double aval = brent(h, 0.0, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
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
    piecewiseSurvivalTime, stratumFraction,
    kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration, nullVariance);


  // obtain results under H0 by matching the maximum information
  if (!fixedFollowup) {
    auto h = [rateRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              kappa1, kappa2, lambda2, gamma1, gamma2,
              accrualDuration, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa1, kappa2, lambda2*rateRatioH0, lambda2,
                  gamma1, gamma2,
                  accrualDuration, aval, fixedFollowup, 0);
                DataFrame nb = DataFrame(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
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
      auto g = [rateRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                kappa1, kappa2, lambda2, gamma1, gamma2,
                fixedFollowup, maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  List na = nbstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    kappa1, kappa2, lambda2*rateRatioH0, lambda2,
                    gamma1, gamma2, aval, 0, fixedFollowup, 0);
                  DataFrame nb = DataFrame(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[18])) - maxInformation;
                };

      accrualDuration = brent(g, 1.0e-6, accrualDuration, 1.0e-6);
      followupTime = 0.0;
      studyDuration = accrualDuration + followupTime;
    }
  } else { // fixed follow-up
    auto h = [rateRatioH0, allocationRatioPlanned,
              accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              kappa1, kappa2, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                List na = nbstat(
                  u0, 1, allocationRatioPlanned,
                  accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa1, kappa2, lambda2*rateRatioH0, lambda2,
                  gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup, 0);
                DataFrame nb = DataFrame(na["resultsUnderH1"]);
                return sum(NumericVector(nb[18])) - maxInformation;
              };

    if (h(followupTime) < 0) { // increase accrual duration
      auto g = [rateRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                kappa1, kappa2, lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval + followupTime);
                  List na = nbstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    kappa1, kappa2, lambda2*rateRatioH0, lambda2,
                    gamma1, gamma2,
                    aval, followupTime, fixedFollowup, 0);
                  DataFrame nb = DataFrame(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[18])) - maxInformation;
                };

      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (g(upper) < 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(g, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (h(0) < 0) { // shorten study duration
      double aval = brent(h, 0.0, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    } else { // decrease accrual duration
      auto g = [rateRatioH0, allocationRatioPlanned,
                accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                kappa1, kappa2, lambda2, gamma1, gamma2,
                followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  List na = nbstat(
                    u0, 1, allocationRatioPlanned,
                    accrualTime, accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    kappa1, kappa2, lambda2*rateRatioH0, lambda2,
                    gamma1, gamma2,
                    aval, followupTime, fixedFollowup, 0);
                  DataFrame nb = DataFrame(na["resultsUnderH1"]);
                  return sum(NumericVector(nb[18])) - maxInformation;
                };

      accrualDuration = brent(g, 1.0e-6, accrualDuration, 1.0e-6);
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
    piecewiseSurvivalTime, stratumFraction,
    kappa1, kappa2, lambda2*rateRatioH0, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration, 0);

  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}


//' @title Power for One-Sample Negative Binomial Rate
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
//' @param lambdaH0 The rate parameter of the negative binomial distribution
//'   under the null hypothesis.
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param kappa The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) of the negative binomial
//'   distribution by stratum.
//' @param lambda The rate parameter of the negative binomial distribution
//'   under the alternative hypothesis by stratum.
//' @param gamma The hazard rate for exponential dropout or a vector of
//'   hazard rates for piecewise exponential dropout by stratum.
//'   Defaults to 0 for no dropout.
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
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{lambdaH0}: The rate parameter of the negative binomial
//'       distribution under the null hypothesis.
//'
//'     - \code{lambda}: The overall rate parameter of the negative binomial
//'       distribution under the alternative hypothesis.
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
//'   \code{stratumFraction}, \code{kappa}, \code{lambda}, \code{gamma},
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
//'           lambdaH0 = 0.125, accrualIntensity = 500,
//'           stratumFraction = c(0.2, 0.8),
//'           kappa = c(3, 5), lambda = c(0.0875, 0.085),
//'           gamma = 0, accrualDuration = 1.25,
//'           followupTime = 2.75, fixedFollowup = FALSE)
//'
//' # Example 2: Fixed follow-up design
//'
//' nbpower1s(kMax = 2, informationRates = c(0.5, 1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           lambdaH0 = 8.4, accrualIntensity = 40,
//'           kappa = 3, lambda = 0.5*8.4,
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
               const std::string typeAlphaSpending = "sfOF",
               const double parameterAlphaSpending = NA_REAL,
               const NumericVector& userAlphaSpending = NA_REAL,
               const NumericVector& futilityBounds = NA_REAL,
               const std::string typeBetaSpending = "none",
               const double parameterBetaSpending = NA_REAL,
               const double lambdaH0 = NA_REAL,
               const NumericVector& accrualTime = 0,
               const NumericVector& accrualIntensity = NA_REAL,
               const NumericVector& piecewiseSurvivalTime = 0,
               const NumericVector& stratumFraction = 1,
               const NumericVector& kappa = NA_REAL,
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

  if (std::isnan(lambdaH0)) {
    stop("lambdaH0 must be provided");
  }

  if (lambdaH0 <= 0) {
    stop("lambdaH0 must be positive");
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

  if (is_true(any(is_na(kappa)))) {
    stop("kappa must be provided");
  }

  if (is_true(any(is_na(lambda)))) {
    stop("lambda must be provided");
  }

  if (is_true(any(kappa < 0))) {
    stop("kappa must be non-negative");
  }

  if (is_true(any(lambda <= 0))) {
    stop("lambda must be positive");
  }

  if (is_true(any(gamma < 0))) {
    stop("gamma must be non-negative");
  }

  if (kappa.size() != 1 && kappa.size() != nstrata) {
    stop("Invalid length for kappa");
  }

  if (lambda.size() != 1 && lambda.size() != nstrata) {
    stop("Invalid length for lambda");
  }

  if (gamma.size() != 1 && gamma.size() != nintervals &&
      gamma.size() != nsi) {
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
  List na;
  DataFrame nb, nc;
  NumericVector time(kMax), rateu(kMax), ratel(kMax);
  NumericVector u0(1, studyDuration1);
  na = nbstat(u0, 1, 1, accrualTime, 2.0*accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              kappa, kappa, lambda, lambda, gamma, gamma,
              accrualDuration, followupTime, fixedFollowup, 0);
  nb = DataFrame(na["resultsUnderH1"]);
  nc = DataFrame(na["resultsUnderH0"]);
  double maxInformation = 2.0*sum(NumericVector(nb[18]));
  double rate = sum(NumericVector(nc[3]));
  double theta1 = -(log(rate) - log(lambdaH0));
  NumericVector theta(kMax, theta1);
  NumericVector I = maxInformation*informationRates1;

  double information1;
  auto f = [accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            kappa, lambda, gamma,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double t)->double {
              NumericVector u0(1, t);
              List na = nbstat(
                u0, 1, 1, accrualTime, 2.0*accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                kappa, kappa, lambda, lambda, gamma, gamma,
                accrualDuration, followupTime, fixedFollowup, 0);
              DataFrame nb = DataFrame(na["resultsUnderH1"]);
              return 2.0*sum(NumericVector(nb[18])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    // match the predicted information to the target
    information1 = std::max(I[i], 0.0);
    time[i] = brent(f, 1.0e-6, studyDuration1, 1.0e-6);
  }
  time[kMax-1] = studyDuration1;

  na = nbstat(time, 1, 1, accrualTime, 2.0*accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              kappa, kappa, lambda, lambda, gamma, gamma,
              accrualDuration, followupTime, fixedFollowup, 0);
  nb = DataFrame(na["resultsUnderH1"]);

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

  rateu = lambdaH0*exp(-criticalValues1/sqrt(I));
  ratel = lambdaH0*exp(-futilityBounds1/sqrt(I));

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
    _["accrualDuration"] = accrualDuration,
    _["followupTime"] = followupTime,
    _["fixedFollowup"] = fixedFollowup,
    _["kMax"] = kMax,
    _["lambdaH0"] = lambdaH0,
    _["lambda"] = rate);

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
    _["kappa"] = kappa,
    _["lambda"] = lambda,
    _["gamma"] = gamma,
    _["spendingTime"] = spendingTime);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "nbpower1s";

  return result;
}


//' @title Sample Size for One-Sample Negative Binomial Rate
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
//' @param lambdaH0 The rate parameter of the negative binomial distribution
//'   under the null hypothesis.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param kappa The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) of the negative binomial
//'   distribution by stratum.
//' @param lambda The rate parameter of the negative binomial distribution
//'   under the alternative hypothesis by stratum.
//' @param gamma The hazard rate for exponential dropout or a vector of
//'   hazard rates for piecewise exponential dropout by stratum.
//'   Defaults to 0 for no dropout.
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
//'                lambdaH0 = 0.125, accrualIntensity = 500,
//'                stratumFraction = c(0.2, 0.8),
//'                kappa = c(3, 5), lambda = c(0.0875, 0.085),
//'                gamma = 0, accrualDuration = 1.25,
//'                followupTime = NA, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains accrual intensity given power, accrual duration, and
//' # follow-up duration for variable follow-up
//'
//' nbsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.5, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                lambdaH0 = 0.125, accrualIntensity = 100,
//'                kappa = 5, lambda = 0.0875,
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
//'                lambdaH0 = 8.4, accrualIntensity = 40,
//'                kappa = 3, lambda = 4.2,
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
                    const std::string typeAlphaSpending = "sfOF",
                    const double parameterAlphaSpending = NA_REAL,
                    const NumericVector& userAlphaSpending = NA_REAL,
                    const NumericVector& futilityBounds = NA_REAL,
                    const std::string typeBetaSpending = "none",
                    const double parameterBetaSpending = NA_REAL,
                    const NumericVector& userBetaSpending = NA_REAL,
                    const double lambdaH0 = NA_REAL,
                    const NumericVector& accrualTime = 0,
                    const NumericVector& accrualIntensity = NA_REAL,
                    const NumericVector& piecewiseSurvivalTime = 0,
                    const NumericVector& stratumFraction = 1,
                    const NumericVector& kappa = NA_REAL,
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
  NumericVector lambdax(nstrata);


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

  if (std::isnan(lambdaH0)) {
    stop("lambdaH0 must be provided");
  }

  if (lambdaH0 <= 0) {
    stop("lambdaH0 must be positive");
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

  if (is_true(any(is_na(kappa)))) {
    stop("kappa must be provided");
  }

  if (is_true(any(is_na(lambda)))) {
    stop("lambda must be provided");
  }

  if (is_true(any(kappa < 0))) {
    stop("kappa must be non-negative");
  }

  if (is_true(any(lambda <= 0))) {
    stop("lambda must be positive");
  }

  if (is_true(any(gamma < 0))) {
    stop("gamma must be non-negative");
  }

  if (kappa.size() != 1 && kappa.size() != nstrata) {
    stop("Invalid length for kappa");
  }

  if (lambda.size() == 1) {
    lambdax = rep(lambda, nstrata);
  } else if (lambda.size() == nstrata) {
    lambdax = lambda;
  } else {
    stop("Invalid length for lambda");
  }

  if (gamma.size() != 1 && gamma.size() != nintervals &&
      gamma.size() != nsi) {
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


  double rate = exp(sum(stratumFraction*log(lambdax)));
  double theta1 = -(log(rate) - log(lambdaH0));
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

  auto f = [accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            kappa, lambda, gamma,
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
                u0, 1, 1, accrualTime, 2.0*accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                kappa, kappa, lambda, lambda, gamma, gamma,
                dur1, dur2, fixedFollowup, 0);
              DataFrame nb = DataFrame(na["resultsUnderH1"]);
              return 2.0*sum(NumericVector(nb[18])) - maxInformation;
            };

  if (unknown == "accrualDuration") {
    accrualDuration = brent(f, interval[0], interval[1], 1.0e-6);
  } else if (unknown == "followupTime") {
    followupTime = brent(f, interval[0], interval[1], 1.0e-6);
  } else if (unknown == "accrualIntensity") {
    double aval = brent(f, interval[0], interval[1], 1.0e-6);
    accrualIntensity1 = aval*accrualIntensity;
  }


  // output the results
  List resultH1, resultH0, result;
  double studyDuration = accrualDuration + followupTime;

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
        auto h = [accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa, lambda, gamma,
                  accrualDuration, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    List na = nbstat(
                      u0, 1, 1, accrualTime, 2.0*accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      kappa, kappa, lambda, lambda, gamma, gamma,
                      accrualDuration, aval, fixedFollowup, 0);
                    DataFrame nb = DataFrame(na["resultsUnderH1"]);
                    return 2.0*sum(NumericVector(nb[18])) - maxInformation;
                  };

        double lower = 0.0, upper = 1.1*followupTime;
        while (h(upper) < 0) {
          lower = upper;
          upper = 2.0*upper;
        }
        followupTime = brent(h, lower, upper, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // adjust study duration to obtain the target maximum information
        auto h = [accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa, lambda, gamma,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    List na = nbstat(
                      u0, 1, 1, accrualTime, 2.0*accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      kappa, kappa, lambda, lambda, gamma, gamma,
                      accrualDuration, followupTime, fixedFollowup, 0);
                    DataFrame nb = DataFrame(na["resultsUnderH1"]);
                    return 2.0*sum(NumericVector(nb[18])) - maxInformation;
                  };

        double aval = brent(h, 0.0, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }

  resultH1 = nbpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending, lambdaH0,
    accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    kappa, lambda, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);


  // obtain results under H0 by matching the maximum information
  NumericVector lamH0(nstrata, lambdaH0);

  if (!fixedFollowup) {
    auto h = [lamH0, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              kappa, gamma,
              accrualDuration, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                List na = nbstat(
                  u0, 1, 1, accrualTime, 2.0*accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa, kappa, lamH0, lamH0, gamma, gamma,
                  accrualDuration, aval, fixedFollowup, 0);
                DataFrame nb = DataFrame(na["resultsUnderH1"]);
                return 2.0*sum(NumericVector(nb[18])) - maxInformation;
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
      auto g = [lamH0, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                kappa, gamma, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  List na = nbstat(
                    u0, 1, 1, accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    kappa, kappa, lamH0, lamH0, gamma, gamma,
                    aval, 0, fixedFollowup, 0);
                  DataFrame nb = DataFrame(na["resultsUnderH1"]);
                  return 2.0*sum(NumericVector(nb[18])) - maxInformation;
                };

      accrualDuration = brent(g, 1.0e-6, accrualDuration, 1.0e-6);
      followupTime = 0.0;
      studyDuration = accrualDuration + followupTime;
    }
  } else { // fixed follow-up
    auto h = [lamH0, accrualTime, accrualIntensity1,
              piecewiseSurvivalTime, stratumFraction,
              kappa, gamma,
              accrualDuration, followupTime, fixedFollowup,
              maxInformation](double aval)->double {
                NumericVector u0(1, accrualDuration + aval);
                List na = nbstat(
                  u0, 1, 1, accrualTime, 2.0*accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa, kappa, lamH0, lamH0, gamma, gamma,
                  accrualDuration, followupTime, fixedFollowup, 0);
                DataFrame nb = DataFrame(na["resultsUnderH1"]);
                return 2.0*sum(NumericVector(nb[18])) - maxInformation;
              };

    if (h(followupTime) < 0) { // increase the accrual duration
      auto g = [lamH0, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                kappa, gamma, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval + followupTime);
                  List na = nbstat(
                    u0, 1, 1, accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    kappa, kappa, lamH0, lamH0, gamma, gamma,
                    aval, followupTime, fixedFollowup, 0);
                  DataFrame nb = DataFrame(na["resultsUnderH1"]);
                  return 2.0*sum(NumericVector(nb[18])) - maxInformation;
                };

      double lower = accrualDuration, upper = 2.0*accrualDuration;
      while (g(upper) < 0) {
        lower = upper;
        upper = 2.0*upper;
      }
      accrualDuration = brent(g, lower, upper, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (h(0) < 0) { // shorten the study duration
      double aval = brent(h, 0.0, followupTime, 1.0e-6);
      studyDuration = accrualDuration + aval;
    } else { // decrease the accrual duration
      auto g = [lamH0, accrualTime, accrualIntensity1,
                piecewiseSurvivalTime, stratumFraction,
                kappa, gamma, followupTime, fixedFollowup,
                maxInformation](double aval)->double {
                  NumericVector u0(1, aval);
                  List na = nbstat(
                    u0, 1, 1, accrualTime, 2.0*accrualIntensity1,
                    piecewiseSurvivalTime, stratumFraction,
                    kappa, kappa, lamH0, lamH0, gamma, gamma,
                    aval, followupTime, fixedFollowup, 0);
                  DataFrame nb = DataFrame(na["resultsUnderH1"]);
                  return 2.0*sum(NumericVector(nb[18])) - maxInformation;
                };

      accrualDuration = brent(g, 1.0e-6, accrualDuration, 1.0e-6);
      studyDuration = accrualDuration;
    }
  }


  // use the same stopping boundaries as under H1
  resultH0 = nbpower1s(
    kMax, informationRates1,
    efficacyStopping1, futilityStopping1, criticalValues1,
    alpha1, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, futilityBounds1,
    typeBetaSpending, parameterBetaSpending, lambdaH0,
    accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    kappa, lamH0, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  result = List::create(
    _["resultsUnderH1"] = resultH1,
    _["resultsUnderH0"] = resultH0);

  return result;
}


//' @title Power for Equivalence in Negative Binomial Rate Ratio
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
//' @inheritParams param_stratumFraction
//' @param kappa1 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the active treatment group
//'   by stratum.
//' @param kappa2 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the control group by stratum.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group by stratum.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group by stratum.
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
//' @return An S3 class \code{nbpowerequiv} object with 4 components:
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
//'     - \code{exposure}: The total exposure.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
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
//'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
//'   \code{kappa1}, \code{kappa2},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   \code{spendingTime}.
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
//'              followupTime = 2.75, fixedFollowup = FALSE)
//'
//' # Example 2: Fixed follow-up design
//' nbpowerequiv(kMax = 2, informationRates = c(0.5, 1),
//'              alpha = 0.05, typeAlphaSpending = "sfOF",
//'              rateRatioLower = 0.5, rateRatioUpper = 2,
//'              accrualIntensity = 220/1.5,
//'              stratumFraction = c(0.2, 0.8),
//'              kappa1 = 3, kappa2 = 3,
//'              lambda1 = c(8.4, 10.2),
//'              lambda2 = c(8.0, 11.5),
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
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& userAlphaSpending = NA_REAL,
                  const double rateRatioLower = NA_REAL,
                  const double rateRatioUpper = NA_REAL,
                  const double allocationRatioPlanned = 1,
                  const NumericVector& accrualTime = 0,
                  const NumericVector& accrualIntensity = NA_REAL,
                  const NumericVector& piecewiseSurvivalTime = 0,
                  const NumericVector& stratumFraction = 1,
                  const NumericVector& kappa1 = NA_REAL,
                  const NumericVector& kappa2 = NA_REAL,
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

  if (std::isnan(rateRatioLower)) {
    stop("rateRatioLower must be provided");
  }

  if (std::isnan(rateRatioUpper)) {
    stop("rateRatioUpper must be provided");
  }

  if (rateRatioLower <= 0) {
    stop("rateRatioLower must be positive");
  }

  if (rateRatioLower >= rateRatioUpper) {
    stop("rateRatio lower must be less than rateRatioUpper");
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

  if (is_true(any(is_na(kappa1)))) {
    stop("kappa1 must be provided");
  }

  if (is_true(any(is_na(kappa2)))) {
    stop("kappa2 must be provided");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(kappa1 < 0))) {
    stop("kappa1 must be non-negative");
  }

  if (is_true(any(kappa2 < 0))) {
    stop("kappa2 must be non-negative");
  }

  if (is_true(any(lambda1 <= 0))) {
    stop("lambda1 must be positive");
  }

  if (is_true(any(lambda2 <= 0))) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (kappa1.size() != 1 && kappa1.size() != nstrata) {
    stop("Invalid length for kappa1");
  }

  if (kappa2.size() != 1 && kappa2.size() != nstrata) {
    stop("Invalid length for kappa2");
  }

  if (lambda1.size() != 1 && lambda1.size() != nstrata) {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() != 1 && lambda2.size() != nstrata) {
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
  List na;
  DataFrame nb, nc;
  NumericVector time(kMax);
  NumericVector u0(1, studyDuration1);
  na = nbstat(u0, 1, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup, 0);
  nb = DataFrame(na["resultsUnderH1"]);

  double maxInformation = sum(NumericVector(nb[18]));
  double rateRatio = sum(NumericVector(nb[14]));
  double theta1 = log(rateRatio);
  NumericVector theta(kMax, theta1);
  NumericVector I = maxInformation*informationRates1;

  double information1;
  auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
            accrualDuration, followupTime, fixedFollowup,
            &information1](double t)->double {
              NumericVector u0(1, t);
              List na = nbstat(
                u0, 1, allocationRatioPlanned,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction,
                kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                accrualDuration, followupTime, fixedFollowup, 0);
              DataFrame nb = DataFrame(na["resultsUnderH1"]);
              return sum(NumericVector(nb[18])) - information1;
            };

  for (int i=0; i<kMax-1; i++) {
    // match the predicted information to the target
    information1 = std::max(I[i], 0.0);
    time[i] = brent(f, 1.0e-6, studyDuration1, 1.0e-6);
  }
  time[kMax-1] = studyDuration1;

  na = nbstat(time, 1, allocationRatioPlanned,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
              accrualDuration, followupTime, fixedFollowup, 0);
  nb = DataFrame(na["resultsUnderH1"]);

  double phi = allocationRatioPlanned/(allocationRatioPlanned+1);
  NumericVector nsubjects = NumericVector(nb[1]);
  NumericVector nsubjects1 = phi*nsubjects;
  NumericVector nsubjects2 = (1-phi)*nsubjects;
  NumericVector nevents = NumericVector(nb[2]);
  NumericVector nevents1 = NumericVector(nb[3]);
  NumericVector nevents2 = NumericVector(nb[4]);
  NumericVector ndropouts = NumericVector(nb[5]);
  NumericVector ndropouts1 = NumericVector(nb[6]);
  NumericVector ndropouts2 = NumericVector(nb[7]);
  NumericVector exposure = NumericVector(nb[11]);
  NumericVector exposure1 = NumericVector(nb[12]);
  NumericVector exposure2 = NumericVector(nb[13]);

  // calculate cumulative rejection probability under H1
  NumericVector theta10 = rep(log(rateRatioLower), kMax);
  NumericVector theta20 = rep(log(rateRatioUpper), kMax);
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

  NumericVector efficacyRateRatioLower = exp(theta10 + b/sqrt(I));
  NumericVector efficacyRateRatioUpper = exp(theta20 - b/sqrt(I));

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
  double expectedExposure = sum(q*exposure);
  double expectedNumberOfEvents1 = sum(q*nevents1);
  double expectedNumberOfDropouts1 = sum(q*ndropouts1);
  double expectedNumberOfSubjects1 = sum(q*nsubjects1);
  double expectedExposure1 = sum(q*exposure1);
  double expectedNumberOfEvents2 = sum(q*nevents2);
  double expectedNumberOfDropouts2 = sum(q*ndropouts2);
  double expectedNumberOfSubjects2 = sum(q*nsubjects2);
  double expectedExposure2 = sum(q*exposure2);
  double expectedStudyDuration = sum(q*time);
  double expectedInformation = sum(q*I);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = alpha,
    _["numberOfEvents"] = (nevents[kMax-1]),
    _["numberOfSubjects"] = (nsubjects[kMax-1]),
    _["exposure"] = (exposure[kMax-1]),
    _["studyDuration"] = (time[kMax-1]),
    _["information"] = maxInformation,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedExposure"] = expectedExposure,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["expectedInformation"] = expectedInformation,
    _["kMax"] = kMax,
    _["rateRatioLower"] = rateRatioLower,
    _["rateRatioUpper"] = rateRatioUpper,
    _["rateRatio"] = rateRatio,
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
    _["exposure"] = exposure,
    _["analysisTime"] = time,
    _["efficacyRateRatioLower"] = efficacyRateRatioLower,
    _["efficacyRateRatioUpper"] = efficacyRateRatioUpper,
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
    _["kappa1"] = kappa1,
    _["kappa2"] = kappa2,
    _["lambda1"] = lambda1,
    _["lambda2"] = lambda2,
    _["gamma1"] = gamma1,
    _["gamma2"] = gamma2,
    _["spendingTime"] = spendingTime);

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


//' @title Sample Size for Equivalence in Negative Binomial Rate Ratio
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
//' @inheritParams param_stratumFraction
//' @param kappa1 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the active treatment group by
//'   stratum.
//' @param kappa2 The dispersion parameter (reciprocal of the shape parameter
//'   of the gamma mixing distribution) for the control group by stratum.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group by stratum.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group by stratum.
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
//'                   stratumFraction = c(0.2, 0.8),
//'                   kappa1 = c(3, 5),
//'                   kappa2 = c(2, 3),
//'                   lambda1 = c(0.125, 0.165),
//'                   lambda2 = c(0.135, 0.175),
//'                   gamma1 = -log(1-0.05),
//'                   gamma2 = -log(1-0.10),
//'                   accrualDuration = 1.25,
//'                   followupTime = NA, fixedFollowup = FALSE)
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
                       const std::string typeAlphaSpending = "sfOF",
                       const double parameterAlphaSpending = NA_REAL,
                       const NumericVector& userAlphaSpending = NA_REAL,
                       const double rateRatioLower = NA_REAL,
                       const double rateRatioUpper = NA_REAL,
                       const double allocationRatioPlanned = 1,
                       const NumericVector& accrualTime = 0,
                       const NumericVector& accrualIntensity = NA_REAL,
                       const NumericVector& piecewiseSurvivalTime = 0,
                       const NumericVector& stratumFraction = 1,
                       const NumericVector& kappa1 = NA_REAL,
                       const NumericVector& kappa2 = NA_REAL,
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
  NumericVector lambda1x(nstrata), lambda2x(nstrata);


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

  if (std::isnan(rateRatioLower)) {
    stop("rateRatioLower must be provided");
  }

  if (std::isnan(rateRatioUpper)) {
    stop("rateRatioUpper must be provided");
  }

  if (rateRatioLower <= 0) {
    stop("rateRatioLower must be positive");
  }

  if (rateRatioLower >= rateRatioUpper) {
    stop("rateRatio lower must be less than rateRatioUpper");
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

  if (is_true(any(is_na(kappa1)))) {
    stop("kappa1 must be provided");
  }

  if (is_true(any(is_na(kappa2)))) {
    stop("kappa2 must be provided");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(kappa1 < 0))) {
    stop("kappa1 must be non-negative");
  }

  if (is_true(any(kappa2 < 0))) {
    stop("kappa2 must be non-negative");
  }

  if (is_true(any(lambda1 <= 0))) {
    stop("lambda1 must be positive");
  }

  if (is_true(any(lambda2 <= 0))) {
    stop("lambda2 must be positive");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (kappa1.size() != 1 && kappa1.size() != nstrata) {
    stop("Invalid length for kappa1");
  }

  if (kappa2.size() != 1 && kappa2.size() != nstrata) {
    stop("Invalid length for kappa2");
  }

  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nstrata) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nstrata) {
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
  double theta10 = log(rateRatioLower), theta20 = log(rateRatioUpper);
  double rateRatio = exp(sum(stratumFraction*log(lambda1x/lambda2x)));
  double theta1 = log(rateRatio);
  NumericVector theta(kMax, theta1);
  double maxInformation;

  List design = getDesignEquiv(
    beta, NA_REAL, theta10, theta20, theta1,
    kMax, informationRates1, criticalValues1,
    alpha, asf, asfpar, userAlphaSpending, spendingTime1);

  DataFrame overallResults = DataFrame(design["overallResults"]);
  maxInformation = overallResults["information"];

  auto f = [allocationRatioPlanned, accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
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
                piecewiseSurvivalTime, stratumFraction,
                kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                dur1, dur2, fixedFollowup, 0);
              DataFrame nb = DataFrame(na["resultsUnderH1"]);
              return sum(NumericVector(nb[18])) - maxInformation;
            };

  if (unknown == "accrualDuration") {
    accrualDuration = brent(f, interval[0], interval[1], 1.0e-6);
  } else if (unknown == "followupTime") {
    followupTime = brent(f, interval[0], interval[1], 1.0e-6);
  } else if (unknown == "accrualIntensity") {
    double aval = brent(f, interval[0], interval[1], 1.0e-6);
    accrualIntensity1 = aval*accrualIntensity;
  }

  double studyDuration = accrualDuration + followupTime;


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
        auto h = [allocationRatioPlanned, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    List na = nbstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, aval, fixedFollowup, 0);
                    DataFrame nb = DataFrame(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[18])) - maxInformation;
                  };

        double lower = 0.0, upper = 1.1*followupTime;
        while (h(upper) < 0) {
          lower = upper;
          upper = 2.0*upper;
        }
        followupTime = brent(h, lower, upper, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // adjust study duration to obtain the target maximum information
        auto h = [allocationRatioPlanned, accrualTime, accrualIntensity1,
                  piecewiseSurvivalTime, stratumFraction,
                  kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                  accrualDuration, followupTime, fixedFollowup,
                  maxInformation](double aval)->double {
                    NumericVector u0(1, accrualDuration + aval);
                    List na = nbstat(
                      u0, 1, allocationRatioPlanned,
                      accrualTime, accrualIntensity1,
                      piecewiseSurvivalTime, stratumFraction,
                      kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
                      accrualDuration, followupTime, fixedFollowup, 0);
                    DataFrame nb = DataFrame(na["resultsUnderH1"]);
                    return sum(NumericVector(nb[18])) - maxInformation;
                  };

        double aval = brent(h, 0.0, followupTime, 1.0e-6);
        studyDuration = accrualDuration + aval;
      }
    }
  }

  List result = nbpowerequiv(
    kMax, informationRates1, criticalValues1,
    alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlphaSpending, rateRatioLower, rateRatioUpper,
    allocationRatioPlanned, accrualTime, accrualIntensity1,
    piecewiseSurvivalTime, stratumFraction,
    kappa1, kappa2, lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  return result;
}
