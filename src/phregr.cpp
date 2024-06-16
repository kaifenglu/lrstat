#include "utilities.h"
using namespace Rcpp;

// define functions in likelihood inference, algorithms adapted from coxph
typedef struct {
  int nused;
  IntegerVector strata;
  NumericVector tstart;
  NumericVector tstop;
  IntegerVector event;
  NumericVector weight;
  NumericMatrix z;
  IntegerVector order1;
  int method;
} coxparams;


// negative log likelihood
double f_nllik_2(int p, double *par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, k, person;

  NumericVector beta(p);
  for (i=0; i<p; i++) beta[i] = par[i];

  double loglik = 0;        // log-likelihood value
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double zbeta;             // linear predictor
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    zbeta = 0;
    for (i=0; i<p; i++) zbeta += beta[i]*param->z(person,i);
    eta[person] = zbeta;
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0;
    }

    dtime = param->tstop[person];
    while (person < param->nused && param->tstop[person] == dtime) {
      // walk through this set of tied times
      risk = param->weight[person]*exp(eta[person]);
      if (param->event[person] == 0) {
        denom += risk;
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
        loglik += param->weight[person]*eta[person];
      }

      person++;

      if (person < param->nused && param->strata[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1 < param->nused; i1++) {
      p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      risk = param->weight[p1]*exp(eta[p1]);
      denom -= risk;
    }

    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        loglik -= deadwt*log(denom);
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          loglik -= meanwt*log(denom);
        }
      }

      // reset for the next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
    }
  }

  return -loglik;
}


// negative score vector
void f_nscore_2(int p, double *par, double *gr, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, k, person;

  NumericVector beta(p);
  for (i=0; i<p; i++) beta[i] = par[i];

  NumericVector u(p);       // score vector
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double zbeta;             // linear predictor
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    zbeta = 0;
    for (i=0; i<p; i++) zbeta += beta[i]*param->z(person,i);
    eta[person] = zbeta;
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0;
      for (i=0; i<p; i++) a[i] = 0;
    }

    dtime = param->tstop[person];
    while (person < param->nused && param->tstop[person] == dtime) {
      // walk through this set of tied times
      risk = param->weight[person]*exp(eta[person]);
      if (param->event[person] == 0) {
        denom += risk;
        for (i=0; i<p; i++) a[i] += risk*param->z(person,i);
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
        for (i=0; i<p; i++) {
          a2[i] += risk*param->z(person,i);
          u[i] += param->weight[person]*param->z(person,i);
        }
      }

      person++;

      if (person < param->nused && param->strata[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1<param->nused; i1++) {
      p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      risk = param->weight[p1]*exp(eta[p1]);
      denom -= risk;
      for (i=0; i<p; i++) a[i] -= risk*param->z(p1,i);
    }

    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          xbar = a[i]/denom;
          u[i] -= deadwt*xbar;
        }
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar = a[i]/denom;
            u[i] -= meanwt*xbar;
          }
        }
      }

      // reset for the next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) a2[i] = 0;
    }
  }

  for (i=0; i<p; i++) gr[i] = -u[i];
}


// observed information matrix
NumericMatrix f_info_2(int p, double *par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, j, k, person;

  NumericVector beta(p);
  for (i=0; i<p; i++) beta[i] = par[i];

  NumericMatrix imat(p,p);  // information matrix
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double zbeta;             // linear predictor
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)
  NumericMatrix cmat(p,p);  // s2(beta,k,t)
  NumericMatrix cmat2(p,p); // sum of w*exp(zbeta)*z*z' for the deaths

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    zbeta = 0;
    for (i=0; i<p; i++) zbeta += beta[i]*param->z(person,i);
    eta[person] = zbeta;
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) {  // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0;
      for (i=0; i<p; i++) {
        a[i] = 0;
        for (j=0; j<p; j++) cmat(i,j) = 0;
      }
    }

    dtime = param->tstop[person];
    while (person < param->nused && param->tstop[person] == dtime) {
      // walk through this set of tied times
      risk = param->weight[person]*exp(eta[person]);
      if (param->event[person] == 0) {
        denom += risk;
        for (i=0; i<p; i++) {
          a[i] += risk*param->z(person,i);
          for (j=0; j<=i; j++) {
            cmat(i,j) += risk*param->z(person,i)*param->z(person,j);
          }
        }
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
        for (i=0; i<p; i++) {
          a2[i] += risk*param->z(person,i);
          for (j=0; j<=i; j++) {
            cmat2(i,j) += risk*param->z(person,i)*param->z(person,j);
          }
        }
      }

      person++;

      if (person < param->nused && param->strata[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1<param->nused; i1++) {
      p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      risk = param->weight[p1]*exp(eta[p1]);
      denom -= risk;
      for (i=0; i<p; i++) {
        a[i] -= risk*param->z(p1,i);
        for (j=0; j<=i; j++) {
          cmat(i,j) -= risk*param->z(p1,i)*param->z(p1,j);
        }
      }
    }

    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          xbar = a[i]/denom;
          for (j=0; j<=i; j++) {
            cmat(i,j) += cmat2(i,j);
            imat(i,j) += deadwt*(cmat(i,j) - xbar*a[j])/denom;
          }
        }
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar = a[i]/denom;
            for (j=0; j<=i; j++) {
              cmat(i,j) += cmat2(i,j)/ndead;
              imat(i,j) += meanwt*(cmat(i,j) - xbar*a[j])/denom;
            }
          }
        }
      }

      // reset for next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) {
        a2[i] = 0;
        for (j=0; j<=i; j++) cmat2(i,j) = 0;
      }
    }
  }

  for (i=0; i<p-1; i++) {
    for (j=i+1; j<p; j++) {
      imat(i,j) = imat(j,i);
    }
  }

  return imat;
}


// score residual matrix
NumericMatrix f_ressco_2(NumericVector beta, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, j, k, person;
  int p = beta.size();

  NumericMatrix resid(param->tstart.size(), p);
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double zbeta;             // linear predictor
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)

  double downwt, hazard, cumhaz = 0;
  NumericVector xhaz(p), mh1(p), mh2(p), mh3(p);

  NumericVector score(param->nused);
  for (person = 0; person < param->nused; person++) {
    zbeta = 0;
    for (i=0; i<p; i++) zbeta += beta[i]*param->z(person,i);
    score[person] = exp(zbeta);
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person];
      // first obs of a new stratum, finish off the prior stratum
      for (; i1 < param->nused && param->order1[i1] < person; i1++) {
        p1 = param->order1[i1];
        for (i=0; i<p; i++) {
          resid(p1,i) -= score[p1]*(param->z(p1,i)*cumhaz - xhaz[i]);
        }
      }
      denom = 0; // reset temporary variables
      cumhaz = 0;
      for (i=0; i<p; i++) {
        a[i] =0;
        xhaz[i] =0;
      }
    }

    dtime = param->tstop[person];
    while (person < param->nused && param->tstop[person] == dtime) {
      // walk through this set of tied times
      // initialize residuals to score[i]*(x[i]*cumhaz - xhaz), before
      // updating cumhaz and xhaz
      for (i=0; i<p; i++) {
        resid(person,i) = score[person]*(param->z(person,i)*cumhaz - xhaz[i]);
      }

      risk = param->weight[person]*score[person];
      if (param->event[person] == 0) {
        denom += risk;
        for (i=0; i<p; i++) a[i] += risk*param->z(person,i);
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
        for (i=0; i<p; i++) {
          a2[i] += risk*param->z(person,i);
        }
      }

      person++;

      if (person < param->nused && param->strata[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1<param->nused; i1++) {
      p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      risk = param->weight[p1]*score[p1];
      denom -= risk;
      for (i=0; i<p; i++) {
        // finish the residual by subtracting score[i]*(x[i]*cumhaz - xhaz)
        resid(p1,i) -= score[p1]*(param->z(p1,i)*cumhaz - xhaz[i]);
        a[i] -= risk*param->z(p1,i);
      }
    }

    // update the cumulative sums at death times
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        hazard = deadwt/denom;
        cumhaz += hazard;
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          xbar = a[i]/denom;
          xhaz[i] += xbar*hazard;
          for (j=person-1; j>=person-ndead; j--) {
            resid(j,i) += param->z(j,i) - xbar;
          }
        }
      } else {
        for (i=0; i<p; i++) {
          mh1[i] = 0;
          mh2[i] = 0;
          mh3[i] = 0;
        }
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          hazard = meanwt/denom;
          cumhaz += hazard;
          downwt = (ndead-k-1.0)/ndead;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar = a[i]/denom;
            xhaz[i] += xbar*hazard;
            mh1[i]  += hazard*downwt;
            mh2[i]  += xbar*hazard*downwt;
            mh3[i]  += xbar/ndead;
          }
        }

        for (j=person-1; j>=person-ndead; j--) {
          for (i=0; i<p; i++) {
            resid(j,i) += (param->z(j,i) - mh3[i]) +
              score[j]*(param->z(j,i)*mh1[i] - mh2[i]);
          }
        }
      }

      // reset for next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) a2[i] = 0;
    }
  }

  // finish those remaining in the final stratum
  for (; i1<param->nused; i1++) {
    p1 = param->order1[i1];
    for (i=0; i<p; i++)
      resid(p1,i) -= score[p1]*(param->z(p1,i)*cumhaz - xhaz[i]);
  }

  return resid;
}


//' @title Proportional hazards regression model
//' @description Obtains the hazard ratio estimates from the proportional
//' hazards regression model with right censored or counting process data.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{rep}: The replication for by-group processing.
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{time}: The follow-up time for right censored data, or
//'     the left end of each interval for counting process data.
//'
//'   * \code{time2}: The right end of each interval for counting process
//'     data only. Intervals are assumed to be open on the left
//'     and closed on the right, and event indicates whether an event
//'     occurred at the right end of each interval.
//'
//'   * \code{event}: The event indicator, normally 1=event, 0=no event.
//'
//'   * \code{covariates}: The values of baseline covariates (and
//'     time-dependent covariates in each interval for counting
//'     process data). This is the full-rank design matrix for the Cox
//'     model, assuming that factor variables have already been
//'     expanded into dummy variables.
//'
//'   * \code{weight}: The weight for each observation.
//'
//'   * \code{id}: The optional subject ID for counting process data
//'     with time-dependent covariates.
//'
//' @param rep The name of the replication variable in the input data.
//' @param stratum The name of the stratum variable in the input data.
//' @param time The name of the time variable or the left end of each
//'   interval for counting process data in the input data.
//' @param time2 The name of the right end of each interval for counting
//'   process data in the input data.
//' @param event The name of the event variable in the input data.
//' @param covariates The vector of names of baseline and time-dependent
//'   covariates in the input data.
//' @param weight The name of the weighting variable in the input data.
//' @param id The name of the id variable in the input data.
//' @param ties The method for handling ties with options including
//'   "breslow" and "efron" (default).
//' @param robust Whether a robust sandwich variance estimate should be
//'   computed. The default is TRUE if there are fractional weights or
//'   there is at least 1 id with >1 event. In the presence of the id
//'   variable, the score residual will be aggregated for each id when
//'   computing the robust sandwich variance estimate.
//'
//' @return A list with the following components:
//'
//' * \code{sumstat}: The data frame of summary statistics of model fit
//'   with the following variables:
//'
//'     - \code{rep}: The replication.
//'
//'     - \code{n}: The number of observations.
//'
//'     - \code{nevents}: The number of events.
//'
//'     - \code{loglik0}: The log-likelihood under null.
//'
//'     - \code{loglik1}: The maximum log-likelihood.
//'
//'     - \code{scoretest}: The score test statistic.
//'
//' * \code{parest}: The data frame of parameter estimates with the
//'   following variables:
//'
//'     - \code{rep}: The replication.
//'
//'     - \code{param}: The name of the covariate for the parameter estimate.
//'
//'     - \code{beta}: The log hazard ratio estimate.
//'
//'     - \code{sebeta}: The standard error of log hazard ratio estimate.
//'
//'     - \code{rsebeta}: The robust standard error of log hazard ratio
//'       estimate if robust variance is requested.
//'
//'     - \code{z}: The Wald test statistic for log hazard ratio. The
//'       \code{rsebeta} will be used if robust variance is requested.
//'
//'     - \code{hazardRatio}: The hazard ratio estimate.
//'
//'     - \code{vbeta}: The covariance matrix for parameter estimates.
//'
//'     - \code{rvbeta}: The robust covariance matrix for parameter
//'       estimates if robust variance is requested.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' library(dplyr)
//'
//' # Example 1 with right-censored data
//' phregr(data = rawdata %>% mutate(treat = 1*(treatmentGroup == 1)),
//'        rep = "iterationNumber", stratum = "stratum",
//'        time = "timeUnderObservation", event = "event",
//'        covariates = "treat")
//'
//' # Example 2 with counting process data and robust variance estimate
//' phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
//'        time = "start", time2 = "stop", event = "event",
//'        covariates = c("rx", "age"), id = "id", robust = 1)
//'
//' @export
// [[Rcpp::export]]
List phregr(const DataFrame data,
            const std::string rep = "rep",
            const std::string stratum = "stratum",
            const std::string time = "time",
            const std::string time2 = "time2",
            const std::string event = "event",
            const StringVector& covariates = "treat",
            const std::string weight = "weight",
            const std::string id = "id",
            const std::string ties = "efron",
            bool robust = 0) {

  int h, i, j, k, n = data.nrows(), p = covariates.size();

  bool has_rep = hasVariable(data, rep);
  bool has_stratum = hasVariable(data, stratum);

  bool has_time = hasVariable(data, time);
  if (!has_time) {
    stop("data must contain the time variable");
  }

  NumericVector timen = data[time];
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each subject");
  }

  bool has_time2 = hasVariable(data, time2);

  NumericVector time2n(n);
  if (has_time2) {
    time2n = data[time2];
    if (is_true(any(time2n <= timen))) {
      stop("time2 must be greater than time for each observation");
    }
  }

  bool has_event = hasVariable(data, event);
  if (!has_event) {
    stop("data must contain the event variable");
  }

  IntegerVector eventn = data[event];
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }

  if (is_true(all(eventn == 0))) {
    stop("at least 1 event is needed to fit the Cox model");
  }

  NumericMatrix zn(n,p);
  for (j=0; j<p; j++) {
    String zj = covariates[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in covariates");
    }
    NumericVector u = data[zj];
    for (i=0; i<n; i++) zn(i,j) = u[i];
  }

  bool has_weight = hasVariable(data, weight);

  NumericVector weightn(n, 1.0);
  if (has_weight) {
    weightn = data[weight];
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }

  bool has_id = hasVariable(data, id);

  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = std::tolower(c);
  });

  int method = meth == "efron" ? 1 : 0;

  // create the numeric rep variable
  IntegerVector repn(n);
  IntegerVector repwi;
  NumericVector repwn;
  StringVector repwc;
  if (!has_rep) {
    repn.fill(1);
  } else {
    if (TYPEOF(data[rep]) == INTSXP) {
      IntegerVector repv = data[rep];
      repwi = unique(repv);
      repwi.sort();
      repn = match(repv, repwi);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      NumericVector repv = data[rep];
      repwn = unique(repv);
      repwn.sort();
      repn = match(repv, repwn);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      StringVector repv = data[rep];
      repwc = unique(repv);
      repwc.sort();
      repn = match(repv, repwc);
    } else {
      stop("incorrect type for the rep variable in the input data");
    }
  }

  // create the numeric stratum variable
  IntegerVector stratumn(n);
  if (!has_stratum) {
    stratumn.fill(1);
  } else {
    if (TYPEOF(data[stratum]) == INTSXP) {
      IntegerVector stratumv = data[stratum];
      IntegerVector stratumwi = unique(stratumv);
      stratumwi.sort();
      stratumn = match(stratumv, stratumwi);
    } else if (TYPEOF(data[stratum]) == REALSXP) {
      NumericVector stratumv = data[stratum];
      NumericVector stratumwn = unique(stratumv);
      stratumwn.sort();
      stratumn = match(stratumv, stratumwn);
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      StringVector stratumv = data[stratum];
      StringVector stratumwc = unique(stratumv);
      stratumwc.sort();
      stratumn = match(stratumv, stratumwc);
    } else {
      stop("incorrect type for the stratum variable in the input data");
    }
  }

  // create the numeric id variable
  IntegerVector idn(n);
  if (!has_id) {
    idn = seq(1,n);
  } else {
    if (TYPEOF(data[id]) == INTSXP) {
      IntegerVector idv = data[id];
      IntegerVector idwi = unique(idv);
      idwi.sort();
      idn = match(idv, idwi);
    } else if (TYPEOF(data[id]) == REALSXP) {
      NumericVector idv = data[id];
      NumericVector idwn = unique(idv);
      idwn.sort();
      idn = match(idv, idwn);
    } else if (TYPEOF(data[id]) == STRSXP) {
      StringVector idv = data[id];
      StringVector idwc = unique(idv);
      idwc.sort();
      idn = match(idv, idwc);
    } else {
      stop("incorrect type for the id variable in the input data");
    }
  }

  // check if there is as least one id with more than 1 event
  IntegerVector idn0 = idn[eventn==1];
  IntegerVector idw0 = unique(idn0);
  bool dup_id = idw0.size() < idn0.size();

  if (R_isnancpp(robust)) {
    if ((has_id && dup_id) || is_true(any(weightn != floor(weightn)))) {
      robust = 1;
    } else {
      robust = 0;
    }
  }

  if (robust && has_time2 && !has_id) {
    stop("id is needed for counting process data with robust variance");
  }

  // create an unique id variable to recover the order of input data set
  IntegerVector uid = seq(1,n);

  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  repn = repn[order];
  stratumn = stratumn[order];
  timen = timen[order];
  time2n = time2n[order];
  eventn = eventn[order];
  weightn = weightn[order];
  idn = idn[order];
  uid = uid[order];

  NumericMatrix z(n,p);
  for (i=0; i<n; i++) {
    for (j=0; j<p; j++) {
      z(i,j) = zn(order[i], j);
    }
  }
  zn = z;

  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (i=1; i<n; i++) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }

  int nreps = idx.size();
  idx.push_back(n);

  // variables in the output data sets
  IntegerVector rep01 = seq(1,nreps);
  IntegerVector nobs(nreps), nevents(nreps);
  NumericMatrix loglik(nreps,2);
  NumericVector scoretest(nreps);

  IntegerVector rep0(nreps*p);
  StringVector par0(nreps*p);
  NumericVector beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p);
  NumericMatrix vbeta0(nreps*p,p), rvbeta0(nreps*p,p);

  for (h=0; h<nreps; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = q1.size();

    IntegerVector stratum1 = stratumn[q1];
    NumericVector time1 = timen[q1];
    NumericVector time21 = time2n[q1];
    IntegerVector event1 = eventn[q1];
    NumericVector weight1 = weightn[q1];
    IntegerVector id1 = idn[q1];
    IntegerVector uid1 = uid[q1];

    NumericMatrix z1(n1,p);
    for (i=0; i<n1; i++) {
      for (j=0; j<p; j++) {
        z1(i,j) = zn(q1[i],j);
      }
    }

    nobs[h] = n1;
    nevents[h] = sum(event1);

    // unify right censored data with counting process data
    NumericVector tstart(n1), tstop(n1);
    if (!has_time2) {
      tstop = time1;
    } else {
      tstart = time1;
      tstop = time21;
    }

    // ignore subjects not at risk for any event time
    double delta = max(tstop) + 1.0; // ensure no overlap between strata
    for (i=0; i<n1; i++) {
      tstart[i] = tstart[i] + stratum1[i]*delta;
      tstop[i] = tstop[i] + stratum1[i]*delta;
    }

    NumericVector etime = tstop[event1==1];
    etime = unique(etime);
    etime.sort();

    IntegerVector index1 = findInterval3(tstart, etime);
    IntegerVector index2 = findInterval3(tstop, etime);
    IntegerVector ignore1(n1);
    for (i=0; i<n1; i++) {
      if ((index1[i] == index2[i]) || is_true(any(is_na(z1(i,_))))) {
        ignore1[i] = 1;
      } else {
        ignore1[i] = 0;
      }
    }

    int nused = n1 - sum(ignore1);

    // sort by descending stopping time within each stratum
    IntegerVector order2 = seq(0, n1-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      return (ignore1[i] < ignore1[j]) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] < stratum1[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstop[i] > tstop[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstop[i] == tstop[j]) && (event1[i] < event1[j]));
    });

    stratum1 = stratum1[order2];
    tstart = tstart[order2];
    tstop = tstop[order2];
    event1 = event1[order2];
    weight1 = weight1[order2];
    id1 = id1[order2];
    uid1 = uid1[order2];
    ignore1 = ignore1[order2];

    NumericMatrix zy(n1,p);
    for (i=0; i<n1; i++) {
      for (j=0; j<p; j++) {
        zy(i,j) = z1(order2[i],j);
      }
    }
    z1 = zy;

    // sort by starting time in descending order within each stratum
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return (ignore1[i] < ignore1[j]) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] < stratum1[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstart[i] > tstart[j]));
    });


    // parameter estimates and standard errors
    coxparams param = {nused, stratum1, tstart, tstop, event1,
                       weight1, z1, order1, method};

    NumericVector b0(p);
    List out = bmini(b0, f_nllik_2, f_nscore_2, &param, 1e-9);
    NumericVector b = out["par"];

    std::vector<double> parb(b.begin(), b.end());
    NumericMatrix infob = f_info_2(p, parb.data(), &param);
    NumericMatrix vb = invsympd(infob);

    NumericVector seb(p);
    for (j=0; j<p; j++) seb[j] = sqrt(vb(j,j));

    for (i=0; i<p; i++) {
      rep0[h*p+i] = h+1;
      par0[h*p+i] = covariates[i];
      beta0[h*p+i] = b[i];
      sebeta0[h*p+i] = seb[i];
      for (j=0; j<p; j++) {
        vbeta0(h*p+i,j) = vb(i,j);
      }
    }

    // log-likelihoods and score test statistic
    std::vector<double> parb0(b0.begin(), b0.end());
    loglik(h,0) = -f_nllik_2(p, parb0.data(), &param);
    loglik(h,1) = -as<double>(out["value"]);

    NumericVector score(p);
    std::vector<double> pars(score.begin(), score.end());
    f_nscore_2(p, parb0.data(), pars.data(), &param);
    for (j=0; j<p; j++) score[j] = -pars[j];

    NumericMatrix infob0 = f_info_2(p, parb0.data(), &param);
    NumericMatrix vb0 = invsympd(infob0);
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        scoretest[h] += score[i]*vb0(i,j)*score[j];
      }
    }

    // robust variance estimates
    if (robust) {
      NumericMatrix ressco = f_ressco_2(b, &param); // score residuals

      int nr; // number of rows in the score residual matrix
      if (!has_id) {
        nr = n1;
      } else { // need to sum up score residuals by id
        IntegerVector order = seq(0, n1-1);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return id1[i] < id1[j];
        });

        IntegerVector id2 = id1[order];
        IntegerVector idx(1,0);
        for (i=1; i<n1; i++) {
          if (id2[i] != id2[i-1]) {
            idx.push_back(i);
          }
        }

        int nids = idx.size();
        idx.push_back(n1);

        NumericMatrix resid(n1,p);
        for (i=0; i<n1; i++) {
          for (j=0; j<p; j++) {
            resid(i,j) += ressco(order[i],j);
          }
        }

        NumericMatrix ressco2(nids,p);
        for (i=0; i<nids; i++) {
          for (j=0; j<p; j++) {
            for (k=idx[i]; k<idx[i+1]; k++) {
              ressco2(i,j) += resid(k,j);
            }
          }
        }

        ressco = ressco2;  // update the score residuals
        nr = nids;
      }

      NumericMatrix D(nr,p); // DFBETA
      for (i=0; i<nr; i++) {
        for (j=0; j<p; j++) {
          for (k=0; k<p; k++) {
            D(i,j) += weight1[i]*ressco(i,k)*vb(k,j);
          }
        }
      }

      NumericMatrix rvb(p,p); // robust variance matrix for betahat
      for (j=0; j<p; j++) {
        for (k=0; k<p; k++) {
          for (i=0; i<nr; i++) {
            rvb(j,k) += D(i,j)*D(i,k);
          }
        }
      }

      NumericVector rseb(p);  // robust standard error for betahat
      for (i=0; i<p; i++) rseb[i] = sqrt(rvb(i,i));

      for (i=0; i<p; i++) {
        rsebeta0[h*p+i] = rseb[i];
        for (j=0; j<p; j++) {
          rvbeta0(h*p+i,j) = rvb(i,j);
        }
      }
    }
  }

  NumericVector hazardRatio0 = exp(beta0);
  NumericVector z0(nreps*p);
  if (!robust) z0 = beta0/sebeta0;
  else z0 = beta0/rsebeta0;

  DataFrame sumstat = List::create(
    _["n"] = nobs,
    _["nevents"] = nevents,
    _["loglik0"] = loglik(_,0),
    _["loglik1"] = loglik(_,1),
    _["scoretest"] = scoretest);

  DataFrame parest;
  if (!robust) {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = sebeta0,
      _["z"] = z0,
      _["hazardRatio"] = hazardRatio0,
      _["vbeta"] = vbeta0);
  } else {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = sebeta0,
      _["rsebeta"] = rsebeta0,
      _["z"] = z0,
      _["hazardRatio"] = hazardRatio0,
      _["vbeta"] = vbeta0,
      _["rvbeta"] = rvbeta0);
  }

  if (has_rep) {
    if (TYPEOF(data[rep]) == INTSXP) {
      sumstat.push_back(repwi[rep01-1], rep);
      parest.push_back(repwi[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      sumstat.push_back(repwn[rep01-1], rep);
      parest.push_back(repwn[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      sumstat.push_back(repwc[rep01-1], rep);
      parest.push_back(repwc[rep0-1], rep);
    }
  }

  List result = List::create(
      _["sumstat"] = sumstat,
      _["parest"] = parest);

  return result;
}
