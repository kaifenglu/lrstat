#include "utilities.h"
using namespace Rcpp;

// define functions in likelihood inference
typedef struct {
  std::string dist;
  IntegerVector strata;
  NumericVector tstart;
  NumericVector tstop;
  IntegerVector status;
  NumericVector weight;
  NumericMatrix z;
} aftparams;


// negative log likelihood
double f_nllik_1(int p, double *par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  double loglik = 0;
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];

    if (param->status[person] == 1) { // event
      double logsig = log(sigma);
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*(u - exp(u) - logsig);
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*(R::dnorm(u, 0, 1, 1) - logsig);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*(R::dlogis(u, 0, 1, 1) - logsig);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*(R::dnorm(u, 0, 1, 1) - logsig);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*(R::dlogis(u, 0, 1, 1) - logsig);
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*log(exp(-exp(v)) - exp(-exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*log(R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*log(R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*log(R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*log(R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*log(1.0 - exp(-exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*R::pnorm(u, 0, 1, 1, 1);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*R::plogis(u, 0, 1, 1, 1);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*R::pnorm(u, 0, 1, 1, 1);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*R::plogis(u, 0, 1, 1, 1);
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*(-exp(v));
      } else if (param->dist == "lognormal") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*R::pnorm(v, 0, 1, 0, 1);
      } else if (param->dist == "loglogistic") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*R::plogis(v, 0, 1, 0, 1);
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*R::pnorm(v, 0, 1, 0, 1);
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*R::plogis(v, 0, 1, 0, 1);
      }
    }
  }

  return -loglik;
}


// negative score vector
void f_nscore_1(int p, double *par, double *gr, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  NumericVector score(p);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];
    NumericVector z = param->z(person, _)/sigma;
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - exp(u));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - exp(u));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*((1 - exp(u))*(-u) - 1);
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(u*u - 1);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c0 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(c0*u - 1);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(u*u - 1);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c0 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(c0*u - 1);
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(exp(v - exp(v)) - exp(u - exp(u)))/
          (exp(-exp(v)) - exp(-exp(u)));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(exp(v - exp(v)) - exp(u - exp(u)))/
          (exp(-exp(v)) - exp(-exp(u)));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(exp(v - exp(v))*v - exp(u - exp(u))*u)/
          (exp(-exp(v)) - exp(-exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(R::dnorm(v, 0, 1, 0)*v - R::dnorm(u, 0, 1, 0)*u)/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(R::dlogis(v, 0, 1, 0) - R::dlogis(u, 0, 1, 0))/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(R::dlogis(v, 0, 1, 0)*v - R::dlogis(u, 0, 1, 0)*u)/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(R::dnorm(v, 0, 1, 0)*v - R::dnorm(u, 0, 1, 0)*u)/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(R::dlogis(v, 0, 1, 0) - R::dlogis(u, 0, 1, 0))/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(R::dlogis(v, 0, 1, 0)*v - R::dlogis(u, 0, 1, 0)*u)/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-exp(u - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-exp(u - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "lognormal") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "loglogistic") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      }
    }
  }

  for (i=0; i<p; i++) gr[i] = -score[i];
}


// observed information matrix
NumericMatrix f_info_1(int p, double *par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, j, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  NumericMatrix imat(p,p);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];
    NumericVector z = param->z(person, _)/sigma;
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*exp(u);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*exp(u);
        double c2 = wt*(exp(u)*u - (1 - exp(u)));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c2 = wt*2*u;
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += wt*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*2*R::dlogis(u, 0, 1, 0);
        double c2 = wt*(2*R::dlogis(u, 0, 1, 0)*u +
                        1 - 2*R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c2 = wt*2*u;
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += wt*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*2*R::dlogis(u, 0, 1, 0);
        double c2 = wt*(2*R::dlogis(u, 0, 1, 0)*u +
                        1 - 2*R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(pow((exp(v - exp(v)) - exp(u - exp(u)))/
                        (exp(-exp(v)) - exp(-exp(u))), 2) +
                          (exp(v - exp(v))*(1 - exp(v)) -
                          exp(u - exp(u))*(1 - exp(u)))/
                            (exp(-exp(v)) - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(pow((exp(v - exp(v)) - exp(u - exp(u)))/
                        (exp(-exp(v)) - exp(-exp(u))), 2) +
                          (exp(v - exp(v))*(1 - exp(v)) -
                          exp(u - exp(u))*(1 - exp(u)))/
                            (exp(-exp(v)) - exp(-exp(u))));
        double c2 = wt*((exp(v - exp(v)) - exp(u - exp(u)))*
                        (exp(v - exp(v))*v - exp(u - exp(u))*u)/
                          pow(exp(-exp(v)) - exp(-exp(u)), 2) +
                            (exp(v - exp(v))*(1 + (1 - exp(v))*v) -
                            exp(u - exp(u))*(1 + (1 - exp(u))*u))/
                              (exp(-exp(v)) - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((exp(v - exp(v))*v - exp(u - exp(u))*u)/
          (exp(-exp(v)) - exp(-exp(u))), 2) +
            (exp(v - exp(v))*(1 + (1 - exp(v))*v)*v
            - exp(u - exp(u))*(1 + (1 - exp(u))*u)*u)/
            (exp(-exp(v)) - exp(-exp(u))));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(pow((R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
                        (R::pnorm(v, 0, 1, 0, 0) -
                          R::pnorm(u, 0, 1, 0, 0)), 2) +
                          (-R::dnorm(v, 0, 1, 0)*v +
                          R::dnorm(u, 0, 1, 0)*u)/
                            (R::pnorm(v, 0, 1, 0, 0) -
                              R::pnorm(u, 0, 1, 0, 0)));
        double c2 = wt*((R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))*
                        (R::dnorm(v, 0, 1, 0)*v - R::dnorm(u, 0, 1, 0)*u)/
                          pow(R::pnorm(v, 0, 1, 0, 0) -
                            R::pnorm(u, 0, 1, 0, 0), 2) +
                            (R::dnorm(v, 0, 1, 0)*(1 - v*v) -
                            R::dnorm(u, 0, 1, 0)*(1 - u*u))/
                              (R::pnorm(v, 0, 1, 0, 0) -
                                R::pnorm(u, 0, 1, 0, 0)));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((R::dnorm(v, 0, 1, 0)*v -
          R::dnorm(u, 0, 1, 0)*u)/
            (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0)), 2) +
              (R::dnorm(v, 0, 1, 0)*(1 - v*v)*v -
              R::dnorm(u, 0, 1, 0)*(1 - u*u)*u)/
                (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0)));
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double q1 = R::plogis(v, 0, 1, 0, 0);
        double q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*(pow((q1*(1-q1) - q2*(1-q2))/(q1-q2), 2) +
                        (q1*(1-q1)*(2*q1-1) - q2*(1-q2)*(2*q2-1))/(q1-q2));
        double c2 = wt*((q1*(1-q1) - q2*(1-q2))*
                        (q1*(1-q1)*v - q2*(1-q2)*u)/pow(q1-q2, 2) +
                        (q1*(1-q1)*(1+(2*q1-1)*v) -
                        q2*(1-q2)*(1+(2*q2-1)*u))/(q1-q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((q1*(1-q1)*v - q2*(1-q2)*u)/(q1-q2), 2) +
          (q1*(1-q1)*(1+(2*q1-1)*v)*v - q2*(1-q2)*(1+(2*q2-1)*u)*u)/(q1-q2));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(pow((R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
                        (R::pnorm(v, 0, 1, 0, 0) -
                          R::pnorm(u, 0, 1, 0, 0)), 2) +
                          (-R::dnorm(v, 0, 1, 0)*v +
                          R::dnorm(u, 0, 1, 0)*u)/
                            (R::pnorm(v, 0, 1, 0, 0) -
                              R::pnorm(u, 0, 1, 0, 0)));
        double c2 = wt*((R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))*
                        (R::dnorm(v, 0, 1, 0)*v - R::dnorm(u, 0, 1, 0)*u)/
                          pow(R::pnorm(v, 0, 1, 0, 0) -
                            R::pnorm(u, 0, 1, 0, 0), 2) +
                            (R::dnorm(v, 0, 1, 0)*(1 - v*v) -
                            R::dnorm(u, 0, 1, 0)*(1 - u*u))/
                              (R::pnorm(v, 0, 1, 0, 0) -
                                R::pnorm(u, 0, 1, 0, 0)));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((R::dnorm(v, 0, 1, 0)*v -
          R::dnorm(u, 0, 1, 0)*u)/
            (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0)), 2) +
              (R::dnorm(v, 0, 1, 0)*(1 - v*v)*v -
              R::dnorm(u, 0, 1, 0)*(1 - u*u)*u)/
                (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0)));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double q1 = R::plogis(v, 0, 1, 0, 0);
        double q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*(pow((q1*(1-q1) - q2*(1-q2))/(q1-q2), 2) +
                        (q1*(1-q1)*(2*q1-1) - q2*(1-q2)*(2*q2-1))/(q1-q2));
        double c2 = wt*((q1*(1-q1) - q2*(1-q2))*
                        (q1*(1-q1)*v - q2*(1-q2)*u)/pow(q1-q2, 2) +
                        (q1*(1-q1)*(1+(2*q1-1)*v) -
                        q2*(1-q2)*(1+(2*q2-1)*u))/(q1-q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((q1*(1-q1)*v - q2*(1-q2)*u)/(q1-q2), 2) +
          (q1*(1-q1)*(1+(2*q1-1)*v)*v - q2*(1-q2)*(1+(2*q2-1)*u)*u)/(q1-q2));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(pow(exp(u - exp(u))/(1 - exp(-exp(u))), 2) -
                        exp(u - exp(u))*(1 - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(pow(exp(u - exp(u))/(1 - exp(-exp(u))), 2) -
                        exp(u - exp(u))*(1 - exp(u))/(1 - exp(-exp(u))));
        double c2 = wt*(pow(exp(u - exp(u))/(1 - exp(-exp(u))), 2)*u -
                        exp(u - exp(u))*(1 + (1 - exp(u))*u)/
                          (1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(pow(R::dnorm(u, 0, 1, 0)/
                        R::pnorm(u, 0, 1, 1, 0), 2) +
                          R::dnorm(u, 0, 1, 0)*u/R::pnorm(u, 0, 1, 1, 0));
        double c2 = wt*(pow(R::dnorm(u, 0, 1, 0)/
                        R::pnorm(u, 0, 1, 1, 0), 2)*u -
                          R::dnorm(u, 0, 1, 0)*(1 - u*u)/
                            R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*q2*(1-q2);
        double c2 = wt*(q2*(1-q2)*u - q2);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(pow(R::dnorm(u, 0, 1, 0)/
                        R::pnorm(u, 0, 1, 1, 0), 2) +
                          R::dnorm(u, 0, 1, 0)*u/R::pnorm(u, 0, 1, 1, 0));
        double c2 = wt*(pow(R::dnorm(u, 0, 1, 0)/
                        R::pnorm(u, 0, 1, 1, 0), 2)*u -
                          R::dnorm(u, 0, 1, 0)*(1 - u*u)/
                            R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*q2*(1-q2);
        double c2 = wt*(q2*(1-q2)*u - q2);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        double c2 = wt*exp(v)*(1+v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
      } else if (param->dist == "lognormal") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(pow(R::dnorm(v, 0, 1, 0)/
                        R::pnorm(v, 0, 1, 0, 0), 2) -
                          R::dnorm(v, 0, 1, 0)*v/R::pnorm(v, 0, 1, 0, 0));
        double c2 = wt*(pow(R::dnorm(v, 0, 1, 0)/
                        R::pnorm(v, 0, 1, 0, 0), 2)*v +
                          R::dnorm(v, 0, 1, 0)*(1 - v*v)/
                            R::pnorm(v, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
      } else if (param->dist == "loglogistic") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double q1 = R::plogis(v, 0, 1, 0, 0);
        double c1 = wt*q1*(1-q1);
        double c2 = wt*(1-q1+q1*(1-q1)*v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(pow(R::dnorm(v, 0, 1, 0)/
                        R::pnorm(v, 0, 1, 0, 0), 2) -
                          R::dnorm(v, 0, 1, 0)*v/R::pnorm(v, 0, 1, 0, 0));
        double c2 = wt*(pow(R::dnorm(v, 0, 1, 0)/
                        R::pnorm(v, 0, 1, 0, 0), 2)*v +
                          R::dnorm(v, 0, 1, 0)*(1 - v*v)/
                            R::pnorm(v, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double q1 = R::plogis(v, 0, 1, 0, 0);
        double c1 = wt*q1*(1-q1);
        double c2 = wt*(1-q1+q1*(1-q1)*v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
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
NumericMatrix f_ressco_1(NumericVector par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int p = par.size();
  int person, i, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  NumericMatrix resid(n, p);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];
    NumericVector z = param->z(person, _)/sigma;
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - exp(u));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - exp(u));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*((1 - exp(u))*(-u) - 1);
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(u*u - 1);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c0 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(c0*u - 1);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(u*u - 1);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c0 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(c0*u - 1);
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(exp(v - exp(v)) - exp(u - exp(u)))/
          (exp(-exp(v)) - exp(-exp(u)));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(exp(v - exp(v)) - exp(u - exp(u)))/
          (exp(-exp(v)) - exp(-exp(u)));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(exp(v - exp(v))*v - exp(u - exp(u))*u)/
          (exp(-exp(v)) - exp(-exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(R::dnorm(v, 0, 1, 0)*v -
          R::dnorm(u, 0, 1, 0)*u)/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(R::dlogis(v, 0, 1, 0) - R::dlogis(u, 0, 1, 0))/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(R::dlogis(v, 0, 1, 0)*v -
          R::dlogis(u, 0, 1, 0)*u)/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(R::dnorm(v, 0, 1, 0)*v -
          R::dnorm(u, 0, 1, 0)*u)/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(R::dlogis(v, 0, 1, 0) - R::dlogis(u, 0, 1, 0))/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(R::dlogis(v, 0, 1, 0)*v -
          R::dlogis(u, 0, 1, 0)*u)/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-exp(u - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-exp(u - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "lognormal") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "loglogistic") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      }
    }
  }

  return resid;
}


//' @title Parametric regression models for failure time data
//' @description Obtains the parameter estimates from parametric
//' regression models with uncensored, right censored, left censored, or
//' interval censored data.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{rep}: The replication for by-group processing.
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{time}: The follow-up time for right censored data, or
//'     the left end of each interval for interval censored data.
//'
//'   * \code{time2}: The right end of each interval for interval
//'     censored data.
//'
//'   * \code{event}: The event indicator, normally 1=event, 0=no event.
//'
//'   * \code{covariates}: The values of baseline covariates.
//'     This is the full-rank design matrix (excluding the intercept)
//'     for the regression model, assuming that factor variables
//'     have already been expanded into dummy variables.
//'     The intercept will be added automatically.
//'
//'   * \code{weight}: The weight for each observation.
//'
//'   * \code{id}: The optional subject ID to group the score residuals
//'     in computing the robust sandwich variance.
//'
//' @param rep The name of the replication variable in the input data.
//' @param stratum The name of the stratum variable in the input data.
//' @param time The name of the time variable or the left end of each
//'   interval for interval censored data in the input data.
//' @param time2 The name of the right end of each interval for
//'   interval censored data in the input data.
//' @param event The name of the event variable in the input data
//'   for right censored data.
//' @param covariates The vector of names of baseline covariates
//'   in the input data.
//' @param weight The name of the weighting variable in the input data.
//' @param id The name of the id variable in the input data.
//' @param dist The assumed distribution for time to event. Options include
//'   "exponential", "weibull", "lognormal", and "loglogistic" to be
//'   modeled on the log-scale, and "normal" and "logistic" to be modeled
//'   on the original scale.
//' @param robust Whether a robust sandwich variance estimate should be
//'   computed. The default is TRUE if there are fractional weights or
//'   there is at least 1 id with >1 event. In the presence of the id
//'   variable, the score residual will be aggregated for each id when
//'   computing the robust sandwich variance estimate.
//'
//' @details There are two ways to specify the model, one for right censored
//' data through the time and event variables, and the other for interval
//' censored data through the time and time2 variables. For the second form,
//' we follow the convention used in SAS PROC LIFEREG:
//'
//' * If lower is not missing, upper is not missing, and lower is equal
//'   to upper, then there is no censoring and the event occurred at
//'   time lower.
//'
//' * If lower is not missing, upper is not missing, and lower < upper,
//'   then the event time is censored within the interval (lower, upper).
//'
//' * If lower is missing, but upper is not missing, then upper will be
//'   used as the left censoring value.
//'
//' * If lower is not missing, but upper is missing, then lower will be
//'   used as the right censoring value.
//'
//' * If lower is not missing, upper is not missing, but lower > upper,
//'   or if both lower and upper are missing, then the observation will
//'   not be used.
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
//'     - \code{beta}: The parameter estimate.
//'
//'     - \code{sebeta}: The standard error of parameter estimate.
//'
//'     - \code{z}: The Wald test statistic.
//'
//'     - \code{expbeta}: The exponentiated parameter.
//'
//'     - \code{vbeta}: The covariance matrix for parameter estimates.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' library(dplyr)
//'
//' # right censored data
//' liferegr(data = rawdata %>% mutate(treat = 1*(treatmentGroup == 1)),
//'          rep = "iterationNumber", stratum = "stratum",
//'          time = "timeUnderObservation", event = "event",
//'          covariates = "treat", dist = "weibull")
//'
//' # tobit regression for left censored data
//' liferegr(data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
//'          time = "time", time2 = "durable",
//'          covariates = c("age", "quant"), dist = "normal")
//'
//' @export
// [[Rcpp::export]]
List liferegr(const DataFrame data,
              const std::string rep = "rep",
              const std::string stratum = "stratum",
              const std::string time = "time",
              const std::string time2 = "time2",
              const std::string event = "event",
              const StringVector& covariates = "treat",
              const std::string weight = "weight",
              const std::string id = "id",
              const std::string dist = "weibull",
              bool robust = 0) {
  std::string dist1 = dist;
  std::for_each(dist1.begin(), dist1.end(), [](char & c) {
    c = std::tolower(c);
  });

  if ((dist1 == "log-logistic") || (dist1 == "llogistic")) {
    dist1 = "loglogistic";
  } else if  ((dist1 == "log-normal") || (dist1 == "lnormal")) {
    dist1 = "lognormal";
  } else if (dist1 == "gaussian") {
    dist1 = "normal";
  }

  if (!((dist1 == "exponential") || (dist1 == "weibull") ||
      (dist1 == "lognormal") || (dist1 == "loglogistic") ||
      (dist1 == "normal") || (dist1 == "logistic"))) {
    std::string str1 = "dist must be exponential, weibull, lognormal,";
    std::string str2 = "loglogistic, normal, or logistic";
    std::string errmsg = str1 + " " + str2;
    stop(errmsg);
  }

  int h, i, j, k, n = data.nrows();
  int nvar = covariates.size() + 1;

  bool has_rep = hasVariable(data, rep);
  bool has_stratum = hasVariable(data, stratum);

  bool has_time = hasVariable(data, time);
  if (!has_time) {
    stop("data must contain the time variable");
  }

  NumericVector timen = data[time];
  for (i=0; i<n; i++) {
    if (!R_isnancpp(timen[i]) && ((dist1 == "exponential") ||
        (dist1 == "weibull") || (dist1 == "lognormal") ||
        (dist1 == "loglogistic")) && (timen[i] <= 0)) {
      std::string str1 = "time must be positive for each subject for the";
      std::string str2 = "distribution";
      std::string errmsg = str1 + " " + dist1 + " " + str2;
      stop(errmsg);
    }
  }

  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    time2n = data[time2];
    for (i=0; i<n; i++) {
      if (!R_isnancpp(time2n[i]) && ((dist1 == "exponential") ||
          (dist1 == "weibull") || (dist1 == "lognormal") ||
          (dist1 == "loglogistic")) && (time2n[i] <= 0)) {
        std::string str1 = "time2 must be positive for each subject for the";
        std::string str2 = "distribution";
        std::string errmsg = str1 + " " + dist1 + " " + str2;
        stop(errmsg);
      }
    }
  }

  bool has_event = hasVariable(data, event);
  if (!has_time2 && !has_event) {
    stop("data must contain the event variable for right censored data");
  }

  IntegerVector eventn(n);
  if (has_event) {
    eventn = data[event];
    if (is_true(any((eventn != 1) & (eventn != 0)))) {
      stop("event must be 1 or 0 for each subject");
    }

    if (is_true(all(eventn == 0))) {
      stop("at least 1 event is needed to fit the parametric model");
    }
  }

  NumericMatrix zn(n,nvar);
  for (i=0; i<n; i++) zn(i,0) = 1; // intercept

  for (j=0; j<nvar-1; j++) {
    String zj = covariates[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in covariates");
    }
    NumericVector u = data[zj];
    for (i=0; i<n; i++) zn(i,j+1) = u[i];
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

  IntegerVector stratumn1 = unique(stratumn);
  int nstrata = stratumn1.size();
  int p = dist1 == "exponential" ? nvar : (nvar+nstrata);

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

  NumericMatrix z(n,nvar);
  for (i=0; i<n; i++) {
    for (j=0; j<nvar; j++) {
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

    NumericMatrix z1(n1,nvar);
    for (i=0; i<n1; i++) {
      for (j=0; j<nvar; j++) {
        z1(i,j) = zn(q1[i],j);
      }
    }

    // unify right censored data with interval censored data
    NumericVector tstart(n1), tstop(n1);
    if (!has_time2) {
      tstart = time1;
      for (i=0; i<n1; i++) {
        tstop[i] = event1[i] == 1 ? tstart[i] : NA_REAL;
      }
    } else {
      tstart = time1;
      tstop = time21;
    }

    IntegerVector status(n1);
    for (i=0; i<n1; i++) {
      if (!R_isnancpp(tstart[i]) && !R_isnancpp(tstop[i]) &&
        (tstart[i] == tstop[i])) {
        status[i] = 1; // event
      } else if (!R_isnancpp(tstart[i]) && !R_isnancpp(tstop[i]) &&
        (tstart[i] < tstop[i])) {
        status[i] = 3; // interval censoring
      } else if (R_isnancpp(tstart[i]) && !R_isnancpp(tstop[i])) {
        status[i] = 2; // left censoring
      } else if (!R_isnancpp(tstart[i]) && R_isnancpp(tstop[i])) {
        status[i] = 0; // right censoring
      } else {
        status[i] = -1; // exclude the observation
      }
    }

    nobs[h] = n1;
    nevents[h] = sum(status == 1);

    // initial parameter values
    NumericVector time0(n1);
    for (i=0; i<n1; i++) {
      if (status[i] == 1) { // event
        time0[i] = tstart[i];
      } else if (status[i] == 3) { // interval censoring
        time0[i] = (tstart[i] + tstop[i])/2;
      } else if (status[i] == 2) { // left censoring
        time0[i] = tstop[i];
      } else if (status[i] == 0) { // right censoring
        time0[i] = tstart[i];
      } else {
        time0[i] = NA_REAL;
      }
    }

    // intercept only model
    LogicalVector sub = !is_na(time0);
    NumericVector y0 = time0[sub];
    if ((dist1 == "exponential") || (dist1 == "weibull") ||
        (dist1 == "lognormal") || (dist1 == "loglogistic")) {
      y0 = log(y0);
    }

    double int0 = mean(y0);
    double logsig0 = log(sd(y0));

    int pint = dist1 == "exponential" ? 1 : (nstrata+1);
    NumericVector bint0(pint);
    if (dist1 == "exponential") {
      bint0[0] = int0;
    } else {
      bint0[0] = int0;
      for (i=0; i<nstrata; i++) {
        bint0[i+1] = logsig0;
      }
    }

    NumericMatrix zi(n1,1);
    for (i=0; i<n1; i++) zi(i,0) = 1;

    // parameter estimates and standard errors for the null model
    aftparams parami = {dist1, stratum1, tstart, tstop, status, weight1, zi};
    List outint = bmini(bint0, f_nllik_1, f_nscore_1, &parami, 1e-9);
    NumericVector bint = outint["par"];
    std::vector<double> parbint(bint.begin(), bint.end());

    NumericVector b0(p);
    if (dist1 == "exponential") {
      b0[0] = bint[0];
    } else {
      b0[0] = bint[0];
      for (i=0; i<nstrata; i++) {
        b0[nvar+i] = bint[i+1];
      }
    }

    // exclude observations with missing covariates
    for (i=0; i<n1; i++) sub[i] = is_false(any(is_na(z1(i,_))));
    int nsub = sum(sub);

    if (nsub < n1) {
      stratum1 = stratum1[sub];
      tstart = tstart[sub];
      tstop = tstop[sub];
      status = status[sub];
      weight1 = weight1[sub];
      id1 = id1[sub];

      NumericMatrix z2(nsub, nvar);
      j = 0;
      for (i = 0; i < n1; i++) {
        if (sub[i]) {
          z2(j, _) = z1(i, _);
          j++;
        }
      }
      z1 = z2;
    }

    // parameter estimates and standard errors for the full model
    aftparams param = {dist1, stratum1, tstart, tstop, status, weight1, z1};
    List out = bmini(b0, f_nllik_1, f_nscore_1, &param, 1e-9);
    NumericVector b = out["par"];

    std::vector<double> parb(b.begin(), b.end());
    NumericMatrix infob = f_info_1(p, parb.data(), &param);
    NumericMatrix vb = invsympd(infob);

    NumericVector seb(p);
    for (j=0; j<p; j++) {
      seb[j] = sqrt(vb(j,j));
    }

    for (i=0; i<p; i++) {
      rep0[h*p+i] = h+1;

      if (i==0) {
        par0[h*p+i] = "(Intercept)";
      } else if (i < nvar) {
        par0[h*p+i] = covariates[i-1];
      } else {
        if (nstrata == 1) {
          par0[h*p+i] = "Log(scale)";
        } else {
          std::string str1 = "Log(scale ";
          std::string str2 = ")";
          par0[h*p+i] = str1 + std::to_string(i-nvar+1) + str2;
        }
      }

      beta0[h*p+i] = b[i];
      sebeta0[h*p+i] = seb[i];
      for (j=0; j<p; j++) {
        vbeta0(h*p+i,j) = vb(i,j);
      }
    }

    // log-likelihoods and score test statistic
    std::vector<double> parb0(b0.begin(), b0.end());
    loglik(h,0) = -f_nllik_1(p, parb0.data(), &param);
    loglik(h,1) = -as<double>(out["value"]);

    NumericVector score(p);
    std::vector<double> pars(score.begin(), score.end());
    f_nscore_1(p, parb0.data(), pars.data(), &param);
    for (j=0; j<p; j++) score[j] = -pars[j];

    NumericMatrix infob0 = f_info_1(p, parb0.data(), &param);
    NumericMatrix vb0 = invsympd(infob0);
    for (i=1; i<nvar; i++) {
      for (j=1; j<nvar; j++) {
        scoretest[h] += score[i]*vb0(i,j)*score[j];
      }
    }

    // robust variance estimates
    if (robust) {
      NumericMatrix ressco = f_ressco_1(b, &param); // score residuals

      int nr; // number of rows in the score residual matrix
      if (!has_id) {
        nr = nsub;
      } else { // need to sum up score residuals by id
        IntegerVector order = seq(0, nsub-1);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return id1[i] < id1[j];
        });

        IntegerVector id2 = id1[order];
        IntegerVector idx(1,0);
        for (i=1; i<nsub; i++) {
          if (id2[i] != id2[i-1]) {
            idx.push_back(i);
          }
        }

        int nids = idx.size();
        idx.push_back(nsub);

        NumericMatrix resid(nsub,p);
        for (i=0; i<nsub; i++) {
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

  NumericVector expbeta0 = exp(beta0);
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
      _["expbeta"] = expbeta0,
      _["vbeta"] = vbeta0);
  } else {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = sebeta0,
      _["rsebeta"] = rsebeta0,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
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
