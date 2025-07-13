#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "utilities.h"
#include "logistic_regression.h"

using namespace Rcpp;


// log likelihood
double f_llik_0(int p, NumericVector par, void *ex) {
  logparams *param = (logparams *) ex;
  int n = param->n;
  int person, i;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  double loglik = 0;
  if (param->link == "logit") {
    for (person = 0; person < n; person++) {
      double f = param->freq[person];
      double w = param->weight[person];
      double y = param->y[person];
      double r = R::plogis(eta[person], 0, 1, 1, 0);
      double v = y*eta[person] + log(1-r);
      loglik += f*w*v;
    }
  } else if (param->link == "probit") {
    for (person = 0; person < n; person++) {
      double f = param->freq[person];
      double w = param->weight[person];
      double y = param->y[person];
      double r = R::pnorm(eta[person], 0, 1, 1, 0);
      double v = y*log(r/(1-r)) + log(1-r);
      loglik += f*w*v;
    }
  } else if (param->link == "cloglog") {
    for (person = 0; person < n; person++) {
      double f = param->freq[person];
      double w = param->weight[person];
      double y = param->y[person];
      double r = 1 - exp(-exp(eta[person]));
      double v = y*log(r/(1-r)) + log(1-r);
      loglik += f*w*v;
    }
  }

  return loglik;
}


// score vector
NumericVector f_score_0(int p, NumericVector par, void *ex) {
  logparams *param = (logparams *) ex;
  int n = param->n;
  int person, i;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector score(p);
  if (param->link == "logit") {
    for (person = 0; person < n; person++) {
      double f = param->freq[person];
      double w = param->weight[person];
      double r = R::plogis(eta[person], 0, 1, 1, 0);
      double v = param->y[person] - r;
      NumericVector z = param->z(person, _);
      for (i=0; i<p; i++) {
        score[i] += f*w*v*z[i];
      }
    }
  } else if (param->link == "probit") {
    for (person = 0; person < n; person++) {
      double f = param->freq[person];
      double w = param->weight[person];
      double r = R::pnorm(eta[person], 0, 1, 1, 0);
      double phi = R::dnorm(eta[person], 0, 1, 0);
      double d = phi/(r*(1-r));
      double v = param->y[person] - r;
      NumericVector z = param->z(person, _);
      for (i=0; i<p; i++) {
        score[i] += f*w*v*d*z[i];
      }
    }
  } else if (param->link == "cloglog") {
    for (person = 0; person < n; person++) {
      double f = param->freq[person];
      double w = param->weight[person];
      double r = 1 - exp(-exp(eta[person]));
      double phi = exp(eta[person] - exp(eta[person]));
      double d = phi/(r*(1-r));
      double v = param->y[person] - r;
      NumericVector z = param->z(person, _);
      for (i=0; i<p; i++) {
        score[i] += f*w*v*d*z[i];
      }
    }
  }

  return score;
}


// expected information matrix
NumericMatrix f_info_0(int p, NumericVector par, void *ex) {
  logparams *param = (logparams *) ex;
  int n = param->n;
  int person, i, j;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericMatrix imat(p,p);
  if (param->link == "logit") {
    for (person = 0; person < n; person++) {
      double f = param->freq[person];
      double w = param->weight[person];
      double v = R::dlogis(eta[person], 0, 1, 0);
      NumericVector z = param->z(person, _);
      for (i=0; i<p; i++) {
        for (j=0; j<=i; j++) {
          imat(i,j) += f*w*v*z[i]*z[j];
        }
      }
    }
  } else if (param->link == "probit") {
    for (person = 0; person < n; person++) {
      double f = param->freq[person];
      double w = param->weight[person];
      double r = R::pnorm(eta[person], 0, 1, 1, 0);
      double phi = R::dnorm(eta[person], 0, 1, 0);
      double v = phi*phi/(r*(1-r));
      NumericVector z = param->z(person, _);
      for (i=0; i<p; i++) {
        for (j=0; j<=i; j++) {
          imat(i,j) += f*w*v*z[i]*z[j];
        }
      }
    }
  } else if (param->link == "cloglog") {
    for (person = 0; person < n; person++) {
      double f = param->freq[person];
      double w = param->weight[person];
      double r = 1 - exp(-exp(eta[person]));
      double phi = exp(eta[person] - exp(eta[person]));
      double v = phi*phi/(r*(1-r));
      NumericVector z = param->z(person, _);
      for (i=0; i<p; i++) {
        for (j=0; j<=i; j++) {
          imat(i,j) += f*w*v*z[i]*z[j];
        }
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


// Firth's bias-reducing penalized log likelihood
double f_pen_llik_0(int p, NumericVector par, void *ex) {
  logparams *param = (logparams *) ex;
  double loglik = f_llik_0(p, par, param);
  NumericMatrix imat = f_info_0(p, par, param);

  // obtain the determinant of imat
  double toler = 1e-12;
  int i = cholesky2(imat, p, toler);

  double v = 0;
  for (i=0; i<p; i++) {
    v += log(imat(i,i));
  }

  return loglik + 0.5*v;
}


// Firth's penalized score vector
NumericVector f_pen_score_0(int p, NumericVector par, void *ex) {
  logparams *param = (logparams *) ex;
  int n = param->n;
  int person, i, j;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector pi(n), d(n), a(n), b(n);
  if (param->link == "logit") {
    for (person = 0; person < n; person++) {
      double r = R::plogis(eta[person], 0, 1, 1, 0);
      pi[person] = r;
      d[person] = 1;
      a[person] = r*(1-r);
      b[person] = 1-2*r;
    }
  } else if (param->link == "probit") {
    for (person = 0; person < n; person++) {
      double r = R::pnorm(eta[person], 0, 1, 1, 0);
      double phi = R::dnorm(eta[person], 0, 1, 0);
      double dphi = -eta[person];
      pi[person] = r;
      d[person] = phi/(r*(1-r));
      a[person] = phi*phi/(r*(1-r));
      b[person] = (2*r-1)*phi/(r*(1-r)) + 2*dphi;
    }
  } else if (param->link == "cloglog") {
    for (person = 0; person < n; person++) {
      double r = 1 - exp(-exp(eta[person]));
      double phi = exp(eta[person] - exp(eta[person]));
      double dphi = 1 - exp(eta[person]);
      pi[person] = r;
      d[person] = phi/(r*(1-r));
      a[person] = phi*phi/(r*(1-r));
      b[person] = (2*r-1)*phi/(r*(1-r)) + 2*dphi;
    }
  }

  NumericMatrix imat(p,p); // X^T W X
  for (person = 0; person < n; person++) {
    double f = param->freq[person];
    double w = param->weight[person];
    NumericVector z = param->z(person, _);
    for (i=0; i<p; i++) {
      for (j=0; j<=i; j++) {
        imat(i,j) += f*w*a[person]*z[i]*z[j];
      }
    }
  }

  for (i=0; i<p-1; i++) {
    for (j=i+1; j<p; j++) {
      imat(i,j) = imat(j,i);
    }
  }

  double toler = 1e-12;
  NumericMatrix var = invsympd(imat, p, toler);

  NumericVector score(p);
  for (person = 0; person < n; person++) {
    double f = param->freq[person];
    double w = param->weight[person];
    NumericVector z = param->z(person, _);
    double h = 0; // diagonal of H = W^{1/2}*X*inverse(X^T W X)*X^T W^{1/2}
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        h += var(i,j)*z[i]*z[j];
      }
    }
    h *= f*w*a[person];

    double resid = param->y[person] - pi[person];
    double u = f*w*resid*d[person] + 0.5*b[person]*h;
    for (i=0; i<p; i++) {
      score[i] += u*z[i];
    }
  }

  return score;
}


// score residual matrix (without firth, weight and frequency)
NumericMatrix f_ressco_0(int p, NumericVector par, void *ex) {
  logparams *param = (logparams *) ex;
  int n = param->n;
  int person, i;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericMatrix resid(n, p);
  if (param->link == "logit") {
    for (person = 0; person < n; person++) {
      double r = R::plogis(eta[person], 0, 1, 1, 0);
      double v = param->y[person] - r;
      NumericVector z = param->z(person, _);
      for (i=0; i<p; i++) {
        resid(person, i) = v*z[i];
      }
    }
  } else if (param->link == "probit") {
    for (person = 0; person < n; person++) {
      double r = R::pnorm(eta[person], 0, 1, 1, 0);
      double phi = R::dnorm(eta[person], 0, 1, 0);
      double d = phi/(r*(1-r));
      double v = param->y[person] - r;
      NumericVector z = param->z(person, _);
      for (i=0; i<p; i++) {
        resid(person, i) = v*d*z[i];
      }
    }
  } else if (param->link == "cloglog") {
    for (person = 0; person < n; person++) {
      double r = 1 - exp(-exp(eta[person]));
      double phi = exp(eta[person] - exp(eta[person]));
      double d = phi/(r*(1-r));
      double v = param->y[person] - r;
      NumericVector z = param->z(person, _);
      for (i=0; i<p; i++) {
        resid(person, i) = v*d*z[i];
      }
    }
  }

  return resid;
}


// underlying optimization algorithm for logisreg
//   colfit: vector of indices of parameters to update
//   ncolfit: number of parameters to update
List logisregloop(int p, NumericVector par, void *ex,
                  int maxiter, double eps, bool firth,
                  IntegerVector colfit, int ncolfit) {
  logparams *param = (logparams *) ex;
  int i, j, iter, halving = 0;
  bool fail;

  double toler = 1e-12;
  NumericVector beta(p), newbeta(p);
  double loglik, newlk = 0;
  NumericVector u(p);
  NumericMatrix imat(p,p);
  NumericVector u1(ncolfit);
  NumericMatrix imat1(ncolfit, ncolfit);


  // initial beta and log likelihood
  for (i=0; i<p; i++) {
    beta[i] = par[i];
  }

  if (!firth) {
    loglik = f_llik_0(p, beta, param);
    u = f_score_0(p, beta, param);
  } else {
    loglik = f_pen_llik_0(p, beta, param);
    u = f_pen_score_0(p, beta, param);
  }

  for (i=0; i<ncolfit; i++) {
    u1[i] = u[colfit[i]];
  }

  imat = f_info_0(p, beta, param);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1(i,j) = imat(colfit[i], colfit[j]);
    }
  }

  i = cholesky2(imat1, ncolfit, toler);

  chsolve2(imat1, ncolfit, u1);

  u.fill(0);
  for (i=0; i<ncolfit; i++) {
    u[colfit[i]] = u1[i];
  }

  for (i=0; i<p; i++) {
    newbeta[i] = beta[i] + u[i];
  }

  for (iter=0; iter<maxiter; iter++) {
    if (!firth) newlk = f_llik_0(p, newbeta, param);
    else newlk = f_pen_llik_0(p, newbeta, param);

    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk) == 1;

    if (!fail && halving == 0 && fabs(1 - (loglik/newlk)) < eps) {
      break;
    }

    if (fail || (newlk < loglik)) {
      // adjust step size if likelihood decreases
      halving++;
      for (i=0; i<p; i++) {
        newbeta[i] = (beta[i] + newbeta[i])/2;
      }
    } else { // update beta normally
      halving = 0;

      for (i=0; i<p; i++) {
        beta[i] = newbeta[i];
      }
      loglik = newlk;

      if (!firth) {
        u = f_score_0(p, beta, param);
      } else {
        u = f_pen_score_0(p, beta, param);
      }

      for (i=0; i<ncolfit; i++) {
        u1[i] = u[colfit[i]];
      }

      imat = f_info_0(p, beta, param);
      for (i=0; i<ncolfit; i++) {
        for (j=0; j<ncolfit; j++) {
          imat1(i,j) = imat(colfit[i], colfit[j]);
        }
      }

      i = cholesky2(imat1, ncolfit, toler);

      chsolve2(imat1, ncolfit, u1);

      u.fill(0);
      for (i=0; i<ncolfit; i++) {
        u[colfit[i]] = u1[i];
      }

      for (i=0; i<p; i++) {
        newbeta[i] = beta[i] + u[i];
      }
    }
  }

  if (iter == maxiter) fail = 1;

  imat = f_info_0(p, newbeta, param);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1(i,j) = imat(colfit[i], colfit[j]);
    }
  }

  NumericMatrix var1 = invsympd(imat1, ncolfit, toler);
  NumericMatrix var(p,p);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      var(colfit[i], colfit[j]) = var1(i,j);
    }
  }

  return List::create(
    Named("coef") = newbeta,
    Named("iter") = iter,
    Named("var") = var,
    Named("loglik") = newlk,
    Named("fail") = fail);
}


// confidence limit of profile likelihood method
//   k: index of the parameter to obtain PL CI
//   which: computation direction (-1 for lower limit or 1 for upper limit)
//   l0: boundary of profile loglik = max loglik - 0.5*qchisq(1-alpha, 1)
// refer to SAS PROC LOGISTIC documentation for Likelihood Ratio-Based
// Confidence Intervals for Parameters
double logisregplloop(int p, NumericVector par, void *ex,
                      int maxiter, double eps, bool firth,
                      int k, int which, double l0) {
  logparams *param = (logparams *) ex;

  int i, j, iter;
  bool fail = 0;

  NumericVector beta(p), newbeta(p);
  double loglik, newlk;
  NumericVector u(p);
  NumericVector delta(p);
  NumericMatrix imat(p,p);
  NumericMatrix v(p,p);

  // initial beta and log likelihood
  for (i=0; i<p; i++) {
    beta[i] = par[i];
  }

  if (!firth) {
    loglik = f_llik_0(p, beta, param);
    u = f_score_0(p, beta, param);
  } else {
    loglik = f_pen_llik_0(p, beta, param);
    u = f_pen_score_0(p, beta, param);
  }

  imat = f_info_0(p, beta, param);

  // Lagrange multiplier method as used in SAS PROC LOGISTIC
  double toler = 1e-12;
  v = invsympd(imat, p, toler);
  v = -1.0*v;

  double w = 0;
  for (i=0; i<p; i++) {
    for (j=0; j<p; j++) {
      w += u[i]*v(i,j)*u[j];
    }
  }

  double underroot = 2*(l0 - loglik + 0.5*w)/v(k,k);
  double lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
  u[k] += lambda;

  delta.fill(0.0);
  for (i=0; i<p; i++) {
    for (j=0; j<p; j++) {
      delta[i] -= v(i,j)*u[j];
    }
  }

  // update beta
  for (i=0; i<p; i++) {
    newbeta[i] = beta[i] + delta[i];
  }

  for (iter=0; iter<maxiter; iter++) {
    if (!firth) newlk = f_llik_0(p, newbeta, param);
    else newlk = f_pen_llik_0(p, newbeta, param);

    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk) == 1;

    if (!fail && fabs(newlk - l0) < eps && w < eps) {
      break;
    }

    for (i=0; i<p; i++) {
      beta[i] = newbeta[i];
    }
    loglik = newlk;

    if (!firth) {
      u = f_score_0(p, beta, param);
    } else {
      u = f_pen_score_0(p, beta, param);
    }

    imat = f_info_0(p, beta, param);

    // Lagrange multiplier method as used in SAS PROC LOGISTIC
    v = invsympd(imat, p, toler);
    v = -1.0*v;

    w = 0;
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        w += u[i]*v(i,j)*u[j];
      }
    }

    underroot = 2*(l0 - loglik + 0.5*w)/v(k,k);
    lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
    u[k] += lambda;

    delta.fill(0.0);
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        delta[i] -= v(i,j)*u[j];
      }
    }

    // update beta
    for (i=0; i<p; i++) {
      newbeta[i] = beta[i] + delta[i];
    }
  }

  if (iter == maxiter) fail = 1;

  if (fail) {
    warning("The algorithm in logisregplloop did not converge");
  }

  return newbeta[k];
}


// [[Rcpp::export]]
List logisregcpp(const DataFrame data,
                 const StringVector& rep = "",
                 const std::string event = "event",
                 const StringVector& covariates = "",
                 const std::string freq = "",
                 const std::string weight = "",
                 const std::string offset = "",
                 const std::string id = "",
                 const std::string link = "logit",
                 const NumericVector& init = NA_REAL,
                 const bool robust = 0,
                 const bool firth = 0,
                 const bool flic = 0,
                 const bool plci = 0,
                 const double alpha = 0.05,
                 const int maxiter = 50,
                 const double eps = 1.0e-9) {

  int h, i, j, k, n = data.nrows();
  int p = static_cast<int>(covariates.size()) + 1;
  if (p == 2 && (covariates[0] == "" || covariates[0] == "none")) {
    p = 1;
  }

  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(1);
  } else {
    List out = bygroup(data, rep);
    has_rep = 1;
    repn = out["index"];
    u_rep = DataFrame(out["lookup"]);
  }

  bool has_event = hasVariable(data, event);
  if (!has_event) {
    stop("data must contain the event variable");
  }
  NumericVector eventn(n);
  if (has_event) {
    NumericVector eventnz = data[event];
    eventn = clone(eventnz);
    if (is_true(any((eventn != 1) & (eventn != 0)))) {
      stop("event must be 1 or 0 for each subject");
    }

    if (is_true(all(eventn == 0))) {
      stop("at least 1 event is needed to fit the parametric model");
    }
  }

  NumericMatrix zn(n,p);
  for (i=0; i<n; i++) {
    zn(i,0) = 1; // intercept
  }
  for (j=0; j<p-1; j++) {
    String zj = covariates[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in covariates");
    }
    NumericVector u = data[zj];
    for (i=0; i<n; i++) {
      zn(i,j+1) = u[i];
    }
  }

  bool has_freq = hasVariable(data, freq);
  NumericVector freqn(n, 1.0);
  if (has_freq) {
    NumericVector freqnz = data[freq];
    freqn = clone(freqnz);
    if (is_true(any((freqn <= 0) | (freqn != floor(freqn))))) {
      stop("freq must be positive integers");
    }
  }

  bool has_weight = hasVariable(data, weight);
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }

  bool has_offset = hasVariable(data, offset);
  NumericVector offsetn(n);
  if (has_offset) {
    NumericVector offsetnz = data[offset];
    offsetn = clone(offsetnz);
  }

  bool has_id = hasVariable(data, id);
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

  std::string link1 = link;
  std::for_each(link1.begin(), link1.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  if (!(link1=="logit" || link1=="probit" || link1=="cloglog")) {
    stop("Invalid value for link");
  }

  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  repn = repn[order];
  eventn = eventn[order];
  freqn = freqn[order];
  weightn = weightn[order];
  offsetn = offsetn[order];
  idn = idn[order];
  zn = subset_matrix_by_row(zn, order);

  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (i=0; i<n; i++) {
    if ((repn[i] == NA_INTEGER) || (eventn[i] == NA_INTEGER) ||
        (std::isnan(freqn[i])) || (std::isnan(weightn[i])) ||
        (std::isnan(offsetn[i])) || (idn[i] == NA_INTEGER)) {
      sub[i] = 0;
    }
    for (j=0; j<p-1; j++) {
      if (std::isnan(zn(i,j+1))) sub[i] = 0;
    }
  }

  order = which(sub);
  repn = repn[order];
  eventn = eventn[order];
  freqn = freqn[order];
  weightn = weightn[order];
  offsetn = offsetn[order];
  idn = idn[order];
  zn = subset_matrix_by_row(zn, order);
  n = sum(sub); // number of nonmissing observations

  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (i=1; i<n; i++) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }

  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);

  // variables in the output data sets
  IntegerVector rep01 = seq(1,nreps);
  NumericVector nobs(nreps), nevents(nreps);
  NumericVector loglik0(nreps), loglik1(nreps);
  NumericVector regloglik0(nreps), regloglik1(nreps); // regular loglik
  IntegerVector niter(nreps);
  LogicalVector fails(nreps);

  IntegerVector rep0(nreps*p);
  StringVector par0(nreps*p);
  NumericVector beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p);
  NumericMatrix vbeta0(nreps*p,p), rvbeta0(nreps*p,p);
  NumericVector lb0(nreps*p), ub0(nreps*p), prob0(nreps*p);
  StringVector clparm0(nreps*p);

  NumericVector linear_predictors(n), fitted_values(n);

  int N=0; // offset for the current replication data set
  for (h=0; h<nreps; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    NumericVector event1 = eventn[q1];
    NumericVector freq1 = freqn[q1];
    NumericVector weight1 = weightn[q1];
    NumericVector offset1 = offsetn[q1];
    IntegerVector id1 = idn[q1];
    NumericMatrix z1 = subset_matrix_by_row(zn, q1);

    // number of trials and number of events accounting for frequencies
    nobs[h] = sum(freq1);
    nevents[h] = sum(freq1*event1);

    // intercept only model
    double num = 0, den = 0;
    for (i=0; i<n1; i++) {
      num += freq1[i]*weight1[i]*event1[i];
      den += freq1[i]*weight1[i];
    }
    if (firth) {
      num += 0.5;
      den += 1.0;
    }

    NumericVector bint0(p);
    bint0[0] = R::qlogis(num/den, 0, 1, 1, 0);

    IntegerVector colfit0(1);
    logparams param = {n1, link1, event1, z1, freq1, weight1, offset1};
    List outint = logisregloop(p, bint0, &param, maxiter, eps, firth,
                               colfit0, 1);

    NumericVector bint = outint["coef"];
    NumericMatrix vbint = as<NumericMatrix>(outint["var"]);

    NumericVector b(p);
    NumericMatrix vb(p,p);
    List out;

    if (p > 1) {
      // parameter estimates and standard errors for the full model
      IntegerVector colfit = seq(0,p-1);
      if (is_false(any(is_na(init))) && init.size() == p) {
        out = logisregloop(p, init, &param, maxiter, eps, firth, colfit, p);
      } else {
        out = logisregloop(p, bint, &param, maxiter, eps, firth, colfit, p);
      }

      bool fail = out["fail"];
      if (fail) warning("The algorithm in logisregr did not converge");

      b = out["coef"];
      vb = as<NumericMatrix>(out["var"]);

      // intercept correction
      if (flic) {
        NumericVector lp(n1);  // linear predictor excluding intercept
        for (int person = 0; person < n1; person++) {
          lp[person] = offset1[person];
          for (i=1; i<p; i++) {
            lp[person] += b[i]*z1(person,i);
          }
        }

        logparams param0 = {n1, link1, event1, z1, freq1, weight1, lp};
        NumericVector bint00(1, bint0[0]);
        outint = logisregloop(1, bint00, &param0, maxiter, eps, 0, 
                              colfit0, 1);
        double a = as<double>(outint["coef"]);
        double va = as<double>(outint["var"]);

        // update the intercept estimate
        b[0] = a;

        // partial derivative of alpha(beta) with respect to beta
        NumericMatrix iflic = f_info_0(p, b, &param);
        NumericVector der(p-1);
        for (i=0; i<p-1; i++) {
          der[i] = -iflic(i+1,0)/iflic(0,0);
        }

        // update the variance of alpha
        vb(0,0) = va;
        for (i=0; i<p-1; i++) {
          for (j=0; j<p-1; j++) {
            vb(0,0) += der[i]*vb(i+1,j+1)*der[j];
          }
        }

        // update the covariance between alpha and beta
        for (i=0; i<p-1; i++) {
          vb(i+1,0) = 0;
          for (j=0; j<p-1; j++) {
            vb(i+1,0) += vb(i+1,j+1)*der[j];
          }
          vb(0,i+1) = vb(i+1,0);
        }
      }
    } else {
      b = bint;
      vb = vbint;
      out = outint;

      if (flic) {
        out = logisregloop(p, bint0, &param, maxiter, eps, 0, colfit0, 1);
        b = out["coef"];
        vb = as<NumericMatrix>(outint["var"]);
      }
    }

    NumericVector seb(p);
    for (j=0; j<p; j++) {
      seb[j] = sqrt(vb(j,j));
    }

    for (i=0; i<p; i++) {
      rep0[h*p+i] = h+1;
      if (i==0) {
        par0[h*p+i] = "(Intercept)";
      } else {
        par0[h*p+i] = covariates[i-1];
      }
      beta0[h*p+i] = b[i];
      sebeta0[h*p+i] = seb[i];
      for (j=0; j<p; j++) {
        vbeta0(h*p+i,j) = vb(i,j);
      }
    }

    // linear predictors and fitted values
    int person;
    NumericVector eta(n1);
    for (person = 0; person < n1; person++) {
      eta[person] = offset1[person];
      for (i=0; i<p; i++) {
        eta[person] += b[i]*z1(person,i);
      }
      linear_predictors[N+person] = eta[person];
    }

    if (link1 == "logit") {
      for (person = 0; person < n1; person++) {
        fitted_values[N+person] = R::plogis(eta[person], 0, 1, 1, 0);
      }
    } else if (link1 == "probit") {
      for (person = 0; person < n1; person++) {
        fitted_values[N+person] = R::pnorm(eta[person], 0, 1, 1, 0);
      }
    } else if (link1 == "cloglog") {
      for (person = 0; person < n1; person++) {
        fitted_values[N+person] = 1 - exp(-exp(eta[person]));
      }
    }

    N += n1;

    niter[h] = out["iter"];
    fails[h] = out["fail"];

    // robust variance estimates
    NumericVector rseb(p);  // robust standard error for betahat
    if (robust) {
      NumericMatrix ressco = f_ressco_0(p, b, &param);

      int nr; // number of rows in the score residual matrix
      NumericVector freqr;
      if (!has_id) {
        for (i=0; i<n1; i++) {
          for (j=0; j<p; j++) {
            ressco(i,j) = weight1[i]*ressco(i,j);
          }
        }
        nr = n1;
        freqr = clone(freq1);
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

        int nids = static_cast<int>(idx.size());
        idx.push_back(n1);

        NumericMatrix resid(n1,p);
        for (i=0; i<n1; i++) {
          for (j=0; j<p; j++) {
            resid(i,j) = ressco(order[i],j);
          }
        }

        NumericVector weight2 = weight1[order];
        NumericVector freq2 = freq1[order];
        NumericVector freqr0(nids); // cluster frequency

        NumericMatrix ressco2(nids,p);
        for (i=0; i<nids; i++) {
          for (j=0; j<p; j++) {
            for (k=idx[i]; k<idx[i+1]; k++) {
              ressco2(i,j) += weight2[k]*resid(k,j);
            }
          }
          freqr0[i] = freq2[idx[i]];
        }

        ressco = ressco2;  // update the score residuals
        nr = nids;
        freqr = freqr0;
      }

      NumericMatrix D(nr,p); // DFBETA
      for (i=0; i<nr; i++) {
        for (j=0; j<p; j++) {
          for (k=0; k<p; k++) {
            D(i,j) += ressco(i,k)*vb(k,j);
          }
        }
      }

      NumericMatrix rvb(p,p); // robust variance matrix for betahat
      for (j=0; j<p; j++) {
        for (k=0; k<p; k++) {
          for (i=0; i<nr; i++) {
            rvb(j,k) += freqr[i]*D(i,j)*D(i,k);
          }
        }
      }

      for (i=0; i<p; i++) {
        rseb[i] = sqrt(rvb(i,i));
      }

      for (i=0; i<p; i++) {
        rsebeta0[h*p+i] = rseb[i];
        for (j=0; j<p; j++) {
          rvbeta0(h*p+i,j) = rvb(i,j);
        }
      }
    }

    // profile likelihood confidence interval for regression coefficients
    NumericVector lb(p), ub(p), prob(p);
    StringVector clparm(p);

    double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
    if (plci) {
      double lmax;
      if (firth) {
        lmax = f_pen_llik_0(p, b, &param);
      } else {
        lmax = f_llik_0(p, b, &param);
      }
      double l0 = lmax - 0.5*R::qchisq(1-alpha, 1, 1, 0);

      if (!(firth && flic)) {
        for (k=0; k<p; k++) {
          lb[k] = logisregplloop(p, b, &param, maxiter, eps, firth, 
                                 k, -1, l0);
          ub[k] = logisregplloop(p, b, &param, maxiter, eps, firth, 
                                 k, 1, l0);

          IntegerVector colfit1(p-1);
          for (i=0; i<k; i++) {
            colfit1[i] = i;
          }
          for (i=k+1; i<p; i++) {
            colfit1[i-1] = i;
          }

          NumericVector b0(p);
          List out0 = logisregloop(p, b0, &param, maxiter, eps, firth,
                                   colfit1, p-1);
          double lmax0 = out0["loglik"];
          prob[k] = R::pchisq(-2*(lmax0 - lmax), 1, 0, 0);
          clparm[k] = "PL";
        }
      } else { // Wald CI for intercept and PL CI for slopes
        if (!robust) {
          lb[0] = b[0] - zcrit*seb[0];
          ub[0] = b[0] + zcrit*seb[0];
          prob[0] = R::pchisq(pow(b[0]/seb[0], 2), 1, 0, 0);
        } else {
          lb[0] = b[0] - zcrit*rseb[0];
          ub[0] = b[0] + zcrit*rseb[0];
          prob[0] = R::pchisq(pow(b[0]/rseb[0], 2), 1, 0, 0);
        }
        clparm[0] = "Wald";

        for (k=1; k<p; k++) {
          lb[k] = logisregplloop(p, b, &param, maxiter, eps, firth, 
                                 k, -1, l0);
          ub[k] = logisregplloop(p, b, &param, maxiter, eps, firth, 
                                 k, 1, l0);

          IntegerVector colfit1(p-1);
          for (i=0; i<k; i++) {
            colfit1[i] = i;
          }
          for (i=k+1; i<p; i++) {
            colfit1[i-1] = i;
          }

          NumericVector b0(p);
          List out0 = logisregloop(p, b0, &param, maxiter, eps, firth, 
                                   colfit1, p-1);
          double lmax0 = out0["loglik"];
          prob[k] = R::pchisq(-2*(lmax0 - lmax), 1, 0, 0);
          clparm[k] = "PL";
        }
      }
    } else {
      for (k=0; k<p; k++) {
        if (!robust) {
          lb[k] = b[k] - zcrit*seb[k];
          ub[k] = b[k] + zcrit*seb[k];
          prob[k] = R::pchisq(pow(b[k]/seb[k], 2), 1, 0, 0);
        } else {
          lb[k] = b[k] - zcrit*rseb[k];
          ub[k] = b[k] + zcrit*rseb[k];
          prob[k] = R::pchisq(pow(b[k]/rseb[k], 2), 1, 0, 0);
        }
        clparm[k] = "Wald";
      }
    }

    for (i=0; i<p; i++) {
      lb0[h*p+i] = lb[i];
      ub0[h*p+i] = ub[i];
      prob0[h*p+i] = prob[i];
      clparm0[h*p+i] = clparm[i];
    }

    // log-likelihoods
    if (firth) {
      loglik0[h] = f_pen_llik_0(p, bint, &param);
      loglik1[h] = f_pen_llik_0(p, b, &param);
      regloglik0[h] = f_llik_0(p, bint, &param);
      regloglik1[h] = f_llik_0(p, b, &param);
    } else {
      loglik0[h] = f_llik_0(p, bint, &param);
      loglik1[h] = f_llik_0(p, b, &param);
    }
  }

  NumericVector expbeta0 = exp(beta0);
  NumericVector z0(nreps*p);
  if (!robust) z0 = beta0/sebeta0;
  else z0 = beta0/rsebeta0;

  DataFrame sumstat = List::create(
    _["n"] = nobs,
    _["nevents"] = nevents,
    _["loglik0"] = loglik0,
    _["loglik1"] = loglik1,
    _["niter"] = niter,
    _["p"] = p,
    _["link"] = link1,
    _["robust"] = robust,
    _["firth"] = firth,
    _["flic"] = flic,
    _["fail"] = fails);

  if (firth) {
    sumstat.push_back(regloglik0, "loglik0_unpenalized");
    sumstat.push_back(regloglik1, "loglik1_unpenalized");
  }

  DataFrame parest;
  if (!robust) {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = sebeta0,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
      _["vbeta"] = vbeta0,
      _["lower"] = lb0,
      _["upper"] = ub0,
      _["p"] = prob0,
      _["method"] = clparm0);
  } else {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = rsebeta0,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
      _["vbeta"] = rvbeta0,
      _["lower"] = lb0,
      _["upper"] = ub0,
      _["p"] = prob0,
      _["method"] = clparm0,
      _["sebeta_naive"] = sebeta0,
      _["vbeta_naive"] = vbeta0);
  }

  DataFrame fitted = DataFrame::create(
    Named("linear_predictors") = linear_predictors,
    Named("fitted_values") = fitted_values);

  if (has_rep) {
    for (i=0; i<p_rep; i++) {
      String s = rep[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector repwi = u_rep[s];
        sumstat.push_back(repwi[rep01-1], s);
        parest.push_back(repwi[rep0-1], s);
        fitted.push_back(repwi[repn-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector repwn = u_rep[s];
        sumstat.push_back(repwn[rep01-1], s);
        parest.push_back(repwn[rep0-1], s);
        fitted.push_back(repwn[repn-1], s);
      } else if (TYPEOF(data[rep]) == STRSXP) {
        StringVector repwc = u_rep[s];
        sumstat.push_back(repwc[rep01-1], s);
        parest.push_back(repwc[rep0-1], s);
        fitted.push_back(repwc[repn-1], s);
      }
    }
  }

  List result = List::create(
    _["sumstat"] = sumstat,
    _["parest"] = parest,
    _["fitted"] = fitted);

  return result;
}
