#include "logistic_regression.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <Rcpp.h>

#include <algorithm> // any_of, fill, for_each, none_of, sort
#include <cctype>    // tolower
#include <climits>   // INT_MIN
#include <cmath>     // exp, fabs, isinf, isnan, log, pow
#include <numeric>   // accumulate, inner_product, iota
#include <stdexcept> // invalid_argument, runtime_error
#include <string>    // string
#include <utility>   // move
#include <vector>    // vector

// structure to hold parameters for logistic regression
struct logparams {
  int n;
  int link_code; // 1: logit, 2: probit, 3: cloglog
  std::vector<int> y;
  FlatMatrix z; // n x p column-major
  std::vector<int> freq;
  std::vector<double> weight;
  std::vector<double> offset;
};

// --------------------------- f_der_0 (log-likelihood, score, information) ----
ListCpp f_der_0(int p, const std::vector<double>& par, void *ex, bool firth) {
  logparams *param = (logparams *) ex;
  const int n = param->n;
  const int link_code = param->link_code;
  const std::vector<int>& yv = param->y;
  const double* zptr = param->z.data_ptr();
  const std::vector<int>& freq = param->freq;
  const std::vector<double>& weight = param->weight;
  std::vector<double> fwvec(n);  // freq * weight per observation
  for (int i = 0; i < n; ++i) {
    fwvec[i] = freq[i] * weight[i];
  }
  
  // compute linear predictor eta efficiently using column-major storage:
  std::vector<double> eta = param->offset; // initialize with offset
  // add contributions of each coefficient times column
  for (int i = 0; i < p; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const double* zi = zptr + i * n;
    for (int r = 0; r < n; ++r) {
      eta[r] += beta * zi[r];
    }
  }
  
  double loglik = 0.0;
  std::vector<double> score(p);
  FlatMatrix imat(p, p); // information matrix in column-major p x p
  
  // Pre-allocate per-observation temporaries
  std::vector<double> rvec(n); // fitted probabilities
  std::vector<double> c1(n);   // contribution for score
  std::vector<double> c2(n);   // contribution for information

  // firth temporaries
  std::vector<double> pi, d, a, b;
  if (firth) {
    pi.assign(n, 0.0);
    d.assign(n, 0.0);
    a.assign(n, 0.0);
    b.assign(n, 0.0);
  }
  
  // 1) single pass over observations: compute r, loglik, c1, c2 + firth temporaries
  for (int person = 0; person < n; ++person) {
    double fw = fwvec[person];
    double y = yv[person];
    double et = eta[person];
    
    if (link_code == 1) {                // logit
      double r = boost_plogis(et);
      rvec[person] = r;
      loglik += fw * (y * et + std::log(1.0 - r));
      c1[person] = fw * (y - r);
      c2[person] = fw * r * (1.0 - r);
      if (firth) {
        pi[person] = r;
        d[person]  = 1.0;
        a[person]  = r * (1.0 - r);
        b[person]  = 1.0 - 2.0 * r;
      }
    } else if (link_code == 2) {         // probit
      double r = boost_pnorm(et);
      double phi = boost_dnorm(et);
      rvec[person] = r;
      loglik += fw * (y * std::log(r / (1.0 - r)) + std::log(1.0 - r));
      c1[person] = fw * (y - r) * (phi / (r * (1.0 - r)));
      c2[person] = fw * (phi * phi / (r * (1.0 - r)));
      if (firth) {
        double dphi = -et;
        pi[person] = r;
        d[person]  = phi / (r * (1.0 - r));
        a[person]  = phi * phi / (r * (1.0 - r));
        b[person]  = (2.0 * r - 1.0) * phi / (r * (1.0 - r)) + 2.0 * dphi;
      }
    } else {                                    // cloglog / extreme
      double r = boost_pextreme(et);
      double phi = boost_dextreme(et);
      rvec[person] = r;
      loglik += fw * (y * std::log(r / (1.0 - r)) + std::log(1.0 - r));
      c1[person] = fw * (y - r) * (phi / (r * (1.0 - r)));
      c2[person] = fw * (phi * phi / (r * (1.0 - r)));
      if (firth) {
        double dphi = 1 - std::exp(et);
        pi[person] = r;
        d[person]  = phi / (r * (1.0 - r));
        a[person]  = phi * phi / (r * (1.0 - r));
        b[person]  = (2.0 * r - 1.0) * phi / (r * (1.0 - r)) + 2.0 * dphi;
      }
    }
  }
  
  // 2) compute score vector by column dot-products: score[j] = sum_i c1[i] * z_{i,j}
  for (int j = 0; j < p; ++j) {
    const double* zcol = zptr + j * n;   // start of column j
    double s = 0.0;
    // inner loop is contiguous load of zcol[i]
    for (int i = 0; i < n; ++i) s += c1[i] * zcol[i];
    score[j] = s;
  }
  
  // 3) compute information matrix (lower triangle) using c2 as per-observation weight
  for (int j = 0; j < p; ++j) {
    const double* zj = zptr + j * n;
    for (int i = j; i < p; ++i) {
    const double* zi = zptr + i * n;
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += c2[k] * zi[k] * zj[k];
      }
      imat(i, j) = sum;
    }
  }
  
  // 4) fill upper triangle as before (keep symmetric)
  for (int j = 1; j < p; ++j) {
    for (int i = 0; i < j; ++i) {
      imat(i, j) = imat(j, i);
    }
  }
  
  if (firth) {
    // make a copy for Cholesky so we preserve imat for output
    FlatMatrix imat0 = imat;
    cholesky2(imat0, p); // in-place Cholesky on imat0
    
    double sumlog = 0.0;
    for (int i = 0; i < p; ++i) sumlog += std::log(imat0(i, i));
    
    double penloglik = loglik + 0.5 * sumlog;
    
    // precompute per-person scalar c[k] = f * w * a[k]
    std::vector<double> c(n);
    for (int person = 0; person < n; ++person) {
      c[person] = fwvec[person] * a[person];
    }
    
    FlatMatrix xwx(p, p); // data initially zeroed by constructor
    
    // compute lower triangle (i >= j) using column-major access
    for (int j = 0; j < p; ++j) {
      const double* zj = zptr + j * n;      // Z(:, j)
      for (int i = j; i < p; ++i) {
      const double* zi = zptr + i * n;        // Z(:, i)
        double sum = 0.0;
        // inner loop reads zi[k] and zj[k] contiguously
        for (int k = 0; k < n; ++k) {
          sum += c[k] * zi[k] * zj[k];
        }
        xwx(i, j) = sum;
      }
    }
    
    // fill upper triangle of xwx
    for (int j = 1; j < p; ++j) {
      for (int i = 0; i < j; ++i) {
        xwx(i, j) = xwx(j, i);
      }
    }
    
    // compute inverse of xwx as FlatMatrix
    FlatMatrix var = invsympd(xwx, p);
    const double* vptr = var.data_ptr();
    
    std::vector<double> h0(n, 0.0);
    for (int k = 0; k < p; ++k) {
      const double* zk = zptr + k * n;   // column k of Z
      const double* vk = vptr + k * p;   // column k of var
      for (int j = 0; j < p; ++j) {
        const double* zj = zptr + j * n; // column j of Z
        const double v = vk[j];          // var(j, k)
        for (int i = 0; i < n; ++i) {
          h0[i] += v * zj[i] * zk[i];
        }
      }
    }
    
    // compute u[r] = f * w * resid * d + 0.5 * b * (f * w * a * h0)
    std::vector<double> u(n);
    for (int i = 0; i < n; ++i) {
      double fw = fwvec[i];
      double resid = yv[i] - pi[i];
      double h_scaled = h0[i] * fw * a[i];
      u[i] = fw * resid * d[i] + 0.5 * b[i] * h_scaled;
    }
    
    // compute g = Z^T * u
    std::vector<double> g(p, 0.0); 
    for (int j = 0; j < p; ++j) { 
      const double* zj = zptr + j * n; // Z(:, j) 
      double sum = 0.0; 
      for (int i = 0; i < n; ++i) sum += u[i] * zj[i]; 
      g[j] = sum;
    }
    
    ListCpp result;
    result.push_back(penloglik, "loglik");
    result.push_back(std::move(g), "score");
    result.push_back(std::move(imat), "imat");
    result.push_back(loglik, "regloglik");
    result.push_back(std::move(score), "regscore");
    return result;
  } else {
    ListCpp result;
    result.push_back(loglik, "loglik");
    result.push_back(std::move(score), "score");
    result.push_back(std::move(imat), "imat");       // FlatMatrix
    return result;
  }
}


// --------------------------- f_ressco_0 (score residuals) --------------------
// Returns an n x p FlatMatrix (column-major), where entry (r, c) equals 
// residual for observation r and covariate c.
FlatMatrix f_ressco_0(int p, const std::vector<double>& par, void *ex) {
  logparams *param = (logparams *) ex;
  const int n = param->n;
  const int link_code = param->link_code;
  const std::vector<int>& yv = param->y;
  const double* zptr = param->z.data_ptr();
  
  // compute eta similarly to f_der_0
  std::vector<double> eta = param->offset;
  for (int i = 0; i < p; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const double* zi = zptr + i * n;
    for (int r = 0; r < n; ++r) {
      eta[r] += beta * zi[r];
    }
  }
  
  FlatMatrix resid(n, p);
  double* rptr = resid.data_ptr(); 
  
  switch (link_code) {
  case 1: { // logit
    std::vector<double> v(n);
    for (int person = 0; person < n; ++person) {
      double r = boost_plogis(eta[person]);
      v[person] = yv[person] - r;
    }
    
    for (int i = 0; i < p; ++i) {
      int off = i * n;
      const double* zcol = zptr + off;
      double* rcol = rptr + off;
      for (int person = 0; person < n; ++person) {
        rcol[person] = v[person] * zcol[person];
      }
    }
    break;
  }
    
  case 2: { // probit
    std::vector<double> vd(n);
    for (int person = 0; person < n; ++person) {
      double r = boost_pnorm(eta[person]);
      double phi = boost_dnorm(eta[person]);
      double d = phi / (r * (1.0 - r));
      double v = yv[person] - r;
      vd[person] = v * d;
    }
    
    for (int i = 0; i < p; ++i) {
      int off = i * n;
      const double* zcol = zptr + off;
      double* rcol = rptr + off;
      for (int person = 0; person < n; ++person) {
        rcol[person] = vd[person] * zcol[person];
      }
    }
    break;
  }
  case 3: {// cloglog / extreme
    std::vector<double> vd(n);
    for (int person = 0; person < n; ++person) {
      double r = boost_pextreme(eta[person]);
      double phi = boost_dextreme(eta[person]);
      double d = phi / (r * (1.0 - r));
      double v = yv[person] - r;
      vd[person] = v * d;
    }
    
    for (int i = 0; i < p; ++i) {
      int off = i * n;
      const double* zcol = zptr + off;
      double* rcol = rptr + off;
      for (int person = 0; person < n; ++person) {
        rcol[person] = vd[person] * zcol[person];
      }
    }
    break;
  }
  }
  return resid;
}

// --------------------------- logisregloop, logisregplloop and logisregcpp ----
ListCpp logisregloop(int p, const std::vector<double>& par, void *ex,
                     int maxiter, double eps, bool firth,
                     const std::vector<int>& colfit, int ncolfit) {
  logparams *param = (logparams *) ex;
  
  int iter = 0, halving = 0;
  bool fail = false;
  
  std::vector<double> beta = par;
  std::vector<double> newbeta(p);
  double loglik = 0.0, newlk = 0.0;
  std::vector<double> u(p);
  std::vector<double> uu1(ncolfit);
  double* u1 = uu1.data();
  FlatMatrix imat(p, p);
  FlatMatrix imat1(ncolfit, ncolfit);
  
  // --- first step ---
  ListCpp der = f_der_0(p, beta, param, firth);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<FlatMatrix>("imat");
  
  for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
  
  for (int j = 0; j < ncolfit; ++j)
    for (int i = 0; i < ncolfit; ++i)
      imat1(i,j) = imat(colfit[i], colfit[j]);
  
  cholesky2(imat1, ncolfit);
  chsolve2(imat1, ncolfit, u1);
  
  std::fill(u.begin(), u.end(), 0.0);
  for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_0(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && std::fabs(1 - (loglik / newlk)) < eps) break;
    
    if (fail || newlk < loglik) {
      ++halving;
      for (int i = 0; i < p; ++i) {
        newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
      }
      continue;
    }
    
    halving = 0;
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    
    for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
    
    for (int j = 0; j < ncolfit; ++j)
      for (int i = 0; i < ncolfit; ++i)
        imat1(i, j) = imat(colfit[i], colfit[j]);
    
    cholesky2(imat1, ncolfit);
    chsolve2(imat1, ncolfit, u1);
    
    std::fill(u.begin(), u.end(), 0.0);
    for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  }
  
  if (iter == maxiter) fail = true;
  
  // final variance assembly
  imat = der.get<FlatMatrix>("imat");
  
  for (int j = 0; j < ncolfit; ++j)
    for (int i = 0; i < ncolfit; ++i)
      imat1(i, j) = imat(colfit[i], colfit[j]);
  
  FlatMatrix var1 = invsympd(imat1, ncolfit);
  FlatMatrix var(p, p);
  for (int j = 0; j < ncolfit; ++j)
    for (int i = 0; i < ncolfit; ++i)
      var(colfit[i], colfit[j]) = var1(i, j);
  
  ListCpp result;
  result.push_back(std::move(newbeta), "coef");
  result.push_back(iter, "iter");
  result.push_back(std::move(var), "var");
  result.push_back(newlk, "loglik");
  result.push_back(fail, "fail");
  
  if (firth) {
    double regloglik = der.get<double>("regloglik");
    result.push_back(regloglik, "regloglik");
  }
  
  return result;
}

// --------------------------- logisregplloop (profile likelihood solver) -----
double logisregplloop(int p, const std::vector<double>& par,
                      void *ex, int maxiter, double eps, bool firth,
                      int k, int direction, double l0) {
  logparams *param = (logparams *) ex;
  int iter;
  bool fail = false;
  
  std::vector<double> beta = par;
  std::vector<double> newbeta(p);
  double loglik = 0.0, newlk = 0.0;
  std::vector<double> u(p);
  std::vector<double> delta(p);
  FlatMatrix imat(p, p);
  FlatMatrix v(p, p);
  
  ListCpp der = f_der_0(p, beta, param, firth);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<FlatMatrix>("imat");
  v = invsympd(imat, p);
  
  // compute w = - u^T v u
  double w = -quadsym(u, v);
  double underroot = -2 * (l0 - loglik + 0.5 * w) / v(k, k);
  double lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
  u[k] += lambda;
  delta = mat_vec_mult(v, u);
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  
  // iterate to convergence 
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_0(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && std::fabs(newlk - l0) < eps && w < eps) break;
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    v = invsympd(imat, p);
    w = -quadsym(u, v);
    underroot = -2 * (l0 - newlk + 0.5 * w) / v(k, k);
    lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
    u[k] += lambda;
    delta = mat_vec_mult(v, u);
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  }
  
  if (iter == maxiter) fail = true;
  if (fail) thread_utils::push_thread_warning("logisregplloop did not converge.");
  
  return newbeta[k];
}

// --------------------------- logisregcpp (high-level API) --------------------
// Convert inputs to FlatMatrix design matrix and use new functions where appropriate.
ListCpp logisregcpp(const DataFrameCpp& data,
                    const std::string& event,
                    const std::vector<std::string>& covariates,
                    const std::string& freq,
                    const std::string& weight,
                    const std::string& offset,
                    const std::string& id,
                    const std::string& link,
                    const std::vector<double>& init,
                    const bool robust,
                    const bool firth,
                    const bool flic,
                    const bool plci,
                    const double alpha,
                    const int maxiter,
                    const double eps) {
  
  int n = data.nrows();
  int p = covariates.size() + 1;
  if (p == 2 && covariates[0].empty()) p = 1;
  
  if (event.empty()) throw std::invalid_argument("event variable is not specified");
  if (!data.containElementNamed(event)) 
    throw std::invalid_argument("data must contain the event variable");
  std::vector<int> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
  } else if (data.int_cols.count(event)) {
    eventn = data.get<int>(event);
  } else if (data.numeric_cols.count(event)) {
    const std::vector<double>& vd = data.get<double>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1) 
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  // construct design matrix zn (n x p) as FlatMatrix column-major
  FlatMatrix zn(n, p);
  // intercept column
  for (int i = 0; i < n; ++i) zn.data[i] = 1.0;
  
  // fill covariate columns (1..p-1)
  for (int j = 0; j < p - 1; ++j) {
    const std::string& zj = covariates[j];
    if (!data.containElementNamed(zj)) 
      throw std::invalid_argument("data must contain the variables in covariates");
    double* zn_col = zn.data_ptr() + (j + 1) * n;
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      std::memcpy(zn_col, vd.data(), n * sizeof(double));
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }
  
  // freq, weight, offset
  std::vector<int> freqn(n, 1);
  if (!freq.empty() && data.containElementNamed(freq)) {
    if (data.int_cols.count(freq)) {
      freqn = data.get<int>(freq);
    } else if (data.numeric_cols.count(freq)) {
      const std::vector<double>& vd = data.get<double>(freq);
      for (int i = 0; i < n; ++i) freqn[i] = static_cast<int>(vd[i]);
    } else throw std::invalid_argument("freq variable must be integer or numeric");
    for (double v : freqn) if (v <= 0) 
      throw std::invalid_argument("freq must be positive integers");
  }
  
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const auto& weighti = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(weighti[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else throw std::invalid_argument("weight variable must be integer or numeric");
    for (double v : weightn) if (v <= 0.0) 
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    if (data.int_cols.count(offset)) {
      const auto& offseti = data.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(offseti[i]);
    } else if (data.numeric_cols.count(offset)) {
      offsetn = data.get<double>(offset);
    } else throw std::invalid_argument("offset variable must be integer or numeric");
  }
  
  // id processing (unchanged semantics)
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(n);
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      auto v = data.get<int>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      auto v = data.get<double>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      auto v = data.get<std::string>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else {
      throw std::invalid_argument(
          "incorrect type for the id variable in the input data");
    }
  }
  
  // link code mapping
  std::string link1 = link;
  for (char &c : link1) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  
  if (link1 == "log-log" || link1 == "loglog") link1 = "cloglog";
  
  int link_code = 0;
  if (link1 == "logit") link_code = 1;
  else if (link1 == "probit") link_code = 2;
  else if (link1 == "cloglog") link_code = 3;
  else throw std::invalid_argument("invalid link: " + link1);
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n, 1);
  for (int i = 0; i < n; ++i) { 
    if (eventn[i] == INT_MIN || freqn[i] == INT_MIN || 
        std::isnan(weightn[i]) || std::isnan(offsetn[i]) || 
        idn[i] == INT_MIN) {
      sub[i] = 0;
      continue;
    }
    for (int j = 0; j < p - 1; ++j) {
      if (std::isnan(zn(i, j+1))) { sub[i] = 0; break; }
    }
  }
  
  std::vector<int> order = which(sub);
  subset_in_place(eventn, order);
  subset_in_place(freqn, order);
  subset_in_place(weightn, order);
  subset_in_place(offsetn, order);
  subset_in_place(idn, order);
  subset_in_place_flatmatrix(zn, order);
  n = order.size();
  if (n == 0) throw std::invalid_argument("no observations without missing values");
  
  // sumstat data set
  int nobs, nevents;
  double loglik0, loglik1;
  double regloglik0, regloglik1; // regular loglik
  int niter;
  bool fail;
  
  // parest data set
  std::vector<std::string> par(p);
  std::vector<double> b(p), seb(p), rseb(p);
  std::vector<double> z(p), expbeta(p);
  FlatMatrix vb(p, p), rvb(p, p);
  std::vector<double> lb(p), ub(p), prob(p);
  std::vector<std::string> clparm(p);
  
  // linear predictor and fitted values for all observations
  std::vector<double> linear_predictors(n), fitted_values(n);
  
  double zcrit = boost_qnorm(1.0 - alpha / 2.0);
  double xcrit = zcrit * zcrit;
  
  // number of trials and number of events accounting for frequencies
  nobs = std::accumulate(freqn.begin(), freqn.end(), 0);
  nevents = std::inner_product(freqn.begin(), freqn.end(), eventn.begin(), 0);
  if (nevents == 0) {
    for (int i = 0; i < p; ++i) {
      par[i] = (i == 0) ? "(Intercept)" : covariates[i-1];
      b[i] = NaN;
      seb[i] = 0;
      rseb[i] = 0;
      z[i] = NaN;
      expbeta[i] = NaN;
      lb[i] = NaN;
      ub[i] = NaN;
      prob[i] = NaN;
      clparm[i] = "Wald";
    }
    
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < p; ++i) {
        vb(i,j) = 0;    
        rvb(i,j) = 0;
      }
    }
    
    for (int person = 0; person < n; ++person) {
      linear_predictors[person] = offsetn[person];
    }
    switch (link_code) {
    case 1:
      for (int person = 0; person < n; ++person)
        fitted_values[person] = boost_plogis(linear_predictors[person]);
      break;
    case 2:
      for (int person = 0; person < n; ++person)
        fitted_values[person] = boost_pnorm(linear_predictors[person]);
      break;
    case 3:
      for (int person = 0; person < n; ++person)
        fitted_values[person] = boost_pextreme(linear_predictors[person]);
      break;
    }
    
    loglik0 = NaN;
    loglik1 = NaN;
    regloglik0 = NaN;
    regloglik1 = NaN;
    niter = 0;
    fail = true;
  } else {
    // intercept only model
    double num = 0, den = 0;
    for (int i = 0; i < n; ++i) {
      num += freqn[i] * weightn[i] * eventn[i];
      den += freqn[i] * weightn[i];
    }
    if (firth) {
      num += 0.5;
      den += 1.0;
    }
    
    std::vector<double> bint0(p);
    bint0[0] = boost_qlogis(num / den);
    
    std::vector<int> colfit0(1);
    logparams param = {n, link_code, eventn, zn, freqn, weightn, offsetn};
    ListCpp outint = logisregloop(p, bint0, &param, maxiter, eps, firth, colfit0, 1);
    
    std::vector<double> bint = outint.get<std::vector<double>>("coef");
    FlatMatrix vbint = outint.get<FlatMatrix>("var");
    
    ListCpp out;
    
    if (p > 1) {
      // parameter estimates and standard errors for the full model
      std::vector<int> colfit = seqcpp(0,p-1);
      if (!init.empty() && static_cast<int>(init.size()) == p && 
          std::none_of(init.begin(), init.end(), [](double val){ 
            return std::isnan(val); })) {
        out = logisregloop(p, init, &param, maxiter, eps, firth, colfit, p);
      } else {
        out = logisregloop(p, bint, &param, maxiter, eps, firth, colfit, p);
      }
      
      fail = out.get<bool>("fail");
      if (fail) {
        thread_utils::push_thread_warning(
          "logisregloop failed to converge for the full model; "
          "continuing with current results.");
      }
      
      b = out.get<std::vector<double>>("coef");
      vb = out.get<FlatMatrix>("var");
      
      // intercept correction
      if (flic) {
        std::vector<double> lp = offsetn;  // linear predictor excluding intercept
        
        // pointers into data for fastest access
        const double* zptr = zn.data_ptr();
        double* lpptr = lp.data();
        
        for (int col = 1; col < p; ++col) {   // skip intercept column 0
          double coef = b[col];
          if (coef == 0.0) continue;                       // skip zero coefficient
          const double* zcol = zptr + col * n;
          // inner loop is contiguous in memory (good locality)
          for (int row = 0; row < n; ++row) {
            lpptr[row] += coef * zcol[row];
          }
        }
        
        logparams param0 = {n, link_code, eventn, zn, freqn, weightn, lp};
        std::vector<double> bint00(1, bint0[0]);
        ListCpp outint0 = logisregloop(1, bint00, &param0, maxiter, 
                                       eps, 0, colfit0, 1);
        double a = outint0.get<std::vector<double>>("coef")[0];
        double va = outint0.get<FlatMatrix>("var")(0,0);
        
        // update the intercept estimate
        b[0] = a;
        
        // partial derivative of alpha(beta) with respect to beta
        ListCpp derint = f_der_0(p, b, &param, 0);
        FlatMatrix iflic = derint.get<FlatMatrix>("imat");
        std::vector<double> der(p-1);
        for (int i = 0; i < p - 1; ++i) {
          der[i] = -iflic(i+1, 0) / iflic(0, 0);
        }
        
        // update variance-covariance matrix via the delta method
        if (p <= 1) { // only intercept
          vb(0,0) = va;
        } else {
          int p1 = p - 1;                  // length of der and of the submatrix
          const double* derp = der.data();
          double* vbptr = vb.data_ptr();    // column-major
          
          // accumulator for vb(1..p-1,0): indexed 0..p1-1 corresponds to rows 1..p-1
          std::vector<double> col0(p1, 0.0);
          
          // Compute col0 = M * der, where M = vb[1..p-1, 1..p-1]
          // columns of M correspond to vb columns 1..p-1; 
          // for column j (0..p1-1) start at
          // src = vbptr + (j+1)*nrows and element M(row = r+1, col = j+1) is src[r+1]
          for (int j = 0; j < p1; ++j) {
            double dj = derp[j];
            if (dj == 0.0) continue;                   // skip zero coeffs
            const double* src_col = vbptr + (j + 1) * p;
            // accumulate rows 1..p-1 -> indices r=0..p1-1 map to src_col[r+1]
            for (int r = 0; r < p1; ++r) {
              col0[r] += src_col[r + 1] * dj;
            }
          }
          
          // write back vb(1..p-1,0) and vb(0,1..p-1) and compute der^T * col0
          double dot = 0.0;
          // vb( row, col ) -> vbptr[ col*nrows + row ]
          for (int r = 0; r < p1; ++r) {
            double v = col0[r];
            vbptr[0 * p + (r + 1)] = v;          // vb(r+1, 0)
            vbptr[(r + 1) * p + 0] = v;          // vb(0, r+1)
            dot += derp[r] * v;
          }
          
          // final scalar update
          vbptr[0] = va + dot;                       // vb(0,0)  
        }
      }
    } else {
      b = bint;
      vb = vbint;
      out = outint;
      
      if (flic) {
        out = logisregloop(p, bint0, &param, maxiter, eps, 0, colfit0, 1);
        b = out.get<std::vector<double>>("coef");
        vb = outint.get<FlatMatrix>("var");
      }
    }
    
    for (int i = 0; i < p; ++i) {
      seb[i] = std::sqrt(vb(i,i));
    }
    
    for (int i = 0; i < p; ++i) {
      par[i] = (i == 0) ? "(Intercept)" : covariates[i-1];
    }
    
    // linear predictors and fitted values
    std::vector<double> eta = offsetn;
    for (int i = 0; i < p; ++i) {
      double beta = b[i];
      if (beta == 0.0) continue;
      const double* zn_col = zn.data_ptr() + i * n;
      for (int r = 0; r < n; ++r) {
        eta[r] += beta * zn_col[r];
      }
    }
    linear_predictors = eta;
    
    switch (link_code) {
    case 1: // logit
      for (int person = 0; person < n; ++person) {
        fitted_values[person] = boost_plogis(eta[person]);
      }
      break;
    case 2: // probit
      for (int person = 0; person < n; ++person) {
        fitted_values[person] = boost_pnorm(eta[person]);
      }
      break;
    case 3: // cloglog
      for (int person = 0; person < n; ++person) {
        fitted_values[person] = boost_pextreme(eta[person]);
      }
      break;
    }
    
    niter = out.get<int>("iter");
    fail = out.get<bool>("fail");
    
    // robust variance estimates
    if (robust) {
      FlatMatrix ressco = f_ressco_0(p, b, &param);
      
      int nr; // number of rows in the score residual matrix
      std::vector<int> freqr;
      if (!has_id) {
        for (int j = 0; j < p; ++j) {
          double* rcol = ressco.data_ptr() + j * n;
          for (int i = 0; i < n; ++i) {
            rcol[i] *= weightn[i];
          }
        }
        
        nr = n;
        freqr = freqn;
      } else { // need to sum up score residuals by id
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return idn[i] < idn[j];
        });
        
        std::vector<int> id1 = subset(idn, order);
        std::vector<int> idx(1,0);
        for (int i = 1; i < n; ++i) {
          if (id1[i] != id1[i-1]) {
            idx.push_back(i);
          }
        }
        
        int nids = static_cast<int>(idx.size());
        idx.push_back(n);
        
        FlatMatrix ressco1(nids, p);
        for (int j = 0; j < p; ++j) {
          const double* rcol = ressco.data_ptr() + j * n;
          double* rcol1 = ressco1.data_ptr() + j * nids;
          for (int i = 0; i < nids; ++i) {
            double sum = 0.0;
            for (int k = idx[i]; k < idx[i+1]; ++k) {
              int row = order[k];
              sum += weightn[row] * rcol[row];
            }
            rcol1[i] = sum;
          }
        }
        
        std::vector<int> freq1(nids); // cluster frequency
        for (int i = 0; i < nids; ++i) {
          freq1[i] = freqn[order[idx[i]]];
        }
        
        // update the score residuals
        ressco = std::move(ressco1);  
        nr = nids;
        freqr = std::move(freq1);
      }
      
      FlatMatrix D = mat_mat_mult(ressco, vb);
      
      const double* Dptr = D.data_ptr();      // Dcol_j starts at Dptr + j*nr_sz
      double* rvbptr     = rvb.data_ptr();    // rvbcol_k starts at rvbptr + k*p
      
      // rvb(j,k) corresponds to rvbptr[k*p + j] because column-major: column k, row j
      for (int j = 0; j < p; ++j) {
        const double* Dj = Dptr + j * nr;          // pointer to D(:,j)
        for (int k = 0; k <= j; ++k) {
          const double* Dk = Dptr + k * nr;        // pointer to D(:,k)
          double sum = 0.0;
          // inner loop reads Dj[i] and Dk[i] contiguously
          for (int i = 0; i < nr; ++i) {
            sum += freqr[i] * Dj[i] * Dk[i];
          }
          // write into rvb(j,k) and mirror
          rvbptr[k * p + j] = sum;                // rvb(j,k)
          if (j != k) rvbptr[j * p + k] = sum;    // rvb(k,j) = rvb(j,k)  (mirror)
        }
      }
      
      for (int i = 0; i < p; ++i) {
        rseb[i] = std::sqrt(rvb(i, i));
      }
    }
    
    // profile likelihood confidence interval for regression coefficients
    if (plci) {
      double lmax = out.get<double>("loglik");
      double l0 = lmax - 0.5 * xcrit;
      
      if (!(firth && flic)) { // PL CI for all parameters
        for (int k = 0; k < p; ++k) {
          lb[k] = logisregplloop(p, b, &param, maxiter, eps, firth, k, -1, l0);
          ub[k] = logisregplloop(p, b, &param, maxiter, eps, firth, k, 1, l0);
          
          std::vector<int> colfit1(p-1);
          for (int i = 0, j = 0; i < p; ++i) {
            if (i == k) continue;
            colfit1[j++] = i;
          }
          
          std::vector<double> b0(p);
          ListCpp out0 = logisregloop(p, b0, &param, maxiter, eps, 
                                      firth, colfit1, p-1);
          double lmax0 = out0.get<double>("loglik");
          prob[k] = boost_pchisq(-2.0 * (lmax0 - lmax), 1, 0);
          clparm[k] = "PL";
        }
      } else { // Wald CI for intercept and PL CI for slopes
        if (!robust) {
          lb[0] = b[0] - zcrit * seb[0];
          ub[0] = b[0] + zcrit * seb[0];
          prob[0] = boost_pchisq(sq(b[0] / seb[0]), 1, 0);
        } else {
          lb[0] = b[0] - zcrit * rseb[0];
          ub[0] = b[0] + zcrit * rseb[0];
          prob[0] = boost_pchisq(sq(b[0] / rseb[0]), 1, 0);
        }
        clparm[0] = "Wald";
        
        for (int k = 1; k < p; ++k) {
          lb[k] = logisregplloop(p, b, &param, maxiter, eps, firth, k, -1, l0);
          ub[k] = logisregplloop(p, b, &param, maxiter, eps, firth, k, 1, l0);
          
          std::vector<int> colfit1(p-1);
          for (int i = 0; i < k; ++i) {
            colfit1[i] = i;
          }
          for (int i = k + 1; i < p; ++i) {
            colfit1[i-1] = i;
          }
          
          std::vector<double> b0(p);
          ListCpp out0 = logisregloop(p, b0, &param, maxiter, eps, 
                                      firth, colfit1, p-1);
          double lmax0 = out0.get<double>("loglik");
          prob[k] = boost_pchisq(-2.0 * (lmax0 - lmax), 1, 0);
          clparm[k] = "PL";
        }
      }
    } else { // Wald confidence interval for all parameters
      for (int k = 0; k < p; ++k) {
        if (!robust) {
          lb[k] = b[k] - zcrit * seb[k];
          ub[k] = b[k] + zcrit * seb[k];
          prob[k] = boost_pchisq(sq(b[k] / seb[k]), 1, 0);
        } else {
          lb[k] = b[k] - zcrit * rseb[k];
          ub[k] = b[k] + zcrit * rseb[k];
          prob[k] = boost_pchisq(sq(b[k] / rseb[k]), 1, 0);
        }
        clparm[k] = "Wald";
      }
    }
    
    // log-likelihoods
    if (p > 0 && firth) {
      loglik0 = outint.get<double>("loglik");
      loglik1 = out.get<double>("loglik");
      regloglik0 = outint.get<double>("regloglik");
      regloglik1 = out.get<double>("regloglik");
    } else {
      loglik0 = outint.get<double>("loglik");
      loglik1 = out.get<double>("loglik");
      regloglik0 = loglik0;
      regloglik1 = loglik1;
    }
    
    // compute exp(beta)
    for (int i = 0; i < p; ++i) {
      expbeta[i] = std::exp(b[i]);
    }
    
    // compute z statistics
    if (robust) {
      for (int i = 0; i < p; ++i) {
        if (rseb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i] / rseb[i];
        }
      }
    } else {
      for (int i = 0; i < p; ++i) {
        if (seb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i] / seb[i];
        }
      }
    }
  }
  
  // prepare the output data sets
  DataFrameCpp sumstat;
  sumstat.push_back(nobs, "n");
  sumstat.push_back(nevents, "nevents");
  sumstat.push_back(loglik0, "loglik0");
  sumstat.push_back(loglik1, "loglik1");
  sumstat.push_back(niter, "niter");
  sumstat.push_back(p, "p");
  sumstat.push_back(link1, "link");
  sumstat.push_back(robust, "robust");
  sumstat.push_back(firth, "firth");
  sumstat.push_back(flic, "flic");
  sumstat.push_back(fail, "fail");
  
  if (p > 0 && firth) {
    sumstat.push_back(regloglik0, "loglik0_unpenalized");
    sumstat.push_back(regloglik1, "loglik1_unpenalized");
  }
  
  std::vector<double> sebeta = robust ? rseb : seb;
  FlatMatrix vbeta = robust ? rvb : vb;
  
  DataFrameCpp parest;
  parest.push_back(std::move(par), "param");
  parest.push_back(std::move(b), "beta");
  parest.push_back(std::move(sebeta), "sebeta");
  parest.push_back(std::move(z), "z");
  parest.push_back(std::move(expbeta), "expbeta");
  parest.push_back(std::move(lb), "lower");
  parest.push_back(std::move(ub), "upper");
  parest.push_back(std::move(prob), "p");
  parest.push_back(std::move(clparm), "method");
  if (robust) parest.push_back(std::move(seb), "sebeta_naive");
  
  DataFrameCpp fitted;
  fitted.push_back(std::move(linear_predictors), "linear_predictors");
  fitted.push_back(std::move(fitted_values), "fitted_values");
  
  ListCpp result;
  result.push_back(std::move(sumstat), "sumstat");
  result.push_back(std::move(parest), "parest");
  result.push_back(std::move(vbeta), "vbeta");
  if (robust) result.push_back(std::move(vb), "vbeta_naive");
  result.push_back(std::move(fitted), "fitted");
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List logisregRcpp(const Rcpp::DataFrame& data,
                        const std::string& event,
                        const std::vector<std::string>& covariates,
                        const std::string& freq,
                        const std::string& weight,
                        const std::string& offset,
                        const std::string& id,
                        const std::string& link,
                        const std::vector<double>& init,
                        const bool robust,
                        const bool firth,
                        const bool flic,
                        const bool plci,
                        const double alpha,
                        const int maxiter,
                        const double eps) {
  
  DataFrameCpp dfcpp = convertRDataFrameToCpp(data);
  
  ListCpp cpp_result = logisregcpp(
    dfcpp, event, covariates, freq, weight, offset, id, link,
    init, robust, firth, flic, plci, alpha, maxiter, eps
  );
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(cpp_result);
}
