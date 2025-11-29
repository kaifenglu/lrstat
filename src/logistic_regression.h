#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

#ifndef __LOGISTIC_REGRESSION__
#define __LOGISTIC_REGRESSION__

struct logparams {
  int n;
  int link_code; // 0: logit, 1: probit, 2: cloglog
  NumericVector y;
  NumericMatrix z;
  NumericVector freq;
  NumericVector weight;
  NumericVector offset;
};

List f_der_0(int p, const NumericVector& par, void *ex, bool firth);

NumericMatrix f_ressco_0(int p, const NumericVector& par, void *ex);

List logisregloop(int p, const NumericVector& par, void *ex,
                  int maxiter, double eps, bool firth,
                  const IntegerVector& colfit, int ncolfit);

double logisregplloop(int p, const NumericVector& par, void *ex,
                      int maxiter, double eps, bool firth,
                      int k, int which, double l0);

List logisregcpp(const DataFrame data,
                 const StringVector& rep,
                 const std::string event,
                 const StringVector& covariates,
                 const std::string freq,
                 const std::string weight,
                 const std::string offset,
                 const std::string id,
                 const std::string link,
                 const NumericVector& init,
                 const bool robust,
                 const bool firth,
                 const bool flic,
                 const bool plci,
                 const double alpha,
                 const int maxiter,
                 const double eps);                      

#endif // __LOGISTIC_REGRESSION__
