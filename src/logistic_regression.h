#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

#ifndef __LOGISTIC_REGRESSION__
#define __LOGISTIC_REGRESSION__

struct logparams {
  int n;
  std::string link;
  NumericVector y;
  NumericMatrix z;
  NumericVector freq;
  NumericVector weight;
  NumericVector offset;
};

double f_llik_0(int p, NumericVector par, void *ex);
NumericVector f_score_0(int p, NumericVector par, void *ex);
NumericMatrix f_info_0(int p, NumericVector par, void *ex);

double f_pen_llik_0(int p, NumericVector par, void *ex);
NumericVector f_pen_score_0(int p, NumericVector par, void *ex);
NumericVector f_firth_score_0(int p, NumericVector par, void *ex);

NumericMatrix f_ressco_0(int p, NumericVector par, void *ex);

List logisregloop(int p, NumericVector par, void *ex,
                  int maxiter, double eps, bool firth, bool bc,
                  IntegerVector colfit, int ncolfit);

double logisregplloop(int p, NumericVector par, void *ex,
                      int maxiter, double eps, bool firth, bool bc,
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
                 const bool robust,
                 const bool firth,
                 const bool bc,
                 const bool flic,
                 const bool plci,
                 const double alpha);

#endif // __LOGISTIC_REGRESSION__
