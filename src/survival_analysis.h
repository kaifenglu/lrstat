#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

#ifndef __SURVIVAL_ANALYSIS__
#define __SURVIVAL_ANALYSIS__

DataFrame kmest(const DataFrame data,
                const std::string rep,
                const std::string stratum,
                const std::string time,
                const std::string event,
                const std::string conftype,
                const double confint);

DataFrame lrtest(const DataFrame data,
                 const std::string rep,
                 const std::string stratum,
                 const std::string treat,
                 const std::string time,
                 const std::string event,
                 const double rho1,
                 const double rho2);

struct aftparams {
  std::string dist;
  IntegerVector strata;
  NumericVector tstart;
  NumericVector tstop;
  IntegerVector status;
  NumericVector weight;
  NumericMatrix z;
};


double f_nllik_1(int p, double *par, void *ex);

void f_nscore_1(int p, double *par, double *gr, void *ex);

NumericMatrix f_info_1(int p, double *par, void *ex);

NumericMatrix f_ressco_1(NumericVector par, void *ex);

List liferegr(const DataFrame data,
              const std::string rep,
              const std::string stratum,
              const std::string time,
              const std::string time2,
              const std::string event,
              const StringVector& covariates,
              const std::string weight,
              const std::string id,
              const std::string dist,
              bool robust);

struct coxparams {
  int nused;
  IntegerVector strata;
  NumericVector tstart;
  NumericVector tstop;
  IntegerVector event;
  NumericVector weight;
  NumericMatrix z;
  IntegerVector order1;
  int method;
};

double f_nllik_2(int p, double *par, void *ex);

void f_nscore_2(int p, double *par, double *gr, void *ex);

NumericMatrix f_info_2(int p, double *par, void *ex);

NumericMatrix f_ressco_2(NumericVector beta, void *ex);

List phregr(const DataFrame data,
            const std::string rep,
            const std::string stratum,
            const std::string time,
            const std::string time2,
            const std::string event,
            const StringVector& covariates,
            const std::string weight,
            const std::string id,
            const std::string ties,
            bool robust);

#endif // __SURVIVAL_ANALYSIS__
