#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

#ifndef __SURVIVAL_ANALYSIS__
#define __SURVIVAL_ANALYSIS__

NumericVector fsurvci(double surv, double sesurv, String ct, double z);

DataFrame kmest(const DataFrame data,
                const StringVector& rep,
                const StringVector& stratum,
                const std::string time,
                const std::string event,
                const std::string conftype,
                const double conflev);

DataFrame kmdiff(const DataFrame data,
                 const StringVector& rep,
                 const StringVector& stratum,
                 const std::string treat,
                 const std::string time,
                 const std::string event,
                 const double milestone,
                 const double survDiffH0,
                 const double conflev);

DataFrame lrtest(const DataFrame data,
                 const StringVector& rep,
                 const StringVector& stratum,
                 const std::string treat,
                 const std::string time,
                 const std::string event,
                 const double rho1,
                 const double rho2);


DataFrame rmest(const DataFrame data,
                const StringVector& rep,
                const StringVector& stratum,
                const std::string time,
                const std::string event,
                const double milestone,
                const double conflev,
                const bool biascorrection);

DataFrame rmdiff(const DataFrame data,
                 const StringVector& rep,
                 const StringVector& stratum,
                 const std::string treat,
                 const std::string time,
                 const std::string event,
                 const double milestone,
                 const double rmstDiffH0,
                 const double conflev,
                 const bool biascorrection);

struct aftparams {
  std::string dist;
  IntegerVector strata;
  NumericVector tstart;
  NumericVector tstop;
  IntegerVector status;
  NumericVector weight;
  NumericVector offset;
  NumericMatrix z;
  int nstrata;
};


double f_llik_1(int p, NumericVector par, void *ex);

NumericVector f_score_1(int p, NumericVector par, void *ex);

NumericMatrix f_info_1(int p, NumericVector par, void *ex);

NumericMatrix f_ressco_1(int p, NumericVector par, void *ex);

NumericMatrix f_jj_1(int p, NumericVector par, void *ex);

List liferegloop(int p, NumericVector par, void *ex,
                 int maxiter, double eps,
                 IntegerVector colfit, int ncolfit);

double liferegplloop(int p, NumericVector par, void *ex,
                     int maxiter, double eps,
                     int k, int which, double l0);

List liferegcpp(const DataFrame data,
                const StringVector& rep,
                const StringVector& stratum,
                const std::string time,
                const std::string time2,
                const std::string event,
                const StringVector& covariates,
                const std::string weight,
                const std::string offset,
                const std::string id,
                const std::string dist,
                const bool robust,
                const bool plci,
                const double alpha);

struct coxparams {
  int nused;
  IntegerVector strata;
  NumericVector tstart;
  NumericVector tstop;
  IntegerVector event;
  NumericVector weight;
  NumericVector offset;
  NumericMatrix z;
  IntegerVector order1;
  int method;
};

double f_llik_2(int p, NumericVector par, void *ex);

NumericVector f_score_2(int p, NumericVector par, void *ex);

NumericMatrix f_info_2(int p, NumericVector par, void *ex);

double f_pen_llik_2(int p, NumericVector par, void *ex);

NumericVector f_pen_score_2(int p, NumericVector par, void *ex);

NumericMatrix f_ressco_2(int p, NumericVector par, void *ex);

NumericMatrix f_jj_2(int p, NumericVector par, void *ex);

List phregloop(int p, NumericVector par, void *ex,
               int maxiter, double eps, bool firth,
               IntegerVector colfit, int ncolfit);

double phregplloop(int p, NumericVector par, void *ex,
                   int maxiter, double eps, bool firth,
                   int k, int which, double l0);

List f_basehaz(int p, NumericVector par, void *ex);

NumericVector f_resmart(int p, NumericVector par, void *ex);

List f_ressch(int p, NumericVector par, void *ex);

List phregcpp(const DataFrame data,
              const StringVector& rep,
              const StringVector& stratum,
              const std::string time,
              const std::string time2,
              const std::string event,
              const StringVector& covariates,
              const std::string weight,
              const std::string offset,
              const std::string id,
              const std::string ties,
              const bool robust,
              const bool est_basehaz,
              const bool est_resid,
              const bool firth,
              const bool plci,
              const double alpha);

DataFrame survfit_phregcpp(const int p,
                           const NumericVector& beta,
                           const NumericMatrix& vbeta,
                           DataFrame basehaz,
                           DataFrame newdata,
                           const StringVector& covariates,
                           const StringVector& stratum,
                           const std::string offset,
                           const std::string id,
                           const std::string tstart,
                           const std::string tstop,
                           const bool sefit,
                           const String conftype,
                           const double conflev);

List residuals_phregcpp(const int p,
                        const NumericVector& beta,
                        DataFrame data,
                        const StringVector& stratum,
                        const std::string time,
                        const std::string time2,
                        const std::string event,
                        const StringVector& covariates,
                        const std::string weight,
                        const std::string offset,
                        const std::string id,
                        const std::string ties,
                        const std::string type);

#endif // __SURVIVAL_ANALYSIS__
