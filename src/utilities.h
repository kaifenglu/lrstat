#include <Rcpp.h>
using namespace Rcpp;

#ifndef __UTILITIES__
#define __UTILITIES__

void set_seed(int seed);

NumericVector stl_sort(const NumericVector& x);

IntegerVector which(const LogicalVector& vector);

IntegerVector findInterval2(NumericVector x,
                            NumericVector breaks);

double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol);

double errorSpentcpp(const double t,
                     const double error,
                     const String sf,
                     const double sfpar);

List exitprobcpp(const NumericVector& b,
                 const NumericVector& a,
                 const NumericVector& theta,
                 const NumericVector& I);

NumericVector ptpwexpcpp(const NumericVector& q,
                         const NumericVector& piecewiseSurvivalTime,
                         const NumericVector& lambda,
                         const double lowerBound);

double qtpwexpcpp1(const double p,
                   const NumericVector& piecewiseSurvivalTime,
                   const NumericVector& lambda,
                   const double lowerBound);

NumericVector qtpwexpcpp(const NumericVector& p,
                         const NumericVector& piecewiseSurvivalTime,
                         const NumericVector& lambda,
                         const double lowerBound);

NumericVector rtpwexpcpp(const int n,
                         const NumericVector& piecewiseSurvivalTime,
                         const NumericVector& lambda,
                         const double lowerBound);

NumericVector getBoundcpp(const int k,
                          const NumericVector& informationRates,
                          const double alpha,
                          const String typeAlphaSpending,
                          const double parameterAlphaSpending,
                          const NumericVector& userAlphaSpending,
                          const NumericVector& spendingTime,
                          const LogicalVector& efficacyStopping);

List getPower(const double alpha,
              const int kMax,
              const NumericVector& b,
              const NumericVector& theta,
              const NumericVector& I,
              const std::string bsf,
              const double bsfpar,
              const NumericVector& st,
              const LogicalVector& futilityStopping,
              const NumericVector& w);

double intnorm(const std::function<double(double)>& f,
               double mu, double sigma, double a, double b);

NumericVector polint(const NumericVector& xa,
                     const NumericVector& ya,
                     int n, double x);

double trapzd(const std::function<double(double)>& f,
              double a, double b, int n);

double qromb(const std::function<double(double)>& f,
             double a, double b);

NumericVector mini(const std::function<double(double)>& f,
                   double x1, double x2, double tol);

#endif // __UTILITIES__
