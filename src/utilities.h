#include <Rcpp.h>
using namespace Rcpp;

#ifndef __UTILITIES__
#define __UTILITIES__

void set_seed(int seed);

NumericVector stl_sort(const NumericVector& x);

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

double qtpwexpcpp(const double probability,
                  const NumericVector& piecewiseSurvivalTime,
                  const NumericVector& lambda,
                  const double lowerBound);


#endif // __UTILITIES__
