#include <Rcpp.h>
using namespace Rcpp;

#ifndef __UTILITIES__
#define __UTILITIES__

void set_seed(int seed);

NumericVector stl_sort(NumericVector x);

IntegerVector findInterval2(NumericVector x,
                            NumericVector breaks);

double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol);

List exitprob(const NumericVector& b,
              NumericVector a,
              NumericVector theta,
              const NumericVector& I);


double qtpwexp(const double probability,
               const NumericVector& piecewiseSurvivalTime,
               const NumericVector& lambda,
               const double lowerBound);

#endif // __UTILITIES__
