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

double errorSpent(const double t, 
                  const double error,
                  const String sf, 
                  const double sfpar);

List exitprob(const NumericVector& b,
              NumericVector a,
              NumericVector theta,
              const NumericVector& I);

double qtpwexp(const double probability,
               const NumericVector& piecewiseSurvivalTime,
               const NumericVector& lambda,
               const double lowerBound);

List updateGraph(const NumericVector& w, 
                 const NumericMatrix& G, 
                 const IntegerVector& I, 
                 const int j);

#endif // __UTILITIES__
