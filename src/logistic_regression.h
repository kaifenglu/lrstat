#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

#ifndef __LOGISTIC_REGRESSION__
#define __LOGISTIC_REGRESSION__

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
