#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;

//' @title Set seed
//' @description Sets the R seed in the cpp program based on set.seed() in R.
//'
//' @param seed The seed to use for generating random numbers.
//' @return No return value, called for side effects.
//'
//' @keywords internal
//'
//' @examples
//' set_seed(123)
//'
//' @export
// [[Rcpp::export]]
void set_seed(int seed) {
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}


//' @title Sort a numeric vector
//' @description Sorts a numeric vector in the cpp program.
//'
//' @param x The numeric vector to sort.
//' @return A vector obtained after sorting the input vector.
//'
//' @keywords internal
//'
//' @examples
//' stl_sort(c(3, 4.2, 1))
//'
//' @export
// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

//' @title Find interval numbers of indices
//' @description The implementation of \code{findInterval()} in R from Advanced
//' R by Hadley Wickham. Given a vector of non-decreasing breakpoints in v,
//' find the interval containing each element of x; i.e., if
//' \code{i <- findInterval2(x,v)}, for each index \code{j} in \code{x},
//' v[i[j]] <= x[j] < v[i[j] + 1]
//' where v[0] := -Inf, v[N+1] := +Inf, and \code{N = length(v)}.
//'
//' @param x The numeric vector of interest.
//' @param v The vector of break points.
//' @return A vector of \code{length(x)} with values in \code{0:N} where
//'   \code{N = length(v)}.
//'
//' @keywords internal
//'
//' @examples
//' x <- 2:18
//' v <- c(5, 10, 15) # create two bins [5,10) and [10,15)
//' cbind(x, findInterval2(x, v))
//'
//' @export
// [[Rcpp::export]]
IntegerVector findInterval2(NumericVector x, NumericVector v) {
  IntegerVector out(x.size());
  
  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;
  
  NumericVector::iterator x_begin=x.begin(), x_end=x.end();
  NumericVector::iterator v_begin=v.begin(), v_end=v.end();
  
  for(it = x_begin, out_it = out.begin(); it != x_end; ++it, ++out_it) {
    pos = std::upper_bound(v_begin, v_end, *it);
    *out_it = std::distance(v_begin, pos);
  }
  
  return out;
}


#include <algorithm>
#define ITMAX 100
#define EPS 3.0e-8
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//' @title Brent's method for root-finding
//' @description Using Brent's method, find the root of a function known to
//' lie between x1 and x2. Program based on the book - Numerical Recipes in C
//' The Art of Scientific Computing - Second Edition, by William H. Press,
//' Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. It mimics
//' the uniroot() function in R.
//'
//' @param f Name of the univariate objective function.
//' @param x1 One end of the interval bracket.
//' @param x2 The other end of the interval bracket.
//' @param tol The tolerance limit for stopping the iteration.
//'
//' @return The root x between x1 and x2 such that f(x) = 0.
//'
//' @examples
//' brent(sin, -1, 1, 0.0001)
//' @export
// [[Rcpp::plugins(cpp11)]]
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol) {
  int iter;
  double a=x1, b=x2, c=x2, d, d1, min1, min2;
  double fa=f(a), fb=f(b), fc, p, q, r, s, tol1, xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    stop("Root must be bracketed in brent");
  }
  
  fc = fb;
  for (iter=1; iter<=ITMAX; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;     // Rename a, b, c and adjust bounding interval d
      fc = fa;
      d = b - a;
      d1 = d;
    }
    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    // Convergence check
    tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
    xm = 0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) {
      return b;
    }
    
    if (fabs(d1) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb/fa; // Attempt inverse quadratic interpolation
      if (a == c) {
        p = 2.0*xm*s;
        q = 1.0-s;
      } else {
        q = fa/fc;
        r = fb/fc;
        p = s*(2.0*xm*q*(q-r) - (b-a)*(r-1.0));
        q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) {
        q = -q;  // Check whether in bounds
      }
      p = fabs(p);
      min1 = 3.0*xm*q - fabs(tol1*q);
      min2 = fabs(d1)*fabs(q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        d1 = d;  // Accept interpolation
        d = p/q;
      } else {  // Interpolation failed, use bisection
        d = xm;
        d1 = d;
      }
    } else {  // Bounds decreasing too slowly, use bisection
      d = xm;
      d1 = d;
    }
    a = b;  // Move last best guess to a
    fa = fb;
    if (fabs(d) > tol1) { // Evaluate new trial root
      b += d;
    } else {
      b += SIGN(tol1, xm);
    }
    fb = f(b);
  }
  stop("Maximum number of iterations exceeded in brent");
  return 0.0; // Never get here
}



//' @title Stagewise exit probabilities
//' @description Obtains the stagewise exit probabilities for both efficacy and
//' futility stopping.
//'
//' @param b Upper boundaries on the z-test statistic scale.
//' @param a Lower boundaries on the z-test statistic scale. Defaults to
//' \code{c(rep(-6.0, kMax-1), b[kMax])} if left unspecified, where
//' \code{kMax = length(b)}.
//' @param theta Stagewise parameter of interest, e.g., \code{-U/V} for
//' weighted log-rank test, where \code{U} is the mean and \code{V} is
//' the variance of the weighted log-rank test score statistic at each stage.
//' For proportional hazards and conventional log-rank test, use the
//' scalar input, \code{theta = -log(HR)}.
//' @param I Stagewise cumulative information, e.g., \code{V}, the variance
//' of the weighted log-rank test score statistic at each stage. For
//' conventional log-rank test, information can be approximated by
//' \code{phi*(1-phi)*D}, where \code{phi} is the probability of being
//' allocated to the active arm, and \code{D} is the total number of events at
//' each stage.
//'
//' @return A list of stagewise exit probabilities: one vector for efficacy
//' stopping probabilities, and the other vector for futility stopping
//' probabilities.
//'
//' @examples
//' exitprob(b = c(3.471, 2.454, 2.004), a = NA,
//'          theta = -log(0.6), I = c(50, 100, 150)/4)
//'
//' @export
// [[Rcpp::export]]
List exitprob(const NumericVector& b,
              NumericVector a,
              NumericVector theta,
              const NumericVector& I) {
  
  // Integer value controlling grid for numerical integration as in
  // Jennison and Turnbull (2000)
  const int r = 18;
  
  // variable declarations
  // kMax is the total number of stages
  // m0, z0, h0 for the previous stage
  // m, z, h for the current stage
  int kMax=b.size(), r1=6*r-1, r2=12*r-3, i0, i1=0, i2=r1-1, i, j,
    m0=r2, m1=r1, m=r2;
  double t, tlower, tupper, xlower, xupper;
  
  NumericVector sqrtI(kMax), thetaSqrtI(kMax), thetaI(kMax), dI(kMax),
  dThetaI(kMax), exitProbUpper(kMax), exitProbLower(kMax),
  shift(r1), x1(r1), x(r1), z0(r2), z(r2), w(r2), h0(r2), h(r2);
  
  // set default parameter values
  if (is_true(any(is_na(a)))) {
    NumericVector tem(kMax);
    for (i=0; i<kMax; i++) {
      if (i<kMax-1) {
        tem[i] = -6.0;
      } else {
        tem[i] = b[i];
      }
    }
    a = tem;
  }
  
  
  
  if (theta.size()==1) {
    theta = rep(theta, kMax);
  }
  
  if (a.size()!=kMax || theta.size()!=kMax || I.size()!=kMax) {
    stop("The input parameters must have the same length");
  }
  
  // edit check
  for (i=0; i<kMax; i++) {
    if (a[i] > b[i]) {
      stop("Lower bounds (a) must be less than upper bounds (b)");
    }
  }
  
  if (I[0] <= 0) {
    stop("Elements of I must be positive");
  } else if (kMax > 1 && is_true(any(diff(I) <= 0))) {
    stop("Elements of I must be increasing");
  }
  
  // constant shifts relative to the means, use floating point computation
  for (i=0; i<r1; i++) {
    if (i < r-1) {
      shift[i] = -3 - 4*log(r/(i+1.0));
    } else if (i < 5*r) {
      shift[i] = -3 + 3*(i+1.0-r)/(2*r);
    } else {
      shift[i] = 3 + 4*log(r/(6*r-i-1.0));
    }
  }
  
  
  // obtain various vectors associated with theta and I
  for (j=0; j<kMax; j++) {
    sqrtI[j] = sqrt(I[j]);
    thetaSqrtI[j] = theta[j]*sqrtI[j];
    thetaI[j] = theta[j]*I[j];
    if (j==0) {
      dI[j] = I[j];
      dThetaI[j] = thetaI[j];
    } else {
      dI[j] = I[j] - I[j-1];
      dThetaI[j] = thetaI[j] - thetaI[j-1];
    }
  }
  
  // loop over stages
  for (j=0; j<kMax; j++) {
    
    // initialize x values
    for (i=0; i<r1; i++) {
      x1[i] = thetaSqrtI[j] + shift[i];
    }
    
    // trim off x values outside (a[j], b[j])
    // trim from below
    if (a[j] >= x1[0]) {
      i1 = 0;
      while (x1[i1] <= a[j]) {
        i1++;
      }
      i1--;
      xlower = a[j]; // lower bound on x
    } else {
      i1 = 0;
      xlower = x1[0];
    }
    
    // trim from above
    if (b[j] <= x1[r1-1]) {
      i2 = r1-1;
      while (x1[i2] >= b[j]) {
        i2--;
      }
      i2++;
      xupper = b[j]; // upper bound on x
    } else {
      i2 = r1-1;
      xupper = x1[r1-1];
    }
    
    // save the trimmed portion to x
    m1 = i2 - i1 + 1;
    x[0] = xlower;
    x[m1-1] = xupper;
    for (i=1; i<m1-1; i++) {
      x[i] = x1[i+i1];
    }
    
    // derive the grid points for z
    m = 2*m1 - 1;
    
    // odd grid points;
    for (i=0; i<m1; i++) {
      z[2*i] = x[i];
    }
    
    // even grid points;
    for (i=0; i<m1-1; i++) {
      z[2*i+1] = (z[2*i] + z[2*i+2])/2;
    }
    
    
    // derive the weights
    w[0] = 1.0/6*(z[2] - z[0]);
    
    for (i0=1; i0<=m1-2; i0++) {
      i = 2*i0;
      w[i] = 1.0/6*(z[i+2] - z[i-2]);
    }
    
    for (i0=1; i0<=m1-1; i0++) {
      i = 2*i0-1;
      w[i] = 4.0/6*(z[i+1] - z[i-1]);
    }
    
    w[m-1] = 1.0/6*(z[m-1] - z[m-3]);
    
    
    // first stage is easy
    if (j==0) {
      // exit probabilities
      exitProbUpper[j] = R::pnorm(-b[j] + thetaSqrtI[j], 0.0, 1.0, 1, 0);
      exitProbLower[j] = R::pnorm(a[j] - thetaSqrtI[j], 0.0, 1.0, 1, 0);
      
      // prepare h0, m0, z0 for the next stage
      if (kMax > 1) {
        for (i=0; i<m; i++) {
          h0[i] = w[i]*R::dnorm(z[i] - thetaSqrtI[j], 0.0, 1.0, 0);
        }
        
        m0 = m;
        z0 = z+0.0; // adding 0.0 to avoid passing by reference
      }
      
    } else {
      // calculate exit probabilities using h0 from the previous stage
      for (i0=0; i0<m0; i0++) {
        tupper = (z0[i0]*sqrtI[j-1] - b[j]*sqrtI[j] + dThetaI[j])/sqrt(dI[j]);
        tlower = (-z0[i0]*sqrtI[j-1] + a[j]*sqrtI[j] - dThetaI[j])/sqrt(dI[j]);
        exitProbUpper[j] += h0[i0]*R::pnorm(tupper, 0.0, 1.0, 1, 0);
        exitProbLower[j] += h0[i0]*R::pnorm(tlower, 0.0, 1.0, 1, 0);
      }
      
      // prepare h0, m0, z0 for the next stage
      if (j < kMax-1) {
        for (i=0; i<m; i++) {
          h[i] = 0;
          for (i0=0; i0<m0; i0++) {
            t = (z[i]*sqrtI[j] - z0[i0]*sqrtI[j-1] - dThetaI[j])/sqrt(dI[j]);
            h[i] += h0[i0]*R::dnorm(t, 0.0, 1.0, 0);
          }
          h[i] *= w[i]*sqrt(I[j]/dI[j]); // factors invariant to i0
        }
        
        h0 = h+0.0; // adding 0.0 to avoid passing by reference
        m0 = m;
        z0 = z+0.0;
      }
    }
    
  }
  
  // return a list of stagewise exit probabilities
  return List::create(Named("exitProbUpper") = exitProbUpper,
                      Named("exitProbLower") = exitProbLower);
  
}



//' @title Quantile function of truncated piecewise exponential distribution
//' @description Obtains the quantile of a piecewise expoenential distribution
//' given that it exceeds a specified lower bound.
//'
//' @param probability The scalar probability corresponding to the quantile.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @param lowerBound The left truncation time point for the survival time.
//' Defaults to 0 for no truncation.
//'
//' @return The quantile x such that
//' P(X > x | X > lowerBound) = 1 - probability.
//'
//' @keywords internal
//'
//' @examples
//' qtpwexp(probability = 0.3, piecewiseSurvivalTime = c(0, 6, 9, 15),
//'         lambda = c(0.025, 0.04, 0.015, 0.007), lowerBound = 0)
//'
//' @export
// [[Rcpp::export]]
double qtpwexp(const double probability,
               const NumericVector& piecewiseSurvivalTime,
               const NumericVector& lambda,
               const double lowerBound) {
  
  int j, j1, m;
  double q, v, v1;
  
  // cumulative hazard from lowerBound until the quantile
  v1 = -log(1 - probability);
  
  // identify the time interval containing the lowerBound
  m = piecewiseSurvivalTime.size();
  for (j=0; j<m; j++) {
    if (piecewiseSurvivalTime[j] > lowerBound) break;
  }
  j1 = (j==0 ? 0 : j-1); // to handle floating point precision
  
  if (j1 == m-1) { // in the last interval
    q = (lambda[j1]==0.0 ? 1.0e+8 : v1/lambda[j1] + lowerBound);
  } else {
    // accumulate the pieces on the cumulative hazard scale
    v = 0;
    for (j=j1; j<m-1; j++) {
      if (j==j1) {
        v += lambda[j]*(piecewiseSurvivalTime[j+1] - lowerBound);
      } else {
        v += lambda[j]*(piecewiseSurvivalTime[j+1] - piecewiseSurvivalTime[j]);
      }
      if (v >= v1) break;
    }
    
    if (j == m-1) { // in the last interval
      q = (lambda[j]==0.0 ? 1.0e+8 :
             (v1 - v)/lambda[j] + piecewiseSurvivalTime[j]);
    } else {
      q = (lambda[j]==0.0 ? 1.0e+8 :
             piecewiseSurvivalTime[j+1] - (v - v1)/lambda[j]);
    }
  }
  
  return q;
}

