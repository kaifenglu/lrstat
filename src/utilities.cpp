#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "utilities.h"

using namespace Rcpp;


// [[Rcpp::export]]
void set_seed(int seed) {
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}


// [[Rcpp::export]]
NumericVector stl_sort(const NumericVector& x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}


// Function to find the indices of all TRUE elements in a logical vector
IntegerVector which(const LogicalVector& vector) {
  IntegerVector true_indices;
  for (int i = 0; i < vector.size(); i++) {
    if (vector[i]) {
      true_indices.push_back(i);
    }
  }
  return true_indices;
}


//' @title Find interval numbers of indices
//' @description The implementation of \code{findInterval()} in R from
//' Advanced R by Hadley Wickham. Given a vector of non-decreasing
//' breakpoints in v, find the interval containing each element of x; i.e.,
//' if \code{i <- findInterval3(x,v)}, for each index \code{j} in \code{x},
//' \code{v[i[j]] <= x[j] < v[i[j] + 1]}, where \code{v[0] := -Inf},
//' \code{v[N+1] := +Inf}, and \code{N = length(v)}.
//'
//' @param x The numeric vector of interest.
//' @param v The vector of break points.
//' @return A vector of \code{length(x)} with values in \code{0:N} where
//'   \code{N = length(v)}.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' x <- 2:18
//' v <- c(5, 10, 15) # create two bins [5,10) and [10,15)
//' cbind(x, findInterval3(x, v))
//'
//' @export
// [[Rcpp::export]]
IntegerVector findInterval3(NumericVector x, NumericVector v) {
  IntegerVector out(x.size());

  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;

  NumericVector::iterator x_begin=x.begin(), x_end=x.end();
  NumericVector::iterator v_begin=v.begin(), v_end=v.end();

  for(it = x_begin, out_it = out.begin(); it != x_end; ++it, ++out_it) {
    pos = std::upper_bound(v_begin, v_end, *it);
    *out_it = static_cast<int>(std::distance(v_begin, pos));
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
//' Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery.
//' It mimics the uniroot() function in R.
//'
//' @param f Name of the univariate objective function.
//' @param x1 One end of the interval bracket.
//' @param x2 The other end of the interval bracket.
//' @param tol The tolerance limit for stopping the iteration.
//'
//' @return The root x between x1 and x2 such that f(x) = 0.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' brent(sin, -1, 1, 0.0001)
//' @export
//'
// [[Rcpp::plugins(cpp11)]]
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol) {
  int iter;
  double a=x1, b=x2, c=x2, d, d1 = 0.0, min1, min2;
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



// [[Rcpp::export]]
double errorSpentcpp(const double t = NA_REAL,
                     const double error = NA_REAL,
                     const String sf = NA_STRING,
                     const double sfpar = NA_REAL) {
  if (error <= 0 || error >= 1) {
    stop("error must be a number between 0 and 1");
  }
  if (t <= 0 || t > 1) {
    stop("t must be a number between 0 and 1");
  }

  std::string asf = sf;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double aval;
  if (asf == "sfp") {
    aval = error*log(1 + (exp(1) - 1)*t);
  } else if (asf == "sfof") {
    aval = R::qnorm(1-error/2, 0, 1, 1, 0);
    aval = 2*(1 - R::pnorm(aval/sqrt(t), 0, 1, 1, 0));
  } else if (asf == "sfkd") {
    if (R_isnancpp(sfpar)) {
      stop("Parameter sfpar is missing for sfKD");
    } else if (sfpar <= 0) {
      stop ("sfpar must be positive for sfKD");
    } else {
      aval = error*pow(t, sfpar);
    }
  } else if (asf == "sfhsd") {
    if (R_isnancpp(sfpar)) {
      stop("Parameter sfpar is missing for sfHSD");
    } else if (sfpar == 0) {
      aval = error*t;
    } else {
      aval = error*(1 - exp(-sfpar*t))/(1 - exp(-sfpar));
    }
  } else {
    stop("Invalid spending function");
  }
  return aval;
}



// [[Rcpp::export]]
List exitprobcpp(const NumericVector& b,
                 const NumericVector& a,
                 const NumericVector& theta,
                 const NumericVector& I) {

  NumericVector a1 = clone(a);
  NumericVector theta1 = clone(theta);
  NumericVector I1 = clone(I);

  // Integer value controlling grid for numerical integration as in
  // Jennison and Turnbull (2000)
  const int r = 18;

  // variable declarations
  // kMax is the total number of stages
  // m0, z0, h0 for the previous stage
  // m, z, h for the current stage
  int kMax=static_cast<int>(b.size());
  int r1=6*r-1, r2=12*r-3, i0, i1=0, i2=r1-1, i, j, m0=r2, m1=r1, m=r2;
  double t, tlower, tupper, xlower, xupper;

  NumericVector sqrtI(kMax), thetaSqrtI(kMax), thetaI(kMax), dI(kMax),
  dThetaI(kMax), exitProbUpper(kMax), exitProbLower(kMax),
  shift(r1), x1(r1), x(r1), z0(r2), z(r2), w(r2), h0(r2), h(r2);

  // set default parameter values
  if (is_false(any(is_na(a)))) {
    if (a.size() != kMax) {
      stop("Invalid length for a");
    }
  } else {
    NumericVector tem(kMax);
    for (i=0; i<kMax; i++) {
      if (i<kMax-1) {
        tem[i] = -6.0;
      } else {
        tem[i] = b[i];
      }
    }
    a1 = tem;
  }

  // edit check of boundaries
  for (i=0; i<kMax; i++) {
    if (a1[i] > b[i]) {
      stop("Lower bounds (a) must be less than upper bounds (b)");
    }
  }


  if (is_false(any(is_na(theta)))) {
    if (theta.size() == 1) {
      theta1 = rep(theta, kMax);
    } else if (theta.size() != kMax) {
      stop("Invalid length for theta");
    }
  } else {
    theta1 = rep(0, kMax);
  }


  if (is_false(any(is_na(I)))) {
    if (I.size() != kMax) {
      stop("Invalid length for I");
    } else if (I[0] <= 0) {
      stop("Elements of I must be positive");
    } else if (kMax > 1 && is_true(any(diff(I) <= 0))) {
      stop("Elements of I must be increasing");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    I1 = NumericVector(tem);
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
    sqrtI[j] = sqrt(I1[j]);
    thetaSqrtI[j] = theta1[j]*sqrtI[j];
    thetaI[j] = theta1[j]*I1[j];
    if (j==0) {
      dI[j] = I1[j];
      dThetaI[j] = thetaI[j];
    } else {
      dI[j] = I1[j] - I1[j-1];
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
    if (a1[j] >= x1[0]) {
      i1 = 0;
      while (x1[i1] <= a1[j]) {
        i1++;
      }
      i1--;
      xlower = a1[j]; // lower bound on x
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
      exitProbLower[j] = R::pnorm(a1[j] - thetaSqrtI[j], 0.0, 1.0, 1, 0);

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
        tupper = (z0[i0]*sqrtI[j-1] - b[j]*sqrtI[j] +
          dThetaI[j])/sqrt(dI[j]);
        tlower = (-z0[i0]*sqrtI[j-1] + a1[j]*sqrtI[j] -
          dThetaI[j])/sqrt(dI[j]);
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
          h[i] *= w[i]*sqrt(I1[j]/dI[j]); // factors invariant to i0
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


// [[Rcpp::export]]
NumericVector ptpwexpcpp(const NumericVector& q,
                         const NumericVector& piecewiseSurvivalTime,
                         const NumericVector& lambda,
                         const double lowerBound,
                         const bool lowertail,
                         const bool logp) {
  int n = static_cast<int>(q.size());
  NumericVector p(n);
  for (int h=0; h<n; h++) {
    if (q[h] <= lowerBound) {
      p[h] = 0;
    } else {
      NumericVector y = NumericVector::create(lowerBound, q[h]);
      IntegerVector i = findInterval3(y, piecewiseSurvivalTime);
      double v;
      if (i[0] == i[1]) {
        v = lambda[i[0]-1]*(q[h] - lowerBound);
      } else {
        v = lambda[i[0]-1]*(piecewiseSurvivalTime[i[0]] - lowerBound);
        for (int j=i[0]; j<i[1]-1; j++) {
          v += lambda[j]*(piecewiseSurvivalTime[j+1] -
            piecewiseSurvivalTime[j]);
        }
        v += lambda[i[1]-1]*(q[h] - piecewiseSurvivalTime[i[1]-1]);
      }
      p[h] = 1 - exp(-v);
    }
  }

  if (!lowertail) p = 1.0 - p;
  if (logp) p = log(p);

  return p;
}


// [[Rcpp::export]]
double qtpwexpcpp1(const double p,
                   const NumericVector& piecewiseSurvivalTime,
                   const NumericVector& lambda,
                   const double lowerBound,
                   const bool lowertail,
                   const bool logp) {
  int j, j1, m = static_cast<int>(piecewiseSurvivalTime.size());
  double q, u = p, v, v1;

  // cumulative hazard from lowerBound until the quantile
  if (logp) u = exp(p);
  if (!lowertail) u = 1.0 - u;

  v1 = -log(1.0 - u);

  // identify the time interval containing the lowerBound
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
        v += lambda[j]*(piecewiseSurvivalTime[j+1] -
          piecewiseSurvivalTime[j]);
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


// [[Rcpp::export]]
NumericVector qtpwexpcpp(const NumericVector& p,
                         const NumericVector& piecewiseSurvivalTime,
                         const NumericVector& lambda,
                         const double lowerBound,
                         const bool lowertail,
                         const bool logp) {
  int n = static_cast<int>(p.size());
  NumericVector q(n);
  for (int h=0; h<n; h++) {
    q[h] = qtpwexpcpp1(p[h], piecewiseSurvivalTime, lambda, lowerBound,
                       lowertail, logp);
  }

  return q;
}


// [[Rcpp::export]]
NumericVector rtpwexpcpp(
    const int n = NA_INTEGER,
    const NumericVector& piecewiseSurvivalTime = NA_REAL,
    const NumericVector& lambda = NA_REAL,
    const double lowerBound = NA_REAL) {

  NumericVector p(n);
  for (int i=0; i<n; i++) {
    p[i] = R::runif(0,1);
  }

  return qtpwexpcpp(p, piecewiseSurvivalTime, lambda, lowerBound, 1, 0);
}


// [[Rcpp::export]]
NumericVector getBoundcpp(
    const int k = NA_INTEGER,
    const NumericVector& informationRates = NA_REAL,
    const double alpha = NA_REAL,
    const String typeAlphaSpending = NA_STRING,
    const double parameterAlphaSpending = NA_REAL,
    const NumericVector& userAlphaSpending = NA_REAL,
    const NumericVector& spendingTime = NA_REAL,
    const LogicalVector& efficacyStopping = NA_LOGICAL) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector spendingTime1 = clone(spendingTime);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);

  if (k == NA_INTEGER) {
    stop("k must be provided");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != k) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (k > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[k-1] > 1) {
      stop("informationRates must not exceed 1");
    }
  } else {
    IntegerVector tem = seq_len(k);
    informationRates1 = NumericVector(tem)/(k+0.0);
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != k) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (k > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[k-1] > 1) {
      stop("spendingTime must not exceed 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != k) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[k-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, k);
  }

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  if (asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < k) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (k > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[k-1] > alpha) {
      stop("userAlphaSpending must not exceed the specified alpha");
    }
  }

  if (asf=="of" || asf=="p" || asf=="wt") {
    IntegerVector tem = seq_len(k);
    NumericVector informationRates2 = NumericVector(tem)/(k+0.0);
    if (max(abs(informationRates1 - informationRates2)) > 1.0e-6) {
      warning("Equal spacing is used for OF, P, and WT boundaries");
    }
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  NumericVector theta(k); // mean values under H0, initialized to zero
  IntegerVector tem = seq_len(k);
  NumericVector I = NumericVector(tem);
  NumericVector t = clone(informationRates1); // info time for test stat
  NumericVector s = clone(spendingTime1); // spending time for alpha-spending
  NumericVector criticalValues(k);

  if (asf == "none") {
    for (int i=0; i<k-1; i++) {
      criticalValues[i] = 6.0;
    }
    criticalValues[k-1] = R::qnorm(1-alpha, 0, 1, 1, 0);
  } else if (asf == "of" || asf == "p" || asf == "wt") {
    double Delta;
    if (asf == "of") {
      Delta = 0;
    } else if (asf == "p") {
      Delta = 0.5;
    } else {
      Delta = asfpar;
    }

    auto f = [k, alpha, Delta, theta, I,
              efficacyStopping1] (double aval)->double {
      NumericVector u(k), l(k);
      for (int i=0; i<k; i++) {
        u[i] = aval*pow((i+1.0)/k, Delta-0.5);
        if (!efficacyStopping1[i]) u[i] = 6.0;
        l[i] = -6.0;
      }

      List probs = exitprobcpp(u, l, theta, I);
      double cpu = sum(NumericVector(probs[0]));
      return cpu - alpha;
    };

    double cwt = brent(f, 0.0, 10.0, 1.0e-6);
    for (int i=0; i<k; i++) {
      criticalValues[i] = cwt*pow((i+1.0)/k, Delta-0.5);
      if (!efficacyStopping1[i]) criticalValues[i] = 6.0;
    }
  } else if (asf == "sfof" || asf == "sfp" || asf == "sfkd" ||
    asf == "sfhsd" || asf == "user") {

    // stage 1
    double cumAlphaSpent;
    if (asf == "user") {
      cumAlphaSpent = userAlphaSpending[0];
    } else {
      cumAlphaSpent = errorSpentcpp(s[0], alpha, asf, asfpar);
    }

    if (!efficacyStopping1[0]) {
      criticalValues[0] = 6.0;
    } else {
      criticalValues[0] = R::qnorm(1 - cumAlphaSpent, 0, 1, 1, 0);
    }


    // lambda expression for finding the critical Values at stage k
    int k1=0;
    auto f = [&k1, &cumAlphaSpent, &criticalValues,
              theta, t](double aval)->double {
                NumericVector u(k1+1), l(k1+1);
                for (int i=0; i<k1; i++) {
                  u[i] = criticalValues[i];
                  l[i] = -6.0;
                }
                u[k1] = aval;
                l[k1] = -6.0;

                IntegerVector idx = Range(0,k1);
                List probs = exitprobcpp(u, l, theta[idx], t[idx]);
                double cpu = sum(NumericVector(probs[0]));
                return cpu - cumAlphaSpent;
              };

    // subsequent stages
    for (k1=1; k1<k; k1++) {
      if (asf == "user") {
        cumAlphaSpent = userAlphaSpending[k1];
      } else {
        cumAlphaSpent = errorSpentcpp(s[k1], alpha, asf, asfpar);
      }

      if (!efficacyStopping1[k1]) {
        criticalValues[k1] = 6.0;
      } else {
        if (f(6.0) > 0) { // no alpha spent at current visit
          criticalValues[k1] = 6.0;
        } else {
          criticalValues[k1] = brent(f, -5.0, 6.0, 1.0e-6);
        }
      }
    }
  } else {
    stop("Invalid value for typeAlphaSpending");
  }

  return criticalValues;
}


// [[Rcpp::export]]
List getPower(const double alpha,
              const int kMax,
              const NumericVector& b,
              const NumericVector& theta,
              const NumericVector& I,
              const std::string bsf,
              const double bsfpar,
              const NumericVector& st,
              const LogicalVector& futilityStopping,
              const NumericVector& w) { // w is the sqrt of variance ratio

  double beta;
  NumericVector a(kMax);
  List probs;
  auto f = [kMax, b, futilityStopping, &a,
            bsf, bsfpar, theta, I, st, w](double beta)->double {
              // initialize futility bound to be updated
              a = NumericVector(kMax);
              double eps;

              // first stage
              int k = 0;
              double cb = errorSpentcpp(st[0], beta, bsf, bsfpar);
              if (!futilityStopping[0]) {
                a[0] = -6.0;
              } else {
                eps = R::pnorm(b[0]*w[0] - theta[0]*sqrt(I[0]), 0, 1, 1, 0)
                      - cb;
                if (eps < 0) return -1.0; // to decrease beta
                a[0] = (R::qnorm(cb, 0, 1, 1, 0) + theta[0]*sqrt(I[0]))/w[0];
              }

              // lambda expression for finding futility bound at stage k
              auto g = [&k, &cb, b, &a, theta, I, w](double aval)->double {
                NumericVector u(k+1), l(k+1);
                for (int i=0; i<k; i++) {
                  u[i] = b[i]*w[i];
                  l[i] = a[i]*w[i];
                }
                u[k] = 6.0;
                l[k] = aval*w[k];

                IntegerVector idx = Range(0,k);
                List probs = exitprobcpp(u, l, theta[idx], I[idx]);
                double cpl = sum(NumericVector(probs[1]));
                return cpl - cb;
              };

              for (k=1; k<kMax; k++) {
                cb = errorSpentcpp(st[k], beta, bsf, bsfpar);

                if (!futilityStopping[k]) {
                  a[k] = -6.0;
                } else {
                  eps = g(b[k]);

                  if (g(-6.0) > 0) { // no beta spent at current visit
                    a[k] = -6.0;
                  } else if (eps > 0) {
                    a[k] = brent(g, -6.0, b[k], 1.0e-6);
                  } else if (k < kMax-1) {
                    return -1.0;
                  }
                }
              }

              return eps;
            };

  double v1 = f(0.0001), v2 = f(1-alpha);

  if (v1 == -1.0 || (v1 < 0 && a[kMax-1] == 0)) {
    stop("Power must be less than 0.9999 to use beta spending");
  } else if (v2 > 0) {
    stop("Power must be greater than alpha to use beta spending");
  } else {
    beta = brent(f, 0.0001, 1-alpha, 1.0e-6);
    a[kMax-1] = b[kMax-1];
    probs = exitprobcpp(b*w, a*w, theta, I);
  }

  List result = List::create(
    _["beta"] = beta,
    _["futilityBounds"] = a,
    _["probs"] = probs);

  return result;
}


//' @title Integration with respect to a normal density
//' @description Integrate a function f(theta) with respect to a normal
//' density of theta.
//'
//' @param f Name of the univariate objective function.
//' @param mu The mean of the normal distribution for theta.
//' @param sigma The standard deviation of the normal distribution for theta.
//' @param a One end of the interval bracket.
//' @param b The other end of the interval bracket.
//'
//' @return The value of the integration:
//'   integrate(function(theta) f(theta)*dnorm(theta, mu, sigma), a, b)/
//'   (pnorm(b, mu, sigma) - pnorm(a, mu, sigma)).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @export
//'
// [[Rcpp::plugins(cpp11)]]
double intnorm(const std::function<double(double)>& f,
               double mu, double sigma, double a, double b) {

  int r=18, r1=6*r-1, r2=12*r-3, i, i0, i1=0, i2=r1-1, m=r2, m1=r1;
  double a1=(a-mu)/sigma , b1=(b-mu)/sigma, xlower, xupper, aval;
  NumericVector x1(r1), x(r1), z(r2), w(r2);

  for (i=0; i<r1; i++) {
    if (i < r-1) {
      x1[i] = -3 - 4*log(r/(i+1.0));
    } else if (i < 5*r) {
      x1[i] = -3 + 3*(i+1.0-r)/(2*r);
    } else {
      x1[i] = 3 + 4*log(r/(6*r-i-1.0));
    }
  }

  // trim off x values outside (a1, b1)
  // trim from below
  if (a1 >= x1[0]) {
    i1 = 0;
    while (x1[i1] <= a1) {
      i1++;
    }
    i1--;
    xlower = a1; // lower bound on x
  } else {
    i1 = 0;
    xlower = x1[0];
  }

  // trim from above
  if (b1 <= x1[r1-1]) {
    i2 = r1-1;
    while (x1[i2] >= b1) {
      i2--;
    }
    i2++;
    xupper = b1; // upper bound on x
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


  // integrate
  aval = 0;
  for (i=0; i<m; i++) {
    aval += w[i]*f(mu + sigma*z[i])*R::dnorm(z[i], 0, 1, 0);
  }

  double denom = R::pnorm(b1, 0, 1, 1, 0) - R::pnorm(a1, 0, 1, 1, 0);

  return aval/denom;
}



#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

//' @title Brent's method for minimization
//' @description Using Brent's method, find the abscissa of the minimum of a
//' function known to lie between x1 and x2. Program based on the book -
//' Numerical Recipes in C The Art of Scientific Computing - Second Edition,
//' by William H. Press, Saul A. Teukolsky, William T. Vetterling, and
//' Brian P. Flannery. It mimics the optimize() function in R.
//'
//' @param f Name of the univariate objective function.
//' @param x1 One end of the interval bracket.
//' @param x2 The other end of the interval bracket.
//' @param tol The tolerance limit for stopping the iteration.
//'
//' @return The abscissa x between x1 and x2 such that f(x) = min f(u).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' mini(sin, 0, 2, 0.0001)
//' @export
//'
// [[Rcpp::plugins(cpp11)]]
NumericVector mini(const std::function<double(double)>& f,
                   double x1, double x2, double tol) {
  int iter;
  double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double d=0.0, e=0.0;

  a=x1; b=x2;
  x=w=v=a+CGOLD*(b-a);
  fw=fv=fx=f(x);
  for (iter=0;iter<ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      return NumericVector::create(x,fx);
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      // if the new step size is at least half of the step before last, or
      // if the new point is outside of (a,b) (to the left or to the right)
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      } else{
        d=p/q; // new step size
        u=x+d; // new point
        if (u-a < tol2 || b-u < tol2) {
          d = SIGN(tol1, xm-x);
        }
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=f(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      } else if (fu <= fv || v == x || v == w) {
        v=u;
        fv=fu;
      }
    }
  }
  stop("Too many iterations in mini");
  return NumericVector::create(x, fx); // Never get here
}


NumericVector quad(integr_fn f, void *ex, double lower, double upper,
                   double tol) {
  double epsabs=tol, epsrel=tol, value, abserr;
  int neval, ier, limit=100, lenw=4*limit, last;

  int *iwork = new int[limit];
  double *work = new double[lenw];

  if (!std::isinf(lower) && !std::isinf(upper)) {
    Rdqags(f, ex, &lower, &upper, &epsabs, &epsrel, &value, &abserr,
           &neval, &ier, &limit, &lenw, &last, iwork, work);
  } else if (!std::isinf(lower) && std::isinf(upper)) {
    double bound = lower;
    int inf = 1;
    Rdqagi(f, ex, &bound, &inf, &epsabs, &epsrel, &value, &abserr,
           &neval, &ier, &limit, &lenw, &last, iwork, work);
  } else if (std::isinf(lower) && !std::isinf(upper)) {
    double bound = upper;
    int inf = -1;
    Rdqagi(f, ex, &bound, &inf, &epsabs, &epsrel, &value, &abserr,
           &neval, &ier, &limit, &lenw, &last, iwork, work);
  } else {
    double bound = 0.0;
    int inf = 2;
    Rdqagi(f, ex, &bound, &inf, &epsabs, &epsrel, &value, &abserr,
           &neval, &ier, &limit, &lenw, &last, iwork, work);
  }

  delete[] iwork;
  delete[] work;

  return NumericVector::create(Named("value") = value,
                               Named("abserr") = abserr,
                               Named("neval") = neval,
                               Named("ier") = ier);
}


// Wrapper function for vmmin
List bmini(NumericVector x0, optimfn fn, optimgr gr = nullptr,
           void *ex = nullptr, double eps = 1e-8) {

  int maxit = 100;
  int trace = 0;
  double abstol = eps, reltol = eps;
  int nREPORT = 10;

  int n = static_cast<int>(x0.size());
  double Fmin;
  int fncount = 0, grcount = 0, fail = 0;
  IntegerVector mask(n, 1);  // All parameters are free

  // Convert NumericVector to standard double array
  std::vector<double> x(x0.begin(), x0.end());

  // Call vmmin function
  vmmin(n, x.data(), &Fmin, fn, gr, maxit, trace,
        mask.begin(), abstol, reltol, nREPORT,
        ex, &fncount, &grcount, &fail);

  // Return results as a list
  return List::create(Named("par") = NumericVector(x.begin(), x.end()),
                      Named("value") = Fmin,
                      Named("fncount") = fncount,
                      Named("grcount") = grcount,
                      Named("fail") = fail);
}


//' @title Number of enrolled subjects
//' @description Obtains the number of subjects enrolled by given calendar
//' times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of enrolled subjects.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_accrualDuration
//'
//' @return A vector of total number of subjects enrolled by the
//' specified calendar times.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Example 1: Uniform enrollment with 20 patients per month for 12 months.
//'
//' accrual(time = 3, accrualTime = 0, accrualIntensity = 20,
//'         accrualDuration = 12)
//'
//'
//' # Example 2: Piecewise accrual, 10 patients per month for the first
//' # 3 months, and 20 patients per month thereafter. Patient recruitment
//' # ends at 12 months for the study.
//'
//' accrual(time = c(2, 9), accrualTime = c(0, 3),
//'         accrualIntensity = c(10, 20), accrualDuration = 12)
//'
//' @export
// [[Rcpp::export]]
NumericVector accrual(const NumericVector& time = NA_REAL,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const double accrualDuration = NA_REAL) {

  int i, j, k = static_cast<int>(time.size());
  NumericVector n(k);

  // up to end of enrollment
  NumericVector t = pmax(pmin(time, accrualDuration), 0.0);

  // identify the time interval containing t
  IntegerVector m = pmax(findInterval3(t, accrualTime), 1);

  // sum up patients enrolled in each interval up to t
  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        n[i] += accrualIntensity[j]*(accrualTime[j+1] - accrualTime[j]);
      } else {
        n[i] += accrualIntensity[j]*(t[i] - accrualTime[j]);
      }
    }
  }

  return n;
}


//' @title Accrual duration to enroll target number of subjects
//' @description Obtains the accrual duration to enroll the target number
//' of subjects.
//'
//' @param nsubjects The vector of target number of subjects.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//'
//' @return A vector of accrual durations.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' getAccrualDurationFromN(nsubjects = c(20, 150), accrualTime = c(0, 3),
//'                         accrualIntensity = c(10, 20))
//'
//' @export
// [[Rcpp::export]]
NumericVector getAccrualDurationFromN(
    const NumericVector& nsubjects = NA_REAL,
    const NumericVector& accrualTime = 0,
    const NumericVector& accrualIntensity = NA_REAL) {
  int i, j, I = static_cast<int>(nsubjects.size());
  int J = static_cast<int>(accrualTime.size());
  NumericVector t(I), p(J);

  p[0] = 0;
  for (j=0; j<J-1; j++) {
    p[j+1] = p[j] + accrualIntensity[j]*(accrualTime[j+1] - accrualTime[j]);
  }

  IntegerVector m = findInterval3(nsubjects, p);

  for (i=0; i<I; i++) {
    j = m[i] - 1;
    t[i] = accrualTime[j] + (nsubjects[i] - p[j])/accrualIntensity[j];
  }

  return t;
}


//' @title Probability of being at risk
//' @description Obtains the probability of being at risk at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//'   probability of being at risk.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of being at risk at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' patrisk(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'         lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
NumericVector patrisk(const NumericVector& time = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda = NA_REAL,
                      const NumericVector& gamma = 0) {

  // identify the time interval containing the specified analysis time
  IntegerVector m = pmax(findInterval3(time, piecewiseSurvivalTime), 1);
  int i, j, k = static_cast<int>(time.size());
  int J = static_cast<int>(lambda.size());

  // hazard for failure or dropout
  NumericVector lg(J);
  if (gamma.size() == 1) {
    lg = lambda + gamma[0];
  } else {
    lg = lambda + gamma;
  }

  NumericVector t = piecewiseSurvivalTime;

  // sum up cumulative hazard up to time
  NumericVector a(k);
  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        a[i] += lg[j]*(t[j+1] - t[j]);
      } else {
        a[i] += lg[j]*(time[i] - t[j]);
      }
    }
  }

  return exp(-a);
}


//' @title Probability of having an event
//' @description Obtains the probability of having an event at given analysis
//' times.
//'
//' @param time A vector of analysis times at which to calculate the
//'   probability of having an event.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of probabilities of having an event at the specified
//' analysis times after enrollment for a patient in a treatment group with
//' specified piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' pevent(time = c(3, 9), piecewiseSurvivalTime = c(0, 6),
//'        lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
NumericVector pevent(const NumericVector& time = NA_REAL,
                     const NumericVector& piecewiseSurvivalTime = 0,
                     const NumericVector& lambda = NA_REAL,
                     const NumericVector& gamma = 0) {

  // identify the time interval containing the specified analysis time
  IntegerVector m = pmax(findInterval3(time, piecewiseSurvivalTime), 1);
  int i, j, k = static_cast<int>(time.size());
  int J = static_cast<int>(lambda.size());

  // hazard for failure or dropout
  NumericVector lg(J);
  if (gamma.size() == 1) {
    lg = lambda + gamma[0];
  } else {
    lg = lambda + gamma;
  }

  // sum up cumulative hazard up to time
  NumericVector t = piecewiseSurvivalTime;
  NumericVector n = patrisk(t, t, lambda, gamma);
  NumericVector a(k);
  double p;

  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        p = lambda[j]/lg[j]*(1 - exp(-lg[j]*(t[j+1] - t[j])));
      } else {
        p = lambda[j]/lg[j]*(1 - exp(-lg[j]*(time[i] - t[j])));
      }
      a[i] += n[j]*p;
    }
  }

  return a;
}


//' @title Integrated event probability over an interval with constant hazard
//' @description Obtains the integration probability of having an event
//' during an interval with constant hazard.
//'
//' @param j The analysis time interval with constant hazard.
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A value for the integrated probability of having an event
//' during an interval with constant hazard for a treatment
//' group with specified piecewise exponential survival and dropout
//' distributions.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' hd(j = 1, t1 = 1, t2 = 3, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double hd(const int j = NA_INTEGER,
          const double t1 = NA_REAL,
          const double t2 = NA_REAL,
          const NumericVector& piecewiseSurvivalTime = 0,
          const NumericVector& lambda = NA_REAL,
          const NumericVector& gamma = 0) {

  int j1 = j-1;

  // lower bound of time interval j for piecewise exponential distribution
  NumericVector t0 = NumericVector::create(piecewiseSurvivalTime[j1]);

  // probability of being at risk at the start of interval j
  NumericVector n0 = patrisk(t0, piecewiseSurvivalTime, lambda, gamma);

  // probability of having an event at the start of interval j
  NumericVector d0 = pevent(t0, piecewiseSurvivalTime, lambda, gamma);


  int J = static_cast<int>(lambda.size());

  // hazard for failure or dropout
  NumericVector lg(J);
  if (gamma.size() == 1) {
    lg = lambda + gamma[0];
  } else {
    lg = lambda + gamma;
  }

  // integration of conditional probability of having an event over (t1,t2)
  // given survival at the start of interval j
  double q1 = (exp(-lg[j1]*(t1-t0[0])) - exp(-lg[j1]*(t2-t0[0])))/lg[j1];
  double q = lambda[j1]/lg[j1] * (t2-t1 - q1);

  // sum up the integration for the already failed and to-be-failed
  return d0[0]*(t2-t1) + n0[0]*q;
}


//' @title Integrated event probability over an interval
//' @description Obtains the integration of the probability of having an
//' event during an interval. The specified analysis time interval can span
//' more than one analysis time interval with constant hazard.
//'
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A value for the integrated probability of having an event
//' during an interval for a treatment group with specified
//' piecewise exponential survival and dropout distributions.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survival with hazard 0.0533 in the first 6
//' # months, and hazard 0.0309 thereafter, and 5% dropout by the end of
//' # 1 year.
//'
//' pd(t1 = 1, t2 = 8, piecewiseSurvivalTime = c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
double pd(const double t1 = NA_REAL,
          const double t2 = NA_REAL,
          const NumericVector& piecewiseSurvivalTime = 0,
          const NumericVector& lambda = NA_REAL,
          const NumericVector& gamma = 0) {

  // identify the analysis time intervals containing t1 and t2
  NumericVector t12 = NumericVector::create(t1, t2);
  IntegerVector j12 = pmax(findInterval3(t12, piecewiseSurvivalTime), 1) - 1;

  NumericVector t = piecewiseSurvivalTime;

  int j, j1=j12[0], j2=j12[1];

  // sum up the integrated event probabilities across analysis time intervals
  double a=0, x;
  for (j=j1; j<=j2; j++) {
    if (j1==j2) {
      x = hd(j+1, t1, t2, t, lambda, gamma);
    } else if (j==j1) {
      x = hd(j+1, t1, t[j+1], t, lambda, gamma);
    } else if (j==j2) {
      x = hd(j+1, t[j], t2, t, lambda, gamma);
    } else {
      x = hd(j+1, t[j], t[j+1], t, lambda, gamma);
    }
    a += x;
  }

  return a;
}


//' @title Number of patients enrolled during an interval and having an event
//' by specified calendar times
//' @description Obtains the number of patients who are enrolled during a
//' specified enrollment time interval and have an event by the specified
//' calendar times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of patients having an event.
//' @param u1 Lower bound of the accrual time interval.
//' @param u2 Upper bound of the accrual time interval.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//' @inheritParams param_gamma
//'
//' @return A vector of number of patients who are enrolled during a
//' specified enrollment time interval and have an event by the specified
//' calendar times for a given treatment group had the enrollment being
//' restricted to the treatment group. By definition, we must have
//' \code{time >= u2}.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, 10 patients per month for the first 3 months, and
//' # 20 patients per month thereafter. Piecewise exponential survival with
//' # hazard 0.0533 in the first 6 months, and hazard 0.0309 thereafter,
//' # and 5% dropout by the end of 1 year.
//'
//' ad(time = c(9, 15), u1 = 1, u2 = 8, accrualTime = c(0, 3),
//'    accrualIntensity = c(10, 20), piecewiseSurvivalTime=c(0, 6),
//'    lambda = c(0.0533, 0.0309), gamma = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
NumericVector ad(const NumericVector& time = NA_REAL,
                 const double u1 = NA_REAL,
                 const double u2 = NA_REAL,
                 const NumericVector& accrualTime = 0,
                 const NumericVector& accrualIntensity = NA_REAL,
                 const NumericVector& piecewiseSurvivalTime = 0,
                 const NumericVector& lambda = NA_REAL,
                 const NumericVector& gamma = 0) {

  // identify the accrual time intervals containing u1 and u2
  NumericVector u12 = NumericVector::create(u1, u2);
  IntegerVector j12 = pmax(findInterval3(u12, accrualTime), 1) - 1;

  NumericVector u = accrualTime;

  int i, j, j1=j12[0], j2=j12[1], k=static_cast<int>(time.size());

  NumericVector a(k);

  // sum up the number of patients with event across accrual time intervals
  double t, x;
  for (i=0; i<k; i++) {
    t = time[i];
    for (j=j1; j<=j2; j++) {
      if (j1==j2) {
        x = pd(t-u2, t-u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j==j1) {
        x = pd(t-u[j+1], t-u1, piecewiseSurvivalTime, lambda, gamma);
      } else if (j==j2) {
        x = pd(t-u2, t-u[j], piecewiseSurvivalTime, lambda, gamma);
      } else {
        x = pd(t-u[j+1], t-u[j], piecewiseSurvivalTime, lambda, gamma);
      }
      a[i] += accrualIntensity[j]*x;
    }
  }

  return a;
}


//' @title Number of subjects at risk
//' @description Obtains the number of subjects at risk at given analysis
//' times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
//'   of patients at risk.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients at risk at the specified
//' analysis times (row) for each treatment group (column).
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' natrisk(time = c(9, 24), allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, minFollowupTime = 18,
//'         maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix natrisk(const NumericVector& time = NA_REAL,
                      const double allocationRatioPlanned = 1,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda1 = NA_REAL,
                      const NumericVector& lambda2 = NA_REAL,
                      const NumericVector& gamma1 = 0,
                      const NumericVector& gamma2 = 0,
                      const double accrualDuration = NA_REAL,
                      const double minFollowupTime = NA_REAL,
                      const double maxFollowupTime = NA_REAL) {

  // truncate the analysis time by the maximum follow-up
  NumericVector t = pmin(time, maxFollowupTime);

  // enrollment time
  NumericVector u = pmin(accrualDuration+minFollowupTime-t, accrualDuration);

  // number of patients enrolled
  NumericVector a = accrual(u, accrualTime, accrualIntensity,
                            accrualDuration);

  // probability of being randomized to the active treatment group
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);

  // number of patients at risk in each treatment group
  int k = static_cast<int>(time.size());
  NumericMatrix n(k, 2);
  n(_, 0) = phi*a*patrisk(t, piecewiseSurvivalTime, lambda1, gamma1);
  n(_, 1) = (1-phi)*a*patrisk(t, piecewiseSurvivalTime, lambda2, gamma2);

  return n;
}


//' @title Number of subjects having an event
//' @description Obtains the number of subjects having an event by given
//' analysis times for each treatment group.
//'
//' @param time A vector of analysis times at which to calculate the number
//'   of patients having an event.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients having an event at the
//' specified analysis times (row) for each treatment group (column).
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' nevent(time = c(9, 24), allocationRatioPlanned = 1,
//'        accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'        piecewiseSurvivalTime = c(0, 6),
//'        lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'        gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'        accrualDuration = 12, minFollowupTime = 18,
//'        maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix nevent(const NumericVector& time = NA_REAL,
                     const double allocationRatioPlanned = 1,
                     const NumericVector& accrualTime = 0,
                     const NumericVector& accrualIntensity = NA_REAL,
                     const NumericVector& piecewiseSurvivalTime = 0,
                     const NumericVector& lambda1 = NA_REAL,
                     const NumericVector& lambda2 = NA_REAL,
                     const NumericVector& gamma1 = 0,
                     const NumericVector& gamma2 = 0,
                     const double accrualDuration = NA_REAL,
                     const double minFollowupTime = NA_REAL,
                     const double maxFollowupTime = NA_REAL) {

  // truncate the analysis time by the maximum follow-up
  NumericVector t = pmin(time, maxFollowupTime);

  // enrollment time
  NumericVector u = pmin(accrualDuration+minFollowupTime-t, accrualDuration);

  // number of patients enrolled
  NumericVector a = accrual(u, accrualTime, accrualIntensity,
                            accrualDuration);

  // probability of being randomized to the active treatment group
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);

  // number of patients having an event in each treatment group
  NumericVector u1(1);
  u1[0] = accrualDuration + minFollowupTime;

  int i, k = static_cast<int>(time.size());
  NumericMatrix d(k, 2);

  NumericVector d1(k), d2(k);
  d1 = a*pevent(t, piecewiseSurvivalTime, lambda1, gamma1);
  d2 = a*pevent(t, piecewiseSurvivalTime, lambda2, gamma2);

  for (i=0; i<k; i++) {
    d(i,0) = phi*(d1[i] + ad(u1, u[i], accrualDuration, accrualTime,
                  accrualIntensity, piecewiseSurvivalTime,
                  lambda1, gamma1)[0]);
    d(i,1) = (1-phi)*(d2[i] + ad(u1, u[i], accrualDuration, accrualTime,
              accrualIntensity, piecewiseSurvivalTime, lambda2, gamma2)[0]);
  }

  return d;
}


//' @title Number of subjects having an event by calendar time
//' @description Obtains the number of subjects having an event by given
//' calendar times for each treatment group.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of patients having an event.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_minFollowupTime
//' @inheritParams param_maxFollowupTime
//'
//' @return A matrix of the number of patients having an event at the
//' specified calendar times (row) for each treatment group (column).
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//' nevent2(time = c(9, 24), allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, minFollowupTime = 18,
//'         maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix nevent2(const NumericVector& time = NA_REAL,
                      const double allocationRatioPlanned = 1,
                      const NumericVector& accrualTime = 0,
                      const NumericVector& accrualIntensity = NA_REAL,
                      const NumericVector& piecewiseSurvivalTime = 0,
                      const NumericVector& lambda1 = NA_REAL,
                      const NumericVector& lambda2 = NA_REAL,
                      const NumericVector& gamma1 = 0,
                      const NumericVector& gamma2 = 0,
                      const double accrualDuration = NA_REAL,
                      const double minFollowupTime = NA_REAL,
                      const double maxFollowupTime = NA_REAL) {

  // truncate the calendar time by study end
  NumericVector t = pmin(time, accrualDuration + minFollowupTime);

  // enrollment time
  NumericVector u = pmin(pmax(t - maxFollowupTime, 0.0), accrualDuration);
  NumericVector w = pmin(t, accrualDuration);

  // number of patients enrolled
  NumericVector a = accrual(u, accrualTime, accrualIntensity,
                            accrualDuration);

  // probability of being randomized to the active treatment group
  double phi = allocationRatioPlanned/(1+allocationRatioPlanned);

  // number of patients having an event in each treatment group
  NumericVector s(1), v(1);
  s[0] = maxFollowupTime;

  int i, k = static_cast<int>(time.size());
  NumericMatrix d(k, 2);

  NumericVector d1(k), d2(k);
  d1 = a*pevent(s, piecewiseSurvivalTime, lambda1, gamma1)[0];
  d2 = a*pevent(s, piecewiseSurvivalTime, lambda2, gamma2)[0];

  for (i=0; i<k; i++) {
    v[0] = t[i];
    d(i,0) = phi*(d1[i] + ad(v, u[i], w[i], accrualTime, accrualIntensity,
                  piecewiseSurvivalTime, lambda1, gamma1)[0]);
    d(i,1) = (1-phi)*(d2[i] + ad(v, u[i], w[i], accrualTime,
              accrualIntensity, piecewiseSurvivalTime, lambda2, gamma2)[0]);
  }

  return d;
}


//' @title Power and sample size for a generic group sequential design
//' @description Obtains the maximum information and stopping boundaries
//' for a generic group sequential design assuming a constant treatment
//' effect, or obtains the power given the maximum information and
//' stopping boundaries.
//'
//' @param beta The type II error.
//' @param IMax The maximum information. Either \code{beta} or \code{IMax}
//'   should be provided while the other one should be missing.
//' @param theta The parameter value.
//' @inheritParams param_kMax
//' @param informationRates The information rates. Fixed prior to the trial.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param varianceRatio The ratio of the variance under H0 to the
//'   variance under H1.
//'
//' @return An S3 class \code{design} object with three components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{attainedAlpha}: The attained significance level, which is
//'       different from the overall significance level in the presence of
//'       futility stopping.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{theta}: The parameter value.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedInformationH1}: The expected information under H1.
//'
//'     - \code{expectedInformationH0}: The expected information under H0.
//'
//'     - \code{drift}: The drift parameter, equal to
//'       \code{theta*sqrt(information)}.
//'
//'     - \code{inflationFactor}: The inflation factor (relative to the
//'       fixed design).
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'
//'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{futilityPerStage}: The probability for futility stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeFutility}: The cumulative probability for futility
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'
//'     - \code{efficacyTheta}: The efficacy boundaries on the parameter
//'       scale.
//'
//'     - \code{futilityTheta}: The futility boundaries on the parameter
//'       scale.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{futilityP}: The futility boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'
//'     - \code{futilityStopping}: Whether to allow futility stopping.
//'
//'     - \code{rejectPerStageH0}: The probability for efficacy stopping
//'       under H0.
//'
//'     - \code{futilityPerStageH0}: The probability for futility stopping
//'       under H0.
//'
//'     - \code{cumulativeRejectionH0}: The cumulative probability for
//'       efficacy stopping under H0.
//'
//'     - \code{cumulativeFutilityH0}: The cumulative probability for
//'       futility stopping under H0.
//'
//' * \code{settings}: A list containing the following input parameters:
//'
//'     - \code{typeAlphaSpending}: The type of alpha spending.
//'
//'     - \code{parameterAlphaSpending}: The parameter value for alpha
//'       spending.
//'
//'     - \code{userAlphaSpending}: The user defined alpha spending.
//'
//'     - \code{typeBetaSpending}: The type of beta spending.
//'
//'     - \code{parameterBetaSpending}: The parameter value for beta
//'       spending.
//'
//'     - \code{userBetaSpending}: The user defined beta spending.
//'
//'     - \code{spendingTime}: The error spending time at each analysis.
//'
//'     - \code{varianceRatio}: The ratio of the variance under H0
//'       to the variance under H1.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Christopher Jennison, Bruce W. Turnbull. Group Sequential Methods with
//' Applications to Clinical Trials. Chapman & Hall/CRC: Boca Raton, 2000,
//' ISBN:0849303168
//'
//' @examples
//'
//' # Example 1: obtain the maximum information given power
//' (design1 <- getDesign(
//'   beta = 0.2, theta = -log(0.7),
//'   kMax = 2, informationRates = c(0.5,1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   typeBetaSpending = "sfP"))
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesign(
//'   IMax = 72.5, theta = -log(0.7),
//'   kMax = 3, informationRates = c(0.5, 0.75, 1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   typeBetaSpending = "sfP"))
//'
//' @export
// [[Rcpp::export]]
List getDesign(const double beta = NA_REAL,
               const double IMax = NA_REAL,
               const double theta = NA_REAL,
               const int kMax = 1,
               const NumericVector& informationRates = NA_REAL,
               const LogicalVector& efficacyStopping = NA_LOGICAL,
               const LogicalVector& futilityStopping = NA_LOGICAL,
               const NumericVector& criticalValues = NA_REAL,
               const double alpha = 0.025,
               const std::string typeAlphaSpending = "sfOF",
               const double parameterAlphaSpending = NA_REAL,
               const NumericVector& userAlphaSpending = NA_REAL,
               const NumericVector& futilityBounds = NA_REAL,
               const std::string typeBetaSpending = "none",
               const double parameterBetaSpending = NA_REAL,
               const NumericVector& userBetaSpending = NA_REAL,
               const NumericVector& spendingTime = NA_REAL,
               const double varianceRatio = 1) {

  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector spendingTime1 = clone(spendingTime);

  double alpha1 = alpha;
  double beta1 = beta;
  double IMax1 = IMax;
  double drift, inflationFactor;

  std::string unknown;

  if (R_isnancpp(beta) && R_isnancpp(IMax)) {
    stop("beta and IMax cannot be both missing");
  }

  if (!R_isnancpp(beta) && !R_isnancpp(IMax)) {
    stop("Only one of beta and IMax should be provided");
  }

  if (!R_isnancpp(IMax)) {
    if (IMax <= 0) {
      stop("IMax must be positive");
    }
    unknown = "beta";
  } else if (!R_isnancpp(beta)) {
    unknown = "IMax";
  }

  if (R_isnancpp(theta)) {
    stop("theta must be provided");
  }

  if (kMax == NA_INTEGER) {
    stop("kMax must be provided");
  }

  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = NumericVector(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
    }
  } else {
    futilityStopping1 = rep(1, kMax);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if ((unknown == "IMax") && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }


  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  }

  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }

    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      stop("futilityBounds and criticalValues must meet at final analysis");
    }
  }

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  if (unknown == "IMax") {
    if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
        bsf=="sfkd" || bsf=="sfhsd" || bsf=="user" || bsf=="none")) {
      stop("Invalid value for typeBetaSpending");
    }
  } else {
    if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
        bsf=="sfkd" || bsf=="sfhsd" || bsf=="none")) {
      stop("Invalid value for typeBetaSpending");
    }
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (unknown=="IMax" && bsf=="user") {
    if (is_true(any(is_na(userBetaSpending)))) {
      stop("userBetaSpending must be specified");
    } else if (userBetaSpending.size() < kMax) {
      stop("Insufficient length of userBetaSpending");
    } else if (userBetaSpending[0] < 0) {
      stop("Elements of userBetaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userBetaSpending) < 0))) {
      stop("Elements of userBetaSpending must be nondecreasing");
    } else if (userBetaSpending[kMax-1] != beta) {
      stop("userBetaSpending must end with specified beta");
    }
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (varianceRatio <= 0) {
    stop("varianceRatio must be positive");
  }


  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1, efficacyStopping1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, l, zero, informationRates1);
  alpha1 = sum(NumericVector(probs[0]));

  bool missingFutilityBounds = is_true(any(is_na(futilityBounds)));
  if (kMax > 1) {
    if (missingFutilityBounds && bsf=="none") {
      futilityBounds1 = rep(-6.0, kMax);
      futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    } else if (!missingFutilityBounds &&
      futilityBounds1.size() == kMax-1) {
      futilityBounds1.push_back(criticalValues1[kMax-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futilityBounds1 = criticalValues1[kMax-1];
    }
  }

  // multiplier for the boundaries under the alternative hypothesis
  NumericVector w = rep(sqrt(varianceRatio), kMax);

  NumericVector t = informationRates1;
  NumericVector st = spendingTime1;
  NumericVector theta1(kMax);
  if (unknown == "IMax") {
    auto f = [beta, kMax, t, futilityStopping1,
              criticalValues1, &futilityBounds1,
              bsf, bsfpar, userBetaSpending, st, w,
              missingFutilityBounds](double aval)->double {

                NumericVector theta1 = rep(aval, kMax);

                // compute stagewise exit probabilities
                if (!missingFutilityBounds || bsf=="none" || kMax==1) {
                  List probs = exitprobcpp(
                    criticalValues1*w, futilityBounds1*w, theta1, t);
                  NumericVector pu = NumericVector(probs[0]);
                  double overallReject = sum(pu);
                  return overallReject - (1-beta);
                } else {
                  // initialize futility bound to be updated
                  futilityBounds1 = NumericVector(kMax);
                  double epsilon;

                  // first stage
                  int k = 0;
                  double cumBetaSpent;
                  if (bsf=="user") {
                    cumBetaSpent = userBetaSpending[0];
                  } else {
                    cumBetaSpent = errorSpentcpp(st[0], beta, bsf, bsfpar);
                  }

                  if (!futilityStopping1[0]) {
                    futilityBounds1[0] = -6.0;
                  } else {
                    epsilon = R::pnorm(criticalValues1[0]*w[0] -
                      theta1[0]*sqrt(t[0]), 0, 1, 1, 0) - cumBetaSpent;
                    if (epsilon < 0) return -1.0;
                    futilityBounds1[0] = (R::qnorm(cumBetaSpent, 0, 1, 1, 0)
                                            + theta1[0]*sqrt(t[0]))/w[0];
                  }


                  // lambda expression for finding futility bound at stage k
                  auto g = [&k, &cumBetaSpent, criticalValues1,
                            &futilityBounds1, theta1, w,
                            t](double aval)->double {
                              NumericVector u(k+1), l(k+1);
                              for (int i=0; i<k; i++) {
                                u[i] = criticalValues1[i]*w[i];
                                l[i] = futilityBounds1[i]*w[i];
                              }
                              u[k] = 6.0;
                              l[k] = aval*w[k];

                              IntegerVector idx = Range(0,k);
                              List probs = exitprobcpp(
                                u, l, theta1[idx], t[idx]);
                              double cpl = sum(NumericVector(probs[1]));
                              return cpl - cumBetaSpent;
                            };


                  for (k=1; k<kMax; k++) {
                    if (bsf == "user") {
                      cumBetaSpent = userBetaSpending[k];
                    } else {
                      cumBetaSpent = errorSpentcpp(st[k], beta, bsf, bsfpar);
                    }

                    if (!futilityStopping1[k]) {
                      futilityBounds1[k] = -6.0;
                    } else {
                      epsilon = g(criticalValues1[k]);

                      if (g(-6.0) > 0) { // no beta spent at current visit
                        futilityBounds1[k] = -6.0;
                      } else if (epsilon > 0) {
                        futilityBounds1[k] = brent(
                          g, -6.0, criticalValues1[k], 1.0e-6);
                      } else if (k < kMax-1) {
                        return -1.0;
                      }
                    }
                  }

                  return epsilon;
                }
              };

    drift = brent(f, 0.0, 6.0, 1.0e-6);
    IMax1 = pow(drift/theta, 2);
    futilityBounds1[kMax-1] = criticalValues1[kMax-1];
    theta1 = rep(drift, kMax);
    probs = exitprobcpp(criticalValues1*w, futilityBounds1*w, theta1, t);
  } else {
    drift = theta*sqrt(IMax1);
    theta1 = rep(drift, kMax);

    if (!missingFutilityBounds || bsf=="none" || kMax==1) {
      probs = exitprobcpp(criticalValues1*w, futilityBounds1*w, theta1, t);
      beta1 = 1 - sum(NumericVector(probs[0]));
    } else {
      List out = getPower(alpha1, kMax, criticalValues1, theta1, t,
                          bsf, bsfpar, st, futilityStopping1, w);

      beta1 = out[0];
      futilityBounds1 = out[1];
      probs = out[2];
    }
  }

  double driftf = R::qnorm(1-alpha1, 0, 1, 1, 0)*w[0] +
    R::qnorm(1-beta1, 0, 1, 1, 0);
  inflationFactor = pow(drift/driftf, 2);


  // output the results
  NumericVector information(kMax);
  NumericVector efficacyTheta(kMax);
  NumericVector futilityTheta(kMax);
  NumericVector efficacyP(kMax);
  NumericVector futilityP(kMax);
  for (int i=0; i<kMax; i++) {
    information[i] = IMax1*informationRates1[i];
    efficacyTheta[i] = criticalValues1[i]/sqrt(information[i])*w[i];
    futilityTheta[i] = futilityBounds1[i]/sqrt(information[i])*w[i];
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
    futilityP[i] = 1 - R::pnorm(futilityBounds1[i], 0, 1, 1, 0);
  }

  // stagewise exit probabilities under H1
  NumericVector pu(kMax), pl(kMax), ptotal(kMax);
  pu = NumericVector(probs[0]);
  pl = NumericVector(probs[1]);
  ptotal = pu + pl;

  double expectedInformationH1 = sum(ptotal*information);

  double overallReject = sum(pu);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);

  // cumulative alpha spent under H0 with non-binding futility
  NumericVector futilityBounds0(kMax, -6.0), theta0(kMax);
  List probs0 = exitprobcpp(criticalValues1, futilityBounds0, theta0, t);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs0[0]));

  // stagewise exit probabilities under H0 with binding futility
  probs0 = exitprobcpp(criticalValues1, futilityBounds1, theta0, t);
  NumericVector pu0(kMax), pl0(kMax), ptotal0(kMax);
  pu0 = NumericVector(probs0[0]);
  pl0 = NumericVector(probs0[1]);
  ptotal0 = pu0 + pl0;

  double expectedInformationH0 = sum(ptotal0*information);

  double overallRejectH0 = sum(pu0);
  NumericVector cpu0 = cumsum(pu0);
  NumericVector cpl0 = cumsum(pl0);

  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      efficacyStopping1[i] = 0;
    }
    if (futilityBounds1[i] == -6) {
      futilityStopping1[i] = 0;
    }
  }

  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["futilityBounds"] = futilityBounds1,
    _["rejectPerStage"] = pu,
    _["futilityPerStage"] = pl,
    _["cumulativeRejection"] = cpu,
    _["cumulativeFutility"] = cpl,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["efficacyTheta"] = efficacyTheta,
    _["futilityTheta"] = futilityTheta,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["information"] = information,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1,
    _["rejectPerStageH0"] = pu0,
    _["futilityPerStageH0"] = pl0,
    _["cumulativeRejectionH0"] = cpu0,
    _["cumulativeFutilityH0"] = cpl0);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["attainedAlpha"] = overallRejectH0,
    _["kMax"] = kMax,
    _["theta"] = theta,
    _["information"] = IMax1,
    _["expectedInformationH1"] = expectedInformationH1,
    _["expectedInformationH0"] = expectedInformationH0,
    _["drift"] = drift,
    _["inflationFactor"] = inflationFactor);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["typeBetaSpending"] = typeBetaSpending,
    _["parameterBetaSpending"] = parameterBetaSpending,
    _["userBetaSpending"] = userBetaSpending,
    _["spendingTime"] = spendingTime,
    _["varianceRatio"] = varianceRatio);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "design";

  return result;
}


//' @title Power and sample size for a generic group sequential equivalence
//' design
//'
//' @description Obtains the maximum information and stopping boundaries
//' for a generic group sequential equivalence design assuming a constant
//' treatment effect, or obtains the power given the maximum information
//' and stopping boundaries.
//'
//' @param beta The type II error.
//' @param IMax The maximum information. Either \code{beta} or \code{IMax}
//'   should be provided while the other one should be missing.
//' @param thetaLower The parameter value at the lower equivalence limit.
//' @param thetaUpper The parameter value at the upper equivalence limit.
//' @param theta The parameter value under the alternative hypothesis.
//' @inheritParams param_kMax
//' @param informationRates The information rates. Fixed prior to the trial.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_criticalValues
//' @param alpha The significance level for each of the two one-sided
//'   tests, e.g., 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param varianceRatioH10 The ratio of the variance under H10 to
//'   the variance under H1.
//' @param varianceRatioH20 The ratio of the variance under H20 to
//'   the variance under H1.
//' @param varianceRatioH12 The ratio of the variance under H10 to
//'   the variance under H20.
//' @param varianceRatioH21 The ratio of the variance under H20 to
//'   the variance under H10.
//'
//' @details
//' Consider the equivalence design with two one-sided hypotheses:
//' \deqn{H_{10}: \theta \leq \theta_{10},}
//' \deqn{H_{20}: \theta \geq \theta_{20}.}
//' We reject \eqn{H_{10}} at or before look \eqn{k} if
//' \deqn{Z_{1j} = (\hat{\theta}_j - \theta_{10})\sqrt{\frac{n_j}{v_{10}}}
//' \geq b_j}
//' for some \eqn{j=1,\ldots,k}, where \eqn{\{b_j:j=1,\ldots,K\}} are the
//' critical values associated with the specified alpha-spending function,
//' and \eqn{v_{10}} is the null variance of
//' \eqn{\hat{\theta}} based on the restricted maximum likelihood (reml)
//' estimate of model parameters subject to the constraint imposed by
//' \eqn{H_{10}} for one sampling unit drawn from \eqn{H_1}. For example,
//' for estimating the risk difference \eqn{\theta = \pi_1 - \pi_2},
//' the asymptotic limits of the
//' reml estimates of \eqn{\pi_1} and \eqn{\pi_2} subject to the constraint
//' imposed by \eqn{H_{10}} are given by
//' \deqn{(\tilde{\pi}_1, \tilde{\pi}_2) = f(\theta_{10}, r, r\pi_1,
//' 1-r, (1-r)\pi_2),}
//' where \eqn{f(\theta_0, n_1, y_1, n_2, y_2)} is the function to obtain
//' the reml of \eqn{\pi_1} and \eqn{\pi_2} subject to the constraint that
//' \eqn{\pi_1-\pi_2 = \theta_0} with observed data
//' \eqn{(n_1, y_1, n_2, y_2)} for the number of subjects and number of
//' responses in the active treatment and control groups,
//' \eqn{r} is the randomization probability for the active treatment
//' group, and \deqn{v_{10} = \frac{\tilde{\pi}_1 (1-\tilde{\pi}_1)}{r} +
//' \frac{\tilde{\pi}_2 (1-\tilde{\pi}_2)}{1-r}.}
//'
//' Let \eqn{I_j = n_j/v_1} denote the information for \eqn{\theta} at the
//' \eqn{j}th look, where
//' \deqn{v_{1} = \frac{\pi_1 (1-\pi_1)}{r} + \frac{\pi_2 (1-\pi_2)}{1-r}}
//' denotes the variance of \eqn{\hat{\theta}} under \eqn{H_1} for one
//' sampling unit. It follows that
//' \deqn{(Z_{1j} \geq b_j) = (Z_j \geq w_{10} b_j +
//' (\theta_{10}-\theta)\sqrt{I_j}),}
//' where \eqn{Z_j = (\hat{\theta}_j - \theta)\sqrt{I_j}}, and
//' \eqn{w_{10} = \sqrt{v_{10}/v_1}}.
//'
//' Similarly, we reject \eqn{H_{20}} at or before look \eqn{k} if
//' \deqn{Z_{2j} = (\hat{\theta}_j - \theta_{20})\sqrt{\frac{n_j}{v_{20}}}
//' \leq -b_j} for some \eqn{j=1,\ldots,k}, where \eqn{v_{20}} is the null
//' variance of \eqn{\hat{\theta}} based on the reml estimate of model
//' parameters subject to the constraint imposed by \eqn{H_{20}} for
//' one sampling unit drawn from \eqn{H_1}. We have
//' \deqn{(Z_{2j} \leq -b_j) = (Z_j \leq -w_{20} b_j +
//' (\theta_{20}-\theta)\sqrt{I_j}),}
//' where \eqn{w_{20} = \sqrt{v_{20}/v_1}}.
//'
//' Let \eqn{l_j = w_{10}b_j + (\theta_{10}-\theta)\sqrt{I_j}},
//' and \eqn{u_j = -w_{20}b_j + (\theta_{20}-\theta)\sqrt{I_j}}.
//' The cumulative probability to reject \eqn{H_0 = H_{10} \cup H_{20}} at
//' or before look \eqn{k} under the alternative hypothesis \eqn{H_1} is
//' given by
//' \deqn{P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j) \cap
//' \cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right) = p_1 + p_2 + p_{12},}
//' where
//' \deqn{p_1 = P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \geq l_j)\right),}
//' \deqn{p_2 = P_\theta\left(\cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \leq u_j)\right),}
//' and
//' \deqn{p_{12} = P_\theta\left(\cup_{j=1}^{k} \{(Z_j \geq l_j) \cup
//' (Z_j \leq u_j)\}\right).}
//' Of note, both \eqn{p_1} and \eqn{p_2} can be evaluated using
//' one-sided exit probabilities for group sequential designs.
//' If there exists \eqn{j\leq k} such that \eqn{l_j \leq u_j}, then
//' \eqn{p_{12} = 1}. Otherwise, \eqn{p_{12}} can be evaluated using
//' two-sided exit probabilities for group sequential designs.
//'
//' To evaluate the type I error of the equivalence trial under
//' \eqn{H_{10}}, we first match the information under \eqn{H_{10}}
//' with the information under \eqn{H_1}. For example, for estimating
//' the risk difference for two independent samples, the sample size
//' \eqn{n_{10}} under \eqn{H_{10}} must satisfy
//' \deqn{\frac{1}{n_{10}}\left(\frac{(\pi_2 + \theta_{10})
//' (1 - \pi_2 - \theta_{10})}{r} + \frac{\pi_2 (1-\pi_2)}{1-r}\right)
//' = \frac{1}{n}\left(\frac{\pi_1(1-\pi_1)}{r} +
//' \frac{\pi_2 (1-\pi_2)}{1-r}\right).}
//' Then we obtain the reml estimates of \eqn{\pi_1} and \eqn{\pi_2}
//' subject to the constraint imposed by \eqn{H_{20}} for one sampling
//' unit drawn from \eqn{H_{10}},
//' \deqn{(\tilde{\pi}_{10}, \tilde{\pi}_{20}) = f(\theta_{20}, r,
//' r(\pi_2 + \theta_{10}), 1-r, (1-r)\pi_2).}
//' Let \eqn{t_j} denote the information fraction at look \eqn{j}.
//' Define \deqn{\tilde{v}_1 = \frac{(\pi_2 + \theta_{10})
//' (1-\pi_2 -\theta_{10})}{r} + \frac{\pi_2 (1-\pi_2)}{1-r},} and
//' \deqn{\tilde{v}_{20} = \frac{\tilde{\pi}_{10}(1-\tilde{\pi}_{10})}{r} +
//' \frac{\tilde{\pi}_{20} (1-\tilde{\pi}_{20})}{1-r}.}
//'
//' The cumulative rejection probability under \eqn{H_{10}} at or before
//' look \eqn{k} is given by
//' \deqn{P_{\theta_{10}}\left(\cup_{j=1}^{k} \{(\hat{\theta}_j - \theta_{10})
//' \sqrt{n_{10} t_j/\tilde{v}_1} \geq b_j\} \cap
//' \cup_{j=1}^{k} \{(\hat{\theta}_j - \theta_{20})
//' \sqrt{n_{10} t_j/\tilde{v}_{20}} \leq -b_j\}\right) =
//' q_1 + q_2 + q_{12},}
//' where
//' \deqn{q_1 = P_{\theta_{10}}\left(\cup_{j=1}^{k}
//' \{(\hat{\theta}_j - \theta_{10})
//' \sqrt{n_{10} t_j/\tilde{v}_1} \geq b_j\}\right) =
//' P_{\theta_{10}}\left(\cup_{j=1}^{k} (Z_j \geq b_j)\right),}
//' \deqn{q_2 = P_{\theta_{10}}\left(\cup_{j=1}^{k}
//' \{(\hat{\theta}_j - \theta_{20})
//' \sqrt{n_{10} t_j/\tilde{v}_{20}} \leq -b_j\}\right) =
//' P_{\theta_{10}}\left(\cup_{j=1}^{k} (Z_j \leq -b_j w_{21} +
//' (\theta_{20} - \theta_{10})\sqrt{I_j})\right),}
//' and
//' \deqn{q_{12} = P_{\theta_{10}}\left(\cup_{j=1}^{k}
//' \{(Z_j \geq b_j) \cup (Z_j \leq -w_{21} b_j +
//' (\theta_{20} - \theta_{10})\sqrt{I_j})\}\right).}
//' Here \eqn{Z_j = (\hat{\theta}_j - \theta_{10}) \sqrt{I_j}}, and
//' \eqn{w_{21} = \sqrt{\tilde{v}_{20}/\tilde{v}_1}}.
//' Of note, \eqn{q_1}, \eqn{q_2}, and \eqn{q_{12}}
//' can be evaluated using group sequential exit probabilities.
//' Similarly, we can define \eqn{\tilde{v}_2}, \eqn{\tilde{v}_{10}},
//' and \eqn{w_{12} = \sqrt{\tilde{v}_{10}/\tilde{v}_2}}, and
//' evaluate the type I error under \eqn{H_{20}}.
//'
//' The variance ratios correspond to
//' \deqn{\text{varianceRatioH10} = v_{10}/v_1,}
//' \deqn{\text{varianceRatioH20} = v_{20}/v_1,}
//' \deqn{\text{varianceRatioH12} = \tilde{v}_{10}/\tilde{v}_2,}
//' \deqn{\text{varianceRatioH21} = \tilde{v}_{20}/\tilde{v}_1.}
//' If the alternative variance is used, then the variance ratios
//' are all equal to 1.
//'
//' @return An S3 class \code{designEquiv} object with three components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{attainedAlphaH10}: The attained significance level under H10.
//'
//'     - \code{attainedAlphaH20}: The attained significance level under H20.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{thetaLower}: The parameter value at the lower equivalence
//'       limit.
//'
//'     - \code{thetaUpper}: The parameter value at the upper equivalence
//'       limit.
//'
//'     - \code{theta}: The parameter value under the alternative hypothesis.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedInformationH1}: The expected information under H1.
//'
//'     - \code{expectedInformationH10}: The expected information under H10.
//'
//'     - \code{expectedInformationH20}: The expected information under H20.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
//'       each of the two one-sided tests.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
//'       the two one-sided tests.
//'
//'     - \code{cumulativeAttainedAlphaH10}: The cumulative probability for
//'       efficacy stopping under H10.
//'
//'     - \code{cumulativeAttainedAlphaH20}: The cumulative probability for
//'       efficacy stopping under H20.
//'
//'     - \code{efficacyThetaLower}: The efficacy boundaries on the
//'       parameter scale for the one-sided null hypothesis at the
//'       lower equivalence limit.
//'
//'     - \code{efficacyThetaUpper}: The efficacy boundaries on the
//'       parameter scale for the one-sided null hypothesis at the
//'       upper equivalence limit.
//'
//'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
//'       each of the two one-sided tests.
//'
//'     - \code{information}: The cumulative information.
//'
//' * \code{settings}: A list containing the following components:
//'
//'     - \code{typeAlphaSpending}: The type of alpha spending.
//'
//'     - \code{parameterAlphaSpending}: The parameter value for alpha
//'       spending.
//'
//'     - \code{userAlphaSpending}: The user defined alpha spending.
//'
//'     - \code{spendingTime}: The error spending time at each analysis.
//'
//'     - \code{varianceRatioH10}: The ratio of the variance under H10 to
//'       the variance under H1.
//'
//'     - \code{varianceRatioH20}: The ratio of the variance under H20 to
//'       the variance under H1.
//'
//'     - \code{varianceRatioH12}: The ratio of the variance under H10 to
//'       the variance under H20.
//'
//'     - \code{varianceRatioH21}: The ratio of the variance under H20 to
//'       the variance under H10.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' # Example 1: obtain the maximum information given power
//' (design1 <- getDesignEquiv(
//'   beta = 0.2, thetaLower = log(0.8), thetaUpper = log(1.25),
//'   kMax = 2, informationRates = c(0.5, 1),
//'   alpha = 0.05, typeAlphaSpending = "sfOF"))
//'
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesignEquiv(
//'   IMax = 72.5, thetaLower = log(0.7), thetaUpper = -log(0.7),
//'   kMax = 3, informationRates = c(0.5, 0.75, 1),
//'   alpha = 0.05, typeAlphaSpending = "sfOF"))
//'
//' @export
// [[Rcpp::export]]
List getDesignEquiv(const double beta = NA_REAL,
                    const double IMax = NA_REAL,
                    const double thetaLower = NA_REAL,
                    const double thetaUpper = NA_REAL,
                    const double theta = 0,
                    const int kMax = 1,
                    const NumericVector& informationRates = NA_REAL,
                    const NumericVector& criticalValues = NA_REAL,
                    const double alpha = 0.05,
                    const std::string typeAlphaSpending = "sfOF",
                    const double parameterAlphaSpending = NA_REAL,
                    const NumericVector& userAlphaSpending = NA_REAL,
                    const NumericVector& spendingTime = NA_REAL,
                    const double varianceRatioH10 = 1,
                    const double varianceRatioH20 = 1,
                    const double varianceRatioH12 = 1,
                    const double varianceRatioH21 = 1) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector spendingTime1 = clone(spendingTime);

  double IMax1 = IMax;

  std::string unknown;

  if (R_isnancpp(beta) && R_isnancpp(IMax)) {
    stop("beta and IMax cannot be both missing");
  }

  if (!R_isnancpp(beta) && !R_isnancpp(IMax)) {
    stop("Only one of beta and IMax should be provided");
  }

  if (!R_isnancpp(IMax)) {
    if (IMax <= 0) {
      stop("IMax must be positive");
    }
    unknown = "beta";
  } else if (!R_isnancpp(beta)) {
    unknown = "IMax";
  }

  if (R_isnancpp(thetaLower)) {
    stop("thetaLower must be provided");
  }

  if (R_isnancpp(thetaUpper)) {
    stop("thetaUpper must be provided");
  }

  if (thetaLower >= theta) {
    stop("thetaLower must be less than theta");
  }

  if (thetaUpper <= theta) {
    stop("thetaUpper must be greater than theta");
  }

  if (kMax == NA_INTEGER) {
    stop("kMax must be provided");
  }

  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = NumericVector(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if ((unknown == "IMax") && (beta >= 1-alpha || beta < 0.0001)) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }


  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    spendingTime1 = clone(informationRates1);
  }

  if (varianceRatioH10 <= 0) {
    stop("varianceRatioH10 must be positive");
  }

  if (varianceRatioH20 <= 0) {
    stop("varianceRatioH20 must be positive");
  }

  if (varianceRatioH12 <= 0) {
    stop("varianceRatioH12 must be positive");
  }

  if (varianceRatioH21 <= 0) {
    stop("varianceRatioH21 must be positive");
  }


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, informationRates1,
                criticalValues, alpha](double aval)->double {
                  NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
                  for (int i=0; i<kMax-1; i++) {
                    u[i] = criticalValues[i];
                  }
                  u[kMax-1] = aval;

                  List probs = exitprobcpp(u, l, zero, informationRates1);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };

      criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      LogicalVector efficacyStopping1(kMax, 1);
      criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
                                    asf, asfpar, userAlphaSpending,
                                    spendingTime1, efficacyStopping1);
    }
  }

  NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);
  List probs = exitprobcpp(criticalValues1, li, zero, informationRates1);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));

  NumericVector efficacyP(kMax);
  for (int i=0; i<kMax; i++) {
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
  }

  double wH10 = sqrt(varianceRatioH10), wH20 = sqrt(varianceRatioH20);
  double wH12 = sqrt(varianceRatioH12), wH21 = sqrt(varianceRatioH21);

  // calculate cumulative rejection probability under H1
  NumericVector t = informationRates1;
  NumericVector b = criticalValues1;
  double deltaLower = thetaLower - theta;
  double deltaUpper = thetaUpper - theta;

  // obtain IMax if needed
  if (unknown == "IMax") {
    auto f = [beta, t, b, wH10, wH20, deltaLower, deltaUpper,
              li, ui, zero](double aval)->double {
                NumericVector I = t*aval;
                NumericVector l = b*wH10 + deltaLower*sqrt(I);
                NumericVector u = -b*wH20 + deltaUpper*sqrt(I);

                List probs1 = exitprobcpp(pmax(l, li), li, zero, I);
                List probs2 = exitprobcpp(ui, pmin(u, ui), zero, I);

                double cpl = sum(NumericVector(probs1[0]));
                double cpu = sum(NumericVector(probs2[1]));

                double power;
                if (is_true(any(l <= u))) {
                  power = cpl + cpu - 1;
                } else {
                  List a = exitprobcpp(l, u, zero, I);
                  double ca = sum(NumericVector(a[0]) + NumericVector(a[1]));
                  power = cpl + cpu - ca;
                }

                return power - (1-beta);
              };

    double z0 = R::qnorm(1-alpha, 0, 1, 1, 0);
    double z1 = R::qnorm(1-beta, 0, 1, 1, 0);
    double IMax10 = pow((z0*wH10 + z1)/deltaLower, 2);
    double IMax20 = pow((z0*wH20 + z1)/deltaUpper, 2);
    double IMaxLower = 0.5*std::min(IMax10, IMax20);
    double IMaxUpper = 1.5*std::max(IMax10, IMax20);
    IMax1 = brent(f, IMaxLower, IMaxUpper, 1.0e-6);
  }

  // obtain cumulative rejection probabilities under H1
  NumericVector I = t*IMax1;
  NumericVector l = b*wH10 + deltaLower*sqrt(I);
  NumericVector u = -b*wH20 + deltaUpper*sqrt(I);

  List probs1 = exitprobcpp(pmax(l, li), li, zero, I);
  List probs2 = exitprobcpp(ui, pmin(u, ui), zero, I);

  NumericVector cpl = cumsum(NumericVector(probs1[0]));
  NumericVector cpu = cumsum(NumericVector(probs2[1]));

  IntegerVector k = which(l >= u);
  NumericVector cp(kMax);
  if (k.size() == 0) {
    cp = cpl + cpu - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    List a = exitprobcpp(l[idx], u[idx], zero[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cp[i] = cpl[i] + cpu[i] - ca[i];
      } else {
        cp[i] = cpl[i] + cpu[i] - 1;
      }
    }
  }

  // incremental exit probabilities under H1
  NumericVector q(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      q[i] = cp[i];
    } else if (i<kMax-1) {
      q[i] = cp[i] - cp[i-1];
    } else {
      q[i] = 1 - cp[i-1];
    }
  }

  NumericVector rejectPerStage(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      rejectPerStage[i] = cp[i];
    } else {
      rejectPerStage[i] = cp[i] - cp[i-1];
    }
  }

  double overallReject = cp[kMax-1];
  double expectedInformationH1 = sum(q*I);

  NumericVector efficacyThetaLower = b/sqrt(I)*wH10 + thetaLower;
  NumericVector efficacyThetaUpper = -b/sqrt(I)*wH20 + thetaUpper;


  // cumulative rejection probability under H10
  NumericVector lH10 = b;
  NumericVector uH10 = -b*wH21 + (thetaUpper - thetaLower)*sqrt(I);
  List probs2H10 = exitprobcpp(ui, pmin(uH10, ui), zero, I);
  NumericVector cpuH10 = cumsum(NumericVector(probs2H10[1]));
  NumericVector cplH10 = cumAlphaSpent;

  IntegerVector kH10 = which(lH10 >= uH10);
  NumericVector cpH10(kMax);
  if (kH10.size() == 0) {
    cpH10 = cplH10 + cpuH10 - 1;
  } else {
    int K = max(kH10);
    IntegerVector idx = Range(0, K);
    List aH10 = exitprobcpp(lH10[idx], uH10[idx], zero[idx], I[idx]);
    NumericVector caH10 = cumsum(NumericVector(aH10[0]) +
      NumericVector(aH10[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH10[i] = cplH10[i] + cpuH10[i] - caH10[i];
      } else {
        cpH10[i] = cplH10[i] + cpuH10[i] - 1;
      }
    }
  }

  // incremental exit probabilities under H10
  NumericVector qH10(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      qH10[i] = cpH10[i];
    } else if (i<kMax-1) {
      qH10[i] = cpH10[i] - cpH10[i-1];
    } else {
      qH10[i] = 1 - cpH10[i-1];
    }
  }

  double attainedAlphaH10 = cpH10[kMax-1];
  double expectedInformationH10 = sum(qH10*I);


  // cumulative rejection probability under H20
  NumericVector lH20 = b*wH12 + (thetaLower - thetaUpper)*sqrt(I);
  NumericVector uH20 = -b;
  List probs1H20 = exitprobcpp(pmax(lH20, li), li, zero, I);
  NumericVector cplH20 = cumsum(NumericVector(probs1H20[0]));
  NumericVector cpuH20 = cumAlphaSpent;

  IntegerVector kH20 = which(lH20 >= uH20);
  NumericVector cpH20(kMax);
  if (kH20.size() == 0) {
    cpH20 = cplH20 + cpuH20 - 1;
  } else {
    int K = max(kH20);
    IntegerVector idx = Range(0, K);
    List aH20 = exitprobcpp(lH20[idx], uH20[idx], zero[idx], I[idx]);
    NumericVector caH20 = cumsum(NumericVector(aH20[0]) +
      NumericVector(aH20[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH20[i] = cplH20[i] + cpuH20[i] - caH20[i];
      } else {
        cpH20[i] = cplH20[i] + cpuH20[i] - 1;
      }
    }
  }

  // incremental exit probabilities under H20
  NumericVector qH20(kMax);
  for (int i=0; i<kMax; i++) {
    if (i==0) {
      qH20[i] = cpH20[i];
    } else if (i<kMax-1) {
      qH20[i] = cpH20[i] - cpH20[i-1];
    } else {
      qH20[i] = 1 - cpH20[i-1];
    }
  }

  double attainedAlphaH20 = cpH20[kMax-1];
  double expectedInformationH20 = sum(qH20*I);


  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["rejectPerStage"] = rejectPerStage,
    _["cumulativeRejection"] = cp,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["cumulativeAttainedAlphaH10"] = cpH10,
    _["cumulativeAttainedAlphaH20"] = cpH20,
    _["efficacyThetaLower"] = efficacyThetaLower,
    _["efficacyThetaUpper"] = efficacyThetaUpper,
    _["efficacyP"] = efficacyP,
    _["information"] = I);

  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = alpha,
    _["attainedAlphaH10"] = attainedAlphaH10,
    _["attainedAlphaH20"] = attainedAlphaH20,
    _["kMax"] = kMax,
    _["thetaLower"] = thetaLower,
    _["thetaUpper"] = thetaUpper,
    _["theta"] = theta,
    _["information"] = IMax1,
    _["expectedInformationH1"] = expectedInformationH1,
    _["expectedInformationH10"] = expectedInformationH10,
    _["expectedInformationH20"] = expectedInformationH20);

  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["spendingTime"] = spendingTime,
    _["varianceRatioH10"] = varianceRatioH10,
    _["varianceRatioH20"] = varianceRatioH20,
    _["varianceRatioH12"] = varianceRatioH12,
    _["varianceRatioH21"] = varianceRatioH21);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "designEquiv";

  return result;
}



//' @title Adaptive design at an interim look
//' @description Obtains the conditional power for specified incremental
//' information given the interim results, parameter value, and
//' data-dependent changes in the error spending function, and the number
//' and spacing of interim looks. Conversely, obtains the incremental
//' information needed to attain a specified conditional power given
//' the interim results, parameter value, and data-dependent changes
//' in the error spending function, and the number and spacing of
//' interim looks.
//'
//' @param betaNew The type II error for the secondary trial.
//' @param INew The maximum information of the secondary trial. Either
//'   \code{betaNew} or \code{INew} should be provided while the other one
//'   should be missing.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of
//'   the primary trial.
//' @param theta The parameter value.
//' @param IMax The maximum information of the primary trial. Must be
//'   provided if \code{futilityBounds} is missing and
//'   \code{typeBetaSpending} is not equal to "none", or
//'   if conditional power calculation is desired.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to true
//'   if left unspecified.
//' @param futilityStopping Indicators of whether futility stopping is
//'   allowed at each stage of the primary trial. Defaults to true
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping for the primary trial.
//' @param alpha The significance level of the primary trial.
//'   Defaults to 0.025.
//' @param typeAlphaSpending The type of alpha spending for the primary
//'   trial. One of the following:
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function,
//'   "user" for user defined spending, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to Delta for "WT", rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param userAlphaSpending The user defined alpha spending for the primary
//'   trial. Cumulative alpha spent up to each stage.
//' @param futilityBounds The lower boundaries on the z-test statistic scale
//'   for futility stopping for the primary trial. Defaults to
//'   \code{rep(-6, kMax-1)} if left unspecified.
//' @param typeBetaSpending The type of beta spending for the primary trial.
//'   One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early futility stopping.
//'   Defaults to "none".
//' @param parameterBetaSpending The parameter value of beta spending
//'   for the primary trial. Corresponds to rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to true
//'   if left unspecified.
//' @param futilityStoppingNew The indicators of whether futility stopping is
//'   allowed at each look of the secondary trial. Defaults to true
//'   if left unspecified.
//' @param typeAlphaSpendingNew The type of alpha spending for the secondary
//'   trial. One of the following:
//'   "OF" for O'Brien-Fleming boundaries,
//'   "P" for Pocock boundaries,
//'   "WT" for Wang & Tsiatis boundaries,
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early efficacy stopping.
//'   Defaults to "sfOF".
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to Delta for "WT", rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param typeBetaSpendingNew The type of beta spending for the secondary
//'   trial. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function,
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function,
//'   "user" for user defined spending, and
//'   "none" for no early futility stopping.
//'   Defaults to "none".
//' @param parameterBetaSpendingNew The parameter value of beta spending
//'   for the secondary trial. Corresponds to rho for "sfKD",
//'   and gamma for "sfHSD".
//' @param userBetaSpendingNew The user defined cumulative beta spending.
//'   Cumulative beta spent up to each stage of the secondary trial.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case, it is the same as
//'   \code{informationRatesNew}.
//' @param varianceRatio The ratio of the variance under H0 to the
//'   variance under H1.
//'
//' @return An \code{adaptDesign} object with two list components:
//'
//' * \code{primaryTrial}: A list of selected information for the primary
//'   trial, including \code{L}, \code{zL}, \code{theta}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
//'   and \code{MullerSchafer}.
//'
//' * \code{secondaryTrial}: A \code{design} object for the secondary trial.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Lu Chi, H. M. James Hung, and Sue-Jane Wang.
//' Modification of sample size in group sequential clinical trials.
//' Biometrics 1999;55:853-857.
//'
//' Hans-Helge Muller and Helmut Schafer.
//' Adaptive group sequential designs for clinical trials:
//' Combining the advantages of adaptive and of
//' classical group sequential approaches.
//' Biometrics 2001;57:886-891.
//'
//' @seealso \code{\link{getDesign}}
//'
//' @examples
//'
//' # original group sequential design with 90% power to detect delta = 6
//' delta = 6
//' sigma = 17
//' n = 282
//' (des1 = getDesign(IMax = n/(4*sigma^2), theta = delta, kMax = 3,
//'                   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'                   parameterAlphaSpending = -4))
//'
//' # interim look results
//' L = 1
//' n1 = n/3
//' delta1 = 4.5
//' sigma1 = 20
//' zL = delta1/sqrt(4/n1*sigma1^2)
//'
//' t = des1$byStageResults$informationRates
//'
//' # conditional power with sample size increase
//' (des2 = adaptDesign(
//'   betaNew = NA, INew = 420/(4*sigma1^2),
//'   L, zL, theta = delta1,
//'   IMax = n/(4*sigma1^2), kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4))
//'
//' # Muller & Schafer (2001) method to design the secondary trial:
//' # 3-look gamma(-2) spending with 84% power at delta = 4.5 and sigma = 20
//' (des2 = adaptDesign(
//'   betaNew = 0.16, INew = NA,
//'   L, zL, theta = delta1,
//'   IMax = n/(4*sigma1^2), kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpendingNew = "sfHSD",
//'   parameterAlphaSpendingNew = -2))
//'
//' # incremental sample size for sigma = 20
//' (nNew = 4*sigma1^2*des2$secondaryTrial$overallResults$information)
//'
//' @export
// [[Rcpp::export]]
List adaptDesign(double betaNew = NA_REAL,
                 double INew = NA_REAL,
                 const int L = NA_INTEGER,
                 const double zL = NA_REAL,
                 const double theta = NA_REAL,
                 const double IMax = NA_REAL,
                 const int kMax = NA_INTEGER,
                 const NumericVector& informationRates = NA_REAL,
                 const LogicalVector& efficacyStopping = NA_LOGICAL,
                 const LogicalVector& futilityStopping = NA_LOGICAL,
                 const NumericVector& criticalValues = NA_REAL,
                 const double alpha = 0.025,
                 const std::string typeAlphaSpending = "sfOF",
                 const double parameterAlphaSpending = NA_REAL,
                 const NumericVector& userAlphaSpending = NA_REAL,
                 const NumericVector& futilityBounds = NA_REAL,
                 const std::string typeBetaSpending = "none",
                 const double parameterBetaSpending = NA_REAL,
                 const NumericVector& spendingTime = NA_REAL,
                 const bool MullerSchafer = 0,
                 const int kNew = NA_INTEGER,
                 const NumericVector& informationRatesNew = NA_REAL,
                 const LogicalVector& efficacyStoppingNew = NA_LOGICAL,
                 const LogicalVector& futilityStoppingNew = NA_LOGICAL,
                 const std::string typeAlphaSpendingNew = "sfOF",
                 const double parameterAlphaSpendingNew = NA_REAL,
                 const std::string typeBetaSpendingNew = "none",
                 const double parameterBetaSpendingNew = NA_REAL,
                 const NumericVector& userBetaSpendingNew = NA_REAL,
                 const NumericVector& spendingTimeNew = NA_REAL,
                 const double varianceRatio = 1) {

  NumericVector t = clone(informationRates);
  LogicalVector es = clone(efficacyStopping);
  LogicalVector fs = clone(futilityStopping);
  NumericVector b = clone(criticalValues);
  NumericVector a = clone(futilityBounds);
  NumericVector st = clone(spendingTime);

  NumericVector tNew = clone(informationRatesNew);
  LogicalVector esNew = clone(efficacyStoppingNew);
  LogicalVector fsNew = clone(futilityStoppingNew);
  NumericVector stNew = clone(spendingTimeNew);

  double alpha1 = alpha;

  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfpar = parameterAlphaSpending;

  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfpar = parameterBetaSpending;

  std::string asfNew = typeAlphaSpendingNew;
  std::for_each(asfNew.begin(), asfNew.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double asfparNew = parameterAlphaSpendingNew;

  std::string bsfNew = typeBetaSpendingNew;
  std::for_each(bsfNew.begin(), bsfNew.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  double bsfparNew = parameterBetaSpendingNew;

  if (R_isnancpp(betaNew) && R_isnancpp(INew)) {
    stop("betaNew and INew cannot be both missing");
  }

  if (!R_isnancpp(betaNew) && !R_isnancpp(INew)) {
    stop("Only one of betaNew and INew should be provided");
  }

  if (!R_isnancpp(betaNew) && betaNew < 0.0001 && betaNew >= 1) {
    stop("betaNew must be greater than or equal to 0.0001 and less than 1");
  }

  if (!R_isnancpp(INew) && INew <= 0) {
    stop("INew must be positive");
  }

  if (L == NA_INTEGER) {
    stop("L must be provided");
  }

  if (L < 1) {
    stop("L must be a positive integer");
  }

  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }

  if (R_isnancpp(theta)) {
    stop("theta must be provided");
  }

  if (kMax == NA_INTEGER) {
    stop("kMax must be provided");
  }

  if (kMax <= L) {
    stop("kMax must be greater than L");
  }

  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    t = NumericVector(tem)/(kMax+0.0);
  }

  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) |
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    es = rep(1, kMax);
  }

  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    } else if (is_false(all((futilityStopping == 1) |
      (futilityStopping == 0)))) {
      stop("Elements of futilityStopping must be 1 or 0");
    }
  } else {
    fs = rep(1, kMax);
  }

  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }

  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && R_isnancpp(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }

  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }

  if (is_true(any(is_na(criticalValues))) && asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }

  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  }

  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }

    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      stop("futilityBounds and criticalValues must meet at final analysis");
    }
  }

  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfof" || bsf=="sfp" ||
      bsf=="sfkd" || bsf=="sfhsd" || bsf=="none")) {
    stop("Invalid value for typeBetaSpending");
  }

  if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
    stop("Missing value for parameterBetaSpending");
  }

  if (bsf=="sfkd" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }

  if (is_false(any(is_na(spendingTime)))) {
    if (spendingTime.size() != kMax) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[kMax-1] != 1) {
      stop("spendingTime must end with 1");
    }
  } else {
    st = clone(t);
  }


  if (MullerSchafer) {
    if (kNew == NA_INTEGER) {
      stop("kNew must be provided");
    }

    if (is_false(any(is_na(informationRatesNew)))) {
      if (informationRatesNew.size() != kNew) {
        stop("Invalid length for informationRatesNew");
      } else if (informationRatesNew[0] <= 0) {
        stop("Elements of informationRatesNew must be positive");
      } else if (kNew > 1 && is_true(any(diff(informationRatesNew) <= 0))) {
        stop("Elements of informationRatesNew must be increasing");
      } else if (informationRatesNew[kNew-1] != 1) {
        stop("informationRatesNew must end with 1");
      }
    } else {
      IntegerVector tem = seq_len(kNew);
      tNew = NumericVector(tem)/(kNew+0.0);
    }

    if (is_false(any(is_na(efficacyStoppingNew)))) {
      if (efficacyStoppingNew.size() != kNew) {
        stop("Invalid length for efficacyStoppingNew");
      } else if (efficacyStoppingNew[kNew-1] != 1) {
        stop("efficacyStoppingNew must end with 1");
      } else if (is_false(all((efficacyStoppingNew == 1) |
        (efficacyStoppingNew == 0)))) {
        stop("Elements of efficacyStoppingNew must be 1 or 0");
      }
    } else {
      esNew = rep(1, kNew);
    }

    if (is_false(any(is_na(futilityStoppingNew)))) {
      if (futilityStoppingNew.size() != kNew) {
        stop("Invalid length for futilityStoppingNew");
      } else if (futilityStoppingNew[kNew-1] != 1) {
        stop("futilityStoppingNew must end with 1");
      } else if (is_false(all((futilityStoppingNew == 1) |
        (futilityStoppingNew == 0)))) {
        stop("Elements of futilityStoppingNew must be 1 or 0");
      }
    } else {
      fsNew = rep(1, kNew);
    }

    if (!(asfNew=="of" || asfNew=="p" || asfNew=="wt" ||
        asfNew=="sfof" || asfNew=="sfp" ||
        asfNew=="sfkd" || asfNew=="sfhsd" || asfNew=="none")) {
      stop("Invalid value for typeAlphaSpendingNew");
    }

    if ((asfNew=="wt" || asfNew=="sfkd" || asfNew=="sfhsd") &&
        R_isnancpp(asfparNew)) {
      stop("Missing value for parameterAlphaSpendingNew");
    }

    if (asfNew=="sfkd" && asfparNew <= 0) {
      stop ("parameterAlphaSpendingNew must be positive for sfKD");
    }

    if (R_isnancpp(INew) && !(bsfNew=="sfof" || bsfNew=="sfp" ||
        bsfNew=="sfkd" || bsfNew=="sfhsd" ||
        bsfNew=="user" || bsfNew=="none")) {
      stop("Invalid value for typeBetaSpendingNew");
    } else if (!(bsfNew=="sfof" || bsfNew=="sfp" || bsfNew=="sfkd" ||
      bsfNew=="sfhsd" || bsfNew=="none")) {
      stop("Invalid value for typeBetaSpendingNew");
    }

    if ((bsfNew=="sfkd" || bsfNew=="sfhsd") && R_isnancpp(bsfparNew)) {
      stop("Missing value for parameterBetaSpendingNew");
    }

    if (bsfNew=="sfkd" && bsfparNew <= 0) {
      stop ("parameterBetaSpendingNew must be positive for sfKD");
    }

    if (R_isnancpp(INew) && bsfNew=="user") {
      if (is_true(any(is_na(userBetaSpendingNew)))) {
        stop("userBetaSpendingNew must be specified");
      } else if (userBetaSpendingNew.size() < kNew) {
        stop("Insufficient length of userBetaSpendingNew");
      } else if (userBetaSpendingNew[0] < 0) {
        stop("Elements of userBetaSpendingNew must be nonnegative");
      } else if (kNew > 1 && is_true(any(diff(userBetaSpendingNew) < 0))) {
        stop("Elements of userBetaSpendingNew must be nondecreasing");
      } else if (userBetaSpendingNew[kNew] != betaNew) {
        stop("userBetaSpendingNew must end with specified betaNew");
      }
    }

    if (is_false(any(is_na(spendingTimeNew)))) {
      if (spendingTimeNew.size() != kNew) {
        stop("Invalid length for spendingTimeNew");
      } else if (spendingTimeNew[0] <= 0) {
        stop("Elements of spendingTimeNew must be positive");
      } else if (kNew > 1 && is_true(any(diff(spendingTimeNew) <= 0))) {
        stop("Elements of spendingTimeNew must be increasing");
      } else if (spendingTimeNew[kNew-1] != 1) {
        stop("spendingTimeNew must end with 1");
      }
    } else {
      stNew = clone(tNew);
    }
  }

  if (varianceRatio <= 0) {
    stop("varianceRatio must be positive");
  }

  NumericVector w = rep(sqrt(varianceRatio), kMax);


  // obtain critical values for the primary trial
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        R_isnancpp(criticalValues[kMax-1])) { // Haybittle & Peto

      auto f = [kMax, t, es, criticalValues, alpha](double aval)->double {
        NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
        for (int i=0; i<kMax-1; i++) {
          u[i] = criticalValues[i];
          if (!es[i]) u[i] = 6.0;
        }
        u[kMax-1] = aval;

        List probs = exitprobcpp(u, l, zero, t);
        double cpu = sum(NumericVector(probs[0]));
        return cpu - alpha;
      };

      b[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
    } else {
      b = getBoundcpp(kMax, t, alpha, asf, asfpar, userAlphaSpending,
                      st, es);
    }
  }

  NumericVector l(kMax, -6.0), zero(kMax);
  List probs = exitprobcpp(b, l, zero, t);
  alpha1 = sum(NumericVector(probs[0]));

  // obtain futility bounds for the primary trial
  if (kMax > 1) {
    if (is_true(any(is_na(futilityBounds))) && bsf=="none") {
      a = rep(-6.0, kMax);
      a[kMax-1] = b[kMax-1];
    } else if (is_false(any(is_na(futilityBounds))) && a.size() == kMax-1) {
      a.push_back(b[kMax-1]);
    }
  } else {
    if (is_true(any(is_na(futilityBounds)))) {
      a = b[kMax-1];
    }
  }

  if (is_true(any(is_na(a)))) {
    if (R_isnancpp(IMax)) {
      stop("IMax must be provided");
    }

    if (IMax <= 0) {
      stop("IMax must be positive");
    }

    NumericVector theta1(kMax, theta);
    List out = getPower(alpha1, kMax, b, theta1, IMax*t, bsf, bsfpar,
                        st, fs, w);
    a = out[1];
  }

  int k1 = kMax - L;
  double alphaNew, conditionalPower, predictivePower;

  NumericVector t1(k1), r1(k1), b1(k1), a1(k1, -6.0), theta0(k1);
  for (int l=0; l<k1; l++) {
    t1[l] = (t[l+L] - t[L-1])/(1 - t[L-1]);
    r1[l] = t[L-1]/t[l+L];
    b1[l] = (b[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
    if (!es[l+L]) b1[l] = 6.0;
  }

  // conditional type I error
  probs = exitprobcpp(b1, a1, theta0, t1);
  alphaNew = sum(NumericVector(probs[0]));

  // conditional power
  for (int l=0; l<k1; l++) {
    a1[l] = (a[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
    if (!fs[l+L]) a1[l] = -6.0;
  }

  if (!R_isnancpp(IMax)) {
    double sigma = 1/sqrt(IMax*t[L-1]);
    double mu = zL*sigma;
    NumericVector theta1(k1, mu);

    NumericVector I1(k1);
    for (int l=0; l<k1; l++) {
      I1[l] = IMax*(t[l+L] - t[L-1]);
    }

    probs = exitprobcpp(b1, a1, theta1, I1);
    conditionalPower = sum(NumericVector(probs[0]));

    // predictive power
    auto f = [k1, b1, a1, I1](double theta)->double {
      NumericVector theta1(k1, theta);
      List probs = exitprobcpp(b1, a1, theta1, I1);
      return sum(NumericVector(probs[0]));
    };

    double lower = mu - 6*sigma, upper = mu + 6*sigma;
    predictivePower = intnorm(f, mu, sigma, lower, upper);
  } else {
    conditionalPower = NA_REAL;
    predictivePower = NA_REAL;
  }

  List des1 = List::create(
    _["L"] = L,
    _["zL"] = zL,
    _["theta"] = theta,
    _["kMax"] = kMax,
    _["informationRates"] = t,
    _["efficacyBounds"] = b,
    _["futilityBounds"] = a,
    _["conditionalAlpha"] = alphaNew,
    _["conditionalPower"] = conditionalPower,
    _["predictivePower"] = predictivePower,
    _["MullerSchafer"] = MullerSchafer);


  List des2;

  if (!MullerSchafer) {
    IntegerVector idx = Range(L, kMax-1);
    LogicalVector esNew = es[idx];
    LogicalVector fsNew = fs[idx];

    des2 = getDesign(betaNew, INew, theta, k1, t1, esNew, fsNew,
                     b1, NA_REAL, typeAlphaSpendingNew,
                     parameterAlphaSpendingNew, 0,
                     a1, typeBetaSpendingNew, parameterBetaSpendingNew,
                     userBetaSpendingNew, stNew, varianceRatio);
  } else {
    if (!R_isnancpp(betaNew) && betaNew >= 1-alphaNew) {
      stop("betaNew must be less than 1 minus conditional type I error");
    }

    NumericVector b1New(kNew, NA_REAL), a1New(kNew, NA_REAL);

    des2 = getDesign(betaNew, INew, theta, kNew, tNew, esNew, fsNew,
                     b1New, alphaNew, typeAlphaSpendingNew,
                     parameterAlphaSpendingNew, 0,
                     a1New, typeBetaSpendingNew, parameterBetaSpendingNew,
                     userBetaSpendingNew, stNew, varianceRatio);
  }

  List result = List::create(
    _["primaryTrial"] = des1,
    _["secondaryTrial"] = des2);

  result.attr("class") = "adaptDesign";

  return result;
}


// [[Rcpp::export]]
bool hasVariable(DataFrame df, std::string varName) {
  StringVector names = df.names();
  for (int i = 0; i < names.size(); i++) {
    if (names[i] == varName) {
      return true;
    }
  }
  return false;
}


// [[Rcpp::export]]
NumericMatrix invsympd(const NumericMatrix& a) {
  int n = a.nrow();
  if (a.ncol() != n) {
    stop("a is not a symmetric matrix");
  }

  int i,j,k;
  double sum;
  NumericMatrix L(n,n);  // A = L*L^T
  for (i=0; i<n; i++) {
    for (j=i; j<n; j++) {
      sum = a(i,j);
      for (k=0; k<i; k++) sum -= L(i,k)*L(j,k);
      if (i == j) {
        if (sum <= 0.0) stop("a is not positve definite");
        L(i,i) = sqrt(sum);
      } else L(j,i) = sum/L(i,i);
    }
  }

  NumericMatrix U(n,n); // U = L^-1
  for (i=0; i<n; i++) {
    U(i,i) = 1.0/L(i,i);
    for (j=i-1; j>=0; j--) {
      sum = 0.0;
      for (k=j+1; k<=i; k++) sum -= U(i,k)*L(k,j);
      U(i,j) = sum/L(j,j);
    }
  }

  NumericMatrix b(n,n); // A^-1 = U^T*U
  for (i=0; i<n; i++) {
    for (j=i; j<n; j++) {
      for (k=j; k<n; k++) {
        b(i,j) += U(k,i)*U(k,j);
      }
      if (i<j) b(j,i) = b(i,j);
    }
  }

  return b;
}


// [[Rcpp::export]]
double quantilecpp(const NumericVector& x, const double p) {
  int n = static_cast<int>(x.size());
  NumericVector y = clone(x);
  y.sort();
  double u = n*p + 1 - p;
  int j = floor(u);
  double g = u - j;
  double result = (1-g)*y[j-1] + g*y[j];
  return result;
}
