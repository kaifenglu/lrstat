#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "utilities.h"

using namespace Rcpp;


void set_seed(int seed) {
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}


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


//' @title Find Interval Numbers of Indices
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

//' @title Brent's Method for Root-Finding
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
    if (std::isnan(sfpar)) {
      stop("Parameter sfpar is missing for sfKD");
    } else if (sfpar <= 0) {
      stop ("sfpar must be positive for sfKD");
    } else {
      aval = error*pow(t, sfpar);
    }
  } else if (asf == "sfhsd") {
    if (std::isnan(sfpar)) {
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

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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


//' @title Integration With Respect to a Normal Density
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

//' @title Brent's Method for Minimization
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


void f_bvnorm(double *x, int n, void *ex) {
  bvnparams *param = (bvnparams *) ex;
  double corr = param->corr;
  double s = sqrt(1 - corr*corr);
  for (int i=0; i<n; i++) {
    double a = (param->a2 - corr*x[i])/s;
    double b = (param->b2 - corr*x[i])/s;
    double t1 = R::dnorm(x[i],0,1,0);
    double t2 = R::pnorm(b,0,1,1,0) - R::pnorm(a,0,1,1,0);
    x[i] = t1*t2;
  }
}


// [[Rcpp::export]]
double pbvnormcpp(NumericVector lower, NumericVector upper, double corr) {
  double result;
  if (corr == 0) {
    double v1 = R::pnorm(upper[0],0,1,1,0) - R::pnorm(lower[0],0,1,1,0);
    double v2 = R::pnorm(upper[1],0,1,1,0) - R::pnorm(lower[1],0,1,1,0);
    result = v1*v2;
  } else {
    double tol = 1.0e-8;
    bvnparams param = {corr, lower[1], upper[1]};
    result = quad(f_bvnorm, &param, lower[0], upper[0], tol)[0];
  }
  return result;
}


// [[Rcpp::export]]
NumericVector hazard_pdcpp(const NumericVector& piecewiseSurvivalTime,
                           const NumericVector& hazard_pfs,
                           const NumericVector& hazard_os,
                           const double corr_pd_os) {
  int n = static_cast<int>(piecewiseSurvivalTime.size());
  int i;
  NumericVector u(n);
  for (i=0; i<n-1; i++) {
    u[i] = piecewiseSurvivalTime[i+1];
  }
  u[n-1] = piecewiseSurvivalTime[n-1] + log(2)/hazard_pfs[n-1];

  NumericVector hazard_pd(n);
  NumericVector t(1), v(0), hazard(0), haz_pfs(0), haz_os(0);
  auto f = [&t, &v, &hazard, &haz_pfs, &haz_os,
            corr_pd_os](double haz)->double {
              NumericVector haz_pd = clone(hazard);
              haz_pd.push_back(haz);
              NumericVector lower(2);
              double a = ptpwexpcpp(t, v, haz_pd, 0, 1, 0)[0];
              double b = ptpwexpcpp(t, v, haz_os, 0, 1, 0)[0];
              lower[0] = R::qnorm(a,0,1,1,0);
              lower[1] = R::qnorm(b,0,1,1,0);
              NumericVector upper(2, R_PosInf);
              double q = pbvnormcpp(lower, upper, corr_pd_os);
              return q - ptpwexpcpp(t, v, haz_pfs, 0, 0, 0)[0];
            };

  double tol = 1e-6;
  for (i=0; i<n; i++) {
    t[0] = u[i];
    v.push_back(piecewiseSurvivalTime[i]);
    haz_pfs.push_back(hazard_pfs[i]);
    haz_os.push_back(hazard_os[i]);
    hazard_pd[i] = brent(f, 0.5*(hazard_pfs[i]-hazard_os[i]),
                         hazard_pfs[i], tol);
    hazard.push_back(hazard_pd[i]);
  }

  return hazard_pd;
}


// [[Rcpp::export]]
NumericVector hazard_subcpp(const NumericVector& piecewiseSurvivalTime,
                            const NumericVector& hazard_itt,
                            const NumericVector& hazard_pos,
                            const double p_pos) {
  int n = static_cast<int>(piecewiseSurvivalTime.size());
  int i;
  NumericVector u(n);
  for (i=0; i<n-1; i++) {
    u[i] = piecewiseSurvivalTime[i+1];
  }
  u[n-1] = piecewiseSurvivalTime[n-1] + log(2)/hazard_itt[n-1];

  NumericVector hazard_neg(n);
  NumericVector t(1), v(0), hazard(0), haz_itt(0), haz_pos(0);
  auto f = [&t, &v, &hazard, &haz_itt, &haz_pos, p_pos](double haz)->double {
    NumericVector haz_neg = clone(hazard);
    haz_neg.push_back(haz);
    double a = ptpwexpcpp(t, v, haz_pos, 0, 1, 0)[0];
    double b = ptpwexpcpp(t, v, haz_neg, 0, 1, 0)[0];
    double q = p_pos*a + (1-p_pos)*b;
    return q - ptpwexpcpp(t, v, haz_itt, 0, 1, 0)[0];
  };

  double tol = 1e-6;
  for (i=0; i<n; i++) {
    t[0] = u[i];
    v.push_back(piecewiseSurvivalTime[i]);
    haz_itt.push_back(hazard_itt[i]);
    haz_pos.push_back(hazard_pos[i]);
    hazard_neg[i] = brent(f, hazard_itt[i]-p_pos*hazard_pos[i],
                          hazard_itt[i]/(1-p_pos), tol);
    hazard.push_back(hazard_neg[i]);
  }

  return hazard_neg;
}


// Wrapper function for vmmin
List bmini(NumericVector x0, optimfn fn, optimgr gr, void *ex, double eps) {
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


//' @title Number of Enrolled Subjects
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


//' @title Accrual Duration to Enroll Target Number of Subjects
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


//' @title Probability of Being at Risk
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
  int J = static_cast<int>(piecewiseSurvivalTime.size());

  // hazard for failure or dropout
  NumericVector lambdax(J), gammax(J);

  if (lambda.size() == 1) {
    lambdax = rep(lambda, J);
  } else if (lambda.size() == J) {
    lambdax = lambda;
  } else {
    stop("Invalid length for lambda");
  }

  if (gamma.size() == 1) {
    gammax = rep(gamma, J);
  } else if (gamma.size() == J) {
    gammax = gamma;
  } else {
    stop("Invalid length for gamma");
  }

  NumericVector lamgam = lambdax + gammax;

  NumericVector t = piecewiseSurvivalTime;

  // sum up cumulative hazard up to time
  NumericVector a(k);
  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        a[i] += lamgam[j]*(t[j+1] - t[j]);
      } else {
        a[i] += lamgam[j]*(time[i] - t[j]);
      }
    }
  }

  return exp(-a);
}


//' @title Probability of Having an Event
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

  int J = static_cast<int>(piecewiseSurvivalTime.size());

  // hazard for failure or dropout
  NumericVector lambdax(J), gammax(J);

  if (lambda.size() == 1) {
    lambdax = rep(lambda, J);
  } else if (lambda.size() == J) {
    lambdax = lambda;
  } else {
    stop("Invalid length for lambda");
  }

  if (gamma.size() == 1) {
    gammax = rep(gamma, J);
  } else if (gamma.size() == J) {
    gammax = gamma;
  } else {
    stop("Invalid length for gamma");
  }

  NumericVector lamgam = lambdax + gammax;

  // sum up cumulative hazard up to time
  NumericVector t = piecewiseSurvivalTime;
  NumericVector n = patrisk(t, t, lambda, gamma);
  NumericVector a(k);
  double p;

  for (i=0; i<k; i++) {
    for (j=0; j<m[i]; j++) {
      if (j<m[i]-1) {
        p = lambda[j]/lamgam[j]*(1 - exp(-lamgam[j]*(t[j+1] - t[j])));
      } else {
        p = lambda[j]/lamgam[j]*(1 - exp(-lamgam[j]*(time[i] - t[j])));
      }
      a[i] += n[j]*p;
    }
  }

  return a;
}


//' @title Integrated Event Probability Over an Interval With Constant Hazard
//' @description Obtains the integrated probability of having an event
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

  int J = static_cast<int>(piecewiseSurvivalTime.size());

  // hazard for failure or dropout
  NumericVector lambdax(J), gammax(J);

  if (lambda.size() == 1) {
    lambdax = rep(lambda, J);
  } else if (lambda.size() == J) {
    lambdax = lambda;
  } else {
    stop("Invalid length for lambda");
  }

  if (gamma.size() == 1) {
    gammax = rep(gamma, J);
  } else if (gamma.size() == J) {
    gammax = gamma;
  } else {
    stop("Invalid length for gamma");
  }

  NumericVector lamgam = lambdax + gammax;

  // integration of conditional probability of having an event over (t1,t2)
  // given survival at the start of interval j
  double q1 = (exp(-lamgam[j1]*(t1-t0[0])) -
               exp(-lamgam[j1]*(t2-t0[0])))/lamgam[j1];
  double q = lambda[j1]/lamgam[j1] * (t2-t1 - q1);

  // sum up the integration for the already failed and to-be-failed
  return d0[0]*(t2-t1) + n0[0]*q;
}


//' @title Integrated Event Probability Over an Interval
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


//' @title Number of Patients Enrolled During an Interval and Having an Event
//' by Specified Calendar Times
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


//' @title Number of Subjects at Risk
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


//' @title Number of Subjects Having an Event
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


//' @title Number of Subjects Having an Event by Calendar Time
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


//' @title Power and Sample Size for a Generic Group Sequential Design
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

  if (std::isnan(beta) && std::isnan(IMax)) {
    stop("beta and IMax cannot be both missing");
  }

  if (!std::isnan(beta) && !std::isnan(IMax)) {
    stop("Only one of beta and IMax should be provided");
  }

  if (!std::isnan(IMax)) {
    if (IMax <= 0) {
      stop("IMax must be positive");
    }
    unknown = "beta";
  } else if (!std::isnan(beta)) {
    unknown = "IMax";
  }

  if (std::isnan(theta)) {
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

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
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

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
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
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

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


//' @title Power and Sample Size for a Generic Group Sequential Equivalence
//' Design
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
//'
//' @details
//' Consider the equivalence design with two one-sided hypotheses:
//' \deqn{H_{10}: \theta \leq \theta_{10},}
//' \deqn{H_{20}: \theta \geq \theta_{20}.}
//' We reject \eqn{H_{10}} at or before look \eqn{k} if
//' \deqn{Z_{1j} = (\hat{\theta}_j - \theta_{10})\sqrt{I_j}
//' \geq b_j}
//' for some \eqn{j=1,\ldots,k}, where \eqn{\{b_j:j=1,\ldots,K\}} are the
//' critical values associated with the specified alpha-spending function,
//' and \eqn{I_j} is the information for \eqn{\theta} (inverse variance of
//' \eqn{\hat{\theta}}) at the
//' \eqn{j}th look. For example,
//' for estimating the risk difference \eqn{\theta = \pi_1 - \pi_2},
//' \deqn{I_j = \left\{\frac{\pi_1 (1-\pi_1)}{n_{1j}} +
//' \frac{\pi_2(1-\pi_2)}{n_{2j}}\right\}^{-1}.}
//' It follows that
//' \deqn{(Z_{1j} \geq b_j) = (Z_j \geq b_j +
//' \theta_{10}\sqrt{I_j}),}
//' where \eqn{Z_j = \hat{\theta}_j \sqrt{I_j}}.
//'
//' Similarly, we reject \eqn{H_{20}} at or before look \eqn{k} if
//' \deqn{Z_{2j} = (\hat{\theta}_j - \theta_{20})\sqrt{I_j}
//' \leq -b_j} for some \eqn{j=1,\ldots,k}. We have
//' \deqn{(Z_{2j} \leq -b_j) = (Z_j \leq - b_j +
//' \theta_{20}\sqrt{I_j}).}
//'
//' Let \eqn{l_j = b_j + \theta_{10}\sqrt{I_j}},
//' and \eqn{u_j = -b_j + \theta_{20}\sqrt{I_j}}.
//' The cumulative probability to reject \eqn{H_0 = H_{10} \cup H_{20}} at
//' or before look \eqn{k} under the alternative hypothesis \eqn{H_1} is
//' given by
//' \deqn{P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j) \cap
//' \cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right) = p_1 + p_2 - p_{12},}
//' where
//' \deqn{p_1 = P_\theta\left(\cup_{j=1}^{k} (Z_{1j} \geq b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \geq l_j)\right),}
//' \deqn{p_2 = P_\theta\left(\cup_{j=1}^{k} (Z_{2j} \leq -b_j)\right)
//' = P_\theta\left(\cup_{j=1}^{k} (Z_j \leq u_j)\right),}
//' and
//' \deqn{p_{12} = P_\theta\left(\cup_{j=1}^{k} (Z_j \geq l_j) \cup
//' (Z_j \leq u_j)\right).}
//' Of note, both \eqn{p_1} and \eqn{p_2} can be evaluated using
//' one-sided exit probabilities for group sequential designs.
//' If there exists \eqn{j\leq k} such that \eqn{l_j \leq u_j}, then
//' \eqn{p_{12} = 1}. Otherwise, \eqn{p_{12}} can be evaluated using
//' two-sided exit probabilities for group sequential designs.
//'
//' Since the equivalent hypothesis is tested using two one-sided tests,
//' the type I error is controlled. To evaluate the attained type I error
//' of the equivalence trial under \eqn{H_{10}} (or \eqn{H_{20}}),
//' we simply fix the control group parameters, update the active
//' treatment group parameters according to the null hypothesis, and
//' use the parameters in the power calculation outlined above.
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
                    const NumericVector& spendingTime = NA_REAL) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector spendingTime1 = clone(spendingTime);

  double IMax1 = IMax;

  std::string unknown;

  if (std::isnan(beta) && std::isnan(IMax)) {
    stop("beta and IMax cannot be both missing");
  }

  if (!std::isnan(beta) && !std::isnan(IMax)) {
    stop("Only one of beta and IMax should be provided");
  }

  if (!std::isnan(IMax)) {
    if (IMax <= 0) {
      stop("IMax must be positive");
    }
    unknown = "beta";
  } else if (!std::isnan(beta)) {
    unknown = "IMax";
  }

  if (std::isnan(thetaLower)) {
    stop("thetaLower must be provided");
  }

  if (std::isnan(thetaUpper)) {
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

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
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

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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


  // obtain criticalValues
  if (is_true(any(is_na(criticalValues)))) {
    if (kMax > 1 && criticalValues.size() == kMax &&
        is_false(any(is_na(head(criticalValues, kMax-1)))) &&
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

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

  // calculate cumulative rejection probability under H1
  NumericVector t = informationRates1;
  NumericVector b = criticalValues1;
  double deltaLower = thetaLower - theta;
  double deltaUpper = thetaUpper - theta;

  // obtain IMax if needed
  if (unknown == "IMax") {
    auto f = [beta, t, b, deltaLower, deltaUpper,
              li, ui, zero](double aval)->double {
                NumericVector I = t*aval;
                NumericVector l = b + deltaLower*sqrt(I);
                NumericVector u = -b + deltaUpper*sqrt(I);

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
    double IMax0 = pow((z0 + z1)/deltaLower, 2);
    double IMaxLower = 0.5*IMax0;
    double IMaxUpper = 1.5*IMax0;
    IMax1 = brent(f, IMaxLower, IMaxUpper, 1.0e-6);
  }

  // obtain cumulative rejection probabilities under H1
  NumericVector I = t*IMax1;
  NumericVector l = b + deltaLower*sqrt(I);
  NumericVector u = -b + deltaUpper*sqrt(I);

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

  NumericVector efficacyThetaLower = b/sqrt(I) + thetaLower;
  NumericVector efficacyThetaUpper = -b/sqrt(I) + thetaUpper;


  NumericVector theta10 = rep(deltaLower, kMax);
  NumericVector theta20 = rep(deltaUpper, kMax);

  // cumulative rejection probability under H10
  List probs2H10 = exitprobcpp(ui, pmin(u, ui), theta10, I);
  NumericVector cplH10 = cumAlphaSpent;
  NumericVector cpuH10 = cumsum(NumericVector(probs2H10[1]));

  NumericVector cpH10(kMax);
  if (k.size() == 0) {
    cpH10 = cplH10 + cpuH10 - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    List a = exitprobcpp(l[idx], u[idx], theta10[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH10[i] = cplH10[i] + cpuH10[i] - ca[i];
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
  List probs1H20 = exitprobcpp(pmax(l, li), li, theta20, I);
  NumericVector cplH20 = cumsum(NumericVector(probs1H20[0]));
  NumericVector cpuH20 = cumAlphaSpent;

  NumericVector cpH20(kMax);
  if (k.size() == 0) {
    cpH20 = cplH20 + cpuH20 - 1;
  } else {
    int K = max(k);
    IntegerVector idx = Range(0, K);
    List a = exitprobcpp(l[idx], u[idx], theta20[idx], I[idx]);
    NumericVector ca = cumsum(NumericVector(a[0]) +
      NumericVector(a[1]));

    for (int i=0; i<kMax; i++) {
      if (i <= K) {
        cpH20[i] = cplH20[i] + cpuH20[i] - ca[i];
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
    _["spendingTime"] = spendingTime);

  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);

  result.attr("class") = "designEquiv";

  return result;
}



//' @title Adaptive Design at an Interim Look
//' @description
//' Calculates the conditional power for specified incremental
//' information, given the interim results, parameter value,
//' data-dependent changes in the error spending function, and
//' the number and spacing of interim looks. Conversely,
//' calculates the incremental information required to attain
//' a specified conditional power, given the interim results,
//' parameter value, data-dependent changes in the error
//' spending function, and the number and spacing of interim looks.
//'
//' @param betaNew The type II error for the secondary trial.
//' @param INew The maximum information of the secondary trial. Either
//'   \code{betaNew} or \code{INew} should be provided, while the other
//'   must be missing.
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
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param futilityStopping Indicators of whether futility stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
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
//' @param userAlphaSpending The user-defined alpha spending for the
//'   primary trial. Represents the cumulative alpha spent up to each stage.
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
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param futilityStoppingNew The indicators of whether futility stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
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
//' @param userBetaSpendingNew The user-defined cumulative beta spending.
//'   Represents the cumulative beta spent up to each stage of the
//'   secondary trial.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
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
//'   L = L, zL = zL, theta = delta1,
//'   IMax = n/(4*sigma1^2), kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4))
//'
//' # Muller & Schafer (2001) method to design the secondary trial:
//' # 3-look gamma(-2) spending with 84% power at delta = 4.5 and sigma = 20
//' (des2 = adaptDesign(
//'   betaNew = 0.16, INew = NA,
//'   L = L, zL = zL, theta = delta1,
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

  if (std::isnan(betaNew) && std::isnan(INew)) {
    stop("betaNew and INew cannot be both missing");
  }

  if (!std::isnan(betaNew) && !std::isnan(INew)) {
    stop("Only one of betaNew and INew should be provided");
  }

  if (!std::isnan(betaNew) && betaNew < 0.0001 && betaNew >= 1) {
    stop("betaNew must be greater than or equal to 0.0001 and less than 1");
  }

  if (!std::isnan(INew) && INew <= 0) {
    stop("INew must be positive");
  }

  if (L == NA_INTEGER) {
    stop("L must be provided");
  }

  if (L < 1) {
    stop("L must be a positive integer");
  }

  if (std::isnan(zL)) {
    stop("zL must be provided");
  }

  if (std::isnan(theta)) {
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

  if (!std::isnan(alpha)) {
    if (alpha < 0.00001 || alpha >= 1) {
      stop("alpha must lie in [0.00001, 1)");
    }
  }

  if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
    stop("alpha must be provided when criticalValues is missing");
  }

  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }

  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
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

  if ((bsf=="sfkd" || bsf=="sfhsd") && std::isnan(bsfpar)) {
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
        std::isnan(asfparNew)) {
      stop("Missing value for parameterAlphaSpendingNew");
    }

    if (asfNew=="sfkd" && asfparNew <= 0) {
      stop ("parameterAlphaSpendingNew must be positive for sfKD");
    }

    if (std::isnan(INew) && !(bsfNew=="sfof" || bsfNew=="sfp" ||
        bsfNew=="sfkd" || bsfNew=="sfhsd" ||
        bsfNew=="user" || bsfNew=="none")) {
      stop("Invalid value for typeBetaSpendingNew");
    } else if (!(bsfNew=="sfof" || bsfNew=="sfp" || bsfNew=="sfkd" ||
      bsfNew=="sfhsd" || bsfNew=="none")) {
      stop("Invalid value for typeBetaSpendingNew");
    }

    if ((bsfNew=="sfkd" || bsfNew=="sfhsd") && std::isnan(bsfparNew)) {
      stop("Missing value for parameterBetaSpendingNew");
    }

    if (bsfNew=="sfkd" && bsfparNew <= 0) {
      stop ("parameterBetaSpendingNew must be positive for sfKD");
    }

    if (std::isnan(INew) && bsfNew=="user") {
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
        std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto

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
    if (std::isnan(IMax)) {
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

  if (!std::isnan(IMax)) {
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
    if (!std::isnan(betaNew) && betaNew >= 1-alphaNew) {
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


double quantilecpp(const NumericVector& x, const double p) {
  int n = static_cast<int>(x.size());
  NumericVector y = clone(x);
  y.sort();
  double u = n*p + 1 - p;
  int j = static_cast<int>(std::floor(u));
  double g = u - j;
  double result = (1-g)*y[j-1] + g*y[j];
  return result;
}


// [[Rcpp::plugins(cpp11)]]
double squantilecpp(const std::function<double(double)>& S, double p) {
  double lower = 0;
  double upper = 1;
  while (S(upper) > p) {
    lower = upper;
    upper = 2*upper;
  }

  auto f = [S, p](double t)->double{
    return S(t) - p;
  };

  return brent(f, lower, upper, 1e-6);
}


IntegerVector c_vectors_i(IntegerVector vec1, IntegerVector vec2) {
  IntegerVector result(vec1.size() + vec2.size());
  std::copy(vec1.begin(), vec1.end(), result.begin());
  std::copy(vec2.begin(), vec2.end(), result.begin() + vec1.size());
  return result;
}


NumericVector c_vectors(NumericVector vec1, NumericVector vec2) {
  NumericVector result(vec1.size() + vec2.size());
  std::copy(vec1.begin(), vec1.end(), result.begin());
  std::copy(vec2.begin(), vec2.end(), result.begin() + vec1.size());
  return result;
}


NumericMatrix subset_matrix_by_row(NumericMatrix a, IntegerVector q) {
  int i, j, n = static_cast<int>(q.size()), p = a.ncol();
  NumericMatrix b(n,p);
  for (j=0; j<p; j++) {
    for (i=0; i<n; i++) {
      b(i,j) = a(q[i],j);
    }
  }
  return b;
}


NumericMatrix c_matrices(NumericMatrix a1, NumericMatrix a2) {
  int h, i, j, n1 = a1.nrow(), n2 = a2.nrow(), p = a1.ncol();
  NumericMatrix b(n1+n2, p);
  for (i=0; i<n1; i++) {
    for (j=0; j<p; j++) {
      b(i,j) = a1(i,j);
    }
  }

  for (i=0; i<n2; i++) {
    h = i+n1;
    for (j=0; j<p; j++) {
      b(h,j) = a2(i,j);
    }
  }

  return b;
}


List bygroup(DataFrame data, const StringVector& variables) {
  int i;
  int n = data.nrows();
  int p = static_cast<int>(variables.size());

  IntegerVector d(p);   // the number of unique values
  List u(p);            // the vector of unique values
  IntegerMatrix x(n,p); // indices of original values in unique values
  for (i=0; i<p; i++) {
    String s = variables[i];
    if (!hasVariable(data, s)) {
      stop("data must contain the variables");
    }

    if (TYPEOF(data[s]) == LGLSXP || TYPEOF(data[s]) == INTSXP) {
      IntegerVector v = data[s];
      IntegerVector w = unique(v);
      w.sort();
      d[i] = static_cast<int>(w.size());
      u[i] = w;
      x(_,i) = match(v,w) - 1;
    } if (TYPEOF(data[s]) == REALSXP) {
      NumericVector v = data[s];
      NumericVector w = unique(v);
      w.sort();
      d[i] = static_cast<int>(w.size());
      u[i] = w;
      x(_,i) = match(v,w) - 1;
    } if (TYPEOF(data[s]) == STRSXP) {
      StringVector v = data[s];
      StringVector w = unique(v);
      w.sort();
      d[i] = static_cast<int>(w.size());
      u[i] = w;
      x(_,i) = match(v,w) - 1;
    }
  }

  int frac = 1;
  int orep = 1;
  for (i=0; i<p; i++) {
    orep = orep*d[i];
  }

  IntegerVector index(n);
  DataFrame lookup;
  for (i=0; i<p; i++) {
    orep = orep/d[i];
    index = index + x(_,i)*orep;

    IntegerVector j = rep(rep_each(seq(0, d[i]-1), orep), frac);
    String s = variables[i];
    if (TYPEOF(data[s]) == LGLSXP || TYPEOF(data[s]) == INTSXP) {
      IntegerVector w = u[i];
      lookup.push_back(w[j],s);
    } else if (TYPEOF(data[s]) == REALSXP) {
      NumericVector w = u[i];
      lookup.push_back(w[j],s);
    } else if (TYPEOF(data[s]) == STRSXP) {
      StringVector w = u[i];
      lookup.push_back(w[j],s);
    }

    frac = frac*d[i];
  }

  return List::create(
    Named("nlevels") = d,
    Named("indices") = x+1,
    Named("lookups") = u,
    Named("index") = index+1,
    Named("lookup") = lookup);
}


// The following three utilities functions are from the survival package
// and are used to compute the Cholesky decomposition of a symmetric
// positive-definite matrix, solve a linear system, and compute the
// inverse of a symmetric positive-definite matrix.

// The matrix A is modified in place. Let A = U' d U, where U is upper
// triangular with unit diagonals, and d is a diagonal matrix.
// The lower triangular part of A is left unchanged, the diagonal part
// is modified to contain d, and the upper triangular part is modified
// to contain U. The toler parameter
// is used to determine the threshold for considering a diagonal
// element as zero. If the diagonal element is less than toler times
// the largest diagonal element, it is considered zero. The function
// returns the rank of the matrix, which is the number of non-zero
// diagonal elements in the Cholesky decomposition.
// [[Rcpp::export]]
int cholesky2(NumericMatrix matrix, int n, double toler) {
  double temp;
  int i, j, k;
  double eps, pivot;
  int rank;
  int nonneg;

  nonneg = 1;
  eps = 0;
  for (i=0; i<n; i++) {
    if (matrix(i,i) > eps) eps = matrix(i,i);
  }
  if (eps==0) eps = toler; // no positive diagonals!
  else eps *= toler;

  rank = 0;
  for (i=0; i<n; i++) {
    pivot = matrix(i,i);
    if (std::isinf(pivot) == 1 || pivot < eps) {
      matrix(i,i) = 0;
      if (pivot < -8*eps) nonneg = -1;
    }
    else  {
      rank++;
      for (j=i+1; j<n; j++) {
        temp = matrix(i,j)/pivot;
        matrix(i,j) = temp;
        matrix(j,j) -= temp*temp*pivot;
        for (k=j+1; k<n; k++) matrix(j,k) -= temp*matrix(i,k);
      }
    }
  }

  return(rank*nonneg);
}

// [[Rcpp::export]]
void chsolve2(NumericMatrix matrix, int n, NumericVector y) {
  int i, j;
  double temp;

  for (i=0; i<n; i++) {
    temp = y[i];
    for (j=0; j<i; j++)
      temp -= y[j]*matrix(j,i);
    y[i] = temp;
  }

  for (i=n-1; i>=0; i--) {
    if (matrix(i,i) == 0) y[i] = 0;
    else {
      temp = y[i]/matrix(i,i);
      for (j=i+1; j<n; j++)
        temp -= y[j]*matrix(i,j);
      y[i] = temp;
    }
  }
}


void chinv2(NumericMatrix matrix, int n) {
  double temp;
  int i, j, k;

  for (i=0; i<n; i++){
    if (matrix(i,i) > 0) {
      matrix(i,i) = 1/matrix(i,i);   // this line inverts D
      for (j=i+1; j<n; j++) {
        matrix(i,j) = -matrix(i,j);
        for (k=0; k<i; k++)     // sweep operator
          matrix(k,j) += matrix(i,j)*matrix(k,i);
      }
    }
  }

  for (i=0; i<n; i++) {
    if (matrix(i,i) == 0) {  // singular row
      for (j=0; j<i; j++) matrix(i,j) = 0;
      for (j=i; j<n; j++) matrix(j,i) = 0;
    }
    else {
      for (j=i+1; j<n; j++) {
        temp = matrix(i,j)*matrix(j,j);
        matrix(j,i) = temp;
        for (k=i; k<j; k++)
          matrix(k,i) += temp*matrix(k,j);
      }
    }
  }
}


NumericMatrix invsympd(NumericMatrix matrix, int n, double toler) {
  int i, j;
  NumericMatrix v = clone(matrix);
  i = cholesky2(v, n, toler);
  chinv2(v, n);
  for (i=1; i<n; i++) {
    for (j=0; j<i; j++) {
      v(j,i) = v(i,j);
    }
  }

  return v;
}


//' @title Split a survival data set at specified cut points
//' @description For a given survival dataset and specified cut times,
//' each record is split into multiple subrecords at each cut time.
//' The resulting dataset is in counting process format, with each
//' subrecord containing a start time, stop time, and event status.
//' This is adapted from the survsplit.c function from the survival package.
//'
//' @param tstart The starting time of the time interval for
//'   counting-process data.
//' @param tstop The stopping time of the time interval for
//'   counting-process data.
//' @param cut The vector of cut points.
//'
//' @return A data frame with the following variables:
//'
//' * \code{row}: The row number of the observation in the input data
//'   (starting from 0).
//'
//' * \code{start}: The starting time of the resulting subrecord.
//'
//' * \code{end}: The ending time of the resulting subrecord.
//'
//' * \code{censor}: Whether the subrecord lies strictly within a record
//'   in the input data (1 for all but the last interval and 0 for the
//'   last interval with cutpoint set equal to tstop).
//'
//' * \code{interval}: The interval number derived from cut (starting
//'   from 0 if the interval lies to the left of the first cutpoint).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @keywords internal
//'
//' @examples
//'
//' survsplit(15, 60, c(10, 30, 40))
//'
//' @export
// [[Rcpp::export]]
DataFrame survsplit(NumericVector tstart,
                    NumericVector tstop,
                    NumericVector cut) {
  int i, j, k, extra;
  int n = static_cast<int>(tstart.size());
  int ncut = static_cast<int>(cut.size());

  // Each cut point strictly within an interval generates an extra line.
  // NA inputs are left alone.
  extra = 0;
  for (i=0; i<n; i++) {
    for (j=0; j<ncut; j++) {
      if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
          cut[j] > tstart[i] && cut[j] < tstop[i]) extra++;
    }
  }

  int n2 = n + extra;
  IntegerVector row(n2), interval(n2);
  NumericVector start(n2), end(n2);
  LogicalVector censor(n2);

  k = 0;
  for (i=0; i<n; i++) {
    if (std::isnan(tstart[i]) || std::isnan(tstop[i])) {
      start[k] = tstart[i];
      end[k] = tstop[i];
      row[k] = i;           // row in the original data
      interval[k] = 1;
      k++;
    } else {
      // find the first cut point after tstart
      for (j=0; j < ncut && cut[j] <= tstart[i]; j++);
      start[k] = tstart[i];
      row[k] = i;
      interval[k] = j;
      for (; j < ncut && cut[j] < tstop[i]; j++) {
        if (cut[j] > tstart[i]) {
          end[k] = cut[j];
          censor[k] = 1;
          k++; // create the next sub-interval
          start[k] = cut[j];
          row[k] = i;
          interval[k] = j+1;
        }
      }
      end[k] = tstop[i]; // finish the last sub-interval
      censor[k] = 0;
      k++;
    }
  }

  DataFrame result = DataFrame::create(
    Named("row") = row,
    Named("start") = start,
    Named("end") = end,
    Named("censor") = censor,
    Named("interval") = interval);

  return result;
}



bool is_sorted(NumericVector x) {
  int n = x.size();

  // Loop through the vector and check if it is sorted
  for (int i = 1; i < n; ++i) {
    if (x[i] < x[i - 1]) {
      return 0;  // Return false if any element is smaller than the previous
    }
  }

  return 1;  // If no violations, the vector is sorted
}


// Householder vector
// Given an n-vector x, this function computes an n-vector v with v(1) = 1
// such that (I - 2*v*t(v)/t(v)*v)*x is zero in all but the first component.
NumericVector house(const NumericVector& x) {
  int n = static_cast<int>(x.size());
  double mu = sqrt(sum(x*x));
  NumericVector v = clone(x);
  if (mu > 0.0) {
    double beta = x[0] + std::copysign(1.0, x[0])*mu;
    for (int i=1; i<n; i++) {
      v[i] /= beta;
    }
  }
  v[0] = 1.0;
  return v;
}


// Householder pre-multiplication
// Given an m-by-n matrix A and a nonzero m-vector v with v(1) = 1,
// the following algorithm overwrites A with P*A where
// P = I - 2*v*t(v)/t(v)*v.
void row_house(NumericMatrix& A, const int i1, const int i2,
               const int j1, const int j2, const NumericVector& v) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow()) {
    stop("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol()) {
    stop("Invalid column indices j1 and j2");
  }

  int i, j, m = i2-i1+1, n = j2-j1+1;
  double beta = -2.0/sum(v*v);
  NumericVector w(n);
  for (j=0; j<n; j++) {
    for (i=0; i<m; i++) {
      w[j] += A(i+i1,j+j1)*v[i];
    }
    w[j] *= beta;
  }

  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      A(i+i1,j+j1) += v[i]*w[j];
    }
  }
}


// Householder post-multiplication
// Given an m-by-n matrix A and a nonzero n-vector v with v(1) = 1,
// the following algorithm overwrites A with A*P where
// P = I - 2*v*t(v)/t(v)*v.
void col_house(NumericMatrix& A, const int i1, const int i2,
               const int j1, const int j2, const NumericVector& v) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow()) {
    stop("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol()) {
    stop("Invalid column indices j1 and j2");
  }

  int i, j, m = i2-i1+1, n = j2-j1+1;
  double beta = -2.0/sum(v*v);
  NumericVector w(m);
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      w[i] += A(i+i1,j+j1)*v[j];
    }
    w[i] *= beta;
  }

  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      A(i+i1,j+j1) += w[i]*v[j];
    }
  }
}


//' @title QR Decomposition of a Matrix
//' @description Computes the QR decomposition of a matrix.
//'
//' @param X A numeric matrix whose QR decomposition is to be computed.
//' @param tol The tolerance for detecting linear dependencies in the
//'   columns of \code{X}.
//'
//' @details
//' This function performs Householder QR with column pivoting:
//' Given an \eqn{m}-by-\eqn{n} matrix \eqn{A} with \eqn{m \geq n},
//' the following algorithm computes \eqn{r = \textrm{rank}(A)} and
//' the factorization \eqn{Q^T A P} equal to
//' \tabular{ccccc}{
//' | \tab \eqn{R_{11}} \tab \eqn{R_{12}} \tab | \tab \eqn{r} \cr
//' | \tab 0 \tab 0 \tab | \tab \eqn{m-r} \cr
//'   \tab \eqn{r} \tab \eqn{n-r} \tab \tab
//' }
//' with \eqn{Q = H_1 \cdots H_r} and \eqn{P = P_1 \cdots P_r}.
//' The upper triangular part of \eqn{A}
//' is overwritten by the upper triangular part of \eqn{R} and
//' components \eqn{(j+1):m} of
//' the \eqn{j}th Householder vector are stored in \eqn{A((j+1):m, j)}.
//' The permutation \eqn{P} is encoded in an integer vector \code{pivot}.
//'
//' @return A list with the following components:
//'
//' * \code{qr}: A matrix with the same dimensions as \code{X}. The upper
//'   triangle contains the \code{R} of the decomposition and the lower
//'   triangle contains Householder vectors (stored in compact form).
//'
//' * \code{rank}: The rank of \code{X} as computed by the decomposition.
//'
//' * \code{pivot}: The column permutation for the pivoting strategy used
//'   during the decomposition.
//'
//' * \code{Q}: The complete \eqn{m}-by-\eqn{m} orthogonal matrix \eqn{Q}.
//'
//' * \code{R}: The complete \eqn{m}-by-\eqn{n} upper triangular
//'   matrix \eqn{R}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Gene N. Golub and Charles F. Van Loan.
//' Matrix Computations, second edition. Baltimore, Maryland:
//' The John Hopkins University Press, 1989, p.235.
//'
//' @examples
//'
//' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
//' h9 <- hilbert(9)
//' qrcpp(h9)
//'
//' @export
// [[Rcpp::export]]
List qrcpp(const NumericMatrix& X, double tol = 1e-12) {
  int i, j, k, l, m = X.nrow(), n = X.ncol();
  NumericMatrix A = clone(X);
  NumericVector c(n);
  for (j=0; j<n; j++) {
    c[j] = sum(A(_,j)*A(_,j));
  }

  double tau = max(c);
  for (k=0; k<n; k++) {
    if (c[k] > tol) break;
  }

  int r = -1;
  IntegerVector piv = seq(0,n-1);
  double u;
  while (tau > tol) {
    r++;

    // exchange column r with column k
    l = piv[r];
    piv[r] = piv[k];
    piv[k] = l;

    for (i=0; i<m; i++) {
      u = A(i,r);
      A(i,r) = A(i,k);
      A(i,k) = u;
    }

    u = c[r];
    c[r] = c[k];
    c[k] = u;

    // find the Householder vector
    NumericVector v(m-r);
    for (i=0; i<m-r; i++) {
      v[i] = A(i+r,r);
    }
    v = house(v);

    // pre-multiply by the Householder matrix
    row_house(A, r, m-1, r, n-1, v);

    // update the sub-diagonal elements of column r
    for (i=1; i<m-r; i++) {
      A(i+r,r) = v[i];
    }

    // go to the next column and update the squared norm
    for (i=r+1; i<n; i++) {
      c[i] -= A(r,i)*A(r,i);
    }

    // identify the pivot column
    if (r < n-1) {
      tau = max(c[Range(r+1,n-1)]);
      for (k=r+1; k<n; k++) {
        if (c[k] > tol) break;
      }
    } else {
      tau = 0.0;
    }
  }

  // recover the Q matrix
  NumericMatrix Q = NumericMatrix::diag(m, 1.0);
  for (k=r; k>=0; k--) {
    NumericVector v(m-k);
    v[0] = 1.0;
    for (i=1; i<m-k; i++) {
      v[i] = A(i+k,k);
    }

    row_house(Q, k, m-1, k, m-1, v);
  }

  // recover the R matrix
  NumericMatrix R(m,n);
  for (j=0; j<n; j++) {
    for (i=0; i<=j; i++) {
      R(i,j) = A(i,j);
    }
  }

  List result = List::create(
    Named("qr") = A,
    Named("rank") = r+1,
    Named("pivot") = piv+1,
    Named("Q") = Q,
    Named("R") = R
  );

  return result;
}


// Given scalars a and b, this function computes
// c = cos(theta) and s = sin(theta) so that
//               |  c   s  |^T | a |  =  | r |
//               | -s   c  |   | b |     | 0 |
NumericVector givens(const double a, const double b) {
  double c, s, tau;

  if (b == 0.0) {
    c = 1.0; s = 0.0;
  } else {
    if (fabs(b) > fabs(a)) {
      double d = -std::copysign(1.0, b);
      tau = -a/b; s = d*1.0/sqrt(1.0 + tau*tau); c = s*tau;
    } else {
      double d = std::copysign(1.0, a);
      tau = -b/a; c = d*1.0/sqrt(1.0 + tau*tau); s = c*tau;
    }
  }

  return NumericVector::create(c,s);
}


// Given A in R^(2xq), c = cos(theta), and s = sin(theta),
// the following algorithm overwrites A with the matrix
//               |  c   s  |^T  A
//               | -s   c  |
void row_rot(NumericMatrix& A, const int i1, const int i2,
             const int j1, const int j2,
             const double c, const double s) {
  if (i1 < 0 || i1 >= i2 || i2 >= A.nrow()) {
    stop("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol()) {
    stop("Invalid column indices j1 and j2");
  }

  int q = j2-j1+1;
  for (int j=0; j<q; j++) {
    double tau1 = A(i1,j+j1);
    double tau2 = A(i2,j+j1);
    A(i1,j+j1) = c*tau1 - s*tau2;
    A(i2,j+j1) = s*tau1 + c*tau2;
  }
}


// Given A in R^(qx2), c = cos(theta), and s = sin(theta),
// the following algorithm overwrites A with the matrix
//               A  |  c   s  |
//                  | -s   c  |
void col_rot(NumericMatrix& A, const int i1, const int i2,
             const int j1, const int j2,
             const double c, const double s) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow()) {
    stop("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 >= j2 || j2 >= A.ncol()) {
    stop("Invalid column indices j1 and j2");
  }

  int q = i2-i1+1;
  for (int i=0; i<q; i++) {
    double tau1 = A(i+i1,j1);
    double tau2 = A(i+i1,j2);
    A(i+i1,j1) = c*tau1 - s*tau2;
    A(i+i1,j2) = s*tau1 + c*tau2;
  }
}


// Householder Bidiagonalization
// Given A in R^(mxn) with m>=n, the following algorithm overwrites
// the upper bidiagonal part of A with the upper bidiagonal part of
// t(U)*A*V = B, where B is upper bidiagonal and U = U_1 ... U_n and
// V = V_1 ... V_{n-2}. The essential part of U_j's Householder vector
// is stored in A((j+1):m, j), while the essential part of V_j's
// Householder vector is stored in A(j, (j+2):n).
List house_bidiag(NumericMatrix& A, const bool outtransform = 1) {
  int i, j, m = A.nrow(), n = A.ncol();
  if (m < n) {
    stop("The input matrix must have number of rows >= number of columns");
  }
  double tol = 1e-12;
  NumericMatrix B(n,n);
  NumericMatrix U = NumericMatrix::diag(m, 1.0);
  NumericMatrix V = NumericMatrix::diag(n, 1.0);

  bool bidiag = 1;
  for (i=0; i<n-2; i++) {
    for (j=i+2; j<n; j++) {
      if (fabs(A(i,j)) > tol) {
        bidiag = 0;
        break;
      }
    }
  }
  for (i=1; i<n; i++) {
    for (j=0; j<i; j++) {
      if (fabs(A(i,j)) > tol) {
        bidiag = 0;
        break;
      }
    }
  }
  for (i=n; i<m-1; i++) {
    for (j=0; j<n; j++) {
      if (fabs(A(i,j)) > tol) {
        bidiag = 0;
        break;
      }
    }
  }

  if (bidiag) {
    B = clone(A);
  } else {
    for (j=0; j<n; j++) {
      NumericVector v(m-j);
      for (i=0; i<m-j; i++) {
        v[i] = A(i+j,j);
      }
      v = house(v);

      row_house(A, j, m-1, j, n-1, v);

      // update the sub-diagonal elements of column j
      for (i=1; i<m-j; i++) {
        A(i+j,j) = v[i];
      }

      if (j < n-2) {
        NumericVector v(n-j-1);
        for (i=0; i<n-j-1; i++) {
          v[i] = A(j,i+j+1);
        }
        v = house(v);

        col_house(A, j, m-1, j+1, n-1, v);

        // update the elements of row j
        for (i=1; i<n-j-1; i++) {
          A(j,i+j+1) = v[i];
        }
      }
    }

    if (outtransform) {
      for (j=n-1; j>=0; j--) {
        NumericVector v(m-j);
        v[0] = 1.0;
        for (i=1; i<m-j; i++) {
          v[i] = A(i+j,j);
        }

        row_house(U, j, m-1, j, m-1, v);
      }

      for (j=n-3; j>=0; j--) {
        NumericVector v(n-j-1);
        v[0] = 1.0;
        for (i=1; i<n-j-1; i++) {
          v[i] = A(j,i+j+1);
        }

        row_house(V, j+1, n-1, j+1, n-1, v);
      }
    }

    for (j=0; j<n; j++) {
      B(j,j) = A(j,j);
      if (j<n-1) {
        B(j,j+1) = A(j,j+1);
      }
    }
  }

  if (outtransform) {
    return List::create(
      Named("B") = B,
      Named("U") = U,
      Named("V") = V
    );
  } else {
    return List::create(
      Named("B") = B
    );
  }
}


// Given a bidiagonal matrix with a zero diagonal, premultiplication
// by a sequence of Givens transformations to zero the entire row
List zero_diagonal(NumericMatrix& B, const int k,
                   const bool outtransform = 1) {
  int j, n = B.nrow();
  if (B.ncol() != n) {
    stop("The input matrix must be a square matrix");
  }
  if (k < 0 || k >= n-1) {
    stop("Invalid value for index k");
  }
  NumericMatrix U = NumericMatrix::diag(n, 1.0);

  for (j=k+1; j<n; j++) {
    NumericVector v = givens(B(k,j), B(j,j));
    double w = v[0];
    v[0] = -v[1]; v[1] = w;
    int j1 = j < n-1 ? j+1 : n-1;
    row_rot(B, k, j, j, j1, v[0], v[1]);
    if (outtransform) col_rot(U, k, j, k, j, v[0], v[1]);
  }

  if (outtransform) {
    return List::create(
      Named("B") = B,
      Named("U") = U
    );
  } else {
    return List::create(
      Named("B") = B
    );
  }
}


// Golub-Kahan SVD Step
// Given a bidiagonal matrix B having no zeros on its diagonal or
// superdiagonal, the following algorithm overwrites B with the
// bidiagonal matrix t(U)*B*V, where U and V are orthogonal and V
// is essentially the orthogonal matrix that would be obtained by
// applying Algorithm 8.2.2 in Golub and Van Loan (1989) to T = t(B)*B.
List svd_step(NumericMatrix& B, const bool outtransform = 1) {
  int k, n = B.ncol();
  NumericMatrix U = NumericMatrix::diag(n, 1.0);
  NumericMatrix V = NumericMatrix::diag(n, 1.0);

  double f1 = B(n-3,n-2), f2 = B(n-2,n-1);
  double d1 = B(n-2,n-2), d2 = B(n-1,n-1);
  double a1 = f1*f1 + d1*d1, a2 = f2*f2 + d2*d2, b1 = f2*d1;
  double d = 0.5*(a1-a2);
  double mu = a2 + d - std::copysign(1.0, d)*sqrt(d*d + b1*b1);
  double y = B(0,0)*B(0,0) - mu;
  double z = B(0,0)*B(0,1);
  NumericVector v(2);
  for (k=0; k<n-1; k++) {
    v = givens(y,z);
    int k1 = k > 0 ? k-1 : 0;
    col_rot(B, k1, k+1, k, k+1, v[0], v[1]);
    if (outtransform) col_rot(V, 0, k+1, k, k+1, v[0], v[1]);

    y = B(k,k);
    z = B(k+1,k);
    v = givens(y,z);
    int k2 = k < n-2 ? k+2 : n-1;
    row_rot(B, k, k+1, k, k2, v[0], v[1]);
    if (outtransform) col_rot(U, 0, k+1, k, k+1, v[0], v[1]);

    if (k < n-2) {
      y = B(k,k+1);
      z = B(k,k+2);
    }
  }

  if (outtransform) {
    return List::create(
      Named("B") = B,
      Named("U") = U,
      Named("V") = V
    );
  } else {
    return List::create(
      Named("B") = B
    );
  }
}


//' @title Singular Value Decomposition of a Matrix
//' @description Computes the singular-value decomposition of a
//' rectangular matrix.
//'
//' @param X A numeric matrix whose SVD decomposition is to be computed.
//' @param outtransform Whether the orthogonal matrices composing of the
//'   left and right singular vectors are to be computed.
//' @param decreasing Whether the singular values should be sorted in
//'   decreasing order and the corresponding singular vectors rearranged
//'   accordingly.
//'
//' @details
//' Given \eqn{A \in R^{m\times n} (m \geq n)}, the following algorithm
//' overwrites \eqn{A} with \eqn{U^T A V = D}, where
//' \eqn{U\in R^{m\times m}} is orthogonal, \eqn{V \in R^{n\times n}} is
//' orthogonal, and \eqn{D \in R^{m\times n}} is diagonal.
//'
//' @return A list with the following components:
//'
//' * \code{d}: A vector containing the singular values of \eqn{X}.
//'
//' * \code{U}: A matrix whose columns contain the left singular vectors
//'   of \eqn{X}.
//'
//' * \code{V}: A matrix whose columns contain the right singular vectors
//'   of \eqn{X}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Gene N. Golub and Charles F. Van Loan.
//' Matrix Computations, second edition. Baltimore, Maryland:
//' The John Hopkins University Press, 1989, p.434.
//'
//' @examples
//'
//' A <- matrix(c(1,0,0,0, 1,2,0,0, 0,1,3,0, 0,0,1,4), 4, 4)
//' svdcpp(A)
//'
//' @export
// [[Rcpp::export]]
List svdcpp(const NumericMatrix& X, const bool outtransform = 1,
            const bool decreasing = 1) {
  int i, j, k, l, m1 = X.nrow(), n1 = X.ncol(), m, n;
  if (m1 >= n1) {
    m = m1; n = n1;
  } else {
    m = n1; n = m1;
  }
  NumericMatrix Y(m,n);
  NumericMatrix U = NumericMatrix::diag(m, 1.0);
  NumericMatrix V = NumericMatrix::diag(n, 1.0);
  if (m1 >= n1) {
    Y = clone(X);
  } else {
    Y = clone(transpose(X));
  }
  double tol = 1e-12;

  List a = house_bidiag(Y, outtransform);
  NumericMatrix B = a["B"];
  if (outtransform) {
    U = as<NumericMatrix>(a["U"]);
    V = as<NumericMatrix>(a["V"]);
  }

  int p, q = 0;
  while (q < n) {
    for (i=1; i<n; i++) {
      if (fabs(B(i-1,i)) <= tol*(fabs(B(i-1,i-1)) + fabs(B(i,i)))) {
        B(i-1,i) = 0.0;
      }
    }

    // find the largest non-negative q and the smallest non-negative p
    // such that
    //               |  B11   0    0   |   p
    //           B = |   0   B22   0   |   n-p-q
    //               |   0    0   B33  |   q
    //                   p  n-p-q  q
    // where B33 is diagonal and B22 has nonzero superdiagonal
    q = n;
    for (i=n-1; i>=1; i--) {
      if (B(i-1,i) != 0.0) {
        q = n-i-1;
        break;
      }
    }

    p = 0;
    for (i=n-q-2; i>=1; i--) {
      if (B(i-1,i) == 0.0) {
        p = i;
        break;
      }
    }

    if (q < n) {
      // if any diagonal entry in B22 is zero, then zero the superdiagonal
      // entry in the same row
      NumericMatrix B22 = B(Range(p,n-q-1), Range(p,n-q-1));
      for (i=0; i<n-p-q-1; i++) {
        if (fabs(B22(i,i)) < tol) {
          List b = zero_diagonal(B22, i, outtransform);
          if (outtransform) {
            NumericMatrix Z = b["U"];
            NumericMatrix W = U(Range(0,m-1), Range(p,n-q-1));
            for (j=0; j<m; j++) {
              for (k=0; k<n-p-q; k++) {
                U(j,k+p) = 0.0;
                for (l=0; l<n-p-q; l++) {
                  U(j,k+p) += W(j,l)*Z(l,k);
                }
              }
            }
          }
        }
      }

      // apply Algorithm 8.3.1 to B22
      List c = svd_step(B22, outtransform);

      // update B22
      for (i=0; i<n-p-q; i++) {
        for (j=0; j<n-p-q; j++) {
          B(i+p,j+p) = B22(i,j);
        }
      }

      if (outtransform) {
        NumericMatrix Z1 = c["U"];
        NumericMatrix W1 = U(Range(0,m-1), Range(p,n-q-1));
        for (i=0; i<m; i++) {
          for (j=0; j<n-p-q; j++) {
            U(i,j+p) = 0.0;
            for (k=0; k<n-p-q; k++) {
              U(i,j+p) += W1(i,k)*Z1(k,j);
            }
          }
        }

        NumericMatrix Z2 = c["V"];
        NumericMatrix W2 = V(Range(0,n-1), Range(p,n-q-1));
        for (i=0; i<n; i++) {
          for (j=0; j<n-p-q; j++) {
            V(i,j+p) = 0.0;
            for (k=0; k<n-p-q; k++) {
              V(i,j+p) += W2(i,k)*Z2(k,j);
            }
          }
        }
      }
    }
  }

  NumericVector d(n);
  for (i=0; i<n; i++) {
    d[i] = B(i,i);
  }

  // ensure the singular values are positive
  for (i=0; i<n; i++) {
    if (d[i] < 0.0) {
      d[i] = -d[i];
      V(_,i) = -V(_,i);
    }
  }

  if (decreasing) {
    // order the singular values from the largest to the smallest
    // and the arrange the associated vectors accordingly
    IntegerVector order = seq(0, n-1);
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return d[i] > d[j];
    });
    d = d[order];
    if (outtransform) {
      NumericMatrix Z = clone(U);
      NumericMatrix W = clone(V);
      for (i=0; i<n; i++) {
        U(_,i) = Z(_,order[i]);
        V(_,i) = W(_,order[i]);
      }
    }
  }

  // switch U and V if m1 < n1
  NumericMatrix U1(m1,m1), V1(n1,n1);
  if (m1 >= n1) {
    U1 = U;
    V1 = V;
  } else {
    U1 = V;
    V1 = U;
  }

  if (outtransform) {
    return List::create(
      Named("d") = d,
      Named("U") = U1,
      Named("V") = V1
    );
  } else {
    return List::create(
      Named("d") = d
    );
  }
}


NumericMatrix rmvnorm(int n, NumericVector mean, NumericMatrix sigma) {
  int i,j,k;
  int p = static_cast<int>(mean.size());
  double toler = 1.818989e-12;
  NumericMatrix v = clone(sigma);
  i = cholesky2(v, p, toler);

  NumericMatrix H(p,p);
  for (i=0; i<p; i++) {
    H(i,i) = sqrt(v(i,i));
    for (j=0; j<i; j++) {
      H(i,j) = v(j,i)*H(j,j);
    }
  }

  NumericMatrix result(n,p);
  NumericVector z(p);
  for (i=0; i<n; i++) {
    for (j=0; j<p; j++) {
      z[j] = R::rnorm(0,1);
    }

    for (j=0; j<p; j++) {
      result(i,j) = mean[j];
      for (k=0; k<p; k++) {
        result(i,j) += H(j,k)*z[k];
      }
    }
  }

  return result;
}


//' @title Converting a decimal to a fraction
//' @description Converts a decimal to a fraction based on the algorithm
//' from http://stackoverflow.com/a/5128558/221955.
//'
//' @param x The fraction in decimal form.
//' @param tol The tolerance level for the conversion error.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' float_to_fraction(5/3)
//'
//' @export
// [[Rcpp::export]]
NumericVector float_to_fraction(const double x, const double tol=0.000001) {
  NumericVector v(2);
  double x1 = x;
  double n = std::floor(x1);
  x1 = x1 - n;
  if (x1 < tol) {
    v[0] = n;
    v[1] = 1;
  } else if (1 - tol < x1) {
    v[0] = n+1;
    v[1] = 1;
  } else {
    // The lower fraction is 0/1
    double lower_n = 0;
    double lower_d = 1;

    // The upper fraction is 1/1
    double upper_n = 1;
    double upper_d = 1;

    bool cond = 1;
    while (cond) {
      // The middle fraction is (lower_n + upper_n) / (lower_d + upper_d)
      double middle_n = lower_n + upper_n;
      double middle_d = lower_d + upper_d;

      // If x + tol < middle
      if (middle_d * (x1 + tol) < middle_n) {
        // middle is our new upper
        upper_n = middle_n;
        upper_d = middle_d;
      } else if (middle_n < (x1 - tol) * middle_d) {
        // Else If middle < x - tol
        // middle is our new lower
        lower_n = middle_n;
        lower_d = middle_d;
      } else {
        // Else middle is our best fraction
        v[0] = n * middle_d + middle_n;
        v[1] = middle_d;
        cond = 0;
      }
    }
  }

  return v;
}
