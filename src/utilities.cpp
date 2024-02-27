#include <Rcpp.h>
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
//' if \code{i <- findInterval2(x,v)}, for each index \code{j} in \code{x},
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



// [[Rcpp::export]]
double errorSpentcpp(const double t, const double error,
                     const String sf, const double sfpar) {
  if (error <= 0 || error >= 1) {
    stop("error must be a number between 0 and 1");
  }
  if (t <= 0 || t > 1) {
    stop("t must be a number between 0 and 1");
  }

  std::string asf = sf;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
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
  int kMax=b.size(), r1=6*r-1, r2=12*r-3, i0, i1=0, i2=r1-1, i, j,
    m0=r2, m1=r1, m=r2;
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
    I1 = as<NumericVector>(tem);
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
                         const double lowerBound) {
  int n = q.size();
  NumericVector p(n);
  for (int h=0; h<n; h++) {
    if (q[h] <= lowerBound) {
      p[h] = 0;
    } else {
      NumericVector y = NumericVector::create(lowerBound, q[h]);
      IntegerVector i = findInterval2(y, piecewiseSurvivalTime);
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

  return p;
}


// [[Rcpp::export]]
double qtpwexpcpp1(const double p,
                   const NumericVector& piecewiseSurvivalTime,
                   const NumericVector& lambda,
                   const double lowerBound) {
  int j, j1, m = piecewiseSurvivalTime.size();
  double q, v, v1;

  // cumulative hazard from lowerBound until the quantile
  v1 = -log(1 - p);

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
                         const double lowerBound) {
  int n = p.size();
  NumericVector q(n);
  for (int h=0; h<n; h++) {
    q[h] = qtpwexpcpp1(p[h], piecewiseSurvivalTime, lambda, lowerBound);
  }

  return q;
}


// [[Rcpp::export]]
NumericVector rtpwexpcpp(const int n,
                         const NumericVector& piecewiseSurvivalTime,
                         const NumericVector& lambda,
                         const double lowerBound) {
  NumericVector p(n);
  for (int i=0; i<n; i++) {
    p[i] = R::runif(0,1);
  }

  return qtpwexpcpp(p, piecewiseSurvivalTime, lambda, lowerBound);
}


// [[Rcpp::export]]
NumericVector getBoundcpp(const int k,
                          const NumericVector& informationRates,
                          const double alpha,
                          const String typeAlphaSpending,
                          const double parameterAlphaSpending,
                          const NumericVector& userAlphaSpending,
                          const NumericVector& spendingTime,
                          const LogicalVector& efficacyStopping) {

  NumericVector informationRates1 = clone(informationRates);
  NumericVector spendingTime1 = clone(spendingTime);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);

  if (R_isnancpp(k)) {
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
    informationRates1 = as<NumericVector>(tem)/(k+0.0);
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
    c = std::tolower(c);
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
    NumericVector informationRates2 = as<NumericVector>(tem)/(k+0.0);
    if (max(abs(informationRates1 - informationRates2)) > 1e-6) {
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
  NumericVector I = as<NumericVector>(tem);
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

    double cwt = brent(f, 0.0, 10.0, 1e-6);
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
                    a[k] = brent(g, -6.0, b[k], 1e-6);
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
    beta = brent(f, 0.0001, 1-alpha, 1e-6);
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


// polynomial interpolation
// [[Rcpp::export]]
NumericVector polint(const NumericVector& xa,
                     const NumericVector& ya,
                     int n, double x) {
  int i, m, ns=0;
  double den, dif, dift, ho, hp, w;
  NumericVector c(n), d(n);

  dif = fabs(x - xa[0]);
  for (i=0; i<n; i++) {
    // find the index ns of the closest table entry
    if ((dift = fabs(x - xa[i])) < dif) {
      ns = i;
      dif = dift;
    }

    // initialize the tableau of c's and d's
    c[i] = ya[i];
    d[i] = ya[i];
  }

  // The initial approximation to y = f(x)
  double y = ya[ns], dy = 0;

  for (m=1; m<n; m++) { // for each column of the tableau
    for (i=0; i<n-m; i++) { // loop through and update current c's and d's
      ho = xa[i] - x;
      hp = xa[i+m] - x;
      w = c[i+1] - d[i];
      if ((den = ho - hp) == 0.0) {
        stop("Error in routine polint");
      } // this error can occur only if two input xa's are identical

      den = w/den;
      d[i] = hp*den; // update the c's and d's
      c[i] = ho*den;
    }

    // After each column in the tableau is completed, decide which
    // correction, c or d, we want to add to our accumulating value of y,
    // i.e., which path to take through the tableau -- forking up or down.
    // We do this in such a way to take the most ''straight line'' route
    // through the tableau to its apex, updating ns accordingly to keep
    // track of where we are.
    y += (dy = (2*ns < (n-m) ? c[ns] : d[--ns]));
  }

  NumericVector results = NumericVector::create(y, dy);
  return results;
}


// trapezoidal rule
// [[Rcpp::plugins(cpp11)]]
double trapzd(const std::function<double(double)>& f,
              double a, double b, int n) {
  double x, tnm, sum, del;
  static double s;
  int it, j;

  if (n == 1) {
    return (s = 0.5*(b-a)*(f(a) + f(b)));
  } else {
    for (it=1, j=1; j<n-1; j++) it <<= 1;
    tnm = it;
    del = (b-a)/tnm; // spacing of the points to be added
    x = a + 0.5*del;
    for (sum=0.0, j=1; j<=it; j++, x+=del) sum += f(x);
    s = 0.5*(s + sum*del); // replace s with its refined value
    return s;
  }
}


#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

//' @title Romberg integration
//' @description Integration by Romberg's method of order 2K,
//' where, e.g., K=2 is Simpson's rule.
//' Program based on the book - Numerical Recipes in C
//' The Art of Scientific Computing - Second Edition, by William H. Press,
//' Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery.
//'
//' @param f Name of the univariate integrand function.
//' @param a One end of the interval.
//' @param b The other end of the interval.
//'
//' @return The integral of the function f from a to b.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' qromb(sin, 0, pi)
//'
//' @export
//'
// [[Rcpp::plugins(cpp11)]]
double qromb(const std::function<double(double)>& f,
             double a, double b) {

  // successive trapezoidal approximations and their relative step sizes
  NumericVector s = NumericVector(JMAXP, NA_REAL);
  NumericVector h = NumericVector(JMAXP+1, NA_REAL);
  NumericVector xa = NumericVector(JMAXP, NA_REAL);
  NumericVector ya = NumericVector(JMAXP, NA_REAL);

  h[0] = 1.0;
  for (int j=0; j<JMAX; j++) {
    s[j] = trapzd(f, a, b, j);
    if (j >= K-1) {
      IntegerVector q = Range(j-K+1, JMAX);
      xa = h[q];
      ya = s[q];
      NumericVector ss = polint(xa, ya, K, 0.0);
      if (fabs(ss[1]) <= EPS*fabs(ss[0])) return ss[0];
    }
    h[j+1] = 0.25*h[j];
  }
  stop("Too many steps in routine qromb");
  return 0.0; // Never get here
}


#include <algorithm>
#define ITMAX 100
#define EPS 3.0e-8
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
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



