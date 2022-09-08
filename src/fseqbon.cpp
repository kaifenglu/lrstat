#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;


//' @title Get efficacy boundaries for group sequential design
//' @description Obtains the efficacy stopping boundaries for a generic group
//' sequential design.
//'
//' @inheritParams param_kMax
//' @inheritParams param_informationRates
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//'
//' @return A numeric vector of critical values.
//'
//' @examples
//'
//' getBound(kMax = 2, informationRates = c(0.5,1),
//'          alpha = 0.025, typeAlphaSpending = "sfOF")
//'
//' @export
// [[Rcpp::export]]
NumericVector getBound(
    const int kMax = NA_INTEGER,
    const NumericVector& informationRates = NA_REAL,
    const double alpha = 0.025,
    const String typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const NumericVector& userAlphaSpending = NA_REAL) {
  
  NumericVector informationRates1 = clone(informationRates);
  
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
    informationRates1 = as<NumericVector>(tem)/(kMax+0.0);
  }
  
  String asf = typeAlphaSpending;
  double asfpar = parameterAlphaSpending;
  
  if (asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegnative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }
  
  
  if (asf=="WT" && R_isnancpp(asfpar)) {
    stop("Missing parameter for WT");
  }
  
  if (asf=="sfKD" && R_isnancpp(asfpar)) {
    stop("Missing parameter for sfKD");
  }
  
  if (asf=="sfHSD" && R_isnancpp(asfpar)) {
    stop("Missing parameter for sfHSD");
  }
  
  if (asf=="sfKD" && asfpar <= 0) {
    stop ("asfpar must be positive for sfKD");
  }
  
  NumericVector theta(kMax); // mean values under H0, initialized to zero
  NumericVector t = clone(informationRates1); // information time
  NumericVector criticalValues(kMax);
  
  if (asf == "none") {
    for (int i=0; i<kMax-1; i++) {
      criticalValues[i] = 6.0;
    }
    criticalValues[kMax-1] = R::qnorm(1-alpha, 0, 1, 1, 0);
  } else if (asf == "OF" || asf == "P" || asf == "WT") {
    double Delta;
    if (asf == "OF") {
      Delta = 0;
    } else if (asf == "P") {
      Delta = 0.5;
    } else {
      Delta = asfpar;
    }
    
    auto f = [kMax, alpha, Delta, theta, t] (double aval)->double {
      NumericVector u(kMax), l(kMax);
      for (int i=0; i<kMax; i++) {
        u[i] = aval*pow(t[i], Delta-0.5);
        l[i] = -6.0;
      }
      
      List probs = exitprob(u, l, theta, t);
      double cpu = sum(NumericVector(probs[0]));
      return cpu - alpha;
    };
    
    double cwt = brent(f, 0, 10, 1e-6);
    for (int i=0; i<kMax; i++) {
      criticalValues[i] = cwt*pow(t[i], Delta-0.5);
    }
  } else if (asf == "sfOF" || asf == "sfP" || asf == "sfKD" ||
    asf == "sfHSD" || asf == "user") {
    
    // stage 1
    double cumAlphaSpent;
    if (asf == "user") {
      cumAlphaSpent = userAlphaSpending[0];
    } else {
      cumAlphaSpent = errorSpent(t[0], alpha, asf, asfpar);
    }
    
    criticalValues[0] = R::qnorm(1 - cumAlphaSpent, 0, 1, 1, 0);
    
    
    // lambda expression for finding the critical Values at stage k
    int k=0;
    auto f = [&k, &cumAlphaSpent, &criticalValues,
              theta, t](double aval)->double {
                NumericVector u(k+1), l(k+1);
                for (int i=0; i<k; i++) {
                  u[i] = criticalValues[i];
                  l[i] = -6.0;
                }
                u[k] = aval;
                l[k] = -6.0;
                
                IntegerVector idx = Range(0,k);
                List probs = exitprob(u, l, theta[idx], t[idx]);
                double cpu = sum(NumericVector(probs[0]));
                return cpu - cumAlphaSpent;
              };
    
    // subsequent stages
    for (k=1; k<kMax; k++) {
      if (asf == "user") {
        cumAlphaSpent = userAlphaSpending[k];
      } else {
        cumAlphaSpent = errorSpent(t[k], alpha, asf, asfpar);
      }
      
      if (f(6) > 0) { // no alpha spent at current visit
        criticalValues[k] = 6.0;
      } else {
        criticalValues[k] = brent(f, 0, 6, 1e-6);
      }
    }
  } else {
    stop("Invalid type of alpha spending");
  }
  
  return criticalValues;
}


//' @title Group sequential trials using graphical approaches
//' @description Obtains the test results for group sequential trials using 
//' graphical approaches based on weighted Bonferroni tests.
//'
//' @param w The vector of initial weights for elementary hypotheses.
//' @param G The initial transition matrix.
//' @param alpha The significance level. Defaults to 0.025.
//' @param asf The vector of alpha spending functions. Each element is one 
//'   of the following: "OF" for O'Brien-Fleming boundaries, 
//'   "P" for Pocock boundaries, "WT" for Wang & Tsiatis boundaries, 
//'   "sfOF" for O'Brien-Fleming type spending function, 
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function, 
//'   "sfHSD" for Hwang, Shi & DeCani spending function, 
//'   and "none" for no early efficacy stopping.
//' @param asfpar The vector of parameter values for the alpha spending 
//'   functions. Each element corresponds to the value of Delta for "WT", 
//'   rho for "sfKD", or gamma for "sfHSD".
//' @param incid The incidence matrix indicating whether the specific 
//'   hypothesis will be tested at the given look.
//' @param maxinfo The vector of targeted maximum information for each 
//'   hypothesis.
//' @param info The matrix of observed information for each hypothesis by look.
//' @param p The matrix of raw p-values for each hypothesis by look. 
//' 
//' @return A logical matrix to indicate whether the specific hypothesis 
//'   is rejected at the given look.
//'
//' @references
//' Willi Maurer and Frank Bretz. Multiple testing in group sequential
//' trials using graphical approaches. Statistics in Biopharmaceutical 
//' Research. 2013;5:311-320. \doi{10.1080/19466315.2013.807748}
//'
//' @examples
//' 
//' # Case study from Maurer & Bretz (2013) 
//' 
//' fseqbon(
//'   w = c(0.5, 0.5, 0, 0),
//'   G = matrix(c(0, 0.5, 0.5, 0,  0.5, 0, 0, 0.5,  
//'                0, 1, 0, 0,  1, 0, 0, 0), 
//'              nrow=4, ncol=4, byrow=TRUE),
//'   alpha = 0.025,
//'   asf = rep("sfOF", 4),
//'   incid = matrix(1, nrow=4, ncol=3),
//'   maxinfo = rep(1,4),
//'   info = matrix(c(rep(1/3, 4), rep(2/3, 4)), nrow=4, ncol=2),
//'   p = matrix(c(0.0062, 0.017, 0.009, 0.13, 
//'                0.0002, 0.0035, 0.002, 0.06), 
//'              nrow=4, ncol=2))
//' 
//'
//' @export
// [[Rcpp::export]]
LogicalMatrix fseqbon(const NumericVector& w, 
                      const NumericMatrix& G, 
                      const double alpha = 0.025,
                      const StringVector& asf = NA_STRING, 
                      const NumericVector& asfpar = NA_REAL, 
                      const LogicalMatrix& incid = NA_LOGICAL, 
                      const NumericVector& maxinfo = NA_REAL, 
                      const NumericMatrix& info = NA_REAL, 
                      const NumericMatrix& p = NA_REAL) {
  
  int i, j, k, l, m = w.size(), K = incid.ncol(), B = info.nrow()/m;
  LogicalMatrix reject(B*m, K);  // rejection indicators

  NumericVector wx(m);    // dynamic weights
  NumericMatrix g(m,m);   // dynamic transition matrix
  NumericMatrix g1(m,m);  // temporary transition matrix 
  
  NumericVector asfpars = clone(asfpar);
  
  if (G.nrow() != m || G.ncol() != m) {
    stop("Invalid dimension for G");
  }
  
  if (is_true(any(w < 0.0))) {
    stop("w must be nonnegative");
  }
  
  if (sum(w) != 1.0) {
    stop("w must sum to 1");
  }
  
  if (is_true(any(G < 0.0))) {
    stop("G must be nonnegative");
  }
  
  if (is_true(any(rowSums(G) > 1.0 + 1.0e-8))) {
    stop("Row sums of G must be less than or equal to 1");
  }
  
  for (i=0; i<m; i++) {
    if (G(i,i) != 0.0) {
      stop("Diagonal elements of G must be equal to 0");
    }
  }
  
  if (asf.size() != m) {
    stop("Invalid length for asf");
  }
  
  if (is_true(any(is_na(asfpar))) && asfpar.size()==1) {
    asfpars = rep(NA_REAL, m);
  }
  
  if (asfpars.size() != m) {
    stop("Invalid length for asfpar");
  }

  for (i=0; i<m; i++) {
    if (!(asf(i)=="OF" || asf(i)=="P" || asf(i)=="WT" || 
        asf(i)=="sfOF" || asf(i)=="sfP" || asf(i)=="sfKD" || 
        asf(i)=="sfHSD" || asf(i)=="none")) {
      stop("Invalid type for asf");
    }
    
    if (asf(i)=="WT" && R_isnancpp(asfpars(i))) {
      stop("Missing parameter for WT");
    }
    
    if (asf(i)=="sfKD" && R_isnancpp(asfpars(i))) {
      stop("Missing parameter for sfKD");
    }
    
    if (asf(i)=="sfHSD" && R_isnancpp(asfpars(i))) {
      stop("Missing parameter for sfHSD");
    }
    
    if (asf(i)=="sfKD" && asfpars(i) <= 0) {
      stop ("asfpar must be positive for sfKD");
    }
  }
  
  if (incid.nrow() != m) {
    stop("Invalid number of rows for incid");
  }
  
  if (maxinfo.size() != m) {
    stop("Invalid length for maxinfo");
  }
  
  if (is_true(any(maxinfo <= 0.0))) {
    stop("maxinfo must be positive");
  }
  
  if (info.nrow() % m != 0) {
    stop("Invalid number of rows for info");
  }
  
  if (info.ncol() > K) {
    stop("Invalid number of columns for info");
  }
  
  if (p.nrow() % m != 0) {
    stop("Invalid number of rows for p");
  }
  
  if (p.ncol() > K) {
    stop("Invalid number of columns for p");
  }
  
  if (info.nrow() != p.nrow()) {
    stop("info and p must have the same number of rows");
  }
  
  if (info.ncol() != p.ncol()) {
    stop("info and p must have the same number of columns");
  }
  
  
  IntegerVector K0 = rowSums(incid); // maximum number or stages
  
  IntegerMatrix idx1(m, K); // study look
  IntegerMatrix idx2(m, K); // hypothesis look
  for (j=0; j<m; j++) {
    l = 0;
    for (k=0; k<K; k++) {
      if (incid(j,k)==1) {
        idx1(j,l) = k;
        idx2(j,k) = l;
        l++;
      } else {
        idx2(j,k) = NA_INTEGER;
      }
    }
    for (k=l; k<K; k++) {
      idx1(j,k) = NA_INTEGER;
    }
  }
  
  
  int step, kj;
  double alphaj, alphastar, asfpar1;
  String asf1;
  int K3, K4 = 0;
  
  NumericMatrix info1(m, K), p1(m, K);
  IntegerVector K1(m), K2(m);
  LogicalVector r(m);  
  
  for (int iter=0; iter<B; iter++) {
    for (j=0; j<m; j++) {
      NumericVector irow = info.row(iter*m+j);
      NumericVector prow = p.row(iter*m+j);
      
      irow = irow[!is_na(irow)]; // exclude missing values
      prow = prow[!is_na(prow)];
      
      if (irow.size() != prow.size()) {
        stop("info and p should have the same number of nonmissing elements");
      }
      
      K1(j) = irow.size();
      if (K1(j) > K0(j)) {
        K1(j) = K0(j); // truncated number of looks for hypothesis j 
      }
      
      K2(j) = idx1(j, K1(j)-1) + 1;  // last study look for hypothesis j
      
      info1(j, _) = irow;
      p1(j, _) = prow;
    }

    K3 = max(K2);
    if (K3 > K4) {
      K4 = K3;
    }
    
    wx = clone(w);  // reset
    g = clone(G); 
    r.fill(0);
    
    for (step=0; step<K3; step++) {  // loop over study look
      for (j=0; j<m; j++) {
        if (!incid(j, step)) {
          reject(iter*m+j, step) = NA_LOGICAL;
        } 
      }
      
      
      for (j=0; j<m; j++) {
        // testable non-rejected hypotheses with positive weights
        if (incid(j, step) && !r(j) && wx(j) > 1.0e-4) { 
          asf1 = Rcpp::String(asf(j));
          asfpar1 = asfpars(j);
          
          // calculate the nominal significance level 
          kj = idx2(j, step);  // index of current look for hypothesis j
          alphaj = wx(j)*alpha;   // alpha allocated to hypothesis j
          
          
          if (kj < K0(j)-1 && info1(j, kj) < maxinfo(j)) {
            NumericVector t(kj+2);  // make the next look the last look
            for (i=0; i<=kj; i++) {
              t(i) = info1(j, i)/maxinfo(j);
            } 
            t(kj+1) = 1;
            
            NumericVector u = getBound(kj+2, t, alphaj, asf1, asfpar1, 0); 
            alphastar = 1 - R::pnorm(u(kj), 0, 1, 1, 0);
          } else {
            NumericVector t(kj+1);  // current look is the last look
            
            // critical values for previous looks are w.r.t. target max info
            for (i=0; i<kj; i++) {
              t(i) = info1(j, i)/maxinfo(j);
            } 
            t(kj) = 1;
            
            NumericVector u = getBound(kj+1, t, alphaj, asf1, asfpar1, 0); 
            
            // update the critical value at the last look if the observed
            // maximum information is different from the target
            if (abs(info1(j, kj) - maxinfo(j)) > 1e-12) {
              NumericVector infoj = info1(j, _);
              auto f = [u, alphaj, infoj, kj](double aval)->double {
                NumericVector u1 = clone(u);
                u1[kj] = aval;
                NumericVector l = rep(-6.0, kj+1);
                NumericVector theta = rep(0.0, kj+1);
                NumericVector I = infoj[Range(0,kj)];
                List p2 = exitprob(u1, l, theta, I);
                return sum(NumericVector(p2[0])) - alphaj;
              };
              
              u(kj) = brent(f, 0, 6, 1e-6);
            }
            
            alphastar = 1 - R::pnorm(u(kj), 0, 1, 1, 0);
          }
          
          
          if (p1(j, kj) < alphastar) {
            // reject Hj
            r(j) = 1;
            reject(iter*m+j, step) = 1;
            
            // update weights
            for (l=0; l<m; l++) {
              if (r(l) == 0) {
                wx(l) = wx(l) + wx(j)*g(j,l);
              }
            }
            wx(j) = 0.0;
            
            // update transition matrix
            g1.fill(0.0);
            for (l=0; l<m; l++) {
              if (r[l] == 0) {
                for (k=0; k<m; k++) {
                  if ((r[k] == 0) && (l != k) && 
                      (g(l,j)*g(j,l) < 1.0 - 1.0e-12)) {
                    g1(l,k) = (g(l,k) + g(l,j)*g(j,k))/(1.0 - g(l,j)*g(j,l));
                  }
                }
              }
            }
            g = clone(g1);
            
            j = -1;  // restart the loop over the index set 
          }
        }
      }
      
      // stop if all hypotheses have been rejected
      if (sum(r) == m) {
        step = K3;
      }
    }
  }
  
  return reject(_, Range(0, K4-1));
}

