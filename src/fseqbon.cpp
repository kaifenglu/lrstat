#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;


//' @title Update graph for graphical approaches
//' @description Updates the weights and transition matrix for graphical 
//' approaches.
//'
//' @param w The current vector of weights for elementary hypotheses.
//' @param G The current transition matrix.
//' @param I The set of indices for yet to be rejected hypotheses.
//' @param j The hypothesis to remove from index set \code{I}.
//'
//' @return A list containing the new vector of weights and the new 
//' transition matrix for the graph, and the new set of indices of yet 
//' to be rejected hypotheses.
//'
//' @examples
//' updateGraph(w = c(0.5, 0.5, 0, 0), 
//'             G = matrix(c(0, 0.5, 0.5, 0,  0.5, 0, 0, 0.5,  
//'                          0, 1, 0, 0,  1, 0, 0, 0), 
//'                        nrow=4, ncol=4, byrow=TRUE), 
//'             I = c(1, 2, 3, 4), 
//'             j = 1)
//'
//' @export
// [[Rcpp::export]]
List updateGraph(const NumericVector& w, const NumericMatrix& G, 
                 const IntegerVector& I, const int j) {
  int k, l, m = w.size();
  
  if (G.nrow() != m || G.ncol() != m) {
    stop("Invalid dimension for G");
  }
  
  if (min(I) < 1 || max(I) > m) {
    stop("Elements of I must be integers between 1 and m");
  }
  
  if (is_true(any(duplicated(I)))) {
    stop("The index set I must not contain duplicates");
  }
  
  if (std::find(I.begin(), I.end(), j) == I.end()) {
    stop("j must be in I");
  }
  
  int j1 = j-1;
  IntegerVector I1 = I-1;
  
  LogicalVector r(m);
  r.fill(1);
  r[I1] = 0;
  r(j1) = 1;
  
  // update weights
  NumericVector wx = clone(w);
  for (l=0; l<m; l++) {
    if (r(l) == 0) {
      wx(l) = wx(l) + wx(j1)*G(j1,l);
    }
  }
  wx(j1) = 0.0;
  
  // update transition matrix
  NumericMatrix g(m,m);
  for (l=0; l<m; l++) {
    if (r[l] == 0) {
      for (k=0; k<m; k++) {
        if ((r[k] == 0) && (l != k) && (G(l,j1)*G(j1,l) < 1.0 - 1.0e-12)) {
          g(l,k) = (G(l,k) + G(l,j1)*G(j1,k))/(1.0 - G(l,j1)*G(j1,l));
        }
      }
    }
  }
  
  List result = List::create(
    _["w"] = wx,
    _["G"] = g,
    _["I"] = I[I!=j]);
  
  return result;
  
}



// [[Rcpp::export]]
NumericMatrix fadjpboncpp(const NumericVector& w, 
                          const NumericMatrix& G, 
                          const NumericMatrix& p) {
  
  int i, j, k, l, m = w.size(), iter, iters = p.nrow(), step;
  double pmax; // running maximum adjusted p-value
  
  NumericMatrix padj(iters,m);  // adjusted p-values
  LogicalVector r(m);  // rejection indicators
  NumericVector pvalues(m);  // raw p-values
  NumericVector q(m);  // ratios of raw p-values over weights
  
  NumericVector wx(m);    // dynamic weights
  NumericMatrix g(m,m);  // dynamic transition matrix
  NumericMatrix g1(m,m);  // temporary transition matrix 
  
  
  if (is_true(any(w < 0.0))) {
    stop("w must be nonnegative");
  }
  
  if (sum(w) != 1.0) {
    stop("w must sum to 1");
  }
  
  if (G.nrow() != m || G.ncol() != m) {
    stop("Invalid dimension for G");
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
  
  if (p.ncol() != m) {
    stop("Invalid number of columns of p");
  }
  
  
  for (iter=0; iter<iters; iter++) {
    wx = clone(w);  // reset
    g = clone(G);    
    r.fill(0);
    pmax = 0.0;
    pvalues = p.row(iter);
    
    for (step=0; step<m; step++) {
      
      // ratios of raw p-values divided by weights
      q.fill(0.0);
      for (i=0; i<m; i++) {
        if (wx(i) > 0.0) {
          q(i) = pvalues(i)/wx(i);
        }
      }
      
      q[q==0.0] = max(q) + 1.0;
      
      
      // identify the hypothesis to reject
      j = which_min(q);
      padj(iter,j) = std::max(std::min(q(j), 1.0), pmax);
      pmax = padj(iter, j);
      r(j) = 1;
      
      
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
            if ((r[k] == 0) && (l != k) && (g(l,j)*g(j,l) < 1.0 - 1.0e-12)) {
              g1(l,k) = (g(l,k) + g(l,j)*g(j,k))/(1.0 - g(l,j)*g(j,l));
            }
          }
        }
      }
      g = clone(g1);
      
    }
  }
  
  return padj;
}


//' @title Weight matrix for all intersection hypotheses
//' @description Obtains the weight matrix for all intersection hypotheses.
//'
//' @param w The vector of weights for elementary hypotheses.
//' @param G The transition matrix.
//'
//' @return The weight matrix starting with the global null hypothesis.
//'
//' @examples
//'
//' w = c(0.5,0.5,0,0)
//' g = matrix(c(0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0), nrow=4, ncol=4, byrow=TRUE)
//' (wgtmat = fwgtmat(w,g))
//'
//' @export
// [[Rcpp::export]]
NumericMatrix fwgtmat(const NumericVector& w,  
                      const NumericMatrix& G) {
  int m = w.size();
  int i, j, k, l;
  int ntests = (1 << m) - 1;
  
  if (is_true(any(w < 0.0))) {
    stop("w must be nonnegative");
  }
  
  if (sum(w) != 1.0) {
    stop("w must sum to 1");
  }
  
  if (G.nrow() != m || G.ncol() != m) {
    stop("Invalid dimension for G");
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
  
  NumericVector wx = clone(w);
  NumericMatrix g = clone(G);
  NumericMatrix wgtmat(ntests, m);
  NumericMatrix gtrmat((ntests+1)/2, m*m); // only need to store first half 
  for (i=0; i<ntests; i++) {
    int number = ntests - i;
    
    // binary representation of elementary hypotheses in the intersection
    IntegerVector cc(m);
    for (j=0; j<m; j++) {
      cc(j) = (number/(1 << (m - 1 - j))) % 2;
    }
    
    
    if (i >= 1) {
      j = which_min(cc);
      
      // indicators for hypotheses not in the super set
      IntegerVector cc1 = 1 - cc;
      cc1(j) = 0;
      
      // index of the super set
      int ip = 0;
      for (k=0; k<m; k++) {
        if (cc1(k)) {
          ip += (1 << (m - 1 - k));
        }
      }
      
      // load the weights from the super set
      for (k=0; k<m; k++) {
        wx(k) = wgtmat(ip, k);
      }
      
      // load the transition matrix from the super set
      for (k=0; k<m; k++) {
        for (l=0; l<m; l++) {
          g(k,l) = gtrmat(ip, k*m+l);
        }
      }
      
      // update the weights
      for (k=0; k<m; k++) {
        if (cc(k)) {
          wx(k) += wx(j)*g(j,k);
        }
      }
      wx(j) = 0;
      
      // update the transition matrix
      NumericMatrix g1(m,m);
      for (k=0; k<m; k++) {
        for (l=0; l<m; l++) {
          if (cc(k) && cc(l) && (k != l) && (g(k,j)*g(j,k) < 1.0 - 1.0e-12)) {
            g1(k,l) = (g(k,l) + g(k,j)*g(j,l))/(1 - g(k,j)*g(j,k));
          }
        }
      }
      g = g1;
      
    }
    
    // save the weights
    for (k=0; k<m; k++) {
      wgtmat(i,k) = wx(k);
    }
    
    // save the transition matrix
    if (i<(ntests+1)/2) {
      for (k=0; k<m; k++) {
        for (l=0; l<m; l++) {
          gtrmat(i, k*m+l) = g(k,l);
        }
      }
    }
  }
  
  return wgtmat;
}



// [[Rcpp::export]]
NumericMatrix fadjpsimcpp(const NumericMatrix& wgtmat,
                          const NumericMatrix& p,
                          const LogicalMatrix& family) {
  
  int ntests = wgtmat.nrow();
  int m = wgtmat.ncol();
  int niters = p.nrow();
  int nfams = family.nrow();
  int i, j, k, l, s, t, iter;
  LogicalMatrix incid(ntests, m);
  NumericMatrix pinter(niters, ntests);
  NumericMatrix padj(niters, m);
  
  if (family.ncol() != m) {
    stop("family must have as many individual hypotheses as columns");
  }
  
  for (j=0; j<m; j++) {
    if (sum(family(_,j)) != 1) {
      stop("Each hypothesis should belong to one or only one family");
    }
  }
  
  for (i=0; i<ntests; i++) {
    int number = ntests - i;
    
    // binary representation of elementary hypotheses in the intersection
    LogicalVector cc(m);
    for (j=0; j<m; j++) {
      cc(j) = (number/(1 << (m - 1 - j))) % 2;
    }
    
    // identify the active families and active hypotheses
    LogicalMatrix family0(nfams, m);
    for (k=0; k<nfams; k++) {
      for (j=0; j<m; j++) {
        family0(k,j) = family(k,j) && cc(j);
      }
    }
    
    int nhyps = sum(family0);
    IntegerVector nhyps0(nfams), hyp(nhyps), fam(nhyps);
    l = 0;
    for (k=0; k<nfams; k++) {
      for (j=0; j<m; j++) {
        if (family0(k,j)) {
          nhyps0(k)++;  // number of active hypotheses in family k
          fam(l) = k;   // family of the l-th active hypothesis
          hyp(l) = j;   // index of the l-th active hypothesis
          l++;
        }
      }
    }
    
    LogicalVector sub = (nhyps0 > 0);
    int nfamil1 = sum(sub); // number of active families
    IntegerVector nhyps1 = nhyps0[sub]; // # of active hypotheses by family
    
    NumericVector w(nhyps);
    for (j=0; j<nhyps; j++) {
      w(j) = wgtmat(i, hyp(j));
    }
    
    
    for (iter=0; iter<niters; iter++) {
      NumericVector pval(nhyps), cw(nhyps);
      for (j=0; j<nhyps; j++) {
        pval(j) = p(iter, hyp(j));
      }
      
      // sort p-values within each family and obtain associated cum weights
      s = 0;
      for (k=0; k<nfamil1; k++) {
        t = nhyps1(k);
        
        // extract p-values and weights in the family
        NumericVector p1(t), w1(t);
        for (j=0; j<t; j++) {
          p1(j) = pval(s+j);
          w1(j) = w(s+j);
        }
        
        
        // obtain the index of sorted p-values with the family
        IntegerVector index = seq_len(t) - 1;
        std::sort(index.begin(), index.end(),
                  [p1](const int&a, const int& b) {
                    return (p1(a) < p1(b));
                  });
        
        // replace original with sorted values
        for (j=0; j<t; j++) {
          pval(s+j) = p1(index(j));
          
          // obtain the cumulative weights within each family
          if (j==0) {
            cw(s+j) = w1(index(j));
          } else {
            cw(s+j) = cw(s+j-1) + w1(index(j));
          }
        }
        
        s += t;
      }
      
      double q = 1;
      for (j=0; j<nhyps; j++) {
        if (cw(j) > 0) {
          q = std::min(q, pval(j)/cw(j));
        }
      }
      
      pinter(iter,i) = q;
    }
    
    incid(i, _) = cc;
  }
  
  // obtain the adjusted p-values for individual hypotheses
  for (iter=0; iter<niters; iter++) {
    for (j=0; j<m; j++) {
      padj(iter,j) = 0;
      for (i=0; i<ntests; i++) {
        if (incid(i,j) && pinter(iter, i) > padj(iter,j)) {
          padj(iter,j) = pinter(iter, i);
        }
      }
    }
  }
  
  return padj;
}




//' @title Get efficacy boundaries for group sequential design
//' @description Obtains the efficacy stopping boundaries for a group
//' sequential design.
//'
//' @param k Look number for the current analysis.
//' @param informationRates Information rates up to the current look. Must be
//'   increasing and less than or equal to 1.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A vector of length \code{k} for the error spending  
//'   time at each analysis. Must be increasing and less than or equal to 1. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRates}.
//'
//' @return A numeric vector of critical values up to the current look.
//'
//' @examples
//'
//' getBound(k = 2, informationRates = c(0.5,1),
//'          alpha = 0.025, typeAlphaSpending = "sfOF")
//'
//' @export
// [[Rcpp::export]]
NumericVector getBound(
    const int k = NA_INTEGER,
    const NumericVector& informationRates = NA_REAL,
    const double alpha = 0.025,
    const String typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const NumericVector& userAlphaSpending = NA_REAL,
    const NumericVector& spendingTime = NA_REAL) {
  
  NumericVector informationRates1 = clone(informationRates);
  NumericVector spendingTime1 = clone(spendingTime);
  
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
  
  
  String asf = typeAlphaSpending;
  double asfpar = parameterAlphaSpending;
  
  if (asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < k) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegnative");
    } else if (k > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[k-1] > alpha) {
      stop("userAlphaSpending must not exceed the specified alpha");
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
  
  NumericVector theta(k); // mean values under H0, initialized to zero
  NumericVector t = clone(informationRates1); // information time for test stat
  NumericVector s = clone(spendingTime1); // spending time for alpha-spending
  NumericVector criticalValues(k);
  
  if (asf == "none") {
    for (int i=0; i<k-1; i++) {
      criticalValues[i] = 6.0;
    }
    criticalValues[k-1] = R::qnorm(1-alpha, 0, 1, 1, 0);
  } else if (asf == "OF" || asf == "P" || asf == "WT") {
    double Delta;
    if (asf == "OF") {
      Delta = 0;
    } else if (asf == "P") {
      Delta = 0.5;
    } else {
      Delta = asfpar;
    }
    
    auto f = [k, alpha, Delta, theta, t] (double aval)->double {
      NumericVector u(k), l(k);
      for (int i=0; i<k; i++) {
        u[i] = aval*pow(t[i], Delta-0.5);
        l[i] = -6.0;
      }
      
      List probs = exitprobcpp(u, l, theta, t);
      double cpu = sum(NumericVector(probs[0]));
      return cpu - alpha;
    };
    
    double cwt = brent(f, 0.0, 10.0, 1e-6);
    for (int i=0; i<k; i++) {
      criticalValues[i] = cwt*pow(t[i], Delta-0.5);
    }
  } else if (asf == "sfOF" || asf == "sfP" || asf == "sfKD" ||
    asf == "sfHSD" || asf == "user") {
    
    // stage 1
    double cumAlphaSpent;
    if (asf == "user") {
      cumAlphaSpent = userAlphaSpending[0];
    } else {
      cumAlphaSpent = errorSpentcpp(s[0], alpha, asf, asfpar);
    }
    
    criticalValues[0] = R::qnorm(1 - cumAlphaSpent, 0, 1, 1, 0);
    
    
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
      
      if (f(6.0) > 0) { // no alpha spent at current visit
        criticalValues[k1] = 6.0;
      } else {
        criticalValues[k1] = brent(f, -5.0, 6.0, 1.0e-6);
      }
    }
  } else {
    stop("Invalid type of alpha spending");
  }
  
  return criticalValues;
}


// [[Rcpp::export]]
NumericVector repeatedPValuecpp(
    const int kMax = NA_INTEGER,
    const String typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const double maxInformation = 1,
    const NumericMatrix& p = NA_REAL,
    const NumericMatrix& information = NA_REAL,
    const NumericMatrix& spendingTime = NA_REAL) {
  
  int iter, i, j, l, L;
  int B = p.nrow(), k = p.ncol();
  
  if (kMax <= 0) {
    stop("kMax must be a positive integer");
  }
  
  String asf = typeAlphaSpending;
  double asfpar = parameterAlphaSpending;
  
  if (!(asf=="OF" || asf=="P" || asf=="WT" || 
      asf=="sfOF" || asf=="sfP" || asf=="sfKD" ||
      asf=="sfHSD" || asf=="none")) {
    stop("Invalid type for typeAlphaSpending");
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
  
  
  if (maxInformation <= 0) {
    stop("maxInformation must be positive");
  }
  
  
  NumericMatrix info(B, k);
  
  if (information.ncol() != k) {
    stop("Invalid number of columns for information");
  } else if (information.nrow() != 1 && information.nrow() != B) {
    stop("Invalid number of rows for information");
  } else if (information.nrow() == 1 && B > 1) {
    for (iter=0; iter<B; iter++) {
      info(iter, _) = information(0, _);
    }
  } else {
    info = information;
  }
  
  for (iter=0; iter<B; iter++) {
    if (info(iter,0) <= 0) {
      stop("Elements of information must be positive");
    } else if (k>1 && is_true(any(diff(info(iter,_)) <= 0))) {
      stop("Elements of information must be increasing over time");
    }
  }
  
  
  NumericMatrix st(B, k);
  
  if (spendingTime.nrow()==1 && spendingTime.ncol()==1 
        && spendingTime(0,0)==0) {
    st.fill(NA_REAL);
  } else if (spendingTime.ncol() != k) {
    stop("Invalid number of columns for spendingTime");
  } else if (spendingTime.nrow() != 1 && spendingTime.nrow() != B) {
    stop("Invalid number of rows for spendingTime");
  } else if (spendingTime.nrow() == 1 && B > 1) {
    for (iter=0; iter<B; iter++) {
      st(iter, _) = spendingTime(0, _);
    }
  } else {
    st = spendingTime;
  }
  
  for (iter=0; iter<B; iter++) {
    if (is_false(all(is_na(st(iter,_))))) {
      if (st(iter,0) <= 0) {
        stop("Elements of spendingTime must be positive");
      } else if (k>1 && is_true(any(diff(st(iter,_)) <= 0))) {
        stop("Elements of spendingTime must be increasing over time");
      } else if (st(iter,k-1) > 1) {
        stop("spendingTime must be less than or equal to 1");
      }
    }
  }
  
  
  NumericMatrix repp(B, k);
  repp.fill(NA_REAL);

  for (iter=0; iter<B; iter++) {
    if (is_true(all(is_na(st(iter,_))))) { // use information rates
      LogicalVector qq = (info(iter,_) >= maxInformation);
      if (is_false(any(qq))) { // all observed info < maxinfo
        L = k-1;
      } else { // find index of first look with observed info >= maxinfo
        L = which_max(qq);
      } 
    } else { // use spending time
      L = k-1;
    }
    
    // information time for forming covariance matrix of test statistics
    NumericVector t1(L+1);
    for (l=0; l<=L; l++) {
      t1(l) = info(iter,l)/info(iter,L);
    }
    
    // spending time for error spending
    NumericVector s1(L+1);
    if (is_true(all(is_na(st(iter,_))))) { // use information rates
      for (l=0; l<=L; l++) {
        if (l == kMax-1 || info(iter,l) >= maxInformation) {
          s1(l) = 1.0;
        } else {
          s1(l) = info(iter,l)/maxInformation;
        }
      }
    } else { // using spending time
      for (l=0; l<=L; l++) {
        s1(l) = st(iter,l);
      }
    }
    
    
    for (i=0; i<=L; i++) {
      NumericVector t(i+1), s(i+1);
      for (j=0; j<=i; j++) {
        t(j) = t1(j);
        s(j) = s1(j);
      }
      
      double pvalue = p(iter,i);
      
      auto f = [i, t, asf, asfpar, s, pvalue](double aval)->double {
        NumericVector u = getBound(i+1, t, aval, asf, asfpar, 0, s);  
        return 1.0 - R::pnorm(u(i), 0, 1, 1, 0) - pvalue;
      };
      
      if (f(0.000001) > 0) {
        repp(iter,i) = 0.000001;
      } else if (f(0.999999) < 0) {
        repp(iter,i) = 0.999999;
      } else {
        repp(iter,i) = brent(f, 0.000001, 0.999999, 1.0e-6);
      }
    }
  }
  
  return repp;
}


// [[Rcpp::export]]
IntegerVector fseqboncpp(
    const NumericVector& w,
    const NumericMatrix& G,
    const double alpha = 0.025,
    const int kMax = NA_INTEGER,
    const StringVector& typeAlphaSpending = NA_STRING, 
    const NumericVector& parameterAlphaSpending = NA_REAL, 
    const LogicalMatrix& incidenceMatrix = NA_LOGICAL, 
    const NumericVector& maxInformation = NA_REAL, 
    const NumericMatrix& p = NA_REAL,
    const NumericMatrix& information = NA_REAL, 
    const NumericMatrix& spendingTime = NA_REAL) {
  
  // alias (shorter variable names)
  StringVector asf = typeAlphaSpending;
  NumericVector asfpar = clone(parameterAlphaSpending);
  LogicalMatrix incid = incidenceMatrix;
  NumericVector maxinfo = maxInformation;
  NumericMatrix rawp = clone(p);
  
  int m = w.size();
  
  if (incid.ncol() != kMax) {
    stop("Invalid number of columns for incidenceMatrix");
  };
  
  if (p.nrow() % m != 0) {
    stop("Invalid number of rows for p");
  }
  
  int k1 = p.ncol();
  
  if (k1 > kMax) {
    stop("Invalid number of columns for p");
  }
  
  
  int B = p.nrow()/m;
  
  int iter, i, j, k, l;
  IntegerVector reject(B*m);  // first look when the hypothesis is rejected

  NumericMatrix info(B*m, k1); // matrix of observed information
  NumericMatrix st(B*m, k1);  // matrix of spending time
  
  NumericVector wx(m);    // dynamic weights
  NumericMatrix g(m,m);   // dynamic transition matrix
  NumericMatrix g1(m,m);  // temporary transition matrix 
  
  
  if (is_true(any(w < 0.0))) {
    stop("w must be nonnegative");
  }
  
  if (sum(w) != 1.0) {
    stop("w must sum to 1");
  }
  
  if (G.nrow() != m || G.ncol() != m) {
    stop("Invalid dimension for G");
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
    stop("Invalid length for typeAlphaSpending");
  }
  
  if (is_true(any(is_na(asfpar))) && asfpar.size()==1) {
    asfpar = rep(NA_REAL, m);
  }
  
  if (asfpar.size() != m) {
    stop("Invalid length for parameterAlphaSpending");
  }
  
  for (i=0; i<m; i++) {
    if (!(asf(i)=="OF" || asf(i)=="P" || asf(i)=="WT" || 
        asf(i)=="sfOF" || asf(i)=="sfP" || asf(i)=="sfKD" || 
        asf(i)=="sfHSD" || asf(i)=="none")) {
      stop("Invalid type for typeAlphaSpending");
    }
    
    if (asf(i)=="WT" && R_isnancpp(asfpar(i))) {
      stop("Missing parameter for WT");
    }
    
    if (asf(i)=="sfKD" && R_isnancpp(asfpar(i))) {
      stop("Missing parameter for sfKD");
    }
    
    if (asf(i)=="sfHSD" && R_isnancpp(asfpar(i))) {
      stop("Missing parameter for sfHSD");
    }
    
    if (asf(i)=="sfKD" && asfpar(i) <= 0) {
      stop ("asfpar must be positive for sfKD");
    }
  }
  
  if (incid.nrow() != m) {
    stop("Invalid number of rows for incidenceMatrix");
  }
  
  if (maxinfo.size() != m) {
    stop("Invalid length for maxInformation");
  }
  
  if (is_true(any(maxinfo <= 0.0))) {
    stop("maxInformation must be positive");
  }
  
  
  if (information.ncol() != k1) {
    stop("information and p must have the same number of columns");
  } else if (information.nrow() != m && information.nrow() != B*m) {
    stop("Invalid number of rows for information");    
  } else if (information.nrow() == m && B > 1) {
    for (iter=0; iter<B; iter++) {
      for (j=0; j<m; j++) {
        info(iter*m+j, _) = information(j, _);  
      }
    }
  } else {
    info = information;
  }
  
  for (iter=0; iter<B*m; iter++) {
    if (info(iter,0) <= 0) {
      stop("Elements of information must be positive");
    } else if (k1>1 && is_true(any(diff(info(iter,_)) <= 0))) {
      stop("Elements of information must be increasing over time");
    }
  }
  
  
  
  if (spendingTime.nrow()==1 && spendingTime.ncol()==1 
        && spendingTime(0,0)==0) {
    st.fill(NA_REAL);
  } else if (spendingTime.ncol() != k1) {
    stop("spendingTime and p must have the same number of columns");
  } else if (spendingTime.nrow() != m && spendingTime.nrow() != B*m) {
    stop("Invalid number of rows for spendingTime");
  } else if (spendingTime.nrow() == m && B > 1) {
    for (iter=0; iter<B; iter++) {
      for (j=0; j<m; j++) {
        st(iter*m+j, _) = spendingTime(j, _);  
      }
    }
  } else {
    st = spendingTime;
  }
  
  
  
  // set information, p values, and spending time to missing at a study look
  // if the hypothesis is not to be tested at the study look
  for (j=0; j<m; j++) {
    for (k=0; k<k1; k++) {
      if (!incid(j,k)) {
        for (iter=0; iter<B; iter++) {
          info(iter*m+j,k) = NA_REAL;
          rawp(iter*m+j,k) = NA_REAL;
          st(iter*m+j,k) = NA_REAL;
        }
      }
    }
  }
  
  
  
  IntegerVector K0 = rowSums(incid); // maximum number of stages
  
  IntegerMatrix idx1(m, kMax); // study look
  IntegerMatrix idx2(m, kMax); // hypothesis look
  for (j=0; j<m; j++) {
    l = 0;
    for (k=0; k<kMax; k++) {
      if (incid(j,k)==1) {
        idx1(j,l) = k;
        idx2(j,k) = l;
        l++;
      } else {
        idx2(j,k) = NA_INTEGER;
      }
    }
    for (k=l; k<kMax; k++) {
      idx1(j,k) = NA_INTEGER;
    }
  }
  
  
  int step;
  double alphastar, asfpar1;
  String asf1;
  int K3, K4 = 0;
  
  NumericMatrix info1(m, k1), p1(m, k1), st1(m, k1), t1(m, k1), s1(m, k1);
  IntegerVector K1(m), K2(m), L(m);
  LogicalVector r(m);

  
  for (iter=0; iter<B; iter++) {
    wx = clone(w);  // reset
    g = clone(G); 
    r.fill(0);

    // extract iteration specific information, p-values, and spending time
    for (j=0; j<m; j++) {
      NumericVector irow = info.row(iter*m+j);
      NumericVector prow = rawp.row(iter*m+j);
      
      irow = irow[!is_na(irow)]; // exclude missing values
      prow = prow[!is_na(prow)];
      
      K1(j) = prow.size();
      K2(j) = idx1(j, K1(j)-1) + 1;  // last study look for hypothesis j
      
      if (irow.size() != prow.size()) {
        stop("information & p must have the same # of nonmissing elements");
      } else if (irow(0) <= 0) {
        stop("Elements of information must be positive");
      } else if (K1(j)>1 && is_true(any(diff(irow) <= 0))) {
        stop("Elements of information must be increasing over time");
      }
      
      info1(j, _) = irow;
      p1(j, _) = prow;
      
      NumericVector strow = st.row(iter*m+j);
      if (is_false(all(is_na(strow)))) {
        strow = strow[!is_na(strow)];
        if (strow.size() != prow.size()) {
          stop("spendingTime & p must have the same # of nonmissing elements");
        } else if (strow(0) <= 0) {
          stop("Elements of spendingTime must be positive");
        } else if (K1(j)>1 && is_true(any(diff(strow) <= 0))) {
          stop("Elements of spendingTime must be increasing over time");
        } else if (strow(K1(j)) > 1) {
          stop("spendingTime must be less than or equal to 1");
        }
        
        st1(j, _) = strow;
      } else {
        st1(j, _) = strow;
      }
      
      // index of the last look for each hypothesis
      if (is_true(all(is_na(st1(j, _))))) { // use information rates
        LogicalVector qq = (irow >= maxinfo(j));
        if (is_false(any(qq))) { // all observed info < maxinfo
          L(j) = K1(j) - 1;
        } else { // find index of first look with observed info >= maxinfo
          L(j) = which_max(qq);
        }
      } else { // use spending time
        L(j) = K1(j) - 1;
      }
      
      
      // information time for forming covariance matrix of test statistics
      for (l=0; l<=L(j); l++) {
        t1(j,l) = irow(l)/irow(L(j));
      }
      
      // spending time for error spending
      if (is_true(all(is_na(st1(j, _))))) { // use information rates
        for (l=0; l<=L(j); l++) {
          if (l == K0(j)-1 || irow(l) >= maxinfo(j)) {
            s1(j,l) = 1.0;
          } else {
            s1(j,l) = irow(l)/maxinfo(j);
          }
        }
      } else { // use spending time
        for (l=0; l<=L(j); l++) {
          s1(j,l) = st1(j,l);
        }
      }
    }
    
    K3 = max(K2);           // maximum look for the iteration 
    K4 = std::max(K3, K4);  // maximum look overall
    
    for (step=0; step<K3; step++) {  // loop over study look
      for (i=0; i<m; i++) {
        // find a hypothesis that can be rejected
        bool done = 1;
        for (j=0; j<m; j++) {
          if (incid(j, step) && wx(j) > 1.0e-4) {
            k = idx2(j, step);
            if (k <= L(j)) {
              NumericVector t(k+1);
              NumericVector s(k+1);
              for (l=0; l<=k; l++) {
                t(l) = t1(j,l);
                s(l) = s1(j,l);
              }
              
              asf1 = Rcpp::String(asf(j));
              asfpar1 = asfpar(j);
              
              double alp = wx(j)*alpha;
              NumericVector u = getBound(k+1, t, alp, asf1, asfpar1, 0, s);
              alphastar = 1.0 - R::pnorm(u(k), 0, 1, 1, 0);
              
              if (p1(j,k) < alphastar) {
                done = 0;
                break;
              }
            }
          }
        }
        
        
        if (!done) {
          // reject Hypothesis j
          r(j) = 1;
          reject(iter*m+j) = step+1;
          
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
        } else { // no more hypothesis to reject at this look
          break;
        }
      }
      
      // stop if all hypotheses have been rejected
      if (sum(r) == m) {
        break;
      }
    }
  }
  
  return reject;
}

