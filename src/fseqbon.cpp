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
//' @param informationRates Information rates up to the current look, must be
//'   increasing and less than or equal to 1.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A vector of length \code{k} for the error spending  
//'   time at each analysis, must be increasing and less than or equal to 1. 
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

//' @title Repeated p-values for group sequential design
//' @description Obtains the repeated p-values for a group sequential design.
//'
//' @param k Look number for the current analysis.
//' @param informationRates Information rates up to the current look, must be
//'   increasing and less than or equal to 1.
//' @param p The raw p-values at look 1 to look \code{k}. 
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @param spendingTime A vector of length \code{k} for the error spending  
//'   time at each analysis, must be increasing and less than or equal to 1. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRates}.
//'
//' @return The repeated p-values at look 1 to look \code{k}.
//'
//' @examples
//'
//' repeatedPValue(k = 3, informationRates = c(529, 700, 800)/800,
//'                p = c(0.2, 0.15, 0.1), typeAlphaSpending = "sfOF",
//'                spendingTime = c(0.6271186, 0.8305085, 1))
//'
//' @export
// [[Rcpp::export]]
NumericVector repeatedPValue(
    const int k = NA_INTEGER,
    const NumericVector& informationRates = NA_REAL,
    const NumericVector& p = NA_REAL,
    const String typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
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
  
  
  NumericVector repp(k);

  for (int i=0; i<k; i++) {
    NumericVector t(i+1), s(i+1);
    for (int j=0; j<=i; j++) {
      t(j) = informationRates1(j);
      s(j) = spendingTime1(j);
    }
    
    double pvalue = p(i);
    auto f = [i, t, asf, asfpar, s, pvalue](double aval)->double {
      NumericVector u = getBound(i+1, t, aval, asf, asfpar, 0, s);  
      return 1.0 - R::pnorm(u(i), 0, 1, 1, 0) - pvalue;
    };
    repp(i) = brent(f, 0.000001, 0.999999, 1.0e-6);
  }
  
  return repp;
}


//' @title Group sequential trials using Bonferroni-based graphical 
//' approaches
//' 
//' @description Obtains the test results for group sequential trials using 
//' graphical approaches based on weighted Bonferroni tests.
//'
//' @param w The vector of initial weights for elementary hypotheses.
//' @param G The initial transition matrix.
//' @param alpha The significance level. Defaults to 0.025.
//' @param typeAlphaSpending The vector of alpha spending functions. 
//'   Each element is one of the following: 
//'   "OF" for O'Brien-Fleming boundaries, 
//'   "P" for Pocock boundaries, "WT" for Wang & Tsiatis boundaries, 
//'   "sfOF" for O'Brien-Fleming type spending function, 
//'   "sfP" for Pocock type spending function,
//'   "sfKD" for Kim & DeMets spending function, 
//'   "sfHSD" for Hwang, Shi & DeCani spending function, 
//'   and "none" for no early efficacy stopping.
//' @param parameterAlphaSpending The vector of parameter values for the 
//'   alpha spending functions. Each element corresponds to the value of 
//'   Delta for "WT", rho for "sfKD", or gamma for "sfHSD".
//' @param incidenceMatrix The incidence matrix indicating whether the 
//'   specific hypothesis will be tested at the given look.
//' @param maxInformation The vector of targeted maximum information for each 
//'   hypothesis.
//' @param information The matrix of observed information for each hypothesis 
//'   by look.
//' @param p The matrix of raw p-values for each hypothesis by look. 
//' @param spendingTime The spending time for alpha spending. Defaults to 
//'   missing, in which case, it is the same as \code{informationRates} 
//'   calculated from \code{information} and \code{maxInformation}.
//' @param repeatedPValueFlag Flag for whether to report the repeated 
//'   p-values. Defaults to FALSE.
//' 
//' @return A list with a component vector to indicate the first look the 
//'   specific hypothesis rejected (0 if the hypothesis is not rejected) and,  
//'   if \code{reppfl = 1}, a component matrix for the repeated p-values for 
//'   each hypothesis over time. 
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
//'   typeAlphaSpending = rep("sfOF", 4),
//'   incidenceMatrix = matrix(1, nrow=4, ncol=3),
//'   maxInformation = rep(1,4),
//'   information = matrix(c(rep(1/3, 4), rep(2/3, 4)), nrow=4, ncol=2),
//'   p = matrix(c(0.0062, 0.017, 0.009, 0.13, 
//'                0.0002, 0.0035, 0.002, 0.06), 
//'              nrow=4, ncol=2),
//'   repeatedPValueFlag = 1)
//' 
//'
//' @export
// [[Rcpp::export]]
List fseqbon(const NumericVector& w, 
             const NumericMatrix& G, 
             const double alpha = 0.025,
             const StringVector& typeAlphaSpending = NA_STRING, 
             const NumericVector& parameterAlphaSpending = NA_REAL, 
             const LogicalMatrix& incidenceMatrix = NA_LOGICAL, 
             const NumericVector& maxInformation = NA_REAL, 
             const NumericMatrix& information = NA_REAL, 
             const NumericMatrix& p = NA_REAL,
             const NumericMatrix& spendingTime = NumericMatrix(0,1,1),
             const bool repeatedPValueFlag = 0) {
  
  // alias shorter variable names
  StringVector asf = typeAlphaSpending;
  NumericVector asfpar = parameterAlphaSpending;
  LogicalMatrix incid = incidenceMatrix;
  NumericVector maxinfo = maxInformation;
  NumericMatrix info = information;
  NumericMatrix sTime = spendingTime;
  bool reppfl = repeatedPValueFlag;
  
  int i, j, k, l, m = w.size(), K = incid.ncol(), B = info.nrow()/m;
  IntegerVector reject(B*m);  // first look when the hypothesis is rejected
  NumericMatrix repp(B*m, K);  // repeated p-values
  repp.fill(NA_REAL);
  
  NumericMatrix st(B*m, K);  // matrix of spending time
  
  NumericVector wx(m);    // dynamic weights
  NumericMatrix g(m,m);   // dynamic transition matrix
  NumericMatrix g1(m,m);  // temporary transition matrix 
  
  
  NumericVector asfpars = clone(asfpar);
  
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
    asfpars = rep(NA_REAL, m);
  }
  
  if (asfpars.size() != m) {
    stop("Invalid length for parameterAlphaSpending");
  }
  
  for (i=0; i<m; i++) {
    if (!(asf(i)=="OF" || asf(i)=="P" || asf(i)=="WT" || 
        asf(i)=="sfOF" || asf(i)=="sfP" || asf(i)=="sfKD" || 
        asf(i)=="sfHSD" || asf(i)=="none")) {
      stop("Invalid type for typeAlphaSpending");
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
    stop("Invalid number of rows for incidenceMatrix");
  }
  
  if (maxinfo.size() != m) {
    stop("Invalid length for maxInformation");
  }
  
  if (is_true(any(maxinfo <= 0.0))) {
    stop("maxInformation must be positive");
  }
  
  if (info.nrow() % m != 0) {
    stop("Invalid number of rows for information");
  }
  
  if (info.ncol() > K) {
    stop("Invalid number of columns for information");
  }
  
  if (p.nrow() % m != 0) {
    stop("Invalid number of rows for p");
  }
  
  if (p.ncol() > K) {
    stop("Invalid number of columns for p");
  }
  
  if (info.nrow() != p.nrow()) {
    stop("information and p must have the same number of rows");
  }
  
  if (info.ncol() != p.ncol()) {
    stop("information and p must have the same number of columns");
  }
  
  
  if (sTime.nrow()==1 && sTime.ncol()==1 && sTime(0,0)==0) {
    st.fill(NA_REAL);
  } else if (sTime.nrow() != p.nrow()) {
    stop("spendingTime and p must have the same number of rows");
  } else if (sTime.ncol() != p.ncol()) {
    stop("spendingTime and p must have the same number of columns");
  } else {
    st = clone(sTime);
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
  double alphastar, asfpar1;
  String asf1;
  int K3, K4 = 0;
  
  NumericMatrix info1(m, K), p1(m, K), st1(m, K);
  IntegerVector K1(m), K2(m);
  LogicalVector r(m);
  LogicalVector done(m);
  LogicalMatrix visited(m, K);
  
  for (int iter=0; iter<B; iter++) {
    wx = clone(w);  // reset
    g = clone(G); 
    r.fill(0);
    done.fill(0);
    visited.fill(0);
    
    // extract iteration specific information, p-values, and spending time
    for (j=0; j<m; j++) {
      NumericVector irow = info.row(iter*m+j);
      NumericVector prow = p.row(iter*m+j);
      
      irow = irow[!is_na(irow)]; // exclude missing values
      prow = prow[!is_na(prow)];
      
      if (irow.size() != prow.size()) {
        stop("information & p must have the same # of nonmissing elements");
      }
      
      K1(j) = irow.size();
      if (K1(j) > K0(j)) {
        K1(j) = K0(j); // truncated number of looks for hypothesis j 
      }
      
      K2(j) = idx1(j, K1(j)-1) + 1;  // last study look for hypothesis j
      
      info1(j, _) = irow;
      p1(j, _) = prow;
      
      NumericVector strow = st.row(iter*m+j);
      if (is_false(all(is_na(strow)))) {
        strow = strow[!is_na(strow)];
        if (strow.size() != prow.size()) {
          stop("spendingTime & p must have the same # of nonmissing elements");
        }
        st1(j, _) = strow;
      } else {
        st1(j, _) = strow;
      }
    }
    
    K3 = max(K2);           // maximum look for the iteration 
    K4 = std::max(K3, K4);  // maximum look overall 
    
    for (step=0; step<K3; step++) {  // loop over study look
      for (j=0; j<m; j++) {
        if (incid(j, step) && !done(j)) { // testable hypotheses
          kj = idx2(j, step);  // index of current look for hypothesis j
          NumericVector t(kj+1);  // information time
          NumericVector s(kj+1);  // spending time
          
          if ((reppfl && !visited(j, kj)) || (!r(j) && wx(j) > 1.0e-4)) {
            for (i=0; i<=kj; i++) {
              t(i) = info1(j, i)/info1(j, kj);
            }
            
            if (is_true(all(is_na(st1(j, _))))) { // use information rates
              if (kj < K0(j)-1 && info1(j, kj) < maxinfo(j)) {
                for (i=0; i<=kj; i++) {
                  s(i) = info1(j, i)/maxinfo(j);
                }
              } else {
                for (i=0; i<=kj; i++) {
                  s(i) = i < kj ? info1(j, i)/maxinfo(j) : 1.0;
                }
              }
            } else { // use spending time
              for (i=0; i<=kj; i++) {
                s(i) = st1(j, i);
              }
            }
          }
          
          
          if (reppfl && !visited(j, kj)) { // generate repeated p-values
            asf1 = Rcpp::String(asf(j));
            asfpar1 = asfpars(j);
            
            double pvalue = p1(j, kj);
            auto f = [kj, t, asf1, asfpar1, s, pvalue](double aval)->double {
              NumericVector u = getBound(kj+1, t, aval, asf1, asfpar1, 0, s);  
              return 1.0 - R::pnorm(u(kj), 0, 1, 1, 0) - pvalue;
            };
            
            double root;
            if (f(0.000001) > 0) {
              root = 0.000001;
            } else if (f(0.999999) < 0) {
              root = 0.999999;
            } else {
              root = brent(f, 0.000001, 0.999999, 1.0e-6);
            }
            repp(iter*m+j, step) = root;
            
            visited(j,kj) = 1;
          }
          
          
          if (!r(j) && wx(j) > 1.0e-4) {
            asf1 = Rcpp::String(asf(j));
            asfpar1 = asfpars(j);
            
            double alphaj = wx(j)*alpha;
            NumericVector u = getBound(kj+1, t, alphaj, asf1, asfpar1, 0, s);  
            alphastar = 1.0 - R::pnorm(u(kj), 0, 1, 1, 0);
            
            if (p1(j, kj) < alphastar) {
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
              
              j = -1;  // restart the loop over the index set 
            }
          }
        }
      }
      
      
      // done with a hypothesis if the max step or info has been reached
      for (j=0; j<m; j++) {
        if (incid(j, step) && !done(j)) {
          kj = idx2(j, step);  // index of current look for hypothesis j
          
          if (is_true(all(is_na(st1(j, _))))) {
            if (kj == K0(j)-1 || info1(j, kj) >= maxinfo(j)) {
              done(j) = 1; // current look is the last look
            }
          } else {
            if (kj == K0(j)-1) {
              done(j) = 1; // current look is the last look
            }
          }
        }
      }
      
      // stop if all hypotheses have been rejected
      if (sum(r) == m) {
        step = K3;
      }
    }
  }
  
  
  List result;
  
  if (reppfl) {
    result = List::create(_["rejectFirstLook"] = reject, 
                          _["repeatedPValues"] = repp(_, Range(0, K4-1)));
  } else {
    result = List::create(_["rejectFirstLook"] = reject);
  }
  
  return result;
}

