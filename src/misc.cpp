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
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
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
  
  LogicalVector r(m,1);
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
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
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
  
  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfpar = parameterAlphaSpending;
  
  if (!(asf=="of" || asf=="p" || asf=="wt" || asf=="sfof" || asf=="sfp" || 
      asf=="sfkd" || asf=="sfhsd" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }
  
  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }
  
  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
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
      LogicalVector x(i+1,1);
      
      double pvalue = p(iter,i);
      
      auto f = [i, t, asf, asfpar, s, x, pvalue](double aval)->double {
        NumericVector u = getBoundcpp(i+1, t, aval, asf, asfpar, 0, s, x);  
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
    std::string asfi = Rcpp::as<std::string>(asf(i));
    std::for_each(asfi.begin(), asfi.end(), [](char & c) {
      c = std::tolower(c);
    });
    
    if (!(asfi=="of" || asfi=="p" || asfi=="wt" || 
        asfi=="sfof" || asfi=="sfp" || asfi=="sfkd" || 
        asfi=="sfhsd" || asfi=="none")) {
      stop("Invalid value for typeAlphaSpending");
    }
    
    if ((asfi=="wt" || asfi=="sfkd" || asfi=="sfhsd") && 
        R_isnancpp(asfpar(i))) {
      stop("Missing value for parameterAlphaSpending");
    }
    
    if (asfi=="sfkd" && asfpar(i) <= 0) {
      stop ("parameterAlphaSpending must be positive for sfKD");
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
          stop("spendingTime & p must have equal # of nonmissing elements");
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
              
              LogicalVector x(k+1,1);
              
              asf1 = Rcpp::String(asf(j));
              asfpar1 = asfpar(j);
              
              double alp = wx(j)*alpha;
              NumericVector u = getBoundcpp(k+1, t, alp, asf1, asfpar1, 0, 
                                            s, x);
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


// [[Rcpp::export]]
NumericMatrix fstp2seqcpp(const NumericMatrix& p,
                          const NumericVector& gamma,
                          const String test = "hochberg",
                          const bool retest = 1) {
  
  std::string test1 = test;
  std::for_each(test1.begin(), test1.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  int nreps = p.nrow();
  int n = p.ncol();
  NumericMatrix padj(nreps, n);
  int m = n/2;
  
  int iter, i, j, s;
  
  NumericMatrix a(nreps, m);
  for (iter=0; iter<nreps; iter++) {
    for (j=0; j<m; j++) {
      double x1 = p(iter, 2*j);
      double x2 = p(iter, 2*j+1);
      
      a(iter,j) = 2*std::max(x1,x2)/(1+gamma[j]);
      if (test1=="holm") {
        a(iter,j) = std::max(a(iter,j), 2*std::min(x1,x2));
      }
    }
  }
  
  NumericMatrix d(m, m);
  for (s=0; s<m; s++) {
    double gmax = 0;
    for (j=s; j<m; j++) {
      if (j>s) gmax = std::max(gmax, gamma[j-1]);
      d(s,j) = std::max((1-gmax)/2, 1e-12);
    }
  }
  
  
  for (iter=0; iter<nreps; iter++) {
    NumericVector a1 = a(iter, _);
    for (i=0; i<m; i++) {
      double t1 = max(a1[Range(0,i)]);
      
      double t2 = 1;
      for (s=0; s<=i; s++) {
        double lhs = 0;
        
        if (s>0) {
          double y = max(a1[Range(0, s-1)]);
          lhs = std::max(lhs, y);
        }
        
        for (j=s; j<=i; j++) {
          lhs = std::max(lhs, p(iter,2*j)/d(s,j));
        }
        
        double rhs = 2*p(iter,2*s+1)/(1+gamma[s]);
        
        if (lhs < rhs) {
          t2 = std::min(t2, lhs);
        }
      }
      
      if (retest && m>1) {
        double t3 = 1;
        for (s=0; s<=std::min(i,m-2); s++) {
          double lhs = 0;
          
          if (s>0) {
            double y = max(a1[Range(0, s-1)]);
            lhs = std::max(lhs, y);
          }
          
          for (j=s; j<m; j++) {
            lhs = std::max(lhs, p(iter,2*j+1)/d(s,j));
          }
          
          for (j=s; j<=i; j++) {
            lhs = std::max(lhs, p(iter,2*j));
          }
          
          double rhs = 2*p(iter,2*s)/(1+gamma[s]);
          
          if (lhs < rhs) {
            t2 = std::min(t3, lhs);
          }
        }
        
        padj(iter,2*i) = std::min(t1, std::min(t2, t3));
      } else {
        padj(iter,2*i) = std::min(t1, t2);
      }
      
      
      t2 = 1;
      for (s=0; s<=i; s++) {
        double lhs = 0;
        
        if (s>0) {
          double y = max(a1[Range(0, s-1)]);
          lhs = std::max(lhs, y);
        }
        
        for (j=s; j<=i; j++) {
          lhs = std::max(lhs, p(iter,2*j+1)/d(s,j));
        }
        
        double rhs = 2*p(iter,2*s)/(1+gamma[s]);
        
        if (lhs < rhs) {
          t2 = std::min(t2, lhs);
        }
      }
      
      if (retest && m>1) {
        double t3 = 1;
        for (s=0; s<=std::min(i,m-2); s++) {
          double lhs = 0;
          
          if (s>0) {
            double y = max(a1[Range(0, s-1)]);
            lhs = std::max(lhs, y);                    
          }
          
          for (j=s; j<m; j++) {
            lhs = std::max(lhs, p(iter,2*j)/d(s,j));
          }
          
          for (j=s; j<=i; j++) {
            lhs = std::max(lhs, p(iter,2*j+1));
          }
          
          double rhs = 2*p(iter,2*s+1)/(1+gamma[s]);
          
          if (lhs < rhs) {
            t3 = std::min(t3, lhs);
          }
        }
        
        padj(iter, 2*i+1) = std::min(t1, std::min(t2, t3));
      } else {
        padj(iter, 2*i+1) = std::min(t1, t2);
      }
    }
  }
  
  return padj;
}


// Function to find the indices of all TRUE elements in a logical vector
IntegerVector which(LogicalVector vector) {
  IntegerVector true_indices;
  for (int i = 0; i < vector.size(); i++) {
    if (vector[i]) {
      true_indices.push_back(i);
    }
  }
  return true_indices;
}


// [[Rcpp::export]]
NumericMatrix fstdmixcpp(const NumericMatrix& p,
                         const LogicalMatrix& family,
                         const LogicalMatrix& serial,
                         const LogicalMatrix& parallel,
                         const NumericVector& gamma,
                         const String test = "hommel",
                         const bool exhaust = 1) {
  
  std::string test1 = test;
  std::for_each(test1.begin(), test1.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  // initialize various quantities
  int nreps = p.nrow();
  int m = p.ncol();
  int ntests = pow(2,m) - 1;
  int nfamily = family.nrow();
  IntegerVector nhyps = rowSums(family);
  
  // to store local p-values for the intersection tests
  NumericMatrix pinter(nreps, ntests);
  
  // incidence matrix for the elementary hypotheses
  NumericMatrix incid(ntests, m);
  
  
  for (int i=0; i<ntests; i++) {
    // expand to binary representations of the intersection hypothesis
    int number = ntests - i;
    LogicalVector cc(m);
    for (int j=0; j<m; j++) {
      cc[j] = (number/(int)std::pow(2, m-1-j)) % 2;
    }
    
    // indicator of active hyp in each family
    LogicalMatrix family0(nfamily, m);
    for (int j=0; j<nfamily; j++) {
      family0(j, _) = family(j, _) * cc;
    }
    
    // number of active hyp in each family
    IntegerVector nhyps0 = rowSums(family0);
    
    
    // determine restricted index set for each family
    LogicalVector cc1(m);
    for (int j=0; j<m; j++) {
      cc1[j] = 1;
    }
    
    for (int j=0; j<m; j++) {
      if (sum(serial(j,_))>0) {
        
        for (int k=0; k<m; k++) {
          if (serial(j,k) && cc[k]) cc1[j] = 0;
        }
        
        for (int k=0; k<m; k++) {
          if (serial(j,k) && !cc1[k]) cc1[j] = 0;
        }
      }
      
      if (sum(parallel(j,_))>0) {
        bool hit = 1;
        for (int k=0; k<m; k++) {
          if (parallel(j,k) && !cc[k]) hit = 0;
        }
        if (hit) cc1[j] = 0;
        
        hit = 1;
        for (int k=0; k<m; k++) {
          if (parallel(j,k) && cc1[k]) hit = 0;
        }
        if (hit) cc1[j] = 0;
      }
    }
    
    cc1 = cc1 * cc;
    
    
    // error rate function divided by alpha
    NumericVector errf(nfamily);
    for (int j=0; j<nfamily; j++) {
      errf[j] = nhyps0[j]>0 ? gamma[j] + (1-gamma[j])*nhyps0[j]/nhyps[j] : 0;
    }
    
    
    // allocated fraction of alpha for each family
    NumericVector coef(nfamily);
    coef[0] = 1;
    for (int j=1; j<nfamily; j++) {
      coef[j] = coef[j-1]*(1 - errf[j-1]);
    }
    
    int kmax = max(which(coef > 0));
    
    
    LogicalMatrix family1(kmax+1, m);
    for (int j=0; j<=kmax; j++) {
      family1(j, _) = family(j, _) * cc1;
    }
    
    // number of active hyp in each family
    IntegerVector nhyps1 = rowSums(family1);
    
    // index of active families
    IntegerVector sub = which(nhyps1 > 0);
    
    // active families;
    int nfamily2 = sub.size();
    
    // active families
    LogicalMatrix family2(nfamily2, m);
    for (int j=0; j<nfamily2; j++) {
      family2(j, _) = family1(sub[j], _);
    }
    
    // number of active hyp in active families
    IntegerVector nhyps2 = nhyps1[sub];
    
    // family indicators for active hypotheses
    IntegerVector fam, hyps2;
    for (int j=0; j<nfamily2; j++) {
      for (int k=0; k<m; k++) {
        if (family2(j,k)) {
          fam.push_back(j);
          hyps2.push_back(k);
        }
      }
    }
    
    // number of elementary hyp in the intersection
    int n = hyps2.size();
    
    
    
    // incidence matrix to ensure active hypotheses are clustered by family
    LogicalMatrix dup(nfamily2, n);
    for (int j=0; j<nfamily2; j++) {
      for (int k=0; k<n; k++) {
        if (fam[k] == j) {
          dup(j,k) = 1;
        } 
      }
    }
    
    
    // relative importance for active families in the intersection
    NumericVector c(n), coef1=coef[sub];
    for (int k=0; k<n; k++) {
      c[k] = 0;
      for (int j=0; j<nfamily2; j++) {
        c[k] += coef1[j] * dup(j,k);
      }
    }
    
    // weights for ordered p-values within each family
    NumericVector w(n);
    
    // truncation parameters for each family
    // use regular nonseparable procedure for last family
    NumericVector gam2 = gamma[sub];
    if (exhaust) gam2[nfamily2-1] = 1;
    
    // Bonferroni part of the weights
    NumericVector tbon(n);
    for (int j=0; j<nfamily2; j++) {
      coef1[j] = (1-gam2[j])/nhyps[sub[j]];
    }
    
    for (int k=0; k<n; k++) {
      tbon[k] = 0;
      for (int j=0; j<nfamily2; j++) {
        tbon[k] += coef1[j] * dup(j,k);
      }
    }
    
    
    // cumulative number of active hypotheses by family
    NumericVector ck(nfamily2+1);
    for (int j=1; j<=nfamily2; j++) {
      ck[j] = ck[j-1] + nhyps2[j-1];
    }
    
    // denominator weight for the ordered p-values in each family
    if (test1 == "hommel") {
      for (int k=0; k<n; k++) {
        // index of the hypothesis within a family
        int l = fam[k];
        int j = (k+1) - ck[l];
        w[k] = j * gam2[l]/nhyps2[l] + tbon[k];
      }
    } else if (test1 == "hochberg") {
      for (int k=0; k<n; k++) {
        int l = fam[k];
        int j = (k+1) - ck[l];
        w[k] = gam2[l]/(nhyps2[l] - j + 1) + tbon[k];
      }
    } else if (test1 == "holm") {
      for (int k=0; k<n; k++) {
        int l = fam[k];
        w[k] = gam2[l]/nhyps2[l] + tbon[k];
      }
    }
    
    
    for (int iter=0; iter<nreps; iter++) {
      // raw p-values
      NumericVector p1(n);
      for (int k=0; k<n; k++) {
        p1[k] = p(iter, hyps2[k]);
      }
      
      // order the p-values within each family 
      NumericVector p2(n);
      for (int j=0; j<nfamily2; j++) {
        Range indices = Range(ck[j], ck[j+1]-1);
        NumericVector p1s = p1[indices];
        p2[indices] = stl_sort(p1s);
      }
      
      NumericVector q = p2 / (w*c);
      double x = min(q);
      
      pinter(iter, i) = x;
    }
    
    incid(i, _) = cc;
  }
  
  // obtain the adjusted p-values for the elementary hypotheses
  NumericMatrix padj(nreps, m);
  for (int j=0; j<m; j++) {
    for (int iter=0; iter<nreps; iter++) {
      padj(iter,j) = 0;
      for (int i=0; i<ntests; i++) {
        if (incid(i,j)) {
          padj(iter,j) = std::max(padj(iter,j), pinter(iter,i));
        }
      }
      padj(iter,j) = std::min(padj(iter,j), 1.0);
    }
  }
  
  // modify the adjusted p-values to conform with the logical restrictions
  for (int j=0; j<m; j++) {
    if (sum(serial(j,_)) > 0) {
      for (int iter=0; iter<nreps; iter++) {
        double pre = 0;
        for (int k=0; k<m; k++) {
          if (serial(j,k)) {
            pre = std::max(pre, padj(iter,k));
          }
        }
        padj(iter,j) = std::max(padj(iter,j), pre);
      }
    }
    
    if (sum(parallel(j,_)) > 0) {
      for (int iter=0; iter<nreps; iter++) {
        double pre = 1;
        for (int k=0; k<m; k++) {
          if (parallel(j,k)) {
            pre = std::min(pre, padj(iter,k));
          }
        }
        padj(iter,j) = std::max(padj(iter,j), pre);
      }
    }
  }
  
  return(padj);
}


// [[Rcpp::export]]
NumericMatrix fmodmixcpp(const NumericMatrix& p,
                         const LogicalMatrix& family,
                         const LogicalMatrix& serial,
                         const LogicalMatrix& parallel,
                         const NumericVector& gamma,
                         const String test = "hommel",
                         const bool exhaust = 1) {
  
  std::string test1 = test;
  std::for_each(test1.begin(), test1.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  // initialize various quantities
  int nreps = p.nrow();
  int m = p.ncol();
  int ntests = pow(2,m) - 1;
  int nfamily = family.nrow();
  IntegerVector nhyps = rowSums(family);
  
  // to store local p-values for the intersection tests
  NumericMatrix pinter(nreps, ntests);
  
  // incidence matrix for the elementary hypotheses
  NumericMatrix incid(ntests, m);
  
  
  for (int i=0; i<ntests; i++) {
    // expand to binary representations of the intersection hypothesis
    int number = ntests - i;
    LogicalVector cc(m);
    for (int j=0; j<m; j++) {
      cc[j] = (number/(int)std::pow(2, m-1-j)) % 2;
    }
    
    // indicator of active hyp in each family
    LogicalMatrix family0(nfamily, m);
    for (int j=0; j<nfamily; j++) {
      family0(j, _) = family(j, _) * cc;
    }
    
    // number of active hyp in each family
    IntegerVector nhyps0 = rowSums(family0);
    
    
    // determine restricted index set for each family
    LogicalVector cc1(m);
    for (int j=0; j<m; j++) {
      cc1[j] = 1;
    }
    
    for (int j=0; j<m; j++) {
      if (sum(serial(j,_))>0) {
        
        for (int k=0; k<m; k++) {
          if (serial(j,k) && cc[k]) cc1[j] = 0;
        }
        
        for (int k=0; k<m; k++) {
          if (serial(j,k) && !cc1[k]) cc1[j] = 0;
        }
      }
      
      if (sum(parallel(j,_))>0) {
        bool hit = 1;
        for (int k=0; k<m; k++) {
          if (parallel(j,k) && !cc[k]) hit = 0;
        }
        if (hit) cc1[j] = 0;
        
        hit = 1;
        for (int k=0; k<m; k++) {
          if (parallel(j,k) && cc1[k]) hit = 0;
        }
        if (hit) cc1[j] = 0;
      }
    }
    
    LogicalVector cc2 = clone(cc1);  // denominator, Nstar
    cc1 = cc1 * cc;  // numerator, Istar
    
    
    // error rate function divided by alpha
    IntegerVector kstar(nfamily), nstar(nfamily);
    NumericVector errf(nfamily);
    for (int j=0; j<nfamily; j++) {
      kstar[j] = sum(family(j,_) * cc1);
      nstar[j] = sum(family(j,_) * cc2);
      errf[j] = kstar[j]>0 ? gamma[j] + (1-gamma[j])*kstar[j]/nstar[j] : 0;
    }
    
    
    // allocated fraction of alpha for each family
    NumericVector coef(nfamily);
    coef[0] = 1;
    for (int j=1; j<nfamily; j++) {
      coef[j] = coef[j-1]*(1 - errf[j-1]);
    }
    
    int kmax = max(which(coef > 0));
    
    
    LogicalMatrix family1(kmax+1, m);
    for (int j=0; j<=kmax; j++) {
      family1(j, _) = family(j, _) * cc1;
    }
    
    // number of active hyp in each family
    IntegerVector nhyps1 = rowSums(family1);
    
    // index of active families
    IntegerVector sub = which(nhyps1 > 0);
    
    // active families;
    int nfamily2 = sub.size();
    
    // active families
    LogicalMatrix family2(nfamily2, m);
    for (int j=0; j<nfamily2; j++) {
      family2(j, _) = family1(sub[j], _);
    }
    
    // number of active hyp in active families
    IntegerVector nhyps2 = nhyps1[sub];
    
    // family indicators for active hypotheses
    IntegerVector fam, hyps2;
    for (int j=0; j<nfamily2; j++) {
      for (int k=0; k<m; k++) {
        if (family2(j,k)) {
          fam.push_back(j);
          hyps2.push_back(k);
        }
      }
    }
    
    // number of elementary hyp in the intersection
    int n = hyps2.size();
    
    
    
    // incidence matrix to ensure active hypotheses are clustered by family
    LogicalMatrix dup(nfamily2, n);
    for (int j=0; j<nfamily2; j++) {
      for (int k=0; k<n; k++) {
        if (fam[k] == j) {
          dup(j,k) = 1;
        } 
      }
    }
    
    
    // relative importance for active families in the intersection
    NumericVector c(n), coef1=coef[sub];
    for (int k=0; k<n; k++) {
      c[k] = 0;
      for (int j=0; j<nfamily2; j++) {
        c[k] += coef1[j] * dup(j,k);
      }
    }
    
    // weights for ordered p-values within each family
    NumericVector w(n);
    
    // truncation parameters for each family
    // use regular nonseparable procedure for last family
    NumericVector gam2 = gamma[sub];
    if (exhaust) gam2[nfamily2-1] = 1;
    
    // Bonferroni part of the weights
    NumericVector tbon(n);
    for (int j=0; j<nfamily2; j++) {
      coef1[j] = (1-gam2[j])/nstar[j];
    }
    
    for (int k=0; k<n; k++) {
      tbon[k] = 0;
      for (int j=0; j<nfamily2; j++) {
        tbon[k] += coef1[j] * dup(j,k);
      }
    }
    
    
    // cumulative number of active hypotheses by family
    NumericVector ck(nfamily2+1);
    for (int j=1; j<=nfamily2; j++) {
      ck[j] = ck[j-1] + nhyps2[j-1];
    }
    
    // denominator weight for the ordered p-values in each family
    if (test1 == "hommel") {
      for (int k=0; k<n; k++) {
        // index of the hypothesis within a family
        int l = fam[k];
        int j = (k+1) - ck[l];
        w[k] = j * gam2[l]/nhyps2[l] + tbon[k];
      }
    } else if (test1 == "hochberg") {
      for (int k=0; k<n; k++) {
        int l = fam[k];
        int j = (k+1) - ck[l];
        w[k] = gam2[l]/(nhyps2[l] - j + 1) + tbon[k];
      }
    } else if (test1 == "holm") {
      for (int k=0; k<n; k++) {
        int l = fam[k];
        w[k] = gam2[l]/nhyps2[l] + tbon[k];
      }
    }
    
    
    for (int iter=0; iter<nreps; iter++) {
      // raw p-values
      NumericVector p1(n);
      for (int k=0; k<n; k++) {
        p1[k] = p(iter, hyps2[k]);
      }
      
      // order the p-values within each family 
      NumericVector p2(n);
      for (int j=0; j<nfamily2; j++) {
        Range indices = Range(ck[j], ck[j+1]-1);
        NumericVector p1s = p1[indices];
        p2[indices] = stl_sort(p1s);
      }
      
      NumericVector q = p2 / (w*c);
      double x = min(q);
      
      pinter(iter, i) = x;
    }
    
    incid(i, _) = cc;
  }
  
  // obtain the adjusted p-values for the elementary hypotheses
  NumericMatrix padj(nreps, m);
  for (int j=0; j<m; j++) {
    for (int iter=0; iter<nreps; iter++) {
      padj(iter,j) = 0;
      for (int i=0; i<ntests; i++) {
        if (incid(i,j)) {
          padj(iter,j) = std::max(padj(iter,j), pinter(iter,i));
        }
      }
      padj(iter,j) = std::min(padj(iter,j), 1.0);
    }
  }
  
  // modify the adjusted p-values to conform with the logical restrictions
  for (int j=0; j<m; j++) {
    if (sum(serial(j,_)) > 0) {
      for (int iter=0; iter<nreps; iter++) {
        double pre = 0;
        for (int k=0; k<m; k++) {
          if (serial(j,k)) {
            pre = std::max(pre, padj(iter,k));
          }
        }
        padj(iter,j) = std::max(padj(iter,j), pre);
      }
    }
    
    if (sum(parallel(j,_)) > 0) {
      for (int iter=0; iter<nreps; iter++) {
        double pre = 1;
        for (int k=0; k<m; k++) {
          if (parallel(j,k)) {
            pre = std::min(pre, padj(iter,k));
          }
        }
        padj(iter,j) = std::max(padj(iter,j), pre);
      }
    }
  }
  
  return(padj);
}


double f_pvalue(const double theta, 
                const NumericVector& b = NA_REAL, 
                const NumericVector& I = NA_REAL,
                const int L = NA_INTEGER, 
                const double zL = NA_REAL) {
  
  NumericVector upper(L), lower(L, -6.0), mu(L, theta), information(L);
  
  for (int l=0; l<L-1; l++) {
    upper[l] = b[l];
  }
  upper[L-1] = zL;
  
  for (int l=0; l<L; l++) {
    information[l] = I[l];
  }
  
  List probs = exitprobcpp(upper, lower, mu, information);
  
  return sum(NumericVector(probs[0]));
}


//' @title Confidence interval after trial termination
//' @description Obtains the p-value, median unbiased point estimate, and 
//' confidence interval after the end of a group sequential trial.
//'
//' @param b The upper boundaries on the Z-test statistic scale
//'   for efficacy stopping.
//' @param I The vector of cumulative information.
//' @param L The interim look.
//' @param zL The Z-test statistic at the interim look.
//'
//' @return A list with the following components: 
//' 
//' * \code{pvalue}: p-value for rejecting the null hypothesis.
//'  
//' * \code{thetahat}: Median unbiased point estimate of the parameter.
//' 
//' * \code{cilevel}: Confidence interval level.
//' 
//' * \code{lower}: Lower bound of confidence interval.
//' 
//' * \code{upper}: Upper bound of confidence interval.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' # group sequential design with 90% power to detect delta = 6
//' delta = 6
//' sigma = 17
//' n = 282
//' (des1 = getDesign(IMax = n/(4*sigma^2), theta = delta, kMax = 3, 
//'                   alpha = 0.05, typeAlphaSpending = "sfHSD", 
//'                   parameterAlphaSpending = -4))
//' 
//' # crossed the boundary at the second look
//' L = 2
//' n1 = n*2/3
//' delta1 = 7
//' sigma1 = 20
//' zL = delta1/sqrt(4/n1*sigma1^2)
//' 
//' # information based on estimated nuisance parameter
//' b = des1$byStageResults$efficacyBounds
//' t = des1$byStageResults$informationRates
//' I = n*t/(4*sigma1^2)
//' 
//' # p-value, point estimate, and confidence interval
//' getCI(b, I, L, zL) 
//' 
//' @export
// [[Rcpp::export]]
List getCI(const NumericVector& b = NA_REAL, 
           const NumericVector& I = NA_REAL,
           const int L = NA_INTEGER, 
           const double zL = NA_REAL) {
  
  if (is_true(any(is_na(b)))) {
    stop("b must be provided");
  }
  
  if (is_true(any(is_na(I)))) {
    stop("I must be provided");
  }
  
  if (b.size() < L) {
    stop("Invalid length for b");
  }
  
  if (I.size() < L) {
    stop("Invalid length for I");
  } else if (I[0] <= 0) {
    stop("Elements of I must be positive");
  } else if (I.size() > 1 && is_true(any(diff(I) <= 0))) {
    stop("Elements of I must be increasing");
  }
  
  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (L < 1) {
    stop("L must be a positive integer");
  }
  
  double pvalue = f_pvalue(0, b, I, L, zL);
  
  int kMax = b.size();
  NumericVector a(kMax, -6.0), zero(kMax); 
  List probs = exitprobcpp(b, a, zero, I);
  double cilevel = 1 - 2*sum(NumericVector(probs[0]));
  
  NumericVector interval(2);
  interval[0] = (zL - 6)/sqrt(I[L-1]);
  interval[1] = (zL + 6)/sqrt(I[L-1]);
  double tol = 0.0001;
  
  auto f = [b, I, L, zL](double theta)->double {
    return f_pvalue(theta, b, I, L, zL) - 0.5;
  };
  double thetahat = brent(f, interval[0], interval[1], tol);
  
  auto f1 = [b, I, L, zL, cilevel](double theta)->double {
    return f_pvalue(theta, b, I, L, zL) - (1-cilevel)/2;
  };
  double lower = brent(f1, interval[0], thetahat, tol);
  
  auto f2 = [b, I, L, zL, cilevel](double theta)->double {
    return f_pvalue(theta, b, I, L, zL) - (1+cilevel)/2;
  };
  double upper = brent(f2, thetahat, interval[1], tol);
  
  List result = List::create(
    _["pvalue"] = pvalue,
    _["thetahat"] = thetahat,
    _["cilevel"] = cilevel,
    _["lower"] = lower,
    _["upper"] = upper);
  
  return result;
}


//' @title Repeated confidence interval for group sequential design
//' @description Obtains the repeated confidence interval 
//' for a group sequential trial.
//'
//' @param IMax The maximum information.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @param spendingTime A vector of length \code{kMax} for the error spending 
//'   time at each analysis. Defaults to missing, in which case, it is the 
//'   same as \code{informationRates}.
//' @param L The interim look.
//' @param zL The Z-test statistic at the interim look.
//'
//' @return A list with the following components: 
//' 
//' * \code{pvalue}: Repeated p-value for rejecting the null hypothesis.
//'  
//' * \code{thetahat}: Point estimate of the parameter.
//' 
//' * \code{cilevel}: Confidence interval level.
//' 
//' * \code{lower}: Lower bound of repeated confidence interval.
//' 
//' * \code{upper}: Upper bound of repeated confidence interval.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' # group sequential design with 90% power to detect delta = 6
//' delta = 6
//' sigma = 17
//' n = 282
//' (des1 = getDesign(IMax = n/(4*sigma^2), theta = delta, kMax = 3, 
//'                   alpha = 0.05, typeAlphaSpending = "sfHSD", 
//'                   parameterAlphaSpending = -4))
//' 
//' # results at the second look
//' L = 2
//' n1 = n*2/3
//' delta1 = 7
//' sigma1 = 20
//' zL = delta1/sqrt(4/n1*sigma1^2)
//' 
//' # information based on estimated nuisance parameter
//' b = des1$byStageResults$efficacyBounds
//' t = des1$byStageResults$informationRates
//' I = n*t/(4*sigma1^2)
//' 
//' # repeated confidence interval
//' getRCI(IMax = n/(4*sigma1^2), kMax = 3, alpha = 0.05, 
//'        typeAlphaSpending = "sfHSD", parameterAlphaSpending = -4, 
//'        L = L, zL = zL)
//' 
//' @export
// [[Rcpp::export]]
List getRCI(const double IMax = NA_REAL,
            const int kMax = NA_REAL,
            const NumericVector& informationRates = NA_REAL,
            const double alpha = 0.025,
            const String typeAlphaSpending = "sfOF",
            const double parameterAlphaSpending = NA_REAL,
            const NumericVector spendingTime = NA_REAL,
            const int L = NA_INTEGER, 
            const double zL = NA_REAL) {
  
  NumericVector t = clone(informationRates);
  NumericVector s = clone(spendingTime);
  
  if (R_isnancpp(IMax)) {
    stop("IMax must be provided");
  }
  
  if (IMax <= 0) {
    stop("IMax must be positive");
  }
  
  if (R_isnancpp(kMax)) {
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
    t = as<NumericVector>(tem)/(kMax+0.0);
  }
  
  if (R_isnancpp(alpha)) {
    stop("alpha must be provided");
  }
  
  if (alpha < 0.00001 || alpha >= 0.5) {
    stop("alpha must lie in [0.00001, 0.5)");
  }
  
  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfpar = parameterAlphaSpending;
  
  if (!(asf=="of" || asf=="p" || asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }
  
  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }
  
  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
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
    s = clone(t);
  }
  
  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (L < 1) {
    stop("L must be a positive integer");
  }
  
  if (kMax < L) {
    stop("kMax must be greater than or equal to L");
  }
  
  
  NumericVector I = IMax*t;
  LogicalVector x(kMax, 1);
  IntegerVector i = Range(0, L-1);
  NumericVector b = getBoundcpp(L, t[i], alpha, asf, asfpar, 0, s[i], x[i]);
  
  // repeated confidence interval
  double lower = (zL - b[L-1])/sqrt(I[L-1]);
  double upper = (zL + b[L-1])/sqrt(I[L-1]);
  
  // point estimate is the lower bound for alpha = 0.5
  NumericVector u = getBoundcpp(L, t[i], 0.5, asf, asfpar, 0, s[i], x[i]);
  double thetahat = (zL - u[L-1])/sqrt(I[L-1]);
  
  // repeated p-value is alpha for which the lower bound of theta is zero
  auto f = [L, zL, t, asf, asfpar, s, x, i](double aval)->double {
    NumericVector u = getBoundcpp(L, t[i], aval, asf, asfpar, 0, s[i], x[i]);  
    return zL - u[L-1];
  };
  
  double pvalue;
  if (f(0.000001) > 0) {
    pvalue = 0.000001;
  } else if (f(0.999999) < 0) {
    pvalue = 0.999999;
  } else {
    pvalue = brent(f, 0.000001, 0.999999, 1.0e-6);
  }
  
  List result = List::create(
    _["pvalue"] = pvalue,
    _["thetahat"] = thetahat,
    _["cilevel"] = 1-2*alpha,
    _["lower"] = lower,
    _["upper"] = upper);
  
  return result;
}


double f_astar(const double theta, 
               const NumericVector& b2, 
               const NumericVector& I2, 
               const int L2, 
               const double zL2) {
  
  NumericVector upper(L2), lower(L2, -6.0), mu(L2, theta), information(L2);
  
  for (int l=0; l<L2-1; l++) {
    upper[l] = b2[l];
  }
  upper[L2-1] = zL2;
  
  for (int l=0; l<L2; l++) {
    information[l] = I2[l];
  }
  
  List probs = exitprobcpp(upper, lower, mu, information);
  return sum(NumericVector(probs[0]));
}


List f_bwimage(const double theta, 
               const int kMax, 
               const NumericVector& b, 
               const NumericVector& I, 
               const int L, 
               const double zL, 
               const NumericVector& b2, 
               const NumericVector& I2, 
               const int L2, 
               const double zL2) {
  
  double astar = f_astar(theta, b2, I2, L2, zL2);
  int k1 = kMax - L;
  
  NumericVector b1(k1), a1(k1, -6.0), mu(k1, theta), I1(k1);
  for (int l=0; l<k1; l++) {
    b1[l] = (b[l+L] - sqrt(I[L-1]/I[l+L])*zL)/sqrt(1 - I[L-1]/I[l+L]);
    I1[l] = I[l+L] - I[L-1];
  }
  
  List probs = exitprobcpp(b1, a1, mu, I1);
  NumericVector pu = NumericVector(probs[0]);
  
  NumericVector p(k1+1);
  p[0] = 0;
  for (int l=0; l<k1; l++) {
    p[l+1] = p[l] + pu[l];
  }
  
  NumericVector astars(1, astar);
  IntegerVector js = findInterval2(astars, p);
  int j = js[0];
  
  double z1j;
  if (j==1) {
    z1j = R::qnorm(1 - astar, 0, 1, 1, 0);
  } else {
    auto f = [j, b1, I1, theta, astar](double z)->double {
      NumericVector upper(j), lower(j, -6.0), mu(j, theta), information(j);
      
      for (int l=0; l<j-1; l++) {
        upper[l] = b1[l];
      }
      upper[j-1] = z;
      
      for (int l=0; l<j; l++) {
        information[l] = I1[l];
      }
      
      List probs = exitprobcpp(upper, lower, mu, information);
      
      return sum(NumericVector(probs[0])) - astar;
    };
    
    z1j = brent(f, -6, 6, 0.0001);
  }
  
  int J = L+j;
  double zJ = sqrt(I[L-1]/I[J-1])*zL + sqrt(1-I[L-1]/I[J-1])*z1j;
  
  List result = List::create(
    _["J"] = J,
    _["zJ"] = zJ);
  
  return result;
}


double f_bwpvalue(const double theta, 
                  const int kMax = NA_INTEGER, 
                  const NumericVector& b = NA_REAL, 
                  const NumericVector& I = NA_REAL,
                  const int L = NA_INTEGER, 
                  const double zL = NA_REAL, 
                  const NumericVector& b2 = NA_REAL, 
                  const NumericVector& I2 = NA_REAL,
                  const int L2 = NA_INTEGER, 
                  const double zL2 = NA_REAL) {
  
  List a = f_bwimage(theta, kMax, b, I, L, zL, b2, I2, L2, zL2);
  
  int J = a[0];
  double zJ = a[1];
  
  NumericVector upper(J), lower(J, -6.0), mu(J, theta), information(J);
  
  for (int l=0; l<J-1; l++) {
    upper[l] = b[l];
  }
  upper[J-1] = zJ;
  
  for (int l=0; l<J; l++) {
    information[l] = I[l];
  }
  
  List probs = exitprobcpp(upper, lower, mu, information);
  
  return sum(NumericVector(probs[0]));
}


//' @title Confidence interval after adaptation
//' @description Obtains the p-value, median unbiased point estimate, and 
//' confidence interval using the backward image method after the end of 
//' an adaptive trial.
//'
//' @param kMax The maximum number of stages for the primary trial.
//' @param b The upper boundaries on the Z-test statistic scale
//'   for efficacy stopping for the primary trial.
//' @param I The vector of cumulative information for the primary trial.
//' @param L The interim look of the primary trial.
//' @param zL The Z-test statistic at the interim look of the primary trial.
//' @param b2 The upper boundaries on the Z-test statistic scale
//'   for efficacy stopping for the secondary trial.
//' @param I2 The vector of cumulative information for the secondary trial.
//' @param L2 The interim look of the secondary trial.
//' @param zL2 The Z-test statistic at the interim look of the 
//'   secondary trial.
//'
//' @return A list with the following components: 
//' 
//' * \code{pvalue}: p-value for rejecting the null hypothesis.
//'  
//' * \code{thetahat}: Median unbiased point estimate of the parameter.
//' 
//' * \code{cilevel}: Confidence interval level.
//' 
//' * \code{lower}: Lower bound of confidence interval.
//' 
//' * \code{upper}: Upper bound of confidence interval.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{adaptDesign}}
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
//' kMax = des1$overallResults$kMax
//' b = des1$byStageResults$efficacyBounds
//' t = des1$byStageResults$informationRates
//' I = n*t/(4*sigma1^2)  # information based on estimated nuisance parameter
//' 
//' # Muller & Schafer (2001) method to design the secondary trial: 
//' # 3-look gamma(-2) spending with 84% power at delta = 4.5, sigma = 17
//' n2 = 300
//' (des2 = adaptDesign(
//'   beta = NA, INew = n2/(4*sigma^2), L, zL, theta = delta1,
//'   kMax = des1$overallResults$kMax, 
//'   informationRates = des1$byStageResults$informationRates,
//'   criticalValues = des1$byStageResults$efficacyBounds,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpending = "sfHSD", 
//'   parameterAlphaSpending = -2))
//'   
//' # termination at the second look of the secondary trial
//' L2 = 2
//' theta2 = 6.6
//' sigma2 = 19.5
//' zL2 = theta2/sqrt(4*sigma2^2/200)
//'  
//' b2 = des2$secondaryTrial$byStageResults$efficacyBounds
//' t2 = des2$secondaryTrial$byStageResults$informationRates
//' I2 = n2*t2/(4*sigma2^2)
//' 
//' # p-value, point estimate, and confidence interval
//' getADCI(kMax, b, I, L, zL, b2, I2, L2, zL2)
//' 
//' @export
// [[Rcpp::export]]
List getADCI(const int kMax = NA_INTEGER, 
             const NumericVector& b = NA_REAL, 
             const NumericVector& I = NA_REAL,
             const int L = NA_INTEGER, 
             const double zL = NA_REAL, 
             const NumericVector& b2 = NA_REAL, 
             const NumericVector& I2 = NA_REAL,
             const int L2 = NA_INTEGER, 
             const double zL2 = NA_REAL) {
  
  if (R_isnancpp(kMax)) {
    stop("kMax must be provided");
  }
  
  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }
  
  if (is_true(any(is_na(b)))) {
    stop("b must be provided");
  }
  
  if (is_true(any(is_na(I)))) {
    stop("I must be provided");
  }
  
  if (I.size() != kMax) {
    stop("Invalid length for I");
  } else if (I[0] <= 0) {
    stop("Elements of I must be positive");
  } else if (I.size() > 1 && is_true(any(diff(I) <= 0))) {
    stop("Elements of I must be increasing");
  }
  
  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (L < 1) {
    stop("L must be a positive integer");
  }
  
  if (kMax <= L) {
    stop("kMax must be greater than L");
  }
  
  if (is_true(any(is_na(b2)))) {
    stop("b2 must be provided");
  }
  
  if (is_true(any(is_na(I2)))) {
    stop("I2 must be provided");
  }
  
  if (I2[0] <= 0) {
    stop("Elements of I2 must be positive");
  } else if (I2.size() > 1 && is_true(any(diff(I2) <= 0))) {
    stop("Elements of I2 must be increasing");
  }
  
  if (R_isnancpp(L2)) {
    stop("L2 must be provided");
  }
  
  if (R_isnancpp(zL2)) {
    stop("zL2 must be provided");
  }
  
  if (L2 < 1) {
    stop("L2 must be a positive integer");
  }
  
  NumericVector a(kMax, -6.0), zero(kMax); 
  List probs = exitprobcpp(b, a, zero, I);
  double cilevel = 1 - 2*sum(NumericVector(probs[0]));
  
  int K = kMax;
  double pvalue = f_bwpvalue(0, K, b, I, L, zL, b2, I2, L2, zL2);
  
  NumericVector interval(2);
  interval[0] = (zL - b[L-1])/sqrt(I[L-1]);
  interval[1] = (zL + b[L-1])/sqrt(I[L-1]);
  double tol = 0.0001;
  
  auto f = [K, b, I, L, zL, b2, I2, L2, zL2](double theta)->double {
    return f_bwpvalue(theta, K, b, I, L, zL, b2, I2, L2, zL2) - 0.5;
  };
  double thetahat = brent(f, interval[0], interval[1], tol);
  
  auto f1 = [K, b, I, L, zL, b2, I2, L2, zL2, cilevel](double theta)->double {
    return f_bwpvalue(theta, K, b, I, L, zL, b2, I2, L2, zL2) - (1-cilevel)/2;
  };
  double lower = brent(f1, interval[0], thetahat, tol);
  
  auto f2 = [K, b, I, L, zL, b2, I2, L2, zL2, cilevel](double theta)->double {
    return f_bwpvalue(theta, K, b, I, L, zL, b2, I2, L2, zL2) - (1+cilevel)/2;
  };
  double upper = brent(f2, thetahat, interval[1], tol);
  
  List result = List::create(
    _["pvalue"] = pvalue,
    _["thetahat"] = thetahat,
    _["cilevel"] = cilevel,
    _["lower"] = lower,
    _["upper"] = upper);
  
  return result;
}


//' @title Repeated confidence interval after adaptation
//' @description Obtains the repeated p-value, point estimate, and 
//' repeated confidence interval for an adaptive group sequential trial.
//' 
//' @param IMax The maximum information of the primary trial.
//' @param kMax The maximum number of stages for the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param alpha The significance level for the primary trial. 
//'   Defaults to 0.025.
//' @param typeAlphaSpending The type of alpha spending for the primary trial. 
//'   One of the following: 
//'   "OF" for O'Brien-Fleming boundaries, 
//'   "P" for Pocock boundaries, 
//'   "WT" for Wang & Tsiatis boundaries, 
//'   "sfOF" for O'Brien-Fleming type spending function, 
//'   "sfP" for Pocock type spending function, 
//'   "sfKD" for Kim & DeMets spending function, 
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and 
//'   "none" for no early efficacy stopping. 
//'   Defaults to "sfOF".
//' @param parameterAlphaSpending The parameter value for the alpha 
//'   spending for the primary trial. 
//'   Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".
//' @param spendingTime A vector of length \code{kMax} for the error 
//'   spending time at each analysis of the primary trial. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRates}.
//' @param L The interim look of the primary trial.
//' @param zL The Z-test statistic at the interim look of the primary trial.
//' @param INew The maximum information for the secondary trial.
//' @param L2 The interim look of the secondary trial.
//' @param zL2 The Z-test statistic at the interim look of the 
//'   secondary trial.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method 
//'   for trial adaptation.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
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
//' @param parameterAlphaSpendingNew The parameter value for the alpha 
//'   spending for the secondary trial. 
//'   Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".
//' @param spendingTimeNew A vector of length \code{kNew} for the error 
//'   spending time at each analysis of the secondary trial. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRatesNew}.
//'
//' @return A list with the following components: 
//' 
//' * \code{pvalue}: Repeated p-value for rejecting the null hypothesis.
//'  
//' * \code{thetahat}: Point estimate of the parameter.
//' 
//' * \code{cilevel}: Confidence interval level.
//' 
//' * \code{lower}: Lower bound of repeated confidence interval.
//' 
//' * \code{upper}: Upper bound of repeated confidence interval.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{adaptDesign}}
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
//' kMax = des1$overallResults$kMax
//' b = des1$byStageResults$efficacyBounds
//' t = des1$byStageResults$informationRates
//' I = n*t/(4*sigma1^2)  # information based on estimated nuisance parameter
//' 
//' # Muller & Schafer (2001) method to design the secondary trial: 
//' # 3-look gamma(-2) spending with 84% power at delta = 4.5, sigma = 17
//' n2 = 300
//' (des2 = adaptDesign(
//'   beta = NA, INew = n2/(4*sigma^2), L, zL, theta = delta1,
//'   kMax = des1$overallResults$kMax, 
//'   informationRates = des1$byStageResults$informationRates,
//'   criticalValues = des1$byStageResults$efficacyBounds,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpending = "sfHSD", 
//'   parameterAlphaSpending = -2))
//'   
//' # termination at the second look of the secondary trial
//' L2 = 2
//' theta2 = 6.6
//' sigma2 = 19.5
//' zL2 = theta2/sqrt(4*sigma2^2/200)
//'  
//' b2 = des2$secondaryTrial$byStageResults$efficacyBounds
//' t2 = des2$secondaryTrial$byStageResults$informationRates
//' I2 = n2*t2/(4*sigma2^2)
//' 
//' # repeated confidence interval
//' getADRCI(IMax = n/(4*sigma1^2), kMax = 3, 
//'          informationRates = t,
//'          alpha = 0.05,  typeAlphaSpending = "sfHSD", 
//'          parameterAlphaSpending = -4,
//'          L = L, zL = zL, 
//'          INew = n2/(4*sigma2^2), 
//'          L2 = L2, zL2 = zL2, 
//'          MullerSchafer = TRUE, 
//'          kNew = 3, informationRatesNew = t2, 
//'          typeAlphaSpendingNew = "sfHSD", 
//'          parameterAlphaSpendingNew = -2)
//' 
//' @export
// [[Rcpp::export]]
List getADRCI(const double IMax = NA_REAL,
              const int kMax = NA_INTEGER,
              const NumericVector& informationRates = NA_REAL, 
              const double alpha = 0.25,
              const String typeAlphaSpending = "sfOF",
              const double parameterAlphaSpending = NA_REAL, 
              const NumericVector& spendingTime = NA_REAL,
              const int L = NA_INTEGER, 
              const double zL = NA_REAL, 
              const double INew = NA_REAL, 
              const int L2 = NA_INTEGER, 
              const double zL2 = NA_REAL, 
              const bool MullerSchafer = 0,
              const int kNew = NA_INTEGER, 
              const NumericVector& informationRatesNew = NA_REAL, 
              const String typeAlphaSpendingNew = "sfOF",
              const double parameterAlphaSpendingNew = NA_REAL, 
              const NumericVector& spendingTimeNew = NA_REAL) {
  
  NumericVector t = clone(informationRates);
  NumericVector t2 = clone(informationRatesNew);
  NumericVector s = clone(spendingTime);
  NumericVector s2 = clone(spendingTimeNew);
  
  if (R_isnancpp(IMax)) {
    stop("IMax must be provided");
  }
  
  if (IMax <= 0) {
    stop("IMax must be positive");
  }
  
  if (R_isnancpp(kMax)) {
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
    t = as<NumericVector>(tem)/(kMax+0.0);
  }
  
  if (R_isnancpp(alpha)) {
    stop("alpha must be provided");
  }
  
  if (alpha < 0.00001 || alpha >= 0.5) {
    stop("alpha must lie in [0.00001, 0.5)");
  }
  
  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfpar = parameterAlphaSpending;
  
  if (!(asf=="of" || asf=="p" || asf=="wt" || asf=="sfof" || asf=="sfp" ||
      asf=="sfkd" || asf=="sfhsd" || asf=="none")) {
    stop("Invalid value for typeAlphaSpending");
  }
  
  if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
    stop("Missing value for parameterAlphaSpending");
  }
  
  if (asf=="sfkd" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
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
    s = clone(t);
  }
  
  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (L < 1) {
    stop("L must be a positive integer");
  }
  
  if (kMax <= L) {
    stop("kMax must be greater than L");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (R_isnancpp(INew)) {
    stop("INew must be provided");
  }
  
  if (INew <= 0) {
    stop("INew must be positive");
  }
  
  if (R_isnancpp(L2)) {
    stop("L2 must be provided");
  }
  
  if (L2 < 1) {
    stop("L2 must be a positive integer");
  }
  
  if (R_isnancpp(zL2)) {
    stop("zL2 must be provided");
  }
  
  std::string asf2 = typeAlphaSpendingNew;
  std::for_each(asf2.begin(), asf2.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfpar2 = parameterAlphaSpendingNew;
  
  if (MullerSchafer) {
    if (R_isnancpp(kNew)) {
      stop("kNew must be provided");
    } else if (kNew < 1) {
      stop("kNew must be a positive integer");
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
      t2 = as<NumericVector>(tem)/(kNew+0.0);
    }
    
    if (!(asf2=="of" || asf2=="p" || asf2=="wt" || 
        asf2=="sfof" || asf2=="sfp" ||
        asf2=="sfkd" || asf2=="sfhsd" || asf2=="none")) {
      stop("Invalid value for typeAlphaSpendingNew");
    }
    
    if ((asf2=="wt" || asf2=="sfkd" || asf2=="sfhsd") && 
        R_isnancpp(asfpar2)) {
      stop("Missing value for parameterAlphaSpendingNew");
    }
    
    if (asf2=="sfkd" && asfpar2 <= 0) {
      stop ("parameterAlphaSpendingNew must be positive for sfKD");
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
      s2 = clone(t2);
    }
  }
  
  NumericVector I = IMax*t;
  LogicalVector x(kMax, 1);
  int J = L+L2;
  IntegerVector i = Range(0, J-1);
  NumericVector b = getBoundcpp(J, t[i], alpha, asf, asfpar, 0, s[i], x[i]);
  
  double lower, upper, thetahat, pvalue;
  if (!MullerSchafer) {
    double I1 = IMax*t[L-1];
    double I2 = INew*(t[J-1] - t[L-1])/(1 - t[L-1]);
    
    double r = t[L-1]/t[J-1];
    double c1 = sqrt(r)*zL + sqrt(1-r)*zL2;
    double c2 = sqrt(r)*sqrt(I1) + sqrt(1-r)*sqrt(I2);
    
    lower = (c1 - b[J-1])/c2;
    upper = (c1 + b[J-1])/c2;
    
    // point estimate is the lower bound for alpha = 0.5
    NumericVector u = getBoundcpp(J, t[i], 0.5, asf, asfpar, 0, s[i], x[i]);
    thetahat = (c1 - u[J-1])/c2;
    
    // repeated p-value is alpha for which the lower bound of theta is zero
    auto f = [J, c1, t, asf, asfpar, s, x, i](double aval)->double {
      NumericVector u = getBoundcpp(J, t[i], aval, asf, asfpar, 0, 
                                    s[i], x[i]);
      return c1 - u[J-1];
    };
    
    if (f(0.000001) > 0) {
      pvalue = 0.000001;
    } else if (f(0.999999) < 0) {
      pvalue = 0.999999;
    } else {
      pvalue = brent(f, 0.000001, 0.999999, 1.0e-6);
    }
  } else {
    double I1 = IMax*t[L-1];
    double I2 = INew*t2[L2-1];
    int k1 = kMax - L;
    
    NumericVector t1(k1), r1(k1), a1(k1, -6.0), theta0(k1);
    for (int l=0; l<k1; l++) {
      t1[l] = (t[l+L] - t[L-1])/(1 - t[L-1]);
      r1[l] = t[L-1]/t[l+L];
    }
    
    NumericVector interval(2);
    interval[0] = (zL - b[L-1])/sqrt(I1);
    interval[1] = (zL + b[L-1])/sqrt(I1);
    double tol = 0.0001;
    
    LogicalVector x2(kNew, 1);
    IntegerVector j = Range(0, L2-1);
    auto f1 = [L, zL, I1, k1, t1, r1, a1, &b, 
               L2, zL2, I2, t2, asf2, asfpar2, s2, x2, j,
               theta0](double theta)->double {
                 
                 // conditional type I error under shifted null
                 double zL1 = zL - theta*sqrt(I1);
                 NumericVector b1(k1);
                 for (int l=0; l<k1; l++) {
                   b1[l] = (b[l+L] - sqrt(r1[l])*zL1)/sqrt(1 - r1[l]);
                 }
                 
                 List probs = exitprobcpp(b1, a1, theta0, t1);
                 double alphaNew = sum(NumericVector(probs[0]));
                 
                 // efficacy boundaries for the secondary trial
                 NumericVector b2 = getBoundcpp(
                   L2, t2[j], alphaNew, asf2, asfpar2, 0, s2[j], x2[j]);
                 
                 return zL2 - theta*sqrt(I2) - b2[L2-1];
               };
    
    lower = brent(f1, interval[0], interval[1], tol);
    
    auto f2 = [L, zL, I1, k1, t1, r1, a1, &b, 
               L2, zL2, I2, t2, asf2, asfpar2, s2, x2, j,
               theta0](double theta)->double {
                 
                 // conditional type I error under shifted null
                 double zL1 = -zL + theta*sqrt(I1);
                 NumericVector b1(k1);
                 for (int l=0; l<k1; l++) {
                   b1[l] = (b[l+L] - sqrt(r1[l])*zL1)/sqrt(1 - r1[l]);
                 }
                 
                 List probs = exitprobcpp(b1, a1, theta0, t1);
                 double alphaNew = sum(NumericVector(probs[0]));
                 
                 
                 // efficacy boundaries for the secondary trial
                 NumericVector b2 = getBoundcpp(
                   L2, t2[j], alphaNew, asf2, asfpar2, 0, s2[j], x2[j]);
                 
                 return -zL2 + theta*sqrt(I2) - b2[L2-1];
               };
    
    upper = brent(f2, interval[0], interval[1], tol);
    
    // point estimate is the lower bound for alpha = 0.5
    b = getBoundcpp(J, t[i], 0.5, asf, asfpar, 0, s[i], x[i]);
    thetahat = brent(f1, lower, upper, tol);
    
    // repeated p-value is alpha for which the lower bound of theta is zero
    auto f = [J, t, asf, asfpar, s, x, i, L, zL, I1, k1, t1, r1, a1,
              L2, zL2, I2, kNew, t2, asf2, asfpar2, s2, x2, j,
              theta0, interval, tol](double aval)->double {
                NumericVector u = getBoundcpp(
                  J, t[i], aval, asf, asfpar, 0, s[i], x[i]);
                
                auto g = [L, zL, I1, k1, t1, r1, a1, &u, 
                          L2, zL2, I2, t2, asf2, asfpar2, s2, x2, j, 
                          theta0](double theta)->double {
                            
                            // conditional type I error under shifted null
                            double zL1 = zL - theta*sqrt(I1);
                            NumericVector b1(k1);
                            for (int l=0; l<k1; l++) {
                              b1[l] = (u[l+L] - sqrt(r1[l])*zL1)/
                                sqrt(1 - r1[l]);
                            }
                            
                            List probs = exitprobcpp(b1, a1, theta0, t1);
                            double alphaNew = sum(NumericVector(probs[0]));
                            
                            // efficacy boundaries for the secondary trial
                            NumericVector b2 = getBoundcpp(
                              L2, t2[j], alphaNew, asf2, asfpar2, 0, 
                              s2[j], x2[j]);
                            
                            return zL2 - theta*sqrt(I2) - b2[L2-1];
                          };
                
                return brent(g, interval[0], interval[1], tol);
              };
    
    if (f(0.000001) >= 0) {
      pvalue = 0.000001;
    } else {
      double left = 0.000001, right = 0.5;
      int count = 0;
      while (f(right) <= 0 && count <= 18) {
        left = right;
        right = (left + 1.0)/2.0;
        count++;
      }
      
      if (count <= 18) {
        pvalue = brent(f, left, right, 1.0e-6);
      } else {
        pvalue = right;
      }
    }
  }
  
  List result = List::create(
    _["pvalue"] = pvalue,
    _["thetahat"] = thetahat,
    _["cilevel"] = 1-2*alpha,
    _["lower"] = lower,
    _["upper"] = upper);
  
  return result;
}


//' @title Conditional power allowing for varying parameter values
//' @description Obtains the conditional power for specified incremental 
//' information given the interim results, parameter values, and 
//' data-dependent changes in the error spending function and the 
//' number and spacing of interim looks. 
//'
//' @param INew The maximum information for the secondary trial.
//' @param L The interim look of the primary trial.
//' @param zL The Z-test statistic at the interim look of the primary trial.
//' @param theta A scalar or a vector of parameter values of 
//'   length \code{1 + kMax - L} if \code{MullerSchafer = FALSE} or 
//'   length \code{1 + kNew} if \code{MullerSchafer = TRUE}.
//' @param kMax The maximum number of stages for the primary trial.
//' @param IMax The maximum information of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param criticalValues The upper boundaries on the Z-test statistic scale
//'   for efficacy stopping for the primary trial.
//' @param futilityBounds The lower boundaries on the Z-test statistic scale
//'   for futility stopping for the primary trial.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method 
//'   for trial adaptation.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is 
//'   allowed at each look of the secondary trial. 
//'   Defaults to true if left unspecified.
//' @param futilityStoppingNew The indicators of whether futility stopping is 
//'   allowed at each look of the secondary trial. 
//'   Defaults to true if left unspecified.
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
//' @param parameterAlphaSpendingNew The parameter value for the alpha 
//'   spending for the secondary trial. 
//'   Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".
//' @param typeBetaSpendingNew The type of beta spending for the secondary 
//'   trial. One of the following: 
//'   "sfOF" for O'Brien-Fleming type spending function, 
//'   "sfP" for Pocock type spending function, 
//'   "sfKD" for Kim & DeMets spending function, 
//'   "sfHSD" for Hwang, Shi & DeCani spending function, 
//'   "user" for user defined spending, and 
//'   "none" for no early futility stopping. 
//'   Defaults to "none".
//' @param parameterBetaSpendingNew The parameter value for the beta 
//'   spending for the secondary trial. 
//'   Corresponds to rho for "sfKD", and gamma for "sfHSD".
//' @param spendingTimeNew A vector of length \code{kNew} for the error 
//'   spending time at each analysis of the secondary trial. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRatesNew}.
//'
//' @return The conditional power given the interim results, parameter 
//' values, and data-dependent design changes.
//' 
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{getDesign}}
//' 
//' @examples
//' 
//' # Conditional power calculation with delayed treatment effect
//' 
//' # Two interim analyses have occurred with 179 and 266 events, respectively
//' # The observed hazard ratio at the second interim look is 0.81
//' 
//' trialsdt = as.Date("2020-03-04")                       # trial start date
//' iadt = c(as.Date("2022-02-01"), as.Date("2022-11-01")) # interim dates
//' mo1 = as.numeric(iadt - trialsdt + 1)/30.4375          # interim months
//' 
//' # Assume a piecewise Poisson enrollment process with a 8-month ramp-up and
//' # 521 patients were enrolled after 17.94 months
//' N = 521                   # total number of patients
//' Ta = 17.94                # enrollment duration
//' Ta1 = 8                   # assumed end of enrollment ramp-up
//' enrate = N / (Ta - Ta1/2) # enrollment rate after ramp-up
//' 
//' # Assume a median survival of 16.7 months for the control group, a 5-month
//' # delay in treatment effect, and a hazard ratio of 0.7 after the delay
//' lam1 = log(2)/16.7  # control group hazard of exponential distribution
//' t1 = 5              # months of delay in treatment effect
//' hr = 0.7            # hazard ratio after delay
//' lam2 = hr*lam1      # treatment group hazard after delay
//' 
//' # Assume an annual dropout rate of 5%
//' pc = 0.05           # annual dropout rate
//' gam = -log(1-pc)/12 # hazard for dropout
//' 
//' 
//' # The original target number of events was 298 and the new target is 3335
//' mo2 <- caltime(
//'   nevents = c(298, 335),
//'   allocationRatioPlanned = 1,
//'   accrualTime = seq(0, Ta1), 
//'   accrualIntensity = enrate*seq(1, Ta1+1)/(Ta1+1),
//'   piecewiseSurvivalTime = c(0, t1),
//'   lambda1 = c(lam1, lam2),
//'   lambda2 = c(lam1, lam1),
//'   gamma1 = gam,
//'   gamma2 = gam,
//'   accrualDuration = Ta,
//'   followupTime = 1000)
//' 
//' # expected number of events and average hazard ratios
//' (lr1 <- lrstat(
//'   time = c(mo1, mo2),
//'   accrualTime = seq(0, Ta1), 
//'   accrualIntensity = enrate*seq(1, Ta1+1)/(Ta1+1),
//'   piecewiseSurvivalTime = c(0, t1),
//'   lambda1 = c(lam1, lam2),
//'   lambda2 = c(lam1, lam1),
//'   gamma1 = gam,
//'   gamma2 = gam,
//'   accrualDuration = Ta,
//'   followupTime = 1000,
//'   predictTarget = 3))
//' 
//' 
//' hr2 = 0.81                    # observed hazard ratio at interim 2
//' z2 = (-log(hr2))*sqrt(266/4)  # corresponding Z-test statistic value
//' 
//' # Assume that the number of events is increased based on unblinded data
//' # Use boundaries based on the original sample size for the CHW statistics
//' b = getBound(k = 3, informationRates = c(179, 266, 298)/298,
//'              alpha = 0.025, typeAlphaSpending = "sfOF")
//' 
//' # expected mean of -log(HR) at interim and final for the new sample size
//' theta = -log(lr1$HR[c(2,4)])
//' 
//' # conditional power for the CHW statistic to cross the boundary at final
//' getCP(INew = (335 - 266)/4, 
//'       L = 2, zL = z2, theta = theta,
//'       kMax = 3, IMax = 298/4, 
//'       informationRates = c(179, 266, 298)/298,
//'       criticalValues = b)
//'
//' @export
// [[Rcpp::export]]
double getCP(double INew = NA_REAL, 
             const int L = NA_INTEGER, 
             const double zL = NA_REAL, 
             const NumericVector& theta = NA_REAL, 
             const int kMax = NA_INTEGER, 
             const double IMax = NA_REAL,
             const NumericVector& informationRates = NA_REAL, 
             const NumericVector& criticalValues = NA_REAL, 
             const NumericVector& futilityBounds = NA_REAL,
             const bool MullerSchafer = 0, 
             const int kNew = NA_INTEGER, 
             const NumericVector& informationRatesNew = NA_REAL, 
             const LogicalVector& efficacyStoppingNew = NA_LOGICAL, 
             const LogicalVector& futilityStoppingNew = NA_LOGICAL,
             const String typeAlphaSpendingNew = "sfOF", 
             const double parameterAlphaSpendingNew = NA_REAL, 
             const String typeBetaSpendingNew = "none", 
             const double parameterBetaSpendingNew = NA_REAL, 
             const NumericVector& spendingTimeNew = NA_REAL) {
  
  NumericVector t = clone(informationRates);
  NumericVector futilityBounds1 = clone(futilityBounds);
  NumericVector tNew = clone(informationRatesNew);
  LogicalVector efficacyStopping1 = clone(efficacyStoppingNew);
  LogicalVector futilityStopping1 = clone(futilityStoppingNew);
  NumericVector spendingTime1 = clone(spendingTimeNew);
  
  std::string asf = typeAlphaSpendingNew;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfpar = parameterAlphaSpendingNew;
  
  std::string bsf = typeBetaSpendingNew;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double bsfpar = parameterBetaSpendingNew;
  
  if (R_isnancpp(INew)) {
    stop("INew must be provided");
  }
  
  if (!R_isnancpp(INew) && INew <= 0) {
    stop("INew must be positive");
  }
  
  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (L <= 0) {
    stop("L must be a positive integer");
  }
  
  if (is_true(any(is_na(theta)))) {
    stop("theta must be provided");
  }
  
  if (R_isnancpp(kMax)) {
    stop("kMax must be provided");
  }
  
  if (kMax <= L) {
    stop("kMax must be greater than L");
  }

  if (R_isnancpp(IMax)) {
    stop("IMax must be provided");
  }
  
  if (IMax <= 0) {
    stop("IMax must be positive");
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
    t = as<NumericVector>(tem)/(kMax+0.0);
  }
  
  if (is_true(any(is_na(criticalValues)))) {
    stop("criticalValues must be provided");
  } else if (criticalValues.size() != kMax) {
    stop("Invalid length for criticalValues");
  }
  
  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  } else {
    futilityBounds1 = rep(-6.0, kMax);
    futilityBounds1[kMax-1] = criticalValues[kMax-1];
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
  
  if (MullerSchafer) {
    if (R_isnancpp(kNew)) {
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
      tNew = as<NumericVector>(tem)/(kNew+0.0);
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
      efficacyStopping1 = rep(1, kNew);
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
      futilityStopping1 = rep(1, kNew);
    }
    
    if (!(asf=="of" || asf=="p" || asf=="wt" || asf=="sfof" || asf=="sfp" ||
        asf=="sfkd" || asf=="sfhsd" || asf=="none")) {
      stop("Invalid value for typeAlphaSpendingNew");
    }
    
    if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && R_isnancpp(asfpar)) {
      stop("Missing value for parameterAlphaSpendingNew");
    }
    
    if (asf=="sfkd" && asfpar <= 0) {
      stop ("parameterAlphaSpendingNew must be positive for sfKD");
    }
    
    if (!(bsf=="sfof" || bsf=="sfp" || bsf=="sfkd" || 
        bsf=="sfhsd" || bsf=="none")) {
      stop("Invalid value for typeBetaSpendingNew");
    }
    
    if ((bsf=="sfkd" || bsf=="sfhsd") && R_isnancpp(bsfpar)) {
      stop("Missing value for parameterBetaSpendingNew");
    }
    
    if (bsf=="sfkd" && bsfpar <= 0) {
      stop ("parameterBetaSpendingNew must be positive for sfKD");
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
      spendingTime1 = clone(tNew);
    }
  }
  
  int k1 = kMax - L;
  NumericVector t1(k1), r1(k1), b1(k1), a1(k1, -6.0), zero(k1);
  for (int l=0; l<k1; l++) {
    t1[l] = (t[l+L] - t[L-1])/(1 - t[L-1]);
    r1[l] = t[L-1]/t[l+L];
    b1[l] = (criticalValues[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
  }
  
  double result;
  if (!MullerSchafer) {
    NumericVector theta1(k1+1);
    if (theta.size() == 1) {
      theta1.fill(theta[0]);
    } else if (theta.size() == k1+1){
      theta1 = clone(theta);
    } else {
      stop("Invalid length for theta"); 
    }
    
    for (int l=0; l<k1; l++) {
      a1[l] = (futilityBounds1[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
    }
    
    NumericVector mu(k1), I2(k1);
    for (int l=0; l<k1; l++) {
      double r = IMax*t[L-1]/(IMax*t[L-1] + INew*t1[l]);
      mu[l] = (theta1[l+1] - r*theta1[0])/(1 - r);
      I2[l] = INew*t1[l];
    }
    
    List probs = exitprobcpp(b1, a1, mu, I2);
    result = sum(NumericVector(probs[0]));
  } else {
    NumericVector theta1(kNew+1);
    if (theta.size() == 1) {
      theta1.fill(theta[0]);
    } else if (theta.size() == kNew+1) {
      theta1 = clone(theta);
    } else {
      stop("Invalid length for theta");
    }
    
    // conditional type I error
    List probs = exitprobcpp(b1, a1, zero, t1);
    double alphaNew = sum(NumericVector(probs[0]));
    
    // obtain the efficacy boundaries for the secondary trial
    NumericVector b2 = getBoundcpp(kNew, tNew, alphaNew, asf, asfpar, 
                                   0, spendingTime1, efficacyStopping1);
    
    // obtain the conditional power
    NumericVector mu(kNew), I2(kNew);
    for (int l=0; l<kNew; l++) {
      double r = IMax*t[L-1]/(IMax*t[L-1] + INew*tNew[l]);
      mu[l] = (theta1[l+1] - r*theta1[0])/(1 - r);
      I2[l] = INew*tNew[l];
    }
    
    List out = getPower(alphaNew, kNew, b2, mu, I2, bsf, bsfpar, 
                        spendingTime1, futilityStopping1);
    result = out[0];
  }
  
  return result;
}

