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
                const int L = NA_INTEGER,
                const double zL = NA_REAL,
                const NumericVector& b = NA_REAL,
                const NumericVector& I = NA_REAL) {
  
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
//' @param L The termination look.
//' @param zL The z-test statistic at the termination look.
//' @param IMax The maximum information of the trial.
//' @param informationRates The information rates up to look \code{L}.
//' @param efficacyStopping Indicators of whether efficacy stopping is 
//'   allowed at each stage up to look \code{L}. 
//'   Defaults to true if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping up to look \code{L}.
//' @inheritParams param_alpha
//' @param typeAlphaSpending The type of alpha spending. 
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
//' @param parameterAlphaSpending The parameter value of alpha spending. 
//'   Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".
//' @param spendingTime The error spending time up to look \code{L}. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRates}.
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
//' @references 
//' Anastasios A. Tsiatis, Gary L. Rosner and Cyrus R. Mehta. 
//' Exact confidence intervals following a group sequential test. 
//' Biometrics 1984;40:797-803.
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
//' # confidence interval
//' getCI(L = L, zL = zL, IMax = n/(4*sigma1^2), 
//'       informationRates = c(1/3, 2/3), alpha = 0.05, 
//'       typeAlphaSpending = "sfHSD", parameterAlphaSpending = -4)
//' 
//' @export
// [[Rcpp::export]]
List getCI(const int L = NA_INTEGER,
           const double zL = NA_REAL,
           const double IMax = NA_REAL,
           const NumericVector& informationRates = NA_REAL,
           const LogicalVector& efficacyStopping = NA_LOGICAL,
           const NumericVector& criticalValues = NA_REAL,
           const double alpha = 0.025,
           const String typeAlphaSpending = "sfOF",
           const double parameterAlphaSpending = NA_REAL,
           const NumericVector& spendingTime = NA_REAL) {
  
  NumericVector t = clone(informationRates);
  LogicalVector es = clone(efficacyStopping);
  NumericVector b = clone(criticalValues);
  NumericVector st = clone(spendingTime);

  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (L < 1) {
    stop("L must be a positive integer");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (R_isnancpp(IMax)) {
    stop("IMax must be provided");
  }
  
  if (IMax <= 0) {
    stop("IMax must be positive");
  }
  
  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != L) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (L > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[L-1] > 1) {
      stop("informationRates must not exceed 1");
    }
  } else {
    stop("informationRates must be provided");
  }
  
  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != L) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[L-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) | 
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    es = rep(1, L);
  }
  
  
  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != L) {
      stop("Invalid length for criticalValues");
    }
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
  
  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
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
    if (spendingTime.size() != L) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (L > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[L-1] > 1) {
      stop("spendingTime must not exceed 1");
    }
  } else {
    st = clone(t);
  }
  
  if (is_true(any(is_na(criticalValues)))) {
    b = getBoundcpp(L, t, alpha, asf, asfpar, 0, st, es);
  }
  
  NumericVector I = IMax*t;
  
  double pvalue = f_pvalue(0, L, zL, b, I);
  
  double cilevel = 1-2*alpha;
  
  NumericVector interval(2);
  interval[0] = (zL - 6)/sqrt(I[L-1]);
  interval[1] = (zL + 6)/sqrt(I[L-1]);
  double tol = 0.0001;
  
  auto f = [L, zL, b, I](double theta)->double {
    return f_pvalue(theta, L, zL, b, I) - 0.5;
  };
  double thetahat = brent(f, interval[0], interval[1], tol);
  
  auto f1 = [L, zL, b, I, cilevel](double theta)->double {
    return f_pvalue(theta, L, zL, b, I) - (1-cilevel)/2;
  };
  double lower = brent(f1, interval[0], thetahat, tol);
  
  auto f2 = [L, zL, b, I, cilevel](double theta)->double {
    return f_pvalue(theta, L, zL, b, I) - (1+cilevel)/2;
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
//' @param L The look of interest.
//' @param zL The z-test statistic at the look.
//' @param IMax The maximum information of the trial.
//' @param informationRates The information rates up to look \code{L}.
//' @param efficacyStopping Indicators of whether efficacy stopping is 
//'   allowed at each stage up to look \code{L}. Defaults to true 
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping up to look \code{L}.
//' @inheritParams param_alpha
//' @param typeAlphaSpending The type of alpha spending. 
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
//' @param parameterAlphaSpending The parameter value of alpha spending. 
//'   Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".
//' @param spendingTime The error spending time up to look \code{L}. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRates}.
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
//' @references 
//' Christopher Jennison and Bruce W. Turnbull. 
//' Interim analyses: the repeated confidence interval approach 
//' (with discussion). 
//' J R Stat Soc Series B. 1989;51:305-361.
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
//' # repeated confidence interval
//' getRCI(L = L, zL = zL, IMax = n/(4*sigma1^2), 
//'        informationRates = c(1/3, 2/3), alpha = 0.05, 
//'        typeAlphaSpending = "sfHSD", parameterAlphaSpending = -4)
//' 
//' @export
// [[Rcpp::export]]
List getRCI(const int L = NA_INTEGER,
            const double zL = NA_REAL,
            const double IMax = NA_REAL,
            const NumericVector& informationRates = NA_REAL,
            const LogicalVector& efficacyStopping = NA_LOGICAL,
            const NumericVector& criticalValues = NA_REAL,
            const double alpha = 0.025,
            const String typeAlphaSpending = "sfOF",
            const double parameterAlphaSpending = NA_REAL,
            const NumericVector& spendingTime = NA_REAL) {
  
  NumericVector t = clone(informationRates);
  LogicalVector es = clone(efficacyStopping);
  NumericVector b = clone(criticalValues);
  NumericVector st = clone(spendingTime);
  
  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (L < 1) {
    stop("L must be a positive integer");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (R_isnancpp(IMax)) {
    stop("IMax must be provided");
  }
  
  if (IMax <= 0) {
    stop("IMax must be positive");
  }
  
  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != L) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (L > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[L-1] > 1) {
      stop("informationRates must not exceed 1");
    }
  } else {
    stop("informationRates must be provided");
  }
  
  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != L) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[L-1] != 1) {
      stop("efficacyStopping must end with 1");
    } else if (is_false(all((efficacyStopping == 1) | 
      (efficacyStopping == 0)))) {
      stop("Elements of efficacyStopping must be 1 or 0");
    }
  } else {
    es = rep(1, L);
  }
  
  
  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != L) {
      stop("Invalid length for criticalValues");
    }
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
  
  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
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
    if (spendingTime.size() != L) {
      stop("Invalid length for spendingTime");
    } else if (spendingTime[0] <= 0) {
      stop("Elements of spendingTime must be positive");
    } else if (L > 1 && is_true(any(diff(spendingTime) <= 0))) {
      stop("Elements of spendingTime must be increasing");
    } else if (spendingTime[L-1] > 1) {
      stop("spendingTime must not exceed 1");
    }
  } else {
    st = clone(t);
  }
  
  if (is_true(any(is_na(criticalValues)))) {
    b = getBoundcpp(L, t, alpha, asf, asfpar, 0, st, es);
  }
  
  NumericVector I = IMax*t;
  
  // repeated confidence interval
  double lower = (zL - b[L-1])/sqrt(I[L-1]);
  double upper = (zL + b[L-1])/sqrt(I[L-1]);
  
  // point estimate is the lower bound for alpha = 0.5
  NumericVector u = getBoundcpp(L, t, 0.5, asf, asfpar, 0, st, es);
  double thetahat = (zL - u[L-1])/sqrt(I[L-1]);
  
  // repeated p-value is alpha for which the lower bound of theta is zero
  auto f = [L, zL, t, asf, asfpar, st, es](double aval)->double {
    NumericVector u = getBoundcpp(L, t, aval, asf, asfpar, 0, st, es);  
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
               const int L2,
               const double zL2,
               const NumericVector& b2,
               const NumericVector& I2) {
  
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
               const int L,
               const double zL,
               const NumericVector& b,
               const NumericVector& I,
               const int L2,
               const double zL2,
               const NumericVector& b2,
               const NumericVector& I2) {
  
  double astar = f_astar(theta, L2, zL2, b2, I2);
  int k1 = kMax - L;
  
  NumericVector b1(k1), a1(k1, -6.0), mu(k1, theta), I1(k1);
  for (int l=0; l<k1; l++) {
    b1[l] = (b[l+L] - sqrt(I[L-1]/I[l+L])*zL)/sqrt(1 - I[L-1]/I[l+L]);
    I1[l] = I[l+L] - I[L-1];
  }
  
  List probs = exitprobcpp(b1, a1, mu, I1);
  NumericVector pu = NumericVector(probs[0]);
  
  // find the interval that contains the rejection probability 
  // in the secondary trial under the shift null
  NumericVector p(k1+1);
  p[0] = 0;
  for (int l=0; l<k1; l++) {
    p[l+1] = p[l] + pu[l];
  }
  
  NumericVector astars(1, astar);
  IntegerVector js = findInterval2(astars, p);
  int j = js[0];
  
  // find the z-test statistic value yielding the rejection probability
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
                  const int L = NA_INTEGER,
                  const double zL = NA_REAL,
                  const NumericVector& b = NA_REAL,
                  const NumericVector& I = NA_REAL,
                  const int L2 = NA_INTEGER,
                  const double zL2 = NA_REAL,
                  const NumericVector& b2 = NA_REAL,
                  const NumericVector& I2 = NA_REAL) {
  
  List bw = f_bwimage(theta, kMax, L, zL, b, I, L2, zL2, b2, I2);
  
  int J = bw[0];
  double zJ = bw[1];
  
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
//' confidence interval after the end of an adaptive trial.
//'
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of 
//'   the primary trial.
//' @param IMax The maximum information of the primary trial.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is 
//'   allowed at each stage of the primary trial. Defaults to true 
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping for the primary trial.
//' @param alpha The significance level of the primary trial. 
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
//' @param parameterAlphaSpending The parameter value of alpha spending 
//'   for the primary trial. Corresponds to Delta for "WT", rho for "sfKD", 
//'   and gamma for "sfHSD".
//' @param spendingTime The error spending time of the primary trial. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRates}.
//' @param L2 The termination look of the secondary trial.
//' @param zL2 The z-test statistic at the termination look of the 
//'   secondary trial.
//' @param INew The maximum information of the secondary trial.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method 
//'   for trial adaptation.
//' @param informationRatesNew The spacing of looks of the secondary trial
//'   up to look \code{L2}.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is 
//'   allowed at each look of the secondary trial up to look \code{L2}. 
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
//' @param parameterAlphaSpendingNew The parameter value of alpha spending 
//'   for the secondary trial. Corresponds to Delta for "WT", 
//'   rho for "sfKD", and gamma for "sfHSD".
//' @param spendingTimeNew The error spending time of the secondary trial
//'   up to look \code{L2}. Defaults to missing, in which case, it is 
//'   the same as \code{informationRatesNew}.
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
//' @references
//' Ping Gao, Lingyun Liu and Cyrus Mehta. 
//' Exact inference for adaptive group sequential designs. 
//' Stat Med. 2013;32(23):3991-4005.
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
//' t = des1$byStageResults$informationRates
//' 
//' # Muller & Schafer (2001) method to design the secondary trial: 
//' des2 = adaptDesign(
//'   betaNew = 0.2, L = L, zL = zL, theta = 5,
//'   kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpendingNew = "sfHSD", 
//'   parameterAlphaSpendingNew = -2)
//' 
//' n2 = ceiling(des2$secondaryTrial$overallResults$maxInformation*4*20^2)
//' ns = round(n2*(1:3)/3)
//'  (des2 = adaptDesign(
//'    INew = n2/(4*20^2), L = L, zL = zL, theta = 5,
//'    kMax = 3, informationRates = t,
//'    alpha = 0.05, typeAlphaSpending = "sfHSD",
//'    parameterAlphaSpending = -4,
//'    MullerSchafer = TRUE,
//'    kNew = 3, informationRatesNew = ns/n2,
//'    typeAlphaSpendingNew = "sfHSD",
//'    parameterAlphaSpendingNew = -2))
//' 
//' # termination at the second look of the secondary trial
//' L2 = 2
//' delta2 = 6.86
//' sigma2 = 21.77
//' zL2 = delta2/sqrt(4/197*sigma2^2)
//' 
//' t2 = des2$secondaryTrial$byStageResults$informationRates[1:L2]
//' 
//' # confidence interval
//' getADCI(L = L, zL = zL,
//'         IMax = n/(4*sigma1^2), kMax = 3,
//'         informationRates = t,
//'         alpha = 0.05, typeAlphaSpending = "sfHSD",
//'         parameterAlphaSpending = -4,
//'         L2 = L2, zL2 = zL2,
//'         INew = n2/(4*sigma2^2),
//'         MullerSchafer = TRUE,
//'         informationRatesNew = t2, 
//'         typeAlphaSpendingNew = "sfHSD",
//'         parameterAlphaSpendingNew = -2)
//' 
//' @export
// [[Rcpp::export]]
List getADCI(const int L = NA_INTEGER,
             const double zL = NA_REAL,
             const double IMax = NA_REAL,
             const int kMax = NA_INTEGER,
             const NumericVector& informationRates = NA_REAL,
             const LogicalVector& efficacyStopping = NA_LOGICAL,
             const NumericVector& criticalValues = NA_REAL,
             const double alpha = 0.25,
             const String typeAlphaSpending = "sfOF",
             const double parameterAlphaSpending = NA_REAL,
             const NumericVector& spendingTime = NA_REAL,
             const int L2 = NA_INTEGER,
             const double zL2 = NA_REAL,
             const double INew = NA_REAL,
             const bool MullerSchafer = 0,
             const NumericVector& informationRatesNew = NA_REAL,
             const LogicalVector& efficacyStoppingNew = NA_LOGICAL,
             const String typeAlphaSpendingNew = "sfOF",
             const double parameterAlphaSpendingNew = NA_REAL,
             const NumericVector& spendingTimeNew = NA_REAL) {
  
  NumericVector t = clone(informationRates);
  LogicalVector es = clone(efficacyStopping);
  NumericVector b = clone(criticalValues);
  NumericVector st = clone(spendingTime);
  NumericVector tNew = clone(informationRatesNew);
  LogicalVector esNew = clone(efficacyStoppingNew);
  NumericVector stNew = clone(spendingTimeNew);
  double alpha1 = alpha;
  
  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfpar = parameterAlphaSpending;
  
  std::string asfNew = typeAlphaSpendingNew;
  std::for_each(asfNew.begin(), asfNew.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfparNew = parameterAlphaSpendingNew;
  
  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (L < 1) {
    stop("L must be a positive integer");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (R_isnancpp(IMax)) {
    stop("IMax must be provided");
  }
  
  if (IMax <= 0) {
    stop("IMax must be positive");
  }
  
  if (R_isnancpp(kMax)) {
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
    t = as<NumericVector>(tem)/(kMax+0.0);
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
  
  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
    
    NumericVector u(kMax), l(kMax, -6.0), theta0(kMax);
    for (int i=0; i<kMax; i++) {
      u[i] = criticalValues[i];
      if (!es[i]) u[i] = 6.0;
    }
    
    List probs = exitprobcpp(u, l, theta0, t);
    alpha1 = sum(NumericVector(probs[0]));
  }
  
  if (!R_isnancpp(alpha1)) {
    if (alpha1 < 0.00001 || alpha1 >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
  } else {
    stop("alpha must be provided for missing criticalValues");
  }
  
  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
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
    st = clone(t);
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
  
  if (R_isnancpp(INew)) {
    stop("INew must be provided");
  }
  
  if (INew <= 0) {
    stop("INew must be positive");
  }
  
  
  if (MullerSchafer) {
    if (is_false(any(is_na(informationRatesNew)))) {
      if (informationRatesNew.size() != L2) {
        stop("Invalid length for informationRatesNew");
      } else if (informationRatesNew[0] <= 0) {
        stop("Elements of informationRatesNew must be positive");
      } else if (L2 > 1 && is_true(any(diff(informationRatesNew) <= 0))) {
        stop("Elements of informationRatesNew must be increasing");
      } else if (informationRatesNew[L2-1] > 1) {
        stop("informationRatesNew must not exceed 1");
      }
    } else {
      stop("informationRatesNew must be provided");
    }
    
    if (is_false(any(is_na(efficacyStoppingNew)))) {
      if (efficacyStoppingNew.size() != L2) {
        stop("Invalid length for efficacyStoppingNew");
      } else if (efficacyStoppingNew[L2-1] != 1) {
        stop("efficacyStoppingNew must end with 1");
      } else if (is_false(all((efficacyStoppingNew == 1) | 
        (efficacyStoppingNew == 0)))) {
        stop("Elements of efficacyStoppingNew must be 1 or 0");
      }
    } else {
      esNew = rep(1, L2);
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
    
    if (is_false(any(is_na(spendingTimeNew)))) {
      if (spendingTimeNew.size() != L2) {
        stop("Invalid length for spendingTimeNew");
      } else if (spendingTimeNew[0] <= 0) {
        stop("Elements of spendingTimeNew must be positive");
      } else if (L2 > 1 && is_true(any(diff(spendingTimeNew) <= 0))) {
        stop("Elements of spendingTimeNew must be increasing");
      } else if (spendingTimeNew[L2-1] > 1) {
        stop("spendingTimeNew must not exceed 1");
      }
    } else {
      stNew = clone(tNew);
    }
  }
  
  // efficacy boundaries for the primary trial
  if (is_true(any(is_na(criticalValues)))) {
    b = getBoundcpp(kMax, t, alpha1, asf, asfpar, 0, st, es);
  }
  
  NumericVector I = IMax*t;
  
  NumericVector b2(L2), I2(L2);
  if (!MullerSchafer) {
    NumericVector t1(L2), r1(L2);
    for (int l=0; l<L2; l++) {
      t1[l] = (t[l+L] - t[L-1])/(1 - t[L-1]);
      r1[l] = t[L-1]/t[l+L];
      b2[l] = (b[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
      if (!es[l+L]) b2[l] = 6.0;
      I2[l] = INew*t1[l];
    }
  } else {
    // conditional type I error
    int k1 = kMax - L;
    NumericVector t1(k1), r1(k1), b1(k1), a1(k1, -6.0), theta0(k1);
    for (int l=0; l<k1; l++) {
      t1[l] = (t[l+L] - t[L-1])/(1 - t[L-1]);
      r1[l] = t[L-1]/t[l+L];
      b1[l] = (b[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
      if (!es[l+L]) b1[l] = 6.0;
    }
    
    List probs = exitprobcpp(b1, a1, theta0, t1);
    double alphaNew = sum(NumericVector(probs[0]));
    
    // efficacy boundaries for the secondary trial
    b2 = getBoundcpp(L2, tNew, alphaNew, asfNew, asfparNew, 0, stNew, esNew);
    
    for (int l=0; l<L2; l++) {
      I2[l] = INew*tNew[l];
    }
  }
  
  double cilevel = 1 - 2*alpha1;
  
  int K = kMax;
  double pvalue = f_bwpvalue(0, K, L, zL, b, I, L2, zL2, b2, I2);
  
  NumericVector interval(2);
  interval[0] = (zL - b[L-1])/sqrt(I[L-1]);
  interval[1] = (zL + b[L-1])/sqrt(I[L-1]);
  double tol = 0.0001;
  
  auto f = [K, L, zL, b, I, L2, zL2, b2, I2](double theta)->double {
    return f_bwpvalue(theta, K, L, zL, b, I, L2, zL2, b2, I2) - 0.5;
  };
  double thetahat = brent(f, interval[0], interval[1], tol);
  
  auto f1 = [K, L, zL, b, I, L2, zL2, b2, I2, cilevel](double theta)->double {
    return f_bwpvalue(theta, K, L, zL, b, I, L2, zL2, b2, I2) - (1-cilevel)/2;
  };
  double lower = brent(f1, interval[0], thetahat, tol);
  
  auto f2 = [K, L, zL, b, I, L2, zL2, b2, I2, cilevel](double theta)->double {
    return f_bwpvalue(theta, K, L, zL, b, I, L2, zL2, b2, I2) - (1+cilevel)/2;
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
//' @description Obtains the repeated p-value, conservative point estimate, 
//' and repeated confidence interval for an adaptive group sequential trial.
//' 
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of 
//'   the primary trial.
//' @param IMax The maximum information of the primary trial.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is 
//'   allowed at each stage of the primary trial. Defaults to true 
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the z-test statistic scale
//'   for efficacy stopping for the primary trial.
//' @param alpha The significance level of the primary trial. 
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
//' @param parameterAlphaSpending The parameter value of alpha spending 
//'   for the primary trial. Corresponds to Delta for "WT", rho for "sfKD", 
//'   and gamma for "sfHSD".
//' @param spendingTime The error spending time of the primary trial. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRates}.
//' @param L2 The look of interest in the secondary trial.
//' @param zL2 The z-test statistic at the look of the secondary trial.
//' @param INew The maximum information of the secondary trial.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method 
//'   for trial adaptation.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is 
//'   allowed at each look of the secondary trial up to look \code{L2}. 
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
//' @param parameterAlphaSpendingNew The parameter value of alpha spending 
//'   for the secondary trial. Corresponds to Delta for "WT", 
//'   rho for "sfKD", and gamma for "sfHSD".
//' @param spendingTimeNew The error spending time of the secondary trial. 
//'   up to look \code{L2}. Defaults to missing, in which case, it is 
//'   the same as \code{informationRatesNew}.
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
//' @references
//' Cyrus R. Mehta, Peter Bauer, Martin Posch and Werner Brannath.
//' Repeated confidence intervals for adaptive group sequential trials.
//' Stat Med. 2007;26:5422–5433.
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
//' t = des1$byStageResults$informationRates
//' 
//' # Muller & Schafer (2001) method to design the secondary trial: 
//' des2 = adaptDesign(
//'   betaNew = 0.2, L = L, zL = zL, theta = 5,
//'   kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, typeAlphaSpendingNew = "sfHSD", 
//'   parameterAlphaSpendingNew = -2)
//' 
//' n2 = ceiling(des2$secondaryTrial$overallResults$maxInformation*4*20^2)
//' ns = round(n2*(1:3)/3)
//' (des2 = adaptDesign(
//'   INew = n2/(4*20^2), L = L, zL = zL, theta = 5,
//'   kMax = 3, informationRates = t,
//'   alpha = 0.05, typeAlphaSpending = "sfHSD",
//'   parameterAlphaSpending = -4,
//'   MullerSchafer = TRUE,
//'   kNew = 3, informationRatesNew = ns/n2,
//'   typeAlphaSpendingNew = "sfHSD",
//'   parameterAlphaSpendingNew = -2))
//' 
//' # termination at the second look of the secondary trial
//' L2 = 2
//' delta2 = 6.86
//' sigma2 = 21.77
//' zL2 = delta2/sqrt(4/197*sigma2^2)
//' 
//' t2 = des2$secondaryTrial$byStageResults$informationRates[1:L2]
//' 
//' # repeated confidence interval
//' getADRCI(L = L, zL = zL,
//'          IMax = n/(4*sigma1^2), kMax = 3,
//'          informationRates = t,
//'          alpha = 0.05, typeAlphaSpending = "sfHSD",
//'          parameterAlphaSpending = -4,
//'          L2 = L2, zL2 = zL2,
//'          INew = n2/(4*sigma2^2),
//'          MullerSchafer = TRUE,
//'          informationRatesNew = t2, 
//'          typeAlphaSpendingNew = "sfHSD",
//'          parameterAlphaSpendingNew = -2)
//' 
//' @export
// [[Rcpp::export]]
List getADRCI(const int L = NA_INTEGER,
              const double zL = NA_REAL,
              const double IMax = NA_REAL,
              const int kMax = NA_INTEGER,
              const NumericVector& informationRates = NA_REAL,
              const LogicalVector& efficacyStopping = NA_LOGICAL,
              const NumericVector& criticalValues = NA_REAL,
              const double alpha = 0.25,
              const String typeAlphaSpending = "sfOF",
              const double parameterAlphaSpending = NA_REAL,
              const NumericVector& spendingTime = NA_REAL,
              const int L2 = NA_INTEGER,
              const double zL2 = NA_REAL,
              const double INew = NA_REAL,
              const bool MullerSchafer = 0,
              const NumericVector& informationRatesNew = NA_REAL,
              const LogicalVector& efficacyStoppingNew = NA_LOGICAL,
              const String typeAlphaSpendingNew = "sfOF",
              const double parameterAlphaSpendingNew = NA_REAL,
              const NumericVector& spendingTimeNew = NA_REAL) {
  
  NumericVector t = clone(informationRates);
  LogicalVector es = clone(efficacyStopping);
  NumericVector b = clone(criticalValues);
  NumericVector st = clone(spendingTime);
  NumericVector tNew = clone(informationRatesNew);
  LogicalVector esNew = clone(efficacyStoppingNew);
  NumericVector stNew = clone(spendingTimeNew);
  double alpha1 = alpha;
  
  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfpar = parameterAlphaSpending;
  
  std::string asfNew = typeAlphaSpendingNew;
  std::for_each(asfNew.begin(), asfNew.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfparNew = parameterAlphaSpendingNew;
  
  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (L < 1) {
    stop("L must be a positive integer");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (R_isnancpp(IMax)) {
    stop("IMax must be provided");
  }
  
  if (IMax <= 0) {
    stop("IMax must be positive");
  }
  
  if (R_isnancpp(kMax)) {
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
    t = as<NumericVector>(tem)/(kMax+0.0);
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
  
  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
    
    NumericVector u(kMax), l(kMax, -6.0), theta0(kMax);
    for (int i=0; i<kMax; i++) {
      u[i] = criticalValues[i];
      if (!es[i]) u[i] = 6.0;
    }
    
    List probs = exitprobcpp(u, l, theta0, t);
    alpha1 = sum(NumericVector(probs[0]));
  }
  
  if (!R_isnancpp(alpha1)) {
    if (alpha1 < 0.00001 || alpha1 >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
  } else {
    stop("alpha must be provided for missing criticalValues");
  }
  
  if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
      asf=="wt" || asf=="sfof" || asf=="sfp" ||
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
    st = clone(t);
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
  
  if (R_isnancpp(INew)) {
    stop("INew must be provided");
  }
  
  if (INew <= 0) {
    stop("INew must be positive");
  }
  
  
  if (MullerSchafer) {
    if (is_false(any(is_na(informationRatesNew)))) {
      if (informationRatesNew.size() != L2) {
        stop("Invalid length for informationRatesNew");
      } else if (informationRatesNew[0] <= 0) {
        stop("Elements of informationRatesNew must be positive");
      } else if (L2 > 1 && is_true(any(diff(informationRatesNew) <= 0))) {
        stop("Elements of informationRatesNew must be increasing");
      } else if (informationRatesNew[L2-1] > 1) {
        stop("informationRatesNew must not exceed 1");
      }
    } else {
      stop("informationRatesNew must be provided");
    }
    
    if (is_false(any(is_na(efficacyStoppingNew)))) {
      if (efficacyStoppingNew.size() != L2) {
        stop("Invalid length for efficacyStoppingNew");
      } else if (efficacyStoppingNew[L2-1] != 1) {
        stop("efficacyStoppingNew must end with 1");
      } else if (is_false(all((efficacyStoppingNew == 1) | 
        (efficacyStoppingNew == 0)))) {
        stop("Elements of efficacyStoppingNew must be 1 or 0");
      }
    } else {
      esNew = rep(1, L2);
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
    
    if (is_false(any(is_na(spendingTimeNew)))) {
      if (spendingTimeNew.size() != L2) {
        stop("Invalid length for spendingTimeNew");
      } else if (spendingTimeNew[0] <= 0) {
        stop("Elements of spendingTimeNew must be positive");
      } else if (L2 > 1 && is_true(any(diff(spendingTimeNew) <= 0))) {
        stop("Elements of spendingTimeNew must be increasing");
      } else if (spendingTimeNew[L2-1] > 1) {
        stop("spendingTimeNew must not exceed 1");
      }
    } else {
      stNew = clone(tNew);
    }
  }
  
  // efficacy boundaries for the primary trial
  if (is_true(any(is_na(criticalValues)))) {
    b = getBoundcpp(kMax, t, alpha1, asf, asfpar, 0, st, es);
  }
  
  NumericVector I = IMax*t;
  
  double lower, upper, thetahat, pvalue;
  if (!MullerSchafer) {
    double I1 = IMax*t[L-1];
    double I2 = INew*(t[L+L2-1] - t[L-1])/(1 - t[L-1]);
    
    double r1 = t[L-1]/t[L+L2-1];
    double c1 = sqrt(r1)*zL + sqrt(1-r1)*zL2;
    double c2 = sqrt(r1)*sqrt(I1) + sqrt(1-r1)*sqrt(I2);
    
    lower = (c1 - b[L+L2-1])/c2;
    upper = (c1 + b[L+L2-1])/c2;
    
    // point estimate is the lower bound for alpha = 0.5
    int J = L+L2;
    IntegerVector i = Range(0, J-1);
    NumericVector u = getBoundcpp(J, t[i], 0.5, asf, asfpar, 0, st[i], es[i]);
    thetahat = (c1 - u[J-1])/c2;
    
    // repeated p-value is alpha for which the lower bound of theta is zero
    auto f = [J, c1, t, asf, asfpar, st, es, i](double aval)->double {
      NumericVector u = getBoundcpp(
        J, t[i], aval, asf, asfpar, 0, st[i], es[i]);
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
    double I2 = INew*tNew[L2-1];
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
    
    // point estimate is the lower bound for alpha = 0.5
    NumericVector u = getBoundcpp(kMax, t, 0.5, asf, asfpar, 0, st, es);
    
    auto f0 = [L, zL, I1, k1, t1, r1, u, es, a1, theta0, 
               L2, zL2, I2, tNew, asfNew, asfparNew, stNew, 
               esNew](double theta)->double {
                 
                 // conditional type I error under shifted null
                 double zL1 = zL - theta*sqrt(I1);
                 NumericVector b1(k1);
                 for (int l=0; l<k1; l++) {
                   b1[l] = (u[l+L] - sqrt(r1[l])*zL1)/sqrt(1 - r1[l]);
                   if (!es[l+L]) b1[l] = 6.0;
                 }
                 
                 List probs = exitprobcpp(b1, a1, theta0, t1);
                 double alphaNew = sum(NumericVector(probs[0]));
                 
                 // efficacy boundaries for the secondary trial
                 NumericVector b2 = getBoundcpp(
                   L2, tNew, alphaNew, asfNew, asfparNew, 0, stNew, esNew);
                 
                 return zL2 - theta*sqrt(I2) - b2[L2-1];
               };
    
    thetahat = brent(f0, interval[0], interval[1], tol);

    auto f1 = [L, zL, I1, k1, t1, r1, b, es, a1, theta0, 
               L2, zL2, I2, tNew, asfNew, asfparNew, stNew, 
               esNew](double theta)->double {
                 
                 // conditional type I error under shifted null
                 double zL1 = zL - theta*sqrt(I1);
                 NumericVector b1(k1);
                 for (int l=0; l<k1; l++) {
                   b1[l] = (b[l+L] - sqrt(r1[l])*zL1)/sqrt(1 - r1[l]);
                   if (!es[l+L]) b1[l] = 6.0;
                 }
                 
                 List probs = exitprobcpp(b1, a1, theta0, t1);
                 double alphaNew = sum(NumericVector(probs[0]));
                 
                 // efficacy boundaries for the secondary trial
                 NumericVector b2 = getBoundcpp(
                   L2, tNew, alphaNew, asfNew, asfparNew, 0, stNew, esNew);
                 
                 return zL2 - theta*sqrt(I2) - b2[L2-1];
               };
    
    lower = brent(f1, interval[0], thetahat, tol);

    auto f2 = [L, zL, I1, k1, t1, r1, b, es, a1, theta0, 
               L2, zL2, I2, tNew, asfNew, asfparNew, stNew, 
               esNew](double theta)->double {
                 
                 // conditional type I error under shifted null
                 double zL1 = -zL + theta*sqrt(I1);
                 NumericVector b1(k1);
                 for (int l=0; l<k1; l++) {
                   b1[l] = (b[l+L] - sqrt(r1[l])*zL1)/sqrt(1 - r1[l]);
                   if (!es[l+L]) b1[l] = 6.0;
                 }
                 
                 List probs = exitprobcpp(b1, a1, theta0, t1);
                 double alphaNew = sum(NumericVector(probs[0]));
                 
                 // efficacy boundaries for the secondary trial
                 NumericVector b2 = getBoundcpp(
                   L2, tNew, alphaNew, asfNew, asfparNew, 0, stNew, esNew);
                 
                 return -zL2 + theta*sqrt(I2) - b2[L2-1];
               };
    
    upper = brent(f2, thetahat, interval[1], tol);

    // repeated p-value is alpha for which the lower bound of theta is zero
    auto f = [kMax, t, asf, asfpar, st, es, 
              L, zL, I1, k1, t1, r1, a1, theta0,
              L2, zL2, I2, tNew, asfNew, asfparNew, stNew, esNew, 
              interval, tol](double aval)->double {
                NumericVector u = getBoundcpp(
                  kMax, t, aval, asf, asfpar, 0, st, es);
                
                auto g = [L, zL, I1, k1, t1, r1, u, es, a1, theta0, 
                          L2, zL2, I2, tNew, asfNew, asfparNew, stNew, 
                          esNew](double theta)->double {
                            
                            // conditional type I error under shifted null
                            double zL1 = zL - theta*sqrt(I1);
                            NumericVector b1(k1);
                            for (int l=0; l<k1; l++) {
                              b1[l] = (u[l+L] - sqrt(r1[l])*zL1)/
                                sqrt(1 - r1[l]);
                              if (!es[l+L]) b1[l] = 6.0;
                            }
                            
                            List probs = exitprobcpp(b1, a1, theta0, t1);
                            double alphaNew = sum(NumericVector(probs[0]));
                            
                            // efficacy boundaries for the secondary trial
                            NumericVector b2 = getBoundcpp(
                              L2, tNew, alphaNew, asfNew, asfparNew, 0, stNew,
                              esNew);
                            
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
//' @param INew The maximum information of the secondary trial.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistic at the interim adaptation look of 
//'   the primary trial.
//' @param theta A scalar or a vector of parameter values of 
//'   length \code{kMax + kMax - L} if \code{MullerSchafer = FALSE} or 
//'   length \code{kMax + kNew} if \code{MullerSchafer = TRUE}.
//' @param IMax The maximum information of the primary trial.
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
//' @param futilityBounds	The lower boundaries on the z-test statistic scale 
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
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and
//'   "none" for no early futility stopping. 
//'   Defaults to "none".
//' @param parameterBetaSpendingNew The parameter value of beta spending 
//'   for the secondary trial. Corresponds to rho for "sfKD", 
//'   and gamma for "sfHSD".
//' @param spendingTimeNew The error spending time of the secondary trial. 
//'   Defaults to missing, in which case, it is the same as 
//'   \code{informationRatesNew}.
//'
//' @return The conditional power given the interim results, parameter 
//' values, and data-dependent design changes.
//' 
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//' 
//' @references 
//' Cyrus R. Mehta and Stuart J. Pocock. 
//' Adaptive increase in sample size when interim results are promising: 
//' A practical guide with examples.
//' Stat Med. 2011;30:3267–3284.
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
//' gam = -log(1-0.05)/12  # hazard for dropout
//' 
//' # The original target number of events was 298 and the new target is 335
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
//' z2 = (-log(hr2))*sqrt(266/4)  # corresponding z-test statistic value
//'  
//' # expected mean of -log(HR) at the original looks and the new final look
//' theta = -log(lr1$HR[c(1,2,3,4)])
//' 
//' # conditional power with sample size increase
//' getCP(INew = (335 - 266)/4, 
//'       L = 2, zL = z2, theta = theta,
//'       IMax = 298/4, kMax = 3, 
//'       informationRates = c(179, 266, 298)/298,
//'       alpha = 0.025, typeAlphaSpending = "sfOF")
//'
//' @export
// [[Rcpp::export]]
double getCP(double INew = NA_REAL,
             const int L = NA_INTEGER,
             const double zL = NA_REAL,
             const NumericVector& theta = NA_REAL,
             const double IMax = NA_REAL,
             const int kMax = NA_INTEGER,
             const NumericVector& informationRates = NA_REAL,
             const LogicalVector& efficacyStopping = NA_LOGICAL,
             const LogicalVector& futilityStopping = NA_LOGICAL,
             const NumericVector& criticalValues = NA_REAL,
             const double alpha = 0.025,
             const String typeAlphaSpending = "sfOF",
             const double parameterAlphaSpending = NA_REAL,
             const NumericVector& userAlphaSpending = NA_REAL,
             const NumericVector& futilityBounds = NA_REAL,
             const String typeBetaSpending = "none",
             const double parameterBetaSpending = NA_REAL,
             const NumericVector& spendingTime = NA_REAL,
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
    c = std::tolower(c);
  });
  
  double asfpar = parameterAlphaSpending;
  
  std::string bsf = typeBetaSpending;
  std::for_each(bsf.begin(), bsf.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double bsfpar = parameterBetaSpending;
  
  std::string asfNew = typeAlphaSpendingNew;
  std::for_each(asfNew.begin(), asfNew.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double asfparNew = parameterAlphaSpendingNew;
  
  std::string bsfNew = typeBetaSpendingNew;
  std::for_each(bsfNew.begin(), bsfNew.end(), [](char & c) {
    c = std::tolower(c);
  });
  
  double bsfparNew = parameterBetaSpendingNew;
  
  if (R_isnancpp(INew)) {
    stop("INew must be provided");
  }
  
  if (INew <= 0) {
    stop("INew must be positive");
  }
  
  if (R_isnancpp(L)) {
    stop("L must be provided");
  }
  
  if (L <= 0) {
    stop("L must be a positive integer");
  }
  
  if (R_isnancpp(zL)) {
    stop("zL must be provided");
  }
  
  if (is_true(any(is_na(theta)))) {
    stop("theta must be provided");
  }
  
  if (R_isnancpp(IMax)) {
    stop("IMax must be provided");
  }
  
  if (IMax <= 0) {
    stop("IMax must be positive");
  }
  
  if (R_isnancpp(kMax)) {
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
    t = as<NumericVector>(tem)/(kMax+0.0);
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
    
    NumericVector u(kMax), l(kMax, -6.0), theta0(kMax);
    for (int i=0; i<kMax; i++) {
      u[i] = criticalValues[i];
      if (!es[i]) u[i] = 6.0;
    }
    
    List probs = exitprobcpp(u, l, theta0, t);
    alpha1 = sum(NumericVector(probs[0]));
  }
  
  if (!R_isnancpp(alpha1)) {
    if (alpha1 < 0.00001 || alpha1 >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
  } else {
    stop("alpha must be provided for missing criticalValues");
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
      stop("Elements of userAlphaSpending must be nonnegnative");
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
    
    if (!(bsfNew=="sfof" || bsfNew=="sfp" || bsfNew=="sfkd" || 
        bsfNew=="sfhsd" || bsfNew=="none")) {
      stop("Invalid value for typeBetaSpendingNew");
    }
    
    if ((bsfNew=="sfkd" || bsfNew=="sfhsd") && R_isnancpp(bsfparNew)) {
      stop("Missing value for parameterBetaSpendingNew");
    }
    
    if (bsfNew=="sfkd" && bsfparNew <= 0) {
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
      stNew = clone(tNew);
    }
  }
  
  
  // obtain critical values for the primary trial
  if (is_true(any(is_na(criticalValues)))) {
    b = getBoundcpp(kMax, t, alpha1, asf, asfpar, userAlphaSpending, st, es);
  }
  
  // obtain futility bounds for the primary trial
  if (kMax > 1) {
    if (is_true(any(is_na(futilityBounds))) && bsf=="none") {
      a = rep(-6.0, kMax);
      a[kMax-1] = b[kMax-1];
    } else if (is_false(any(is_na(futilityBounds))) && 
      a.size() == kMax-1) {
      a.push_back(b[kMax-1]);
    }
  } else {
    if (is_true(any(is_na(futilityBounds)))) {
      a = b[kMax-1];
    }
  }
  
  if (is_true(any(is_na(a)))) {
    NumericVector theta1(kMax);
    if (theta.size() == 1) {
      theta1.fill(theta[0]);
    } else if (theta.size() > kMax) {
      IntegerVector idx = Range(0, kMax-1);
      theta1 = theta[idx];
    } else {
      stop("Invalid length for theta");
    }
    
    List out = getPower(alpha1, kMax, b, theta1, IMax*t, bsf, bsfpar, st, fs);
    a = out[1];
  }
  

  int k1 = kMax - L;
  NumericVector t1(k1), r1(k1), b1(k1), a1(k1, -6.0), zero(k1);
  for (int l=0; l<k1; l++) {
    t1[l] = (t[l+L] - t[L-1])/(1 - t[L-1]);
    r1[l] = t[L-1]/t[l+L];
    b1[l] = (b[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
    if (!es[l+L]) b1[l] = 6.0;
  }
  
  double result;
  if (!MullerSchafer) {
    NumericVector theta1(k1+1);
    if (theta.size() == 1) {
      theta1.fill(theta[0]);
    } else if (theta.size() == kMax+k1){
      theta1[0] = theta[L-1];
      for (int l=0; l<k1; l++) {
        theta1[l+1] = theta[kMax+l];
      }
    } else {
      stop("Invalid length for theta"); 
    }
    
    for (int l=0; l<k1; l++) {
      a1[l] = (a[l+L] - sqrt(r1[l])*zL)/sqrt(1 - r1[l]);
      if (!fs[l+L]) a1[l] = -6.0;
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
    } else if (theta.size() == kMax+kNew) {
      theta1[0] = theta[L-1];
      for (int l=0; l<kNew; l++) {
        theta1[l+1] = theta[kMax+l];
      }
    } else {
      stop("Invalid length for theta");
    }
    
    // obtain conditional type I error
    List probs = exitprobcpp(b1, a1, zero, t1);
    double alphaNew = sum(NumericVector(probs[0]));
    
    // obtain efficacy boundaries for the secondary trial
    NumericVector b2 = getBoundcpp(
      kNew, tNew, alphaNew, asfNew, asfparNew, 0, stNew, esNew);
    
    // obtain conditional power
    NumericVector mu(kNew), I2(kNew);
    for (int l=0; l<kNew; l++) {
      double r = IMax*t[L-1]/(IMax*t[L-1] + INew*tNew[l]);
      mu[l] = (theta1[l+1] - r*theta1[0])/(1 - r);
      I2[l] = INew*tNew[l];
    }
    
    List out = getPower(
      alphaNew, kNew, b2, mu, I2, bsfNew, bsfparNew, stNew, fsNew);
    
    result = out[0];
  }
  
  return result;
}
