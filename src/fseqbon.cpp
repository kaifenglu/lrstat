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
  
  if (!(asf=="of" || asf=="p" || asf=="wt" || 
      asf=="sfof" || asf=="sfp" || asf=="sfkd" ||
      asf=="sfhsd" || asf=="none")) {
    stop("Invalid type for typeAlphaSpending");
  }
  
  if (asf=="wt" && R_isnancpp(asfpar)) {
    stop("Missing parameter for WT");
  }
  
  if (asf=="sfkd" && R_isnancpp(asfpar)) {
    stop("Missing parameter for sfKD");
  }
  
  if (asf=="sfhsd" && R_isnancpp(asfpar)) {
    stop("Missing parameter for sfHSD");
  }
  
  if (asf=="sfkd" && asfpar <= 0) {
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
      LogicalVector x(i+1);
      x.fill(1);
      
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
      stop("Invalid type for typeAlphaSpending");
    }
    
    if (asfi=="wt" && R_isnancpp(asfpar(i))) {
      stop("Missing parameter for WT");
    }
    
    if (asfi=="sfkd" && R_isnancpp(asfpar(i))) {
      stop("Missing parameter for sfKD");
    }
    
    if (asfi=="sfhsd" && R_isnancpp(asfpar(i))) {
      stop("Missing parameter for sfHSD");
    }
    
    if (asfi=="sfkd" && asfpar(i) <= 0) {
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
              
              LogicalVector x(k+1);
              x.fill(1);
              
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

