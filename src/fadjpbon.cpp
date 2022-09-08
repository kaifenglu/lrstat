#include <Rcpp.h>

using namespace Rcpp;


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




