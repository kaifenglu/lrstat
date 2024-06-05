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
//' @return A list containing the new vector of weights, the new
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
List updateGraph(const NumericVector& w = NA_REAL,
                 const NumericMatrix& G = NA_REAL,
                 const IntegerVector& I = NA_INTEGER,
                 const int j = NA_INTEGER) {
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
//' g = matrix(c(0,0,1,0, 0,0,0,1, 0,1,0,0, 1,0,0,0),
//'            nrow=4, ncol=4, byrow=TRUE)
//' (wgtmat = fwgtmat(w,g))
//'
//' @export
// [[Rcpp::export]]
NumericMatrix fwgtmat(const NumericVector& w = NA_REAL,
                      const NumericMatrix& G = NA_REAL) {
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
          if (cc(k) && cc(l) && (k != l) && (g(k,j)*g(j,k) < 1.0 - 1e-12)) {
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
    const std::string typeAlphaSpending = "sfOF",
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
  std::string asf1;
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

              asf1 = as<std::string>(asf(j));
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
                          const std::string test = "hochberg",
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


// [[Rcpp::export]]
NumericMatrix fstdmixcpp(const NumericMatrix& p,
                         const LogicalMatrix& family,
                         const LogicalMatrix& serial,
                         const LogicalMatrix& parallel,
                         const NumericVector& gamma,
                         const std::string test = "hommel",
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
                         const std::string test = "hommel",
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
//' @return A data frame with the following components:
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
DataFrame getCI(const int L = NA_INTEGER,
                const double zL = NA_REAL,
                const double IMax = NA_REAL,
                const NumericVector& informationRates = NA_REAL,
                const LogicalVector& efficacyStopping = NA_LOGICAL,
                const NumericVector& criticalValues = NA_REAL,
                const double alpha = 0.025,
                const std::string typeAlphaSpending = "sfOF",
                const double parameterAlphaSpending = NA_REAL,
                const NumericVector& spendingTime = NA_REAL) {

  NumericVector t = clone(informationRates);
  LogicalVector es = clone(efficacyStopping);
  NumericVector b = clone(criticalValues);
  NumericVector st = clone(spendingTime);

  if (L == NA_INTEGER) {
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
  double tol = 1.0e-6;

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

  DataFrame result = DataFrame::create(
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
//' @return A data frame with the following components:
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
DataFrame getRCI(const int L = NA_INTEGER,
                 const double zL = NA_REAL,
                 const double IMax = NA_REAL,
                 const NumericVector& informationRates = NA_REAL,
                 const LogicalVector& efficacyStopping = NA_LOGICAL,
                 const NumericVector& criticalValues = NA_REAL,
                 const double alpha = 0.025,
                 const std::string typeAlphaSpending = "sfOF",
                 const double parameterAlphaSpending = NA_REAL,
                 const NumericVector& spendingTime = NA_REAL) {

  NumericVector t = clone(informationRates);
  LogicalVector es = clone(efficacyStopping);
  NumericVector b = clone(criticalValues);
  NumericVector st = clone(spendingTime);

  if (L == NA_INTEGER) {
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

  DataFrame result = DataFrame::create(
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
  IntegerVector js = findInterval3(astars, p);
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

    z1j = brent(f, -6.0, 6.0, 1.0e-6);
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
//' @param typeAlphaSpending The type of alpha spending for the primary
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
//' @return A data frame with the following variables:
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
//' n2 = ceiling(des2$secondaryTrial$overallResults$information*4*20^2)
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
DataFrame getADCI(const int L = NA_INTEGER,
                  const double zL = NA_REAL,
                  const double IMax = NA_REAL,
                  const int kMax = NA_INTEGER,
                  const NumericVector& informationRates = NA_REAL,
                  const LogicalVector& efficacyStopping = NA_LOGICAL,
                  const NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.25,
                  const std::string typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const NumericVector& spendingTime = NA_REAL,
                  const int L2 = NA_INTEGER,
                  const double zL2 = NA_REAL,
                  const double INew = NA_REAL,
                  const bool MullerSchafer = 0,
                  const NumericVector& informationRatesNew = NA_REAL,
                  const LogicalVector& efficacyStoppingNew = NA_LOGICAL,
                  const std::string typeAlphaSpendingNew = "sfOF",
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

  if (L == NA_INTEGER) {
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

  if (L2 == NA_INTEGER) {
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
  double pvalue = f_bwpvalue(0,K,L,zL,b,I,L2,zL2,b2,I2);

  NumericVector interval(2);
  interval[0] = (zL - b[L-1])/sqrt(I[L-1]);
  interval[1] = (zL + b[L-1])/sqrt(I[L-1]);
  double tol = 1.0e-6;

  auto f = [K,L,zL,b,I,L2,zL2,b2,I2](double theta)->double {
    return f_bwpvalue(theta,K,L,zL,b,I,L2,zL2,b2,I2) - 0.5;
  };
  double thetahat = brent(f, interval[0], interval[1], tol);

  auto f1 = [K,L,zL,b,I,L2,zL2,b2,I2,cilevel](double theta)->double {
    return f_bwpvalue(theta,K,L,zL,b,I,L2,zL2,b2,I2) - (1-cilevel)/2;
  };
  double lower = brent(f1, interval[0], thetahat, tol);

  auto f2 = [K,L,zL,b,I,L2,zL2,b2,I2,cilevel](double theta)->double {
    return f_bwpvalue(theta,K,L,zL,b,I,L2,zL2,b2,I2) - (1+cilevel)/2;
  };
  double upper = brent(f2, thetahat, interval[1], tol);

  DataFrame result = DataFrame::create(
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
//' @param typeAlphaSpending The type of alpha spending for the primary
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
//' @return A data frame with the following variables:
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
//' n2 = ceiling(des2$secondaryTrial$overallResults$information*4*20^2)
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
DataFrame getADRCI(const int L = NA_INTEGER,
                   const double zL = NA_REAL,
                   const double IMax = NA_REAL,
                   const int kMax = NA_INTEGER,
                   const NumericVector& informationRates = NA_REAL,
                   const LogicalVector& efficacyStopping = NA_LOGICAL,
                   const NumericVector& criticalValues = NA_REAL,
                   const double alpha = 0.025,
                   const std::string typeAlphaSpending = "sfOF",
                   const double parameterAlphaSpending = NA_REAL,
                   const NumericVector& spendingTime = NA_REAL,
                   const int L2 = NA_INTEGER,
                   const double zL2 = NA_REAL,
                   const double INew = NA_REAL,
                   const bool MullerSchafer = 0,
                   const NumericVector& informationRatesNew = NA_REAL,
                   const LogicalVector& efficacyStoppingNew = NA_LOGICAL,
                   const std::string typeAlphaSpendingNew = "sfOF",
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

  if (L == NA_INTEGER) {
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

  if (L2 == NA_INTEGER) {
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
    NumericVector u = getBoundcpp(J, t[i], 0.5, asf, asfpar, 0,
                                  st[i], es[i]);
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
    double tol = 1.0e-6;

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
                              L2, tNew, alphaNew, asfNew, asfparNew, 0,
                              stNew, esNew);

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

  DataFrame result = DataFrame::create(
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
//' data-dependent changes in the error spending function, as well as the
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
//' @param varianceRatio The ratio of the variance under H0 to the variance
//'   under H1.
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
//' # Two interim analyses have occurred with 179 and 266 events,
//' # respectively. The observed hazard ratio at the second interim
//' # look is 0.81.
//'
//' trialsdt = as.Date("2020-03-04")                       # trial start date
//' iadt = c(as.Date("2022-02-01"), as.Date("2022-11-01")) # interim dates
//' mo1 = as.numeric(iadt - trialsdt + 1)/30.4375          # interim months
//'
//' # Assume a piecewise Poisson enrollment process with a 8-month ramp-up
//' # and 521 patients were enrolled after 17.94 months
//' N = 521                   # total number of patients
//' Ta = 17.94                # enrollment duration
//' Ta1 = 8                   # assumed end of enrollment ramp-up
//' enrate = N / (Ta - Ta1/2) # enrollment rate after ramp-up
//'
//' # Assume a median survival of 16.7 months for the control group, a
//' # 5-month delay in treatment effect, and a hazard ratio of 0.7 after
//' # the delay
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

  if (L == NA_INTEGER) {
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

  if (varianceRatio <= 0) {
    stop("varianceRatio must be positive");
  }

  NumericVector w = rep(sqrt(varianceRatio), kMax);


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

    List out = getPower(alpha1, kMax, b, theta1, IMax*t, bsf, bsfpar,
                        st, fs, w);
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
      alphaNew, kNew, b2, mu, I2, bsfNew, bsfparNew, stNew, fsNew, w);

    result = out[0];
  }

  return result;
}


// [[Rcpp::export]]
NumericMatrix ftrunccpp(const NumericMatrix& p,
                        const std::string test,
                        const double gamma) {
  std::string test1 = test;
  std::for_each(test1.begin(), test1.end(), [](char & c) {
    c = std::tolower(c);
  });

  int niters = p.nrow();
  int m = p.ncol();
  int ntests = pow(2, m) - 1;
  int i, j, k, iter;
  LogicalMatrix incid(ntests, m);
  NumericMatrix pinter(niters, ntests);
  NumericMatrix padj(niters, m);

  for (i=0; i<ntests; i++) {
    int number = ntests - i;

    // binary representation of elementary hypotheses in the intersection
    LogicalVector cc(m);
    for (j=0; j<m; j++) {
      cc(j) = (number/(1 << (m - 1 - j))) % 2;
    }

    incid(i, _) = cc;

    // hypotheses included in the intersection
    k = sum(cc);
    IntegerVector hyp = which(cc);

    for (iter=0; iter<niters; iter++) {
      NumericVector p1(k);
      for (j=0; j<k; j++) {
        p1(j) = p(iter, hyp(j));
      }

      NumericVector p2 = stl_sort(p1);

      double q = 1.0, w;
      for (j=0; j<k; j++) {
        if (test == "hommel") {
          w = (j+1)*gamma/k + (1-gamma)/m;
        } else if (test == "hochberg") {
          w = gamma/(k - j) + (1-gamma)/m;
        } else { // holm
          w = gamma/k + (1-gamma)/m;
        }

        q = std::min(q, p2(j)/w);
      }

      pinter(iter, i) = q;
    }
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



// probability of failure to reject H0 for Simon's two-stage design
double pfutile(double p, int n1, int n2, int r1, int r) {
  double aval = R::pbinom(r1, n1, p, 1, 0);
  for (int x=r1+1; x<=std::min(n1,r); x++) {
    aval += R::dbinom(x, n1, p, 0)*R::pbinom(r-x, n2, p, 1, 0);
  }
  return aval;
}


//' @title Simon's two-stage design
//' @description Obtains Simon's two-stage minimax, admissible, and
//' optimal designs.
//'
//' @param alpha Type I error rate (one-sided).
//' @param beta Type II error rate (1-power).
//' @param piH0 Response probability under the null hypothesis.
//' @param pi Response probability under the alternative hypothesis.
//' @param n_max Upper limit for sample size, defaults to 110.
//'
//' @return A data frame containing the following variables:
//'
//' * \code{piH0}: Response probability under the null hypothesis.
//'
//' * \code{pi}: Response probability under the alternative hypothesis.
//'
//' * \code{alpha}: The specified one-sided significance level.
//'
//' * \code{beta}: The specified type II error.
//'
//' * \code{n}: Total sample size.
//'
//' * \code{n1}: Stage 1 sample size.
//'
//' * \code{r1}: Futility boundary for stage 1.
//'
//' * \code{r}: Futility boundary for stage 2.
//'
//' * \code{EN0}: Expected sample size under the null hypothesis.
//'
//' * \code{attainedAlpha}: Attained type 1 error.
//'
//' * \code{power}: Attained power.
//'
//' * \code{PET0}: Probability of early stopping under the null hypothesis.
//'
//' * \code{w_lower}: Lower bound of the interval for \code{w}.
//'
//' * \code{w_upper}: Upper bound of the interval for \code{w}.
//'
//' * \code{design}: Description of the design, e.g., minimax, admissible,
//'   or optimal.
//'
//' Here \code{w} is the weight in the objective function:
//' \code{w*n + (1-w)*EN0}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' simon2stage(0.05, 0.15, 0.1, 0.3)
//'
//' @export
// [[Rcpp::export]]
DataFrame simon2stage(const double alpha = NA_REAL,
                      const double beta = NA_REAL,
                      const double piH0 = NA_REAL,
                      const double pi = NA_REAL,
                      const int n_max = 110) {

  int n1, n, r1, r;

  if (R_isnancpp(alpha)) {
    stop("alpha must be provided");
  }

  if (alpha < 0.00001 || alpha >= 1) {
    stop("alpha must lie in [0.00001, 1)");
  }

  if (R_isnancpp(beta)) {
    stop("beta must be provided");
  }

  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }

  if (R_isnancpp(piH0)) {
    stop("piH0 must be provided");
  }

  if (piH0 <= 0 || piH0 >= 1) {
    stop("piH0 must lie between 0 and 1");
  }

  if (R_isnancpp(pi)) {
    stop("pi must be provided");
  }

  if (pi <= piH0 || pi >= 1) {
    stop("pi must lie between piH0 and 1");
  }

  double p = (piH0 + pi)/2;
  double z1 = R::qnorm(1-alpha, 0, 1, 1, 0);
  double z2 = R::qnorm(1-beta, 0, 1, 1, 0);
  int n_min = std::floor(p*(1-p)*pow((z1 + z2)/(pi - piH0), 2));
  int n_max1 = std::ceil(2.0*n_min);

  int n_lower = std::floor(0.5*n_min);
  int n_upper = std::min(n_max, n_max1);
  int m = n_upper - n_lower + 1;

  double alphastar, betastar;

  IntegerVector nx = IntegerVector(m, NA_INTEGER);
  IntegerVector n1x = IntegerVector(m, NA_INTEGER);
  IntegerVector r1x = IntegerVector(m, NA_INTEGER);
  IntegerVector rx = IntegerVector(m, NA_INTEGER);
  NumericVector en0x = NumericVector(m, NA_REAL);
  NumericVector pet0x = NumericVector(m, NA_REAL);
  NumericVector alphax = NumericVector(m, NA_REAL);
  NumericVector powerx = NumericVector(m, NA_REAL);

  int i = 0;
  for (n=n_lower; n<=n_upper; n++) {
    // identify the design satisfying the type 1 and type 2 error constraints
    // while having the minimum expected sample size under H0 for the given n
    double EN0 = n;
    List optimal;

    bool exist = 0; // initialize to FALSE
    for (n1=1; n1<=n-1; n1++) {
      int n2=n-n1;
      for (r1=0; r1<=n1; r1++) {
        // find the maximum value of r satisfying the type 2 error constraint
        for (r=r1+n2; r>=r1; r--) {
          betastar = pfutile(pi, n1, n2, r1, r);
          if (betastar <= beta) break;
        }

        if (r>=r1) {
          // check the type 1 error constraint
          alphastar = 1 - pfutile(piH0, n1, n2, r1, r);
          if (alphastar <= alpha) {
            exist = 1;
            // calculate the expected sample size under H0
            double pet0 = R::pbinom(r1, n1, piH0, 1, 0);
            double en0 = n1 + (1-pet0)*n2;

            if (en0 < EN0) {
              EN0 = en0;
              optimal = List::create(
                _["n"] = n,
                _["n1"] = n1,
                _["r1"] = r1,
                _["r"] = r,
                _["en0"] = en0,
                _["pet0"] = pet0,
                _["alpha"] = alphastar,
                _["power"] = 1-betastar
              );
            }
          }
        }
      }
    }

    // if there is a design satisfying the type 1 and type 2 constraints
    if (exist) {
      nx[i] = optimal["n"];
      n1x[i] = optimal["n1"];
      r1x[i] = optimal["r1"];
      rx[i] = optimal["r"];
      en0x[i] = optimal["en0"];
      pet0x[i] = optimal["pet0"];
      alphax[i] = optimal["alpha"];
      powerx[i] = optimal["power"];
      i++;
    }
  }

  // only keep nonmissing records
  LogicalVector sub = !is_na(nx);
  en0x = en0x[sub];

  // identify the optimal design
  int I = which_min(en0x);
  IntegerVector q = Range(0, I);

  // only keep candidate designs between minimax and optimal designs
  nx = nx[q];
  n1x = n1x[q];
  r1x = r1x[q];
  rx = rx[q];
  en0x = en0x[q];
  pet0x = pet0x[q];
  alphax = alphax[q];
  powerx = powerx[q];

  // only keep the admissible designs on the convex hull
  IntegerVector u = IntegerVector::create(0);

  // gift wrapping algorithm (Jarvis March)
  i = 0;
  while (i < I) {
    IntegerVector v = Range(i+1, I);
    NumericVector slope = (as<NumericVector>(en0x[v]) - en0x[i])/
      (as<NumericVector>(nx[v]) - nx[i]);
    int j = which_min(slope);
    i += j + 1; // update the current point
    u.push_back(i);
  }

  nx = nx[u];
  n1x = n1x[u];
  r1x = r1x[u];
  rx = rx[u];
  en0x = en0x[u];
  pet0x = pet0x[u];
  alphax = alphax[u];
  powerx = powerx[u];

  // derive the interval for w
  m = u.size();
  NumericVector w1(m), w2(m);
  for (i=0; i<m; i++) {
    if (i<m-1) {
      double slope = (en0x[i+1] - en0x[i])/(nx[i+1] - nx[i]);
      double w = slope/(slope - 1);

      if (i==0) { // minimax design
        w1[i] = w;
        w2[i] = 1;
      } else { // admissible designs
        w1[i] = w;
        w2[i] = w1[i-1];
      }
    } else {  // optimal design
      w1[i] = 0;
      w2[i] = w1[i-1];
    }
  }

  CharacterVector design(m);
  for (i=0; i<m; i++) {
    if (w1[i] == 0) {
      design[i] = "Optimal";
    } else if (w2[i] == 1) {
      design[i] = "Minimax";
    } else {
      design[i] = "Admissible";
    }
  }

  DataFrame results = DataFrame::create(
    _["piH0"] = rep(piH0, m),
    _["pi"] = rep(pi, m),
    _["alpha"] = rep(alpha, m),
    _["beta"] = rep(beta, m),
    _["n"] = nx,
    _["n1"] = n1x,
    _["r1"] = r1x,
    _["r"] = rx,
    _["EN0"] = en0x,
    _["attainedAlpha"] = alphax,
    _["attainedPower"] = powerx,
    _["PET0"] = pet0x,
    _["w_lower"] = w1,
    _["w_upper"] = w2,
    _["design"] = design);

  return results;
}


//' @title Power for binomial one-sample exact test
//' @description Obtains the power for binomial one-sample exact test.
//'
//' @param n The sample size.
//' @param piH0 The response probability under the null hypothesis.
//' @param pi The response probability under the alternative hypothesis.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//'
//' @return A data frame containing the critical value of the number of
//' responses for rejecting the null hypothesis, the attained type I
//' error, the power for the exact test, the sample size, the
//' response probabilities under the null and alternative hypotheses,
//' and the direction of the alternative.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' powerOnePropExact(n = 110, piH0 = 0.15, pi = 0.25, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
DataFrame powerOnePropExact(const int n = NA_INTEGER,
                            const double piH0 = NA_REAL,
                            const double pi = NA_REAL,
                            const double alpha = 0.025) {

  // find the critical value r to meet the type I error condition
  int r;
  double attainedAlpha = 0, power;
  bool directionUpper = pi > piH0;
  if (directionUpper) {
    r = R::qbinom(1-alpha, n, piH0, 1, 0) + 1;
    attainedAlpha = R::pbinom(r-1, n, piH0, 0, 0);
    power = R::pbinom(r-1, n, pi, 0, 0);
  } else {
    int rstar = R::qbinom(alpha, n, piH0, 1, 0);
    r = R::pbinom(rstar, n, piH0, 1, 0) <= alpha ? rstar : rstar - 1;
    attainedAlpha = R::pbinom(r, n, piH0, 1, 0);
    power = R::pbinom(r, n, pi, 1, 0);
  }

  DataFrame result = DataFrame::create(
    _["alpha"] = alpha,
    _["attainedAlpha"] = attainedAlpha,
    _["power"] = power,
    _["n"] = n,
    _["piH0"] = piH0,
    _["pi"] = pi,
    _["r"] = r);

  return result;
}


//' @title Sample size for binomial one-sample exact test
//' @description Obtains the sample size for binomial one-sample exact test.
//'
//' @param beta The type II error.
//' @param piH0 The response probability under the null hypothesis.
//' @param pi The response probability under the alternative hypothesis.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//'
//' @return A data frame containing the critical value of the number of
//' responses for rejecting the null hypothesis, the attained type I
//' error, the power for the exact test, the sample size, the
//' response probabilities under the null and alternative hypotheses,
//' and the direction of the alternative.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' samplesizeOnePropExact(beta = 0.2, piH0 = 0.15, pi = 0.25, alpha = 0.025)
//'
//' @export
// [[Rcpp::export]]
DataFrame samplesizeOnePropExact(const double beta = 0.2,
                                 const double piH0 = NA_REAL,
                                 const double pi = NA_REAL,
                                 const double alpha = 0.025) {

  double z1 = R::qnorm(1-alpha, 0, 1, 1, 0);
  double z2 = R::qnorm(1-beta, 0, 1, 1, 0);
  double tau = pow((z1 + z2)/(pi - piH0), 2);
  int n_lower = std::floor(0.5*tau*std::min(piH0*(1-piH0), pi*(1-pi)));
  int n_upper = std::ceil(2*tau*std::max(piH0*(1-piH0), pi*(1-pi)));

  int n;
  DataFrame a, b;
  for (n=n_lower; n<=n_upper; n++) {
    a = powerOnePropExact(n, piH0, pi, alpha);

    // ensure that the power exceeds 1-beta for the next 10 sample sizes
    if (as<double>(a["power"]) >= 1-beta) {
      int i;
      bool okay = 1;
      for (i=1; i<=10; i++) {
        b = powerOnePropExact(n+i, piH0, pi, alpha);
        if (as<double>(b["power"]) < 1-beta) {
          okay = 0;
          break;
        }
      }

      if (okay) {
        break;
      } else {
        n += i;
      }
    }
  }

  return a;
}


//' @title Power for Poisson one-sample exact test
//' @description Obtains the power for Poisson one-sample exact test.
//'
//' @param n The sample size.
//' @param lambdaH0 The Poisson rate under the null hypothesis.
//' @param lambda The Poisson rate under the alternative hypothesis.
//' @param D The average exposure per subject.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//'
//' @return A data frame containing the critical value of the number of
//' events for rejecting the null hypothesis, the attained type I
//' error, the power for the exact test, the sample size,
//' the Poisson rates under the null and alternative hypotheses,
//' the average exposure, and the direction of the alternative.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' powerOneRateExact(n = 525, lambdaH0 = 0.049, lambda = 0.012,
//'                   D = 0.5, alpha = 0.025)
//'
//' @export
// [[Rcpp::export]]
DataFrame powerOneRateExact(const int n = NA_INTEGER,
                            const double lambdaH0 = NA_REAL,
                            const double lambda = NA_REAL,
                            const double D = 1,
                            const double alpha = 0.025) {

  // find the critical value r to meet the type I error condition
  int r;
  double attainedAlpha = 0, power;
  if (lambda > lambdaH0) {
    r = R::qpois(1-alpha, n*lambdaH0*D, 1, 0) + 1;
    attainedAlpha = R::ppois(r-1, n*lambdaH0*D, 0, 0);
    power = R::ppois(r-1, n*lambda*D, 0, 0);
  } else {
    int rstar = R::qpois(alpha, n*lambdaH0*D, 1, 0);
    r = R::ppois(rstar, n*lambdaH0*D, 1, 0) <= alpha ? rstar : rstar - 1;
    attainedAlpha = R::ppois(r, n*lambdaH0*D, 1, 0);
    power = R::ppois(r, n*lambda*D, 1, 0);
  }

  DataFrame result = DataFrame::create(
    _["alpha"] = alpha,
    _["attainedAlpha"] = attainedAlpha,
    _["power"] = power,
    _["n"] = n,
    _["lambdaH0"] = lambdaH0,
    _["lambda"] = lambda,
    _["D"] = D,
    _["r"] = r);

  return result;
}


//' @title Sample size for Poisson one-sample exact test
//' @description Obtains the sample size for Poisson one-sample exact test.
//'
//' @param beta The type II error.
//' @param lambdaH0 The Poisson rate under the null hypothesis.
//' @param lambda The Poisson rate under the alternative hypothesis.
//' @param D The average exposure per subject.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//'
//' @return A data frame containing the critical value of the number of
//' events for rejecting the null hypothesis, the attained type I
//' error, the power for the exact test, the sample size,
//' the Poisson rates under the null and alternative hypotheses,
//' the average exposure, and the direction of the alternative.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' samplesizeOneRateExact(beta = 0.2, lambdaH0 = 0.2, lambda = 0.3,
//'                        D = 1, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
DataFrame samplesizeOneRateExact(const double beta = 0.2,
                                 const double lambdaH0 = NA_REAL,
                                 const double lambda = NA_REAL,
                                 const double D = 1,
                                 const double alpha = 0.025) {

  double z1 = R::qnorm(1-alpha, 0, 1, 1, 0);
  double z2 = R::qnorm(1-beta, 0, 1, 1, 0);
  double tau = pow((z1 + z2)/log(lambda/lambdaH0), 2);
  int n_lower = std::floor(0.5*tau/(std::max(lambdaH0, lambda)*D));
  int n_upper = std::ceil(2*tau/(std::min(lambdaH0, lambda)*D));

  int n;
  DataFrame a, b;
  for (n=n_lower; n<=n_upper; n++) {
    a = powerOneRateExact(n, lambdaH0, lambda, D, alpha);

    // ensure that the power exceeds 1-beta for the next 10 sample sizes
    if (as<double>(a["power"]) >= 1-beta) {
      int i;
      bool okay = 1;
      for (i=1; i<=10; i++) {
        b = powerOneRateExact(n+i, lambdaH0, lambda, D, alpha);
        if (as<double>(b["power"]) < 1-beta) {
          okay = 0;
          break;
        }
      }

      if (okay) {
        break;
      } else {
        n += i;
      }
    }
  }

  return a;
}


//' @title Power for Fisher's exact test for two proportions
//' @description Obtains the power given sample size for Fisher's exact
//' test for two proportions.
//'
//' @param n The total sample size.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The two-sided significance level. Defaults to 0.05.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The two-sided significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @keywords internal
//'
//' @examples
//'
//' (design1 <- powerFisherExact(
//'   n = 136, pi1 = 0.25, pi2 = 0.05, alpha = 0.05))
//'
//' @export
// [[Rcpp::export]]
DataFrame powerFisherExact(const int n = NA_INTEGER,
                           const double pi1 = NA_REAL,
                           const double pi2 = NA_REAL,
                           const double allocationRatioPlanned = 1,
                           const double alpha = 0.05) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  int n1 = std::round(n*r);
  int n2 = n - n1;
  int i, j;

  double power = 0;
  for (int m=0; m<=n; m++) {
    int lower = std::max(0, m-n2);
    int upper = std::min(n1, m);
    int k = upper - lower + 1;
    IntegerVector x = seq(lower, upper);
    IntegerVector y = m - x;

    // conditional density under the null
    NumericVector p0 = dhyper(x, n1, n2, m, 0);

    // sort p0 in ascending order
    IntegerVector order = seq(0, k-1);
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return p0[i] < p0[j];
    });

    NumericVector p0sorted = p0[order];

    // tally the cumulative probabilities for unique values of p0
    IntegerVector idx(1,0);
    for (i=1; i<k; i++) {
      if (p0sorted[i] != p0sorted[i-1]) {
        idx.push_back(i);
      }
    }
    int k1 = idx.size();

    idx.push_back(k);

    double s = 0;
    NumericVector cp0(k1);
    for (i=0; i<k1; i++) {
      for (j=idx[i]; j<idx[i+1]; j++) {
        s += p0sorted[j];
      }
      cp0[i] = s;
    }

    // reject H0 only if the smallest table probability is below alpha
    if (cp0[0] <= alpha) {
      // find the critical value table
      for (i=1; i<k1; i++) {
        if (cp0[i] > alpha) break;
      }

      // unconditional density under the alternative
      NumericVector p1 = dbinom(x, n1, pi1, 0)*dbinom(y, n2, pi2, 0);
      NumericVector p1sorted = p1[order];

      for (j=0; j<idx[i]; j++) {
        power += p1sorted[j];
      }
    }
  }

  DataFrame result = DataFrame::create(
    _["alpha"] = alpha,
    _["power"] = power,
    _["n"] = n,
    _["pi1"] = pi1,
    _["pi2"] = pi2,
    _["allocationRatioPlanned"] = allocationRatioPlanned);

  return result;
}


//' @title Sample size for Fisher's exact test for two proportions
//' @description Obtains the sample size given power for Fisher's exact
//' test for two proportions.
//'
//' @param beta The type II error.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The two-sided significance level. Defaults to 0.05.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The two-sided significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @keywords internal
//'
//' @examples
//'
//' (design1 <- samplesizeFisherExact(
//'   beta = 0.1, pi1 = 0.25, pi2 = 0.05, alpha = 0.05))
//'
//' @export
// [[Rcpp::export]]
DataFrame samplesizeFisherExact(const double beta = NA_REAL,
                                const double pi1 = NA_REAL,
                                const double pi2 = NA_REAL,
                                const double allocationRatioPlanned = 1,
                                const double alpha = 0.05) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  double v1 = 1/(4*r*(1-r));

  auto f = [pi1, pi2, r, v1, alpha, beta](double n)->double {
    double delta = asin(sqrt(pi1 + 1/(2*n*r)*((pi1 < pi2) - (pi1 > pi2)))) -
      asin(sqrt(pi2 + 1/(2*n*(1-r))*((pi2 < pi1) - (pi2 > pi1))));
    return R::pnorm(fabs(delta)*sqrt(n/v1) - R::qnorm(1-alpha, 0, 1, 1, 0),
                    0, 1, 1, 0) - (1-beta);
  };

  double delta = asin(sqrt(pi1)) - asin(sqrt(pi2));
  double n0 = pow(R::qnorm(1-alpha, 0, 1, 1, 0) +
                  R::qnorm(1-beta, 0, 1, 1, 0), 2)*v1/pow(delta, 2);
  double n1 = brent(f, 0.5*n0, 2*n0, 1.0e-6);

  int n_lower = std::floor(n1), n_upper = std::ceil(10*n1), n;

  DataFrame a, b;
  a = powerFisherExact(n_lower, pi1, pi2, allocationRatioPlanned, alpha);
  while (as<double>(a["power"]) >= 1-beta) {
    a = powerFisherExact(--n_lower, pi1, pi2, allocationRatioPlanned, alpha);
  }

  for (n=n_lower+1; n<=n_upper; n++) {
    a = powerFisherExact(n, pi1, pi2, allocationRatioPlanned, alpha);

    // ensure that the power exceeds 1-beta for the next 10 sample sizes
    if (as<double>(a["power"]) >= 1-beta) {
      int i;
      bool okay = 1;
      for (i=1; i<=10; i++) {
        b = powerFisherExact(n+i, pi1, pi2, allocationRatioPlanned, alpha);
        if (as<double>(b["power"]) < 1-beta) {
          okay = 0;
          break;
        }
      }

      if (okay) {
        break;
      } else {
        n += i;
      }
    }

  }

  return a;
}


//' @title REML estimates of individual proportions with specified risk
//' difference
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified risk difference.
//'
//' @param riskDiffH0 The specified risk difference.
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the response probabilities for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlRiskDiff(riskDiffH0 = 0.1, n1 = 10, y1 = 4, n2 = 20, y2 = 0)
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector remlRiskDiff(const double riskDiffH0 = 0.0,
                           const double n1 = NA_REAL,
                           const double y1 = NA_REAL,
                           const double n2 = NA_REAL,
                           const double y2 = NA_REAL) {
  double n = n1 + n2;
  double y = y1 + y2;

  double p1, p2;
  if (fabs(riskDiffH0) < 1.0e-8) {
    p1 = y/n;
    p2 = p1;
  } else {
    double L3 = n;
    double L2 = (n1 + 2*n2)*riskDiffH0 - n - y;
    double L1 = (n2*riskDiffH0 - n - 2*y2)*riskDiffH0 + y;
    double L0 = y2*riskDiffH0*(1 - riskDiffH0);

    double q = pow(L2,3)/pow(3*L3,3) - L1*L2/(6*pow(L3,2)) + L0/(2*L3);
    double p = (2*(q>0) - 1)*sqrt(pow(L2,2)/pow(3*L3,2) - L1/(3*L3));

    double pi = 2*acos(0.0);
    double a = 1/3.0*(pi + acos(std::min(std::max(q/pow(p,3), -1.0), 1.0)));
    p2 = 2*p*cos(a) - L2/(3*L3);
    p1 = p2 + riskDiffH0;
  }

  return NumericVector::create(p1, p2);
}


// [[Rcpp::export]]
DataFrame remlRiskDiff2(const double riskDiffH0 = 0.0,
                        const NumericVector& n1 = NA_REAL,
                        const NumericVector& y1 = NA_REAL,
                        const NumericVector& n2 = NA_REAL,
                        const NumericVector& y2 = NA_REAL) {
  int k = n1.size();
  NumericVector p1(k), p2(k);
  for (int i=0; i<k; i++) {
    NumericVector a = remlRiskDiff(riskDiffH0, n1[i], y1[i], n2[i], y2[i]);
    p1[i] = a[0];
    p2[i] = a[1];
  }

  return DataFrame::create(_["p1"] = p1, _["p2"] = p2);
}


//' @title Miettinen-Nurminen score test statistic for two-sample risk
//' difference
//' @description Obtains the Miettinen-Nurminen score test statistic
//' for two-sample risk difference possibly with stratification.
//'
//' @param riskDiffH0 The risk difference under the null hypothesis.
//'   Defaults to 0.
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatRiskDiff(riskDiffH0 = 0, n1 = c(10,10), y1 = c(4,3),
//'               n2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
 // [[Rcpp::export]]
double zstatRiskDiff(const double riskDiffH0 = 0.0,
                     const NumericVector& n1 = NA_REAL,
                     const NumericVector& y1 = NA_REAL,
                     const NumericVector& n2 = NA_REAL,
                     const NumericVector& y2 = NA_REAL) {
  List mr = remlRiskDiff2(riskDiffH0, n1, y1, n2, y2);
  NumericVector p1 = as<NumericVector>(mr["p1"]);
  NumericVector p2 = as<NumericVector>(mr["p2"]);

  int i, k = n1.size();
  NumericVector n = n1 + n2;
  NumericVector w(k), md(k), mv(k);
  for (i=0; i<k; i++) {
    w[i] = n1[i]*n2[i]/n[i];
    md[i] = y1[i]/n1[i] - y2[i]/n2[i] - riskDiffH0;
    mv[i] = p1[i]*(1-p1[i])/n1[i] + p2[i]*(1-p2[i])/n2[i];
    mv[i] = std::max(mv[i]*n[i]/(n[i]-1.0), 1.0e-8);
  }

  w = w/sum(w);
  return sum(w*md)/sqrt(sum(pow(w,2)*mv));
}


//' @title Miettinen-Nurminen score confidence interval for
//' two-sample risk difference
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample risk difference possibly with
//' stratification.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input sample size
//'   and number of responses for each treatment group.
//'   It has the following variables:
//'
//'     - \code{n1}: The sample size for the active treatment group.
//'
//'     - \code{y1}: The number of responses for the active treatment group.
//'
//'     - \code{n2}: The sample size for the control group.
//'
//'     - \code{y2}: The number of responses for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for risk difference. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnRiskDiffCI(n1 = c(10,10), y1 = c(4,3), n2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
List mnRiskDiffCI(const NumericVector& n1 = NA_REAL,
                  const NumericVector& y1 = NA_REAL,
                  const NumericVector& n2 = NA_REAL,
                  const NumericVector& y2 = NA_REAL,
                  const double cilevel = 0.95) {
  if (is_true(any(n1 <= 0))) {
    stop("Each element of n1 must be a positive integer");
  }

  if (is_true(any((y1 < 0) | (y1 > n1)))) {
    stop("Each element of y1 must be an integer between 0 and n1");
  }

  if (is_true(any(n2 <= 0))) {
    stop("Each element of n2 must be a positive integer");
  }

  if (is_true(any((y2 < 0) | (y2 > n2)))) {
    stop("Each element of y2 must be an integer between 0 and n2");
  }

  if (cilevel <= 0 || cilevel >= 1) {
    stop("cilevel must lie between 0 and 1");
  }

  int i, k = n1.size();
  NumericVector n = n1 + n2;
  NumericVector w(k);
  for (i=0; i<k; i++) {
    w[i] = n1[i]*n2[i]/n[i];
  }
  w = w/sum(w);

  double estimate = 0;
  for (i=0; i<k; i++) {
    estimate += w[i]*(y1[i]/n1[i] - y2[i]/n2[i]);
  }

  double b = R::qnorm((1 + cilevel)/2, 0, 1, 1, 0);

  auto f1 = [n1, y1, n2, y2, b](double riskDiffH0) {
    return zstatRiskDiff(riskDiffH0, n1, y1, n2, y2) - b;
  };

  double lower = brent(f1, -1.0, estimate, 1.0e-6);

  auto f2 = [n1, y1, n2, y2, b](double riskDiffH0) {
    return zstatRiskDiff(riskDiffH0, n1, y1, n2, y2) + b;
  };

  double upper = brent(f2, estimate, 1.0, 1.0e-6);

  DataFrame data = DataFrame::create(
    _["n1"] = n1,
    _["y1"] = y1,
    _["n2"] = n2,
    _["y2"] = y2);

  DataFrame estimates = DataFrame::create(
    _["scale"] = "risk difference",
    _["estimate"] = estimate,
    _["lower"] = lower,
    _["upper"] = upper,
    _["cilevel"] = cilevel);

  List des = List::create(
    _["data"] = data,
    _["estimates"] = estimates);

  des.attr("class") = "estimateCI";

  return des;
}


//' @title REML estimates of individual proportions with specified risk
//' ratio
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified risk ratio.
//'
//' @param riskRatioH0 The specified risk ratio.
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the response probabilities for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlRiskRatio(riskRatioH0 = 1.2, n1 = 10, y1 = 4, n2 = 20, y2 = 2)
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector remlRiskRatio(const double riskRatioH0 = 1.0,
                            const double n1 = NA_REAL,
                            const double y1 = NA_REAL,
                            const double n2 = NA_REAL,
                            const double y2 = NA_REAL) {
  double n = n1 + n2;
  double y = y1 + y2;

  double p1, p2;
  if (fabs(riskRatioH0 - 1) < 1.0e-8) {
    p1 = y/n;
    p2 = p1;
  } else {
    double a = n*riskRatioH0;
    double b = -(n1*riskRatioH0 + y1 + n2 + y2*riskRatioH0);
    double c = y;

    p2 = (-b - sqrt(pow(b,2) - 4*a*c))/(2*a);
    p1 = p2*riskRatioH0;
  }

  return NumericVector::create(p1, p2);
}


// [[Rcpp::export]]
DataFrame remlRiskRatio2(const double riskRatioH0 = 1.0,
                         const NumericVector& n1 = NA_REAL,
                         const NumericVector& y1 = NA_REAL,
                         const NumericVector& n2 = NA_REAL,
                         const NumericVector& y2 = NA_REAL) {
  int k = n1.size();
  NumericVector p1(k), p2(k);
  for (int i=0; i<k; i++) {
    NumericVector a = remlRiskRatio(riskRatioH0, n1[i], y1[i], n2[i], y2[i]);
    p1[i] = a[0];
    p2[i] = a[1];
  }

  return DataFrame::create(_["p1"] = p1, _["p2"] = p2);
}


//' @title Miettinen-Nurminen score test statistic for two-sample risk ratio
//' @description Obtains the Miettinen-Nurminen score test statistic for
//' two-sample risk ratio possibly with stratification.
//'
//' @param riskRatioH0 The risk ratio under the null hypothesis.
//'   Defaults to 1.
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatRiskRatio(riskRatioH0 = 1, n1 = c(10,10), y1 = c(4,3),
//'                n2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
double zstatRiskRatio(const double riskRatioH0 = 1.0,
                      const NumericVector& n1 = NA_REAL,
                      const NumericVector& y1 = NA_REAL,
                      const NumericVector& n2 = NA_REAL,
                      const NumericVector& y2 = NA_REAL) {
  List mr = remlRiskRatio2(riskRatioH0, n1, y1, n2, y2);
  NumericVector p1 = as<NumericVector>(mr["p1"]);
  NumericVector p2 = as<NumericVector>(mr["p2"]);

  int i, k = n1.size();
  NumericVector n = n1 + n2;
  NumericVector w(k), md(k), mv(k);
  for (i=0; i<k; i++) {
    w[i] = n1[i]*n2[i]/n[i];
    md[i] = y1[i]/n1[i] - y2[i]/n2[i]*riskRatioH0;
    mv[i] = p1[i]*(1-p1[i])/n1[i] + pow(riskRatioH0,2)*p2[i]*(1-p2[i])/n2[i];
    mv[i] = std::max(mv[i]*n[i]/(n[i]-1.0), 1.0e-8);
  }

  w = w/sum(w);
  return sum(w*md)/sqrt(sum(pow(w,2)*mv));
}


//' @title Miettinen-Nurminen score confidence interval for
//' two-sample risk ratio
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample risk ratio possibly with
//' stratification.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input sample size
//'   and number of responses for each treatment group.
//'   It has the following variables:
//'
//'     - \code{n1}: The sample size for the active treatment group.
//'
//'     - \code{y1}: The number of responses for the active treatment group.
//'
//'     - \code{n2}: The sample size for the control group.
//'
//'     - \code{y2}: The number of responses for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for risk ratio. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnRiskRatioCI(n1 = c(10,10), y1 = c(4,3), n2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
List mnRiskRatioCI(const NumericVector& n1 = NA_REAL,
                   const NumericVector& y1 = NA_REAL,
                   const NumericVector& n2 = NA_REAL,
                   const NumericVector& y2 = NA_REAL,
                   const double cilevel = 0.95) {
  if (is_true(any(n1 <= 0))) {
    stop("Each element of n1 must be a positive integer");
  }

  if (is_true(any((y1 < 0) | (y1 > n1)))) {
    stop("Each element of y1 must be an integer between 0 and n1");
  }

  if (is_true(any(n2 <= 0))) {
    stop("Each element of n2 must be a positive integer");
  }

  if (is_true(any((y2 < 0) | (y2 > n2)))) {
    stop("Each element of y2 must be an integer between 0 and n2");
  }

  if (cilevel <= 0 || cilevel >= 1) {
    stop("cilevel must lie between 0 and 1");
  }

  int i, k = n1.size();
  NumericVector n = n1 + n2;
  NumericVector w(k);
  for (i=0; i<k; i++) {
    w[i] = n1[i]*n2[i]/n[i];
  }
  w = w/sum(w);

  double b = R::qnorm((1 + cilevel)/2, 0, 1, 1, 0);

  double estimate, lower, upper;

  if (is_true(any((y1 == 0) & (y2 == 0)))) {
    estimate = NA_REAL;
    lower = NA_REAL;
    upper = NA_REAL;
  } else if (is_true(all(y1 == 0))) {
    estimate = 0;
    lower = 0;

    auto f2 = [n1, y1, n2, y2, b](double riskRatioH0) {
      return zstatRiskRatio(riskRatioH0, n1, y1, n2, y2) + b;
    };

    upper = brent(f2, 0.001, 1000.0, 1.0e-6);
  } else if (is_true(all(y2 == 0))) {
    estimate = R_PosInf;
    upper = R_PosInf;

    auto f1 = [n1, y1, n2, y2, b](double riskRatioH0) {
      return zstatRiskRatio(riskRatioH0, n1, y1, n2, y2) - b;
    };

    lower = brent(f1, 0.001, 1000.0, 1.0e-6);
  } else {
    double p1 = 0, p2 = 0;
    for (i=0; i<k; i++) {
      p1 += w[i]*y1[i]/n1[i];
      p2 += w[i]*y2[i]/n2[i];
    }

    estimate = p1/p2;

    auto f1 = [n1, y1, n2, y2, b](double riskRatioH0) {
      return zstatRiskRatio(riskRatioH0, n1, y1, n2, y2) - b;
    };

    lower = brent(f1, 0.001, estimate, 1.0e-6);

    auto f2 = [n1, y1, n2, y2, b](double riskRatioH0) {
      return zstatRiskRatio(riskRatioH0, n1, y1, n2, y2) + b;
    };

    upper = brent(f2, estimate, 1000.0, 1.0e-6);
  }

  DataFrame data = DataFrame::create(
    _["n1"] = n1,
    _["y1"] = y1,
    _["n2"] = n2,
    _["y2"] = y2);

  DataFrame estimates = DataFrame::create(
    _["scale"] = "risk ratio",
    _["estimate"] = estimate,
    _["lower"] = lower,
    _["upper"] = upper,
    _["cilevel"] = cilevel);

  List des = List::create(
    _["data"] = data,
    _["estimates"] = estimates);

  des.attr("class") = "estimateCI";

  return des;
}


//' @title REML estimates of individual proportions with specified odds
//' ratio
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified odds ratio.
//'
//' @param oddsRatioH0 The specified odds ratio.
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the response probabilities for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlOddsRatio(oddsRatioH0 = 1.25, n1 = 10, y1 = 4, n2 = 20, y2 = 2)
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector remlOddsRatio(const double oddsRatioH0 = 1.0,
                            const double n1 = NA_REAL,
                            const double y1 = NA_REAL,
                            const double n2 = NA_REAL,
                            const double y2 = NA_REAL) {
  double n = n1 + n2;
  double y = y1 + y2;

  double p1, p2;
  if (fabs(oddsRatioH0 - 1) < 1.0e-8) {
    p1 = y/n;
    p2 = p1;
  } else {
    double a = n2*(oddsRatioH0 - 1);
    double b = n1*oddsRatioH0 + n2 - y*(oddsRatioH0 - 1);
    double c = -y;
    p2 = (-b + sqrt(pow(b,2) - 4*a*c))/(2*a);
    p1 = p2*oddsRatioH0/(1 + p2*(oddsRatioH0 - 1));
  }

  return NumericVector::create(p1, p2);
}


// [[Rcpp::export]]
DataFrame remlOddsRatio2(const double oddsRatioH0 = 1.0,
                         const NumericVector& n1 = NA_REAL,
                         const NumericVector& y1 = NA_REAL,
                         const NumericVector& n2 = NA_REAL,
                         const NumericVector& y2 = NA_REAL) {
  int k = n1.size();
  NumericVector p1(k), p2(k);
  for (int i=0; i<k; i++) {
    NumericVector a = remlOddsRatio(oddsRatioH0, n1[i], y1[i], n2[i], y2[i]);
    p1[i] = a[0];
    p2[i] = a[1];
  }

  return DataFrame::create(_["p1"] = p1, _["p2"] = p2);
}


//' @title Miettinen-Nurminen score test statistic for two-sample odds ratio
//' @description Obtains the Miettinen-Nurminen score test statistic for
//' two-sample odds ratio possibly with stratification.
//'
//' @param oddsRatioH0 The odds ratio under the null hypothesis.
//'   Defaults to 1.
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatOddsRatio(oddsRatioH0 = 1, n1 = c(10,10), y1 = c(4,3),
//'                n2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
double zstatOddsRatio(const double oddsRatioH0 = 1.0,
                      const NumericVector& n1 = NA_REAL,
                      const NumericVector& y1 = NA_REAL,
                      const NumericVector& n2 = NA_REAL,
                      const NumericVector& y2 = NA_REAL) {
  List mr = remlOddsRatio2(oddsRatioH0, n1, y1, n2, y2);
  NumericVector p1 = as<NumericVector>(mr["p1"]);
  NumericVector p2 = as<NumericVector>(mr["p2"]);

  int i, k = n1.size();
  NumericVector n = n1 + n2;
  NumericVector w(k), md(k), mv(k);
  for (i=0; i<k; i++) {
    w[i] = n1[i]*n2[i]/n[i];
    md[i] = (y1[i]/n1[i] - p1[i])/(p1[i]*(1-p1[i])) -
      (y2[i]/n2[i] - p2[i])/(p2[i]*(1-p2[i]));
    mv[i] = 1/(n1[i]*p1[i]*(1-p1[i])) + 1/(n2[i]*p2[i]*(1-p2[i]));
    mv[i] = std::max(mv[i]*n[i]/(n[i]-1.0), 1.0e-8);
  }

  w = w/sum(w);
  return sum(w*md)/sqrt(sum(pow(w,2)*mv));
}


//' @title Miettinen-Nurminen score confidence interval for
//' two-sample odds ratio
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample odds ratio possibly with
//' stratification.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel sample size weights are used for stratified
//' samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input sample size
//'   and number of responses for each treatment group.
//'   It has the following variables:
//'
//'     - \code{n1}: The sample size for the active treatment group.
//'
//'     - \code{y1}: The number of responses for the active treatment group.
//'
//'     - \code{n2}: The sample size for the control group.
//'
//'     - \code{y2}: The number of responses for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for odds ratio. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnOddsRatioCI(n1 = c(10,10), y1 = c(4,3), n2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
List mnOddsRatioCI(const NumericVector& n1 = NA_REAL,
                   const NumericVector& y1 = NA_REAL,
                   const NumericVector& n2 = NA_REAL,
                   const NumericVector& y2 = NA_REAL,
                   const double cilevel = 0.95) {
  if (is_true(any(n1 <= 0))) {
    stop("Each element of n1 must be a positive integer");
  }

  if (is_true(any((y1 < 0) | (y1 > n1)))) {
    stop("Each element of y1 must be an integer between 0 and n1");
  }

  if (is_true(any(n2 <= 0))) {
    stop("Each element of n2 must be a positive integer");
  }

  if (is_true(any((y2 < 0) | (y2 > n2)))) {
    stop("Each element of y2 must be an integer between 0 and n2");
  }

  if (cilevel <= 0 || cilevel >= 1) {
    stop("cilevel must lie between 0 and 1");
  }

  int i, k = n1.size();
  NumericVector n = n1 + n2;
  NumericVector w(k);
  for (i=0; i<k; i++) {
    w[i] = n1[i]*n2[i]/n[i];
  }
  w = w/sum(w);

  double b = R::qnorm((1 + cilevel)/2, 0, 1, 1, 0);

  double estimate, lower, upper;

  if (is_true(any((y1 == 0) & (y2 == 0))) ||
      is_true(any((y1 == n1) & (y2 == n2)))) {
    estimate = NA_REAL;
    lower = NA_REAL;
    upper = NA_REAL;
  } else if (is_true(all(y1 == 0))) {
    estimate = 0;
    lower = 0;

    auto f2 = [n1, y1, n2, y2, b](double oddsRatioH0) {
      return zstatOddsRatio(oddsRatioH0, n1, y1, n2, y2) + b;
    };

    upper = brent(f2, 0.001, 1000.0, 1.0e-6);
  } else if (is_true(all(y2 == 0))) {
    estimate = R_PosInf;
    upper = R_PosInf;

    auto f1 = [n1, y1, n2, y2, b](double oddsRatioH0) {
      return zstatOddsRatio(oddsRatioH0, n1, y1, n2, y2) - b;
    };

    lower = brent(f1, 0.001, 1000.0, 1.0e-6);
  } else {
    auto f0 = [n1, y1, n2, y2](double oddsRatioH0) {
      return zstatOddsRatio(oddsRatioH0, n1, y1, n2, y2);
    };

    estimate = brent(f0, 0.001, 1000.0, 1.0e-6);

    auto f1 = [n1, y1, n2, y2, b](double oddsRatioH0) {
      return zstatOddsRatio(oddsRatioH0, n1, y1, n2, y2) - b;
    };

    lower = brent(f1, 0.001, estimate, 1e-6);

    auto f2 = [n1, y1, n2, y2, b](double oddsRatioH0) {
      return zstatOddsRatio(oddsRatioH0, n1, y1, n2, y2) + b;
    };

    upper = brent(f2, estimate, 1000.0, 1.0e-6);
  }

  DataFrame data = DataFrame::create(
    _["n1"] = n1,
    _["y1"] = y1,
    _["n2"] = n2,
    _["y2"] = y2);

  DataFrame estimates = DataFrame::create(
    _["scale"] = "odds ratio",
    _["estimate"] = estimate,
    _["lower"] = lower,
    _["upper"] = upper,
    _["cilevel"] = cilevel);

  List des = List::create(
    _["data"] = data,
    _["estimates"] = estimates);

  des.attr("class") = "estimateCI";

  return des;
}


//' @title REML estimates of individual rates with specified rate
//' difference
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified rate difference.
//'
//' @param rateDiffH0 The specified rate difference.
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the incidence rates for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlRateDiff(rateDiffH0 = 0.1, t1 = 10, y1 = 4, t2 = 20, y2 = 2)
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector remlRateDiff(const double rateDiffH0 = 0.0,
                           const double t1 = NA_REAL,
                           const double y1 = NA_REAL,
                           const double t2 = NA_REAL,
                           const double y2 = NA_REAL) {
  double t = t1 + t2;
  double y = y1 + y2;

  double r1, r2;
  if (fabs(rateDiffH0) < 1.0e-8) {
    r1 = y/t;
    r2 = r1;
  } else {
    double a = t;
    double b = t*rateDiffH0 - y;
    double c = -y2*rateDiffH0;

    r2 = (-b + sqrt(pow(b,2) - 4*a*c))/(2*a);
    r1 = r2 + rateDiffH0;
  }

  return NumericVector::create(r1, r2);
}


// [[Rcpp::export]]
DataFrame remlRateDiff2(const double rateDiffH0 = 0.0,
                        const NumericVector& t1 = NA_REAL,
                        const NumericVector& y1 = NA_REAL,
                        const NumericVector& t2 = NA_REAL,
                        const NumericVector& y2 = NA_REAL) {
  int k = t1.size();
  NumericVector r1(k), r2(k);
  for (int i=0; i<k; i++) {
    NumericVector a = remlRateDiff(rateDiffH0, t1[i], y1[i], t2[i], y2[i]);
    r1[i] = a[0];
    r2[i] = a[1];
  }

  return DataFrame::create(_["r1"] = r1, _["r2"] = r2);
}


//' @title Miettinen-Nurminen score test statistic for two-sample rate
//' difference
//' @description Obtains the Miettinen-Nurminen score test statistic for
//' two-sample rate difference possibly with stratification.
//'
//' @param rateDiffH0 The rate difference under the null hypothesis.
//'   Defaults to 0.
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//'
//' @details
//' The Mantel-Haenszel weights are used for stratified samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatRateDiff(rateDiffH0 = 0, t1 = c(10,10), y1 = c(4,3),
//'               t2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
double zstatRateDiff(const double rateDiffH0 = 0.0,
                     const NumericVector& t1 = NA_REAL,
                     const NumericVector& y1 = NA_REAL,
                     const NumericVector& t2 = NA_REAL,
                     const NumericVector& y2 = NA_REAL) {
  List mr = remlRateDiff2(rateDiffH0, t1, y1, t2, y2);
  NumericVector r1 = as<NumericVector>(mr["r1"]);
  NumericVector r2 = as<NumericVector>(mr["r2"]);

  int i, k = t1.size();
  NumericVector t = t1 + t2;
  NumericVector w(k), md(k), mv(k);
  for (i=0; i<k; i++) {
    w[i] = t1[i]*t2[i]/t[i];
    md[i] = y1[i]/t1[i] - y2[i]/t2[i] - rateDiffH0;
    mv[i] = r1[i]/t1[i] + r2[i]/t2[i];
    mv[i] = std::max(mv[i], 1.0e-8);
  }

  w = w/sum(w);
  return sum(w*md)/sqrt(sum(pow(w,2)*mv));
}


//' @title Miettinen-Nurminen score confidence interval for
//' two-sample rate difference
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample rate difference possibly with
//' stratification.
//'
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel weights are used for stratified samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input exposure
//'   and number of events for each treatment group.
//'   It has the following variables:
//'
//'     - \code{t1}: The exposure for the active treatment group.
//'
//'     - \code{y1}: The number of events for the active treatment group.
//'
//'     - \code{t2}: The exposure for the control group.
//'
//'     - \code{y2}: The number of events for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for rate difference. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnRateDiffCI(t1 = c(10,10), y1 = c(4,3), t2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
List mnRateDiffCI(const NumericVector& t1 = NA_REAL,
                  const NumericVector& y1 = NA_REAL,
                  const NumericVector& t2 = NA_REAL,
                  const NumericVector& y2 = NA_REAL,
                  const double cilevel = 0.95) {
  if (is_true(any(t1 <= 0))) {
    stop("Each element of t1 must be positive");
  }

  if (is_true(any(y1 < 0))) {
    stop("Each element of y1 must be a nonnegative integer");
  }

  if (is_true(any(t2 <= 0))) {
    stop("Each element of t2 must be positive");
  }

  if (is_true(any(y2 < 0))) {
    stop("Each element of y2 must be a nonnegative integer");
  }

  if (cilevel <= 0 || cilevel >= 1) {
    stop("cilevel must lie between 0 and 1");
  }


  int i, k = t1.size();
  NumericVector t = t1 + t2;
  NumericVector w = t1*t2/t;
  w = w/sum(w);

  double estimate = 0;
  for (i=0; i<k; i++) {
    estimate += w[i]*(y1[i]/t1[i] - y2[i]/t2[i]);
  }

  double b = R::qnorm((1 + cilevel)/2, 0, 1, 1, 0);

  auto f1 = [t1, y1, t2, y2, b](double rateDiffH0) {
    return zstatRateDiff(rateDiffH0, t1, y1, t2, y2) - b;
  };

  double lower = brent(f1, -1000.0, estimate, 1.0e-6);

  auto f2 = [t1, y1, t2, y2, b](double rateDiffH0) {
    return zstatRateDiff(rateDiffH0, t1, y1, t2, y2) + b;
  };

  double upper = brent(f2, estimate, 1000.0, 1.0e-6);

  DataFrame data = DataFrame::create(
    _["t1"] = t1,
    _["y1"] = y1,
    _["t2"] = t2,
    _["y2"] = y2);

  DataFrame estimates = DataFrame::create(
    _["scale"] = "rate difference",
    _["estimate"] = estimate,
    _["lower"] = lower,
    _["upper"] = upper,
    _["cilevel"] = cilevel);

  List des = List::create(
    _["data"] = data,
    _["estimates"] = estimates);

  des.attr("class") = "estimateCI";

  return des;
}


//' @title REML estimates of individual rates with specified rate
//' ratio
//' @description Obtains the restricted maximum likelihood estimates of
//' individual proportions with specified rate ratio.
//'
//' @param rateRatioH0 The specified rate ratio.
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//'
//' @return A vector of the restricted maximum likelihood estimates
//' of the incidence rates for the two treatment groups.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' remlRateRatio(rateRatioH0 = 1.1, t1 = 10, y1 = 4, t2 = 20, y2 = 2)
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector remlRateRatio(const double rateRatioH0 = 1.0,
                            const double t1 = NA_REAL,
                            const double y1 = NA_REAL,
                            const double t2 = NA_REAL,
                            const double y2 = NA_REAL) {
  double r2 = (y1 + y2)/(t1*rateRatioH0 + t2);
  double r1 = r2*rateRatioH0;

  return NumericVector::create(r1, r2);
}


// [[Rcpp::export]]
DataFrame remlRateRatio2(const double rateRatioH0 = 1.0,
                         const NumericVector& t1 = NA_REAL,
                         const NumericVector& y1 = NA_REAL,
                         const NumericVector& t2 = NA_REAL,
                         const NumericVector& y2 = NA_REAL) {
  int k = t1.size();
  NumericVector r1(k), r2(k);
  for (int i=0; i<k; i++) {
    NumericVector a = remlRateRatio(rateRatioH0, t1[i], y1[i], t2[i], y2[i]);
    r1[i] = a[0];
    r2[i] = a[1];
  }

  return DataFrame::create(_["r1"] = r1, _["r2"] = r2);
}


//' @title Miettinen-Nurminen score test statistic for two-sample rate ratio
//' @description Obtains the Miettinen-Nurminen score test statistic for
//' two-sample rate ratio possibly with stratification.
//'
//' @param rateRatioH0 The rate ratio under the null hypothesis.
//'   Defaults to 1.
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//'
//' @details
//' The Mantel-Haenszel weights are used for stratified samples.
//'
//' @return The value of the score test statistic.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' zstatRateRatio(rateRatioH0 = 1, t1 = c(10,10), y1 = c(4,3),
//'                t2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
double zstatRateRatio(const double rateRatioH0 = 1.0,
                      const NumericVector& t1 = NA_REAL,
                      const NumericVector& y1 = NA_REAL,
                      const NumericVector& t2 = NA_REAL,
                      const NumericVector& y2 = NA_REAL) {
  List mr = remlRateRatio2(rateRatioH0, t1, y1, t2, y2);
  NumericVector r1 = as<NumericVector>(mr["r1"]);
  NumericVector r2 = as<NumericVector>(mr["r2"]);

  NumericVector t = t1 + t2;
  NumericVector w = t1*t2/t;
  w = w/sum(w);

  int i, k = t1.size();
  NumericVector md(k), mv(k);
  for (i=0; i<k; i++) {
    md[i] = y1[i]/t1[i] - (y2[i]/t2[i])*rateRatioH0;
    mv[i] = r1[i]/t1[i] + pow(rateRatioH0,2)*r2[i]/t2[i];
    mv[i] = std::max(mv[i], 1.0e-8);
  }

  return sum(w*md)/sqrt(sum(pow(w,2)*mv));
}


//' @title Miettinen-Nurminen score confidence interval for
//' two-sample rate ratio
//' @description Obtains the Miettinen-Nurminen score confidence
//' interval for two-sample rate ratio possibly with
//' stratification.
//'
//' @param t1 The exposure for the active treatment group.
//' @param y1 The number of events for the active treatment group.
//' @param t2 The exposure for the control group.
//' @param y2 The number of events for the control group.
//' @param cilevel The confidence interval level.
//'
//' @details
//' The Mantel-Haenszel weights are used for stratified samples.
//'
//' @return A list with two components:
//'
//' * \code{data} A data frame containing the input exposure
//'   and number of events for each treatment group.
//'   It has the following variables:
//'
//'     - \code{t1}: The exposure for the active treatment group.
//'
//'     - \code{y1}: The number of events for the active treatment group.
//'
//'     - \code{t2}: The exposure for the control group.
//'
//'     - \code{y2}: The number of events for the control group.
//'
//' * \code{estimates}: A data frame containing the point estimate
//'   and confidence interval for rate ratio. It has the following
//'   variables:
//'
//'     - \code{scale}: The scale of treatment effect.
//'
//'     - \code{estimate}: The point estimate.
//'
//'     - \code{lower}: The lower limit of the confidence interval.
//'
//'     - \code{upper}: The upper limit of the confidence interval.
//'
//'     - \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' mnRateRatioCI(t1 = c(10,10), y1 = c(4,3), t2 = c(20,10), y2 = c(2,0))
//'
//' @export
//'
// [[Rcpp::export]]
List mnRateRatioCI(const NumericVector& t1,
                   const NumericVector& y1,
                   const NumericVector& t2,
                   const NumericVector& y2,
                   const double cilevel = 0.95) {
  if (is_true(any(t1 <= 0))) {
    stop("Each element of t1 must be positive");
  }

  if (is_true(any(y1 < 0))) {
    stop("Each element of y1 must be a nonnegative integer");
  }

  if (is_true(any(t2 <= 0))) {
    stop("Each element of t2 must be positive");
  }

  if (is_true(any(y2 < 0))) {
    stop("Each element of y2 must be a nonnegative integer");
  }

  if (cilevel <= 0 || cilevel >= 1) {
    stop("cilevel must lie between 0 and 1");
  }


  int i, k = t1.size();
  NumericVector t = t1 + t2;
  NumericVector w = t1*t2/t;
  w = w/sum(w);

  double b = R::qnorm((1 + cilevel)/2, 0, 1, 1, 0);

  double estimate, lower, upper;
  if (is_true(any((y1 == 0) & (y2 == 0)))) {
    estimate = NA_REAL;
    lower = NA_REAL;
    upper = NA_REAL;
  } else if (is_true(all(y1 == 0))) {
    estimate = 0;
    lower = 0;

    auto f2 = [t1, y1, t2, y2, b](double rateRatioH0) {
      return zstatRateRatio(rateRatioH0, t1, y1, t2, y2) + b;
    };

    upper = brent(f2, 0.001, 1000.0, 1.0e-6);
  } else if (is_true(all(y2 == 0))) {
    estimate = NA_REAL;
    upper = NA_REAL;

    auto f1 = [t1, y1, t2, y2, b](double rateRatioH0) {
      return zstatRateRatio(rateRatioH0, t1, y1, t2, y2) - b;
    };

    lower = brent(f1, 0.001, 1000.0, 1.0e-6);
  } else {
    double r1=0, r2=0;
    for (i=0; i<k; i++) {
      r1 += w[i]*y1[i]/t1[i];
      r2 += w[i]*y2[i]/t2[i];
    }
    estimate = r1/r2;

    auto f1 = [t1, y1, t2, y2, b](double rateRatioH0) {
      return zstatRateRatio(rateRatioH0, t1, y1, t2, y2) - b;
    };

    lower = brent(f1, 0.001, estimate, 1.0e-6);

    auto f2 = [t1, y1, t2, y2, b](double rateRatioH0) {
      return zstatRateRatio(rateRatioH0, t1, y1, t2, y2) + b;
    };

    upper = brent(f2, estimate, 1000.0, 1.0e-6);
  }

  DataFrame data = DataFrame::create(
    _["t1"] = t1,
    _["y1"] = y1,
    _["t2"] = t2,
    _["y2"] = y2);

  DataFrame estimates = DataFrame::create(
    _["scale"] = "rate ratio",
    _["estimate"] = estimate,
    _["lower"] = lower,
    _["upper"] = upper,
    _["cilevel"] = cilevel);

  List des = List::create(
    _["data"] = data,
    _["estimates"] = estimates);

  des.attr("class") = "estimateCI";

  return des;
}


//' @title Power for exact unconditional test of risk difference
//' @description Obtains the power given sample size for exact unconditional
//' test of risk difference.
//'
//' @param n The total sample size.
//' @param riskDiffH0 The risk difference under the null hypothesis.
//'   Defaults to 0.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified one-sided significance level.
//'
//' * \code{attainedAlpha}: The attained one-sided significance level if
//'   requested.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskDiffH0}: The risk difference under the null hypothesis.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskDiffBound}: The critical value on the scale of
//'   score test statistic for risk difference.
//'
//' * \code{pi2star}: The response probability in the control group
//'   at which the critical value of the test statistic is attained.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' powerRiskDiffExact(n = 68, pi1 = 0.6, pi2 = 0.25, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
DataFrame powerRiskDiffExact(
    const int n = NA_INTEGER,
    const double riskDiffH0 = 0,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.025,
    const bool calculateAttainedAlpha = 1) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  int n1 = std::round(n*r);
  int n2 = n - n1;
  IntegerVector x1 = seq(0, n1);
  IntegerVector x2 = seq(0, n2);
  int i, k = (n1+1)*(n2+1);

  // The set of values for the test statistic
  NumericVector T(k);
  i = 0;
  for (int y1=0; y1<=n1; y1++) {
    for (int y2=0; y2<=n2; y2++) {
      NumericVector a = remlRiskDiff(riskDiffH0, n1, y1, n2, y2);
      double p1 = a[0], p2 = a[1];
      double p1hat = y1/(n1+0.0), p2hat = y2/(n2+0.0);
      double md = p1hat - p2hat - riskDiffH0;
      double mv = std::max(p1*(1-p1)/n1 + p2*(1-p2)/n2, 1.0e-8);
      T[i++] = md/sqrt(mv);
    }
  }

  // sort T in ascending order
  IntegerVector order = seq(0, k-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return T[i] < T[j];
  });

  NumericVector Tsorted = T[order];

  // identify the locations of the unique values of T
  IntegerVector idx(1,0);
  for (i=1; i<k; i++) {
    if (Tsorted[i] != Tsorted[i-1]) {
      idx.push_back(i);
    }
  }

  NumericVector Tunique = Tsorted[idx];

  int k1 = idx.size();  // how many unique values of T
  idx.push_back(k);     // add the upper bound to the index

  // whether higher values represent better responses
  bool directionUpper = pi1 - pi2 > riskDiffH0;

  // obtain the critical value of the test statistic
  auto f = [n1, n2, x1, x2, riskDiffH0, directionUpper, Tunique,
            k, k1, idx, order, alpha](double p2)->double {
              // table probabilities under H0
              NumericVector q1 = dbinom(x1, n1, p2 + riskDiffH0, 0);
              NumericVector q2 = dbinom(x2, n2, p2, 0);
              NumericVector q(k);

              int i, j, l;
              i = 0;
              for (int y1=0; y1<=n1; y1++) {
                for (int y2=0; y2<=n2; y2++) {
                  q[i++] = q1[y1]*q2[y2];
                }
              }

              NumericVector qsorted = q[order];
              double aval;

              if (directionUpper) { // reject H0 if T >= t
                double s = 0;
                l = k-1;  // flattened index
                for (i=k1-1; i>=0; i--) {  // move down from the largest T
                  for (j=idx[i+1]-1; j>=idx[i]; j--) {
                    s += qsorted[l--];
                  }
                  if (s > alpha) break;
                }

                // backtrack to obtain the critical value of T
                if (i == k1-1) {
                  aval = Tunique[k1-1] + 1; // impossible to reject H0
                } else {
                  aval = Tunique[i+1];
                }
              } else { // reject H0 if T <= t
                double s = 0;
                l = 0;
                for (i=0; i<k1; i++) {  // move up from the smallest T
                  for (j=idx[i]; j<idx[i+1]; j++) {
                    s += qsorted[l++];
                  }
                  if (s > alpha) break;
                }

                // backtrack to obtain the critical value of T
                if (i == 0) {
                  aval = Tunique[0] - 1; // impossible to reject H0
                } else {
                  aval = Tunique[i-1];
                }
              }

              return directionUpper ? -aval : aval;
            };

  // find the critical value independent of pi2
  double pi2lower = std::max(0.0, -riskDiffH0);
  double pi2upper = std::min(1.0, 1-riskDiffH0);

  int K = 100;
  NumericVector a(K), b(K);
  double delta = (pi2upper - pi2lower)/K;
  double tol = 1.0e-6;
  for (i=0; i<K; i++) {
    double lower = pi2lower + i*delta;
    double upper = lower + delta;
    NumericVector out = mini(f, lower, upper, tol);
    a[i] = out[0];
    b[i] = out[1];
  }

  i = which_min(b);
  double pi2star = a[i];
  double t = directionUpper ? -b[i] : b[i];

  // calculate attained alpha and power
  auto g = [n1, n2, x1, x2, T, t,
            directionUpper](double p1, double p2)->double{
              NumericVector q1 = dbinom(x1, n1, p1, 0);
              NumericVector q2 = dbinom(x2, n2, p2, 0);

              double preject = 0;
              int i = 0;
              for (int y1=0; y1<=n1; y1++) {
                for (int y2=0; y2<=n2; y2++) {
                  if ((directionUpper && (T[i] >= t)) ||
                      ((!directionUpper) && (T[i] <= t))) {
                    preject += q1[y1]*q2[y2];
                  }
                  i++;
                }
              }

              return preject;
            };

  double power = g(pi1, pi2);

  DataFrame result;

  if (calculateAttainedAlpha) {
    auto h = [n1, n2, x1, x2, T, t, directionUpper,
              riskDiffH0](double p2)->double {
                double p1 = p2 + riskDiffH0;
                NumericVector q1 = dbinom(x1, n1, p1, 0);
                NumericVector q2 = dbinom(x2, n2, p2, 0);

                double preject = 0;
                int i = 0;
                for (int y1=0; y1<=n1; y1++) {
                  for (int y2=0; y2<=n2; y2++) {
                    if ((directionUpper && (T[i] >= t)) ||
                        ((!directionUpper) && (T[i] <= t))) {
                      preject += q1[y1]*q2[y2];
                    }
                    i++;
                  }
                }

                return -preject;
    };

    for (i=0; i<K; i++) {
      double lower = pi2lower + i*delta;
      double upper = lower + delta;
      NumericVector out = mini(h, lower, upper, tol);
      a[i] = out[0];
      b[i] = out[1];
    }

    i = which_min(b);
    double attainedAlpha = -b[i];

    result = DataFrame::create(
      _["alpha"] = alpha,
      _["attainedAlpha"] = attainedAlpha,
      _["power"] = power,
      _["n"] = n,
      _["riskDiffH0"] = riskDiffH0,
      _["pi1"] = pi1,
      _["pi2"] = pi2,
      _["allocationRatioPlanned"] = allocationRatioPlanned,
      _["zstatRiskDiffBound"] = t,
      _["pi2star"] = pi2star);
  } else {
    result = DataFrame::create(
      _["alpha"] = alpha,
      _["power"] = power,
      _["n"] = n,
      _["riskDiffH0"] = riskDiffH0,
      _["pi1"] = pi1,
      _["pi2"] = pi2,
      _["allocationRatioPlanned"] = allocationRatioPlanned,
      _["zstatRiskDiffBound"] = t,
      _["pi2star"] = pi2star);
  }

  return result;
}


//' @title Sample size for exact unconditional test of risk difference
//' @description Obtains the sample size given power for exact unconditional
//' test of risk difference.
//'
//' @param beta The type II error.
//' @param riskDiffH0 The risk difference under the null hypothesis.
//'   Defaults to 0.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The one-sided significance level.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified one-sided significance level.
//'
//' * \code{attainedAlpha}: The attained one-sided significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskDiffH0}: The risk difference under the null hypothesis.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskDiffBound}: The critical value on the scale of
//'   score test statistic for risk difference.
//'
//' * \code{pi2star}: The response probability in the control group
//'   at which the critical value of the test statistic is attained.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' samplesizeRiskDiffExact(beta = 0.2, riskDiffH0 = -0.2,
//'                         pi1 = 0.8, pi2 = 0.8, alpha = 0.025)
//'
//' @export
// [[Rcpp::export]]
DataFrame samplesizeRiskDiffExact(
    const double beta = NA_REAL,
    const double riskDiffH0 = 0,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.025) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  NumericVector mr = remlRiskDiff(riskDiffH0, r, r*pi1, 1-r, (1-r)*pi2);
  double p1 = mr[0], p2 = mr[1];

  double v0 = p1*(1-p1)/r + p2*(1-p2)/(1-r);
  double v1 = pi1*(1-pi1)/r + pi2*(1-pi2)/(1-r);
  double z0 = R::qnorm(1-alpha, 0, 1, 1, 0);
  double z1 = R::qnorm(1-beta, 0, 1, 1, 0);
  double theta = fabs(pi1 - pi2 - riskDiffH0);
  double n0 = pow(z0*sqrt(v0) + z1*sqrt(v1), 2)/pow(theta, 2);

  int n_lower = std::floor(n0), n_upper = std::ceil(10*n0), n;

  DataFrame a, b;
  a = powerRiskDiffExact(n0, riskDiffH0, pi1, pi2,
                         allocationRatioPlanned, alpha, 0);

  while (as<double>(a["power"]) >= 1-beta) {
    a = powerRiskDiffExact(--n_lower, riskDiffH0, pi1, pi2,
                           allocationRatioPlanned, alpha, 0);
  }

  for (n=n_lower+1; n<=n_upper; n++) {
    a = powerRiskDiffExact(n, riskDiffH0, pi1, pi2,
                           allocationRatioPlanned, alpha, 0);

    // ensure that the power exceeds 1-beta for the next 10 sample sizes
    if (as<double>(a["power"]) >= 1-beta) {
      int i;
      bool okay = 1;
      for (i=1; i<=10; i++) {
        b = powerRiskDiffExact(n+i, riskDiffH0, pi1, pi2,
                               allocationRatioPlanned, alpha, 0);
        if (as<double>(b["power"]) < 1-beta) {
          okay = 0;
          break;
        }
      }

      if (okay) {
        break;
      } else {
        n += i;
      }
    }
  }

  // calculate attained alpha;
  a = powerRiskDiffExact(n, riskDiffH0, pi1, pi2,
                         allocationRatioPlanned, alpha, 1);

  return a;
}


//' @title Power for exact unconditional test of risk ratio
//' @description Obtains the power given sample size for exact unconditional
//' test of risk ratio.
//'
//' @param n The total sample size.
//' @param riskRatioH0 The risk ratio under the null hypothesis.
//'   Defaults to 1.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The one-sided significance level. Defaults to 0.025.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified one-sided significance level.
//'
//' * \code{attainedAlpha}: The attained one-sided significance level if
//'   requested.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskRatioH0}: The risk ratio under the null hypothesis.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskRatioBound}: The critical value on the scale of
//'   score test statistic for risk ratio.
//'
//' * \code{pi2star}: The response probability in the control group
//'   at which the critical value of the test statistic is attained.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' powerRiskRatioExact(n = 68, pi1 = 0.6, pi2 = 0.25, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
DataFrame powerRiskRatioExact(
    const int n = NA_INTEGER,
    const double riskRatioH0 = 1,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.025,
    const bool calculateAttainedAlpha = 1) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  int n1 = std::round(n*r);
  int n2 = n - n1;
  IntegerVector x1 = seq(0, n1);
  IntegerVector x2 = seq(0, n2);
  int i, k = (n1+1)*(n2+1);

  // The set of values for the test statistic
  NumericVector T(k);
  i = 0;
  for (int y1=0; y1<=n1; y1++) {
    for (int y2=0; y2<=n2; y2++) {
      NumericVector a = remlRiskRatio(riskRatioH0, n1, y1, n2, y2);
      double p1 = a[0], p2 = a[1];
      double p1hat = y1/(n1+0.0), p2hat = y2/(n2+0.0);
      double md = p1hat - p2hat*riskRatioH0;
      double mv = std::max(p1*(1-p1)/n1 +
                           pow(riskRatioH0,2)*p2*(1-p2)/n2, 1.0e-8);
      T[i++] = md/sqrt(mv);
    }
  }

  // sort T in ascending order
  IntegerVector order = seq(0, k-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return T[i] < T[j];
  });

  NumericVector Tsorted = T[order];

  // identify the locations of the unique values of T
  IntegerVector idx(1,0);
  for (i=1; i<k; i++) {
    if (Tsorted[i] != Tsorted[i-1]) {
      idx.push_back(i);
    }
  }

  NumericVector Tunique = Tsorted[idx];

  int k1 = idx.size();  // how many unique values of T
  idx.push_back(k);     // add the upper bound to the index

  // whether higher values represent better responses
  bool directionUpper = pi1/pi2 > riskRatioH0;

  // obtain the critical value of the test statistic
  auto f = [n1, n2, x1, x2, riskRatioH0, directionUpper, Tunique,
            k, k1, idx, order, alpha](double p2)->double {
              // table probabilities under H0
              NumericVector q1 = dbinom(x1, n1, p2*riskRatioH0, 0);
              NumericVector q2 = dbinom(x2, n2, p2, 0);
              NumericVector q(k);

              int i, j, l;
              i = 0;
              for (int y1=0; y1<=n1; y1++) {
                for (int y2=0; y2<=n2; y2++) {
                  q[i++] = q1[y1]*q2[y2];
                }
              }

              NumericVector qsorted = q[order];
              double aval;

              if (directionUpper) { // reject H0 if T >= t
                double s = 0;
                l = k-1;  // flattened index
                for (i=k1-1; i>=0; i--) {  // move down from the largest T
                  for (j=idx[i+1]-1; j>=idx[i]; j--) {
                    s += qsorted[l--];
                  }
                  if (s > alpha) break;
                }

                // backtrack to obtain the critical value of T
                if (i == k1-1) {
                  aval = Tunique[k1-1] + 1; // impossible to reject H0
                } else {
                  aval = Tunique[i+1];
                }
              } else { // reject H0 if T <= t
                double s = 0;
                l = 0;
                for (i=0; i<k1; i++) {  // move up from the smallest T
                  for (j=idx[i]; j<idx[i+1]; j++) {
                    s += qsorted[l++];
                  }
                  if (s > alpha) break;
                }

                // backtrack to obtain the critical value of T
                if (i == 0) {
                  aval = Tunique[0] - 1; // impossible to reject H0
                } else {
                  aval = Tunique[i-1];
                }
              }

              return directionUpper ? -aval : aval;
            };

  // find the critical value independent of pi2
  double pi2lower = 0;
  double pi2upper = std::min(1.0, 1/riskRatioH0);

  int K = 100;
  NumericVector a(K), b(K);
  double delta = (pi2upper - pi2lower)/K;
  double tol = 1.0e-6;
  for (i=0; i<K; i++) {
    double lower = pi2lower + i*delta;
    double upper = lower + delta;
    NumericVector out = mini(f, lower, upper, tol);
    a[i] = out[0];
    b[i] = out[1];
  }

  i = which_min(b);
  double pi2star = a[i];
  double t = directionUpper ? -b[i] : b[i];

  // calculate attained alpha and power
  auto g = [n1, n2, x1, x2, T, t,
            directionUpper](double p1, double p2)->double{
              NumericVector q1 = dbinom(x1, n1, p1, 0);
              NumericVector q2 = dbinom(x2, n2, p2, 0);

              double preject = 0;
              int i = 0;
              for (int y1=0; y1<=n1; y1++) {
                for (int y2=0; y2<=n2; y2++) {
                  if ((directionUpper && (T[i] >= t)) ||
                      ((!directionUpper) && (T[i] <= t))) {
                    preject += q1[y1]*q2[y2];
                  }
                  i++;
                }
              }

              return preject;
            };

  double power = g(pi1, pi2);

  DataFrame result;

  if (calculateAttainedAlpha) {
    auto h = [n1, n2, x1, x2, T, t, directionUpper,
              riskRatioH0](double p2)->double {
                double p1 = p2*riskRatioH0;
                NumericVector q1 = dbinom(x1, n1, p1, 0);
                NumericVector q2 = dbinom(x2, n2, p2, 0);

                double preject = 0;
                int i = 0;
                for (int y1=0; y1<=n1; y1++) {
                  for (int y2=0; y2<=n2; y2++) {
                    if ((directionUpper && (T[i] >= t)) ||
                        ((!directionUpper) && (T[i] <= t))) {
                      preject += q1[y1]*q2[y2];
                    }
                    i++;
                  }
                }

                return -preject;
              };

    for (i=0; i<K; i++) {
      double lower = pi2lower + i*delta;
      double upper = lower + delta;
      NumericVector out = mini(h, lower, upper, tol);
      a[i] = out[0];
      b[i] = out[1];
    }

    i = which_min(b);
    double attainedAlpha = -b[i];

    result = DataFrame::create(
      _["alpha"] = alpha,
      _["attainedAlpha"] = attainedAlpha,
      _["power"] = power,
      _["n"] = n,
      _["riskRatioH0"] = riskRatioH0,
      _["pi1"] = pi1,
      _["pi2"] = pi2,
      _["allocationRatioPlanned"] = allocationRatioPlanned,
      _["zstatRiskRatioBound"] = t,
      _["pi2star"] = pi2star);
  } else {
    result = DataFrame::create(
      _["alpha"] = alpha,
      _["power"] = power,
      _["n"] = n,
      _["riskRatioH0"] = riskRatioH0,
      _["pi1"] = pi1,
      _["pi2"] = pi2,
      _["allocationRatioPlanned"] = allocationRatioPlanned,
      _["zstatRiskRatioBound"] = t,
      _["pi2star"] = pi2star);
  }

  return result;
}


//' @title Sample size for exact unconditional test of risk ratio
//' @description Obtains the sample size given power for exact unconditional
//' test of risk ratio.
//'
//' @param beta The type II error.
//' @param riskRatioH0 The risk ratio under the null hypothesis.
//'   Defaults to 1.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The one-sided significance level.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified one-sided significance level.
//'
//' * \code{attainedAlpha}: The attained one-sided significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskRatioH0}: The risk ratio under the null hypothesis.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskRatioBound}: The critical value on the scale of
//'   score test statistic for risk ratio.
//'
//' * \code{pi2star}: The response probability in the control group
//'   at which the critical value of the test statistic is attained.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' samplesizeRiskRatioExact(beta = 0.2, riskRatioH0 = 0.833,
//'                          pi1 = 0.9, pi2 = 0.9, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
DataFrame samplesizeRiskRatioExact(
    const double beta = NA_REAL,
    const double riskRatioH0 = 1,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.025) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  NumericVector mr = remlRiskRatio(riskRatioH0, r, r*pi1, 1-r, (1-r)*pi2);
  double p1 = mr[0], p2 = mr[1];

  double v0 = (1-p1)/(r*p1) + (1-p2)/((1-r)*p2);
  double v1 = (1-pi1)/(r*pi1) + (1-pi2)/((1-r)*pi2);
  double z0 = R::qnorm(1-alpha, 0, 1, 1, 0);
  double z1 = R::qnorm(1-beta, 0, 1, 1, 0);
  double theta = fabs(log(pi1/pi2)- log(riskRatioH0));
  double n0 = pow(z0*sqrt(v0) + z1*sqrt(v1), 2)/pow(theta, 2);

  int n_lower = std::floor(n0), n_upper = std::ceil(10*n0), n;

  DataFrame a, b;
  a = powerRiskRatioExact(n0, riskRatioH0, pi1, pi2,
                          allocationRatioPlanned, alpha, 0);

  while (as<double>(a["power"]) >= 1-beta) {
    a = powerRiskRatioExact(--n_lower, riskRatioH0, pi1, pi2,
                            allocationRatioPlanned, alpha, 0);
  }

  for (n=n_lower+1; n<=n_upper; n++) {
    a = powerRiskRatioExact(n, riskRatioH0, pi1, pi2,
                            allocationRatioPlanned, alpha, 0);

    // ensure that the power exceeds 1-beta for the next 10 sample sizes
    if (as<double>(a["power"]) >= 1-beta) {
      int i;
      bool okay = 1;
      for (i=1; i<=10; i++) {
        b = powerRiskRatioExact(n+i, riskRatioH0, pi1, pi2,
                                allocationRatioPlanned, alpha, 0);
        if (as<double>(b["power"]) < 1-beta) {
          okay = 0;
          break;
        }
      }

      if (okay) {
        break;
      } else {
        n += i;
      }
    }
  }

  a = powerRiskRatioExact(n, riskRatioH0, pi1, pi2,
                          allocationRatioPlanned, alpha, 1);

  return a;
}


//' @title Power for exact unconditional test of equivalence in risk
//' difference
//' @description Obtains the power given sample size for exact unconditional
//' test of equivalence in risk difference.
//'
//' @param n The total sample size.
//' @param riskDiffLower The lower equivalence limit of risk difference.
//' @param riskDiffUpper The upper equivalence limit of risk difference.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified significance level for each of the two
//'   one-sided tests.
//'
//' * \code{attainedAlphaH10}: The attained significance level under H10
//'   if requested.
//'
//' * \code{attainedAlphaH20}: The attained significance level under H20
//'   if requested.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskDiffLower}: The lower equivalence limit of risk difference.
//'
//' * \code{riskDiffUpper}: The upper equivalence limit of risk difference.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{riskDiff}: The risk difference.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskDiffLower}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   lower equivalence limit.
//'
//' * \code{zstatRiskDiffUpper}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   upper equivalence limit.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' powerRiskDiffExactEquiv(
//'   n = 200, riskDiffLower = -0.2, riskDiffUpper = 0.2,
//'   pi1 = 0.775, pi2 = 0.775, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
DataFrame powerRiskDiffExactEquiv(
    const int n = NA_INTEGER,
    const double riskDiffLower = NA_REAL,
    const double riskDiffUpper = NA_REAL,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05,
    const bool calculateAttainedAlpha = 1) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  int n1 = std::round(n*r);
  int n2 = n - n1;
  IntegerVector x1 = seq(0, n1);
  IntegerVector x2 = seq(0, n2);
  int i, k = (n1+1)*(n2+1);

  // The set of values for the test statistic for the lower limit
  NumericVector T1(k);
  i = 0;
  for (int y1=0; y1<=n1; y1++) {
    for (int y2=0; y2<=n2; y2++) {
      NumericVector a = remlRiskDiff(riskDiffLower, n1, y1, n2, y2);
      double p1 = a[0], p2 = a[1];
      double p1hat = y1/(n1+0.0), p2hat = y2/(n2+0.0);
      double md = p1hat - p2hat - riskDiffLower;
      double mv = std::max(p1*(1-p1)/n1 + p2*(1-p2)/n2, 1.0e-8);
      T1[i++] = md/sqrt(mv);
    }
  }

  // sort T1 in ascending order
  IntegerVector order1 = seq(0, k-1);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    return T1[i] < T1[j];
  });

  NumericVector T1sorted = T1[order1];

  // identify the locations of the unique values of T1
  IntegerVector idx1(1,0);
  for (i=1; i<k; i++) {
    if (T1sorted[i] != T1sorted[i-1]) {
      idx1.push_back(i);
    }
  }

  NumericVector T1unique = T1sorted[idx1];

  int k1 = idx1.size();  // how many unique values of T1
  idx1.push_back(k);     // add the upper bound to the index

  // obtain the critical value of the test statistic T1
  auto f1 = [n1, n2, x1, x2, riskDiffLower, T1unique,
             k, k1, idx1, order1, alpha](double p2)->double {
               // table probabilities under H10
               NumericVector q1 = dbinom(x1, n1, p2 + riskDiffLower, 0);
               NumericVector q2 = dbinom(x2, n2, p2, 0);
               NumericVector q(k);

               int i, j, l;
               i = 0;
               for (int y1=0; y1<=n1; y1++) {
                 for (int y2=0; y2<=n2; y2++) {
                   q[i++] = q1[y1]*q2[y2];
                 }
               }

               NumericVector qsorted = q[order1];
               double aval;

               double s = 0;
               l = k-1;  // flattened index
               for (i=k1-1; i>=0; i--) {  // move down from the largest T1
                 for (j=idx1[i+1]-1; j>=idx1[i]; j--) {
                   s += qsorted[l--];
                 }
                 if (s > alpha) break;
               }

               // backtrack to obtain the critical value of T1
               if (i == k1-1) {
                 aval = T1unique[k1-1] + 1; // impossible to reject H10
               } else {
                 aval = T1unique[i+1];
               }

               return -aval;
             };

  // find the critical value for testing H10 independent of pi2
  double pi2lower1 = std::max(0.0, -riskDiffLower);
  double pi2upper1 = std::min(1.0, 1-riskDiffLower);

  int K = 100;
  NumericVector a1(K), b1(K);
  double delta1 = (pi2upper1 - pi2lower1)/K;
  double tol = 1.0e-6;
  for (i=0; i<K; i++) {
    double lower = pi2lower1 + i*delta1;
    double upper = lower + delta1;
    NumericVector out = mini(f1, lower, upper, tol);
    a1[i] = out[0];
    b1[i] = out[1];
  }

  i = which_min(b1);
  double t1 = -b1[i];


  // The set of values for the test statistic for the upper limit
  NumericVector T2(k);
  i = 0;
  for (int y1=0; y1<=n1; y1++) {
    for (int y2=0; y2<=n2; y2++) {
      NumericVector a = remlRiskDiff(riskDiffUpper, n1, y1, n2, y2);
      double p1 = a[0], p2 = a[1];
      double p1hat = y1/(n1+0.0), p2hat = y2/(n2+0.0);
      double md = p1hat - p2hat - riskDiffUpper;
      double mv = std::max(p1*(1-p1)/n1 + p2*(1-p2)/n2, 1.0e-8);
      T2[i++] = md/sqrt(mv);
    }
  }

  // find the critical value independent of pi2
  double pi2lower2 = std::max(0.0, -riskDiffUpper);
  double pi2upper2 = std::min(1.0, 1-riskDiffUpper);
  double t2;
  NumericVector a2(K);

  if (riskDiffLower == -riskDiffUpper) {
    t2 = -t1;
  } else {
    // sort T2 in ascending order
    IntegerVector order2 = seq(0, k-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      return T2[i] < T2[j];
    });

    NumericVector T2sorted = T2[order2];

    // identify the locations of the unique values of T2
    IntegerVector idx2(1,0);
    for (i=1; i<k; i++) {
      if (T2sorted[i] != T2sorted[i-1]) {
        idx2.push_back(i);
      }
    }

    NumericVector T2unique = T2sorted[idx2];

    int k2 = idx2.size();  // how many unique values of T2
    idx2.push_back(k);     // add the upper bound to the index

    // obtain the critical value of the test statistic
    auto f2 = [n1, n2, x1, x2, riskDiffUpper, T2unique,
               k, k2, idx2, order2, alpha](double p2)->double {
                 // table probabilities under H20
                 NumericVector q1 = dbinom(x1, n1, p2 + riskDiffUpper, 0);
                 NumericVector q2 = dbinom(x2, n2, p2, 0);
                 NumericVector q(k);

                 int i, j, l;
                 i = 0;
                 for (int y1=0; y1<=n1; y1++) {
                   for (int y2=0; y2<=n2; y2++) {
                     q[i++] = q1[y1]*q2[y2];
                   }
                 }

                 NumericVector qsorted = q[order2];
                 double aval;

                 double s = 0;
                 l = 0;
                 for (i=0; i<k2; i++) {  // move up from the smallest T2
                   for (j=idx2[i]; j<idx2[i+1]; j++) {
                     s += qsorted[l++];
                   }
                   if (s > alpha) break;
                 }

                 // backtrack to obtain the critical value of T2
                 if (i == 0) {
                   aval = T2unique[0] - 1; // impossible to reject H20
                 } else {
                   aval = T2unique[i-1];
                 }

                 return aval;
               };

    // find the critical value independent of pi2
    NumericVector b2(K);
    double delta2 = (pi2upper2 - pi2lower2)/K;
    for (i=0; i<K; i++) {
      double lower = pi2lower2 + i*delta2;
      double upper = lower + delta2;
      NumericVector out = mini(f2, lower, upper, tol);
      a2[i] = out[0];
      b2[i] = out[1];
    }

    i = which_min(b2);
    t2 = b2[i];
  }


  // calculate attained alpha
  auto g = [n1, n2, x1, x2, T1, T2, t1, t2](double p1, double p2)->double {
    NumericVector q1 = dbinom(x1, n1, p1, 0);
    NumericVector q2 = dbinom(x2, n2, p2, 0);

    double preject = 0;
    int i = 0;
    for (int y1=0; y1<=n1; y1++) {
      for (int y2=0; y2<=n2; y2++) {
        if ((T1[i] >= t1) && (T2[i] <= t2)) {
          preject += q1[y1]*q2[y2];
        }
        i++;
      }
    }

    return preject;
  };

  double power = g(pi1, pi2);

  DataFrame result;

  if (calculateAttainedAlpha) {
    double riskDiffH0;

    auto h = [n1, n2, x1, x2, T1, T2, t1, t2,
              &riskDiffH0](double p2)->double {
                double p1 = p2 + riskDiffH0;
                NumericVector q1 = dbinom(x1, n1, p1, 0);
                NumericVector q2 = dbinom(x2, n2, p2, 0);

                double preject = 0;
                int i = 0;
                for (int y1=0; y1<=n1; y1++) {
                  for (int y2=0; y2<=n2; y2++) {
                    if ((T1[i] >= t1) && (T2[i] <= t2)) {
                      preject += q1[y1]*q2[y2];
                    }
                    i++;
                  }
                }

                return -preject;
              };

    riskDiffH0 = riskDiffLower;
    for (i=0; i<K; i++) {
      double lower = pi2lower1 + i*delta1;
      double upper = lower + delta1;
      NumericVector out = mini(h, lower, upper, tol);
      b1[i] = out[1];
    }

    i = which_min(b1);
    double attainedAlphaH10 = -b1[i];

    double attainedAlphaH20;
    if (riskDiffLower == -riskDiffUpper) {
      attainedAlphaH20 = attainedAlphaH10;
    } else {
      NumericVector b2(K);
      double delta2 = (pi2upper2 - pi2lower2)/K;

      riskDiffH0 = riskDiffUpper;
      for (i=0; i<K; i++) {
        double lower = pi2lower2 + i*delta2;
        double upper = lower + delta2;
        NumericVector out = mini(h, lower, upper, tol);
        b2[i] = out[1];
      }

      i = which_min(b2);
      attainedAlphaH20 = -b2[i];
    }

    result = DataFrame::create(
      _["alpha"] = alpha,
      _["attainedAlphaH10"] = attainedAlphaH10,
      _["attainedAlphaH20"] = attainedAlphaH20,
      _["power"] = power,
      _["n"] = n,
      _["riskDiffLower"] = riskDiffLower,
      _["riskDiffUpper"] = riskDiffUpper,
      _["pi1"] = pi1,
      _["pi2"] = pi2,
      _["riskDiff"] = pi1 - pi2,
      _["allocationRatioPlanned"] = allocationRatioPlanned,
      _["zstatRiskDiffLower"] = t1,
      _["zstatRiskDiffUpper"] = t2);
  } else {
    result = DataFrame::create(
      _["alpha"] = alpha,
      _["power"] = power,
      _["n"] = n,
      _["riskDiffLower"] = riskDiffLower,
      _["riskDiffUpper"] = riskDiffUpper,
      _["pi1"] = pi1,
      _["pi2"] = pi2,
      _["riskDiff"] = pi1 - pi2,
      _["allocationRatioPlanned"] = allocationRatioPlanned,
      _["zstatRiskDiffLower"] = t1,
      _["zstatRiskDiffUpper"] = t2);
  }

  return result;
}


//' @title Sample size for exact unconditional test of equivalence in risk
//' difference
//' @description Obtains the sample size given power for exact unconditional
//' test of equivalence in risk difference.
//'
//' @param beta The type II error.
//' @param riskDiffLower The lower equivalence limit of risk difference.
//' @param riskDiffUpper The upper equivalence limit of risk difference.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified significance level for each of the two
//'   one-sided tests.
//'
//' * \code{attainedAlpha}: The attained significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskDiffLower}: The lower equivalence limit of risk difference.
//'
//' * \code{riskDiffUpper}: The upper equivalence limit of risk difference.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{riskDiff}: The risk difference.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskDiffLower}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   lower equivalence limit.
//'
//' * \code{zstatRiskDiffUpper}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   upper equivalence limit.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' samplesizeRiskDiffExactEquiv(
//'   beta = 0.2, riskDiffLower = -0.3, riskDiffUpper = 0.3,
//'   pi1 = 0.85, pi2 = 0.85, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
DataFrame samplesizeRiskDiffExactEquiv(
    const double beta = NA_REAL,
    const double riskDiffLower = NA_REAL,
    const double riskDiffUpper = NA_REAL,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  double theta = std::min(pi1 - pi2 - riskDiffLower,
                          riskDiffUpper - (pi1 - pi2));

  double v1 = pi1*(1-pi1)/r + pi2*(1-pi2)/(1-r);
  double z0 = R::qnorm(1-alpha, 0, 1, 1, 0);
  double z1 = R::qnorm(1-beta, 0, 1, 1, 0);
  double n0 = pow(z0 + z1, 2)*v1/pow(theta, 2);

  int n_lower = std::floor(n0), n_upper = std::ceil(10*n0), n;

  DataFrame a, b;
  a = powerRiskDiffExactEquiv(n0, riskDiffLower, riskDiffUpper,
                              pi1, pi2, allocationRatioPlanned, alpha, 0);

  while (as<double>(a["power"]) >= 1-beta) {
    a = powerRiskDiffExactEquiv(--n_lower, riskDiffLower, riskDiffUpper,
                                pi1, pi2, allocationRatioPlanned, alpha, 0);
  }

  for (n=n_lower+1; n<=n_upper; n++) {
    a = powerRiskDiffExactEquiv(n, riskDiffLower, riskDiffUpper,
                                pi1, pi2, allocationRatioPlanned, alpha, 0);

    // ensure that the power exceeds 1-beta for the next 5 sample sizes
    if (as<double>(a["power"]) >= 1-beta) {
      int i;
      bool okay = 1;
      for (i=1; i<=5; i++) {
        b = powerRiskDiffExactEquiv(n+i, riskDiffLower, riskDiffUpper,
                                    pi1, pi2, allocationRatioPlanned, alpha, 0);
        if (as<double>(b["power"]) < 1-beta) {
          okay = 0;
          break;
        }
      }

      if (okay) {
        break;
      } else {
        n += i;
      }
    }
  }

  a = powerRiskDiffExactEquiv(n, riskDiffLower, riskDiffUpper,
                              pi1, pi2, allocationRatioPlanned, alpha, 1);

  return a;
}


//' @title Power for exact unconditional test of equivalence in risk
//' ratio
//' @description Obtains the power given sample size for exact unconditional
//' test of equivalence in risk ratio.
//'
//' @param n The total sample size.
//' @param riskRatioLower The lower equivalence limit of risk ratio.
//' @param riskRatioUpper The upper equivalence limit of risk ratio.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @param calculateAttainedAlpha Whether to calculate the attained alpha.
//'
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified significance level for each of the two
//'   one-sided tests.
//'
//' * \code{attainedAlphaH10}: The attained significance level under H10
//'   if requested.
//'
//' * \code{attainedAlphaH20}: The attained significance level under H20
//'   if requested.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskRatioLower}: The lower equivalence limit of risk ratio.
//'
//' * \code{riskRatioUpper}: The upper equivalence limit of risk ratio.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{riskRatio}: The risk ratio.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskRatioLower}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   lower equivalence limit.
//'
//' * \code{zstatRiskRatioUpper}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   upper equivalence limit.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' powerRiskRatioExactEquiv(
//'   n = 200, riskRatioLower = 0.8, riskRatioUpper = 1.25,
//'   pi1 = 0.775, pi2 = 0.775, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
DataFrame powerRiskRatioExactEquiv(
    const int n = NA_INTEGER,
    const double riskRatioLower = NA_REAL,
    const double riskRatioUpper = NA_REAL,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05,
    const bool calculateAttainedAlpha = 1) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  int n1 = std::round(n*r);
  int n2 = n - n1;
  IntegerVector x1 = seq(0, n1);
  IntegerVector x2 = seq(0, n2);
  int i, k = (n1+1)*(n2+1);

  // The set of values for the test statistic for the lower limit
  NumericVector T1(k);
  i = 0;
  for (int y1=0; y1<=n1; y1++) {
    for (int y2=0; y2<=n2; y2++) {
      NumericVector a = remlRiskRatio(riskRatioLower, n1, y1, n2, y2);
      double p1 = a[0], p2 = a[1];
      double p1hat = y1/(n1+0.0), p2hat = y2/(n2+0.0);
      double md = p1hat - p2hat*riskRatioLower;
      double mv = std::max(p1*(1-p1)/n1 +
                           pow(riskRatioLower,2)*p2*(1-p2)/n2, 1.0e-8);
      T1[i++] = md/sqrt(mv);
    }
  }

  // sort T1 in ascending order
  IntegerVector order1 = seq(0, k-1);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    return T1[i] < T1[j];
  });

  NumericVector T1sorted = T1[order1];

  // identify the locations of the unique values of T1
  IntegerVector idx1(1,0);
  for (i=1; i<k; i++) {
    if (T1sorted[i] != T1sorted[i-1]) {
      idx1.push_back(i);
    }
  }

  NumericVector T1unique = T1sorted[idx1];

  int k1 = idx1.size();  // how many unique values of T1
  idx1.push_back(k);     // add the upper bound to the index

  // obtain the critical value of the test statistic T1
  auto f1 = [n1, n2, x1, x2, riskRatioLower, T1unique,
             k, k1, idx1, order1, alpha](double p2)->double {
               // table probabilities under H10
               NumericVector q1 = dbinom(x1, n1, p2*riskRatioLower, 0);
               NumericVector q2 = dbinom(x2, n2, p2, 0);
               NumericVector q(k);

               int i, j, l;
               i = 0;
               for (int y1=0; y1<=n1; y1++) {
                 for (int y2=0; y2<=n2; y2++) {
                   q[i++] = q1[y1]*q2[y2];
                 }
               }

               NumericVector qsorted = q[order1];
               double aval;

               double s = 0;
               l = k-1;  // flattened index
               for (i=k1-1; i>=0; i--) {  // move down from the largest T1
                 for (j=idx1[i+1]-1; j>=idx1[i]; j--) {
                   s += qsorted[l--];
                 }
                 if (s > alpha) break;
               }

               // backtrack to obtain the critical value of T1
               if (i == k1-1) {
                 aval = T1unique[k1-1] + 1; // impossible to reject H10
               } else {
                 aval = T1unique[i+1];
               }

               return -aval;
             };

  // find the critical value for testing H10 independent of pi2
  double pi2lower1 = 0;
  double pi2upper1 = std::min(1.0, 1/riskRatioLower);

  int K = 100;
  NumericVector a1(K), b1(K);
  double delta1 = (pi2upper1 - pi2lower1)/K;
  double tol = 1.0e-6;
  for (i=0; i<K; i++) {
    double lower = pi2lower1 + i*delta1;
    double upper = lower + delta1;
    NumericVector out = mini(f1, lower, upper, tol);
    a1[i] = out[0];
    b1[i] = out[1];
  }

  i = which_min(b1);
  double t1 = -b1[i];


  // The set of values for the test statistic for the upper limit
  NumericVector T2(k);
  i = 0;
  for (int y1=0; y1<=n1; y1++) {
    for (int y2=0; y2<=n2; y2++) {
      NumericVector a = remlRiskRatio(riskRatioUpper, n1, y1, n2, y2);
      double p1 = a[0], p2 = a[1];
      double p1hat = y1/(n1+0.0), p2hat = y2/(n2+0.0);
      double md = p1hat - p2hat*riskRatioUpper;
      double mv = std::max(p1*(1-p1)/n1 +
                           pow(riskRatioUpper,2)*p2*(1-p2)/n2, 1.0e-8);
      T2[i++] = md/sqrt(mv);
    }
  }

  // find the critical value independent of pi2
  double pi2lower2 = 0;
  double pi2upper2 =  std::min(1.0, 1/riskRatioUpper);
  double t2;
  NumericVector a2(K);

  if (fabs(riskRatioLower*riskRatioUpper - 1) < 1.0e-8) {
    t2 = -t1;
  } else {
    // sort T2 in ascending order
    IntegerVector order2 = seq(0, k-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      return T2[i] < T2[j];
    });

    NumericVector T2sorted = T2[order2];

    // identify the locations of the unique values of T2
    IntegerVector idx2(1,0);
    for (i=1; i<k; i++) {
      if (T2sorted[i] != T2sorted[i-1]) {
        idx2.push_back(i);
      }
    }

    NumericVector T2unique = T2sorted[idx2];

    int k2 = idx2.size();  // how many unique values of T2
    idx2.push_back(k);     // add the upper bound to the index

    // obtain the critical value of the test statistic
    auto f2 = [n1, n2, x1, x2, riskRatioUpper, T2unique,
               k, k2, idx2, order2, alpha](double p2)->double {
                 // table probabilities under H20
                 NumericVector q1 = dbinom(x1, n1, p2*riskRatioUpper, 0);
                 NumericVector q2 = dbinom(x2, n2, p2, 0);
                 NumericVector q(k);

                 int i, j, l;
                 i = 0;
                 for (int y1=0; y1<=n1; y1++) {
                   for (int y2=0; y2<=n2; y2++) {
                     q[i++] = q1[y1]*q2[y2];
                   }
                 }

                 NumericVector qsorted = q[order2];

                 double s = 0;
                 l = 0;
                 for (i=0; i<k2; i++) {  // move up from the smallest T2
                   for (j=idx2[i]; j<idx2[i+1]; j++) {
                     s += qsorted[l++];
                   }
                   if (s > alpha) break;
                 }

                 // backtrack to obtain the critical value of T2
                 if (i == 0) {
                   return T2unique[0] - 1; // impossible to reject H20
                 } else {
                   return T2unique[i-1];
                 }
               };

    // find the critical value independent of pi2
    NumericVector b2(K);
    double delta2 = (pi2upper2 - pi2lower2)/K;
    for (i=0; i<K; i++) {
      double lower = pi2lower2 + i*delta2;
      double upper = lower + delta2;
      NumericVector out = mini(f2, lower, upper, tol);
      a2[i] = out[0];
      b2[i] = out[1];
    }

    i = which_min(b2);
    t2 = b2[i];
  }


  // calculate attained alpha
  auto g = [n1, n2, x1, x2, T1, T2, t1, t2](double p1, double p2)->double {
    NumericVector q1 = dbinom(x1, n1, p1, 0);
    NumericVector q2 = dbinom(x2, n2, p2, 0);

    double preject = 0;
    int i = 0;
    for (int y1=0; y1<=n1; y1++) {
      for (int y2=0; y2<=n2; y2++) {
        if ((T1[i] >= t1) && (T2[i] <= t2)) {
          preject += q1[y1]*q2[y2];
        }
        i++;
      }
    }

    return preject;
  };

  double power = g(pi1, pi2);

  DataFrame result;

  if (calculateAttainedAlpha) {
    double riskRatioH0;

    auto h = [n1, n2, x1, x2, T1, T2, t1, t2,
              &riskRatioH0](double p2)->double {
                double p1 = p2*riskRatioH0;
                NumericVector q1 = dbinom(x1, n1, p1, 0);
                NumericVector q2 = dbinom(x2, n2, p2, 0);

                double preject = 0;
                int i = 0;
                for (int y1=0; y1<=n1; y1++) {
                  for (int y2=0; y2<=n2; y2++) {
                    if ((T1[i] >= t1) && (T2[i] <= t2)) {
                      preject += q1[y1]*q2[y2];
                    }
                    i++;
                  }
                }

                return -preject;
              };

    riskRatioH0 = riskRatioLower;
    for (i=0; i<K; i++) {
      double lower = pi2lower1 + i*delta1;
      double upper = lower + delta1;
      NumericVector out = mini(h, lower, upper, tol);
      b1[i] = out[1];
    }

    i = which_min(b1);
    double attainedAlphaH10 = -b1[i];

    double attainedAlphaH20;
    if (fabs(riskRatioLower*riskRatioUpper - 1) < 1.0e-8) {
      attainedAlphaH20 = attainedAlphaH10;
    } else {
      NumericVector b2(K);
      double delta2 = (pi2upper2 - pi2lower2)/K;

      riskRatioH0 = riskRatioUpper;
      for (i=0; i<K; i++) {
        double lower = pi2lower2 + i*delta2;
        double upper = lower + delta2;
        NumericVector out = mini(h, lower, upper, tol);
        b2[i] = out[1];
      }

      i = which_min(b2);
      attainedAlphaH20 = -b2[i];
    }

    result = DataFrame::create(
      _["alpha"] = alpha,
      _["attainedAlphaH10"] = attainedAlphaH10,
      _["attainedAlphaH20"] = attainedAlphaH20,
      _["power"] = power,
      _["n"] = n,
      _["riskRatioLower"] = riskRatioLower,
      _["riskRatioUpper"] = riskRatioUpper,
      _["pi1"] = pi1,
      _["pi2"] = pi2,
      _["riskRatio"] = pi1/pi2,
      _["allocationRatioPlanned"] = allocationRatioPlanned,
      _["zstatRiskRatioLower"] = t1,
      _["zstatRiskRatioUpper"] = t2);
  } else {
    result = DataFrame::create(
      _["alpha"] = alpha,
      _["power"] = power,
      _["n"] = n,
      _["riskRatioLower"] = riskRatioLower,
      _["riskRatioUpper"] = riskRatioUpper,
      _["pi1"] = pi1,
      _["pi2"] = pi2,
      _["riskRatio"] = pi1/pi2,
      _["allocationRatioPlanned"] = allocationRatioPlanned,
      _["zstatRiskRatioLower"] = t1,
      _["zstatRiskRatioUpper"] = t2);
  }

  return result;
}


//' @title Sample size for exact unconditional test of equivalence in risk
//' ratio
//' @description Obtains the sample size given power for exact unconditional
//' test of equivalence in risk ratio.
//'
//' @param beta The type II error.
//' @param riskRatioLower The lower equivalence limit of risk ratio.
//' @param riskRatioUpper The upper equivalence limit of risk ratio.
//' @param pi1 The assumed probability for the active treatment group.
//' @param pi2 The assumed probability for the control group.
//' @param allocationRatioPlanned Allocation ratio for the active treatment
//'   versus control. Defaults to 1 for equal randomization.
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @return A data frame with the following variables:
//'
//' * \code{alpha}: The specified significance level for each of the two
//'   one-sided tests.
//'
//' * \code{attainedAlpha}: The attained significance level.
//'
//' * \code{power}: The power.
//'
//' * \code{n}: The sample size.
//'
//' * \code{riskRatioLower}: The lower equivalence limit of risk ratio.
//'
//' * \code{riskRatioUpper}: The upper equivalence limit of risk ratio.
//'
//' * \code{pi1}: The assumed probability for the active treatment group.
//'
//' * \code{pi2}: The assumed probability for the control group.
//'
//' * \code{riskRatio}: The risk ratio.
//'
//' * \code{allocationRatioPlanned}: Allocation ratio for the active
//'   treatment versus control.
//'
//' * \code{zstatRiskRatioLower}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   lower equivalence limit.
//'
//' * \code{zstatRiskRatioUpper}: The efficacy boundaries on the
//'   z-test statistic scale for the one-sided null hypothesis on the
//'   upper equivalence limit.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' samplesizeRiskRatioExactEquiv(
//'   beta = 0.2, riskRatioLower = 0.7, riskRatioUpper = 1/0.7,
//'   pi1 = 0.85, pi2 = 0.85, alpha = 0.05)
//'
//' @export
// [[Rcpp::export]]
DataFrame samplesizeRiskRatioExactEquiv(
    const double beta = NA_REAL,
    const double riskRatioLower = NA_REAL,
    const double riskRatioUpper = NA_REAL,
    const double pi1 = NA_REAL,
    const double pi2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const double alpha = 0.05) {

  double r = allocationRatioPlanned/(1 + allocationRatioPlanned);
  double theta = std::min(log(pi1/pi2) - log(riskRatioLower),
                          log(riskRatioUpper) - log(pi1/pi2));

  double v1 = (1-pi1)/(r*pi1) + (1-pi2)/((1-r)*pi2);
  double z0 = R::qnorm(1-alpha, 0, 1, 1, 0);
  double z1 = R::qnorm(1-beta, 0, 1, 1, 0);
  double n0 = pow(z0 + z1, 2)*v1/pow(theta, 2);

  int n_lower = std::floor(n0), n_upper = std::ceil(10*n0), n;

  DataFrame a, b;
  a = powerRiskRatioExactEquiv(n0, riskRatioLower, riskRatioUpper,
                               pi1, pi2, allocationRatioPlanned, alpha, 0);

  while (as<double>(a["power"]) >= 1-beta) {
    a = powerRiskRatioExactEquiv(--n_lower, riskRatioLower, riskRatioUpper,
                                 pi1, pi2, allocationRatioPlanned, alpha, 0);
  }

  for (n=n_lower+1; n<=n_upper; n++) {
    a = powerRiskRatioExactEquiv(n, riskRatioLower, riskRatioUpper,
                                 pi1, pi2, allocationRatioPlanned, alpha, 0);

    // ensure that the power exceeds 1-beta for the next 5 sample sizes
    if (as<double>(a["power"]) >= 1-beta) {
      int i;
      bool okay = 1;
      for (i=1; i<=5; i++) {
        b = powerRiskRatioExactEquiv(n+i, riskRatioLower, riskRatioUpper,
                                     pi1, pi2, allocationRatioPlanned,
                                     alpha, 0);
        if (as<double>(b["power"]) < 1-beta) {
          okay = 0;
          break;
        }
      }

      if (okay) {
        break;
      } else {
        n += i;
      }
    }
  }

  a = powerRiskRatioExactEquiv(n, riskRatioLower, riskRatioUpper,
                               pi1, pi2, allocationRatioPlanned, alpha, 1);

  return a;
}



//' @title P-value for exact unconditional test of risk difference
//' @description Obtains the p-value for exact unconditional
//' test of risk difference.
//'
//' @param riskDiffH0 The risk difference under the null hypothesis.
//'   Defaults to 0.
//' @param directionUpper Whether larger values represent better
//'   responses.
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//'
//' @return A data frame containing the following variables:
//'
//' * \code{riskDiffH0}: The risk difference under the null hypothesis.
//'
//' * \code{directionUpper}: Whether larger values represent better
//'   responses.
//'
//' * \code{riskDiff}: The observed risk difference.
//'
//' * \code{zstat}: The observed value of the Z test statistic.
//'
//' * \code{pvalue}: The one-sided p-value for the unconditional exact test.
//'
//' * \code{pi2star}: The value of pi2 that yields the p-value.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' riskDiffExactPValue(riskDiffH0 = 0, directionUpper = 1,
//'                     n1 = 68, y1 = 2, n2 = 65, y2 = 1)
//'
//' @export
// [[Rcpp::export]]
DataFrame riskDiffExactPValue(
    const double riskDiffH0 = 0,
    bool directionUpper = 1,
    const int n1 = NA_INTEGER,
    const int y1 = NA_INTEGER,
    const int n2 = NA_INTEGER,
    const int y2 = NA_INTEGER) {

  int i, k = (n1+1)*(n2+1);

  // The set of values for the test statistic
  NumericVector T(k);
  i = 0;
  for (int x1=0; x1<=n1; x1++) {
    for (int x2=0; x2<=n2; x2++) {
      NumericVector a = remlRiskDiff(riskDiffH0, n1, x1, n2, x2);
      double p1 = a[0], p2 = a[1];
      double p1hat = x1/(n1+0.0), p2hat = x2/(n2+0.0);
      double md = p1hat - p2hat - riskDiffH0;
      double mv = std::max(p1*(1-p1)/n1 + p2*(1-p2)/n2, 1.0e-8);
      T[i++] = md/sqrt(mv);
    }
  }

  double Tobs = T[y1*(n2+1)+y2];
  double riskDiff = y1/(n1+0.0) - y2/(n2+0.0);

  IntegerVector x1 = seq(0, n1);
  IntegerVector x2 = seq(0, n2);

  auto f = [n1, n2, x1, x2, riskDiffH0, directionUpper,
            T, k, Tobs](double p2)->double {
              NumericVector q1 = dbinom(x1, n1, p2 + riskDiffH0, 0);
              NumericVector q2 = dbinom(x2, n2, p2, 0);
              NumericVector q(k);

              int i = 0;
              for (int y1=0; y1<=n1; y1++) {
                for (int y2=0; y2<=n2; y2++) {
                  q[i++] = q1[y1]*q2[y2];
                }
              }

              double s = 0.0;
              for (i=0; i<k; i++) {
                if ((2*directionUpper - 1)*(T[i] - Tobs) >= 0.0) {
                  s += q[i];
                }
              }

              return -s;
            };


  double pi2lower = std::max(0.0, -riskDiffH0);
  double pi2upper = std::min(1.0, 1-riskDiffH0);

  int K = 100;
  NumericVector a(K), b(K);
  double delta = (pi2upper - pi2lower)/K;
  double tol = 1.0e-6;
  for (i=0; i<K; i++) {
    double lower = pi2lower + i*delta;
    double upper = lower + delta;
    NumericVector out = mini(f, lower, upper, tol);
    a[i] = out[0];
    b[i] = out[1];
  }

  i = which_min(b);
  double pi2star = a[i];
  double pvalue = -b[i];

  DataFrame result = DataFrame::create(
    _["riskDiffH0"] = riskDiffH0,
    _["directionUpper"] = directionUpper,
    _["riskDiff"] = riskDiff,
    _["zstat"] = Tobs,
    _["pvalue"] = pvalue,
    _["pi2star"] = pi2star);

  return result;
}



//' @title Exact unconditional confidence interval for risk difference
//' @description Obtains the exact unconditional confidence interval for
//' risk difference based on the standardized score statistic.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @return A data frame containing the following variables:
//'
//' * \code{scale}: The scale of treatment effect.
//'
//' * \code{estimate}: The point estimate.
//'
//' * \code{lower}: The lower limit of the confidence interval.
//'
//' * \code{upper}: The upper limit of the confidence interval.
//'
//' * \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' riskDiffExactCI(n1 = 68, y1 = 2, n2 = 65, y2 = 1, cilevel = 0.95)
//'
//' @export
// [[Rcpp::export]]
DataFrame riskDiffExactCI(
    const int n1 = NA_INTEGER,
    const int y1 = NA_INTEGER,
    const int n2 = NA_INTEGER,
    const int y2 = NA_INTEGER,
    const double cilevel = 0.95) {

  double estimate = y1/(n1+0.0) - y2/(n2+0.0);
  double alpha = 1.0 - cilevel; // two-sided alpha

  auto f1 = [n1, y1, n2, y2, alpha](double riskDiff)->double {
    DataFrame df = riskDiffExactPValue(riskDiff, 1, n1, y1, n2, y2);
    double pvalue = df["pvalue"];
    return pvalue - alpha/2;
  };

  auto f2 = [n1, y1, n2, y2, alpha](double riskDiff)->double {
    DataFrame df = riskDiffExactPValue(riskDiff, 0, n1, y1, n2, y2);
    double pvalue = df["pvalue"];
    return pvalue - alpha/2;
  };

  double lower = brent(f1, -1.0, estimate, 1.0e-6);
  double upper = brent(f2, estimate, 1.0, 1.0e-6);

  DataFrame result = DataFrame::create(
    _["scale"] = "risk difference",
    _["estimate"] = estimate,
    _["lower"] = lower,
    _["upper"] = upper,
    _["cilevel"] = cilevel);

  return result;
}



//' @title P-value for exact unconditional test of risk ratio
//' @description Obtains the p-value for exact unconditional
//' test of risk ratio.
//'
//' @param riskRatioH0 The risk ratio under the null hypothesis.
//'   Defaults to 1.
//' @param directionUpper Whether larger values represent better
//'   responses.
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//'
//' @return A data frame containing the following variables:
//'
//' * \code{riskRatioH0}: The risk ratio under the null hypothesis.
//'
//' * \code{directionUpper}: Whether larger values represent better
//'   responses.
//'
//' * \code{riskRatio}: The observed risk ratio.
//'
//' * \code{zstat}: The observed value of the Z test statistic.
//'
//' * \code{pvalue}: The one-sided p-value for the unconditional exact test.
//'
//' * \code{pi2star}: The value of pi2 that yields the p-value.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' riskRatioExactPValue(riskRatioH0 = 1, directionUpper = 1,
//'                      n1 = 68, y1 = 2, n2 = 65, y2 = 1)
//'
//' @export
// [[Rcpp::export]]
DataFrame riskRatioExactPValue(
    const double riskRatioH0 = 1,
    bool directionUpper = 1,
    const int n1 = NA_INTEGER,
    const int y1 = NA_INTEGER,
    const int n2 = NA_INTEGER,
    const int y2 = NA_INTEGER) {

  int i, k = (n1+1)*(n2+1);

  // The set of values for the test statistic
  NumericVector T(k);
  i = 0;
  for (int x1=0; x1<=n1; x1++) {
    for (int x2=0; x2<=n2; x2++) {
      NumericVector a = remlRiskRatio(riskRatioH0, n1, x1, n2, x2);
      double p1 = a[0], p2 = a[1];
      double p1hat = x1/(n1+0.0), p2hat = x2/(n2+0.0);
      double md = p1hat - p2hat*riskRatioH0;
      double mv = std::max(p1*(1-p1)/n1 +
                           pow(riskRatioH0, 2)*p2*(1-p2)/n2, 1.0e-8);
      T[i++] = md/sqrt(mv);
    }
  }

  double Tobs = T[y1*(n2+1)+y2];
  double riskRatio = y2 > 0 ? (y1/(n1+0.0)) / (y2/(n2+0.0)) : R_PosInf;

  IntegerVector x1 = seq(0, n1);
  IntegerVector x2 = seq(0, n2);

  auto f = [n1, n2, x1, x2, riskRatioH0, directionUpper,
            T, k, Tobs](double p2)->double {
              NumericVector q1 = dbinom(x1, n1, p2*riskRatioH0, 0);
              NumericVector q2 = dbinom(x2, n2, p2, 0);
              NumericVector q(k);

              int i = 0;
              for (int y1=0; y1<=n1; y1++) {
                for (int y2=0; y2<=n2; y2++) {
                  q[i++] = q1[y1]*q2[y2];
                }
              }

              double s = 0.0;
              for (i=0; i<k; i++) {
                if ((2*directionUpper - 1)*(T[i] - Tobs) >= 0.0) {
                  s += q[i];
                }
              }

              return -s;
            };


  double pi2lower = 0.0;
  double pi2upper = std::min(1.0, 1.0/riskRatioH0);

  int K = 100;
  NumericVector a(K), b(K);
  double delta = (pi2upper - pi2lower)/K;
  double tol = 1.0e-6;
  for (i=0; i<K; i++) {
    double lower = pi2lower + i*delta;
    double upper = lower + delta;
    NumericVector out = mini(f, lower, upper, tol);
    a[i] = out[0];
    b[i] = out[1];
  }

  i = which_min(b);
  double pi2star = a[i];
  double pvalue = -b[i];

  DataFrame result = DataFrame::create(
    _["riskRatioH0"] = riskRatioH0,
    _["directionUpper"] = directionUpper,
    _["riskRatio"] = riskRatio,
    _["zstat"] = Tobs,
    _["pvalue"] = pvalue,
    _["pi2star"] = pi2star);

  return result;
}


//' @title Exact unconditional confidence interval for risk ratio
//' @description Obtains the exact unconditional confidence interval for
//' risk ratio based on the standardized score statistic.
//'
//' @param n1 The sample size for the active treatment group.
//' @param y1 The number of responses for the active treatment group.
//' @param n2 The sample size for the control group.
//' @param y2 The number of responses for the control group.
//' @param cilevel The confidence interval level.
//'
//' @return A data frame containing the following variables:
//'
//' * \code{scale}: The scale of treatment effect.
//'
//' * \code{estimate}: The point estimate.
//'
//' * \code{lower}: The lower limit of the confidence interval.
//'
//' * \code{upper}: The upper limit of the confidence interval.
//'
//' * \code{cilevel}: The confidence interval level.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' riskRatioExactCI(n1 = 68, y1 = 2, n2 = 65, y2 = 1, cilevel = 0.95)
//'
//' @export
// [[Rcpp::export]]
DataFrame riskRatioExactCI(
    const int n1 = NA_INTEGER,
    const int y1 = NA_INTEGER,
    const int n2 = NA_INTEGER,
    const int y2 = NA_INTEGER,
    const double cilevel = 0.95) {

  double estimate = y2 > 0 ? (y1/(n1+0.0)) / (y2/(n2+0.0)) : R_PosInf;
  double alpha = 1.0 - cilevel; // two-sided alpha

  auto f1 = [n1, y1, n2, y2, alpha](double riskRatio)->double {
    DataFrame df = riskRatioExactPValue(riskRatio, 1, n1, y1, n2, y2);
    double pvalue = df["pvalue"];
    return pvalue - alpha/2;
  };

  auto f2 = [n1, y1, n2, y2, alpha](double riskRatio)->double {
    DataFrame df = riskRatioExactPValue(riskRatio, 0, n1, y1, n2, y2);
    double pvalue = df["pvalue"];
    return pvalue - alpha/2;
  };

  double lower = brent(f1, 0.001, estimate, 1.0e-6);
  double upper = brent(f2, estimate, 1000.0, 1.0e-6);

  DataFrame result = DataFrame::create(
    _["scale"] = "risk difference",
    _["estimate"] = estimate,
    _["lower"] = lower,
    _["upper"] = upper,
    _["cilevel"] = cilevel);

  return result;
}


