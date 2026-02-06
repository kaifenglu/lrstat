#include <Rcpp.h>

#include "utilities.h"
#include "dataframe_list.h"

#include <algorithm>   // min_element, max_element, find, sort, adjacent_find,
// any_of, fill, distance, min, for_each, transform, tolower
#include <cctype>      // tolower
#include <cmath>       // fabs, isnan
#include <numeric>     // accumulate
#include <stdexcept>   // invalid_argument
#include <string>      // string
#include <vector>      // vector


ListCpp updateGraphcpp(const std::vector<double>& w,
                       const FlatMatrix& G,
                       const std::vector<int>& I,
                       const int j) {
  int m = static_cast<int>(w.size());

  // Check matrix dimension
  if (G.nrow != m || G.ncol != m) {
    throw std::invalid_argument("Invalid dimension for G.");
  }

  // Check if elements of I are between 1 and m
  int minI = *std::min_element(I.begin(), I.end());
  int maxI = *std::max_element(I.begin(), I.end());
  if (minI < 1 || maxI > m) {
    throw std::invalid_argument("Elements of I must be integers between 1 and m.");
  }

  // Check for duplicates in I
  std::vector<int> I_sorted = I;
  std::sort(I_sorted.begin(), I_sorted.end());
  auto it = std::adjacent_find(I_sorted.begin(), I_sorted.end());
  if (it != I_sorted.end()) {
    throw std::invalid_argument("The index set I must not contain duplicates.");
  }

  // Check if j is in I
  if (std::find(I.begin(), I.end(), j) == I.end()) {
    throw std::invalid_argument("j must be in I.");
  }

  // Convert to 0-based indexing
  int j1 = j - 1;
  std::vector<int> I1(I.size());
  for (std::size_t i = 0; i < I.size(); ++i) {
    I1[i] = I[i] - 1;
  }

  // Create r vector: 1 for all, 0 for elements in I1, then set r[j1] = 1
  std::vector<unsigned char> r(m, 1);
  for (int idx : I1) {
    r[idx] = 0;
  }
  r[j1] = 1;

  // Update weights
  std::vector<double> wx = w;  // copy w
  for (int l = 0; l < m; ++l) {
    if (r[l] == 0) {
      wx[l] = wx[l] + wx[j1] * G(j1, l);
    }
  }
  wx[j1] = 0.0;

  // Update transition matrix
  FlatMatrix g(m, m);
  for (int k = 0; k < m; ++k) {
    for (int l = 0; l < m; ++l) {
      if (r[l] == 0) {
        if ((r[k] == 0) && (l != k) && (G(l, j1) * G(j1, l) < 1.0 - 1.0e-12)) {
          g(l, k) = (G(l, k) + G(l, j1) * G(j1, k)) / (1.0 - G(l, j1) * G(j1, l));
        }
      }
    }
  }

  // Create I_new: elements of I that are not equal to j
  std::vector<int> I_new;
  I_new.reserve(I.size() - 1);
  for (int val : I) {
    if (val != j) {
      I_new.push_back(val);
    }
  }

  // Create result list
  ListCpp result;
  result.push_back(std::move(wx), "w");
  result.push_back(std::move(g), "G");
  result.push_back(std::move(I_new), "I");

  return result;
}

//' @title Update Graph for Graphical Approaches
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
Rcpp::List updateGraph(const Rcpp::NumericVector& w,
                       const Rcpp::NumericMatrix& G,
                       const Rcpp::IntegerVector& I,
                       const int j) {
  auto w1 = Rcpp::as<std::vector<double>>(w);
  auto I1 = Rcpp::as<std::vector<int>>(I);
  FlatMatrix G1 = flatmatrix_from_Rmatrix(G);
  ListCpp result = updateGraphcpp(w1, G1, I1, j);
  return Rcpp::wrap(result);
}


FlatMatrix fadjpboncpp1(const std::vector<double>& w,
                        const FlatMatrix& G,
                        const FlatMatrix& p) {

  int m = static_cast<int>(w.size());
  int iters = p.nrow;
  double pmax; // running maximum adjusted p-value

  FlatMatrix padj(iters, m);  // adjusted p-values
  std::vector<unsigned char> r(m);  // rejection indicators
  std::vector<double> pvalues(m);  // raw p-values
  std::vector<double> q(m);  // ratios of raw p-values over weights

  std::vector<double> wx(m);    // dynamic weights
  FlatMatrix g(m, m);  // dynamic transition matrix
  FlatMatrix g1(m, m);  // temporary transition matrix

  // Validation: w must be nonnegative
  if (std::any_of(w.begin(), w.end(), [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("w must be nonnegative");
  }

  // Validation: w must sum to 1
  double sum_w = std::accumulate(w.begin(), w.end(), 0.0);
  if (std::fabs(sum_w - 1.0) > 1.0e-12) {
    throw std::invalid_argument("w must sum to 1");
  }

  // Validation: G dimension
  if (G.nrow != m || G.ncol != m) {
    throw std::invalid_argument("Invalid dimension for G");
  }

  // Validation: G must be nonnegative
  for (int i = 0; i < m * m; ++i) {
    if (G.data[i] < 0.0) {
      throw std::invalid_argument("G must be nonnegative");
    }
  }

  // Validation: Row sums of G must be <= 1
  for (int i = 0; i < m; ++i) {
    double rowsum = 0.0;
    for (int j = 0; j < m; ++j) {
      rowsum += G(i, j);
    }
    if (rowsum > 1.0 + 1.0e-8) {
      throw std::invalid_argument("Row sums of G must be less than or equal to 1");
    }
  }

  // Validation: Diagonal elements of G must be 0
  for (int i = 0; i < m; ++i) {
    if (G(i, i) != 0.0) {
      throw std::invalid_argument("Diagonal elements of G must be equal to 0");
    }
  }

  // Validation: p must have m columns
  if (p.ncol != m) {
    throw std::invalid_argument("Invalid number of columns of p");
  }

  // Main algorithm
  for (int iter = 0; iter < iters; ++iter) {
    // Reset for this iteration
    wx = w;  // copy w
    g = G;   // copy G
    std::fill(r.begin(), r.end(), 0);
    pmax = 0.0;

    // Extract row iter from p
    for (int i = 0; i < m; ++i) {
      pvalues[i] = p(iter, i);
    }

    for (int step = 0; step < m; ++step) {

      // Compute ratios of raw p-values divided by weights
      std::fill(q.begin(), q.end(), 0.0);
      for (int i = 0; i < m; ++i) {
        if (wx[i] > 0.0) {
          q[i] = pvalues[i] / wx[i];
        }
      }

      // Replace q[i] == 0.0 with max(q) + 1.0
      double max_q = *std::max_element(q.begin(), q.end());
      for (int i = 0; i < m; ++i) {
        if (q[i] == 0.0) {
          q[i] = max_q + 1.0;
        }
      }

      // Identify the hypothesis to reject (minimum q)
      int j = static_cast<int>(std::distance(q.begin(),
                                             std::min_element(q.begin(), q.end())));
      padj(iter, j) = std::max(std::min(q[j], 1.0), pmax);
      pmax = padj(iter, j);
      r[j] = 1;

      // Update weights
      for (int l = 0; l < m; ++l) {
        if (r[l] == 0) {
          wx[l] = wx[l] + wx[j] * g(j, l);
        }
      }
      wx[j] = 0.0;

      // Update transition matrix
      std::fill(g1.data.begin(), g1.data.end(), 0.0);

      for (int l = 0; l < m; ++l) {
        if (r[l] == 0) {
          for (int k = 0; k < m; ++k) {
            if ((r[k] == 0) && (l != k) && (g(l, j) * g(j, l) < 1.0 - 1.0e-12)) {
              g1(l, k) = (g(l, k) + g(l, j) * g(j, k)) / (1.0 - g(l, j) * g(j, l));
            }
          }
        }
      }
      g = g1;  // copy g1 to g
    }
  }

  return padj;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix fadjpboncpp(const Rcpp::NumericVector& w,
                                const Rcpp::NumericMatrix& G,
                                const Rcpp::NumericMatrix& p) {
  auto w1 = Rcpp::as<std::vector<double>>(w);
  FlatMatrix G1 = flatmatrix_from_Rmatrix(G);
  FlatMatrix p1 = flatmatrix_from_Rmatrix(p);
  FlatMatrix padj1 = fadjpboncpp1(w1, G1, p1);
  return Rcpp::wrap(padj1);
}


FlatMatrix fwgtmatcpp(const std::vector<double>& w,
                      const FlatMatrix& G) {
  int m = static_cast<int>(w.size());
  int ntests = (1 << m) - 1;

  // Validation: w must be nonnegative
  if (std::any_of(w.begin(), w.end(), [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("w must be nonnegative");
  }

  // Validation: w must sum to 1
  double sum_w = std::accumulate(w.begin(), w.end(), 0.0);
  if (std::fabs(sum_w - 1.0) > 1.0e-12) {
    throw std::invalid_argument("w must sum to 1");
  }

  // Validation: G dimension
  if (G.nrow != m || G.ncol != m) {
    throw std::invalid_argument("Invalid dimension for G");
  }

  // Validation: G must be nonnegative
  for (int i = 0; i < m * m; ++i) {
    if (G.data[i] < 0.0) {
      throw std::invalid_argument("G must be nonnegative");
    }
  }

  // Validation: Row sums of G must be <= 1
  for (int i = 0; i < m; ++i) {
    double rowsum = 0.0;
    for (int j = 0; j < m; ++j) {
      rowsum += G(i, j);
    }
    if (rowsum > 1.0 + 1.0e-8) {
      throw std::invalid_argument("Row sums of G must be less than or equal to 1");
    }
  }

  // Validation: Diagonal elements of G must be 0
  for (int i = 0; i < m; ++i) {
    if (G(i, i) != 0.0) {
      throw std::invalid_argument("Diagonal elements of G must be equal to 0");
    }
  }

  std::vector<double> wx = w;  // copy w
  FlatMatrix g = G;            // copy G
  FlatMatrix wgtmat(ntests, m);
  FlatMatrix gtrmat((ntests + 1) / 2, m * m); // only need to store first half

  for (int i = 0; i < ntests; ++i) {
    int number = ntests - i;

    // Binary representation of elementary hypotheses in the intersection
    std::vector<int> cc(m);
    for (int j = 0; j < m; ++j) {
      cc[j] = (number / (1 << (m - 1 - j))) % 2;
    }

    if (i >= 1) {
      // Find index of minimum element in cc
      int j = static_cast<int>(std::distance(cc.begin(),
                                             std::min_element(cc.begin(), cc.end())));

      // Indicators for hypotheses not in the super set: cc1 = 1 - cc
      std::vector<int> cc1(m);
      for (int k = 0; k < m; ++k) {
        cc1[k] = 1 - cc[k];
      }
      cc1[j] = 0;

      // Index of the super set
      int ip = 0;
      for (int k = 0; k < m; ++k) {
        if (cc1[k]) {
          ip += (1 << (m - 1 - k));
        }
      }

      // Load the weights from the super set
      for (int k = 0; k < m; ++k) {
        wx[k] = wgtmat(ip, k);
      }

      // Load the transition matrix from the super set
      for (int l = 0; l < m; ++l) {
        for (int k = 0; k < m; ++k) {
          g(k, l) = gtrmat(ip, k * m + l);
        }
      }

      // Update the weights
      for (int k = 0; k < m; ++k) {
        if (cc[k]) {
          wx[k] += wx[j] * g(j, k);
        }
      }
      wx[j] = 0;

      // Update the transition matrix
      FlatMatrix g1(m, m);
      std::fill(g1.data.begin(), g1.data.end(), 0.0);
      for (int l = 0; l < m; ++l) {
        for (int k = 0; k < m; ++k) {
          if (cc[k] && cc[l] && (k != l) && (g(k, j) * g(j, k) < 1.0 - 1e-12)) {
            g1(k, l) = (g(k, l) + g(k, j) * g(j, l)) / (1 - g(k, j) * g(j, k));
          }
        }
      }
      g = g1;  // copy g1 to g
    }

    // Save the weights
    for (int k = 0; k < m; ++k) {
      wgtmat(i, k) = wx[k];
    }

    // Save the transition matrix
    if (i < (ntests + 1) / 2) {
      for (int l = 0; l < m; ++l) {
        for (int k = 0; k < m; ++k) {
          gtrmat(i, k * m + l) = g(k, l);
        }
      }
    }
  }

  return wgtmat;
}

//' @title Weight Matrix for All Intersection Hypotheses
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
Rcpp::NumericMatrix fwgtmat(const Rcpp::NumericVector& w,
                            const Rcpp::NumericMatrix& G) {
  auto w1 = Rcpp::as<std::vector<double>>(w);
  FlatMatrix G1 = flatmatrix_from_Rmatrix(G);
  FlatMatrix wgtmat1 = fwgtmatcpp(w1, G1);
  return Rcpp::wrap(wgtmat1);
}


FlatMatrix fadjpsimcpp1(const FlatMatrix& wgtmat,
                        const FlatMatrix& p,
                        const BoolMatrix& family) {

  int ntests = wgtmat.nrow;
  int m = wgtmat.ncol;
  int niters = p.nrow;
  int nfams = family.nrow;
  BoolMatrix incid(ntests, m);  // Changed to BoolMatrix
  FlatMatrix pinter(niters, ntests);

  // Validation: family must have m columns
  if (family.ncol != m) {
    throw std::invalid_argument(
        "family must have as many individual hypotheses as columns");
  }

  // Validation: Each hypothesis should belong to one and only one family
  for (int j = 0; j < m; ++j) {
    int colsum = 0;
    for (int k = 0; k < nfams; ++k) {
      if (family(k, j) != 0) colsum++;
    }
    if (colsum != 1) {
      throw std::invalid_argument(
          "Each hypothesis should belong to one or only one family");
    }
  }

  for (int i = 0; i < ntests; ++i) {
    int number = ntests - i;

    // Binary representation of elementary hypotheses in the intersection
    std::vector<unsigned char> cc(m);
    for (int j = 0; j < m; ++j) {
      cc[j] = (number / (1 << (m - 1 - j))) % 2;
    }

    // Identify the active families and active hypotheses
    BoolMatrix family0(nfams, m);  // Changed to BoolMatrix
    for (int j = 0; j < m; ++j) {
      for (int k = 0; k < nfams; ++k) {
        family0(k, j) = (family(k, j) != 0 && cc[j]) ? 1 : 0;
      }
    }

    // Count total active hypotheses
    int nhyps = 0;
    for (int j = 0; j < m; ++j) {
      for (int k = 0; k < nfams; ++k) {
        if (family0(k, j) != 0) nhyps++;
      }
    }

    std::vector<int> nhyps0(nfams, 0);
    std::vector<int> hyp, fam;
    hyp.reserve(nhyps);
    fam.reserve(nhyps);

    for (int j = 0; j < m; ++j) {
      for (int k = 0; k < nfams; ++k) {
        if (family0(k, j) != 0) {
          nhyps0[k]++;  // number of active hypotheses in family k
          fam.push_back(k);   // family of the l-th active hypothesis
          hyp.push_back(j);   // index of the l-th active hypothesis
        }
      }
    }

    // Create subset of active families
    std::vector<unsigned char> sub(nfams);
    for (int k = 0; k < nfams; ++k) {
      sub[k] = (nhyps0[k] > 0) ? 1 : 0;
    }

    int nfamil1 = 0;  // number of active families
    for (int k = 0; k < nfams; ++k) {
      if (sub[k]) nfamil1++;
    }

    // nhyps1: # of active hypotheses by family (only for active families)
    std::vector<int> nhyps1;
    nhyps1.reserve(nfamil1);
    for (int k = 0; k < nfams; ++k) {
      if (sub[k]) {
        nhyps1.push_back(nhyps0[k]);
      }
    }

    // Extract weights
    std::vector<double> w(nhyps);
    for (int j = 0; j < nhyps; ++j) {
      w[j] = wgtmat(i, hyp[j]);
    }

    for (int iter = 0; iter < niters; ++iter) {
      std::vector<double> pval(nhyps), cw(nhyps);
      for (int j = 0; j < nhyps; ++j) {
        pval[j] = p(iter, hyp[j]);
      }

      // Sort p-values within each family and obtain associated cum weights
      int s = 0;
      for (int k = 0; k < nfamil1; ++k) {
        int t = nhyps1[k];

        // Extract p-values and weights in the family
        std::vector<double> p1(t), w1(t);
        for (int j = 0; j < t; ++j) {
          p1[j] = pval[s + j];
          w1[j] = w[s + j];
        }

        // Obtain the index of sorted p-values within the family
        std::vector<int> index = seqcpp(0, t-1);
        std::sort(index.begin(), index.end(),
                  [p1](const int& a, const int& b) {
                    return p1[a] < p1[b];
                  });

        // Replace original with sorted values
        for (int j = 0; j < t; ++j) {
          pval[s + j] = p1[index[j]];

          // Obtain the cumulative weights within each family
          if (j == 0) {
            cw[s + j] = w1[index[j]];
          } else {
            cw[s + j] = cw[s + j - 1] + w1[index[j]];
          }
        }

        s += t;
      }

      double q = 1.0;
      for (int j = 0; j < nhyps; ++j) {
        if (cw[j] > 0.0) {
          q = std::min(q, pval[j] / cw[j]);
        }
      }

      pinter(iter, i) = q;
    }

    // Save cc to row i of incid
    for (int j = 0; j < m; ++j) {
      incid(i, j) = cc[j];
    }
  }

  // Obtain the adjusted p-values for individual hypotheses
  FlatMatrix padj(niters, m);
  for (int j = 0; j < m; ++j) {
      for (int i = 0; i < ntests; ++i) {
        for (int iter = 0; iter < niters; ++iter) {
          if (incid(i, j) != 0 && pinter(iter, i) > padj(iter, j)) {
          padj(iter, j) = pinter(iter, i);
        }
      }
    }
  }

  return padj;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix fadjpsimcpp(const Rcpp::NumericMatrix& wgtmat,
                                const Rcpp::NumericMatrix& p,
                                const Rcpp::LogicalMatrix& family) {
  FlatMatrix wgtmat1 = flatmatrix_from_Rmatrix(wgtmat);
  FlatMatrix p1 = flatmatrix_from_Rmatrix(p);
  BoolMatrix family1 = boolmatrix_from_Rmatrix(family);
  FlatMatrix padj1 = fadjpsimcpp1(wgtmat1, p1, family1);
  return Rcpp::wrap(padj1);
}


FlatMatrix repeatedPValuecpp1(
    const int kMax,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const double maxInformation,
    const FlatMatrix& p,
    const FlatMatrix& information,
    const FlatMatrix& spendingTime) {

  int B = p.nrow;
  int k = p.ncol;

  // Validation: kMax
  if (kMax <= 0) {
    throw std::invalid_argument("kMax must be a positive integer");
  }

  // Convert typeAlphaSpending to lowercase
  std::string asf = typeAlphaSpending;
  std::for_each(asf.begin(), asf.end(), [](char& c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  // Validation: typeAlphaSpending
  if (!(asf == "of" || asf == "p" || asf == "wt" || asf == "sfof" ||
      asf == "sfp" || asf == "sfkd" || asf == "sfhsd" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }

  // Validation: parameterAlphaSpending
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument("parameterAlphaSpending must be positive for sfKD");
  }

  // Validation: maxInformation
  if (maxInformation <= 0) {
    throw std::invalid_argument("maxInformation must be positive");
  }

  // Process information matrix
  FlatMatrix info(B, k);
  if (information.ncol != k) {
    throw std::invalid_argument("Invalid number of columns for information");
  } else if (information.nrow != 1 && information.nrow != B) {
    throw std::invalid_argument("Invalid number of rows for information");
  } else if (information.nrow == 1 && B > 1) {
    // Replicate single row to all B rows
    for (int j = 0; j < k; ++j) {
      std::fill_n(info.data.data() + j * B, B, information(0, j));
    }
  } else {
    info = information;
  }

  // Validation: information elements
  for (int iter = 0; iter < B; ++iter) {
    if (info(iter, 0) <= 0) {
      throw std::invalid_argument("information must be positive");
    }

    // Check if elements are increasing
    if (k > 1) {
      for (int j = 1; j < k; ++j) {
        if (info(iter, j) - info(iter, j - 1) <= 0) {
          throw std::invalid_argument("information must be increasing over time");
        }
      }
    }
  }

  // Process spendingTime matrix
  FlatMatrix spendTime(B, k);
  if (spendingTime.nrow == 1 && spendingTime.ncol == 1 && spendingTime(0, 0) == 0) {
    std::fill(spendTime.data.begin(), spendTime.data.end(), NaN);
  } else if (spendingTime.ncol != k) {
    throw std::invalid_argument("Invalid number of columns for spendingTime");
  } else if (spendingTime.nrow != 1 && spendingTime.nrow != B) {
    throw std::invalid_argument("Invalid number of rows for spendingTime");
  } else if (spendingTime.nrow == 1 && B > 1) {
    // Replicate single row to all B rows
    for (int j = 0; j < k; ++j) {
      std::fill_n(spendTime.data.data() + j * B, B, spendingTime(0, j));
    }
  } else {
    spendTime = spendingTime;
  }

  // Validation: spendingTime elements
  for (int iter = 0; iter < B; ++iter) {
    // Check if all elements in row are NA
    bool all_st_na = true;
    for (int j = 0; j < k; ++j) {
      if (!std::isnan(spendTime(iter, j))) {
        all_st_na = false; break;
      }
    }

    if (!all_st_na) {
      if (spendTime(iter, 0) <= 0) {
        throw std::invalid_argument("spendingTime must be positive");
      }

      // Check if elements are increasing
      if (k > 1) {
        for (int j = 1; j < k; ++j) {
          if (spendTime(iter, j) - spendTime(iter, j - 1) <= 0) {
            throw std::invalid_argument("spendingTime must be increasing over time");
          }
        }
      }

      if (spendTime(iter, k - 1) > 1) {
        throw std::invalid_argument("spendingTime must be less than or equal to 1");
      }
    }
  }

  // Initialize result matrix with NaN
  FlatMatrix repp(B, k);
  std::fill(repp.data.begin(), repp.data.end(), NaN);

  // Main computation loop
  for (int iter = 0; iter < B; ++iter) {
    int L;

    // Check if all spendingTime values are NA
    bool all_st_na = true;
    for (int j = 0; j < k; ++j) {
      if (!std::isnan(spendTime(iter, j))) {
        all_st_na = false; break;
      }
    }

    if (all_st_na) {  // use information rates
      // Find if any info >= maxInformation
      bool any_gte_max = false;
      int first_gte_idx = -1;
      for (int j = 0; j < k; ++j) {
        if (info(iter, j) >= maxInformation) {
          any_gte_max = true;
          if (first_gte_idx == -1) first_gte_idx = j;
          break;
        }
      }

      if (!any_gte_max) {  // all observed info < maxinfo
        L = k - 1;
      } else {  // find index of first look with observed info >= maxinfo
        L = first_gte_idx;
      }
    } else {  // use spending time
      L = k - 1;
    }

    // Information time for forming covariance matrix of test statistics
    std::vector<double> t1(L + 1);
    for (int l = 0; l <= L; ++l) {
      t1[l] = info(iter, l) / info(iter, L);
    }

    // Spending time for error spending
    std::vector<double> s1(L + 1);
    if (all_st_na) {  // use information rates
      for (int l = 0; l <= L; ++l) {
        if (l == kMax - 1 || info(iter, l) >= maxInformation) {
          s1[l] = 1.0;
        } else {
          s1[l] = info(iter, l) / maxInformation;
        }
      }
    } else {  // using spending time
      for (int l = 0; l <= L; ++l) {
        s1[l] = spendTime(iter, l);
      }
    }

    // Compute repeated p-values
    std::vector<double> t, s;
    std::vector<unsigned char> x;
    t.reserve(L + 1);  // Pre-allocate maximum size
    s.reserve(L + 1);
    x.reserve(L + 1);
    for (int i = 0; i <= L; ++i) {
      t.assign(t1.begin(), t1.begin() + i + 1);
      s.assign(s1.begin(), s1.begin() + i + 1);
      x.assign(i + 1, 1);
      double pvalue = p(iter, i);

      // Lambda function for root finding
      auto f = [&](double a)->double {
        auto u = getBoundcpp(i+1, t, a, asf, parameterAlphaSpending, {NaN}, s, x);
        return 1.0 - boost_pnorm(u[i]) - pvalue;
      };

      // Find root
      if (f(0.000001) > 0) {
        repp(iter, i) = 0.000001;
      } else if (f(0.999999) < 0) {
        repp(iter, i) = 0.999999;
      } else {
        repp(iter, i) = brent(f, 0.000001, 0.999999, 1.0e-6);
      }
    }
  }

  return repp;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix repeatedPValuecpp(
    const int kMax,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const double maxInformation,
    const Rcpp::NumericMatrix& p ,
    const Rcpp::NumericMatrix& information,
    const Rcpp::NumericMatrix& spendingTime) {
  FlatMatrix p1 = flatmatrix_from_Rmatrix(p);
  FlatMatrix information1 = flatmatrix_from_Rmatrix(information);
  FlatMatrix spendingTime1 = flatmatrix_from_Rmatrix(spendingTime);
  FlatMatrix repp1 = repeatedPValuecpp1(
    kMax, typeAlphaSpending, parameterAlphaSpending,
    maxInformation, p1, information1, spendingTime1);
  return Rcpp::wrap(repp1);
}
