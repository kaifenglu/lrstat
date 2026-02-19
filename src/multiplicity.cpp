// [[Rcpp::depends(RcppParallel)]]
#include "generic_design.h"
#include "utilities.h"
#include "dataframe_list.h"

#include <algorithm>     // any_of, distance, fill, max_element, min, sort
#include <cctype>        // tolower
#include <cmath>         // fabs, isnan
#include <cstring>       // memcpy
#include <numeric>       // accumulate
#include <stdexcept>     // invalid_argument
#include <string>        // string
#include <vector>        // vector

#include <Rcpp.h>
#include <RcppParallel.h>


// Helper to update graph for graphical approaches
ListCpp updateGraphcpp(const std::vector<double>& w,
                       const FlatMatrix& G,
                       const std::vector<int>& I,
                       const int j) {

  const int m = static_cast<int>(w.size());

  // Validation: w must be nonnegative
  if (std::any_of(w.begin(), w.end(), [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("w must be nonnegative");
  }

  // Validation: w must sum to 1
  double sum_w = std::accumulate(w.begin(), w.end(), 0.0);
  if (std::fabs(sum_w - 1.0) > 1e-12) {
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
  std::vector<double> rowsum(m, 0.0);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      rowsum[i] += G(i, j);
    }
  }
  for (int i = 0; i < m; ++i) {
    if (rowsum[i] > 1.0 + 1e-8) {
      throw std::invalid_argument("Row sums of G must be less than or equal to 1");
    }
  }

  // Validation: Diagonal elements of G must be 0
  for (int i = 0; i < m; ++i) {
    if (G(i, i) != 0.0) {
      throw std::invalid_argument("Diagonal elements of G must be equal to 0");
    }
  }

  // Check validity of elements of I, create zero-based indexing, and remove j
  std::vector<unsigned char> seen(m, 0); // 0/1 marks presence
  std::vector<int> I_zero; I_zero.reserve(I.size()); // zero-based indexing of I
  std::vector<int> I_new; I_new.reserve(I.size()); // new I after removing j
  std::vector<int> I_zero_new; I_zero_new.reserve(I.size());
  for (int val : I) {
    if (val < 1 || val > m) {
      throw std::invalid_argument("Elements of I must be integers between 1 and m.");
    }
    int idx = val - 1;
    if (seen[idx]) {
      throw std::invalid_argument("The index set I must not contain duplicates.");
    }
    seen[idx] = 1;
    I_zero.push_back(idx);
    if (val != j) {
      I_new.push_back(val); // keep 1-based for return
      I_zero_new.push_back(idx);
    }
  }

  // Check that j belongs to I
  int j1 = j - 1;
  if (!seen[j1]) throw std::invalid_argument("j must be in I.");

  // Update weights
  std::vector<double> wx = w;  // copy w
  const double wxj = wx[j1];
  for (int idx : I_zero_new) {
    wx[idx] += wxj * G(j1, idx);
  }
  wx[j1] = 0.0;

  // Update transition matrix
  std::vector<double> denom(m);
  for (int l : I_zero_new) {
    denom[l] =  1.0 - G(l, j1) * G(j1, l);
  }

  FlatMatrix g(m, m);
  for (int k : I_zero_new) {
    double g_jk = G(j1, k);
    for (int l : I_zero_new) {
      if (l == k) continue;
      double dl = denom[l];
      if (dl > 1e-12) {
        g(l, k) = (G(l, k) + G(l, j1) * g_jk) / dl;
      }
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
  auto G1 = flatmatrix_from_Rmatrix(G);
  auto cpp_result = updateGraphcpp(w1, G1, I1, j);
  return Rcpp::wrap(cpp_result);
}


// Helper to compute adjusted p-values for Bonferroni-based graphical approaches
FlatMatrix fadjpboncpp1(const std::vector<double>& w,
                        const FlatMatrix& G,
                        const FlatMatrix& p) {

  const int m = static_cast<int>(w.size());
  const int iters = p.nrow;

  // Validation: w must be nonnegative
  if (std::any_of(w.begin(), w.end(), [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("w must be nonnegative");
  }

  // Validation: w must sum to 1
  double sum_w = std::accumulate(w.begin(), w.end(), 0.0);
  if (std::fabs(sum_w - 1.0) > 1e-12) {
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
  std::vector<double> rowsum(m, 0.0);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      rowsum[i] += G(i, j);
    }
  }
  for (int i = 0; i < m; ++i) {
    if (rowsum[i] > 1.0 + 1e-8) {
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


  // Output and reusable buffers
  FlatMatrix padj(iters, m);  // adjusted p-values
  std::vector<unsigned char> r(m); // rejection indicators
  std::vector<double> pvalues(m);  // raw p-values
  std::vector<double> wx(m);    // dynamic weights
  std::vector<double> denom(m); // denom reused
  std::vector<int> active;      // currently active hypotheses
  active.reserve(m);
  std::vector<int> pos(m, -1); // position map: pos[idx] = index in active or -1
  FlatMatrix g(m, m), g1(m, m);  // dynamic transition matrices

  // Main outer loop: per-iteration (each row of p)
  for (int iter = 0; iter < iters; ++iter) {
    // Reset for this iteration
    wx = w;  // copy w
    g = G;   // copy G
    std::fill_n(r.data(), r.size(), 0);

    // Extract row iter from p into pvalues
    for (int i = 0; i < m; ++i) {
      pvalues[i] = p(iter, i);
    }

    // build initial active list and pos map
    active.clear();
    for (int idx = 0; idx < m; ++idx) {
      active.push_back(idx);
      pos[idx] = idx;
    }

    double pmax = 0.0;

    // Reuse g1 contents: ensure it's zeroed before first use
    g1.fill(0.0);

    // Steps: reject one hypothesis at a time
    for (int step = 0; step < m; ++step) {

      // Find min ratio among active indices (ignoring zero weights)
      double min_q = POS_INF;
      int j = -1;
      for (int i : active) {
        double wxi = wx[i];
        if (wxi > 0.0) {
          double qval = pvalues[i] / wxi;
          if (qval < min_q) {
            min_q = qval;
            j = i;
          }
        }
      }
      if (j == -1 || !std::isfinite(min_q)) {
        // no eligible hypotheses left
        break;
      }

      // Accept and record adjusted p-value (bounded by [pmax, 1.0])
      double qbounded = std::min(min_q, 1.0);
      if (qbounded < pmax) qbounded = pmax;
      pmax = qbounded;
      padj(iter, j) = pmax;

      // Mark rejection and remove j from active (swap-remove)
      r[j] = 1;
      int remove_pos = pos[j];
      int last_pos = static_cast<int>(active.size()) - 1;
      if (remove_pos != last_pos) {
        int swapped_idx = active[last_pos];
        active[remove_pos] = swapped_idx;
        pos[swapped_idx] = remove_pos;
      }
      active.pop_back();
      pos[j] = -1;

      // Update weights wx for remaining active indices
      double wxj = wx[j];
      for (int l : active) {
        wx[l] += wxj * g(j, l);
      }
      wx[j] = 0.0;

      // Precompute denom for surviving l entries only
      for (int l : active) {
        denom[l] =  1.0 - g(l, j) * g(j, l);
      }

      // Build new g1 columns only for remaining k in active and l in active
      for (int k : active) {
        double g_jk = g(j, k);
        for (int l : active) {
          if (l == k) continue;
          double dl = denom[l];
          if (dl > 1e-12) {
            g1(l, k) = (g(l, k) + g(l, j) * g_jk) / dl;
          }
        }
      }

      std::swap(g.data, g1.data);  // copy g1 to g for next iteration
    }
  }

  return padj;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix fadjpboncpp(const Rcpp::NumericVector& w,
                                const Rcpp::NumericMatrix& G,
                                const Rcpp::NumericMatrix& p) {
  auto w1 = Rcpp::as<std::vector<double>>(w);
  auto G1 = flatmatrix_from_Rmatrix(G);
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto padj1 = fadjpboncpp1(w1, G1, p1);
  return Rcpp::wrap(padj1);
}


// Helper to compute the full weight matrix for graphical approaches
FlatMatrix fwgtmatcpp(const std::vector<double>& w,
                      const FlatMatrix& G) {
  int m = static_cast<int>(w.size());
  int ntests = (1 << m) - 1;
  const int gtr_nrow = (ntests + 1) / 2;
  const int gtr_ncol = m * m;

  // Validation: w must be nonnegative
  if (std::any_of(w.begin(), w.end(), [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("w must be nonnegative");
  }

  // Validation: w must sum to 1
  double sum_w = std::accumulate(w.begin(), w.end(), 0.0);
  if (std::fabs(sum_w - 1.0) > 1e-12) {
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
  std::vector<double> rowsum(m, 0.0);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      rowsum[i] += G(i, j);
    }
  }
  for (int i = 0; i < m; ++i) {
    if (rowsum[i] > 1.0 + 1e-8) {
      throw std::invalid_argument("Row sums of G must be less than or equal to 1");
    }
  }

  // Validation: Diagonal elements of G must be 0
  for (int i = 0; i < m; ++i) {
    if (G(i, i) != 0.0) {
      throw std::invalid_argument("Diagonal elements of G must be equal to 0");
    }
  }

  // Preallocations and initial copies
  std::vector<double> wx = w;   // dynamic weights, reused
  FlatMatrix g = G;             // mutable transition matrix
  FlatMatrix g1(m, m); // temp transition matrix, reused (zeroed when needed)
  FlatMatrix gtrmat(gtr_nrow, gtr_ncol); // store first half of transition matrix
  FlatMatrix wgtmat(ntests, m); // output

  std::vector<int> active; // indices of active hypotheses in intersection
  active.reserve(m);

  // bitmask for m bits, used to flip bits and find super set
  const int mask = (1 << m) - 1;
  std::vector<double> denom(m); // temp denominator for transition matrix update

  for (int i = 0; i < ntests; ++i) {
    if (i >= 1) {
      int number = ntests - i;  // original mapping

      // Build list of active indices
      active.clear();
      int j = 0; // index of minimum active hypothesis
      bool found_zero = false;
      for (int k = 0; k < m; ++k) {
        int bit = (number >> (m - 1 - k)) & 1;
        if (bit) active.push_back(k);
        else if (!found_zero) { j = k; found_zero = true; }
      }
      if (!found_zero) j = 0;

      // index of the super set, with j-th bit set to 0 and others flipped
      int ip = ((~number) & mask) & ~(1 << (m - 1 - j));

      // Load the weights from the super set
      for (int k = 0; k < m; ++k) {
        wx[k] = wgtmat(ip, k);
      }

      // Load the transition matrix from the super set
      for (int k = 0; k < m; ++k) {
        for (int l = 0; l < m; ++l) {
          g(l, k) = gtrmat(ip, k * m + l);
        }
      }

      // Update the weights
      double wxj = wx[j];
      if (wxj != 0.0) {
        for (int k : active) {
          wx[k] += wxj * g(j, k);
        }
        wx[j] = 0.0;
      }

      // Update the transition matrix
      for (int l : active) {
        denom[l] =  1.0 - g(l, j) * g(j, l);
      }

      g1.fill(0.0); // ensure g1 is zeroed before use
      for (int k : active) {
        double g_jk = g(j, k);
        for (int l : active) {
          if (l == k) continue;
          double dl = denom[l];
          if (dl > 1e-12) {
            g1(l, k) = (g(l, k) + g(l, j) * g_jk) / dl;
          }
        }
      }

      std::swap(g.data, g1.data);  // copy g1 to g for next iteration
    }

    // Save the weights
    for (int k = 0; k < m; ++k) {
      wgtmat(i, k) = wx[k];
    }

    // Save the transition matrix
    if (i < gtr_nrow) {
      for (int k = 0; k < m; ++k) {
        for (int l = 0; l < m; ++l) {
          gtrmat(i, k * m + l) = g(l, k);
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
  auto G1 = flatmatrix_from_Rmatrix(G);
  auto wgtmat1 = fwgtmatcpp(w1, G1);
  return Rcpp::wrap(wgtmat1);
}


// Helper to compute adjusted p-values for Simes-based graphical approaches
FlatMatrix fadjpsimcpp1(const FlatMatrix& wgtmat,
                        const FlatMatrix& p,
                        const BoolMatrix& family) {

  const int ntests = wgtmat.nrow;
  const int m = wgtmat.ncol;
  const int niters = p.nrow;
  const int nfams = family.nrow;

  if (family.ncol != m) {
    throw std::invalid_argument(
        "family must have as many individual hypotheses as columns");
  }

  // Precompute hyps in each family (static)
  std::vector<std::vector<int>> hyps_in_family(nfams);
  for (int j = 0; j < m; ++j) {
    int found = 0;
    for (int k = 0; k < nfams; ++k) {
      if (family(k, j)) {
        hyps_in_family[k].push_back(j);
        ++found;
      }
    }
    if (found != 1) {
      throw std::invalid_argument(
          "Each hypothesis should belong to one or only one family");
    }
  }

  // Precompute, for each subset i, the active hypotheses (hyp) and
  // the family block sizes nhyps1
  std::vector<std::vector<int>> subset_hyp(ntests);
  std::vector<std::vector<int>> subset_nhyps1(ntests);

  for (int i = 0; i < ntests; ++i) {
    int number = ntests - i; // MSB-first bit convention as original

    std::vector<int> nhyps0(nfams, 0);
    std::vector<int>& hyp = subset_hyp[i];
    hyp.clear();
    for (int j = 0; j < m; ++j) {
      int bit = (number >> (m - 1 - j)) & 1;
      if (!bit) continue;
      for (int k = 0; k < nfams; ++k) {
        if (family(k, j)) {
          nhyps0[k]++;
          hyp.push_back(j);
          break; // a hypothesis belongs to exactly one family
        }
      }
    }
    std::vector<int>& nh1 = subset_nhyps1[i];
    nh1.clear();
    for (int k = 0; k < nfams; ++k) {
      if (nhyps0[k] > 0) nh1.push_back(nhyps0[k]);
    }
  }

  // Output matrix initialized to 0
  FlatMatrix padj(niters, m);

  // Reusable buffers sized up to m
  std::vector<double> wbuf; wbuf.reserve(m);
  std::vector<double> pbuf; pbuf.reserve(m);
  std::vector<double> cw;   cw.reserve(m);
  std::vector<int> idx;     idx.reserve(m);
  std::vector<double> pinter_col(niters);

  // Main loop over subsets
  for (int i = 0; i < ntests; ++i) {
    const std::vector<int>& hyp = subset_hyp[i];
    const std::vector<int>& nhyps1 = subset_nhyps1[i];
    const int nhyps = static_cast<int>(hyp.size());

    // Extract weights wx for this subset from wgtmat row i (column-major)
    wbuf.assign(nhyps, 0.0);
    for (int t = 0; t < nhyps; ++t) {
      int col = hyp[t];
      wbuf[t] = wgtmat(i, col);
    }

    // For each iteration (row of p), compute pinter(iter, i)
    for (int iter = 0; iter < niters; ++iter) {
      // Extract p-values for active hypotheses in the same hyp order
      pbuf.resize(nhyps);
      cw.assign(nhyps, 0.0);
      for (int t = 0; t < nhyps; ++t) {
        int col = hyp[t];
        pbuf[t] = p(iter, col);
      }

      // Sort p-values within each active family block
      // (nhyps1 gives block sizes in family order).
      int s = 0;
      for (int block = 0; block < (int)nhyps1.size(); ++block) {
        int t = nhyps1[block];
        // snapshot original block to avoid in-place overwrite corruption
        // p1 and w1 are small (<= m) and allocated on the heap but
        // re-used across iterations
        std::vector<double> p1(t), w1(t); // within family block
        for (int u = 0; u < t; ++u) {
          p1[u] = pbuf[s + u];
          w1[u] = wbuf[s + u];
        }

        // obtain index ordering of the snapshot
        idx = seqcpp(0, t-1);
        std::sort(idx.begin(), idx.end(), [&p1](int a, int b) {
                    return p1[a] < p1[b];
                  });

        // copy sorted p and compute cumulative weights from the snapshot
        double cum = 0.0;
        for (int j = 0; j < t; ++j) {
          int src = idx[j];
          double pv = p1[src];
          double wv = w1[src];
          cum += wv;
          pbuf[s + j] = pv;    // write sorted p-values into pbuf
          cw[s + j] = cum;     // cumulative weights
        }
        s += t;
      }

      // compute q = min_j pbuf[j] / cw[j] ignoring cw==0
      double q = 1.0;
      for (int j = 0; j < nhyps; ++j) {
        double cj = cw[j];
        if (cj > 0.0) {
          double ratio = pbuf[j] / cj;
          if (ratio < q) q = ratio;
        }
      }
      pinter_col[iter] = q;
    } // end iter loop

    // Update padj columns for active hypotheses (each hyp[t] gets max over subsets)
    for (int t = 0; t < nhyps; ++t) {
      int col = hyp[t];
      double* padj_col = &padj.data[col * niters]; // contiguous column
      for (int iter = 0; iter < niters; ++iter) {
        double v = pinter_col[iter];
        if (v > padj_col[iter]) padj_col[iter] = v;
      }
    }
  } // end subsets

  return padj;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix fadjpsimcpp(const Rcpp::NumericMatrix& wgtmat,
                                const Rcpp::NumericMatrix& p,
                                const Rcpp::LogicalMatrix& family) {
  auto wgtmat1 = flatmatrix_from_Rmatrix(wgtmat);
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto family1 = boolmatrix_from_Rmatrix(family);
  auto padj1 = fadjpsimcpp1(wgtmat1, p1, family1);
  return Rcpp::wrap(padj1);
}


// Helper to compute repeated p-values for alpha spending approaches
FlatMatrix repeatedPValuecpp1(
    const int kMax,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const double maxInformation,
    const FlatMatrix& p,
    const FlatMatrix& information,
    const FlatMatrix& spendingTime) {

  // Validation: kMax must be a positive integer
  if (kMax <= 0) {
    throw std::invalid_argument("kMax must be a positive integer");
  }

  // Convert typeAlphaSpending to lowercase
  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

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

  // Validation: dimensions of p and information
  int B = p.nrow;
  int k1 = p.ncol;

  if (k1 > kMax) {
    throw std::invalid_argument("Number of columns in p must not exceed kMax");
  }

  // Process information matrix with possible broadcasting
  FlatMatrix info(B, k1);
  if (information.ncol != k1) {
    throw std::invalid_argument("Invalid number of columns for information");
  }
  if (information.nrow != 1 && information.nrow != B) {
    throw std::invalid_argument("Invalid number of rows for information");
  }

  const double* in_ptr = information.data_ptr();
  double* info_ptr = info.data_ptr();
  if (information.nrow == B) {
    // shapes match -> fast block copy
    std::memcpy(info_ptr, in_ptr, B * k1 * sizeof(double));
  } else {
    // information.nrow == 1 -> broadcast the single row to all B rows
    // For each column k, value = in_ptr[k], destination column is info_ptr + k*B
    for (int k = 0; k < k1; ++k) {
      double v = in_ptr[k];                    // source scalar for column k
      double* dst_col = info_ptr + k * B;      // column-major base
      std::fill_n(dst_col, B, v);              // efficient bulk write
    }
  }

  // Process spendingTime matrix with possible broadcasting and NA handling
  FlatMatrix spendTime(B, k1);
  if (spendingTime.nrow == 1 && spendingTime.ncol == 1 && spendingTime(0, 0) == 0) {
    spendTime.fill(NaN);
  } else {
    if (spendingTime.ncol != k1) {
      throw std::invalid_argument("Invalid number of columns for spendingTime");
    }
    if (spendingTime.nrow != 1 && spendingTime.nrow != B) {
      throw std::invalid_argument("Invalid number of rows for spendingTime");
    }

    const double* sp_ptr = spendingTime.data_ptr();
    double* spend_ptr = spendTime.data_ptr();
    if (spendingTime.nrow == B) {
      // shapes match -> fast block copy
      std::memcpy(spend_ptr, sp_ptr, B * k1 * sizeof(double));
    } else {
      // spendingTime.nrow == 1 -> broadcast the single row to all B rows
      // For each column k, value = sp_ptr[k], destination column is spend_ptr + k*B
      for (int k = 0; k < k1; ++k) {
        double v = sp_ptr[k];                    // source scalar for column k
        double* dst_col = spend_ptr + k * B;     // column-major base
        std::fill_n(dst_col, B, v);              // efficient bulk write
      }
    }
  }


  auto f = [&](const int b)-> std::vector<double> {
    std::vector<double> p_vec(k1), i_vec(k1), s_vec(k1);
    for (int k = 0; k < k1; ++k) {
      p_vec[k] = p(b, k);
      i_vec[k] = info(b, k);
      s_vec[k] = spendTime(b, k);
    }

    // Validation: p-values must be between 0 and 1 and not NA
    if (std::any_of(p_vec.begin(), p_vec.end(),
                    [](double v){ return std::isnan(v); })) {
      throw std::invalid_argument("p-values must be provided");
    }
    if (std::any_of(p_vec.begin(), p_vec.end(),
                    [](double v){ return v < 0 || v > 1; })) {
      throw std::invalid_argument("p-values must be between 0 and 1");
    }

    // Validation: information must be positive, increasing, and not NA
    if (std::any_of(i_vec.begin(), i_vec.end(),
                    [](double v){ return std::isnan(v); })) {
      throw std::invalid_argument("information must be provided");
    }
    if (i_vec[0] <= 0.0) {
      throw std::invalid_argument("information must be positive");
    }
    if (any_nonincreasing(i_vec)) {
      throw std::invalid_argument("information must be increasing over time");
    }

    // Validation: spendingTime must be positive, increasing, and not exceeding 1
    bool all_na = std::all_of(s_vec.data(), s_vec.data() + k1,
                              [](double v){ return std::isnan(v); });

    if (!all_na) {
      if (std::any_of(s_vec.data(), s_vec.data() + k1,
                      [](double v){ return std::isnan(v); })) {
        throw std::invalid_argument("spendingTime must be provided");
      }
      if (s_vec[0] <= 0.0) {
        throw std::invalid_argument("spendingTime must be positive");
      }
      if (any_nonincreasing(s_vec)) {
        throw std::invalid_argument("spendingTime must be increasing over time");
      }
       if (s_vec[k1 - 1] > 1.0) {
        throw std::invalid_argument("spendingTime must not exceed 1");
      }
    }


    // Determine L based on maxInformation and spendingTime
    int L;
    if (all_na) {  // use information rates
      // Find if any information >= maxInformation
      auto it = std::lower_bound(i_vec.begin(), i_vec.end(), maxInformation);
      if (it == i_vec.end()) { // none >= maxInformation
        L = k1;
      } else {
        L = static_cast<int>(std::distance(i_vec.begin(), it)) + 1;
      }
    } else {  // use spending time
      L = k1;
    }

    std::vector<double> p1(k1), t1(k1), s1(k1);
    std::memcpy(p1.data(), p_vec.data(), L * sizeof(double));

    // Information time for forming covariance matrix of test statistics
    double info_L = i_vec[L-1];
    for (int l = 0; l < L; ++l) {
      t1[l] = i_vec[l] / info_L;
    }

    // Spending time for error spending
    if (all_na) {  // use information rates
      for (int l = 0; l < L; ++l) {
        if (l == kMax - 1 || i_vec[l] >= maxInformation) {
          s1[l] = 1.0; // the last look is at or beyond maxInformation
        } else {
          s1[l] = i_vec[l] / maxInformation;
        }
      }
    } else {  // using spending time
      std::memcpy(s1.data(), s_vec.data(), L * sizeof(double));
    }

    // Compute repeated p-values
    std::vector<double> empty_user;
    std::vector<double> repp(k1);
    for (int i = 0; i < L; ++i) {
      double pvalue = p1[i];
      std::vector<double> t0(t1.begin(), t1.begin() + i + 1);
      std::vector<double> s0(s1.begin(), s1.begin() + i + 1);
      std::vector<unsigned char> x(i + 1, 1);

      BoundCacheAlpha cache(i + 1, t0, asf, parameterAlphaSpending,
                            empty_user, s0, x, 64, 12);

      // Lambda function for root finding to solve for the repeated p-value at step i
      auto f = [&](double a)->double {
        auto u = cache.get(a);
        return 1.0 - boost_pnorm(u[i]) - pvalue;
      };

      // Find root in (0, 1) using brent's method, with checks at the endpoints
      double fl = f(0.000001);
      if (fl > 0) {
        repp[i] = 0.000001;
      } else {
        double fh = f(0.999999);
        if (fh < 0) {
          repp[i] = 0.999999;
        } else {
          repp[i] = brent(f, 0.000001, 0.999999, 1e-6);
        }
      }
    }

    return repp;
  };


  struct SimulationWorker : public RcppParallel::Worker {
    // references to read-only inputs (no mutation)
    const int k1;
    const std::string& asf;
    const double parameterAlphaSpending;
    const int kMax;
    const double maxInformation;
    const FlatMatrix& p;
    const FlatMatrix& info;
    const FlatMatrix& spendTime;

    // function f and other params that f needs are captured from outer scope
    // capture them by reference here so worker can call f(...)
    std::function<std::vector<double>(const int)> f;

    // result references (each iteration writes unique index into these)
    FlatMatrix& repp_out; // B by k1 matrix to store repeated p-values

    // constructor
    SimulationWorker(const int k1_,
                     const std::string& asf_,
                     const double parameterAlphaSpending_,
                     const int kMax_,
                     const double maxInformation_,
                     const FlatMatrix& p_,
                     const FlatMatrix& info_,
                     const FlatMatrix& spendTime_,
                     decltype(f) f_,
                     FlatMatrix& repp_out_) :

      k1(k1_), asf(asf_), parameterAlphaSpending(parameterAlphaSpending_),
      kMax(kMax_), maxInformation(maxInformation_),
      p(p_), info(info_), spendTime(spendTime_),
      f(std::move(f_)),
      repp_out(repp_out_) {}

    // operator() processes a range of bootstrap iterations [begin, end)
    void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t b = begin; b < end; ++b) {
        // call the (thread-safe) per-iteration function f
        std::vector<double> out = f(static_cast<int>(b));

        // write results
        for (int k = 0; k < k1; ++k) {
          repp_out(b, k) = out[k];
        }
      } // end for b
    } // end operator()
  }; // end BootstrapWorker

  FlatMatrix repp_mat(B, k1); // B by k1 matrix to store repeated values

  // Instantiate the Worker with references to inputs and outputs
  SimulationWorker worker(
      k1, asf, parameterAlphaSpending, kMax,
      maxInformation, p, info, spendTime,
      // bind f into std::function (capture the f we already have)
      std::function<std::vector<double>(const int)>(f),
      repp_mat
  );

  // Run the parallel loop over iterations
  RcppParallel::parallelFor(0, B, worker, 1 /*grain size*/);

  return worker.repp_out;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix repeatedPValuecpp(
    const int kMax,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const double maxInformation,
    const Rcpp::NumericMatrix& p,
    const Rcpp::NumericMatrix& information,
    const Rcpp::NumericMatrix& spendingTime) {

  auto p1 = flatmatrix_from_Rmatrix(p);
  auto information1 = flatmatrix_from_Rmatrix(information);
  auto spendingTime1 = flatmatrix_from_Rmatrix(spendingTime);

  auto repp = repeatedPValuecpp1(
    kMax, typeAlphaSpending, parameterAlphaSpending,
    maxInformation, p1, information1, spendingTime1);
  return Rcpp::wrap(repp);
}


// Helper to compute the first rejection step for Bonferroni-based graphical
// approaches in group sequential trials
IntMatrix fseqboncpp1(
    const std::vector<double>& w,
    const FlatMatrix& G,
    const double alpha,
    const int kMax,
    const std::vector<std::string>& typeAlphaSpending,
    const std::vector<double>& parameterAlphaSpending,
    const std::vector<double>& maxInformation,
    const BoolMatrix& incidenceMatrix,
    const int k1,
    const FlatMatrix& p,
    const FlatMatrix& information,
    const FlatMatrix& spendingTime) {

  int m = static_cast<int>(w.size());

  // Validation: w must be nonnegative
  if (std::any_of(w.begin(), w.end(), [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("w must be nonnegative");
  }

  // Validation: w must sum to 1
  double sum_w = std::accumulate(w.begin(), w.end(), 0.0);
  if (std::fabs(sum_w - 1.0) > 1e-12) {
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
  std::vector<double> rowsum(m, 0.0);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      rowsum[i] += G(i, j);
    }
  }
  for (int i = 0; i < m; ++i) {
    if (rowsum[i] > 1.0 + 1e-8) {
      throw std::invalid_argument("Row sums of G must be less than or equal to 1");
    }
  }

  // Validation: Diagonal elements of G must be 0
  for (int i = 0; i < m; ++i) {
    if (G(i, i) != 0.0) {
      throw std::invalid_argument("Diagonal elements of G must be equal to 0");
    }
  }

  // Validation: alpha must be in (0, 1)
  if (alpha <= 0.0 || alpha >= 1.0) {
    throw std::invalid_argument("alpha must be in (0, 1)");
  }

  // Validation: kMax must be a positive integer
  if (kMax <= 0) {
    throw std::invalid_argument("kMax must be a positive integer");
  }

  // Validation: alpha spending type and parameter values
  std::vector<std::string> asf = typeAlphaSpending;
  std::vector<double> asfpar = parameterAlphaSpending;

  if (asf.size() == 1) asf.resize(m, asf[0]);
  if (asf.size() != static_cast<std::size_t>(m)) {
    throw std::invalid_argument("Invalid length for typeAlphaSpending");
  }

  if (asfpar.size() == 1) asfpar.resize(m, asfpar[0]);
  if (asfpar.size() != static_cast<std::size_t>(m)) {
    throw std::invalid_argument("Invalid length for parameterAlphaSpending");
  }

  for (int i = 0; i < m; ++i) {
    std::string& asfi = asf[i];
    for (char &c : asfi) {
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (!(asfi == "of" || asfi == "p" || asfi == "wt" ||
        asfi == "sfof" || asfi == "sfp" || asfi == "sfkd" ||
        asfi == "sfhsd" || asfi == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpending");
    }

    if ((asfi == "wt" || asfi == "sfkd" || asfi == "sfhsd") && std::isnan(asfpar[i])) {
      throw std::invalid_argument("Missing value for parameterAlphaSpending");
    }

    if (asfi == "sfkd" && asfpar[i] <= 0) {
      throw std::invalid_argument("parameterAlphaSpending must be positive for sfKD");
    }
  }

  // Validation: maxInformation must be positive for each hypothesis
  if (maxInformation.size() != static_cast<std::size_t>(m)) {
    throw std::invalid_argument("Invalid length for maxInformation");
  }
  if (std::any_of(maxInformation.begin(), maxInformation.end(),
                  [](double val) { return val <= 0.0; })) {
    throw std::invalid_argument("maxInformation must be positive");
  }

  // Validation: incidenceMatrix dimension
  if (incidenceMatrix.ncol != m) {
    throw std::invalid_argument("Invalid number of columns for incidenceMatrix");
  }
  if (incidenceMatrix.nrow != kMax) {
    throw std::invalid_argument("Invalid number of rows for incidenceMatrix");
  }

  // Validation: k1 number of looks at interim analysis
  if (k1 <= 0) {
    throw std::invalid_argument("k1 must be a positive integer");
  }
  if (k1 > kMax) {
    throw std::invalid_argument("k1 must be less than or equal to kMax");
  }


  // Validation: p matrix dimension
  if (p.ncol != m) {
    throw std::invalid_argument("Invalid number of columns for p");
  }
  if (p.nrow % k1 != 0) {
    throw std::invalid_argument("p number of rows must be a multiple of k1");
  }

  int B = p.nrow / k1;

  // Validation: information matrix dimension
  if (information.ncol != m) {
    throw std::invalid_argument("Invalid number of columns for information");
  }
  if (information.nrow != k1 && information.nrow != B * k1) {
    throw std::invalid_argument("Invalid number of rows for information");
  }

  FlatMatrix info(B * k1, m);
  if (information.nrow == k1 && B > 1) {
    // replicate information for B iterations
    const double* src_ptr = information.data_ptr();
    double* dst_ptr = info.data_ptr();
    for (int j = 0; j < m; ++ j) {
      for (int b = 0; b < B; ++b) {
        std::memcpy(dst_ptr + (b * k1) + j * (B * k1),
                    src_ptr + j * k1,
                    k1 * sizeof(double));
      }
    }
  } else {
    info = information;
  }


  // Handle spending time
  FlatMatrix spendTime(B * k1, m);
  if (spendingTime.nrow == 1 && spendingTime.ncol == 1 && spendingTime(0, 0) == 0) {
    spendTime.fill(NaN);
  } else {
    if (spendingTime.ncol != m) {
      throw std::invalid_argument("Invalid number of columns for spendingTime");
    }
    if (spendingTime.nrow != k1 && spendingTime.nrow != B * k1) {
      throw std::invalid_argument("Invalid number of rows for spendingTime");
    }

    if (spendingTime.nrow == k1 && B > 1) {
      // replicate spendingTime for B iterations
      const double* src_ptr = spendingTime.data_ptr();
      double* dst_ptr = spendTime.data_ptr();
      for (int j = 0; j < m; ++ j) {
        for (int b = 0; b < B; ++b) {
          std::memcpy(dst_ptr + (b * k1) + j * (B * k1),
                      src_ptr + j * k1,
                      k1 * sizeof(double));
        }
      }
    } else {
      spendTime = spendingTime;
    }
  }


  std::vector<int> K0(m, 0); // max number of testable looks for hypothesis j
  for (int j = 0; j < m; ++j) {
    for (int k = 0; k < kMax; ++k) {
      if (incidenceMatrix(k, j)) K0[j]++;
    }
  }

  //--- study look index and testable look index, and number of testable looks
  IntMatrix idx1(k1, m); // study look index for each testable look
  IntMatrix idx2(m, k1); // testable look index for each study look
  idx1.fill(-1); // initialize with -1 for non-existent test look
  idx2.fill(-1); // initialize with -1 for non-tested study look
  std::vector<int> K1(m); // number of testable looks for hypothesis j at interim
  std::vector<int> K2(m); // last study look for each hypothesis at interim
  for (int j = 0; j < m; ++j) {
    int l = 0; // index for testable looks of hypothesis j
    for (int k = 0; k < k1; ++k) {
      if (incidenceMatrix(k, j)) {
        idx1(l, j) = k; // study look index for the l-th testable look of hypothesis j
        idx2(j, k) = l; // testable look index for the k-th study look of hypothesis j
        ++l;
      }
    }
    K1[j] = l; // number of testable looks for hypothesis j (1-based)
    if (l > 0) {
      K2[j] = idx1(l - 1, j) + 1;  // last study look for hypothesis j (1-based)
    }
  }

  // Implement the logic to determine which hypotheses are rejected at iteration b
  // based on the p-values, information, spending time, and the graphical procedure.
  auto f = [&](const int b)->std::vector<int> {
    std::vector<int> L(m);

    FlatMatrix p1(k1, m), t1(k1, m), s1(k1, m);
    p1.fill(NaN); t1.fill(NaN); s1.fill(NaN);
    double* p1_ptr = p1.data_ptr();
    double* t1_ptr = t1.data_ptr();
    double* s1_ptr = s1.data_ptr();

    // Extract p-values, information, and spending time vectors for each hypothesis
    std::vector<double> p_vec; p_vec.reserve(k1);
    std::vector<double> i_vec; i_vec.reserve(k1);
    std::vector<double> s_vec; s_vec.reserve(k1);
    for (int j = 0; j < m; ++j) { // loop over hypotheses
      int Kj = K1[j]; // number of testable looks for hypothesis j at interim
      if (Kj == 0) continue; // no testable look, skip to next hypothesis
      double maxinfoj = maxInformation[j]; // maxInformation for hypothesis j

      p_vec.resize(Kj);
      i_vec.resize(Kj);
      s_vec.resize(Kj);
      for (int l = 0; l < Kj; ++l) {
        int k = idx1(l, j); // study look index for the l-th testable look

        // row index in p, info, and spendTime for iteration b and study look k
        int h = b * k1 + k;
        if (std::isnan(p(h, j))) {
          throw std::invalid_argument("p must be provided at each testable look");
        }
        if (std::isnan(info(h, j))) {
          throw std::invalid_argument(
              "information must be provided at each testable look");
        }
        p_vec[l] = p(h, j);
        i_vec[l] = info(h, j);
        s_vec[l] = spendTime(h, j);
      }

      // Validate p values
      if (std::any_of(p_vec.begin(), p_vec.end(),
                      [](double v){ return v < 0.0 || v > 1.0; })) {
        throw std::invalid_argument("p must lie between 0 and 1");
      }

      // Validate information
      if (i_vec[0] <= 0.0) {
        throw std::invalid_argument("information must be positive");
      }
      if (any_nonincreasing(i_vec)) {
        throw std::invalid_argument("information must be increasing over time");
      }

      // Validate spending time
      bool all_na = std::all_of(s_vec.begin(), s_vec.end(),
                                [](double v){ return std::isnan(v); });
      if (!all_na) {
        if (std::any_of(s_vec.begin(), s_vec.end(),
                        [](double v){ return std::isnan(v); })) {
          throw std::invalid_argument(
              "spendingTime must be provided at each testable look");
        }
        if (s_vec[0] <= 0.0) {
          throw std::invalid_argument("spendingTime must be positive");
        }
        if (any_nonincreasing(s_vec)) {
          throw std::invalid_argument("spendingTime must be increasing over time");
        }
        if (s_vec[Kj - 1] > 1.0) {
          throw std::invalid_argument("spendingTime must be less than or equal to 1");
        }
      }

      // Determine L[j]: number of looks to consider for hypothesis j at interim
      if (all_na) { // will use information rates, no need to check s_vec further
        auto it = std::lower_bound(i_vec.begin(), i_vec.end(), maxinfoj);
        if (it == i_vec.end()) { // none >= maxInformation
          L[j] = Kj;
        } else { // first >= maxInformation, consider looks up to & including this one
          L[j] = static_cast<int>(std::distance(i_vec.begin(), it)) + 1;
        }
      } else { // will use spending time, consider all testable looks
        L[j] = Kj;
      }
      int Lj = L[j];

      // Copy p values in one block (p_vec contiguous)
      double* p1_col = p1_ptr + j * k1; // column j begins at offset j * k1
      double* t1_col = t1_ptr + j * k1;
      double* s1_col = s1_ptr + j * k1;

      std::memcpy(p1_col, p_vec.data(), Lj * sizeof(double));

      // Information time for forming covariance matrix of test statistics
      double info_L = i_vec[Lj - 1];
      for (int l = 0; l < Lj; ++l) {
        t1_col[l] = i_vec[l] / info_L;
      }

      // Spending time for error spending
      if (all_na) { // use information rates
        for (int l = 0; l < Lj; ++l) {
          if (l == K0[j] - 1 || i_vec[l] >= maxinfoj) {
            s1_col[l] = 1.0; // the last testable look is at or beyond maxInformation
          } else {
            s1_col[l] = i_vec[l] / maxinfoj;
          }
        }
      } else { // use spending time
        std::memcpy(s1_col, s_vec.data(), Lj * sizeof(double));
      }
    }

    int num_rejected = 0; // number of hypotheses rejected so far
    std::vector<int> reject(m, 0); // first look when the hypothesis is rejected
    std::vector<double> wx = w; // current weights for hypotheses, updated in-place
    FlatMatrix g = G; // current transition matrix, updated in-place
    FlatMatrix g1(m, m); // temporary transition matrix for update
    std::vector<double> user; // empty userAlphaSpending for getBoundcpp

    std::vector<int> active(m); // currently active hypotheses
    std::iota(active.begin(), active.end(), 0); // initialize with 0,1,...,m-1
    std::vector<int> pos(m); // position map: pos[idx] = index in active or -1
    std::iota(pos.begin(), pos.end(), 0); // initialize with 0,1,...,m-1

    std::vector<double> denom(m); // temp denominator for transition matrix update
    std::vector<double> u_vec; u_vec.reserve(k1); // upper bound from getBoundcpp

    // temp vectors reuse across hypotheses
    std::vector<double> t;
    std::vector<double> s;
    std::vector<unsigned char> x(k1, 1); // efficacyStopping for getBoundcpp
    std::vector<double> w_pre(k1); // previous weights

    // cached per-hypothesis previous upper bounds (u_pre[j] for hypothesis j)
    std::vector<std::vector<double>> u_pre(m);
    for (int j = 0; j < m; ++j) u_pre[j].resize(k1);

    // Preallocated helper vectors reused for "resuse" branch
    std::vector<double> l_vec(k1, -6.0);
    std::vector<double> theta_vec(k1, 0.0);

    // pointers to idx2
    const int* idx2_ptr = idx2.data_ptr(); // m x k1

    int K3 = *std::max_element(K2.begin(), K2.end());

    for (int step = 0; step < K3; ++step) {  // loop over study look
      const int* idx2_col = idx2_ptr + step * m;

      //Try to find a hypothesis that can be rejected at this step
      for ([[maybe_unused]] int i : active) {
        bool found_reject = false;
        int found_j = -1;

        // scan all hypotheses j to find a rejectable one
        for (int j : active) {
          if (wx[j] < 1e-8) continue;  // weight too small or already rejected
          int l = idx2_col[j];         // testable look index (0-based)
          if (l < 0) continue;         // not testable at this study look
          if (l >= L[j]) continue;     // beyond L[j] considered at interim
          int n = l + 1; // # of testable looks for hypothesis j at this study look

          double alpha1 = wx[j] * alpha;
          const std::string& asf1 = asf[j];
          double asfpar1 = asfpar[j];

          // Compute upper bound
          t = flatmatrix_get_column(t1, j);
          s = flatmatrix_get_column(s1, j);
          if (wx[j] != w_pre[j]) {
            // weights changed, compute full u_vec
            u_vec = getBoundcpp(n, t, alpha1, asf1, asfpar1, user, s, x);
          } else {
            // reuse previous u_pre[j] prefix, only solve for last element
            u_vec.resize(n);
            if (l > 0) std::memcpy(u_vec.data(), u_pre[j].data(), l * sizeof(double));
            // compute cumulative alpha for this n
            double cumAlpha = errorSpentcpp(s[l], alpha1, asf1, asfpar1);

            // small lambda that only sets last element
            auto g = [&](double aval)->double {
              u_vec[l] = aval;
              ListCpp probs = exitprobcpp(u_vec, l_vec, theta_vec, t);
              auto v = probs.get<std::vector<double>>("exitProbUpper");
              double cpu = std::accumulate(v.begin(), v.end(), 0.0);
              return cpu - cumAlpha;
            };

            double g_6 = g(6.0);
            if (g_6 > 0.0) { // no alpha spent at current visit
              u_vec[l] = 6.0;
            } else {
              auto g_for_brent = [&](double aval)->double {
                if (aval == 6.0) return g_6; // avoid recomputation at 6.0
                return g(aval);
              };
              u_vec[l] = brent(g_for_brent, -5.0, 6.0, 1e-6);
            }
          }

          // cache computed u_vec for hypothesis j
          std::memcpy(u_pre[j].data(), u_vec.data(), n * sizeof(double));

          // test rejection
          double alphastar = 1.0 - boost_pnorm(u_vec[l]);
          if (p1(l, j) < alphastar) {
            found_reject = true;
            found_j = j;
            break; // stop scanning j, we'll process rejection
          }
        } // end scan j

        if (!found_reject) { // no more rejections at this study look
          // remember current weights for next study look to check if weights changed
          w_pre = wx;
          break;
        }

        // Process rejection of hypothesis found_j
        int j = found_j;
        reject[j] = step + 1;
        ++num_rejected;

        // Remove j from active set
        int remove_pos = pos[j];
        int last_pos = static_cast<int>(active.size()) - 1;
        if (remove_pos != last_pos) {
          int swapped_idx = active[last_pos];
          active[remove_pos] = swapped_idx;
          pos[swapped_idx] = remove_pos;
        }
        active.pop_back();
        pos[j] = -1;

        // Update weights in-place for remaining active hypotheses
        w_pre = wx; // store current weights before update
        double wxj = wx[j];
        for (int l : active) {
          wx[l] += wxj * g(j, l);
        }
        wx[j] = 0.0; // weight of rejected hypothesis becomes 0

        // Update transition matrix for active hypotheses
        for (int l : active) {
          denom[l] =  1.0 - g(l, j) * g(j, l);
        }

        g1.fill(0.0); // reset g1 to 0 before filling
        for (int k : active) {
          double g_jk = g(j, k);
          for (int l : active) {
            if (l == k) continue;
            double dl = denom[l];
            if (dl > 1e-12) {
              g1(l, k) = (g(l, k) + g(l, j) * g_jk) / dl;
            }
          }
        }

        std::swap(g.data, g1.data);  // copy g1 to g for next iteration

        // Stop if all hypotheses rejected
        if (num_rejected == m) break;
      } // end attempt loop

      if (num_rejected == m) break;
    } // end study look loop

    return reject;
  };


  struct SimulationWorker : public RcppParallel::Worker {
    // references to read-only inputs (no mutation)
    const int m;
    const int k1;
    const std::vector<double>& w;
    const FlatMatrix& G;
    const double alpha;
    const std::vector<std::string>& asf;
    const std::vector<double>& asfpar;
    const std::vector<int>& K0;
    const std::vector<int>& K1;
    const std::vector<int>& K2;
    const IntMatrix& idx1;
    const IntMatrix& idx2;
    const std::vector<double>& maxInformation;
    const FlatMatrix& p;
    const FlatMatrix& info;
    const FlatMatrix& spendTime;

    // function f and other params that f needs are captured from outer scope
    // capture them by reference here so worker can call f(...)
    std::function<std::vector<int>(const int)> f;

    // result references (each iteration writes unique index into these)
    IntMatrix& reject_out; // B by m matrix to store rejection results

    // constructor
    SimulationWorker(const int m_,
                     const int k1_,
                     const std::vector<double>& w_,
                     const FlatMatrix& G_,
                     const double alpha_,
                     const std::vector<std::string>& asf_,
                     const std::vector<double>& asfpar_,
                     const std::vector<int>& K0_,
                     const std::vector<int>& K1_,
                     const std::vector<int>& K2_,
                     const IntMatrix& idx1_,
                     const IntMatrix& idx2_,
                     const std::vector<double>& maxInformation_,
                     const FlatMatrix& p_,
                     const FlatMatrix& info_,
                     const FlatMatrix& spendTime_,
                     decltype(f) f_,
                     IntMatrix& reject_out_) :

      m(m_), k1(k1_), w(w_), G(G_), alpha(alpha_),
      asf(asf_), asfpar(asfpar_), K0(K0_), K1(K1_), K2(K2_),
      idx1(idx1_), idx2(idx2_), maxInformation(maxInformation_),
      p(p_), info(info_), spendTime(spendTime_),
      f(std::move(f_)),
      reject_out(reject_out_) {}


    // operator() processes a range of bootstrap iterations [begin, end)
    void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t b = begin; b < end; ++b) {
        // call the (thread-safe) per-iteration function f
        std::vector<int> out = f(static_cast<int>(b));

        // write results
        for (int j = 0; j < m; ++j) {
          reject_out(b, j) = out[j];
        }
      } // end for b
    } // end operator()
  }; // end BootstrapWorker

  IntMatrix reject_mat(B, m); // B by m matrix to store rejection results

  // Instantiate the Worker with references to inputs and outputs
  SimulationWorker worker(
      m, k1, w, G, alpha, asf, asfpar, K0, K1, K2,
      idx1, idx2, maxInformation, p, info, spendTime,
      // bind f into std::function (capture the f we already have)
      std::function<std::vector<int>(const int)>(f),
      reject_mat
  );

  // Run the parallel loop over iterations
  RcppParallel::parallelFor(0, B, worker, 1 /*grain size*/);

  return worker.reject_out;
}


// [[Rcpp::export]]
Rcpp::IntegerMatrix fseqboncpp(
    const Rcpp::NumericVector& w,
    const Rcpp::NumericMatrix& G,
    const double alpha,
    const int kMax,
    const Rcpp::StringVector& typeAlphaSpending,
    const Rcpp::NumericVector& parameterAlphaSpending,
    const Rcpp::NumericVector& maxInformation,
    const Rcpp::LogicalMatrix& incidenceMatrix,
    const int k1,
    const Rcpp::NumericMatrix& p,
    const Rcpp::NumericMatrix& information,
    const Rcpp::NumericMatrix& spendingTime) {
  auto w1 = Rcpp::as<std::vector<double>>(w);
  auto G1 = flatmatrix_from_Rmatrix(G);
  auto asf1 = Rcpp::as<std::vector<std::string>>(typeAlphaSpending);
  auto asfpar1 = Rcpp::as<std::vector<double>>(parameterAlphaSpending);
  auto maxInfo1 = Rcpp::as<std::vector<double>>(maxInformation);
  auto incid1 = boolmatrix_from_Rmatrix(incidenceMatrix);
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto info1 = flatmatrix_from_Rmatrix(information);
  auto spendTime1 = flatmatrix_from_Rmatrix(spendingTime);
  auto reject1 = fseqboncpp1(w1, G1, alpha, kMax, asf1, asfpar1, maxInfo1,
                             incid1, k1, p1, info1, spendTime1);
  return Rcpp::wrap(reject1);
}



// Converts step-down p-values to sequential adjusted p-values
FlatMatrix fstp2seqcpp1(
    const FlatMatrix& p,
    const std::vector<double>& gamma,
    const std::string& test = "hochberg",
    const bool retest = true) {

  // Validate dimensions
  int nreps = p.nrow;
  int n = p.ncol;
  int m = n / 2;

  // Normalize test string to lowercase
  std::string test1 = test;
  for (char &c : test1) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (test1 != "hochberg" && test1 != "holm") {
    throw std::invalid_argument("test must be 'hochberg' or 'holm'");
  }

  // Validate p-values in [0, 1]
  if (std::any_of(p.data.begin(), p.data.end(),
                  [](double v){ return v < 0.0 || v > 1.0; })) {
    throw std::invalid_argument("p-values must be between 0 and 1");
  }

  bool is_holm = (test1 == "holm");

  // Step 1: Compute matrix a (nreps × m) from pairs of p-values
  FlatMatrix a(nreps, m);
  for (int j = 0; j < m; ++j) {
    double gamma_j = gamma[j];
    double inv_factor = 1.0 / (1.0 + gamma_j);

    for (int iter = 0; iter < nreps; ++iter) {
      double x1 = p(iter, 2*j);
      double x2 = p(iter, 2*j+1);
      double val = 2.0 * std::max(x1, x2) * inv_factor;
      if (is_holm) {
        val = std::max(val, 2.0 * std::min(x1, x2));
      }
      a(iter, j) = val;
    }
  }

  // Step 2: Compute matrix d (m × m) - lower triangular discount factors
  FlatMatrix d(m, m);
  for (int s = 0; s < m; ++s) {
    double gmax = 0.0;
    for (int j = s; j < m; ++j) {
      if (j > s) {
        gmax = std::max(gmax, gamma[j - 1]);
      }
      d(s, j) = std::max((1.0 - gmax) * 0.5, 1e-12);
    }
  }

  // Step 3: Compute adjusted p-values
  FlatMatrix padj(nreps, n);

  // Preallocate reusable vectors
  std::vector<double> a_row(m);        // current row of a
  std::vector<double> a_cummax(m);     // cumulative max
  for (int iter = 0; iter < nreps; ++iter) {
    // Extract row iter from a and compute cumulative max
    for (int j = 0; j < m; ++j) {
      a_row[j] = a(iter, j);
    }

    a_cummax[0] = a_row[0];
    for (int j = 1; j < m; ++j) {
      a_cummax[j] = std::max(a_cummax[j - 1], a_row[j]);
    }

    // For each hypothesis index i
    for (int i = 0; i < m; ++i) {
      double t1 = a_cummax[i]; // max of a[0..i]

      // --- Compute padj(iter, 2*i) (even index) ---
      double t2_even = 1.0;
      for (int s = 0; s <= i; ++s) {
        double lhs = (s > 0) ? a_cummax[s - 1] : 0.0;

        // max over j ∈ [s, i] of p(iter, 2*j) / d(s, j)
        for (int j = s; j <= i; ++j) {
          lhs = std::max(lhs, p(iter, 2*j) / d(s, j));
        }

        double rhs = 2.0 * p(iter, 2*s+1) / (1.0 + gamma[s]);
        if (lhs < rhs) {
          t2_even = std::min(t2_even, lhs);
        }
      }

      double result_even;
      if (retest && m > 1) {
        double t3_even = 1.0;
        int s_max = std::min(i, m - 2);

        for (int s = 0; s <= s_max; ++s) {
          double lhs = (s > 0) ? a_cummax[s - 1] : 0.0;

          // max over j ∈ [s, m-1] of p(iter, 2*j+1) / d(s, j)
          for (int j = s; j < m; ++j) {
            lhs = std::max(lhs, p(iter, 2*j+1) / d(s, j));
          }

          // max over j ∈ [s, i] of p(iter, 2*j)
          for (int j = s; j <= i; ++j) {
            lhs = std::max(lhs, p(iter, 2*j));
          }

          double rhs = 2.0 * p(iter, 2*s) / (1.0 + gamma[s]);
          if (lhs < rhs) {
            t3_even = std::min(t3_even, lhs);
          }
        }

        result_even = std::min({t1, t2_even, t3_even});
      } else {
        result_even = std::min(t1, t2_even);
      }

      padj(iter, 2*i) = result_even;

      // --- Compute padj(iter, 2*i+1) (odd index) ---
      double t2_odd = 1.0;
      for (int s = 0; s <= i; ++s) {
        double lhs = (s > 0) ? a_cummax[s - 1] : 0.0;

        // max over j ∈ [s, i] of p(iter, 2*j+1) / d(s, j)
        for (int j = s; j <= i; ++j) {
          lhs = std::max(lhs, p(iter, 2*j+1) / d(s, j));
        }

        double rhs = 2.0 * p(iter, 2*s) / (1.0 + gamma[s]);
        if (lhs < rhs) {
          t2_odd = std::min(t2_odd, lhs);
        }
      }

      double result_odd;
      if (retest && m > 1) {
        double t3_odd = 1.0;
        int s_max = std::min(i, m - 2);

        for (int s = 0; s <= s_max; ++s) {
          double lhs = (s > 0) ? a_cummax[s - 1] : 0.0;

          // max over j ∈ [s, m-1] of p(iter, 2*j) / d(s, j)
          for (int j = s; j < m; ++j) {
            lhs = std::max(lhs, p(iter, 2*j) / d(s, j));
          }

          // max over j ∈ [s, i] of p(iter, 2*j+1)
          for (int j = s; j <= i; ++j) {
            lhs = std::max(lhs, p(iter, 2*j+1));
          }

          double rhs = 2.0 * p(iter, 2*s+1) / (1.0 + gamma[s]);
          if (lhs < rhs) {
            t3_odd = std::min(t3_odd, lhs);
          }
        }

        result_odd = std::min({t1, t2_odd, t3_odd});
      } else {
        result_odd = std::min(t1, t2_odd);
      }

      padj(iter, 2*i+1) = result_odd;
    }
  }

  return padj;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix fstp2seqcpp(
    const Rcpp::NumericMatrix& p,
    const Rcpp::NumericVector& gamma,
    const std::string& test = "hochberg",
    const bool retest = true) {

  auto p1 = flatmatrix_from_Rmatrix(p);
  auto gamma1 = Rcpp::as<std::vector<double>>(gamma);
  auto padj1 = fstp2seqcpp1(p1, gamma1, test, retest);
  return Rcpp::wrap(padj1);
}


// Helper: compute row sums for BoolMatrix
std::vector<int> boolmatrix_rowsums(const BoolMatrix& mat) {
  int nrow = mat.nrow;
  int ncol = mat.ncol;
  std::vector<int> sums(nrow, 0);
  const unsigned char* data = mat.data_ptr();

  for (int j = 0; j < ncol; ++j) {
    const unsigned char* col = data + j * nrow;
    for (int i = 0; i < nrow; ++i) {
      sums[i] += col[i];
    }
  }
  return sums;
}


// Helper for adjusted p-values for standard mixture gatekeeping procedures
FlatMatrix fstdmixcpp1(
    const FlatMatrix& p,
    const BoolMatrix& family,
    const BoolMatrix& serial,
    const BoolMatrix& parallel,
    const std::vector<double>& gamma,
    const std::string& test = "hommel",
    const bool exhaust = true) {

  // Validate inputs
  int nreps = p.nrow;
  int m = p.ncol;

  if (family.ncol != m) {
    throw std::invalid_argument("family must have m columns");
  }

  if (serial.nrow != m || serial.ncol != m) {
    throw std::invalid_argument("serial must be m x m matrix");
  }

  if (parallel.nrow != m || parallel.ncol != m) {
    throw std::invalid_argument("parallel must be m x m matrix");
  }

  int nfamily = family.nrow;

  if (gamma.size() != static_cast<std::size_t>(nfamily)) {
    throw std::invalid_argument("gamma length must equal nrow(family)");
  }

  // Validate gamma in [0, 1]
  if (std::any_of(gamma.begin(), gamma.end(),
                  [](double g){ return g < 0.0 || g > 1.0; })) {
    throw std::invalid_argument("gamma values must be between 0 and 1");
  }

  // Validate p-values in [0, 1]
  if (std::any_of(p.data.begin(), p.data.end(),
                  [](double v){ return v < 0.0 || v > 1.0; })) {
    throw std::invalid_argument("p-values must be between 0 and 1");
  }

  // Normalize test string
  std::string test1 = test;
  for (char &c : test1) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (test1 != "hommel" && test1 != "hochberg" && test1 != "holm") {
    throw std::invalid_argument("test must be 'hommel', 'hochberg', or 'holm'");
  }

  // Compute number of hypotheses per family
  std::vector<int> nhyps = boolmatrix_rowsums(family);

  // Number of intersection tests
  int ntests = (1 << m) - 1; // 2^m - 1

  // Matrix to store local p-values for intersection tests
  FlatMatrix pinter(nreps, ntests);

  // Incidence matrix for elementary hypotheses
  BoolMatrix incid(ntests, m);

  // Process each intersection hypothesis
  for (int i = 0; i < ntests; ++i) {
    int number = ntests - i;

    // Binary representation of intersection hypothesis
    std::vector<unsigned char> cc(m);
    for (int j = 0; j < m; ++j) {
      cc[j] = (number >> (m - 1 - j)) & 1;
    }

    // Active hypotheses in each family, I1, ..., I_nfamily
    BoolMatrix family0(nfamily, m);
    for (int j = 0; j < m; ++j) {
      for (int h = 0; h < nfamily; ++h) {
        family0(h, j) = family(h, j) & cc[j];
      }
    }

    // Number of active hypotheses by family, k1, ..., k_nfamily
    std::vector<int> nhyps0 = boolmatrix_rowsums(family0);

    // Determine restricted index set for each family
    std::vector<unsigned char> cc1(m, 1);

    for (int j = 0; j < m; ++j) {
      // Check serial constraints
      int serial_sum = 0;
      for (int k = 0; k < m; ++k) {
        serial_sum += serial(j, k);
      }

      if (serial_sum > 0) {
        // if any serial predecessor is accepted, remove j
        for (int k = 0; k < m; ++k) {
          if (serial(j, k) && cc[k]) {
            cc1[j] = 0;
            break;
          }
        }

        // if any serial predecessor is not testable, remove j
        for (int k = 0; k < m; ++k) {
          if (serial(j, k) && !cc1[k]) {
            cc1[j] = 0;
            break;
          }
        }
      }

      // Check parallel constraints
      int parallel_sum = 0;
      for (int k = 0; k < m; ++k) {
        parallel_sum += parallel(j, k);
      }

      if (parallel_sum > 0) {
        // if none of the parallel predecessors are rejected, remove j
        bool hit = true;
        for (int k = 0; k < m; ++k) {
          if (parallel(j, k) && !cc[k]) {
            hit = false;
            break;
          }
        }
        if (hit) cc1[j] = 0;

        // if none of the parallel predecessors are testable, remove j
        hit = true;
        for (int k = 0; k < m; ++k) {
          if (parallel(j, k) && cc1[k]) {
            hit = false;
            break;
          }
        }
        if (hit) cc1[j] = 0;
      }
    }

    // Apply intersection
    for (int j = 0; j < m; ++j) {
      cc1[j] = cc1[j] & cc[j];
    }

    // Error rate function divided by alpha
    std::vector<double> errf(nfamily);
    for (int j = 0; j < nfamily; ++j) {
      errf[j] = nhyps0[j] > 0 ?
      gamma[j] + (1.0 - gamma[j]) * nhyps0[j] / nhyps[j] : 0.0;
    }

    // Allocated fraction of alpha for each family
    std::vector<double> coef(nfamily);
    coef[0] = 1.0;
    for (int j = 1; j < nfamily; ++j) {
      coef[j] = coef[j - 1] * (1.0 - errf[j - 1]);
    }

    int kmax = 0; // last family with positive allocation (1-based index)
    for (int j = nfamily - 1; j >= 0; --j) {
      if (coef[j] > 0.0) {
        kmax = j + 1;
        break;
      }
    }

    // Families up to kmax, I_1*, ..., I_kmax*
    BoolMatrix family1(kmax, m);
    for (int j = 0; j < kmax; ++j) {
      for (int k = 0; k < m; ++k) {
        family1(j, k) = family(j, k) & cc1[k];
      }
    }

    // Number of testable hypotheses by family, k_1*, ..., k_kmax*
    std::vector<int> nhyps1 = boolmatrix_rowsums(family1);

    // Indices of active families
    std::vector<int> sub;
    sub.reserve(kmax);
    for (int j = 0; j < kmax; ++j) {
      if (nhyps1[j] > 0) {
        sub.push_back(j);
      }
    }

    int nfamily2 = static_cast<int>(sub.size());

    // Subset of active families after removing those without testable hypotheses
    BoolMatrix family2(nfamily2, m);
    std::vector<int> nhyps2(nfamily2);
    for (int j = 0; j < nfamily2; ++j) {
      for (int k = 0; k < m; ++k) {
        family2(j, k) = family1(sub[j], k);
      }
      nhyps2[j] = nhyps1[sub[j]];
    }

    // family indices and hypothesis indices for testable hypotheses
    std::vector<int> fam, hyps2;
    for (int j = 0; j < nfamily2; ++j) {
      for (int k = 0; k < m; ++k) {
        if (family2(j, k)) {
          fam.push_back(j);
          hyps2.push_back(k);
        }
      }
    }

    // total number of testable hypotheses across active families
    int n = static_cast<int>(hyps2.size());

    // Relative importance for active families
    std::vector<double> coef1(nfamily2);
    for (int j = 0; j < nfamily2; ++j) {
      coef1[j] = coef[sub[j]];
    }

    // Broadcasted family weight for each testable hypothesis
    std::vector<double> c(n);
    for (int k = 0; k < n; ++k) {
      c[k] = coef1[fam[k]];
    }

    // Truncation parameters
    std::vector<double> gam2(nfamily2);
    for (int j = 0; j < nfamily2; ++j) {
      gam2[j] = gamma[sub[j]];
    }
    if (exhaust) gam2[nfamily2 - 1] = 1.0;

    // Bonferroni part of weights
    std::vector<double> coef2(nfamily2);
    for (int j = 0; j < nfamily2; ++j) {
      coef2[j] = (1.0 - gam2[j]) / nhyps[sub[j]];
    }

    // Broadcasted Bonferroni part of weights to each testable hypothesis
    std::vector<double> tbon(n);
    for (int k = 0; k < n; ++k) {
      tbon[k] = coef2[fam[k]];
    }

    // Cumulative count of hypotheses before the current family
    std::vector<int> ck(nfamily2 + 1, 0);
    for (int j = 1; j <= nfamily2; ++j) {
      ck[j] = ck[j - 1] + nhyps2[j - 1];
    }

    // Compute weights
    std::vector<double> w(n);
    if (test1 == "hommel") {
      for (int k = 0; k < n; ++k) {
        int l = fam[k];
        int j = (k + 1) - ck[l];
        w[k] = j * gam2[l] / nhyps2[l] + tbon[k];
      }
    } else if (test1 == "hochberg") {
      for (int k = 0; k < n; ++k) {
        int l = fam[k];
        int j = (k + 1) - ck[l];
        w[k] = gam2[l] / (nhyps2[l] - j + 1) + tbon[k];
      }
    } else { // holm
      for (int k = 0; k < n; ++k) {
        int l = fam[k];
        w[k] = gam2[l] / nhyps2[l] + tbon[k];
      }
    }

    // Process each replication
    std::vector<double> p1(n), p1s(n), p2(n);
    for (int iter = 0; iter < nreps; ++iter) {
      // Extract raw p-values
      for (int k = 0; k < n; ++k) {
        p1[k] = p(iter, hyps2[k]);
      }

      // Sort p-values within each family
      for (int j = 0; j < nfamily2; ++j) {
        int start = ck[j];
        int end = ck[j + 1];
        int len = end - start;
        p1s = subset(p1, start, end);
        std::sort(p1s.begin(), p1s.end());
        std::memcpy(p2.data() + start, p1s.data(), len * sizeof(double));
      }

      // Compute minimum ratio
      double min_val = 1.0;
      for (int k = 0; k < n; ++k) {
        double ratio = p2[k] / (w[k] * c[k]);
        min_val = std::min(min_val, ratio);
      }

      pinter(iter, i) = min_val;
    }

    // Store incidence
    for (int j = 0; j < m; ++j) {
      incid(i, j) = cc[j];
    }
  }

  // Compute adjusted p-values for elementary hypotheses
  FlatMatrix padj(nreps, m);
  for (int j = 0; j < m; ++j) {
    for (int iter = 0; iter < nreps; ++iter) {
      double max_p = 0.0;
      for (int i = 0; i < ntests; ++i) {
        if (incid(i, j)) {
          max_p = std::max(max_p, pinter(iter, i));
        }
      }
      padj(iter, j) = std::min(max_p, 1.0);
    }
  }

  // Apply logical restrictions (serial constraints)
  for (int j = 0; j < m; ++j) {
    int serial_sum = 0;
    for (int k = 0; k < m; ++k) {
      serial_sum += serial(j, k);
    }

    if (serial_sum > 0) {
      for (int iter = 0; iter < nreps; ++iter) {
        double pre = 0.0;
        for (int k = 0; k < m; ++k) {
          if (serial(j, k)) {
            pre = std::max(pre, padj(iter, k));
          }
        }
        padj(iter, j) = std::max(padj(iter, j), pre);
      }
    }
  }

  // Apply logical restrictions (parallel constraints)
  for (int j = 0; j < m; ++j) {
    int parallel_sum = 0;
    for (int k = 0; k < m; ++k) {
      parallel_sum += parallel(j, k);
    }

    if (parallel_sum > 0) {
      for (int iter = 0; iter < nreps; ++iter) {
        double pre = 1.0;
        for (int k = 0; k < m; ++k) {
          if (parallel(j, k)) {
            pre = std::min(pre, padj(iter, k));
          }
        }
        padj(iter, j) = std::max(padj(iter, j), pre);
      }
    }
  }

  return padj;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix fstdmixcpp(
    const Rcpp::NumericMatrix& p,
    const Rcpp::LogicalMatrix& family,
    const Rcpp::LogicalMatrix& serial,
    const Rcpp::LogicalMatrix& parallel,
    const Rcpp::NumericVector& gamma,
    const std::string& test = "hommel",
    const bool exhaust = true) {
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto family1 = boolmatrix_from_Rmatrix(family);
  auto serial1 = boolmatrix_from_Rmatrix(serial);
  auto parallel1 = boolmatrix_from_Rmatrix(parallel);
  auto gamma1 = Rcpp::as<std::vector<double>>(gamma);
  auto padj1 = fstdmixcpp1(p1, family1, serial1, parallel1, gamma1, test, exhaust);
  return Rcpp::wrap(padj1);
}


// Helper for adjusted p-values for modified mixture gatekeeping procedures
FlatMatrix fmodmixcpp1(
    const FlatMatrix& p,
    const BoolMatrix& family,
    const BoolMatrix& serial,
    const BoolMatrix& parallel,
    const std::vector<double>& gamma,
    const std::string& test = "hommel",
    const bool exhaust = true) {

  // Validate inputs
  int nreps = p.nrow;
  int m = p.ncol;

  if (family.ncol != m) {
    throw std::invalid_argument("family must have m columns");
  }

  if (serial.nrow != m || serial.ncol != m) {
    throw std::invalid_argument("serial must be m x m matrix");
  }

  if (parallel.nrow != m || parallel.ncol != m) {
    throw std::invalid_argument("parallel must be m x m matrix");
  }

  int nfamily = family.nrow;

  if (gamma.size() != static_cast<std::size_t>(nfamily)) {
    throw std::invalid_argument("gamma length must equal nrow(family)");
  }

  // Validate gamma in [0, 1]
  if (std::any_of(gamma.begin(), gamma.end(),
                  [](double g){ return g < 0.0 || g > 1.0; })) {
    throw std::invalid_argument("gamma values must be between 0 and 1");
  }

  // Validate p-values in [0, 1]
  if (std::any_of(p.data.begin(), p.data.end(),
                  [](double v){ return v < 0.0 || v > 1.0; })) {
    throw std::invalid_argument("p-values must be between 0 and 1");
  }

  // Normalize test string
  std::string test1 = test;
  for (char &c : test1) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (test1 != "hommel" && test1 != "hochberg" && test1 != "holm") {
    throw std::invalid_argument("test must be 'hommel', 'hochberg', or 'holm'");
  }

  // Compute number of hypotheses per family
  std::vector<int> nhyps = boolmatrix_rowsums(family);

  // Number of intersection tests
  int ntests = (1 << m) - 1; // 2^m - 1

  // Matrix to store local p-values for intersection tests
  FlatMatrix pinter(nreps, ntests);

  // Incidence matrix for elementary hypotheses
  BoolMatrix incid(ntests, m);

  // Process each intersection hypothesis
  for (int i = 0; i < ntests; ++i) {
    int number = ntests - i;

    // Binary representation of intersection hypothesis
    std::vector<unsigned char> cc(m);
    for (int j = 0; j < m; ++j) {
      cc[j] = (number >> (m - 1 - j)) & 1;
    }

    // Active hypotheses in each family, I1, ..., I_nfamily
    BoolMatrix family0(nfamily, m);
    for (int j = 0; j < m; ++j) {
      for (int h = 0; h < nfamily; ++h) {
        family0(h, j) = family(h, j) & cc[j];
      }
    }

    // Number of active hypotheses by family, k1, ..., k_nfamily
    std::vector<int> nhyps0 = boolmatrix_rowsums(family0);

    // Determine restricted index set for each family
    std::vector<unsigned char> cc1(m, 1);

    for (int j = 0; j < m; ++j) {
      // Check serial constraints
      int serial_sum = 0;
      for (int k = 0; k < m; ++k) {
        serial_sum += serial(j, k);
      }

      if (serial_sum > 0) {
        // if any serial predecessor is accepted, remove j
        for (int k = 0; k < m; ++k) {
          if (serial(j, k) && cc[k]) {
            cc1[j] = 0;
            break;
          }
        }

        // if any serial predecessor is not testable, remove j
        for (int k = 0; k < m; ++k) {
          if (serial(j, k) && !cc1[k]) {
            cc1[j] = 0;
            break;
          }
        }
      }

      // Check parallel constraints
      int parallel_sum = 0;
      for (int k = 0; k < m; ++k) {
        parallel_sum += parallel(j, k);
      }

      if (parallel_sum > 0) {
        // if none of the parallel predecessors are rejected, remove j
        bool hit = true;
        for (int k = 0; k < m; ++k) {
          if (parallel(j, k) && !cc[k]) {
            hit = false;
            break;
          }
        }
        if (hit) cc1[j] = 0;

        // if none of the parallel predecessors are testable, remove j
        hit = true;
        for (int k = 0; k < m; ++k) {
          if (parallel(j, k) && cc1[k]) {
            hit = false;
            break;
          }
        }
        if (hit) cc1[j] = 0;
      }
    }

    std::vector<unsigned char> cc2 = cc1; // for nstar calculation

    // Apply intersection, for kstar calculation
    for (int j = 0; j < m; ++j) {
      cc1[j] = cc1[j] & cc[j];
    }

    // Compute kstar and nstar for error rate function
    std::vector<int> kstar(nfamily), nstar(nfamily);
    std::vector<double> errf(nfamily);
    for (int j = 0; j < nfamily; ++j) {
      kstar[j] = 0;
      nstar[j] = 0;

      for (int k = 0; k < m; ++k) {
        if (family(j, k) && cc1[k]) kstar[j]++;
        if (family(j, k) && cc2[k]) nstar[j]++;
      }

      // KEY DIFFERENCE from stdmix: uses kstar and nstar instead of nhyps0 and nhyps
      errf[j] = kstar[j] > 0 ?
      gamma[j] + (1.0 - gamma[j]) * kstar[j] / nstar[j] : 0.0;
    }

    // Allocated fraction of alpha for each family
    std::vector<double> coef(nfamily);
    coef[0] = 1.0;
    for (int j = 1; j < nfamily; ++j) {
      coef[j] = coef[j - 1] * (1.0 - errf[j - 1]);
    }

    int kmax = 0; // last family with positive allocation (1-based index)
    for (int j = nfamily - 1; j >= 0; --j) {
      if (coef[j] > 0.0) {
        kmax = j + 1;
        break;
      }
    }

    // Families up to kmax, I_1*, ..., I_kmax*
    BoolMatrix family1(kmax, m);
    for (int j = 0; j < kmax; ++j) {
      for (int k = 0; k < m; ++k) {
        family1(j, k) = family(j, k) & cc1[k];
      }
    }

    // Number of testable hypotheses by family, k_1*, ..., k_kmax*
    std::vector<int> nhyps1 = boolmatrix_rowsums(family1);

    // Indices of active families
    std::vector<int> sub;
    sub.reserve(kmax);
    for (int j = 0; j < kmax; ++j) {
      if (nhyps1[j] > 0) {
        sub.push_back(j);
      }
    }

    int nfamily2 = static_cast<int>(sub.size());

    // Subset of active families after removing those without testable hypotheses
    BoolMatrix family2(nfamily2, m);
    std::vector<int> nhyps2(nfamily2);
    for (int j = 0; j < nfamily2; ++j) {
      for (int k = 0; k < m; ++k) {
        family2(j, k) = family1(sub[j], k);
      }
      nhyps2[j] = nhyps1[sub[j]];
    }

    // family indices and hypothesis indices for testable hypotheses
    std::vector<int> fam, hyps2;
    for (int j = 0; j < nfamily2; ++j) {
      for (int k = 0; k < m; ++k) {
        if (family2(j, k)) {
          fam.push_back(j);
          hyps2.push_back(k);
        }
      }
    }

    // total number of testable hypotheses across active families
    int n = static_cast<int>(hyps2.size());

    // Relative importance for active families
    std::vector<double> coef1(nfamily2);
    for (int j = 0; j < nfamily2; ++j) {
      coef1[j] = coef[sub[j]];
    }

    // Broadcasted family weight for each testable hypothesis
    std::vector<double> c(n);
    for (int k = 0; k < n; ++k) {
      c[k] = coef1[fam[k]];
    }

    // Truncation parameters
    std::vector<double> gam2(nfamily2);
    for (int j = 0; j < nfamily2; ++j) {
      gam2[j] = gamma[sub[j]];
    }
    if (exhaust) gam2[nfamily2 - 1] = 1.0;

    // Bonferroni part of weights
    // (KEY DIFFERENCE from stdmix: uses nstar instead of nhyps)
    std::vector<double> coef2(nfamily2);
    for (int j = 0; j < nfamily2; ++j) {
      coef2[j] = (1.0 - gam2[j]) / nstar[sub[j]];
    }

    // Broadcasted Bonferroni part of weights to each testable hypothesis
    std::vector<double> tbon(n);
    for (int k = 0; k < n; ++k) {
      tbon[k] = coef2[fam[k]];
    }

    // Cumulative count of hypotheses before the current family
    std::vector<int> ck(nfamily2 + 1, 0);
    for (int j = 1; j <= nfamily2; ++j) {
      ck[j] = ck[j - 1] + nhyps2[j - 1];
    }

    // Compute weights
    std::vector<double> w(n);
    if (test1 == "hommel") {
      for (int k = 0; k < n; ++k) {
        int l = fam[k];
        int j = (k + 1) - ck[l];
        w[k] = j * gam2[l] / nhyps2[l] + tbon[k];
      }
    } else if (test1 == "hochberg") {
      for (int k = 0; k < n; ++k) {
        int l = fam[k];
        int j = (k + 1) - ck[l];
        w[k] = gam2[l] / (nhyps2[l] - j + 1) + tbon[k];
      }
    } else { // holm
      for (int k = 0; k < n; ++k) {
        int l = fam[k];
        w[k] = gam2[l] / nhyps2[l] + tbon[k];
      }
    }

    // Process each replication
    std::vector<double> p1(n), p1s(n), p2(n);
    for (int iter = 0; iter < nreps; ++iter) {
      // Extract raw p-values
      for (int k = 0; k < n; ++k) {
        p1[k] = p(iter, hyps2[k]);
      }

      // Sort p-values within each family
      for (int j = 0; j < nfamily2; ++j) {
        int start = ck[j];
        int end = ck[j + 1];
        int len = end - start;
        p1s = subset(p1, start, end);
        std::sort(p1s.begin(), p1s.end());
        std::memcpy(p2.data() + start, p1s.data(), len * sizeof(double));
      }

      // Compute minimum ratio
      double min_val = 1.0;
      for (int k = 0; k < n; ++k) {
        double ratio = p2[k] / (w[k] * c[k]);
        min_val = std::min(min_val, ratio);
      }

      pinter(iter, i) = min_val;
    }

    // Store incidence
    for (int j = 0; j < m; ++j) {
      incid(i, j) = cc[j];
    }
  }

  // Compute adjusted p-values for elementary hypotheses
  FlatMatrix padj(nreps, m);
  for (int j = 0; j < m; ++j) {
    for (int iter = 0; iter < nreps; ++iter) {
      double max_p = 0.0;
      for (int i = 0; i < ntests; ++i) {
        if (incid(i, j)) {
          max_p = std::max(max_p, pinter(iter, i));
        }
      }
      padj(iter, j) = std::min(max_p, 1.0);
    }
  }

  // Apply logical restrictions (serial constraints)
  for (int j = 0; j < m; ++j) {
    int serial_sum = 0;
    for (int k = 0; k < m; ++k) {
      serial_sum += serial(j, k);
    }

    if (serial_sum > 0) {
      for (int iter = 0; iter < nreps; ++iter) {
        double pre = 0.0;
        for (int k = 0; k < m; ++k) {
          if (serial(j, k)) {
            pre = std::max(pre, padj(iter, k));
          }
        }
        padj(iter, j) = std::max(padj(iter, j), pre);
      }
    }
  }

  // Apply logical restrictions (parallel constraints)
  for (int j = 0; j < m; ++j) {
    int parallel_sum = 0;
    for (int k = 0; k < m; ++k) {
      parallel_sum += parallel(j, k);
    }

    if (parallel_sum > 0) {
      for (int iter = 0; iter < nreps; ++iter) {
        double pre = 1.0;
        for (int k = 0; k < m; ++k) {
          if (parallel(j, k)) {
            pre = std::min(pre, padj(iter, k));
          }
        }
        padj(iter, j) = std::max(padj(iter, j), pre);
      }
    }
  }

  return padj;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix fmodmixcpp(
    const Rcpp::NumericMatrix& p,
    const Rcpp::LogicalMatrix& family,
    const Rcpp::LogicalMatrix& serial,
    const Rcpp::LogicalMatrix& parallel,
    const Rcpp::NumericVector& gamma,
    const std::string& test = "hommel",
    const bool exhaust = true) {
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto family1 = boolmatrix_from_Rmatrix(family);
  auto serial1 = boolmatrix_from_Rmatrix(serial);
  auto parallel1 = boolmatrix_from_Rmatrix(parallel);
  auto gamma1 = Rcpp::as<std::vector<double>>(gamma);
  auto padj1 = fmodmixcpp1(p1, family1, serial1, parallel1, gamma1, test, exhaust);
  return Rcpp::wrap(padj1);
}


// Helper to compute truncated adjusted p-values for multiple testing
FlatMatrix ftrunccpp1(
    const FlatMatrix& p,
    const std::string& test,
    const double gamma) {

  // Validate inputs
  int niters = p.nrow;
  int m = p.ncol;

  if (gamma < 0.0 || gamma > 1.0) {
    throw std::invalid_argument("gamma must be between 0 and 1");
  }

  // Validate p-values in [0, 1]
  if (std::any_of(p.data.begin(), p.data.end(),
                  [](double v){ return v < 0.0 || v > 1.0; })) {
    throw std::invalid_argument("p-values must be between 0 and 1");
  }

  // Normalize test string
  std::string test1 = test;
  for (char &c : test1) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (test1 != "hommel" && test1 != "hochberg" && test1 != "holm") {
    throw std::invalid_argument("test must be 'hommel', 'hochberg', or 'holm'");
  }

  int ntests = (1 << m) - 1; // 2^m - 1

  // Incidence matrix (ntests x m) - which hypotheses are in each intersection
  BoolMatrix incid(ntests, m);

  // Intersection p-values (niters x ntests)
  FlatMatrix pinter(niters, ntests);

  // Precompute constant factors
  double tbon = (1.0 - gamma) / m;

  // Preallocate reusable vectors
  std::vector<int> hyp_indices; hyp_indices.reserve(m);
  std::vector<double> p1; p1.reserve(m);

  // Process each intersection hypothesis
  for (int i = 0; i < ntests; ++i) {
    int number = ntests - i;

    // Binary representation of elementary hypotheses in intersection
    std::vector<unsigned char> cc(m);
    int k = 0; // count of hypotheses in this intersection
    hyp_indices.clear();

    for (int j = 0; j < m; ++j) {
      cc[j] = (number >> (m - 1 - j)) & 1;
      if (cc[j]) {
        hyp_indices.push_back(j);
        ++k;
      }
    }

    // Store incidence (column-major: column j base = incid_data + j * ntests)
    for (int j = 0; j < m; ++j) {
      incid(i, j) = cc[j];
    }

    // Precompute weights if they're constant across iterations
    std::vector<double> weights(k);
    if (test1 == "hommel") {
      for (int j = 0; j < k; ++j) {
        weights[j] = (j + 1) * gamma / k + tbon;
      }
    } else if (test1 == "hochberg") {
      for (int j = 0; j < k; ++j) {
        weights[j] = gamma / (k - j) + tbon;
      }
    } else { // holm
      double w_holm = gamma / k + tbon;
      std::fill_n(weights.data(), weights.size(), w_holm);
    }

    // Process each iteration/replication
    for (int iter = 0; iter < niters; ++iter) {
      // Extract p-values for hypotheses in this intersection
      p1.resize(k);
      for (int j = 0; j < k; ++j) {
        p1[j] = p(iter, hyp_indices[j]);
      }

      // Sort p-values
      std::sort(p1.begin(), p1.end());

      // Compute minimum ratio
      double q = 1.0;
      for (int j = 0; j < k; ++j) {
        double ratio = p1[j] / weights[j];
        q = std::min(q, ratio);
      }

      pinter(iter, i) = q;
    }
  }

  // Compute adjusted p-values for individual hypotheses
  FlatMatrix padj(niters, m);

  // For each hypothesis j, find max pinter over intersections containing j
  for (int j = 0; j < m; ++j) {
    for (int iter = 0; iter < niters; ++iter) {
      double max_p = 0.0;
      for (int i = 0; i < ntests; ++i) {
        if (incid(i, j)) {
          max_p = std::max(max_p, pinter(iter, i));
        }
      }

      padj(iter, j) = max_p;
    }
  }

  return padj;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix ftrunccpp(
    const Rcpp::NumericMatrix& p,
    const std::string& test,
    const double gamma) {
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto padj1 = ftrunccpp1(p1, test, gamma);
  return Rcpp::wrap(padj1);
}
