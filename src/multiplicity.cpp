#include "multiplicity.h"
#include "dataframe_list.h"
#include "generic_design.h"
#include "mvnormr.h"
#include "utilities.h"

#include <algorithm> // any_of, distance, fill, max_element, min, sort
#include <cctype>    // tolower
#include <cmath>     // fabs, isnan
#include <cstring>   // memcpy
#include <numeric>   // accumulate
#include <stdexcept> // invalid_argument
#include <string>    // string
#include <vector>    // vector

#include <Rcpp.h>
#include <RcppParallel.h>

using std::size_t;

// Helper to update graph for graphical approaches
Graph updateGraphcpp(const std::vector<double> &w, const FlatMatrix &G,
                     const std::vector<size_t> &I, const size_t j) {

  size_t m = w.size();

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
  for (size_t i = 0; i < m * m; ++i) {
    if (G.data[i] < 0.0) {
      throw std::invalid_argument("G must be nonnegative");
    }
  }

  // Validation: Row sums of G must be <= 1
  std::vector<double> rowsum(m, 0.0);
  for (size_t j = 0; j < m; ++j) {
    for (size_t i = 0; i < m; ++i) {
      rowsum[i] += G(i, j);
    }
  }
  for (size_t i = 0; i < m; ++i) {
    if (rowsum[i] > 1.0 + 1e-8) {
      throw std::invalid_argument(
          "Row sums of G must be less than or equal to 1");
    }
  }

  // Validation: Diagonal elements of G must be 0
  for (size_t i = 0; i < m; ++i) {
    if (G(i, i) != 0.0) {
      throw std::invalid_argument("Diagonal elements of G must be equal to 0");
    }
  }

  // Check validity of elements of I, create zero-based indexing, and remove j
  std::vector<unsigned char> seen(m, 0); // 0/1 marks presence
  std::vector<size_t> I_new;
  I_new.reserve(I.size()); // new I after removing j
  for (size_t i : I) {
    if (seen[i]) {
      throw std::invalid_argument(
          "The index set I must not contain duplicates.");
    }
    seen[i] = 1;
    if (i != j)
      I_new.push_back(i); // new I excludes index j
  }

  // Check that j belongs to I
  if (!seen[j])
    throw std::invalid_argument("j must be in I.");

  // Update weights
  const double w_j = w[j];
  std::vector<double> wx = w; // copy w
  for (size_t i : I_new) {
    wx[i] += w_j * G(j, i);
  }
  wx[j] = 0.0;

  // Update transition matrix
  std::vector<double> denom(m);
  for (size_t l : I_new) {
    denom[l] = 1.0 - G(l, j) * G(j, l);
  }

  FlatMatrix g(m, m);
  for (size_t k : I_new) {
    double g_jk = G(j, k);
    for (size_t l : I_new) {
      if (l == k)
        continue;
      double dl = denom[l];
      if (dl > 1e-12) {
        g(l, k) = (G(l, k) + G(l, j) * g_jk) / dl;
      }
    }
  }

  return Graph{std::move(wx), std::move(g), std::move(I_new)};
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
Rcpp::List updateGraph(const std::vector<double> &w,
                       const Rcpp::NumericMatrix &G, const std::vector<int> &I,
                       const int j) {

  int m = w.size();
  auto I1 = validateConvertIdx(I, m, "I");
  if (j <= 0)
    throw std::invalid_argument("j must be positive");
  size_t j1 = j - 1;
  auto G1 = flatmatrix_from_Rmatrix(G);
  auto out = updateGraphcpp(w, G1, I1, j1);
  for (size_t &i : out.I)
    ++i;
  ListCpp result;
  result.push_back(std::move(out.w), "w");
  result.push_back(std::move(out.G), "G");
  result.push_back(std::move(out.I), "I");
  return Rcpp::wrap(result);
}

// Helper to compute the default weight matrix for graphical approaches
WeightMatrix fDefaultWgtmatcpp(size_t m) {
  const size_t ntests = (size_t{1} << m) - 1;
  WeightMatrix out;
  out.inthyp.resize(ntests, m);
  out.wgtmat.resize(ntests, m);

  for (size_t i = 0; i < ntests; ++i) {
    const size_t number = ntests - i;
    size_t nactive = 0;

    for (size_t k = 0; k < m; ++k) {
      const int bit = static_cast<int>((number >> (m - 1 - k)) & size_t{1});
      out.inthyp(i, k) = bit;
      nactive += static_cast<size_t>(bit);
    }

    const double weight = 1.0 / static_cast<double>(nactive);
    for (size_t k = 0; k < m; ++k) {
      if (out.inthyp(i, k)) {
        out.wgtmat(i, k) = weight;
      }
    }
  }

  return out;
}

BoolMatrix fDefaultFamilycpp(const size_t m) {
  BoolMatrix family(1, m);
  family.fill(1);
  return family;
}

FlatMatrix fDefaultCorrcpp(const BoolMatrix &family) {
  const size_t m = family.ncol;
  FlatMatrix corr(m, m);
  corr.fill(NaN);

  for (size_t j = 0; j < m; ++j) {
    corr(j, j) = 1.0;
    for (size_t i = 0; i < j; ++i) {
      bool same_family = false;
      for (size_t h = 0; h < family.nrow; ++h) {
        if (family(h, i) && family(h, j)) {
          same_family = true;
          break;
        }
      }
      if (same_family) {
        corr(i, j) = 0.5;
        corr(j, i) = 0.5;
      }
    }
  }

  return corr;
}

// Helper to compute the full weight matrix for graphical approaches
WeightMatrix fwgtmatcpp(const std::vector<double> &w, const FlatMatrix &G) {
  size_t m = w.size();
  size_t ntests = (1 << m) - 1;
  size_t gtr_nrow = (ntests + 1) / 2;
  size_t gtr_ncol = m * m;

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
  for (size_t i = 0; i < m * m; ++i) {
    if (G.data[i] < 0.0) {
      throw std::invalid_argument("G must be nonnegative");
    }
  }

  // Validation: Row sums of G must be <= 1
  std::vector<double> rowsum(m, 0.0);
  for (size_t j = 0; j < m; ++j) {
    for (size_t i = 0; i < m; ++i) {
      rowsum[i] += G(i, j);
    }
  }
  for (size_t i = 0; i < m; ++i) {
    if (rowsum[i] > 1.0 + 1e-8) {
      throw std::invalid_argument(
          "Row sums of G must be less than or equal to 1");
    }
  }

  // Validation: Diagonal elements of G must be 0
  for (size_t i = 0; i < m; ++i) {
    if (G(i, i) != 0.0) {
      throw std::invalid_argument("Diagonal elements of G must be equal to 0");
    }
  }

  // Preallocations and initial copies
  std::vector<double> wx = w; // dynamic weights, reused
  FlatMatrix g = G;           // mutable transition matrix
  FlatMatrix g1(m, m); // temp transition matrix, reused (zeroed when needed)
  FlatMatrix gtrmat(gtr_nrow, gtr_ncol); // first half of transition matrix
  IntMatrix inthyp(ntests, m);           // store intersection hypotheses
  FlatMatrix wgtmat(ntests, m); // store weights of elementary hypotheses

  std::vector<size_t> active; // indices of active hypotheses in intersection
  active.reserve(m);

  // bitmask for m bits, used to flip bits and find super set
  const int mask = (1 << m) - 1;
  std::vector<double> denom(m); // temp denominator for transition matrix update

  for (size_t i = 0; i < ntests; ++i) {
    if (i >= 1) {
      int number = ntests - i; // original mapping

      // Build list of active indices
      active.clear();
      size_t j = 0; // index of minimum active hypothesis
      bool found_zero = false;
      for (size_t k = 0; k < m; ++k) {
        int bit = (number >> (m - 1 - k)) & 1;
        inthyp(i, k) = bit;
        if (bit)
          active.push_back(k);
        else if (!found_zero) {
          j = k;
          found_zero = true;
        }
      }
      if (!found_zero)
        j = 0;

      // index of the super set, with j-th bit set to 0 and others flipped
      size_t ip = ((~number) & mask) & ~(1 << (m - 1 - j));

      // Load the weights from the super set
      for (size_t k = 0; k < m; ++k) {
        wx[k] = wgtmat(ip, k);
      }

      // Load the transition matrix from the super set
      for (size_t k = 0; k < m; ++k) {
        for (size_t l = 0; l < m; ++l) {
          g(l, k) = gtrmat(ip, k * m + l);
        }
      }

      // Update the weights
      double wxj = wx[j];
      if (wxj != 0.0) {
        for (size_t k : active) {
          wx[k] += wxj * g(j, k);
        }
        wx[j] = 0.0;
      }

      // Update the transition matrix
      for (size_t l : active) {
        denom[l] = 1.0 - g(l, j) * g(j, l);
      }

      g1.fill(0.0); // ensure g1 is zeroed before use
      for (size_t k : active) {
        double g_jk = g(j, k);
        for (size_t l : active) {
          if (l == k)
            continue;
          double dl = denom[l];
          if (dl > 1e-12) {
            g1(l, k) = (g(l, k) + g(l, j) * g_jk) / dl;
          }
        }
      }

      std::swap(g.data, g1.data); // copy g1 to g for next iteration
    } else {
      for (size_t k = 0; k < m; ++k) {
        inthyp(i, k) = 1;
      }
    }

    // Save the weights
    for (size_t k = 0; k < m; ++k) {
      wgtmat(i, k) = wx[k];
    }

    // Save the transition matrix
    if (i < gtr_nrow) {
      for (size_t k = 0; k < m; ++k) {
        for (size_t l = 0; l < m; ++l) {
          gtrmat(i, k * m + l) = g(l, k);
        }
      }
    }
  }

  return WeightMatrix{inthyp, wgtmat};
}

//' @title Default Weight Matrix for All Intersection Hypotheses
//' @description Obtains the default weight matrix for all intersection
//' hypotheses, assigning equal weights to the elementary hypotheses within
//' each intersection hypothesis.
//'
//' @param m The number of elementary hypotheses.
//'
//' @return A list with the following components:
//' * \code{inthyp}: The indicator matrix for the intersection hypotheses.
//' * \code{wgtmat}: The default weight matrix for the elementary hypotheses.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' fDefaultWgtmat(3)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fDefaultWgtmat(const size_t m) {
  auto out = fDefaultWgtmatcpp(m);
  ListCpp result;
  result.push_back(out.inthyp, "inthyp");
  result.push_back(out.wgtmat, "wgtmat");
  return Rcpp::wrap(result);
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
//' w <- c(0.5,0.5,0,0)
//' G <- matrix(c(0,0,1,0, 0,0,0,1, 0,1,0,0, 1,0,0,0),
//'             nrow=4, ncol=4, byrow=TRUE)
//' (wgtmat <- fwgtmat(w,G))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fwgtmat(const Rcpp::NumericVector &w, const Rcpp::NumericMatrix &G) {
  auto w1 = Rcpp::as<std::vector<double>>(w);
  auto G1 = flatmatrix_from_Rmatrix(G);
  auto out = fwgtmatcpp(w1, G1);
  ListCpp result;
  result.push_back(out.inthyp, "inthyp");
  result.push_back(out.wgtmat, "wgtmat");
  return Rcpp::wrap(result);
}

// Helper to compute adjusted p-values for Bonferroni-based graphical approaches
AdjustedPValues fadjpboncpp(const FlatMatrix &p, const WeightMatrix &wgtmat) {

  size_t ntests = wgtmat.inthyp.nrow;
  size_t m = p.ncol;
  size_t niters = p.nrow;

  // Output matrix initialized to 0
  FlatMatrix padj(niters, m);
  FlatMatrix pinter(niters, ntests);
  pinter.fill(1.0); // initialize to 1.0 for min operation

  // Main loop over subsets
  for (size_t i = 0; i < ntests; ++i) {
    double *pinter_col = &pinter.data[i * niters]; // contiguous column

    for (size_t j = 0; j < m; ++j) {
      double w = wgtmat.wgtmat(i, j);
      if (w == 0 || wgtmat.inthyp(i, j) == 0)
        continue;
      for (size_t iter = 0; iter < niters; ++iter) {
        double v = p(iter, j) / w;
        if (v < pinter_col[iter])
          pinter_col[iter] = v;
      }
    }

    // Update padj columns for active hypotheses (each hyp[t] gets max over
    // subsets)
    for (size_t j = 0; j < m; ++j) {
      if (wgtmat.inthyp(i, j) == 0)
        continue;
      double *padj_col = &padj.data[j * niters]; // contiguous column
      for (size_t iter = 0; iter < niters; ++iter) {
        double v = pinter_col[iter];
        if (v > padj_col[iter])
          padj_col[iter] = v;
      }
    }
  } // end subsets

  return AdjustedPValues{wgtmat.inthyp, pinter, padj};
}

AdjustedPValues fadjpboncpp(const FlatMatrix &p) {
  return fadjpboncpp(p, fDefaultWgtmatcpp(p.ncol));
}

// [[Rcpp::export]]
Rcpp::List fadjpbonRcpp(const Rcpp::NumericMatrix &p,
                        const Rcpp::Nullable<Rcpp::List> &wgtmat = R_NilValue) {
  auto p1 = flatmatrix_from_Rmatrix(p);
  if (wgtmat.isNull()) {
    auto out = fadjpboncpp(p1);
    ListCpp result;
    result.push_back(std::move(out.inthyp), "inthyp");
    result.push_back(std::move(out.pinter), "pinter");
    result.push_back(std::move(out.padj), "padj");
    return Rcpp::wrap(result);
  }
  Rcpp::List wgtmat_list(wgtmat.get());
  ListPtr wgtmat_ptr = listcpp_from_rlist(wgtmat_list);
  WeightMatrix wgt_pair;
  wgt_pair.inthyp = wgtmat_ptr->get<IntMatrix>("inthyp");
  wgt_pair.wgtmat = wgtmat_ptr->get<FlatMatrix>("wgtmat");
  auto out = fadjpboncpp(p1, wgt_pair);
  ListCpp result;
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.pinter), "pinter");
  result.push_back(std::move(out.padj), "padj");
  return Rcpp::wrap(result);
}

// Helper to compute adjusted p-values for Simes-based graphical approaches
AdjustedPValues fadjpsimcpp(const FlatMatrix &p, const WeightMatrix &wgtmat,
                            const BoolMatrix &family) {

  size_t ntests = wgtmat.inthyp.nrow;
  size_t m = p.ncol;
  size_t niters = p.nrow;
  size_t nfams = family.nrow;

  if (family.ncol != m) {
    throw std::invalid_argument(
        "family must have as many elementary hypotheses as columns");
  }

  // Output matrix initialized to 0
  FlatMatrix padj(niters, m);
  FlatMatrix pinter(niters, ntests);

  // Reusable buffers sized up to m
  std::vector<double> wbuf;
  wbuf.reserve(m);
  std::vector<double> pbuf;
  pbuf.reserve(m);
  std::vector<double> cw;
  cw.reserve(m);
  std::vector<size_t> idx;
  idx.reserve(m);
  std::vector<double> p1;
  p1.reserve(m);
  std::vector<double> w1;
  w1.reserve(m);

  // Main loop over subsets
  for (size_t i = 0; i < ntests; ++i) {
    std::vector<size_t> nhyps0(nfams, 0);
    std::vector<size_t> hyp;
    for (size_t j = 0; j < m; ++j) {
      if (!wgtmat.inthyp(i, j))
        continue;
      for (size_t k = 0; k < nfams; ++k) {
        if (family(k, j)) {
          nhyps0[k]++;
          hyp.push_back(j);
          break; // a hypothesis belongs to exactly one family
        }
      }
    }
    std::vector<size_t> nhyps1;
    for (size_t k = 0; k < nfams; ++k) {
      if (nhyps0[k] > 0)
        nhyps1.push_back(nhyps0[k]);
    }

    size_t nhyps = hyp.size();

    // Extract weights wx for this subset from wgtmat row i (column-major)
    wbuf.assign(nhyps, 0.0);
    for (size_t t = 0; t < nhyps; ++t) {
      size_t col = hyp[t];
      wbuf[t] = wgtmat.wgtmat(i, col);
    }

    double *pinter_col = &pinter.data[i * niters]; // contiguous column
    // For each iteration (row of p), compute pinter(iter, i)
    for (size_t iter = 0; iter < niters; ++iter) {
      // Extract p-values for active hypotheses in the same hyp order
      pbuf.resize(nhyps);
      cw.assign(nhyps, 0.0);
      for (size_t t = 0; t < nhyps; ++t) {
        size_t col = hyp[t];
        pbuf[t] = p(iter, col);
      }

      // Sort p-values within each active family block
      // (nhyps1 gives block sizes in family order).
      size_t s = 0;
      for (size_t block = 0; block < nhyps1.size(); ++block) {
        size_t t = nhyps1[block];
        // snapshot original block to avoid in-place overwrite corruption
        // p1 and w1 are reusable buffers (reserved to m, resized per block)
        p1.resize(t);
        w1.resize(t); // within family block
        for (size_t u = 0; u < t; ++u) {
          p1[u] = pbuf[s + u];
          w1[u] = wbuf[s + u];
        }

        // obtain index ordering of the snapshot
        idx = seqcpp(0, t - 1);
        std::sort(idx.begin(), idx.end(),
                  [&p1](size_t a, size_t b) { return p1[a] < p1[b]; });

        // copy sorted p and compute cumulative weights from the snapshot
        double cum = 0.0;
        for (size_t j = 0; j < t; ++j) {
          size_t src = idx[j];
          double pv = p1[src];
          double wv = w1[src];
          cum += wv;
          pbuf[s + j] = pv; // write sorted p-values into pbuf
          cw[s + j] = cum;  // cumulative weights
        }
        s += t;
      }

      // compute q = min_j pbuf[j] / cw[j] ignoring cw==0
      double q = 1.0;
      for (size_t j = 0; j < nhyps; ++j) {
        double cj = cw[j];
        if (cj > 0.0) {
          double ratio = pbuf[j] / cj;
          if (ratio < q)
            q = ratio;
        }
      }
      pinter_col[iter] = q;
    } // end iter loop

    // Update padj columns for active hypotheses (each hyp[t] gets max over
    // subsets)
    for (size_t t = 0; t < nhyps; ++t) {
      size_t col = hyp[t];
      double *padj_col = &padj.data[col * niters]; // contiguous column
      for (size_t iter = 0; iter < niters; ++iter) {
        double v = pinter_col[iter];
        if (v > padj_col[iter])
          padj_col[iter] = v;
      }
    }
  } // end subsets

  return AdjustedPValues{wgtmat.inthyp, pinter, padj};
}

AdjustedPValues fadjpsimcpp(const FlatMatrix &p, const BoolMatrix &family) {
  return fadjpsimcpp(p, fDefaultWgtmatcpp(p.ncol), family);
}

AdjustedPValues fadjpsimcpp(const FlatMatrix &p, const WeightMatrix &wgtmat) {
  return fadjpsimcpp(p, wgtmat, fDefaultFamilycpp(p.ncol));
}

// [[Rcpp::export]]
Rcpp::List
fadjpsimRcpp(const Rcpp::NumericMatrix &p,
             const Rcpp::Nullable<Rcpp::List> &wgtmat = R_NilValue,
             const Rcpp::Nullable<Rcpp::LogicalMatrix> &family = R_NilValue) {
  auto p1 = flatmatrix_from_Rmatrix(p);
  BoolMatrix family1 = fDefaultFamilycpp(p1.ncol);
  if (!family.isNull()) {
    Rcpp::LogicalMatrix family_matrix(family.get());
    family1 = boolmatrix_from_Rmatrix(family_matrix);
  }
  if (wgtmat.isNull()) {
    auto out = fadjpsimcpp(p1, family1);
    ListCpp result;
    result.push_back(std::move(out.inthyp), "inthyp");
    result.push_back(std::move(out.pinter), "pinter");
    result.push_back(std::move(out.padj), "padj");
    return Rcpp::wrap(result);
  }
  Rcpp::List wgtmat_list(wgtmat.get());
  ListPtr wgtmat_ptr = listcpp_from_rlist(wgtmat_list);
  WeightMatrix wgt_pair;
  wgt_pair.inthyp = wgtmat_ptr->get<IntMatrix>("inthyp");
  wgt_pair.wgtmat = wgtmat_ptr->get<FlatMatrix>("wgtmat");
  auto out = fadjpsimcpp(p1, wgt_pair, family1);
  ListCpp result;
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.pinter), "pinter");
  result.push_back(std::move(out.padj), "padj");
  return Rcpp::wrap(result);
}

// Helper to compute adjusted p-values for Dunnett-based graphical approaches
AdjustedPValues fadjpduncpp(const FlatMatrix &p, const WeightMatrix &wgtmat,
                            const BoolMatrix &family, const FlatMatrix &corr) {

  size_t ntests = wgtmat.inthyp.nrow;
  size_t m = p.ncol;
  size_t niters = p.nrow;
  size_t nfams = family.nrow;

  if (family.ncol != m) {
    throw std::invalid_argument(
        "family must have as many elementary hypotheses as columns");
  }

  // Output matrix initialized to 0
  FlatMatrix padj(niters, m);
  FlatMatrix pinter(niters, ntests);

  std::vector<size_t> J, J_h;
  std::vector<std::vector<size_t>> J_h_list(nfams);
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> mean;
  FlatMatrix sigma;
  std::vector<std::vector<double>> lower_list(nfams);
  std::vector<std::vector<double>> mean_list(nfams);
  std::vector<FlatMatrix> sigma_list(nfams);

  // Main loop over subsets
  for (size_t i = 0; i < ntests; ++i) {
    double *pinter_col = &pinter.data[i * niters]; // contiguous column

    J.clear();
    for (size_t j = 0; j < m; ++j) {
      if (wgtmat.wgtmat(i, j) > 0.0) {
        J.push_back(j);
      }
    }

    for (size_t h = 0; h < nfams; ++h) {
      J_h.clear();
      for (size_t j : J) {
        if (family(h, j))
          J_h.push_back(j);
      }
      J_h_list[h] = J_h;

      size_t k = J_h.size();
      if (k > 0) {
        lower.resize(k);
        lower.assign(k, -POS_INF);
        lower_list[h] = lower;

        mean.resize(k);
        mean.assign(k, 0.0);
        mean_list[h] = mean;

        sigma.resize(k, k);
        for (size_t t1 = 0; t1 < k; ++t1) {
          for (size_t t2 = 0; t2 < k; ++t2) {
            size_t j1 = J_h[t1];
            size_t j2 = J_h[t2];
            sigma(t1, t2) = corr(j1, j2);
          }
        }
        sigma_list[h] = sigma;
      }
    }

    // For each iteration (row of p), compute pinter(iter, i)
    for (size_t iter = 0; iter < niters; ++iter) {
      double aval = 1.0;
      for (size_t h = 0; h < nfams; ++h) {
        J_h = J_h_list[h];
        size_t k = J_h.size();
        if (k > 0) {
          double sumw = 0.0;
          double q = 1.0;
          for (size_t j : J_h) {
            sumw += wgtmat.wgtmat(i, j);
            double ratio = p(iter, j) / wgtmat.wgtmat(i, j);
            if (ratio < q)
              q = ratio;
          }

          const auto &lowerref = lower_list[h];
          const auto &meanref = mean_list[h];
          const auto &sigmaref = sigma_list[h];

          upper.resize(k);
          for (size_t t = 0; t < k; ++t) {
            size_t j = J_h[t];
            upper[t] = boost_qnorm(1.0 - wgtmat.wgtmat(i, j) * q);
          }

          double v1 = pmvnormcpp(lowerref, upper, meanref, sigmaref).prob;
          double v2 = (1.0 - v1) / sumw;
          if (v2 < aval)
            aval = v2;
        }
      }

      pinter_col[iter] = aval;
    } // end iter loop

    // Update padj columns for active hypotheses
    for (size_t j = 0; j < m; ++j) {
      if (!wgtmat.inthyp(i, j))
        continue;
      double *padj_col = &padj.data[j * niters]; // contiguous column
      for (size_t iter = 0; iter < niters; ++iter) {
        double v = pinter_col[iter];
        if (v > padj_col[iter])
          padj_col[iter] = v;
      }
    }
  } // end subsets

  return AdjustedPValues{wgtmat.inthyp, pinter, padj};
}

AdjustedPValues fadjpduncpp(const FlatMatrix &p, const BoolMatrix &family,
                            const FlatMatrix &corr) {
  return fadjpduncpp(p, fDefaultWgtmatcpp(p.ncol), family, corr);
}

AdjustedPValues fadjpduncpp(const FlatMatrix &p, const WeightMatrix &wgtmat,
                            const BoolMatrix &family) {
  return fadjpduncpp(p, wgtmat, family, fDefaultCorrcpp(family));
}

AdjustedPValues fadjpduncpp(const FlatMatrix &p, const BoolMatrix &family) {
  return fadjpduncpp(p, fDefaultWgtmatcpp(p.ncol), family,
                     fDefaultCorrcpp(family));
}

AdjustedPValues fadjpduncpp(const FlatMatrix &p, const WeightMatrix &wgtmat) {
  const BoolMatrix family = fDefaultFamilycpp(p.ncol);
  return fadjpduncpp(p, wgtmat, family, fDefaultCorrcpp(family));
}

AdjustedPValues fadjpduncpp(const FlatMatrix &p) {
  return fadjpduncpp(p, fDefaultWgtmatcpp(p.ncol));
}

// [[Rcpp::export]]
Rcpp::List
fadjpdunRcpp(const Rcpp::NumericMatrix &p,
             const Rcpp::Nullable<Rcpp::List> &wgtmat = R_NilValue,
             const Rcpp::Nullable<Rcpp::LogicalMatrix> &family = R_NilValue,
             const Rcpp::Nullable<Rcpp::NumericMatrix> &corr = R_NilValue) {
  auto p1 = flatmatrix_from_Rmatrix(p);
  BoolMatrix family1 = fDefaultFamilycpp(p1.ncol);
  if (!family.isNull()) {
    Rcpp::LogicalMatrix family_matrix(family.get());
    family1 = boolmatrix_from_Rmatrix(family_matrix);
  }
  FlatMatrix corr1 = fDefaultCorrcpp(family1);
  if (!corr.isNull()) {
    Rcpp::NumericMatrix corr_matrix(corr.get());
    corr1 = flatmatrix_from_Rmatrix(corr_matrix);
  }
  if (wgtmat.isNull()) {
    auto out = fadjpduncpp(p1, family1, corr1);
    ListCpp result;
    result.push_back(std::move(out.inthyp), "inthyp");
    result.push_back(std::move(out.pinter), "pinter");
    result.push_back(std::move(out.padj), "padj");
    return Rcpp::wrap(result);
  }
  Rcpp::List wgtmat_list(wgtmat.get());
  ListPtr wgtmat_ptr = listcpp_from_rlist(wgtmat_list);
  WeightMatrix wgt_pair;
  wgt_pair.inthyp = wgtmat_ptr->get<IntMatrix>("inthyp");
  wgt_pair.wgtmat = wgtmat_ptr->get<FlatMatrix>("wgtmat");
  auto out = fadjpduncpp(p1, wgt_pair, family1, corr1);
  ListCpp result;
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.pinter), "pinter");
  result.push_back(std::move(out.padj), "padj");
  return Rcpp::wrap(result);
}

// Helper to compute repeated p-values for alpha spending approaches
FlatMatrix repeatedPValuecpp(const size_t kMax,
                             const std::string &typeAlphaSpending,
                             const double parameterAlphaSpending,
                             const double maxInformation, const FlatMatrix &p,
                             const FlatMatrix &information,
                             const FlatMatrix &spendingTime) {

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
    throw std::invalid_argument(
        "parameterAlphaSpending must be positive for sfKD");
  }

  // Validation: maxInformation
  if (maxInformation <= 0) {
    throw std::invalid_argument("maxInformation must be positive");
  }

  // Validation: dimensions of p and information
  size_t B = p.nrow;
  size_t k1 = p.ncol;

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

  const double *in_ptr = information.data_ptr();
  double *info_ptr = info.data_ptr();
  if (information.nrow == B) {
    // shapes match -> fast block copy
    std::memcpy(info_ptr, in_ptr, B * k1 * sizeof(double));
  } else {
    // information.nrow == 1 -> broadcast the single row to all B rows
    // For each column k, value = in_ptr[k], destination column is info_ptr +
    // k*B
    for (size_t k = 0; k < k1; ++k) {
      double v = in_ptr[k];               // source scalar for column k
      double *dst_col = info_ptr + k * B; // column-major base
      std::fill_n(dst_col, B, v);         // efficient bulk write
    }
  }

  // Process spendingTime matrix with possible broadcasting and NA handling
  FlatMatrix spendTime(B, k1);
  if (spendingTime.nrow == 1 && spendingTime.ncol == 1 &&
      spendingTime(0, 0) == 0) {
    spendTime.fill(NaN);
  } else {
    if (spendingTime.ncol != k1) {
      throw std::invalid_argument("Invalid number of columns for spendingTime");
    }
    if (spendingTime.nrow != 1 && spendingTime.nrow != B) {
      throw std::invalid_argument("Invalid number of rows for spendingTime");
    }

    const double *sp_ptr = spendingTime.data_ptr();
    double *spend_ptr = spendTime.data_ptr();
    if (spendingTime.nrow == B) {
      // shapes match -> fast block copy
      std::memcpy(spend_ptr, sp_ptr, B * k1 * sizeof(double));
    } else {
      // spendingTime.nrow == 1 -> broadcast the single row to all B rows
      // For each column k, value = sp_ptr[k], destination column is spend_ptr +
      // k*B
      for (size_t k = 0; k < k1; ++k) {
        double v = sp_ptr[k];                // source scalar for column k
        double *dst_col = spend_ptr + k * B; // column-major base
        std::fill_n(dst_col, B, v);          // efficient bulk write
      }
    }
  }

  auto f = [&](const size_t b) -> std::vector<double> {
    std::vector<double> p_vec(k1), i_vec(k1), s_vec(k1);
    for (size_t k = 0; k < k1; ++k) {
      p_vec[k] = p(b, k);
      i_vec[k] = info(b, k);
      s_vec[k] = spendTime(b, k);
    }

    // Validation: p-values must be between 0 and 1 and not NA
    if (std::any_of(p_vec.begin(), p_vec.end(),
                    [](double v) { return std::isnan(v); })) {
      throw std::invalid_argument("p-values must be provided");
    }
    if (std::any_of(p_vec.begin(), p_vec.end(),
                    [](double v) { return v < 0 || v > 1; })) {
      throw std::invalid_argument("p-values must be between 0 and 1");
    }

    // Validation: information must be positive, increasing, and not NA
    if (std::any_of(i_vec.begin(), i_vec.end(),
                    [](double v) { return std::isnan(v); })) {
      throw std::invalid_argument("information must be provided");
    }
    if (i_vec[0] <= 0.0) {
      throw std::invalid_argument("information must be positive");
    }
    if (any_nonincreasing(i_vec)) {
      throw std::invalid_argument("information must be increasing over time");
    }

    // Validation: spendingTime must be positive, increasing, and not exceeding
    // 1
    bool all_na = std::all_of(s_vec.data(), s_vec.data() + k1,
                              [](double v) { return std::isnan(v); });

    if (!all_na) {
      if (std::any_of(s_vec.data(), s_vec.data() + k1,
                      [](double v) { return std::isnan(v); })) {
        throw std::invalid_argument("spendingTime must be provided");
      }
      if (s_vec[0] <= 0.0) {
        throw std::invalid_argument("spendingTime must be positive");
      }
      if (any_nonincreasing(s_vec)) {
        throw std::invalid_argument(
            "spendingTime must be increasing over time");
      }
      if (s_vec[k1 - 1] > 1.0) {
        throw std::invalid_argument("spendingTime must not exceed 1");
      }
    }

    // Determine L based on maxInformation and spendingTime
    size_t L;
    if (all_na) { // use information rates
      // Find if any information >= maxInformation
      auto it = std::lower_bound(i_vec.begin(), i_vec.end(), maxInformation);
      if (it == i_vec.end()) { // none >= maxInformation
        L = k1;
      } else {
        L = static_cast<size_t>(std::distance(i_vec.begin(), it)) + 1;
      }
    } else { // use spending time
      L = k1;
    }

    std::vector<double> p1(k1), t1(k1), s1(k1);
    std::memcpy(p1.data(), p_vec.data(), L * sizeof(double));

    // Information time for forming covariance matrix of test statistics
    double info_L = i_vec[L - 1];
    for (size_t l = 0; l < L; ++l) {
      t1[l] = i_vec[l] / info_L;
    }

    // Spending time for error spending
    if (all_na) { // use information rates
      for (size_t l = 0; l < L; ++l) {
        if (l == kMax - 1 || i_vec[l] >= maxInformation) {
          s1[l] = 1.0; // the last look is at or beyond maxInformation
        } else {
          s1[l] = i_vec[l] / maxInformation;
        }
      }
    } else { // using spending time
      std::memcpy(s1.data(), s_vec.data(), L * sizeof(double));
    }

    // Compute repeated p-values
    std::vector<double> empty_user;
    std::vector<double> repp(k1);
    for (size_t i = 0; i < L; ++i) {
      double pvalue = p1[i];
      std::vector<double> t0(t1.begin(), t1.begin() + i + 1);
      std::vector<double> s0(s1.begin(), s1.begin() + i + 1);
      std::vector<unsigned char> x(i + 1, 1);

      BoundCacheAlpha cache(i + 1, t0, asf, parameterAlphaSpending, empty_user,
                            s0, x, 64, 12);

      // Lambda function for root finding to solve for the repeated p-value at
      // step i
      auto g = [&](double a) -> double {
        auto u = cache.get(a);
        return 1.0 - boost_pnorm(u[i]) - pvalue;
      };

      // Find root in (0, 1) using brent's method, with checks at the endpoints
      double gl = g(0.000001);
      if (gl > 0) {
        repp[i] = 0.000001;
      } else {
        double gh = g(0.999999);
        if (gh < 0) {
          repp[i] = 0.999999;
        } else {
          repp[i] = brent(g, 0.000001, 0.999999, 1e-6);
        }
      }
    }

    return repp;
  };

  struct SimulationWorker : public RcppParallel::Worker {
    // references to read-only inputs (no mutation)
    const size_t k1;
    const std::string &asf;
    const double parameterAlphaSpending;
    const size_t kMax;
    const double maxInformation;
    const FlatMatrix &p;
    const FlatMatrix &info;
    const FlatMatrix &spendTime;

    // function f and other params that f needs are captured from outer scope
    // capture them by reference here so worker can call f(...)
    std::function<std::vector<double>(const size_t)> f;

    // result references (each iteration writes unique index into these)
    FlatMatrix &repp_out; // B by k1 matrix to store repeated p-values

    // constructor
    SimulationWorker(const size_t k1_, const std::string &asf_,
                     const double parameterAlphaSpending_, const size_t kMax_,
                     const double maxInformation_, const FlatMatrix &p_,
                     const FlatMatrix &info_, const FlatMatrix &spendTime_,
                     decltype(f) f_, FlatMatrix &repp_out_)
        :

          k1(k1_), asf(asf_), parameterAlphaSpending(parameterAlphaSpending_),
          kMax(kMax_), maxInformation(maxInformation_), p(p_), info(info_),
          spendTime(spendTime_), f(std::move(f_)), repp_out(repp_out_) {}

    // operator() processes a range of bootstrap iterations [begin, end)
    void operator()(size_t begin, size_t end) {
      for (size_t b = begin; b < end; ++b) {
        // call the (thread-safe) per-iteration function f
        std::vector<double> out = f(b);

        // write results
        for (size_t k = 0; k < k1; ++k) {
          repp_out(b, k) = out[k];
        }
      } // end for b
    } // end operator()
  }; // end BootstrapWorker

  FlatMatrix repp_mat(B, k1); // B by k1 matrix to store repeated values

  // Instantiate the Worker with references to inputs and outputs
  SimulationWorker worker(
      k1, asf, parameterAlphaSpending, kMax, maxInformation, p, info, spendTime,
      // bind f into std::function (capture the f we already have)
      std::function<std::vector<double>(const size_t)>(f), repp_mat);

  // Run the parallel loop over iterations
  RcppParallel::parallelFor(0, B, worker, 1 /*grain size*/);

  return worker.repp_out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix
repeatedPValueRcpp(const int kMax, const std::string &typeAlphaSpending,
                   const double parameterAlphaSpending,
                   const double maxInformation, const Rcpp::NumericMatrix &p,
                   const Rcpp::NumericMatrix &information,
                   const Rcpp::NumericMatrix &spendingTime) {

  auto p1 = flatmatrix_from_Rmatrix(p);
  auto information1 = flatmatrix_from_Rmatrix(information);
  auto spendingTime1 = flatmatrix_from_Rmatrix(spendingTime);

  auto repp = repeatedPValuecpp(static_cast<size_t>(kMax), typeAlphaSpending,
                                parameterAlphaSpending, maxInformation, p1,
                                information1, spendingTime1);
  return Rcpp::wrap(repp);
}

// Compute the first rejection study-look for each hypothesis under a
// sequential Bonferroni graphical procedure.
//
// High-level flow:
// 1) Validate/standardize inputs and expand per-hypothesis settings.
// 2) For each simulation replicate b, build per-hypothesis test trajectories
//    up to interim look k1 (p-values, information rates, spending times).
// 3) At each study look, iteratively search for rejectable hypotheses,
//    update graph weights/transitions after each rejection, and continue
//    until no further rejection is possible at that look.
//
// Return value is a B x m integer matrix where each entry is the 1-based
// study look of first rejection (0 if not rejected by look k1).
IntMatrix fseqboncpp(const std::vector<double> &w, const FlatMatrix &G,
                     const double alpha, const size_t kMax,
                     const std::vector<std::string> &typeAlphaSpending,
                     const std::vector<double> &parameterAlphaSpending,
                     const std::vector<double> &maxInformation,
                     const BoolMatrix &incidenceMatrix, const size_t k1,
                     const FlatMatrix &p, const FlatMatrix &information,
                     const FlatMatrix &spendingTime, const bool lookback) {

  size_t m = w.size();

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
  for (size_t i = 0; i < m * m; ++i) {
    if (G.data[i] < 0.0) {
      throw std::invalid_argument("G must be nonnegative");
    }
  }

  // Validation: Row sums of G must be <= 1
  std::vector<double> rowsum(m, 0.0);
  for (size_t j = 0; j < m; ++j) {
    for (size_t i = 0; i < m; ++i) {
      rowsum[i] += G(i, j);
    }
  }
  for (size_t i = 0; i < m; ++i) {
    if (rowsum[i] > 1.0 + 1e-8) {
      throw std::invalid_argument(
          "Row sums of G must be less than or equal to 1");
    }
  }

  // Validation: Diagonal elements of G must be 0
  for (size_t i = 0; i < m; ++i) {
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

  if (asf.size() == 1)
    asf.resize(m, asf[0]);
  if (asf.size() != m) {
    throw std::invalid_argument("Invalid length for typeAlphaSpending");
  }

  if (asfpar.size() == 1)
    asfpar.resize(m, asfpar[0]);
  if (asfpar.size() != m) {
    throw std::invalid_argument("Invalid length for parameterAlphaSpending");
  }

  for (size_t i = 0; i < m; ++i) {
    std::string &asfi = asf[i];
    for (char &c : asfi) {
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (!(asfi == "sfof" || asfi == "sfp" || asfi == "sfkd" ||
          asfi == "sfhsd" || asfi == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpending");
    }

    if ((asfi == "sfkd" || asfi == "sfhsd") && std::isnan(asfpar[i])) {
      throw std::invalid_argument("Missing value for parameterAlphaSpending");
    }

    if (asfi == "sfkd" && asfpar[i] <= 0) {
      throw std::invalid_argument(
          "parameterAlphaSpending must be positive for sfKD");
    }
  }

  // Validation: maxInformation must be positive for each hypothesis
  if (maxInformation.size() != m) {
    throw std::invalid_argument("Invalid length for maxInformation");
  }
  if (std::any_of(maxInformation.begin(), maxInformation.end(),
                  [](double val) { return val <= 0.0; })) {
    throw std::invalid_argument("maxInformation must be positive");
  }

  // Validation: incidenceMatrix dimension
  if (incidenceMatrix.ncol != m) {
    throw std::invalid_argument(
        "Invalid number of columns for incidenceMatrix");
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

  size_t B = p.nrow / k1;

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
    const double *src_ptr = information.data_ptr();
    double *dst_ptr = info.data_ptr();
    for (size_t j = 0; j < m; ++j) {
      for (size_t b = 0; b < B; ++b) {
        std::memcpy(dst_ptr + (b * k1) + j * (B * k1), src_ptr + j * k1,
                    k1 * sizeof(double));
      }
    }
  } else {
    info = information;
  }

  // Build spending-time matrix aligned to B * k1 rows.
  // If spendingTime == matrix(0,1,1), treat as missing and derive spending
  // from information rates later.
  FlatMatrix spendTime(B * k1, m);
  if (spendingTime.nrow == 1 && spendingTime.ncol == 1 &&
      spendingTime(0, 0) == 0) {
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
      const double *src_ptr = spendingTime.data_ptr();
      double *dst_ptr = spendTime.data_ptr();
      for (size_t j = 0; j < m; ++j) {
        for (size_t b = 0; b < B; ++b) {
          std::memcpy(dst_ptr + (b * k1) + j * (B * k1), src_ptr + j * k1,
                      k1 * sizeof(double));
        }
      }
    } else {
      spendTime = spendingTime;
    }
  }

  std::vector<size_t> K0(m, 0); // max number of testable looks for hypothesis j
  for (size_t j = 0; j < m; ++j) {
    for (size_t k = 0; k < kMax; ++k) {
      if (incidenceMatrix(k, j))
        K0[j]++;
    }
  }

  // Precompute study-look <-> testable-look mappings at interim k1.
  // idx1(l, j): study look index (0-based) for hypothesis j at testable look l.
  // idx2(k, j): testable look index (0-based) for hypothesis j at study look k,
  //             or -1 when hypothesis j is not testable at study look k.
  IntMatrix idx1(k1, m);     // study look index for each testable look
  IntMatrix idx2(k1, m);     // testable look index for each study look
  idx1.fill(-1);             // initialize with -1 for non-existent test look
  idx2.fill(-1);             // initialize with -1 for non-tested study look
  std::vector<size_t> K1(m); // last study look for each hypothesis at interim
  std::vector<size_t> K2(m); // testable looks for each hypothesis at interim
  for (size_t j = 0; j < m; ++j) {
    size_t l = 0; // index for testable looks of hypothesis j
    for (size_t k = 0; k < k1; ++k) {
      if (incidenceMatrix(k, j)) {
        idx1(l, j) = k; // study look index for the l-th testable look of hyp j
        idx2(k, j) = l; // testable look index for the k-th study look of hyp j
        ++l;
      }
    }
    K2[j] = l; // number of testable looks for hypothesis j (1-based)
    if (l > 0) {
      K1[j] = idx1(l - 1, j) + 1; // last study look for hypothesis j (1-based)
    }
  }

  // Per-replicate rejection engine used by the parallel worker.
  auto f = [&](const size_t b) -> std::vector<int> {
    std::vector<size_t> L(m);

    FlatMatrix p1(k1, m), t1(k1, m), s1(k1, m);
    p1.fill(NaN);
    t1.fill(NaN);
    s1.fill(NaN);
    double *p1_ptr = p1.data_ptr();
    double *t1_ptr = t1.data_ptr();
    double *s1_ptr = s1.data_ptr();

    // Extract p-values, information, and spending time vectors for each
    // hypothesis
    std::vector<double> p_vec;
    p_vec.reserve(k1);
    std::vector<double> i_vec;
    i_vec.reserve(k1);
    std::vector<double> s_vec;
    s_vec.reserve(k1);
    for (size_t j = 0; j < m; ++j) { // loop over hypotheses
      size_t Kj = K2[j]; // number of testable looks for hypothesis j at interim
      if (Kj == 0)
        continue; // no testable look, skip to next hypothesis
      double maxinfoj = maxInformation[j]; // maxInformation for hypothesis j

      p_vec.resize(Kj);
      i_vec.resize(Kj);
      s_vec.resize(Kj);
      for (size_t l = 0; l < Kj; ++l) {
        size_t k = idx1(l, j); // study look index for the l-th testable look

        // row index in p, info, and spendTime for iteration b and study look k
        size_t h = b * k1 + k;
        if (std::isnan(p(h, j))) {
          throw std::invalid_argument(
              "p must be provided at each testable look");
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
                      [](double v) { return v < 0.0 || v > 1.0; })) {
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
                                [](double v) { return std::isnan(v); });
      if (!all_na) {
        if (std::any_of(s_vec.begin(), s_vec.end(),
                        [](double v) { return std::isnan(v); })) {
          throw std::invalid_argument(
              "spendingTime must be provided at each testable look");
        }
        if (s_vec[0] <= 0.0) {
          throw std::invalid_argument("spendingTime must be positive");
        }
        if (any_nonincreasing(s_vec)) {
          throw std::invalid_argument(
              "spendingTime must be increasing over time");
        }
        if (s_vec[Kj - 1] > 1.0) {
          throw std::invalid_argument(
              "spendingTime must be less than or equal to 1");
        }
      }

      // Determine L[j]: number of looks to consider for hypothesis j at interim
      if (all_na) { // use information rates, no need to check s_vec further
        auto it = std::lower_bound(i_vec.begin(), i_vec.end(), maxinfoj);
        if (it == i_vec.end()) { // none >= maxInformation
          L[j] = Kj;
        } else { // first >= max info, consider looks up to & including this one
          L[j] = static_cast<size_t>(std::distance(i_vec.begin(), it)) + 1;
        }
      } else { // will use spending time, consider all testable looks
        L[j] = Kj;
      }
      size_t Lj = L[j];

      // Copy p values in one block (p_vec contiguous)
      double *p1_col = p1_ptr + j * k1; // column j begins at offset j * k1
      double *t1_col = t1_ptr + j * k1;
      double *s1_col = s1_ptr + j * k1;

      std::memcpy(p1_col, p_vec.data(), Lj * sizeof(double));

      // Information time for forming covariance matrix of test statistics
      double info_L = i_vec[Lj - 1];
      for (size_t l = 0; l < Lj; ++l) {
        t1_col[l] = i_vec[l] / info_L;
      }

      // Spending time for error spending
      if (all_na) { // use information rates
        for (size_t l = 0; l < Lj; ++l) {
          if (l == K0[j] - 1 || i_vec[l] >= maxinfoj) {
            s1_col[l] = 1.0; // the last testable look is at or beyond max info
          } else {
            s1_col[l] = i_vec[l] / maxinfoj;
          }
        }
      } else { // use spending time
        std::memcpy(s1_col, s_vec.data(), Lj * sizeof(double));
      }
    }

    size_t num_rejected = 0;       // number of hypotheses rejected so far
    std::vector<int> reject(m, 0); // first look when the hypothesis is rejected
    std::vector<double> wx = w;    // current weights, updated in-place
    FlatMatrix g = G;         // current transition matrix, updated in-place
    FlatMatrix g1(m, m);      // temporary transition matrix for update
    std::vector<double> user; // empty userAlphaSpending for getBoundcpp

    std::vector<size_t> active(m);              // currently active hypotheses
    std::iota(active.begin(), active.end(), 0); // initialize with 0,1,...,m-1
    std::vector<int> pos(m); // position map: pos[idx] = index in active or -1
    std::iota(pos.begin(), pos.end(), 0); // initialize with 0,1,...,m-1

    std::vector<double> denom(m); // temp denom for transition matrix update
    std::vector<double> u_vec;
    u_vec.reserve(k1); // upper bound from getBoundcpp

    // temp vectors reuse across hypotheses
    // Cache full t/s columns once per replicate. They are reused many times
    // in the inner rejection search and do not change during replicate b.
    std::vector<std::vector<double>> t_cols(m);
    std::vector<std::vector<double>> s_cols(m);
    for (size_t j = 0; j < m; ++j) {
      t_cols[j] = flatmatrix_get_column(t1, j);
      s_cols[j] = flatmatrix_get_column(s1, j);
    }

    std::vector<unsigned char> x(k1, 1); // efficacyStopping for getBoundcpp
    std::vector<double> w_pre(m);        // previous weights by hypothesis index

    // cached per-hypothesis previous upper bounds (u_pre[j] for hypothesis j)
    FlatMatrix u_pre(k1, m);
    double *u_pre_ptr = u_pre.data_ptr();

    // Preallocated helper vectors reused for "resuse" branch
    std::vector<double> l_vec(k1, -8.0);
    std::vector<double> theta_vec(k1, 0.0);

    ExitProbResult probs; // to store results from exitprobcpp

    size_t K3 = *std::max_element(K1.begin(), K1.end());

    for (size_t step = 0; step < K3; ++step) { // loop over study look

      // Try to find all hypotheses that can be rejected at this step
      for ([[maybe_unused]] size_t i : active) {
        bool found_reject = false;
        int found_j = -1;
        size_t step_j = 0; // study look (1-based) for rejection of hypothesis j

        if (!lookback) {
          // no lookback, only check hypotheses with testable look at current
          // step scan all hypotheses j to find a rejectable one
          for (size_t j : active) {
            if (wx[j] < 1e-8)
              continue;                // weight too small or already rejected
            int l_int = idx2(step, j); // testable look index (0-based) for
            // hypothesis j at current study look
            if (l_int < 0)
              continue; // not testable at this study look
            size_t l = static_cast<size_t>(l_int);
            if (l >= L[j])
              continue;       // beyond L[j] considered at interim
            size_t n = l + 1; // testable looks for hyp j at this study look

            double alpha1 = wx[j] * alpha;
            const std::string &asf1 = asf[j];
            double asfpar1 = asfpar[j];

            // Compute upper bound
            const std::vector<double> &t = t_cols[j];
            const std::vector<double> &s = s_cols[j];
            if (wx[j] != w_pre[j]) {
              // weights changed, compute full u_vec
              u_vec = getBoundcpp(n, t, alpha1, asf1, asfpar1, user, s, x);
            } else {
              // reuse previous u_pre[j] prefix, only solve for last element
              u_vec.resize(n);
              if (l > 0)
                std::memcpy(u_vec.data(), u_pre_ptr + j * k1,
                            l * sizeof(double));
              // compute cumulative alpha for this n
              double cumAlpha = errorSpentcpp(s[l], alpha1, asf1, asfpar1);

              // small lambda that only sets last element
              auto g = [&](double aval) -> double {
                u_vec[l] = aval;
                probs = exitprobcpp(u_vec, l_vec, theta_vec, t);
                double cpu = std::accumulate(probs.exitProbUpper.begin(),
                                             probs.exitProbUpper.end(), 0.0);
                return cpu - cumAlpha;
              };

              double g_8 = g(8.0);
              if (g_8 > 0.0) { // no alpha spent at current visit
                u_vec[l] = 8.0;
              } else {
                auto g_for_brent = [&](double aval) -> double {
                  if (aval == 8.0)
                    return g_8; // avoid recomputation at 8.0
                  return g(aval);
                };
                u_vec[l] = brent(g_for_brent, -5.0, 8.0, 1e-6);
              }
            }

            // cache computed u_vec for hypothesis j
            std::memcpy(u_pre_ptr + j * k1, u_vec.data(), n * sizeof(double));

            double alphastar = 1.0 - boost_pnorm(u_vec[l]);
            if (p1(l, j) <= alphastar) {
              found_reject = true;
              found_j = j;
              step_j = step + 1; // current study look (1-based)
              break;             // stop scanning j, we'll process rejection
            }
          } // end scan j
        } else {
          // lookback, need to check all active hypotheses at each step
          // scan all hypotheses j to find a rejectable one
          for (size_t j : active) {
            if (wx[j] < 1e-8)
              continue; // weight too small or already rejected

            // find the number of testable looks for hypothesis j at this study
            // look
            int l_int = -1;
            for (size_t ll = 0; ll <= step; ++ll) {
              int l_int_ll = idx2(ll, j); // testable look at study look ll
              if (l_int_ll >= 0) {
                // update the last testable look up to current step
                l_int = l_int_ll;
              }
            }

            if (l_int < 0)
              continue; // not testable up to this study look
            size_t l = static_cast<size_t>(l_int);
            if (l >= L[j])
              continue;       // beyond L[j] considered at interim
            size_t n = l + 1; // testable looks for hyp j up to this study look

            double alpha1 = wx[j] * alpha;
            const std::string &asf1 = asf[j];
            double asfpar1 = asfpar[j];

            // Compute upper bound
            const std::vector<double> &t = t_cols[j];
            const std::vector<double> &s = s_cols[j];
            if (wx[j] != w_pre[j]) {
              // weights changed, compute full u_vec
              u_vec = getBoundcpp(n, t, alpha1, asf1, asfpar1, user, s, x);
              // cache computed u_vec for hypothesis j
              std::memcpy(u_pre_ptr + j * k1, u_vec.data(), n * sizeof(double));
            } else if (idx2(step, j) >= 0) {
              // weights unchanged and testable at current study look,
              // reuse previous u_pre[j] prefix, only solve for last element
              u_vec.resize(n);
              if (l > 0)
                std::memcpy(u_vec.data(), u_pre_ptr + j * k1,
                            l * sizeof(double));
              // compute cumulative alpha for this n
              double cumAlpha = errorSpentcpp(s[l], alpha1, asf1, asfpar1);

              // small lambda that only sets last element
              auto g = [&](double aval) -> double {
                u_vec[l] = aval;
                probs = exitprobcpp(u_vec, l_vec, theta_vec, t);
                double cpu = std::accumulate(probs.exitProbUpper.begin(),
                                             probs.exitProbUpper.end(), 0.0);
                return cpu - cumAlpha;
              };

              double g_8 = g(8.0);
              if (g_8 > 0.0) { // no alpha spent at current visit
                u_vec[l] = 8.0;
              } else {
                auto g_for_brent = [&](double aval) -> double {
                  if (aval == 8.0)
                    return g_8; // avoid recomputation at 8.0
                  return g(aval);
                };
                u_vec[l] = brent(g_for_brent, -5.0, 8.0, 1e-6);
              }
              // cache computed u_vec for hypothesis j
              std::memcpy(u_pre_ptr + j * k1, u_vec.data(), n * sizeof(double));
            } else {
              // no weight change and not testable at current study look,
              // no need to test for rejection at current look
              continue;
            }

            // test rejection
            if (wx[j] == w_pre[j]) {
              // weight unchanged since last check, only need to check current
              // look l
              double alphastar = 1.0 - boost_pnorm(u_vec[l]);
              if (p1(l, j) <= alphastar) {
                found_reject = true;
                found_j = j;
                step_j = step + 1; // current study look (1-based)
                break;             // stop scanning j, we'll process rejection
              }
            } else {
              // weight changed, need to check all previous looks up to l
              bool reject_j = false;
              for (size_t ll = 0; ll <= l; ++ll) {
                double alphastar = 1.0 - boost_pnorm(u_vec[ll]);
                if (p1(ll, j) <= alphastar) {
                  reject_j = true;
                  step_j = idx1(ll, j) + 1; // study look for this testable look
                  break;
                }
              }
              if (reject_j) {
                found_reject = true;
                found_j = j;
                break; // stop scanning j, we'll process rejection
              }
            }
          } // end scan j
        }

        if (!found_reject) { // no more rejections at this study look
          // remember current weights for next study look to check if weights
          // changed
          w_pre = wx;
          break;
        }

        // Process rejection of hypothesis found_j
        size_t j = found_j;
        reject[j] = step_j;
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
        for (size_t l : active) {
          wx[l] += wxj * g(j, l);
        }
        wx[j] = 0.0; // weight of rejected hypothesis becomes 0

        // Update transition matrix for active hypotheses
        for (size_t l : active) {
          denom[l] = 1.0 - g(l, j) * g(j, l);
        }

        g1.fill(0.0); // reset g1 to 0 before filling
        for (size_t k : active) {
          double g_jk = g(j, k);
          for (size_t l : active) {
            if (l == k)
              continue;
            double dl = denom[l];
            if (dl > 1e-12) {
              g1(l, k) = (g(l, k) + g(l, j) * g_jk) / dl;
            }
          }
        }

        std::swap(g.data, g1.data); // copy g1 to g for next iteration

        // Stop if all hypotheses rejected
        if (num_rejected == m)
          break;
      } // end attempt loop

      if (num_rejected == m)
        break;
    } // end study look loop

    return reject;
  };

  struct SimulationWorker : public RcppParallel::Worker {
    // references to read-only inputs (no mutation)
    const size_t m;
    const size_t k1;
    const std::vector<double> &w;
    const FlatMatrix &G;
    const double alpha;
    const std::vector<std::string> &asf;
    const std::vector<double> &asfpar;
    const std::vector<size_t> &K0;
    const std::vector<size_t> &K1;
    const std::vector<size_t> &K2;
    const IntMatrix &idx1;
    const IntMatrix &idx2;
    const std::vector<double> &maxInformation;
    const FlatMatrix &p;
    const FlatMatrix &info;
    const FlatMatrix &spendTime;
    const bool lookback;

    // function f and other params that f needs are captured from outer scope
    // capture them by reference here so worker can call f(...)
    std::function<std::vector<int>(const size_t)> f;

    // result references (each iteration writes unique index into these)
    IntMatrix &reject_out; // B by m matrix to store rejection results

    // constructor
    SimulationWorker(const size_t m_, const size_t k1_,
                     const std::vector<double> &w_, const FlatMatrix &G_,
                     const double alpha_, const std::vector<std::string> &asf_,
                     const std::vector<double> &asfpar_,
                     const std::vector<size_t> &K0_,
                     const std::vector<size_t> &K1_,
                     const std::vector<size_t> &K2_, const IntMatrix &idx1_,
                     const IntMatrix &idx2_,
                     const std::vector<double> &maxInformation_,
                     const FlatMatrix &p_, const FlatMatrix &info_,
                     const FlatMatrix &spendTime_, const bool lookback_,
                     decltype(f) f_, IntMatrix &reject_out_)
        :

          m(m_), k1(k1_), w(w_), G(G_), alpha(alpha_), asf(asf_),
          asfpar(asfpar_), K0(K0_), K1(K1_), K2(K2_), idx1(idx1_), idx2(idx2_),
          maxInformation(maxInformation_), p(p_), info(info_),
          spendTime(spendTime_), lookback(lookback_), f(std::move(f_)),
          reject_out(reject_out_) {}

    // operator() processes a range of bootstrap iterations [begin, end)
    void operator()(size_t begin, size_t end) {
      for (size_t b = begin; b < end; ++b) {
        // call the (thread-safe) per-iteration function f
        std::vector<int> out = f(b);

        // write results
        for (size_t j = 0; j < m; ++j) {
          reject_out(b, j) = out[j];
        }
      } // end for b
    } // end operator()
  }; // end BootstrapWorker

  IntMatrix reject_mat(B, m); // B by m matrix to store rejection results

  // Fast path: avoid parallel scheduling overhead for a single replicate.
  if (B == 1) {
    std::vector<int> out = f(0);
    for (size_t j = 0; j < m; ++j) {
      reject_mat(0, j) = out[j];
    }
    return reject_mat;
  }

  // Instantiate the Worker with references to inputs and outputs
  SimulationWorker worker(
      m, k1, w, G, alpha, asf, asfpar, K0, K1, K2, idx1, idx2, maxInformation,
      p, info, spendTime, lookback,
      // bind f into std::function (capture the f we already have)
      std::function<std::vector<int>(const size_t)>(f), reject_mat);

  // Run the parallel loop over iterations
  RcppParallel::parallelFor(0, B, worker, 1 /*grain size*/);

  return worker.reject_out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix
fseqbonRcpp(const Rcpp::NumericVector &w, const Rcpp::NumericMatrix &G,
            const double alpha, const int kMax,
            const Rcpp::StringVector &typeAlphaSpending,
            const Rcpp::NumericVector &parameterAlphaSpending,
            const Rcpp::NumericVector &maxInformation,
            const Rcpp::LogicalMatrix &incidenceMatrix, const int k1,
            const Rcpp::NumericMatrix &p,
            const Rcpp::NumericMatrix &information,
            const Rcpp::NumericMatrix &spendingTime, const bool lookback) {

  auto w1 = Rcpp::as<std::vector<double>>(w);
  auto G1 = flatmatrix_from_Rmatrix(G);
  auto asf1 = Rcpp::as<std::vector<std::string>>(typeAlphaSpending);
  auto asfpar1 = Rcpp::as<std::vector<double>>(parameterAlphaSpending);
  auto maxInfo1 = Rcpp::as<std::vector<double>>(maxInformation);
  auto incid1 = boolmatrix_from_Rmatrix(incidenceMatrix);
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto info1 = flatmatrix_from_Rmatrix(information);
  auto spendTime1 = flatmatrix_from_Rmatrix(spendingTime);
  auto reject1 = fseqboncpp(w1, G1, alpha, static_cast<size_t>(kMax), asf1,
                            asfpar1, maxInfo1, incid1, static_cast<size_t>(k1),
                            p1, info1, spendTime1, lookback);
  return Rcpp::wrap(reject1);
}

// Converts step-down p-values to sequential adjusted p-values
FlatMatrix fstp2seqcpp(const FlatMatrix &p, const std::vector<double> &gamma,
                       const std::string &test, const bool retest) {

  // Validate dimensions
  size_t nreps = p.nrow;
  size_t n = p.ncol;
  size_t m = n / 2;

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
                  [](double v) { return v < 0.0 || v > 1.0; })) {
    throw std::invalid_argument("p-values must be between 0 and 1");
  }

  bool is_holm = (test1 == "holm");

  // Step 1: Compute matrix a (nreps × m) from pairs of p-values
  FlatMatrix a(nreps, m);
  for (size_t j = 0; j < m; ++j) {
    double gamma_j = gamma[j];
    double inv_factor = 1.0 / (1.0 + gamma_j);

    for (size_t iter = 0; iter < nreps; ++iter) {
      double x1 = p(iter, 2 * j);
      double x2 = p(iter, 2 * j + 1);
      double val = 2.0 * std::max(x1, x2) * inv_factor;
      if (is_holm) {
        val = std::max(val, 2.0 * std::min(x1, x2));
      }
      a(iter, j) = val;
    }
  }

  // Step 2: Compute matrix d (m × m) - lower triangular discount factors
  FlatMatrix d(m, m);
  for (size_t s = 0; s < m; ++s) {
    double gmax = 0.0;
    for (size_t j = s; j < m; ++j) {
      if (j > s) {
        gmax = std::max(gmax, gamma[j - 1]);
      }
      d(s, j) = std::max((1.0 - gmax) * 0.5, 1e-12);
    }
  }

  // Step 3: Compute adjusted p-values
  FlatMatrix padj(nreps, n);

  // Preallocate reusable vectors
  std::vector<double> a_row(m);    // current row of a
  std::vector<double> a_cummax(m); // cumulative max
  for (size_t iter = 0; iter < nreps; ++iter) {
    // Extract row iter from a and compute cumulative max
    for (size_t j = 0; j < m; ++j) {
      a_row[j] = a(iter, j);
    }

    a_cummax[0] = a_row[0];
    for (size_t j = 1; j < m; ++j) {
      a_cummax[j] = std::max(a_cummax[j - 1], a_row[j]);
    }

    // For each hypothesis index i
    for (size_t i = 0; i < m; ++i) {
      double t1 = a_cummax[i]; // max of a[0..i]

      // --- Compute padj(iter, 2*i) (even index) ---
      double t2_even = 1.0;
      for (size_t s = 0; s <= i; ++s) {
        double lhs = (s > 0) ? a_cummax[s - 1] : 0.0;

        // max over j ∈ [s, i] of p(iter, 2*j) / d(s, j)
        for (size_t j = s; j <= i; ++j) {
          lhs = std::max(lhs, p(iter, 2 * j) / d(s, j));
        }

        double rhs = 2.0 * p(iter, 2 * s + 1) / (1.0 + gamma[s]);
        if (lhs < rhs) {
          t2_even = std::min(t2_even, lhs);
        }
      }

      double result_even;
      if (retest && m > 1) {
        double t3_even = 1.0;
        size_t s_max = std::min(i, m - 2);

        for (size_t s = 0; s <= s_max; ++s) {
          double lhs = (s > 0) ? a_cummax[s - 1] : 0.0;

          // max over j ∈ [s, m-1] of p(iter, 2*j+1) / d(s, j)
          for (size_t j = s; j < m; ++j) {
            lhs = std::max(lhs, p(iter, 2 * j + 1) / d(s, j));
          }

          // max over j ∈ [s, i] of p(iter, 2*j)
          for (size_t j = s; j <= i; ++j) {
            lhs = std::max(lhs, p(iter, 2 * j));
          }

          double rhs = 2.0 * p(iter, 2 * s) / (1.0 + gamma[s]);
          if (lhs < rhs) {
            t3_even = std::min(t3_even, lhs);
          }
        }

        result_even = std::min({t1, t2_even, t3_even});
      } else {
        result_even = std::min(t1, t2_even);
      }

      padj(iter, 2 * i) = result_even;

      // --- Compute padj(iter, 2*i+1) (odd index) ---
      double t2_odd = 1.0;
      for (size_t s = 0; s <= i; ++s) {
        double lhs = (s > 0) ? a_cummax[s - 1] : 0.0;

        // max over j ∈ [s, i] of p(iter, 2*j+1) / d(s, j)
        for (size_t j = s; j <= i; ++j) {
          lhs = std::max(lhs, p(iter, 2 * j + 1) / d(s, j));
        }

        double rhs = 2.0 * p(iter, 2 * s) / (1.0 + gamma[s]);
        if (lhs < rhs) {
          t2_odd = std::min(t2_odd, lhs);
        }
      }

      double result_odd;
      if (retest && m > 1) {
        double t3_odd = 1.0;
        size_t s_max = std::min(i, m - 2);

        for (size_t s = 0; s <= s_max; ++s) {
          double lhs = (s > 0) ? a_cummax[s - 1] : 0.0;

          // max over j ∈ [s, m-1] of p(iter, 2*j) / d(s, j)
          for (size_t j = s; j < m; ++j) {
            lhs = std::max(lhs, p(iter, 2 * j) / d(s, j));
          }

          // max over j ∈ [s, i] of p(iter, 2*j+1)
          for (size_t j = s; j <= i; ++j) {
            lhs = std::max(lhs, p(iter, 2 * j + 1));
          }

          double rhs = 2.0 * p(iter, 2 * s + 1) / (1.0 + gamma[s]);
          if (lhs < rhs) {
            t3_odd = std::min(t3_odd, lhs);
          }
        }

        result_odd = std::min({t1, t2_odd, t3_odd});
      } else {
        result_odd = std::min(t1, t2_odd);
      }

      padj(iter, 2 * i + 1) = result_odd;
    }
  }

  return padj;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix fstp2seqRcpp(const Rcpp::NumericMatrix &p,
                                 const Rcpp::NumericVector &gamma,
                                 const std::string &test = "hochberg",
                                 const bool retest = true) {

  auto p1 = flatmatrix_from_Rmatrix(p);
  auto gamma1 = Rcpp::as<std::vector<double>>(gamma);
  auto padj1 = fstp2seqcpp(p1, gamma1, test, retest);
  return Rcpp::wrap(padj1);
}

// Helper: compute row sums for BoolMatrix
std::vector<size_t> boolmatrix_rowsums(const BoolMatrix &mat) {
  size_t nrow = mat.nrow;
  size_t ncol = mat.ncol;
  std::vector<size_t> sums(nrow, 0);
  const unsigned char *data = mat.data_ptr();

  for (size_t j = 0; j < ncol; ++j) {
    const unsigned char *col = data + j * nrow;
    for (size_t i = 0; i < nrow; ++i) {
      sums[i] += col[i];
    }
  }
  return sums;
}

// Helper for adjusted p-values for standard mixture gatekeeping procedures
AdjustedPValues fstdmixcpp(const FlatMatrix &p, const BoolMatrix &family,
                           const BoolMatrix &serial, const BoolMatrix &parallel,
                           const std::vector<double> &gamma,
                           const std::string &test, const bool exhaust) {

  // Validate inputs
  size_t nreps = p.nrow;
  size_t m = p.ncol;

  if (family.ncol != m) {
    throw std::invalid_argument("family must have m columns");
  }

  if (serial.nrow != m || serial.ncol != m) {
    throw std::invalid_argument("serial must be m x m matrix");
  }

  if (parallel.nrow != m || parallel.ncol != m) {
    throw std::invalid_argument("parallel must be m x m matrix");
  }

  size_t nfamily = family.nrow;

  if (gamma.size() != static_cast<size_t>(nfamily)) {
    throw std::invalid_argument("gamma length must equal nrow(family)");
  }

  // Validate gamma in [0, 1]
  if (std::any_of(gamma.begin(), gamma.end(),
                  [](double g) { return g < 0.0 || g > 1.0; })) {
    throw std::invalid_argument("gamma values must be between 0 and 1");
  }

  // Validate p-values in [0, 1]
  if (std::any_of(p.data.begin(), p.data.end(),
                  [](double v) { return v < 0.0 || v > 1.0; })) {
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
  std::vector<size_t> nhyps = boolmatrix_rowsums(family);

  // Number of intersection tests
  size_t ntests = (1 << m) - 1; // 2^m - 1

  // Matrix to store local p-values for intersection tests
  FlatMatrix pinter(nreps, ntests);

  // Incidence matrix for elementary hypotheses
  IntMatrix incid(ntests, m);

  // Process each intersection hypothesis
  for (size_t i = 0; i < ntests; ++i) {
    size_t number = ntests - i;

    // Binary representation of intersection hypothesis
    std::vector<int> cc(m);
    for (size_t j = 0; j < m; ++j) {
      cc[j] = (number >> (m - 1 - j)) & 1;
    }

    // Active hypotheses in each family, I1, ..., I_nfamily
    BoolMatrix family0(nfamily, m);
    for (size_t j = 0; j < m; ++j) {
      for (size_t h = 0; h < nfamily; ++h) {
        family0(h, j) = family(h, j) & cc[j];
      }
    }

    // Number of active hypotheses by family, k1, ..., k_nfamily
    std::vector<size_t> nhyps0 = boolmatrix_rowsums(family0);

    // Determine restricted index set for each family
    std::vector<int> cc1(m, 1);

    for (size_t j = 0; j < m; ++j) {
      // Check serial constraints
      int serial_sum = 0;
      for (size_t k = 0; k < m; ++k) {
        serial_sum += serial(j, k);
      }

      if (serial_sum > 0) {
        // if any serial predecessor is accepted, remove j
        for (size_t k = 0; k < m; ++k) {
          if (serial(j, k) && cc[k]) {
            cc1[j] = 0;
            break;
          }
        }

        // if any serial predecessor is not testable, remove j
        for (size_t k = 0; k < m; ++k) {
          if (serial(j, k) && !cc1[k]) {
            cc1[j] = 0;
            break;
          }
        }
      }

      // Check parallel constraints
      int parallel_sum = 0;
      for (size_t k = 0; k < m; ++k) {
        parallel_sum += parallel(j, k);
      }

      if (parallel_sum > 0) {
        // if none of the parallel predecessors are rejected, remove j
        bool hit = true;
        for (size_t k = 0; k < m; ++k) {
          if (parallel(j, k) && !cc[k]) {
            hit = false;
            break;
          }
        }
        if (hit)
          cc1[j] = 0;

        // if none of the parallel predecessors are testable, remove j
        hit = true;
        for (size_t k = 0; k < m; ++k) {
          if (parallel(j, k) && cc1[k]) {
            hit = false;
            break;
          }
        }
        if (hit)
          cc1[j] = 0;
      }
    }

    // Apply intersection
    for (size_t j = 0; j < m; ++j) {
      cc1[j] = cc1[j] & cc[j];
    }

    // Error rate function divided by alpha
    std::vector<double> errf(nfamily);
    for (size_t j = 0; j < nfamily; ++j) {
      errf[j] = nhyps0[j] > 0
                    ? gamma[j] + (1.0 - gamma[j]) * nhyps0[j] / nhyps[j]
                    : 0.0;
    }

    // Allocated fraction of alpha for each family
    std::vector<double> coef(nfamily);
    coef[0] = 1.0;
    for (size_t j = 1; j < nfamily; ++j) {
      coef[j] = coef[j - 1] * (1.0 - errf[j - 1]);
    }

    size_t kmax = 0; // last family with positive allocation (1-based index)
    for (size_t j = nfamily; j-- > 0;) {
      if (coef[j] > 0.0) {
        kmax = j + 1;
        break;
      }
    }

    // Families up to kmax, I_1*, ..., I_kmax*
    BoolMatrix family1(kmax, m);
    for (size_t j = 0; j < kmax; ++j) {
      for (size_t k = 0; k < m; ++k) {
        family1(j, k) = family(j, k) & cc1[k];
      }
    }

    // Number of testable hypotheses by family, k_1*, ..., k_kmax*
    std::vector<size_t> nhyps1 = boolmatrix_rowsums(family1);

    // Indices of active families
    std::vector<size_t> sub;
    sub.reserve(kmax);
    for (size_t j = 0; j < kmax; ++j) {
      if (nhyps1[j] > 0) {
        sub.push_back(j);
      }
    }

    size_t nfamily2 = sub.size();

    // Subset of active families after removing those without testable
    // hypotheses
    BoolMatrix family2(nfamily2, m);
    std::vector<size_t> nhyps2(nfamily2);
    for (size_t j = 0; j < nfamily2; ++j) {
      for (size_t k = 0; k < m; ++k) {
        family2(j, k) = family1(sub[j], k);
      }
      nhyps2[j] = nhyps1[sub[j]];
    }

    // family indices and hypothesis indices for testable hypotheses
    std::vector<size_t> fam, hyps2;
    for (size_t j = 0; j < nfamily2; ++j) {
      for (size_t k = 0; k < m; ++k) {
        if (family2(j, k)) {
          fam.push_back(j);
          hyps2.push_back(k);
        }
      }
    }

    // total number of testable hypotheses across active families
    size_t n = hyps2.size();

    // Relative importance for active families
    std::vector<double> coef1(nfamily2);
    for (size_t j = 0; j < nfamily2; ++j) {
      coef1[j] = coef[sub[j]];
    }

    // Broadcasted family weight for each testable hypothesis
    std::vector<double> c(n);
    for (size_t k = 0; k < n; ++k) {
      c[k] = coef1[fam[k]];
    }

    // Truncation parameters
    std::vector<double> gam2(nfamily2);
    for (size_t j = 0; j < nfamily2; ++j) {
      gam2[j] = gamma[sub[j]];
    }
    if (exhaust)
      gam2[nfamily2 - 1] = 1.0;

    // Bonferroni part of weights
    std::vector<double> coef2(nfamily2);
    for (size_t j = 0; j < nfamily2; ++j) {
      coef2[j] = (1.0 - gam2[j]) / nhyps[sub[j]];
    }

    // Broadcasted Bonferroni part of weights to each testable hypothesis
    std::vector<double> tbon(n);
    for (size_t k = 0; k < n; ++k) {
      tbon[k] = coef2[fam[k]];
    }

    // Cumulative count of hypotheses before the current family
    std::vector<size_t> ck(nfamily2 + 1, 0);
    for (size_t j = 1; j <= nfamily2; ++j) {
      ck[j] = ck[j - 1] + nhyps2[j - 1];
    }

    // Compute weights
    std::vector<double> w(n);
    if (test1 == "hommel") {
      for (size_t k = 0; k < n; ++k) {
        size_t l = fam[k];
        size_t j = (k + 1) - ck[l];
        w[k] = j * gam2[l] / nhyps2[l] + tbon[k];
      }
    } else if (test1 == "hochberg") {
      for (size_t k = 0; k < n; ++k) {
        size_t l = fam[k];
        size_t j = (k + 1) - ck[l];
        w[k] = gam2[l] / (nhyps2[l] - j + 1) + tbon[k];
      }
    } else { // holm
      for (size_t k = 0; k < n; ++k) {
        size_t l = fam[k];
        w[k] = gam2[l] / nhyps2[l] + tbon[k];
      }
    }

    // Process each replication
    std::vector<double> p1(n), p1s(n), p2(n);
    for (size_t iter = 0; iter < nreps; ++iter) {
      // Extract raw p-values
      for (size_t k = 0; k < n; ++k) {
        p1[k] = p(iter, hyps2[k]);
      }

      // Sort p-values within each family
      for (size_t j = 0; j < nfamily2; ++j) {
        size_t start = ck[j];
        size_t end = ck[j + 1];
        size_t len = end - start;
        p1s = subset(p1, start, end);
        std::sort(p1s.begin(), p1s.end());
        std::memcpy(p2.data() + start, p1s.data(), len * sizeof(double));
      }

      // Compute minimum ratio
      double min_val = 1.0;
      for (size_t k = 0; k < n; ++k) {
        double ratio = p2[k] / (w[k] * c[k]);
        min_val = std::min(min_val, ratio);
      }

      pinter(iter, i) = min_val;
    }

    // Store incidence
    for (size_t j = 0; j < m; ++j) {
      incid(i, j) = cc[j];
    }
  }

  // Compute adjusted p-values for elementary hypotheses
  FlatMatrix padj(nreps, m);
  for (size_t j = 0; j < m; ++j) {
    for (size_t iter = 0; iter < nreps; ++iter) {
      double max_p = 0.0;
      for (size_t i = 0; i < ntests; ++i) {
        if (incid(i, j)) {
          max_p = std::max(max_p, pinter(iter, i));
        }
      }
      padj(iter, j) = std::min(max_p, 1.0);
    }
  }

  // Apply logical restrictions (serial constraints)
  for (size_t j = 0; j < m; ++j) {
    int serial_sum = 0;
    for (size_t k = 0; k < m; ++k) {
      serial_sum += serial(j, k);
    }

    if (serial_sum > 0) {
      for (size_t iter = 0; iter < nreps; ++iter) {
        double pre = 0.0;
        for (size_t k = 0; k < m; ++k) {
          if (serial(j, k)) {
            pre = std::max(pre, padj(iter, k));
          }
        }
        padj(iter, j) = std::max(padj(iter, j), pre);
      }
    }
  }

  // Apply logical restrictions (parallel constraints)
  for (size_t j = 0; j < m; ++j) {
    int parallel_sum = 0;
    for (size_t k = 0; k < m; ++k) {
      parallel_sum += parallel(j, k);
    }

    if (parallel_sum > 0) {
      for (size_t iter = 0; iter < nreps; ++iter) {
        double pre = 1.0;
        for (size_t k = 0; k < m; ++k) {
          if (parallel(j, k)) {
            pre = std::min(pre, padj(iter, k));
          }
        }
        padj(iter, j) = std::max(padj(iter, j), pre);
      }
    }
  }

  return AdjustedPValues{incid, pinter, padj};
}

// [[Rcpp::export]]
Rcpp::List fstdmixRcpp(const Rcpp::NumericMatrix &p,
                       const Rcpp::LogicalMatrix &family,
                       const Rcpp::LogicalMatrix &serial,
                       const Rcpp::LogicalMatrix &parallel,
                       const Rcpp::NumericVector &gamma,
                       const std::string &test = "hommel",
                       const bool exhaust = true) {
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto family1 = boolmatrix_from_Rmatrix(family);
  auto serial1 = boolmatrix_from_Rmatrix(serial);
  auto parallel1 = boolmatrix_from_Rmatrix(parallel);
  auto gamma1 = Rcpp::as<std::vector<double>>(gamma);
  auto out = fstdmixcpp(p1, family1, serial1, parallel1, gamma1, test, exhaust);
  ListCpp result;
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.pinter), "pinter");
  result.push_back(std::move(out.padj), "padj");
  return Rcpp::wrap(result);
}

// Helper for adjusted p-values for modified mixture gatekeeping procedures
AdjustedPValues fmodmixcpp(const FlatMatrix &p, const BoolMatrix &family,
                           const BoolMatrix &serial, const BoolMatrix &parallel,
                           const std::vector<double> &gamma,
                           const std::string &test, const bool exhaust) {

  // Validate inputs
  size_t nreps = p.nrow;
  size_t m = p.ncol;

  if (family.ncol != m) {
    throw std::invalid_argument("family must have m columns");
  }

  if (serial.nrow != m || serial.ncol != m) {
    throw std::invalid_argument("serial must be m x m matrix");
  }

  if (parallel.nrow != m || parallel.ncol != m) {
    throw std::invalid_argument("parallel must be m x m matrix");
  }

  size_t nfamily = family.nrow;

  if (gamma.size() != static_cast<size_t>(nfamily)) {
    throw std::invalid_argument("gamma length must equal nrow(family)");
  }

  // Validate gamma in [0, 1]
  if (std::any_of(gamma.begin(), gamma.end(),
                  [](double g) { return g < 0.0 || g > 1.0; })) {
    throw std::invalid_argument("gamma values must be between 0 and 1");
  }

  // Validate p-values in [0, 1]
  if (std::any_of(p.data.begin(), p.data.end(),
                  [](double v) { return v < 0.0 || v > 1.0; })) {
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
  std::vector<size_t> nhyps = boolmatrix_rowsums(family);

  // Number of intersection tests
  size_t ntests = (1 << m) - 1; // 2^m - 1

  // Matrix to store local p-values for intersection tests
  FlatMatrix pinter(nreps, ntests);

  // Incidence matrix for elementary hypotheses
  IntMatrix incid(ntests, m);

  // Process each intersection hypothesis
  for (size_t i = 0; i < ntests; ++i) {
    size_t number = ntests - i;

    // Binary representation of intersection hypothesis
    std::vector<int> cc(m);
    for (size_t j = 0; j < m; ++j) {
      cc[j] = (number >> (m - 1 - j)) & 1;
    }

    // Active hypotheses in each family, I1, ..., I_nfamily
    BoolMatrix family0(nfamily, m);
    for (size_t j = 0; j < m; ++j) {
      for (size_t h = 0; h < nfamily; ++h) {
        family0(h, j) = family(h, j) & cc[j];
      }
    }

    // Number of active hypotheses by family, k1, ..., k_nfamily
    std::vector<size_t> nhyps0 = boolmatrix_rowsums(family0);

    // Determine restricted index set for each family
    std::vector<int> cc1(m, 1);

    for (size_t j = 0; j < m; ++j) {
      // Check serial constraints
      int serial_sum = 0;
      for (size_t k = 0; k < m; ++k) {
        serial_sum += serial(j, k);
      }

      if (serial_sum > 0) {
        // if any serial predecessor is accepted, remove j
        for (size_t k = 0; k < m; ++k) {
          if (serial(j, k) && cc[k]) {
            cc1[j] = 0;
            break;
          }
        }

        // if any serial predecessor is not testable, remove j
        for (size_t k = 0; k < m; ++k) {
          if (serial(j, k) && !cc1[k]) {
            cc1[j] = 0;
            break;
          }
        }
      }

      // Check parallel constraints
      int parallel_sum = 0;
      for (size_t k = 0; k < m; ++k) {
        parallel_sum += parallel(j, k);
      }

      if (parallel_sum > 0) {
        // if none of the parallel predecessors are rejected, remove j
        bool hit = true;
        for (size_t k = 0; k < m; ++k) {
          if (parallel(j, k) && !cc[k]) {
            hit = false;
            break;
          }
        }
        if (hit)
          cc1[j] = 0;

        // if none of the parallel predecessors are testable, remove j
        hit = true;
        for (size_t k = 0; k < m; ++k) {
          if (parallel(j, k) && cc1[k]) {
            hit = false;
            break;
          }
        }
        if (hit)
          cc1[j] = 0;
      }
    }

    std::vector<int> cc2 = cc1; // for nstar calculation

    // Apply intersection, for kstar calculation
    for (size_t j = 0; j < m; ++j) {
      cc1[j] = cc1[j] & cc[j];
    }

    // Compute kstar and nstar for error rate function
    std::vector<size_t> kstar(nfamily), nstar(nfamily);
    std::vector<double> errf(nfamily);
    for (size_t j = 0; j < nfamily; ++j) {
      kstar[j] = 0;
      nstar[j] = 0;

      for (size_t k = 0; k < m; ++k) {
        if (family(j, k) && cc1[k])
          kstar[j]++;
        if (family(j, k) && cc2[k])
          nstar[j]++;
      }

      // DIFFERENCE from stdmix: uses kstar and nstar instead of nhyps0 and
      // nhyps
      errf[j] = kstar[j] > 0 ? gamma[j] + (1.0 - gamma[j]) * kstar[j] / nstar[j]
                             : 0.0;
    }

    // Allocated fraction of alpha for each family
    std::vector<double> coef(nfamily);
    coef[0] = 1.0;
    for (size_t j = 1; j < nfamily; ++j) {
      coef[j] = coef[j - 1] * (1.0 - errf[j - 1]);
    }

    size_t kmax = 0; // last family with positive allocation (1-based index)
    for (size_t j = nfamily; j-- > 0;) {
      if (coef[j] > 0.0) {
        kmax = j + 1;
        break;
      }
    }

    // Families up to kmax, I_1*, ..., I_kmax*
    BoolMatrix family1(kmax, m);
    for (size_t j = 0; j < kmax; ++j) {
      for (size_t k = 0; k < m; ++k) {
        family1(j, k) = family(j, k) & cc1[k];
      }
    }

    // Number of testable hypotheses by family, k_1*, ..., k_kmax*
    std::vector<size_t> nhyps1 = boolmatrix_rowsums(family1);

    // Indices of active families
    std::vector<size_t> sub;
    sub.reserve(kmax);
    for (size_t j = 0; j < kmax; ++j) {
      if (nhyps1[j] > 0) {
        sub.push_back(j);
      }
    }

    size_t nfamily2 = sub.size();

    // Subset of active families after removing those without testable
    // hypotheses
    BoolMatrix family2(nfamily2, m);
    std::vector<size_t> nhyps2(nfamily2);
    for (size_t j = 0; j < nfamily2; ++j) {
      for (size_t k = 0; k < m; ++k) {
        family2(j, k) = family1(sub[j], k);
      }
      nhyps2[j] = nhyps1[sub[j]];
    }

    // family indices and hypothesis indices for testable hypotheses
    std::vector<size_t> fam, hyps2;
    for (size_t j = 0; j < nfamily2; ++j) {
      for (size_t k = 0; k < m; ++k) {
        if (family2(j, k)) {
          fam.push_back(j);
          hyps2.push_back(k);
        }
      }
    }

    // total number of testable hypotheses across active families
    size_t n = hyps2.size();

    // Relative importance for active families
    std::vector<double> coef1(nfamily2);
    for (size_t j = 0; j < nfamily2; ++j) {
      coef1[j] = coef[sub[j]];
    }

    // Broadcasted family weight for each testable hypothesis
    std::vector<double> c(n);
    for (size_t k = 0; k < n; ++k) {
      c[k] = coef1[fam[k]];
    }

    // Truncation parameters
    std::vector<double> gam2(nfamily2);
    for (size_t j = 0; j < nfamily2; ++j) {
      gam2[j] = gamma[sub[j]];
    }
    if (exhaust)
      gam2[nfamily2 - 1] = 1.0;

    // Bonferroni part of weights
    // (KEY DIFFERENCE from stdmix: uses nstar instead of nhyps)
    std::vector<double> coef2(nfamily2);
    for (size_t j = 0; j < nfamily2; ++j) {
      coef2[j] = (1.0 - gam2[j]) / nstar[sub[j]];
    }

    // Broadcasted Bonferroni part of weights to each testable hypothesis
    std::vector<double> tbon(n);
    for (size_t k = 0; k < n; ++k) {
      tbon[k] = coef2[fam[k]];
    }

    // Cumulative count of hypotheses before the current family
    std::vector<size_t> ck(nfamily2 + 1, 0);
    for (size_t j = 1; j <= nfamily2; ++j) {
      ck[j] = ck[j - 1] + nhyps2[j - 1];
    }

    // Compute weights
    std::vector<double> w(n);
    if (test1 == "hommel") {
      for (size_t k = 0; k < n; ++k) {
        size_t l = fam[k];
        size_t j = (k + 1) - ck[l];
        w[k] = j * gam2[l] / nhyps2[l] + tbon[k];
      }
    } else if (test1 == "hochberg") {
      for (size_t k = 0; k < n; ++k) {
        size_t l = fam[k];
        size_t j = (k + 1) - ck[l];
        w[k] = gam2[l] / (nhyps2[l] - j + 1) + tbon[k];
      }
    } else { // holm
      for (size_t k = 0; k < n; ++k) {
        size_t l = fam[k];
        w[k] = gam2[l] / nhyps2[l] + tbon[k];
      }
    }

    // Process each replication
    std::vector<double> p1(n), p1s(n), p2(n);
    for (size_t iter = 0; iter < nreps; ++iter) {
      // Extract raw p-values
      for (size_t k = 0; k < n; ++k) {
        p1[k] = p(iter, hyps2[k]);
      }

      // Sort p-values within each family
      for (size_t j = 0; j < nfamily2; ++j) {
        size_t start = ck[j];
        size_t end = ck[j + 1];
        size_t len = end - start;
        p1s = subset(p1, start, end);
        std::sort(p1s.begin(), p1s.end());
        std::memcpy(p2.data() + start, p1s.data(), len * sizeof(double));
      }

      // Compute minimum ratio
      double min_val = 1.0;
      for (size_t k = 0; k < n; ++k) {
        double ratio = p2[k] / (w[k] * c[k]);
        min_val = std::min(min_val, ratio);
      }

      pinter(iter, i) = min_val;
    }

    // Store incidence
    for (size_t j = 0; j < m; ++j) {
      incid(i, j) = cc[j];
    }
  }

  // Compute adjusted p-values for elementary hypotheses
  FlatMatrix padj(nreps, m);
  for (size_t j = 0; j < m; ++j) {
    for (size_t iter = 0; iter < nreps; ++iter) {
      double max_p = 0.0;
      for (size_t i = 0; i < ntests; ++i) {
        if (incid(i, j)) {
          max_p = std::max(max_p, pinter(iter, i));
        }
      }
      padj(iter, j) = std::min(max_p, 1.0);
    }
  }

  // Apply logical restrictions (serial constraints)
  for (size_t j = 0; j < m; ++j) {
    int serial_sum = 0;
    for (size_t k = 0; k < m; ++k) {
      serial_sum += serial(j, k);
    }

    if (serial_sum > 0) {
      for (size_t iter = 0; iter < nreps; ++iter) {
        double pre = 0.0;
        for (size_t k = 0; k < m; ++k) {
          if (serial(j, k)) {
            pre = std::max(pre, padj(iter, k));
          }
        }
        padj(iter, j) = std::max(padj(iter, j), pre);
      }
    }
  }

  // Apply logical restrictions (parallel constraints)
  for (size_t j = 0; j < m; ++j) {
    int parallel_sum = 0;
    for (size_t k = 0; k < m; ++k) {
      parallel_sum += parallel(j, k);
    }

    if (parallel_sum > 0) {
      for (size_t iter = 0; iter < nreps; ++iter) {
        double pre = 1.0;
        for (size_t k = 0; k < m; ++k) {
          if (parallel(j, k)) {
            pre = std::min(pre, padj(iter, k));
          }
        }
        padj(iter, j) = std::max(padj(iter, j), pre);
      }
    }
  }

  return AdjustedPValues{incid, pinter, padj};
}

// [[Rcpp::export]]
Rcpp::List fmodmixRcpp(const Rcpp::NumericMatrix &p,
                       const Rcpp::LogicalMatrix &family,
                       const Rcpp::LogicalMatrix &serial,
                       const Rcpp::LogicalMatrix &parallel,
                       const Rcpp::NumericVector &gamma,
                       const std::string &test = "hommel",
                       const bool exhaust = true) {
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto family1 = boolmatrix_from_Rmatrix(family);
  auto serial1 = boolmatrix_from_Rmatrix(serial);
  auto parallel1 = boolmatrix_from_Rmatrix(parallel);
  auto gamma1 = Rcpp::as<std::vector<double>>(gamma);
  auto out = fmodmixcpp(p1, family1, serial1, parallel1, gamma1, test, exhaust);
  ListCpp result;
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.pinter), "pinter");
  result.push_back(std::move(out.padj), "padj");
  return Rcpp::wrap(result);
}

// Helper to compute truncated adjusted p-values for multiple testing
AdjustedPValues ftrunccpp(const FlatMatrix &p, const std::string &test,
                          const double gamma) {

  // Validate inputs
  size_t niters = p.nrow;
  size_t m = p.ncol;

  if (gamma < 0.0 || gamma > 1.0) {
    throw std::invalid_argument("gamma must be between 0 and 1");
  }

  // Validate p-values in [0, 1]
  if (std::any_of(p.data.begin(), p.data.end(),
                  [](double v) { return v < 0.0 || v > 1.0; })) {
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

  size_t ntests = (1 << m) - 1; // 2^m - 1

  // Incidence matrix (ntests x m) - which hypotheses are in each intersection
  IntMatrix incid(ntests, m);

  // Intersection p-values (niters x ntests)
  FlatMatrix pinter(niters, ntests);

  // Precompute constant factors
  double tbon = (1.0 - gamma) / m;

  // Preallocate reusable vectors
  std::vector<size_t> hyp_indices;
  hyp_indices.reserve(m);
  std::vector<double> p1;
  p1.reserve(m);

  // Process each intersection hypothesis
  for (size_t i = 0; i < ntests; ++i) {
    size_t number = ntests - i;

    // Binary representation of elementary hypotheses in intersection
    std::vector<int> cc(m);
    size_t k = 0; // count of hypotheses in this intersection
    hyp_indices.clear();

    for (size_t j = 0; j < m; ++j) {
      cc[j] = (number >> (m - 1 - j)) & 1;
      if (cc[j]) {
        hyp_indices.push_back(j);
        ++k;
      }
    }

    // Store incidence (column-major: column j base = incid_data + j * ntests)
    for (size_t j = 0; j < m; ++j) {
      incid(i, j) = cc[j];
    }

    // Precompute weights if they're constant across iterations
    std::vector<double> weights(k);
    if (test1 == "hommel") {
      for (size_t j = 0; j < k; ++j) {
        weights[j] = (j + 1) * gamma / k + tbon;
      }
    } else if (test1 == "hochberg") {
      for (size_t j = 0; j < k; ++j) {
        weights[j] = gamma / (k - j) + tbon;
      }
    } else { // holm
      double w_holm = gamma / k + tbon;
      std::fill(weights.begin(), weights.end(), w_holm);
    }

    // Process each iteration/replication
    for (size_t iter = 0; iter < niters; ++iter) {
      // Extract p-values for hypotheses in this intersection
      p1.resize(k);
      for (size_t j = 0; j < k; ++j) {
        p1[j] = p(iter, hyp_indices[j]);
      }

      // Sort p-values
      std::sort(p1.begin(), p1.end());

      // Compute minimum ratio
      double q = 1.0;
      for (size_t j = 0; j < k; ++j) {
        double ratio = p1[j] / weights[j];
        q = std::min(q, ratio);
      }

      pinter(iter, i) = q;
    }
  }

  // Compute adjusted p-values for elementary hypotheses
  FlatMatrix padj(niters, m);

  // For each hypothesis j, find max pinter over intersections containing j
  for (size_t j = 0; j < m; ++j) {
    for (size_t iter = 0; iter < niters; ++iter) {
      double max_p = 0.0;
      for (size_t i = 0; i < ntests; ++i) {
        if (incid(i, j)) {
          max_p = std::max(max_p, pinter(iter, i));
        }
      }

      padj(iter, j) = max_p;
    }
  }

  return AdjustedPValues{incid, pinter, padj};
}

// [[Rcpp::export]]
Rcpp::List ftruncRcpp(const Rcpp::NumericMatrix &p,
                      const std::string &test = "hommel",
                      const double gamma = 1.0) {
  auto p1 = flatmatrix_from_Rmatrix(p);
  auto out = ftrunccpp(p1, test, gamma);
  ListCpp result;
  result.push_back(std::move(out.inthyp), "inthyp");
  result.push_back(std::move(out.pinter), "pinter");
  result.push_back(std::move(out.padj), "padj");
  return Rcpp::wrap(result);
}
