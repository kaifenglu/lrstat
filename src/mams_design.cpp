#include "dataframe_list.h"
#include "thread_utils.h"
#include "utilities.h"
#include "generic_design.h"
#include "mvnormr.h"

#include <Rcpp.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using std::size_t;


ListCpp exitprob_mams_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t kMax,
    const FlatMatrix& b,
    const FlatMatrix& a,
    const std::vector<double>& I) {

  if (M < 1) throw std::invalid_argument("M should be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r should be positive");
  if (theta.size() != M) throw std::invalid_argument("theta should have length M");
  if (kMax < 1) throw std::invalid_argument("kMax should be at least 1");
  if (!none_na(b.data)) throw std::invalid_argument("b must be provided");
  if (!none_na(a.data)) throw std::invalid_argument("a must be provided");
  if (b.nrow != M) throw std::invalid_argument("b should have M rows");
  if (b.ncol < kMax) throw std::invalid_argument("Insufficient columns for b");
  if (a.nrow != M) throw std::invalid_argument("a should have M rows");
  if (a.ncol < kMax) throw std::invalid_argument("Insufficient columns for a");
  for (size_t k = 0; k < kMax; ++k) {
    for (size_t m = 0; m < M; ++m) {
      if (b(m, k) < a(m, k)) {
        throw std::invalid_argument("Each b(m, k) must be >= a(m, k)");
      }
    }
  }

  std::vector<double> Ivec; Ivec.reserve(kMax);
  if (none_na(I)) {
    if (I.size() < kMax) throw std::invalid_argument("Insufficient length for I");
    Ivec.assign(I.begin(), I.begin() + kMax);
    if (Ivec[0] <= 0.0) throw std::invalid_argument("I must be positive");
    if (any_nonincreasing(Ivec)) throw std::invalid_argument("I must be increasing");
  } else {
    Ivec.resize(kMax);
    std::iota(Ivec.begin(), Ivec.end(), 1.0);
  }

  double rho = corr_known ? r / (r + 1.0) : 0;

  std::vector<double> sqrtI(kMax);
  for (size_t k = 0; k < kMax; ++k) {
    sqrtI[k] = std::sqrt(Ivec[k]);
  }

  auto lowest_bit_index = [](std::uint64_t x) -> int {
    // x must be non-zero. trailing-zero count by shifting.
    int pos = 0;
    while ((x & 1ULL) == 0ULL) { x >>= 1; ++pos; }
    return pos;
  };


  std::vector<double> exitProbUpper(kMax), exitProbLower(kMax);

  if (rho == 0 || M == 1) {
    double T = Ivec.back();
    std::vector<double> s(kMax);
    for (size_t k = 0; k < kMax; ++k) {
      s[k] = Ivec[k] / T;
    }

    std::size_t fullN = M * kMax;
    std::vector<double> upper(fullN);
    std::vector<double> lower(fullN, -8.0);
    for (size_t k = 0; k < kMax; ++k) {
      for (size_t m = 0; m < M; ++m) {
        size_t idx = k * M + m;
        upper[idx] = b(m, k) - theta[m] * sqrtI[k];
        if (a(m, k) > -8.0) lower[idx] = a(m, k) - theta[m] * sqrtI[k];
      }
    }

    ListCpp probs;
    std::vector<double> v;
    std::vector<double> upper1, lower1, c;
    std::vector<double> upper1a, upper2;
    std::vector<double> lo(kMax, -8.0);
    std::vector<double> zero(kMax, 0.0);
    for (size_t k = 1; k <= kMax; ++k) {
      std::size_t nk = M * k;

      upper1.resize(nk); lower1.resize(nk);
      std::memcpy(upper1.data(), upper.data(), nk * sizeof(double));
      std::memcpy(lower1.data(), lower.data(), nk * sizeof(double));

      upper1a.resize(k);
      upper2.resize(k - 1);

      // --- start of signed sums using Gray code ---

      // working buffer (contiguous), initialized to b for first k - 1 looks
      c.resize(nk);
      size_t nk2 = (k - 1) * M;
      if (k > 1) std::memcpy(c.data(), upper1.data(), nk2 * sizeof(double));

      // collect active indices among 0..k-2 where a[i] != -8
      std::vector<std::size_t> active;
      for (std::size_t i = 0; i + 1 < k; ++i) {
        if (a(0, i) == -8.0) {
          // forced to b[i] (already set)
        } else {
          active.push_back(i);
          // initialize to a[i] for starting gray == 0 (all a)
          size_t j = i * M; // start of the block for index i
          std::memcpy(c.data() + j, lower1.data() + j, M * sizeof(double));
        }
      }

      const std::size_t m = active.size();
      const std::uint64_t totalComb = (m == 0) ? 1ULL : (1ULL << m);
      double vupper1 = 0.0, vupper2 = 0.0, vlower = 0.0;

      // initial sign: choose a for every active index, each contributes -1 => (-1)^m
      double sign = ( (m % 2) ? -1.0 : 1.0 );

      std::uint64_t prevGray = 0;

      for (std::uint64_t i = 0; i < totalComb; ++i) {
        const std::uint64_t gray = i ^ (i >> 1);
        if (i == 0) {
          // c already initialized for gray == 0
        } else {
          const std::uint64_t changed = prevGray ^ gray; // exactly one bit set
          const int bitPos = lowest_bit_index(changed);
          const std::size_t idx = active[bitPos];
          const bool choseB = ((gray >> bitPos) & 1ULL) != 0ULL;
          size_t j = idx * M; // start of the block for index idx
          if (choseB) { // copy b[idx] into c
            std::memcpy(c.data() + j, upper1.data() + j, M * sizeof(double));
          } else { // copy a[idx] into c
            std::memcpy(c.data() + j, lower1.data() + j, M * sizeof(double));
          }
          // flip sign when one factor toggles between -1 and +1
          sign = -sign;
        }

        double pupper1 = 1.0, pupper2 = 1.0, plower = 1.0;
        for (size_t m = 0; m < M; ++m) {
          if (k > 1) {
            for(size_t j = 0; j + 1 < k; ++j) {
              upper2[j] = c[j * M + m];
              if (upper2[j] < -8.0) upper2[j] = -8.0;
            }
            probs = exitprobcpp(upper2, lo, zero, Ivec);
            v = probs.get<std::vector<double>>("exitProbUpper");
            pupper1 *= (1.0 - std::accumulate(v.begin(), v.end(), 0.0));
            std::memcpy(upper1a.data(), upper2.data(), (k - 1) * sizeof(double));
          }

          upper1a[k - 1] = upper1[(k - 1) * M + m];
          if (upper1a[k - 1] < -8.0) upper1a[k - 1] = -8.0;
          probs = exitprobcpp(upper1a, lo, zero, Ivec);
          v = probs.get<std::vector<double>>("exitProbUpper");
          pupper2 *= (1.0 - std::accumulate(v.begin(), v.end(), 0.0));

          upper1a[k - 1] = lower1[(k - 1) * M + m];
          if (upper1a[k - 1] < -8.0) upper1a[k - 1] = -8.0;
          probs = exitprobcpp(upper1a, lo, zero, Ivec);
          v = probs.get<std::vector<double>>("exitProbUpper");
          plower *= (1.0 - std::accumulate(v.begin(), v.end(), 0.0));
        }

        if (k > 1) vupper1 += sign * pupper1;
        vupper2 += sign * pupper2;
        vlower += sign * plower;
        prevGray = gray;
      }

      if (k > 1) exitProbUpper[k - 1] = vupper1 - vupper2;
      else exitProbUpper[k - 1] = 1.0 - vupper2;

      exitProbLower[k - 1] = vlower;
    }
  } else if (kMax == 1) {
    double sqrtrho = std::sqrt(rho);
    double sqrt1mr = std::sqrt(1.0 - rho);
    double sqrtI0 = sqrtI[0];
    auto f1 = [&](double z0) {
      double value = 1.0;
      for (size_t m = 0; m < M; ++m) {
        double upper = (b(m, 0) - theta[m] * sqrtI0 - sqrtrho * z0) / sqrt1mr;
        double mass = boost_pnorm(upper);
        if (mass <= 0.0) return 0.0;
        value *= mass;
      }
      return value * boost_dnorm(z0);
    };

    std::vector<double> breaks = { -8.0, 0.0, 8.0 };
    double p1 = integrate3(f1, breaks, 1e-6);
    exitProbUpper[0] = 1.0 - p1;

    auto f2 = [&](double z0) {
      double value = 1.0;
      for (size_t m = 0; m < M; ++m) {
        double upper = (a(m, 0) - theta[m] * sqrtI0 - sqrtrho * z0) / sqrt1mr;
        double mass = boost_pnorm(upper);
        if (mass <= 0.0) return 0.0;
        value *= mass;
      }
      return value * boost_dnorm(z0);
    };

    double p2 = integrate3(f2, breaks, 1e-6);
    exitProbLower[0] = p2;
  } else {
    double T = Ivec.back();
    double sqrtT = sqrtI.back();
    std::vector<double> s(kMax);
    for (size_t k = 0; k < kMax; ++k) {
      s[k] = Ivec[k] / T;
    }

    std::size_t fullN = M * kMax;
    std::vector<double> upper(fullN);
    std::vector<double> lower(fullN, -8.0);
    for (size_t k = 0; k < kMax; ++k) {
      for (size_t m = 0; m < M; ++m) {
        size_t idx = k * M + m;
        upper[idx] = b(m, k) * sqrt(s[k]) - theta[m] * s[k] * sqrtT;
        if (a(m, k) > -8.0)
          lower[idx] = a(m, k) * sqrt(s[k]) - theta[m] * s[k] * sqrtT;
      }
    }

    std::vector<double> zero(fullN, 0.0);

    // covariance matrix for scaled Brownian motion
    FlatMatrix sigma(fullN, fullN);

    for (size_t k2 = 0; k2 < kMax; ++k2) {
      for (size_t k1 = 0; k1 < kMax; ++k1) {
        size_t k = std::min(k1, k2);
        double val_off = rho * s[k];    // value for off-diagonal (m1 != m2)
        double val_diag = s[k];         // diagonal (m1 == m2)
        // block row base and block col base:
        std::size_t row_base = k1 * M;
        std::size_t col_base = k2 * M;
        // for each column (m2) inside the MxM block, fill contiguous column segment:
        for (size_t m2 = 0; m2 < M; ++m2) {
          // pointer to start of the column (col_index) at row row_base
          double* dst_col = sigma.data_ptr() + (col_base + m2) * fullN + row_base;
          // fill M entries in that column with off-diagonal value
          std::fill_n(dst_col, M, val_off);
          // overwrite the diagonal element (at offset m2 inside this block)
          dst_col[m2] = val_diag;
        }
      }
    }


    std::vector<double> upper1, lower1, c, zero1, lower1m;
    std::vector<double> upper2, lower2, zero2;
    FlatMatrix sigma1, sigma2;
    for (size_t k = 1; k <= kMax; ++k) {
      std::size_t nk = M * k;

      upper1.resize(nk); lower1.resize(nk); zero1.resize(nk); lower1m.resize(nk);
      std::memcpy(upper1.data(), upper.data(), nk * sizeof(double));
      std::memcpy(lower1.data(), lower.data(), nk * sizeof(double));
      std::memcpy(zero1.data(), zero.data(), nk * sizeof(double));
      std::fill_n(lower1m.data(), nk, -8.0);

      sigma1.resize(nk, nk);               // reuses capacity if already large enough
      // copy top-left into sigma1 (column-major): one memcpy per column
      for (std::size_t j = 0; j < nk; ++j) {
        const double* src = sigma.data_ptr() + j * fullN; // start of col c in sigma
        double* dst = sigma1.data_ptr() + j * nk;         // start of col c in sigma1
        std::memcpy(dst, src, nk * sizeof(double));
      }

      size_t nk2 = (k - 1) * M;
      if (k > 1) {
        upper2.resize(nk2); lower2.resize(nk2); zero2.resize(nk2);
        std::fill_n(lower2.data(), nk2, -8.0);
        std::memcpy(zero2.data(), zero1.data(), nk2 * sizeof(double));

        sigma2.resize(nk2, nk2);
        // copy top-left into sigma2 (column-major): one memcpy per column
        for (std::size_t j = 0; j < nk2; ++j) {
          const double* src = sigma1.data_ptr() + j * nk;
          double* dst = sigma2.data_ptr() + j * nk2;
          std::memcpy(dst, src, nk2 * sizeof(double));
        }
      }

      // --- start of signed sums using Gray code ---

      // working buffer (contiguous), initialized to b for first k - 1 looks
      c.resize(nk);
      if (k > 1) std::memcpy(c.data(), upper1.data(), nk2 * sizeof(double));

      // collect active indices among 0..k-2 where a[i] != -8
      std::vector<std::size_t> active;
      for (std::size_t i = 0; i + 1 < k; ++i) {
        if (a(0, i) == -8.0) {
          // forced to b[i] (already set)
        } else {
          active.push_back(i);
          // initialize to a[i] for starting gray == 0 (all a)
          size_t j = i * M; // start of the block for index i
          std::memcpy(c.data() + j, lower1.data() + j, M * sizeof(double));
        }
      }

      const std::size_t m = active.size();
      const std::uint64_t totalComb = (m == 0) ? 1ULL : (1ULL << m);
      double vupper1 = 0.0, vupper2 = 0.0, vlower = 0.0;

      // initial sign: choose a for every active index, each contributes -1 => (-1)^m
      double sign = ( (m % 2) ? -1.0 : 1.0 );

      std::uint64_t prevGray = 0;
      for (std::uint64_t i = 0; i < totalComb; ++i) {
        const std::uint64_t gray = i ^ (i >> 1);
        if (i == 0) {
          // c already initialized for gray == 0
        } else {
          const std::uint64_t changed = prevGray ^ gray; // exactly one bit set
          const int bitPos = lowest_bit_index(changed);
          const std::size_t idx = active[bitPos];
          const bool choseB = ((gray >> bitPos) & 1ULL) != 0ULL;
          size_t j = idx * M; // start of the block for index idx
          if (choseB) { // copy b[idx] into c
            std::memcpy(c.data() + j, upper1.data() + j, M * sizeof(double));
          } else { // copy a[idx] into c
            std::memcpy(c.data() + j, lower1.data() + j, M * sizeof(double));
          }
          // flip sign when one factor toggles between -1 and +1
          sign = -sign;
        }

        if (k > 1) {
          std::memcpy(upper2.data(), c.data(), nk2 * sizeof(double));
          for (size_t j = 0; j < nk2; ++j) {
            if (upper2[j] < -8.0) upper2[j] = -8.0;
          }
          auto rupper1 = pmvnormcpp(lower2, upper2, zero2, sigma2,
                                    1024, 16384, 8, 1e-4, 0.0, 314159, true);
          vupper1 += sign * rupper1.prob;
        }

        std::memcpy(c.data() + nk2, upper1.data() + nk2, M * sizeof(double));
        for (size_t j = 0; j < M; ++j) {
          if (c[nk2 + j] < -8.0) c[nk2 + j] = -8.0;
        }
        auto rupper2 = pmvnormcpp(lower1m, c, zero1, sigma1,
                                  1024, 16384, 8, 1e-4, 0.0, 314159, true);
        vupper2 += sign * rupper2.prob;

        std::memcpy(c.data() + nk2, lower1.data() + nk2, M * sizeof(double));
        for (size_t j = 0; j < M; ++j) {
          if (c[nk2 + j] < -8.0) c[nk2 + j] = -8.0;
        }

        auto rlower = pmvnormcpp(lower1m, c, zero1, sigma1,
                                 1024, 16384, 8, 1e-4, 0.0, 314159, true);
        vlower += sign * rlower.prob;

        prevGray = gray;
      }

      if (k > 1) exitProbUpper[k - 1] = vupper1 - vupper2;
      else exitProbUpper[k - 1] = 1.0 - vupper2;

      exitProbLower[k - 1] = vlower;
    }
  }

  ListCpp exitProb;
  exitProb.push_back(std::move(exitProbUpper), "exitProbUpper");
  exitProb.push_back(std::move(exitProbLower), "exitProbLower");
  return exitProb;
}


// overload with default a matrix of -8.0
ListCpp exitprob_mams_cpp(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const bool corr_known,
    const size_t kMax,
    const FlatMatrix& b,
    const std::vector<double>& I) {
  if (!none_na(b.data)) throw std::invalid_argument("b must be provided");
  double bmin = *std::min_element(b.data_ptr(), b.data_ptr() + b.nrow * b.ncol);
  double amin = std::min(-8.0, bmin);
  FlatMatrix a(M, kMax); a.fill(amin);
  return exitprob_mams_cpp(M, r, theta, corr_known, kMax, b, a, I);
}


// helper function
FlatMatrix as_boundary_matrix(SEXP x, const char* arg, size_t M, size_t K,
                              bool allow_default = false) {
  if (x == R_NilValue) {
    if (allow_default) {
      FlatMatrix out(M, K);
      out.fill(-8.0);
      return out;
    }
    throw std::invalid_argument(std::string(arg) + " must be provided");
  }

  if (Rf_isMatrix(x)) {
    Rcpp::NumericMatrix mat(x);
    if (static_cast<size_t>(mat.nrow()) != M ||
        static_cast<size_t>(mat.ncol()) != K) {
      throw std::invalid_argument(std::string(arg) +
                                  " matrix must have dimension M x kMax");
    }
    return flatmatrix_from_Rmatrix(mat);
  }

  if (Rf_isNumeric(x)) {
    Rcpp::NumericVector vec(x);
    FlatMatrix out(M, K);

    if (vec.size() == 1) {
      out.fill(vec[0]);
      return out;
    }

    if (static_cast<size_t>(vec.size()) == K) {
      for (size_t k = 0; k < K; ++k) {
        for (size_t m = 0; m < M; ++m) {
          out(m, k) = vec[k];
        }
      }
      return out;
    }

    if (static_cast<size_t>(vec.size()) == M * K) {
      for (size_t k = 0; k < K; ++k) {
        for (size_t m = 0; m < M; ++m) {
          out(m, k) = vec[m + M * k];
        }
      }
      return out;
    }

    throw std::invalid_argument(
        std::string(arg) + " vector must have length 1, kMax, or M*kMax");
  }

  throw std::invalid_argument(
      std::string(arg) + " must be either a numeric vector or numeric matrix");
}

//' @title Exit Probabilities for a Multi-Arm Multi-Stage Design
//' @description Computes the exit (rejection) probabilities for a multi-arm
//' multi-stage design.
//'
//' @param M Number of active treatment arms.
//' @param r Randomization ratio of each active arm to the common control.
//' @param theta A vector of length \eqn{M} representing the true treatment
//'   effects for each active arm versus the common control.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param kMax Number of sequential looks.
//' @param b A vector of efficacy boundaries for the max-Z statistics.
//' @param a A vector of futility boundaries for the max-Z statistics.
//' @param I A vector of information levels for any active arm versus the
//'   common control.
//'
//' @details
//' The function assumes a multivariate normal distribution for the Wald
//' statistics and all active arms share the same information level.
//'
//' @return A vector \code{exitProb} of length \code{kMax} containing the
//' probability of rejection at each look.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @examples
//'
//' # Setup: 2 active arms vs control and 3 sequential looks.
//'
//' # Information levels: equal spacing over 3 looks based on a maximum of
//' # 95 patients per arm, SD = 1.0
//' I <- 95 / (2 * 1.0^2) * seq(1, 3)/3
//'
//' # O'Brien-Fleming critical values
//' b <- c(3.886562, 2.748214, 2.243907)
//'
//' # Type I error under the global null hypothesis
//' p0 <- exitprob_mams(M = 2, theta = c(0, 0), kMax = 3, b = b, I = I)
//' cumsum(p0$exitProbUpper)
//'
//' # Power under alternative: Treatment effects of 0.3 and 0.5
//' p1 <- exitprob_mams(M = 2, theta = c(0.3, 0.5), kMax = 3, b = b, I = I)
//' cumsum(p1$exitProbUpper)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List exitprob_mams(
    const int M = NA_INTEGER,
    const double r = 1,
    const Rcpp::NumericVector& theta = NA_REAL,
    const bool corr_known = true,
    const int kMax = NA_INTEGER,
    SEXP b = R_NilValue,
    SEXP a = R_NilValue,
    SEXP I = R_NilValue) {

  if (M == NA_INTEGER || M < 1) {
    throw std::invalid_argument("M must be a positive integer");
  }
  if (kMax == NA_INTEGER || kMax < 1) {
    throw std::invalid_argument("kMax must be a positive integer");
  }

  size_t M_ = static_cast<size_t>(M);
  size_t K = static_cast<size_t>(kMax);

  std::vector<double> thetaVec(theta.begin(), theta.end());

  std::vector<double> IVec;
  if (I == R_NilValue) {
    IVec = std::vector<double>(1, NA_REAL);
  } else {
    Rcpp::NumericVector Iv(I);
    IVec.assign(Iv.begin(), Iv.end());
  }

  FlatMatrix bMat = as_boundary_matrix(b, "b", M_, K, false);
  FlatMatrix aMat = as_boundary_matrix(a, "a", M_, K, true);

  auto probs = exitprob_mams_cpp(M_, r, thetaVec, corr_known, K, bMat, aMat, IVec);
  return Rcpp::wrap(probs);
}



std::vector<double> getBound_mams_cpp(
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t k,
    const std::vector<double>& informationRates,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const std::vector<unsigned char>& efficacyStopping) {

  if (M < 1) throw std::invalid_argument("M should be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r should be positive");
  if (k < 1) throw std::invalid_argument("k should be at least 1");
  if (alpha < 0.00001 || alpha >= 1) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }

  // infoRates: if missing create 1/k, 2/k, ..., k/k
  size_t kMax = k;
  std::vector<double> infoRates;
  if (none_na(informationRates)) {
    kMax = informationRates.size();
    if (kMax < k)
      throw std::invalid_argument("Insufficient length for informationRates");

    infoRates = informationRates;
    if (infoRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(infoRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (infoRates.back() > 1.0)
      throw std::invalid_argument("informationRates must not exceed 1");
  } else {
    infoRates.resize(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(kMax);
    }
  }

  // spendTime: default to infoRates if missing
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() < kMax)
      throw std::invalid_argument("Insufficient length for spendingTime");

    spendTime = spendingTime;
    if (spendTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendTime.back() > 1.0)
      throw std::invalid_argument("spendingTime must not exceed 1");
  } else {
    spendTime = infoRates;
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() < kMax)
      throw std::invalid_argument("Insufficient length for efficacyStopping");

    effStopping = efficacyStopping;
    if (effStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
  } else {
    effStopping.assign(kMax, 1);
  }

  // asf (alpha spending function) to lower-case
  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  // userAlphaSpending checks when asf == "user"
  if (asf == "user") {
    if (!none_na(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be specified");
    if (userAlphaSpending.size() < kMax)
      throw std::invalid_argument("Insufficient length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending.back() > alpha)
      throw std::invalid_argument("userAlphaSpending must not exceed alpha");
  }

  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }

  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  if (asf == "of" || asf == "p" || asf == "wt" || asf == "none") {
    if (infoRates.back() != 1.0) {
      throw std::invalid_argument(
          "informationRates must end with 1 for OF, P, WT, or NONE");
    }
    if (spendTime.back() != 1.0) {
      throw std::invalid_argument(
          "spendingTime must end with 1 for OF, P, WT, or NONE");
    }
  }

  std::vector<double> zero(M, 0.0);
  std::vector<double> criticalValues(kMax, 8.0);
  FlatMatrix b(M, kMax);
  b.fill(8.0);
  ListCpp probs;
  std::vector<double> v;

  if (asf == "none") {
    double* colptr = b.data_ptr() + (kMax - 1) * M;
    auto f = [&](double x)->double {
      std::fill_n(colptr, M, x);
      probs = exitprob_mams_cpp(M, r, zero, corr_known, kMax, b, infoRates);
      v = probs.get<std::vector<double>>("exitProbUpper");
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    criticalValues[kMax-1] = brent(f, 0.0, 8.0, 1e-6);
    return subset(criticalValues, 0, k);
  }

  if (asf == "of" || asf == "p" || asf == "wt") {
    double Delta;
    if (asf == "of") Delta = 0.0;
    else if (asf == "p") Delta = 0.5;
    else Delta = parameterAlphaSpending; // parameterAlphaSpending holds delta for WT

    // for a given multiplier, compute cumulative upper exit probability - alpha
    std::vector<double> u0(kMax);
    for (size_t i = 0; i < kMax; ++i) {
      u0[i] = std::pow(infoRates[i], Delta - 0.5);
    }

    auto f = [&](double x)->double {
      for (size_t i = 0; i < kMax; ++i) {
        if (!effStopping[i]) continue;
        double val = x * u0[i];
        std::fill_n(b.data_ptr() + i * M, M, val); // contiguous write of M doubles
      }

      probs = exitprob_mams_cpp(M, r, zero, corr_known, kMax, b, infoRates);
      v = probs.get<std::vector<double>>("exitProbUpper");
      double cpu = std::accumulate(v.begin(), v.end(), 0.0);
      return cpu - alpha;
    };

    double cwt = brent(f, 0.0, 10.0, 1e-6);
    for (size_t i = 0; i < kMax; ++i) {
      if (effStopping[i]) criticalValues[i] = cwt * u0[i];
    }
    return subset(criticalValues, 0, k);
  }

  if (asf == "sfof" || asf == "sfp" || asf == "sfkd" ||
      asf == "sfhsd" || asf == "user") {
    // stage 1
    double cumAlpha;
    if (asf == "user") cumAlpha = userAlphaSpending[0];
    else cumAlpha = errorSpentcpp(spendTime[0], alpha, asf, parameterAlphaSpending);

    if (!effStopping[0]) criticalValues[0] = 8.0;
    else {
      std::vector<double> mean(M, 0.0);
      FlatMatrix sigma(M, M);
      for (size_t j = 0; j < M; ++j) {
        for (size_t i = 0; i < M; ++i) {
          sigma(i, j) = (i == j) ? 1.0 : r / (r + 1.0);
        }
      }
      criticalValues[0] = qmvnormcpp(1.0 - cumAlpha, mean, sigma,
                                     1024, 16384, 8, 1e-4, 0.0, 314159, true);
    }

    // Preallocate reusable buffers used by the root-finding lambda
    // subsequent stages
    for (size_t k1 = 1; k1 < kMax; ++k1) {
      if (!effStopping[k1]) continue;

      // determine cumulative alpha at this stage
      if (asf == "user") cumAlpha = userAlphaSpending[k1];
      else cumAlpha = errorSpentcpp(spendTime[k1], alpha, asf,
                                    parameterAlphaSpending);

      // Fill columns k1-1 with their known critical values
      double val = criticalValues[k1 - 1];
      double* col_ptr = b.data_ptr() + (k1 - 1) * M;   // column-major contiguous
      std::fill_n(col_ptr, M, val);              // fast contiguous write

      // Define lambda that only sets the last column of b
      double* last_col = b.data_ptr() + k1 * M;
      auto f = [&](double x)->double {
        // set the last column to the current candidate critical value
        std::fill_n(last_col, M, x);  // fill last column fast
        probs = exitprob_mams_cpp(M, r, zero, corr_known, k1 + 1, b, infoRates);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - cumAlpha;
      };

      double f_8 = f(8.0);
      if (f_8 > 0.0) { // no alpha spent at current visit
        criticalValues[k1] = 8.0;
      } else {
        auto f_for_brent = [&](double x)->double {
          if (x == 8.0) return f_8; // avoid recomputation at 8.0
          return f(x);
        };

        criticalValues[k1] = brent(f_for_brent, 0.0, 8.0, 1e-6);
      }
    }

    return subset(criticalValues, 0, k);
  }

  throw std::invalid_argument("Invalid value for typeAlphaSpending");
}


//' @title Efficacy Boundaries for a Multi-Arm Multi-Stage Design
//' @description Calculates the efficacy stopping boundaries for a multi-arm
//' multi-stage design.
//'
//' @param M Number of active treatment arms.
//' @param r Randomization ratio of each active arm to the common control.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is assumed.
//' @param k The index of the current look.
//' @param informationRates A numeric vector of information rates up to the
//'   current look. Values must be strictly increasing and \eqn{\le 1}.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param spendingTime A numeric vector of length \eqn{k} specifying the
//'   error spending time at each analysis. Values must be strictly increasing
//'   and \eqn{\le 1}. If omitted, defaults to \code{informationRates}.
//' @inheritParams param_efficacyStopping
//'
//' @details
//' The function determines critical values by solving for the boundary that
//' satisfies the alpha-spending requirement.
//'
//' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
//' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
//' and \code{spendingTime} must be of full length \code{kMax}, and
//' \code{informationRates} and \code{spendingTime} must end with 1.
//'
//' @return A numeric vector of length \eqn{k} containing the critical
//' values (on the standard normal Z-scale) for each analysis up to the
//' current look.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @examples
//'
//' # Determine O'Brien-Fleming boundaries for a TSSSD with
//' # 2 active arms and 3 looks.
//' getBound_mams(M = 2, k = 3, informationRates = seq(1, 3)/3,
//'               alpha = 0.025, typeAlphaSpending = "OF")
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector getBound_mams(
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int k = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL) {

  std::vector<double> infoRates(informationRates.begin(), informationRates.end());
  std::vector<double> userAlpha(userAlphaSpending.begin(), userAlphaSpending.end());
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());
  auto effStopping = convertLogicalVector(efficacyStopping);

  auto result = getBound_mams_cpp(
    static_cast<size_t>(M), r, corr_known, static_cast<size_t>(k),
    infoRates, alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlpha, spendTime, effStopping
  );

  return Rcpp::wrap(result);
}

// The following function computes the conditional power and futility
// bounds of a secondary trial given the interim result of a primary
// trial and the beta-spending function of the secondary trial, where
// IL and zL are the information and Wald statistics at the interim look
// of the primary trial, M, r, theta, alpha, kMax, bsf, bsfpar, st,
// futStopping are parameters of the secondary trial, critValues,
// futBounds and I are boundaries and information of integrated trial.
ListCpp getPower_mams(
    const size_t M,
    const double r,
    const std::vector<double>& theta,
    const double alpha,
    const size_t kMax,
    const std::vector<double>& critValues,
    const std::vector<double>& I,
    const std::string& bsf,
    const double bsfpar,
    const std::vector<double>& st,
    const std::vector<unsigned char>& futStopping,
    const double IL,
    const std::vector<double>& zL) {

  std::vector<double> I2(kMax), sqrtI2(kMax), sqrtI(kMax);
  for (size_t k = 0; k < kMax; ++k) {
    I2[k] = I[k] - IL;
    sqrtI2[k] = std::sqrt(I2[k]);
    sqrtI[k] = std::sqrt(I[k]);
  }

  double sqrtIL = std::sqrt(IL);
  std::vector<double> zscaled(M);
  for (size_t m = 0; m < M; ++m) {
    zscaled[m] = zL[m] * sqrtIL;
  }

  FlatMatrix a(M, kMax);
  a.fill(-8.0);

  FlatMatrix b(M, kMax);
  for (size_t k = 0; k < kMax; ++k) {
    for (size_t m = 0; m < M; ++m) {
      b(m, k) = (critValues[k] * sqrtI[k] - zscaled[m]) / sqrtI2[k];
    }
  }

  std::vector<double> lo(M, -8.0), hi(M);
  std::vector<double> mu0(M);
  for (size_t m = 0; m < M; ++m) {
    mu0[m] = theta[m] * sqrtI2[0];
  }

  FlatMatrix sigma(M, M);
  sigma.fill(r / (r + 1.0));
  for (size_t m = 0; m < M; ++m) {
    sigma(m, m) = 1.0;
  }

  // reusable buffers for prefixes
  ListCpp probs;
  std::vector<double> v;
  std::vector<double> futBounds(kMax);
  auto f = [&](double x) -> double {
    // reset futility bounds
    std::fill(futBounds.begin(), futBounds.end(), -8.0);
    double eps = 0.0, cb = 0.0;

    // first stage
    if (futStopping[0]) {
      cb = errorSpentcpp(st[0], x, bsf, bsfpar);

      auto g = [&](double aval) -> double {
        for (size_t m = 0; m < M; ++m) {
          hi[m] = (aval * sqrtI[0] - zscaled[m]) / sqrtI2[0];
        }

        auto q = pmvnormcpp(lo, hi, mu0, sigma, 1024, 16384, 8, 1e-4, 0.0, 314159);
        return q.prob - cb;
      };

      eps = g(critValues[0]);
      if (eps < 0.0) return -1.0; // to decrease drift
      futBounds[0] = brent(g, -8.0, critValues[0], 1e-6);
    }

    // subsequent stages
    for (size_t k = 1; k < kMax; ++k) {
      if (futStopping[k]) {
        cb = errorSpentcpp(st[k], x, bsf, bsfpar);

        for (size_t m = 0; m < M; ++m) {
          a(m, k-1) = (futBounds[k-1] * sqrtI[k-1] - zscaled[m]) / sqrtI2[k-1];
          if (a(m, k-1) > b(m, k-1)) a(m, k-1) = b(m, k-1);
        }

        // lambda expression for finding futility bound at stage k
        // it is an increasing function in aval, and we want to find
        // the root where it crosses 0
        auto g = [&](double aval) -> double {
          for (size_t m = 0; m < M; ++m) {
            a(m, k) = (aval * sqrtI[k] - zscaled[m]) / sqrtI2[k];
            if (a(m, k) > b(m, k)) a(m, k) = b(m, k);
          }
          probs = exitprob_mams_cpp(M, r, theta, true, k + 1, b, a, I2);
          v = probs.get<std::vector<double>>("exitProbLower");
          double cpl = std::accumulate(v.begin(), v.end(), 0.0);
          return cpl - cb;
        };

        double bk = critValues[k];
        eps = g(bk);
        double g_minus8 = g(-8.0);

        if (g_minus8 > 0.0) { // no beta spent at current visit
          futBounds[k] = -8.0;
        } else if (eps > 0.0) {
          auto g_for_brent = [&](double aval)->double {
            if (aval == -8.0) return g_minus8;  // avoid recomputation at 8.0
            if (aval == bk) return eps;         // avoid recomputation at b[k]
            return g(aval);
          };

          futBounds[k] = brent(g_for_brent, -8.0, bk, 1e-6);
        } else if (k < kMax-1) {
          return -1.0; // to decrease beta
        } // else it is the final look, a[k] = b[k], so we need eps = g(bk) = 0
        // since eps < 0, it can be used to decrease beta more accurately
      }
    }

    return eps;
  };

  double v1 = f(0.0001), v2 = f(1.0 - alpha);
  double beta = 0.0;
  if (v1 == -1.0 || (v1 < 0.0 && futBounds[kMax-1] == -8.0)) {
    throw std::invalid_argument("Power must be less than 0.9999");
  } else if (v2 > 0.0) {
    throw std::invalid_argument("Power must be greater than alpha");
  } else {
    auto f_for_brent = [&](double x)->double {
      if (x == 0.0001) return v1;  // avoid recomputation at 0.0001
      if (x == 1.0 - alpha) return v2;  // avoid recomputation at 1.0 - alpha
      return f(x);
    };

    beta = brent(f_for_brent, 0.0001, 1.0 - alpha, 1e-6);
    futBounds[kMax-1] = critValues[kMax-1];
    std::memcpy(a.data_ptr() + (kMax - 1) * M, b.data_ptr() + (kMax - 1) * M,
                M * sizeof(double));
    probs = exitprob_mams_cpp(M, r, theta, true, kMax, b, a, I2);
  }

  ListCpp result;
  result.push_back(1.0 - beta, "power");
  result.push_back(std::move(futBounds), "futilityBounds");
  result.push_back(std::move(probs), "probs");
  return result;
}


ListCpp getDesign_mams_cpp(
    const double beta,
    const double IMax,
    const std::vector<double>& theta,
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<unsigned char>& futilityStopping,
    const FlatMatrix& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& futilityBounds,
    const std::vector<double>& futilityCP,
    const std::vector<double>& futilityTheta,
    const std::string& typeBetaSpending,
    const double parameterBetaSpending,
    const std::vector<double>& userBetaSpending,
    const std::vector<double>& spendingTime) {

  // ----------- Input Validation ----------- //
  if (std::isnan(beta) && std::isnan(IMax))
    throw std::invalid_argument("beta and IMax cannot be missing simultaneously");

  if (!std::isnan(beta) && !std::isnan(IMax))
    throw std::invalid_argument("Only one of beta and IMax should be provided");

  if (!std::isnan(IMax) && IMax <= 0)
    throw std::invalid_argument("IMax must be positive");

  if (M < 1) throw std::invalid_argument("M must be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r must be positive");
  if (!none_na(theta)) throw std::invalid_argument("theta must be provided");
  if (theta.size() != M) throw std::invalid_argument("theta should have length M");
  if (kMax < 1) throw std::invalid_argument("kMax must be at least 1");

  // Alpha and Beta must be within valid ranges
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(beta) && (beta >= 1 - alpha || beta < 0.0001)) {
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");
  }

  std::string unknown = std::isnan(beta) ? "beta" : "IMax";

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (informationRates.size() != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates.back() != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (size_t i = 0; i < kMax; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(kMax);
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping.back() != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues.data);
  bool missingFutilityBounds = !none_na(futilityBounds)
    && !none_na(futilityCP) && !none_na(futilityTheta);

  if (!missingCriticalValues && criticalValues.nrow != kMax) {
    throw std::invalid_argument("Invalid number of rows for criticalValues");
  }
  if (!missingCriticalValues && criticalValues.ncol != M) {
    throw std::invalid_argument("Invalid number of columns for criticalValues");
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "user" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  if (missingCriticalValues && asf == "user") {
    if (!none_na(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be specified");
    if (userAlphaSpending.size() != kMax)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending.back() != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (none_na(futilityBounds)) {
      if (futilityBounds.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityBounds");
      }
    } else if (none_na(futilityCP)) {
      if (futilityCP.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityCP");
      }
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (futilityCP[i] < 0.0 || futilityCP[i] > 1.0) {
          throw std::invalid_argument("futilityCP must lie in [0, 1]");
        }
      }
    } else {
      if (futilityTheta.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityTheta");
      }
    }
  }

  std::string bsf = typeBetaSpending;
  for (char &c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (unknown == "IMax") {
    if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
        bsf == "sfkd" || bsf == "sfhsd" || bsf == "user" || bsf == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpending");
    }
  } else {
    if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
        bsf == "sfkd" || bsf == "sfhsd" || bsf == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpending");
    }
  }

  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument ("parameterBetaSpending must be positive for sfKD");
  }

  if (unknown == "IMax" && bsf == "user") {
    if (!none_na(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be specified");
    if (userBetaSpending.size() != kMax)
      throw std::invalid_argument("Invalid length of userBetaSpending");
    if (userBetaSpending[0] < 0.0)
      throw std::invalid_argument("userBetaSpending must be nonnegative");
    if (any_nonincreasing(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be nondecreasing");
    if (userBetaSpending.back() != beta)
      throw std::invalid_argument("userBetaSpending must end with specified beta");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[kMax-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // ----------- End of Input Validation ----------- //
  ListCpp probs;
  std::vector<double> v;

  // set up efficacy bounds
  // construct by-level bounds
  std::vector<int> level(M * kMax);
  std::vector<int> stage(M * kMax);
  std::vector<double> efficacyBounds(M * kMax);
  for (size_t M1 = M; M1 > 0; --M1) {
    for (size_t k = 0; k < kMax; ++k) {
      size_t row = (M - M1) * kMax + k;
      level[row] = M1;
      stage[row] = k + 1;
    }
  }

  bool correctDim = (criticalValues.nrow == kMax && criticalValues.ncol == M);
  if (missingCriticalValues) {
    std::vector<double> cut(kMax);
    for (size_t M1 = M ; M1 > 0; --M1) {
      bool haybittle = false;
      if (correctDim) {
        cut = flatmatrix_get_column(criticalValues, M - M1);

        if (kMax > 1 && cut.size() == kMax) {
          bool hasNaN = false;
          for (size_t i = 0; i < kMax - 1; ++i) {
            if (std::isnan(cut[i])) { hasNaN = true; break; }
          }
          if (!hasNaN && std::isnan(cut[kMax-1])) haybittle = true;
        }
      }

      if (haybittle) { // Haybittle & Peto
        std::vector<double> zero1(M1, 0.0);
        FlatMatrix b1(M1, kMax);

        for (size_t i = 0; i < kMax - 1; ++i) {
          if (!effStopping[i]) cut[i] = 8.0;
          std::fill_n(b1.data_ptr() + i * M1, M1, cut[i]);
        }

        double* last_col = b1.data_ptr() + (kMax - 1) * M1;
        auto f = [&](double x)->double {
          std::fill_n(last_col, M1, x);
          probs = exitprob_mams_cpp(M1, r, zero1, corr_known, kMax, b1, infoRates);
          v = probs.get<std::vector<double>>("exitProbUpper");
          double cpu = std::accumulate(v.begin(), v.end(), 0.0);
          return cpu - alpha;
        };

        cut[kMax - 1] = brent(f, 0.0, 8.0, 1e-6);
      } else {
        cut = getBound_mams_cpp(
          M1, r, corr_known, kMax, infoRates, alpha, asf,
          parameterAlphaSpending, userAlphaSpending, spendTime, effStopping);
      }

      std::copy_n(cut.data(), kMax, efficacyBounds.data() + (M - M1) * kMax);
    }
  } else {
    efficacyBounds = criticalValues.data; // copy from input matrix
  }

  std::vector<double> critValues(kMax);
  std::memcpy(critValues.data(), efficacyBounds.data(), kMax * sizeof(double));

  std::vector<double> zero(M, 0.0);
  FlatMatrix a(M, kMax); a.fill(-8.0);
  FlatMatrix b(M, kMax);
  for (size_t i = 0; i < kMax; ++i) {
    std::fill_n(b.data_ptr() + i * M, M, critValues[i]);
  }
  probs = exitprob_mams_cpp(M, r, zero, corr_known, kMax, b, infoRates);
  v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(kMax);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent[kMax-1];

  FlatMatrix sigma(M, M);
  sigma.fill(r / (r + 1.0));
  for (size_t m = 0; m < M; ++m) {
    sigma(m, m) = 1.0;
  }

  // set up futility bounds
  std::vector<double> futBounds(kMax, NaN);
  if (kMax > 1) {
    if (missingFutilityBounds && bsf == "none") {
      std::fill_n(futBounds.begin(), kMax - 1, -8.0);
      futBounds[kMax-1] = critValues[kMax-1];
    } else if (!missingFutilityBounds) {
      if (none_na(futilityBounds)) {
        for (size_t i = 0; i < kMax - 1; ++i) {
          if (futilityBounds[i] > critValues[i]) {
            throw std::invalid_argument(
                "futilityBounds must lie below criticalValues");
          }
        }
        std::copy_n(futilityBounds.begin(), kMax-1, futBounds.begin());
        futBounds[kMax-1] = critValues[kMax-1];
      } else if (none_na(futilityCP)) {
        double c2 = critValues[kMax - 1];
        for (size_t i = 0; i < kMax - 1; ++i) {
          double q = qmvnormcpp(1 - futilityCP[i], zero, sigma,
                                1024, 16384, 8, 1e-4, 0.0, 314159);
          futBounds[i] = std::sqrt(infoRates[i]) *
            (c2 - std::sqrt(1 - infoRates[i]) * q);
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityCP values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax-1] = critValues[kMax-1];
      } else if (!std::isnan(IMax)) {
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i] * IMax) * futilityTheta[i];
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityTheta values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax-1] = critValues[kMax-1];
      }
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }


  double IMax1 = IMax;
  std::vector<double> information(kMax);
  if (unknown == "IMax") {
    double maxtheta = *std::max_element(theta.begin(), theta.end());
    std::vector<double> mu0(M);
    std::vector<double> lo(M, -8.0), hi(M, critValues[0]);

    auto f = [&](double x)->double {
      double maxInformation = sq(x / maxtheta);
      for (size_t i = 0; i < kMax; ++i) {
        information[i] = infoRates[i] * maxInformation;
      }
      double sqrtI0 = std::sqrt(information[0]);

      // compute stagewise exit probabilities
      if (!missingFutilityBounds || bsf == "none" || kMax == 1) {
        if (!none_na(futilityBounds) && !none_na(futilityCP)
              && none_na(futilityTheta)) {
          for (size_t i = 0; i < kMax - 1; ++i) {
            futBounds[i] = std::sqrt(information[i]) * futilityTheta[i];
            if (futBounds[i] > critValues[i]) return -1.0; // to decrease drift
          }
          futBounds[kMax-1] = critValues[kMax-1];
        }

        for (size_t i = 0; i < kMax; ++i) {
          std::fill_n(a.data_ptr() + i * M, M, futBounds[i]);
        }

        probs = exitprob_mams_cpp(M, r, theta, true, kMax, b, a, information);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
        return (1.0 - overallReject) - beta;
      } else {
        // initialize futility bound to be updated
        std::fill(futBounds.begin(), futBounds.end(), -8.0);
        a.fill(-8.0);
        double eps = 0.0, cb = 0.0;

        // first stage
        if (futStopping[0]) {
          cb = (bsf == "user") ? userBetaSpending[0] :
          errorSpentcpp(spendTime[0], beta, bsf, parameterBetaSpending);

          for (size_t m = 0; m < M; ++m) {
            mu0[m] = theta[m] * sqrtI0;
          }
          auto q = pmvnormcpp(lo, hi, mu0, sigma,
                              1024, 16384, 8, 1e-4, 0.0, 314159);
          eps = q.prob - cb;
          if (eps < 0.0) return -1.0; // to decrease drift
          futBounds[0] = qmvnormcpp(cb, mu0, sigma,
                                    1024, 16384, 8, 1e-4, 0.0, 314159);
        }

        // subsequent stages
        for (size_t k = 1; k < kMax; ++k) {
          if (futStopping[k]) {
            cb = (bsf == "user") ? userBetaSpending[k] :
            errorSpentcpp(spendTime[k], beta, bsf, parameterBetaSpending);

            std::fill_n(a.data_ptr() + (k - 1) * M, M, futBounds[k - 1]);

            // lambda expression for finding futility bound at stage k
            auto g = [&](double aval)->double {
              std::fill_n(a.data_ptr() + k * M, M, aval);
              probs = exitprob_mams_cpp(M, r, theta, true, k + 1, b, a, information);
              v = probs.get<std::vector<double>>("exitProbLower");
              double cpl = std::accumulate(v.begin(), v.end(), 0.0);
              return cpl - cb;
            };

            double bk = critValues[k];
            eps = g(bk);
            double g_minus8 = g(-8.0);

            if (g_minus8 > 0.0) { // no beta spent at current visit
              futBounds[k] = -8.0;
            } else if (eps > 0.0) {
              auto g_for_brent = [&](double aval)->double {
                if (aval == -8.0) return g_minus8;  // avoid recomputation at 8.0
                if (aval == bk) return eps;         // avoid recomputation at b[k]
                return g(aval);
              };

              futBounds[k] = brent(g_for_brent, -8.0, bk, 1e-6);
            } else if (k < kMax-1) {
              return -1.0;
            }
          }
        }

        return eps;
      }
    };

    double drift = brent(f, 0.001, 8.0, 1e-6);
    IMax1 = sq(drift / maxtheta);
    futBounds[kMax-1] = critValues[kMax-1];
    std::fill_n(a.data_ptr() + (kMax - 1) * M, M, futBounds[kMax - 1]);
    probs = exitprob_mams_cpp(M, r, theta, true, kMax, b, a, information);
  } else {
    for (size_t i = 0; i < kMax; ++i) information[i] = infoRates[i] * IMax1;

    if (!missingFutilityBounds || bsf == "none" || kMax == 1) {
      for (size_t i = 0; i < kMax; ++i) {
        std::fill_n(a.data_ptr() + i * M, M, futBounds[i]);
      }
      probs = exitprob_mams_cpp(M, r, theta, true, kMax, b, a, information);
    } else { // beta-spending
      ListCpp out = getPower_mams(M, r, theta, alpha1, kMax, critValues,
                                  information, bsf, parameterBetaSpending,
                                  spendTime, futStopping, 0.0, zero);
      futBounds = out.get<std::vector<double>>("futilityBounds");
      for (size_t i = 0; i < kMax; ++i) {
        std::fill_n(a.data_ptr() + i * M, M, futBounds[i]);
      }
      probs = out.get_list("probs");
    }
  }


  // output the results
  std::vector<double> efficacyTheta(kMax);
  std::vector<double> futTheta(kMax);
  std::vector<double> efficacyP(kMax);
  std::vector<double> futilityP(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    efficacyTheta[i] = critValues[i] / std::sqrt(information[i]);
    futTheta[i] = futBounds[i] / std::sqrt(information[i]);
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
    futilityP[i] = 1.0 - boost_pnorm(futBounds[i]);
  }

  // stagewise exit probabilities under H1
  auto pu = probs.get<std::vector<double>>("exitProbUpper");
  auto pl = probs.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpu(kMax), cpl(kMax);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());
  double overallReject = cpu.back();
  std::vector<double> ptotal(kMax);
  for (size_t i = 0; i < kMax; ++i) ptotal[i] = pu[i] + pl[i];

  double multiplier = (M * r + 1.0) / (r + 1.0);
  std::vector<double> informationOverall(kMax);
  for (size_t i = 0; i < kMax; ++i) {
    informationOverall[i] = information[i] * multiplier;
  }
  double IMaxOverall = informationOverall.back();

  double expectedInformationH1 = std::inner_product(
    ptotal.begin(), ptotal.end(), information.begin(), 0.0);
  double expectedInformationOverallH1 = std::inner_product(
    ptotal.begin(), ptotal.end(), informationOverall.begin(), 0.0);

  // stagewise exit probabilities under H0 with binding futility
  ListCpp probsH0 = exitprob_mams_cpp(M, r, zero, true, kMax, b, a, infoRates);
  auto puH0 = probsH0.get<std::vector<double>>("exitProbUpper");
  auto plH0 = probsH0.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpuH0(kMax), cplH0(kMax);
  std::partial_sum(puH0.begin(), puH0.end(), cpuH0.begin());
  std::partial_sum(plH0.begin(), plH0.end(), cplH0.begin());
  double overallRejectH0 = cpuH0[kMax-1];
  std::vector<double> ptotalH0(kMax);
  for (size_t i = 0; i < kMax; ++i) ptotalH0[i] = puH0[i] + plH0[i];
  double expectedInformationH0 = std::inner_product(
    ptotalH0.begin(), ptotalH0.end(), information.begin(), 0.0);
  double expectedInformationOverallH0 = std::inner_product(
    ptotalH0.begin(), ptotalH0.end(), informationOverall.begin(), 0.0);

  for (size_t i = 0; i < kMax; ++i) {
    if (critValues[i] == 8) effStopping[i] = 0;
    if (futBounds[i] == -8) futStopping[i] = 0;
  }

  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha1, "alpha");
  overallResults.push_back(overallRejectH0, "attainedAlpha");
  overallResults.push_back(M, "M");
  overallResults.push_back(r, "r");
  overallResults.push_back(corr_known, "corr_known");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(IMax1, "information");
  overallResults.push_back(expectedInformationH1, "expectedInformationH1");
  overallResults.push_back(expectedInformationH0, "expectedInformationH0");
  overallResults.push_back(IMaxOverall, "informationOverall");
  overallResults.push_back(expectedInformationOverallH1,
                           "expectedInformationOverallH1");
  overallResults.push_back(expectedInformationOverallH0,
                           "expectedInformationOverallH0");

  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(futBounds), "futilityBounds");
  byStageResults.push_back(std::move(pu), "rejectPerStage");
  byStageResults.push_back(std::move(pl), "futilityPerStage");
  byStageResults.push_back(std::move(cpu), "cumulativeRejection");
  byStageResults.push_back(std::move(cpl), "cumulativeFutility");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(efficacyTheta), "efficacyTheta");
  byStageResults.push_back(std::move(futTheta), "futilityTheta");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(futilityP), "futilityP");
  byStageResults.push_back(std::move(information), "information");
  byStageResults.push_back(std::move(informationOverall), "informationOverall");
  byStageResults.push_back(std::move(effStopping), "efficacyStopping");
  byStageResults.push_back(std::move(futStopping), "futilityStopping");
  byStageResults.push_back(std::move(puH0), "rejectPerStageH0");
  byStageResults.push_back(std::move(plH0), "futilityPerStageH0");
  byStageResults.push_back(std::move(cpuH0), "cumulativeRejectionH0");
  byStageResults.push_back(std::move(cplH0), "cumulativeFutilityH0");

  DataFrameCpp byLevelBounds;
  byLevelBounds.push_back(std::move(level), "level");
  byLevelBounds.push_back(std::move(stage), "stage");
  byLevelBounds.push_back(std::move(efficacyBounds), "efficacyBounds");

  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(typeBetaSpending, "typeBetaSpending");
  settings.push_back(parameterBetaSpending, "parameterBetaSpending");
  settings.push_back(userBetaSpending, "userBetaSpending");
  settings.push_back(spendingTime, "spendingTime");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(byLevelBounds), "byLevelBounds");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");
  return result;
}


//' @title Power and Sample Size for a Multi-Arm Multi-Stage Design
//' @description Computes either the maximum information and stopping
//' boundaries for a multi-arm multi-stage design, or
//' the achieved power when the maximum information and stopping boundaries
//' are provided.
//'
//' @param beta Type II error rate. Provide either \code{beta} or \code{IMax};
//'   the other should be missing.
//' @param IMax Maximum information for any active arm versus the common
//'   control. Provide either \code{IMax} or \code{beta}; the other should
//'   be missing.
//' @param theta A vector of length \eqn{M} representing the true treatment
//'   effects for each active arm versus the common control. The global null
//'   is \eqn{\theta_i = 0} for all \eqn{i}, and alternatives are one-sided:
//'   \eqn{\theta_i > 0} for at least one \eqn{i = 1, \ldots, M}.
//' @param M Number of active treatment arms.
//' @param r Randomization ratio of each active arm to the common control.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param kMax Number of sequential looks.
//' @param informationRates A numeric vector of information rates fixed
//'   before the trial. If unspecified, defaults to \eqn{(1:kMax) / kMax}.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @param criticalValues The matrix of by-level upper boundaries on the
//'   max z-test statistic scale for efficacy stopping.
//'   The first column is for level \code{M}, the second column is for
//'   level \code{M - 1}, and so on, with the last column for level 1.
//'   If left unspecified, the critical values will be computed based
//'   on the specified alpha spending function.
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param futilityBounds A numeric vector of length \code{kMax - 1}
//'   specifying the futility boundaries on the max z-test statistic scale
//'   for futility stopping.
//' @param futilityCP A numeric vector of length \code{kMax - 1} specifying
//'   the futility boundaries on the conditional power scale for futility
//'   stopping.
//' @param futilityTheta A numeric vector of length \code{kMax - 1} specifying
//'   the futility boundaries on the parameter scale for futility stopping.
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param spendingTime A numeric vector of length \code{kMax} specifying the
//'   error spending time at each analysis. Values must be strictly increasing
//'   and ends at 1. If omitted, defaults to \code{informationRates}.
//'
//' @return An S3 object of class \code{mams} with the following components:
//'
//' * \code{overallResults}: A data frame containing:
//'     - \code{overallReject}: Overall probability of rejecting the global
//'       null hypothesis.
//'     - \code{alpha}: Overall significance level.
//'     - \code{attainedAlpha}: The attained significance level, which is
//'       different from the overall significance level in the presence of
//'       futility stopping.
//'     - \code{M}: Number of active arms.
//'     - \code{r}: Randomization ratio per active arm versus control.
//'     - \code{corr_known}: Whether the correlation among Wald statistics
//'       was assumed known.
//'     - \code{kMax}: Number of stages.
//'     - \code{information}: Maximum information for any active arm versus
//'       control.
//'     - \code{expectedInformationH1}: The expected information under H1.
//'     - \code{expectedInformationH0}: The expected information under H0.
//'
//' * \code{byStageResults}: A data frame containing:
//'     - \code{informationRates}: Information rates at each analysis.
//'     - \code{efficacyBounds}: Efficacy boundaries on the max Z-scale.
//'     - \code{futilityBounds}: Futility boundaries on the max Z-scale.
//'     - \code{rejectPerStage}: Probability of efficacy stopping at each stage.
//'     - \code{futilityPerStage}: Probability of futility stopping at each stage.
//'     - \code{cumulativeRejection}: Cumulative probability of efficacy stopping.
//'     - \code{cumulativeFutility}: Cumulative probability of futility stopping.
//'     - \code{cumulativeAlphaSpent}: Cumulative alpha spent.
//'     - \code{efficacyTheta}: Efficacy boundaries on the parameter scale.
//'     - \code{futilityTheta}: Futility boundaries on the parameter scale.
//'     - \code{efficacyP}: Efficacy boundaries on the p-value scale.
//'     - \code{futilityP}: Futility boundaries on the p-value scale.
//'     - \code{information}: Cumulative information for any active arm versus
//'       control at each analysis.
//'     - \code{efficacyStopping}: Indicator of whether efficacy stopping
//'       is permitted at each stage.
//'     - \code{futilityStopping}: Indicator of whether futility stopping
//'       is permitted at each stage.
//'     - \code{rejectPerStageH0}: The probability for efficacy stopping
//'       under H0.
//'     - \code{futilityPerStageH0}: The probability for futility stopping
//'       under H0.
//'     - \code{cumulativeRejectionH0}: The cumulative probability for
//'       efficacy stopping under H0.
//'     - \code{cumulativeFutilityH0}: The cumulative probability for
//'       futility stopping under H0.
//'
//' * \code{settings}: A list of input settings:
//'     - \code{typeAlphaSpending}: The type of alpha spending.
//'     - \code{parameterAlphaSpending}: The parameter value for the chosen
//'       alpha spending function.
//'     - \code{userAlphaSpending}: The user-specified alpha spending values.
//'     - \code{typeBetaSpending}: The type of beta spending.
//'     - \code{parameterBetaSpending}: The parameter value for the chosen
//'       beta spending function.
//'     - \code{userBetaSpending}: The user-specified beta spending values.
//'     - \code{spendingTime}: The error-spending time at each analysis.
//'
//' * \code{byLevelBounds}: A data frame containing the efficacy boundaries
//'   for each level of testing (i.e., number of active arms remaining) and
//'   each stage. Columns include:
//'     - \code{level}: Number of active arms remaining (1 to \eqn{M}).
//'     - \code{stage}: Stage index (1 to \code{kMax}).
//'     - \code{efficacyBounds}: Efficacy boundaries on the max Z-scale
//'       for the given level and stage.
//'
//' @details If \code{corr_known} is \code{FALSE}, critical boundaries are
//' computed assuming independence among the Wald statistics in each stage
//' (a conservative assumption). Power calculations, however, use the
//' correlation implied by the randomization ratio \eqn{r}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @examples
//'
//' # Example 1: obtain the maximum information given power
//' (design1 <- getDesign_mams(
//'   beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   kMax = 3, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' # Example 2: obtain power given the maximum information
//' (design2 <- getDesign_mams(
//'   IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
//'   kMax = 3, informationRates = seq(1, 3)/3,
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getDesign_mams(
    const double beta = NA_REAL,
    const double IMax = NA_REAL,
    const Rcpp::NumericVector& theta = NA_REAL,
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::Nullable<Rcpp::NumericMatrix> criticalValues = R_NilValue,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBounds = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityCP = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityTheta = R_NilValue,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const Rcpp::NumericVector& userBetaSpending = NA_REAL,
    const Rcpp::NumericVector& spendingTime = NA_REAL) {

  std::vector<double> thetaVec(theta.begin(), theta.end());
  std::vector<double> infoRates(informationRates.begin(), informationRates.end());
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  std::vector<double> userAlpha(userAlphaSpending.begin(), userAlphaSpending.end());
  std::vector<double> userBeta(userBetaSpending.begin(), userBetaSpending.end());
  std::vector<double> spendTime(spendingTime.begin(), spendingTime.end());

  // Handle optional matrix safely
  FlatMatrix critValues;
  if (criticalValues.isNotNull()) {
    Rcpp::NumericMatrix cm(criticalValues); // unwrap
    critValues = flatmatrix_from_Rmatrix(cm);
  } else {
    critValues = FlatMatrix(1, 1);
    critValues(0, 0) = std::numeric_limits<double>::quiet_NaN(); // placeholder
  }

  std::vector<double> futBounds, futCP, futTheta;
  if (futilityBounds.isNotNull()) {
    futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  } else {
    futBounds = std::vector<double>(1, NaN);
  }

  if (futilityCP.isNotNull()) {
    futCP = Rcpp::as<std::vector<double>>(futilityCP);
  } else {
    futCP = std::vector<double>(1, NaN);
  }

  if (futilityTheta.isNotNull()) {
    futTheta = Rcpp::as<std::vector<double>>(futilityTheta);
  } else {
    futTheta = std::vector<double>(1, NaN);
  }

  auto cpp_result = getDesign_mams_cpp(
    beta, IMax, thetaVec, static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(kMax), infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlpha, futBounds, futCP, futTheta, typeBetaSpending,
    parameterBetaSpending, userBeta, spendTime
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "mams";
  return result;
}


ListCpp adaptDesign_mams_cpp(
    double betaNew,
    double INew,
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t L,
    const std::vector<double>& zL,
    const std::vector<double>& theta,
    const double IMax,
    const size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<unsigned char>& futilityStopping,
    const FlatMatrix& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& futilityBounds,
    const std::vector<double>& futilityCP,
    const std::vector<double>& futilityTheta,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
    const size_t MNew,
    const std::vector<int>& selected,
    const double rNew,
    const size_t kNew,
    const std::vector<double>& informationRatesNew,
    const std::vector<unsigned char>& efficacyStoppingNew,
    const std::vector<unsigned char>& futilityStoppingNew,
    const std::string& typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double>& futilityBoundsInt,
    const std::vector<double>& futilityCPInt,
    const std::vector<double>& futilityThetaInt,
    const std::string& typeBetaSpendingNew,
    const double parameterBetaSpendingNew,
    const std::vector<double>& userBetaSpendingNew,
    const std::vector<double>& spendingTimeNew) {

  // ----------- Start of Input Validation ----------- //
  if (std::isnan(betaNew) && std::isnan(INew))
    throw std::invalid_argument("betaNew and INew cannot be missing simultaneously");
  if (!std::isnan(betaNew) && !std::isnan(INew))
    throw std::invalid_argument("Only one of betaNew and INew should be provided");
  if (!std::isnan(INew) && INew <= 0.0)
    throw std::invalid_argument("INew must be positive");

  if (M < 1) throw std::invalid_argument("M must be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r must be positive");
  if (L < 1) throw std::invalid_argument("L must be at least 1");
  if (!none_na(zL)) throw std::invalid_argument("zL must be provided");
  if (zL.size() != M) throw std::invalid_argument("Invalid length for zL");
  if (!none_na(theta)) throw std::invalid_argument("theta must be provided");
  if (theta.size() != M) throw std::invalid_argument("Invalid length for theta");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  if (!std::isnan(betaNew) && (betaNew < 0.0001 || betaNew >= 1)) {
    throw std::invalid_argument("betaNew must lie in [0.0001, 1)");
  }

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (informationRates.size() != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates.back() != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (size_t i = 0; i < kMax; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(kMax);
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping.back() != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues.data);
  bool missingFutilityBounds = !none_na(futilityBounds)
    && !none_na(futilityCP) && !none_na(futilityTheta);

  if (!missingCriticalValues && criticalValues.nrow != kMax) {
    throw std::invalid_argument("Invalid number of rows for criticalValues");
  }
  if (!missingCriticalValues && criticalValues.ncol != M) {
    throw std::invalid_argument("Invalid number of columns for criticalValues");
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "user" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  if (missingCriticalValues && asf == "user") {
    if (!none_na(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be specified");
    if (userAlphaSpending.size() != kMax)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending.back() != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (none_na(futilityBounds)) {
      if (futilityBounds.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityBounds");
      }
    } else if (none_na(futilityCP)) {
      if (futilityCP.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityCP");
      }
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (futilityCP[i] < 0.0 || futilityCP[i] > 1.0) {
          throw std::invalid_argument("futilityCP must lie in [0, 1]");
        }
      }
    } else {
      if (futilityTheta.size() < kMax - 1) {
        throw std::invalid_argument("Insufficient length for futilityTheta");
      }
    }
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime.back() != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // ----------- New Design Input Validation ----------- //
  std::vector<size_t> selectedNew;
  std::vector<double> infoRatesNew;
  std::vector<unsigned char> effStoppingNew;
  std::vector<unsigned char> futStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  std::string bsfNew = typeBetaSpendingNew;
  for (char &c : bsfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  bool missingFutilityBoundsInt = !none_na(futilityBoundsInt)
    && !none_na(futilityCPInt) && !none_na(futilityThetaInt);

  if (MNew < 1) throw std::invalid_argument("MNew must be at least 1");
  for (auto i : selected) {
    if (i < 1 || i > static_cast<int>(M)) {
      throw std::invalid_argument("Invalid value in selected");
    }
    selectedNew.push_back(static_cast<size_t>(i - 1));
  }
  selectedNew = unique_sorted(selectedNew);
  if (selectedNew.size() != MNew)
    throw std::invalid_argument("Length of selected does not match MNew");

  if (rNew <= 0.0) throw std::invalid_argument("rNew must be positive");

  size_t k1 = kMax - L;
  std::vector<double> s1(k1);
  for (size_t i = 0; i < k1; ++i) {
    s1[i] = (infoRates[L + i] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
  }

  if (!MullerSchafer) {
    infoRatesNew = s1;
    effStoppingNew.resize(k1);
    std::memcpy(effStoppingNew.data(), effStopping.data() + L,
                k1 * sizeof(unsigned char));

    if (none_na(futilityStoppingNew)) {
      if (futilityStoppingNew.size() != k1)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew.back() != 1)
        throw std::invalid_argument("futilityStoppingNew must end with 1");
      futStoppingNew = futilityStoppingNew; // copy
    } else {
      futStoppingNew.resize(k1);
      std::memcpy(futStoppingNew.data(), futStopping.data() + L,
                  k1 * sizeof(unsigned char));
    }

    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() != k1)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew.back() != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = s1;
    }

  } else {
    if (kNew < 1) throw std::invalid_argument("kNew must be at least 1");

    // informationRatesNew: default to (1:kNew)/kNew if missing
    infoRatesNew.resize(kNew);
    if (none_na(informationRatesNew)) {
      if (informationRatesNew.size() != kNew)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew.back() != 1.0)
        throw std::invalid_argument("informationRatesNew must end with 1");
      infoRatesNew = informationRatesNew; // copy
    } else {
      for (size_t i = 0; i < kNew; ++i)
        infoRatesNew[i] = static_cast<double>(i+1) / static_cast<double>(kNew);
    }

    // effStoppingNew: default to all 1s if missing
    if (none_na(efficacyStoppingNew)) {
      if (efficacyStoppingNew.size() != kNew)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew.back() != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew; // copy
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    // futStoppingNew: default to all 1s if missing
    if (none_na(futilityStoppingNew)) {
      if (futilityStoppingNew.size() != kNew)
        throw std::invalid_argument("Invalid length for futilityStoppingNew");
      if (futilityStoppingNew.back() != 1)
        throw std::invalid_argument("futilityStoppingNew must end with 1");
      futStoppingNew = futilityStoppingNew; // copy
    } else {
      futStoppingNew.assign(kNew, 1);
    }

    if (!(asfNew == "of" || asfNew == "sfof" || asfNew == "sfp" ||
        asfNew == "sfkd" || asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }
    if ((asfNew == "sfkd" || asfNew == "sfhsd") &&
        std::isnan(parameterAlphaSpendingNew)) {
      throw std::invalid_argument("Missing value for parameterAlphaSpendingNew");
    }
    if (asfNew == "sfkd" && parameterAlphaSpendingNew <= 0.0) {
      throw std::invalid_argument (
          "parameterAlphaSpendingNew must be positive for sfKD");
    }

    // spendingTimeNew: default to informationRatesNew if missing
    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() != kNew)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew.back() != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = infoRatesNew;
    }
  }

  size_t k2 = MullerSchafer ? kNew : k1;

  if (!missingFutilityBoundsInt) {
    if (none_na(futilityBoundsInt)) {
      if (futilityBoundsInt.size() < k2 - 1) {
        throw std::invalid_argument("Insufficient length for futilityBoundsInt");
      }
    } else if (none_na(futilityCPInt)) {
      if (futilityCPInt.size() < k2 - 1) {
        throw std::invalid_argument("Insufficient length for futilityCPInt");
      }
      for (size_t i = 0; i < k2 - 1; ++i) {
        if (futilityCPInt[i] < 0.0 || futilityCPInt[i] > 1.0) {
          throw std::invalid_argument("futilityCPInt must lie in [0, 1]");
        }
      }
    } else {
      if (futilityThetaInt.size() < k2 - 1) {
        throw std::invalid_argument("Insufficient length for futilityThetaInt");
      }
    }
  }

  if (std::isnan(INew)) {
    if (missingFutilityBoundsInt && !(bsfNew == "sfof" || bsfNew == "sfp" ||
        bsfNew == "sfkd" || bsfNew == "sfhsd" || bsfNew == "user" ||
        bsfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpendingNew");
    }
  } else {
    if (missingFutilityBoundsInt && !(bsfNew == "sfof" || bsfNew == "sfp" ||
        bsfNew == "sfkd" || bsfNew == "sfhsd" || bsfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeBetaSpendingNew");
    }
  }

  if ((bsfNew == "sfkd" || bsfNew == "sfhsd") &&
      std::isnan(parameterBetaSpendingNew)) {
    throw std::invalid_argument("Missing value for parameterBetaSpendingNew");
  }
  if (bsfNew == "sfkd" && parameterBetaSpendingNew <= 0.0) {
    throw std::invalid_argument(
        "parameterBetaSpendingNew must be positive for sfKD");
  }

  if (std::isnan(INew) && bsfNew == "user") {
    if (!none_na(userBetaSpendingNew))
      throw std::invalid_argument("userBetaSpendingNew must be specified");
    if (userBetaSpendingNew.size() != k2)
      throw std::invalid_argument("Invalid length of userBetaSpendingNew");
    if (userBetaSpendingNew[0] < 0.0)
      throw std::invalid_argument("userBetaSpendingNew must be nonnegative");
    if (any_nonincreasing(userBetaSpendingNew))
      throw std::invalid_argument("userBetaSpendingNew must be nondecreasing");
    if (userBetaSpendingNew.back() != betaNew)
      throw std::invalid_argument(
          "userBetaSpendingNew must end with specified betaNew");
  }


  // ----------- End of Input Validation ----------- //
  ListCpp probs;
  std::vector<double> v;

  // obtain critical values for the primary trial
  // by level critical values for the integrated trial
  // construct by-level bounds for the primary trial
  std::vector<int> level1(M * kMax);
  std::vector<int> stage1(M * kMax);
  std::vector<double> efficacyBounds1(M * kMax);
  for (size_t M1 = M; M1 > 0; --M1) {
    for (size_t k = 0; k < kMax; ++k) {
      size_t row = (M - M1) * kMax + k;
      level1[row] = M1;
      stage1[row] = k + 1;
    }
  }

  bool correctDim = (criticalValues.nrow == kMax && criticalValues.ncol == M);
  if (missingCriticalValues) {
    std::vector<double> cut(kMax);
    for (size_t M1 = M ; M1 > 0; --M1) {
      bool haybittle = false;
      if (correctDim) {
        cut = flatmatrix_get_column(criticalValues, M - M1);

        if (kMax > 1 && cut.size() == kMax) {
          bool hasNaN = false;
          for (size_t i = 0; i < kMax - 1; ++i) {
            if (std::isnan(cut[i])) { hasNaN = true; break; }
          }
          if (!hasNaN && std::isnan(cut[kMax-1])) haybittle = true;
        }
      }

      if (haybittle) { // Haybittle & Peto
        std::vector<double> zero1(M1, 0.0);
        FlatMatrix b1(M1, kMax);

        for (size_t i = 0; i < kMax - 1; ++i) {
          if (!effStopping[i]) cut[i] = 8.0;
          std::fill_n(b1.data_ptr() + i * M1, M1, cut[i]);
        }

        double* last_col = b1.data_ptr() + (kMax - 1) * M1;
        auto f = [&](double x)->double {
          std::fill_n(last_col, M1, x);
          probs = exitprob_mams_cpp(M1, r, zero1, corr_known, kMax, b1, infoRates);
          v = probs.get<std::vector<double>>("exitProbUpper");
          double cpu = std::accumulate(v.begin(), v.end(), 0.0);
          return cpu - alpha;
        };

        cut[kMax - 1] = brent(f, 0.0, 8.0, 1e-6);
      } else {
        cut = getBound_mams_cpp(
          M1, r, corr_known, kMax, infoRates, alpha, asf,
          parameterAlphaSpending, userAlphaSpending, spendTime, effStopping);
      }

      std::copy_n(cut.data(), kMax, efficacyBounds1.data() + (M - M1) * kMax);
    }
  } else {
    efficacyBounds1 = criticalValues.data; // copy from input matrix
  }

  std::vector<double> critValues(kMax);
  std::memcpy(critValues.data(), efficacyBounds1.data(), kMax * sizeof(double));

  std::vector<double> zero(M, 0.0);
  FlatMatrix sigma(M, M);
  sigma.fill(r / (1.0 + r));
  for (size_t i = 0; i < M; ++i) sigma(i, i) = 1.0;

  FlatMatrix a(M, kMax); a.fill(-8.0);
  FlatMatrix b(M, kMax);
  for (size_t i = 0; i < kMax; ++i) {
    std::fill_n(b.data_ptr() + i * M, M, critValues[i]);
  }
  probs = exitprob_mams_cpp(M, r, zero, corr_known, kMax, b, infoRates);
  v = probs.get<std::vector<double>>("exitProbUpper");
  double p0 = std::accumulate(v.begin(), v.end(), 0.0);
  double alpha1 = missingCriticalValues ? alpha : p0;

  // obtain futility bounds for the primary trial
  std::vector<double> futBounds(kMax);
  if (kMax > 1) {
    if (missingFutilityBounds) {
      std::fill_n(futBounds.begin(), kMax - 1, -8.0);
      futBounds[kMax-1] = critValues[kMax-1];
    } else if (!missingFutilityBounds) {
      if (none_na(futilityBounds)) {
        for (size_t i = 0; i < kMax - 1; ++i) {
          if (futilityBounds[i] > critValues[i]) {
            throw std::invalid_argument(
                "futilityBounds must lie below criticalValues");
          }
        }
        std::copy_n(futilityBounds.begin(), kMax-1, futBounds.begin());
        futBounds[kMax-1] = critValues[kMax-1];
      } else if (none_na(futilityCP)) {
        double c2 = critValues[kMax - 1];
        for (size_t i = 0; i < kMax - 1; ++i) {
          double q = qmvnormcpp(1 - futilityCP[i], zero, sigma,
                                1024, 16384, 8, 1e-4, 0.0, 314159);
          futBounds[i] = std::sqrt(infoRates[i]) *
            (c2 - std::sqrt(1 - infoRates[i]) * q);
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityCP values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax-1] = critValues[kMax-1];
      } else {
        for (size_t i = 0; i < kMax - 1; ++i) {
          futBounds[i] = std::sqrt(infoRates[i] * IMax) * futilityTheta[i];
          if (futBounds[i] > critValues[i]) {
            throw std::invalid_argument("futilityTheta values are too large to "
                                          "be compatible with criticalValues");
          }
        }
        futBounds[kMax-1] = critValues[kMax-1];
      }
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // information for the primary trial
  std::vector<double> information1(kMax);
  for (size_t i = 0; i < kMax; ++i) information1[i] = infoRates[i] * IMax;

  for (size_t i = 0; i < kMax; ++i) {
    std::fill_n(a.data_ptr() + i * M, M, futBounds[i]);
  }

  // compute conditional alpha, conditional power
  std::vector<double> r1(k1);
  FlatMatrix b1(M, k1);
  FlatMatrix a1(M, k1);
  a1.fill(-8.0);
  for (size_t i = 0; i < k1; ++i) {
    r1[i] = infoRates[L - 1] / infoRates[L + i];
    double* bptr = b1.data_ptr() + i * M; // start of column i of b1
    if (effStopping[L + i]) {
      double cut = critValues[L + i];
      r1[i] = infoRates[L - 1] / infoRates[L + i];
      double sqrt_r1 = std::sqrt(r1[i]);
      double denom = std::sqrt(1.0 - r1[i]);
      // write contiguous column
      for (size_t m = 0; m < M; ++m) {
        bptr[m] = (cut - zL[m] * sqrt_r1) / denom;
      }
    } else {
      std::fill_n(bptr, M, 8.0);
    }
  }

  // secondary trial information for original design
  double INew1 = IMax * (1.0 - infoRates[L - 1]);
  std::vector<double> I1(k1);
  for (size_t i = 0; i < k1; ++i) {
    I1[i] = INew1 * s1[i];
  }

  // conditional type I error
  probs = exitprob_mams_cpp(M, r, zero, corr_known, k1, b1, I1);
  auto v0 = probs.get<std::vector<double>>("exitProbUpper");
  double alphaNew = std::accumulate(v0.begin(), v0.end(), 0.0);

  // conditional power
  for (size_t i = 0; i < k1; ++i) {
    double* aptr = a1.data_ptr() + i * M; // start of column i of a1
    if (futStopping[L + i]) {
      double cut = futBounds[L + i];
      double sqrt_r1 = std::sqrt(r1[i]);
      double denom = std::sqrt(1.0 - r1[i]);
      // write contiguous column
      for (size_t m = 0; m < M; ++m) {
        aptr[m] = (cut - zL[m] * sqrt_r1) / denom;
      }
    } else {
      std::fill_n(aptr, M, -8.0);
    }
  }

  probs = exitprob_mams_cpp(M, r, theta, true, k1, b1, a1, I1);
  auto v1 = probs.get<std::vector<double>>("exitProbUpper");
  double conditionalPower = std::accumulate(v1.begin(), v1.end(), 0.0);

  // secondary trial design
  double IL = information1[L - 1];
  double sqrtIL = std::sqrt(IL);
  std::vector<double> theta2 = subset(theta, selectedNew);
  std::vector<double> zero2(MNew, 0.0);
  FlatMatrix sigma2(MNew, MNew);
  sigma2.fill(rNew / (rNew + 1.0));
  for (size_t i = 0; i < MNew; ++i) sigma2(i, i) = 1.0;

  std::string asf2;
  double asfpar2;
  std::vector<double> cpu0(k2);
  if (!MullerSchafer) {
    asf2 = "user";
    asfpar2 = NaN;
    std::partial_sum(v0.begin(), v0.end(), cpu0.begin());
  } else {
    asf2 = asfNew;
    asfpar2 = parameterAlphaSpendingNew;
    if (asf2 != "none" && asf2 != "of") {
      for (size_t i = 0; i < k2; ++i) {
        cpu0[i] = errorSpentcpp(spendTimeNew[i], alphaNew, asf2, asfpar2);
      }
    }
  }

  std::vector<double> I2(k2); // information levels for secondary trial
  std::vector<double> Ic(k2); // information levels for integrated trial
  std::vector<double> sqrtI2(k2), sqrtIc(k2);
  std::vector<double> critValues2(k2, 8.0); // for integrated trial
  std::vector<double> futBounds2(k2, -8.0); // for integrated trial
  FlatMatrix b2(MNew, k2); b2.fill(8.0);  // for secondary trial
  FlatMatrix a2(MNew, k2); a2.fill(-8.0); // for secondary trial

  std::vector<double> zscaled(MNew);
  for (size_t j = 0; j < MNew; ++j) zscaled[j] = zL[selectedNew[j]] * sqrtIL;

  double maxtheta = *std::max_element(theta.begin(), theta.end());

  if (std::isnan(betaNew)) { // INew is provided, power calculation problem
    for (size_t i = 0; i < k2; ++i) {
      I2[i] = INew * infoRatesNew[i];
      Ic[i] = I2[i] + IL;
      sqrtI2[i] = std::sqrt(I2[i]);
      sqrtIc[i] = std::sqrt(Ic[i]);
    }

    // first obtain the efficacy bounds for the secondary trial
    if (asf2 == "of") {
      auto g = [&b2, &I2, &sqrtI2, &sqrtIc, &zscaled, &zero2,
                &effStoppingNew, &probs, &v,
                k2, alphaNew, MNew, rNew, corr_known]
      (double aval)->double {
        double col_const = aval * sqrtIc[k2 - 1];
        for (size_t i = 0; i < k2; ++i) {
          double* colptr = b2.data_ptr() + i * MNew;
          if (effStoppingNew[i]) {
            double denom = sqrtI2[i];
            for (size_t j = 0; j < MNew; ++j) {
              colptr[j] = (col_const - zscaled[j]) / denom;
            }
          } else {
            std::fill_n(colptr, MNew, 8.0);
          }
        }

        probs = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, b2, I2);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - alphaNew;
      };

      double cof = brent(g, 0.0, 8.0, 1e-6);
      double col_const = cof * sqrtIc[k2 - 1];
      for (size_t i = 0; i < k2; ++i) {
        double* colptr = b2.data_ptr() + i * MNew;
        if (effStoppingNew[i]) {
          critValues2[i] = col_const / sqrtIc[i];
          double denom = sqrtI2[i];
          for (size_t j = 0; j < MNew; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }
        } else {
          critValues2[i] = 8.0;
          std::fill_n(colptr, MNew, 8.0);
        }
      }
    } else if (asf2 == "none") {
      for (size_t i = 0; i < k2 - 1; ++i) critValues2[i] = 8.0;
      double denom = sqrtI2[k2 - 1];

      auto g = [&b2, &I2, &sqrtIc, &zscaled, &zero2, &probs, &v,
                denom, k2, alphaNew, MNew, rNew, corr_known]
      (double aval)->double {
        double* colptr = b2.data_ptr() + (k2 - 1) * MNew;
        double col_const = aval * sqrtIc[k2 - 1];
        for (size_t j = 0; j < MNew; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }

        probs = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, b2, I2);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - alphaNew;
      };

      double cof = brent(g, 0.0, 8.0, 1e-6);
      critValues2[k2 - 1] = cof;
      double* colptr = b2.data_ptr() + (k2 - 1) * MNew;
      double col_const = cof * sqrtIc[k2 - 1];
      for (size_t j = 0; j < MNew; ++j) {
        colptr[j] = (col_const - zscaled[j]) / denom;
      }
    } else {
      for (size_t i = 0; i < k2; ++i) {
        if (!effStoppingNew[i]) continue;
        double denom = sqrtI2[i];

        auto g = [&b2, &I2, &sqrtIc, &zscaled, &cpu0, &zero2, &probs, &v,
                  denom, i, MNew, rNew, corr_known]
        (double aval)->double {
          double col_const = aval * sqrtIc[i];
          double* colptr = b2.data_ptr() + i * MNew;
          // update critical values of the secondary trial at current look
          for (size_t j = 0; j < MNew; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }

          probs = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, i + 1, b2, I2);
          v = probs.get<std::vector<double>>("exitProbUpper");
          double p0 = std::accumulate(v.begin(), v.end(), 0.0);
          return p0 - cpu0[i];
        };

        double cof = brent(g, 0.0, 8.0, 1e-6);
        double col_const = cof * sqrtIc[i];
        critValues2[i] = cof;
        double* colptr = b2.data_ptr() + i * MNew;
        for (size_t j = 0; j < MNew; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }
      }
    }

    // now compute futility bounds if needed
    if (k2 > 1) {
      if (missingFutilityBoundsInt && bsfNew == "none") {
        std::fill_n(futBounds2.begin(), k2 - 1, -8.0);
        futBounds2[k2-1] = critValues2[k2-1];
      } else if (!missingFutilityBoundsInt) {
        if (none_na(futilityBoundsInt)) {
          for (size_t i = 0; i < k2 - 1; ++i) {
            if (futilityBoundsInt[i] > critValues2[i]) {
              throw std::invalid_argument(
                  "futilityBoundsInt must lie below critical values "
                  "for the integrated trial");
            }
          }
          std::copy_n(futilityBoundsInt.begin(), k2-1, futBounds2.begin());
          futBounds2[k2-1] = critValues2[k2-1];
        } else if (none_na(futilityCPInt)) {
          double c2 = critValues2[k2 - 1];
          for (size_t i = 0; i < k2 - 1; ++i) {
            double q = qmvnormcpp(1 - futilityCPInt[i], zero2, sigma2,
                                  1024, 16384, 8, 1e-4, 0.0, 314159);
            double sc = Ic[i] / Ic[k2 - 1];
            futBounds2[i] = std::sqrt(sc) * (c2 - std::sqrt(1 - sc) * q);
            if (futBounds2[i] > critValues2[i]) {
              throw std::invalid_argument(
                  "futilityCPInt values are too large to be compatible with "
                  "critical values for the integrated trial");
            }
          }
          futBounds2[k2-1] = critValues2[k2-1];
        } else {
          for (size_t i = 0; i < k2 - 1; ++i) {
            futBounds2[i] = std::sqrt(Ic[i]) * futilityThetaInt[i];
            if (futBounds2[i] > critValues2[i]) {
              throw std::invalid_argument(
                  "futilityThetaInt values are too large to be compatible with "
                  "critical values for the integrated trial");
            }
          }
          futBounds2[k2-1] = critValues2[k2-1];
        }
      }
    } else {
      if (missingFutilityBoundsInt) {
        futBounds2 = critValues2;
      }
    }

    if (missingFutilityBoundsInt && bsfNew != "none" && k2 > 1) { // beta-spending
      ListCpp out = getPower_mams(
        MNew, rNew, theta2, alphaNew, k2, critValues2, Ic,
        bsfNew, parameterBetaSpendingNew, spendTimeNew, futStoppingNew, IL, zL);
      futBounds2 = out.get<std::vector<double>>("futilityBounds");
    }

    // update the actual futility bounds of the secondary trial
    for (size_t i = 0; i < k2; ++i) {
      for (size_t m = 0; m < MNew; ++m) {
        a2(m, i) = (futBounds2[i] * sqrtIc[i] - zscaled[m]) / sqrtI2[i];
      }
    }

  } else { // unknown information, sample size calculation problem
    // obtain required max information for the secondary trial given target power
    std::vector<double> lo(MNew, -8.0), hi(MNew), mu0(MNew);

    auto f = [&critValues2, &futBounds2, &b2, &a2, &I2, &Ic,
              &sqrtI2, &sqrtIc, &zscaled, &cpu0, &zero2, &theta2, &infoRatesNew,
              &effStoppingNew, &futStoppingNew, &spendTimeNew,
              &lo, &hi, &mu0, &sigma2, &userBetaSpendingNew,
              &futilityBoundsInt, &futilityCPInt, &futilityThetaInt, &probs, &v,
              betaNew, k2, asf2, alphaNew, IL, MNew, rNew, corr_known,
              missingFutilityBoundsInt, bsfNew, parameterBetaSpendingNew, maxtheta]
    (double x)->double {
      double Inew = sq(x / maxtheta);
      for (size_t i = 0; i < k2; ++i) {
        I2[i] = Inew * infoRatesNew[i];
        Ic[i] = I2[i] + IL;
        sqrtI2[i] = std::sqrt(I2[i]);
        sqrtIc[i] = std::sqrt(Ic[i]);
      }

      for (size_t m = 0; m < MNew; ++m) {
        mu0[m] = theta2[m] * sqrtI2[0];
      }

      // first obtain the efficacy bounds for the secondary trial
      if (asf2 == "of") {
        auto g = [&b2, &I2, &sqrtI2, &sqrtIc, &zscaled, &zero2,
                  &effStoppingNew, &probs, &v,
                  k2, alphaNew, MNew, rNew, corr_known]
        (double aval)->double {
          double col_const = aval * sqrtIc[k2 - 1];
          for (size_t i = 0; i < k2; ++i) {
            double* colptr = b2.data_ptr() + i * MNew;
            if (effStoppingNew[i]) {
              double denom = sqrtI2[i];
              for (size_t j = 0; j < MNew; ++j) {
                colptr[j] = (col_const - zscaled[j]) / denom;
                if (colptr[j] < -8.0) colptr[j] = -8.0;
              }
            } else {
              std::fill_n(colptr, MNew, 8.0);
            }
          }

          probs = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, b2, I2);
          v = probs.get<std::vector<double>>("exitProbUpper");
          double p0 = std::accumulate(v.begin(), v.end(), 0.0);
          return p0 - alphaNew;
        };

        double cof = brent(g, 0.0, 8.0, 1e-6);
        double col_const = cof * sqrtIc[k2 - 1];
        for (size_t i = 0; i < k2; ++i) {
          double* colptr = b2.data_ptr() + i * MNew;
          if (effStoppingNew[i]) {
            critValues2[i] = col_const / sqrtIc[i];
            double denom = sqrtI2[i];
            for (size_t j = 0; j < MNew; ++j) {
              colptr[j] = (col_const - zscaled[j]) / denom;
              if (colptr[j] < -8.0) colptr[j] = -8.0;
            }
          } else {
            critValues2[i] = 8.0;
            std::fill_n(colptr, MNew, 8.0);
          }
        }
      } else if (asf2 == "none") {
        for (size_t i = 0; i < k2 - 1; ++i) critValues2[i] = 8.0;
        double denom = sqrtI2[k2 - 1];

        auto g = [&b2, &I2, &sqrtIc, &zscaled, &zero2, &probs, &v,
                  denom, k2, alphaNew, MNew, rNew, corr_known]
        (double aval)->double {
          double* colptr = b2.data_ptr() + (k2 - 1) * MNew;
          double col_const = aval * sqrtIc[k2 - 1];
          for (size_t j = 0; j < MNew; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
            if (colptr[j] < -8.0) colptr[j] = -8.0;
          }

          probs = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, b2, I2);
          v = probs.get<std::vector<double>>("exitProbUpper");
          double p0 = std::accumulate(v.begin(), v.end(), 0.0);
          return p0 - alphaNew;
        };

        double cof = brent(g, 0.0, 8.0, 1e-6);
        critValues2[k2 - 1] = cof;
        double* colptr = b2.data_ptr() + (k2 - 1) * MNew;
        double col_const = cof * sqrtIc[k2 - 1];
        for (size_t j = 0; j < MNew; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }
      } else {
        for (size_t i = 0; i < k2; ++i) {
          if (!effStoppingNew[i]) continue;
          double denom = sqrtI2[i];

          auto g = [&b2, &I2, &sqrtIc, &zscaled, &cpu0, &zero2, &probs, &v,
                    denom, i, MNew, rNew, corr_known]
          (double aval)->double {
            double col_const = aval * sqrtIc[i];
            double* colptr = b2.data_ptr() + i * MNew;
            // update critical values of the secondary trial at current look
            for (size_t j = 0; j < MNew; ++j) {
              colptr[j] = (col_const - zscaled[j]) / denom;
              if (colptr[j] < -8.0) colptr[j] = -8.0;
            }

            probs = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, i + 1, b2, I2);
            v = probs.get<std::vector<double>>("exitProbUpper");
            double p0 = std::accumulate(v.begin(), v.end(), 0.0);
            return p0 - cpu0[i];
          };

          double cof = brent(g, 0.0, 8.0, 1e-6);
          double col_const = cof * sqrtIc[i];
          critValues2[i] = cof;
          double* colptr = b2.data_ptr() + i * MNew;
          for (size_t j = 0; j < MNew; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }
        }
      }

      // now compute futility bounds if needed
      if (k2 > 1) {
        if (missingFutilityBoundsInt && bsfNew == "none") {
          std::fill_n(futBounds2.begin(), k2 - 1, -8.0);
          futBounds2[k2-1] = critValues2[k2-1];
        } else if (!missingFutilityBoundsInt) {
          if (none_na(futilityBoundsInt)) {
            for (size_t i = 0; i < k2 - 1; ++i) {
              if (futilityBoundsInt[i] > critValues2[i]) {
                throw std::invalid_argument(
                    "futilityBoundsInt must lie below critical values "
                    "for the integrated trial");
              }
            }
            std::copy_n(futilityBoundsInt.begin(), k2-1, futBounds2.begin());
            futBounds2[k2-1] = critValues2[k2-1];
          } else if (none_na(futilityCPInt)) {
            double c2 = critValues2[k2 - 1];
            for (size_t i = 0; i < k2 - 1; ++i) {
              double q = qmvnormcpp(1 - futilityCPInt[i], zero2, sigma2,
                                    1024, 16384, 8, 1e-4, 0.0, 314159);
              double sc = Ic[i] / Ic[k2 - 1];
              futBounds2[i] = std::sqrt(sc) * (c2 - std::sqrt(1 - sc) * q);
              if (futBounds2[i] > critValues2[i]) {
                throw std::invalid_argument(
                    "futilityCPInt values are too large to be compatible with "
                    "critical values for the integrated trial");
              }
            }
            futBounds2[k2-1] = critValues2[k2-1];
          } else {
            for (size_t i = 0; i < k2 - 1; ++i) {
              futBounds2[i] = std::sqrt(Ic[i]) * futilityThetaInt[i];
              if (futBounds2[i] > critValues2[i]) {
                throw std::invalid_argument(
                    "futilityThetaInt values are too large to be compatible with "
                    "critical values for the integrated trial");
              }
            }
            futBounds2[k2-1] = critValues2[k2-1];
          }
        }
      } else {
        if (missingFutilityBoundsInt) {
          futBounds2 = critValues2;
        }
      }

      if (!missingFutilityBoundsInt || bsfNew == "none" || k2 == 1) {
        for (size_t i = 0; i < k2; ++i) {
          for (size_t m = 0; m < MNew; ++m) {
            a2(m, i) = (futBounds2[i] * sqrtIc[i] - zscaled[m]) / sqrtI2[i];
            if (a2(m, i) > b2(m, i)) a2(m, i) = b2(m, i);
          }
        }

        probs = exitprob_mams_cpp(MNew, rNew, theta2, true, k2, b2, a2, I2);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
        return (1.0 - overallReject) - betaNew;
      } else {
        // initialize futility bound to be updated
        std::fill(futBounds2.begin(), futBounds2.end(), -8.0);
        double eps = 0.0, cb = 0.0;

        // first stage
        if (futStoppingNew[0]) {
          cb = (bsfNew == "user") ? userBetaSpendingNew[0] :
          errorSpentcpp(spendTimeNew[0], betaNew, bsfNew,
                        parameterBetaSpendingNew);

          auto g = [&](double aval) -> double {
            for (size_t m = 0; m < MNew; ++m) {
              hi[m] = (aval * sqrtIc[0] - zscaled[m]) / sqrtI2[0];
            }

            auto q = pmvnormcpp(lo, hi, mu0, sigma2,
                                1024, 16384, 8, 1e-4, 0.0, 314159);
            return q.prob - cb;
          };

          eps = g(critValues2[0]);
          if (eps < 0.0) return -1.0; // to decrease drift
          futBounds2[0] = brent(g, -8.0, critValues2[0], 1e-6);
        }

        // subsequent stages
        for (size_t k = 1; k < k2; ++k) {
          if (futStoppingNew[k]) {
            cb = (bsfNew == "user") ? userBetaSpendingNew[k] :
            errorSpentcpp(spendTimeNew[k], betaNew, bsfNew,
                          parameterBetaSpendingNew);

            for (size_t m = 0; m < MNew; ++m) {
              a2(m, k-1) = (futBounds2[k-1] *sqrtIc[k-1] - zscaled[m]) /sqrtI2[k-1];
              if (a2(m, k-1) > b2(m, k-1)) a2(m, k-1) = b2(m, k-1);
            }

            // lambda expression for finding futility bound at stage k
            auto g = [&](double aval)->double {
              for (size_t m = 0; m < MNew; ++m) {
                a2(m, k) = (aval * sqrtIc[k] - zscaled[m]) / sqrtI2[k];
                if (a2(m, k) > b2(m, k)) a2(m, k) = b2(m, k);
              }
              probs = exitprob_mams_cpp(MNew, rNew, theta2, true, k + 1, b2, a2, I2);
              v = probs.get<std::vector<double>>("exitProbLower");
              double cpl = std::accumulate(v.begin(), v.end(), 0.0);
              return cpl - cb;
            };

            double bk = critValues2[k];
            eps = g(bk);
            double g_minus8 = g(-8.0);
            if (g_minus8 > 0.0) { // no beta spent at current visit
              futBounds2[k] = -8.0;
            } else if (eps > 0.0) {
              auto g_for_brent = [&](double aval)->double {
                if (aval == -8.0) return g_minus8;  // avoid recomputation at 8.0
                if (aval == bk) return eps;         // avoid recomputation at b[k]
                return g(aval);
              };

              futBounds2[k] = brent(g_for_brent, -8.0, bk, 1e-6);
            } else if (k < k2-1) {
              return -1.0;
            }
          }
        }

        return eps;
      }
    };

    double drift = brent(f, 0.001, 8.0, 1e-6);
    INew = sq(drift / maxtheta);

    futBounds2[k2-1] = critValues2[k2-1];
    std::memcpy(a2.data_ptr() + (k2-1) * MNew, b2.data_ptr() + (k2-1) * MNew,
                MNew * sizeof(double));
  }

  // compute conditional power of the secondary trial
  probs = exitprob_mams_cpp(MNew, rNew, theta2, true, k2, b2, a2, I2);
  v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cpu1(k2);
  std::partial_sum(v.begin(), v.end(), cpu1.begin());
  double p2 = cpu1.back();

  std::vector<int> hypothesis2(MNew * k2);
  std::vector<int> stage2(MNew * k2);
  std::vector<double> effBounds2_long(MNew * k2);
  std::vector<double> futBounds2_long(k2);
  for (size_t j = 0; j < MNew; ++j) {
    for (size_t i = 0; i < k2; ++i) {
      size_t row = j * k2 + i;
      hypothesis2[row] = selectedNew[j] + 1;
      stage2[row] = i + 1;
      effBounds2_long[row] = b2(j, i);
      futBounds2_long[row] = a2(j, i);
    }
  }

  // integrated trial information and critical values
  size_t kc = L + k2;
  std::vector<double> Ic_full(kc);
  std::copy_n(information1.data(), L, Ic_full.data());
  std::copy_n(Ic.data(), k2, Ic_full.data() + L);

  double IMaxc = Ic_full.back();
  std::vector<double> infoRates_full(kc);
  for (size_t i = 0; i < kc; ++i) infoRates_full[i] = Ic_full[i] / IMaxc;

  std::vector<double> critValues_full(kc);
  std::copy_n(critValues.data(), L, critValues_full.data());
  std::copy_n(critValues2.data(), k2, critValues_full.data() + L);

  std::vector<double> futBounds_full(kc);
  std::copy_n(futBounds.data(), L, futBounds_full.data());
  std::copy_n(futBounds2.data(), k2, futBounds_full.data() + L);

  // critical values for each intersection hypothesis
  size_t ntests = (1 << MNew) - 1; // 2^MNew - 1
  std::vector<std::string> intersectHyp(ntests * kc);
  std::vector<int> stage(ntests * kc);
  std::vector<double> efficacyBounds(ntests * kc);

  std::vector<int> idx_in_selected(M, -1);
  for (size_t t = 0; t < selectedNew.size(); ++t)
    idx_in_selected[selectedNew[t]] = static_cast<int>(t);

  for (size_t h = 0; h < ntests; ++h) {
    size_t number = ntests - h;

    std::vector<size_t> primary;
    std::vector<size_t> selectedNew2;
    primary.reserve(M); selectedNew2.reserve(MNew);

    // Count M2 by popcount
    size_t M2 = 0;
    for (size_t b = 0; b < MNew; ++b) if ((number >> (MNew - 1 - b)) & 1u) ++M2;
    size_t M1 = M2 + (M - MNew);

    for (size_t j = 0; j < M; ++j) {
      int pos = idx_in_selected[j];
      if (pos >= 0) {
        // arm j is among selectedNew; check corresponding bit in 'number'
        if ((number >> (MNew - 1 - pos)) & 1u) {
          primary.push_back(j);
          selectedNew2.push_back(j);
        }
      } else {
        // non-selected arm always included
        primary.push_back(j);
      }
    }

    std::vector<double> critValues(kMax);
    std::copy_n(efficacyBounds1.data() + (M - M1) * kMax, kMax, critValues.data());

    // compute conditional alpha for MAMS with the given M1 hypotheses
    FlatMatrix b1(M1, k1); b1.fill(8.0);
    for (size_t i = 0; i < k1; ++i) {
      double col_const = critValues[L + i];
      double r1 = infoRates[L - 1] / infoRates[L + i];
      double sqrt_r1 = std::sqrt(r1);
      double denom = std::sqrt(1.0 - r1);
      double* colptr = b1.data_ptr() + i * M1; // start of column i
      if (effStopping[L + i]) {
        for (size_t j = 0; j < M1; ++j) {
          colptr[j] = (col_const - zL[primary[j]] * sqrt_r1) / denom;
        }
      } else {
        std::fill_n(colptr, M1, 8.0);
      }
    }

    // conditional type I error
    std::vector<double> zero1(M1, 0.0);
    probs = exitprob_mams_cpp(M1, r, zero1, corr_known, k1, b1, I1);
    v0 = probs.get<std::vector<double>>("exitProbUpper");
    double alphaNew = std::accumulate(v0.begin(), v0.end(), 0.0);

    std::vector<double> cpu0(k2);
    if (!MullerSchafer) {
      std::partial_sum(v0.begin(), v0.end(), cpu0.begin());
    } else {
      if (asf2 != "none" && asf2 != "of") {
        for (size_t i = 0; i < k2; ++i) {
          cpu0[i] = errorSpentcpp(spendTimeNew[i], alphaNew, asf2, asfpar2);
        }
      }
    }

    std::vector<double> zero2(M2, 0.0);
    std::vector<double> critValues2(k2, 8.0);
    FlatMatrix b2(M2, k2); b2.fill(8.0);

    std::vector<double> zscaled(M2);
    for (size_t j = 0; j < M2; ++j) zscaled[j] = zL[selectedNew2[j]] * sqrtIL;

    if (asf2 == "of") {
      auto g = [&b2, &I2, &sqrtI2, &sqrtIc, &zscaled, &zero2,
                &effStoppingNew, &probs, &v,
                k2, alphaNew, M2, rNew, corr_known]
      (double aval)->double {
        double col_const = aval * sqrtIc[k2 - 1];
        for (size_t i = 0; i < k2; ++i) {
          double* colptr = b2.data_ptr() + i * M2;
          if (effStoppingNew[i]) {
            double denom = sqrtI2[i];
            for (size_t j = 0; j < M2; ++j) {
              colptr[j] = (col_const - zscaled[j]) / denom;
            }
          } else {
            std::fill_n(colptr, M2, 8.0);
          }
        }

        probs = exitprob_mams_cpp(M2, rNew, zero2, corr_known, k2, b2, I2);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - alphaNew;
      };

      double cof = brent(g, 0.0, 8.0, 1e-6);
      double col_const = cof * sqrtIc[k2 - 1];
      for (size_t i = 0; i < k2; ++i) {
        if (effStoppingNew[i]) {
          critValues2[i] = col_const / sqrtIc[i];
        } else {
          critValues2[i] = 8.0;
        }
      }
    } else if (asf2 == "none") {
      for (size_t i = 0; i < k2 - 1; ++i) critValues2[i] = 8.0;
      double denom = sqrtI2[k2 - 1];

      auto g = [&b2, &I2, &sqrtIc, &zscaled, &zero2, &probs, &v,
                denom, k2, alphaNew, M2, rNew, corr_known]
      (double aval)->double {
        double col_const = aval * sqrtIc[k2 - 1];
        double* colptr = b2.data_ptr() + (k2 - 1) * M2;
        for (size_t j = 0; j < M2; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }

        probs = exitprob_mams_cpp(M2, rNew, zero2, corr_known, k2, b2, I2);
        v = probs.get<std::vector<double>>("exitProbUpper");
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - alphaNew;
      };

      double cof = brent(g, 0.0, 8.0, 1e-6);
      critValues2[k2 - 1] = cof;
    } else {
      for (size_t i = 0; i < k2; ++i) {
        if (!effStoppingNew[i]) continue;
        double denom = sqrtI2[i];

        auto g = [&b2, &I2, &sqrtIc, &zscaled, &cpu0, &zero2, &probs, &v,
                  denom, i, M2, rNew, corr_known]
        (double aval)->double {
          double col_const = aval * sqrtIc[i];
          double* colptr = b2.data_ptr() + i * M2;
          // update critical values of the secondary trial at current look
          for (size_t j = 0; j < M2; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }

          probs = exitprob_mams_cpp(M2, rNew, zero2, corr_known, i+1, b2, I2);
          v = probs.get<std::vector<double>>("exitProbUpper");
          double p0 = std::accumulate(v.begin(), v.end(), 0.0);
          return p0 - cpu0[i];
        };

        critValues2[i] = brent(g, 0.0, 8.0, 1e-6);
      }
    }

    std::vector<double> critValues_full(kc);
    std::copy_n(critValues.data(), L, critValues_full.data());
    std::copy_n(critValues2.data(), k2, critValues_full.data() + L);
    for (size_t i = 0; i < kc; ++i) {
      size_t row = h * kc + i;

      std::string s;
      s.reserve(3 + selectedNew2.size() * 4); // rough reservation
      for (size_t t = 0; t < selectedNew2.size(); ++t) {
        if (t) { s.append(", "); }
        s.append(std::to_string(selectedNew2[t] + 1));
      }
      intersectHyp[row] = std::move(s);

      stage[row] = i + 1;
      efficacyBounds[row] = critValues_full[i];
    }
  }

  ListCpp des1;
  des1.push_back(M, "M");
  des1.push_back(r, "r");
  des1.push_back(corr_known, "corr_known");
  des1.push_back(L, "L");
  des1.push_back(zL, "zL");
  des1.push_back(theta, "theta");
  des1.push_back(IMax, "maxInformation");
  des1.push_back(kMax, "kMax");
  des1.push_back(std::move(infoRates), "informationRates");
  des1.push_back(std::move(critValues), "efficacyBounds");
  des1.push_back(std::move(futBounds), "futilityBounds");
  des1.push_back(std::move(information1), "information");
  des1.push_back(alpha1, "alpha");
  des1.push_back(alphaNew, "conditionalAlpha");
  des1.push_back(conditionalPower, "conditionalPower");
  des1.push_back(MullerSchafer, "MullerSchafer");

  DataFrameCpp byLevelBounds1;
  byLevelBounds1.push_back(std::move(level1), "level");
  byLevelBounds1.push_back(std::move(stage1), "stage");
  byLevelBounds1.push_back(std::move(efficacyBounds1), "efficacyBounds");
  des1.push_back(std::move(byLevelBounds1), "byLevelBounds");

  ListCpp des2;
  des2.push_back(p2, "overallReject");
  des2.push_back(alphaNew, "alpha");
  des2.push_back(MNew, "M");
  des2.push_back(rNew, "r");
  des2.push_back(selected, "selected");
  des2.push_back(corr_known, "corr_known");
  des2.push_back(k2, "kMax");
  des2.push_back(INew, "maxInformation");
  des2.push_back(std::move(infoRatesNew), "informationRates");
  des2.push_back(std::move(cpu1), "cumulativeRejection");
  des2.push_back(std::move(cpu0), "cumulativeAlphaSpent");
  des2.push_back(std::move(I2), "information");
  des2.push_back(asf2, "typeAlphaSpending");
  des2.push_back(asfpar2, "parameterAlphaSpending");
  des2.push_back(bsfNew, "typeBetaSpending");
  des2.push_back(parameterBetaSpendingNew, "parameterBetaSpending");
  des2.push_back(userBetaSpendingNew, "userBetaSpending");
  des2.push_back(spendTimeNew, "spendingTime");

  DataFrameCpp byHypBounds;
  byHypBounds.push_back(std::move(hypothesis2), "hypothesis");
  byHypBounds.push_back(std::move(stage2), "stage");
  byHypBounds.push_back(std::move(effBounds2_long), "efficacyBounds");
  byHypBounds.push_back(std::move(futBounds2_long), "futilityBounds");
  des2.push_back(std::move(byHypBounds), "byHypothesisBounds");

  ListCpp des3;
  des3.push_back(M, "M");
  des3.push_back(r, "r");
  des3.push_back(corr_known, "corr_known");
  des3.push_back(MNew, "MNew");
  des3.push_back(rNew, "rNew");
  des3.push_back(selected, "selected");
  des3.push_back(L, "L");
  des3.push_back(zL, "zL");
  des3.push_back(theta, "theta");
  des3.push_back(IMaxc, "maxInformation");
  des3.push_back(kc, "kMax");
  des3.push_back(std::move(infoRates_full), "informationRates");
  des3.push_back(std::move(critValues_full), "efficacyBounds");
  des3.push_back(std::move(futBounds_full), "futilityBounds");
  des3.push_back(std::move(Ic_full), "information");

  DataFrameCpp byInterBounds;
  byInterBounds.push_back(std::move(intersectHyp), "intersectionHypothesis");
  byInterBounds.push_back(std::move(stage), "stage");
  byInterBounds.push_back(std::move(efficacyBounds), "efficacyBounds");
  des3.push_back(std::move(byInterBounds), "byIntersectionBounds");

  ListCpp result;
  result.push_back(std::move(des1), "primaryTrial");
  result.push_back(std::move(des2), "secondaryTrial");
  result.push_back(std::move(des3), "integratedTrial");
  return result;
}


//' @title Adaptive Multi-Arm Multi-Stage Design
//' @description
//' Calculates the conditional power for specified incremental
//' information, given the interim results, parameter value,
//' data-dependent changes in treatment selection,
//' the error spending function, and
//' the number and spacing of interim looks. Conversely,
//' calculates the incremental information required to attain
//' a specified conditional power, given the interim results,
//' parameter value, data-dependent changes in treatment selection,
//' the error spending function, and the number and spacing of interim looks.
//'
//' @param betaNew The type II error for the secondary trial.
//' @param INew The maximum information for any active arm versus the common
//'   control in the secondary trial. Either
//'   \code{betaNew} or \code{INew} should be provided, while the other
//'   must be missing.
//' @param M Number of active treatment arms in the primary trial.
//' @param r Randomization ratio of each active arm to the common control
//'   in the primary trial.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is assumed.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistics at the interim adaptation look of
//'   the primary trial.
//' @param theta A vector of length \eqn{M} representing the assumed treatment
//'   effects for each active arm versus the common control. The global null
//'   is \eqn{\theta_i = 0} for all \eqn{i}, and alternatives are one-sided:
//'   \eqn{\theta_i > 0} for at least one \eqn{i = 1, \ldots, M}.
//' @param IMax Maximum information for any active arm versus the common
//'   control for the primary trial. Must be provided.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param futilityStopping Indicators of whether futility stopping is
//'   allowed at each stage of the primary trial.
//' @param criticalValues The matrix of by-level upper boundaries on the
//'   max z-test statistic scale for efficacy stopping for the primary trial.
//'   The first column is for level \code{M}, the second column is for
//'   level \code{M - 1}, and so on, with the last column for level 1.
//'   If left unspecified, the critical values will be computed based
//'   on the specified alpha spending function.
//' @param alpha The significance level of the primary trial.
//'   Defaults to 0.025.
//' @param typeAlphaSpending The type of alpha spending for the primary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function,
//'   \code{"user"} for user defined spending, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param userAlphaSpending The user-defined alpha spending for the
//'   primary trial. Represents the cumulative alpha spent up to each stage.
//' @param futilityBounds The futility boundaries on the max-z statistic
//'   scale for the primary trial. Defaults to \code{rep(-8, kMax-1)}
//'   if left unspecified.
//' @param futilityCP The conditional power-based futility bounds for the
//'   primary trial.
//' @param futilityTheta The parameter value-based futility bounds for the
//'   primary trial.
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param MNew Number of active treatment arms in the secondary trial.
//' @param selected The indices of the selected active treatment arms for
//'   the secondary trial.
//' @param rNew Randomization ratio of each active arm to the common control
//'   in the secondary trial.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param futilityStoppingNew The indicators of whether futility stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param typeAlphaSpendingNew The type of alpha spending for the secondary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param futilityBoundsInt The futility boundaries on the max-z statistic
//'   scale for new stages of the integrated trial.
//' @param futilityCPInt The conditional power-based futility bounds for
//'   new stages of the integrated trial.
//' @param futilityThetaInt The parameter value-based futility bounds for the
//'   new stages of the integrated trial.
//' @param typeBetaSpendingNew The type of beta spending for the secondary
//'   trial. One of the following:
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function,
//'   \code{"user"} for user defined spending, and
//'   \code{"none"} for no early futility stopping.
//'   Defaults to \code{"none"}.
//' @param parameterBetaSpendingNew The parameter value of beta spending
//'   for the secondary trial. Corresponds to \eqn{\rho} for \code{"sfKD"},
//'   and \eqn{\gamma} for \code{"sfHSD"}.
//' @param userBetaSpendingNew The user-defined cumulative beta spending.
//'   Represents the cumulative beta spent up to each stage of the
//'   secondary trial.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRatesNew}.
//'
//' @return An \code{adaptDesign_mams} object with three list components:
//'
//' * \code{primaryTrial}: A list of selected information for the primary
//'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{L},
//'   \code{zL}, \code{theta}, \code{maxInformation}, \code{kMax},
//'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
//'   \code{information}, \code{alpha}, \code{conditionalAlpha},
//'   \code{conditionalPower}, \code{MullerSchafer}, and \code{byLevelBounds},
//'   where \code{byLevelBounds} is a data frame with columns \code{level},
//'   \code{stage}, and \code{efficacyBounds}, representing the efficacy
//'   bounds for each combination of the number of active arms and
//'   the stage of analysis in the primary trial.
//'
//' * \code{secondaryTrial}: A list of selected information for the secondary
//'   trial, including \code{overallReject}, \code{alpha}, \code{M}, \code{r},
//'   \code{selected}, \code{corr_known}, \code{kMax}, \code{maxInformation},
//'   \code{informationRates}, \code{cumulativeRejection},
//'   \code{cumulativeAlphaSpent}, \code{information},
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{typeBetaSpending}, \code{parameterBetaSpending},
//'   \code{userBetaSpending}, \code{spendingTime}, and
//'   \code{byHypothesisBounds}, where \code{byHypothesisBounds} is a
//'   data frame with columns \code{hypothesis}, \code{stage},
//'   \code{efficacyBounds}, and \code{futilityBounds}, representing
//'   the efficacy and futility bounds for each hypothesis and each
//'   stage of analysis in the secondary trial.
//'
//' * \code{integratedTrial}: A list of selected information for the integrated
//'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{MNew},
//'   \code{rNew}, \code{selected}, \code{L}, \code{zL}, \code{theta},
//'   \code{maxInformation}, \code{kMax}, \code{informationRates},
//'   \code{efficacyBounds}, \code{futilityBounds}, \code{information},
//'   and \code{byIntersectionBounds}, where \code{byIntersectionBounds} is
//'   a data frame with columns \code{intersectionHypothesis}, \code{stage},
//'   and \code{efficacyBounds}, representing the efficacy bounds for
//'   each intersection hypothesis and each stage of analysis in the
//'   integrated trial.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @seealso \code{\link{getDesign_mams}}
//'
//' @examples
//'
//' # Two active treatment arms are compared with a common control in a
//' # two-look time-to-event design using O'Brien–Fleming–type alpha spending.
//' # Suppose each active arm has a true hazard ratio of 0.75 versus control,
//' # and the total number of events across all three arms at the final analysis
//' # is 486. This corresponds to approximately 324 events for each active arm
//' # versus the common control. Under these assumptions, the trial has about
//' # 80% power to detect the treatment effect in at least one active arm.
//'
//' (des1 <- getDesign_mams(
//'   IMax = 324 / 4, theta = c(-log(0.75), -log(0.75)),
//'   M = 2, r = 1, kMax = 2, informationRates = c(1/2, 1),
//'   alpha = 0.025, typeAlphaSpending = "OF"))
//'
//' # Now assume that, at the interim analysis, the observed hazard ratios for
//' # the two active arms versus control are 0.91 and 0.78, respectively. Using
//' # the rule “drop any arm with an observed hazard ratio > 0.9”, arm 1 is
//' # dropped. We then aim to achieve 80% conditional power to detect a hazard
//' # ratio of 0.78 for the remaining arm at the final look. The analysis below
//' # indicates that the required total number of events for arm 2 versus control
//' # at the final analysis should be increased from 324 to 535.
//'
//' (des2 <- adaptDesign_mams(
//'   betaNew = 0.2, M = 2, r = 1, corr_known = FALSE,
//'   L = 1, zL = c(-log(0.91), -log(0.78)) * sqrt(324 / 4 / 2),
//'   theta = c(-log(0.91), -log(0.78)),
//'   IMax = 324 / 4, kMax = 2, informationRates = c(1/2, 1),
//'   alpha = 0.025, typeAlphaSpending = "OF",
//'   MNew = 1, selected = 2, rNew = 1))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List adaptDesign_mams(
    double betaNew = NA_REAL,
    double INew = NA_REAL,
    const int M = NA_INTEGER,
    const double r = 1,
    const bool corr_known = true,
    const int L = NA_INTEGER,
    const Rcpp::NumericVector& zL = NA_REAL,
    const Rcpp::NumericVector& theta = NA_REAL,
    const double IMax = NA_REAL,
    const int kMax = NA_INTEGER,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::Nullable<Rcpp::NumericMatrix> criticalValues = R_NilValue,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBounds = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityCP = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityTheta = R_NilValue,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool MullerSchafer = false,
    const int MNew = NA_INTEGER,
    const Rcpp::IntegerVector& selected = NA_INTEGER,
    const double rNew = 1,
    const int kNew = NA_INTEGER,
    const Rcpp::NumericVector& informationRatesNew = NA_REAL,
    const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStoppingNew = NA_LOGICAL,
    const std::string& typeAlphaSpendingNew = "sfOF",
    const double parameterAlphaSpendingNew = NA_REAL,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBoundsInt = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityCPInt = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityThetaInt = R_NilValue,
    const std::string& typeBetaSpendingNew = "none",
    const double parameterBetaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& userBetaSpendingNew = NA_REAL,
    const Rcpp::NumericVector& spendingTimeNew = NA_REAL) {

  auto zLVec = Rcpp::as<std::vector<double>>(zL);
  auto thetaVec = Rcpp::as<std::vector<double>>(theta);
  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto selectedNew = Rcpp::as<std::vector<int>>(selected);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto futStoppingNew = convertLogicalVector(futilityStoppingNew);
  auto userBetaNew = Rcpp::as<std::vector<double>>(userBetaSpendingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  // Handle optional matrix safely
  FlatMatrix critValues;
  if (criticalValues.isNotNull()) {
    Rcpp::NumericMatrix cm(criticalValues); // unwrap
    critValues = flatmatrix_from_Rmatrix(cm);
  } else {
    critValues = FlatMatrix(1, 1);
    critValues(0, 0) = std::numeric_limits<double>::quiet_NaN(); // placeholder
  }

  std::vector<double> futBounds, futCP, futTheta;
  if (futilityBounds.isNotNull()) {
    futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  } else {
    futBounds = std::vector<double>(1, NaN);
  }

  if (futilityCP.isNotNull()) {
    futCP = Rcpp::as<std::vector<double>>(futilityCP);
  } else {
    futCP = std::vector<double>(1, NaN);
  }

  if (futilityTheta.isNotNull()) {
    futTheta = Rcpp::as<std::vector<double>>(futilityTheta);
  } else {
    futTheta = std::vector<double>(1, NaN);
  }

  std::vector<double> futBoundsInt, futCPInt, futThetaInt;
  if (futilityBoundsInt.isNotNull()) {
    futBoundsInt = Rcpp::as<std::vector<double>>(futilityBoundsInt);
  } else {
    futBoundsInt = std::vector<double>(1, NaN);
  }

  if (futilityCPInt.isNotNull()) {
    futCPInt = Rcpp::as<std::vector<double>>(futilityCPInt);
  } else {
    futCPInt = std::vector<double>(1, NaN);
  }

  if (futilityThetaInt.isNotNull()) {
    futThetaInt = Rcpp::as<std::vector<double>>(futilityThetaInt);
  } else {
    futThetaInt = std::vector<double>(1, NaN);
  }

  auto cpp_result = adaptDesign_mams_cpp(
    betaNew, INew,  static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(L), zLVec, thetaVec, IMax,
    static_cast<size_t>(kMax), infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending, parameterAlphaSpending,
    userAlpha, futBounds, futCP, futTheta, spendTime,
    MullerSchafer, static_cast<size_t>(MNew), selectedNew, rNew,
    static_cast<size_t>(kNew), infoRatesNew,
    effStoppingNew, futStoppingNew,
    typeAlphaSpendingNew, parameterAlphaSpendingNew,
    futBoundsInt, futCPInt, futThetaInt, typeBetaSpendingNew,
    parameterBetaSpendingNew, userBetaNew, spendTimeNew
  );

  Rcpp::List result = Rcpp::wrap(cpp_result);
  result.attr("class") = "adaptDesign_mams";
  return result;
}

