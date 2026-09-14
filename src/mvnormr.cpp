#include "mvnormr.h"
#include "dataframe_list.h"
#include "utilities.h"

#include <Rcpp.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using std::size_t;

GHaltonScramble make_scramble(size_t J, std::uint64_t seed_rep,
                              const unsigned int *primes_,
                              const unsigned int *permTN2_) {

  if (J > ghaltonMaxDim)
    throw std::invalid_argument("J exceeds ghaltonMaxDim");

  boost::random::mt19937_64 rng(seed_rep);

  GHaltonScramble S;
  S.J = J;
  S.base.resize(J);
  S.f.resize(J);
  S.inv_base.resize(J);
  S.shift.resize(J * 32);

  for (size_t j = 0; j < J; ++j) {
    const unsigned int b = primes_[j];
    S.base[j] = b;
    S.f[j] = permTN2_[j];
    S.inv_base[j] = 1.0 / static_cast<double>(b);

    boost::random::uniform_int_distribution<unsigned int> Udigit(0u, b - 1u);
    const size_t off = j * 32;
    for (size_t k = 0; k < 32; ++k) {
      S.shift[off + k] = Udigit(rng);
    }
  }

  return S;
}

void ghalton_incremental_init(GHaltonIncremental &G, const GHaltonScramble &S,
                              size_t begin) {

  const size_t dim = (S.J > 0 ? (S.J - 1) : 0);
  G.dim = dim;
  G.i = begin;

  G.digits.assign(dim, std::array<unsigned int, 32>{});
  G.u.assign(dim, 0.0);
  G.invpow.resize(dim * 32);

  // Precompute inv powers per dimension
  for (size_t j = 0; j < dim; ++j) {
    double p = S.inv_base[j]; // 1/b
    const size_t off = j * 32;
    for (size_t k = 0; k < 32; ++k) {
      G.invpow[off + k] = p; // 1/b^(k+1)
      p *= S.inv_base[j];
    }
  }

  // Initialize digits at index 'begin' and compute u[j]
  for (size_t j = 0; j < dim; ++j) {
    const unsigned int b = S.base[j];

    std::uint64_t tmp = begin;
    size_t k = 0;
    while (tmp > 0u && k < 32) {
      G.digits[j][k] = static_cast<unsigned int>(tmp % b);
      tmp /= b;
      ++k;
    }

    double uj = 0.0;
    for (size_t kk = 0; kk < 32; ++kk) {
      const unsigned int dig = scrambled_digit(S, j, kk, G.digits[j][kk]);
      uj += static_cast<double>(dig) * G.invpow[j * 32 + kk];
    }
    G.u[j] = uj;
  }
}

void ghalton_incremental_next(GHaltonIncremental &G, const GHaltonScramble &S) {
  // advance from i to i+1
  G.i++;

  for (size_t j = 0; j < G.dim; ++j) {
    const unsigned int b = S.base[j];
    const size_t off = j * 32;

    // carry propagation in base b
    size_t k = 0;
    while (k < 32) {
      const unsigned int old_digit = G.digits[j][k];
      const unsigned int new_digit = old_digit + 1;

      if (new_digit < b) {
        const unsigned int old_s = scrambled_digit(S, j, k, old_digit);
        const unsigned int new_s = scrambled_digit(S, j, k, new_digit);

        G.digits[j][k] = new_digit;
        G.u[j] += (static_cast<double>(new_s) - static_cast<double>(old_s)) *
                  G.invpow[off + k];
        break;
      } else {
        // wrap to 0 and carry
        const unsigned int old_s = scrambled_digit(S, j, k, old_digit);
        const unsigned int new_s = scrambled_digit(S, j, k, 0u);

        G.digits[j][k] = 0u;
        G.u[j] += (static_cast<double>(new_s) - static_cast<double>(old_s)) *
                  G.invpow[off + k];
        ++k;
      }
    }
  }
}

// sequential conditioning transform contribution for one sample
double mvn_genz_contribution(size_t J, const std::vector<double> &lower_std,
                             const std::vector<double> &upper_std,
                             const std::vector<double> &C, const double *w,
                             double *y) {

  size_t start = 0;
  double p = 1.0;

  for (size_t j = 0; j < J; ++j) {
    double m = 0.0;
    for (size_t k = 0; k < j; ++k) {
      m += C[start + k] * y[k];
    }
    start += j;

    const double a = lower_std[j] - m;
    const double b = upper_std[j] - m;

    double d = pnorm_fast(a);
    double e = pnorm_fast(b);

    const double emd = e - d;
    if (emd <= 0.0)
      return 0.0;
    p *= emd;

    if (j < J - 1) {
      double pp = d + w[j] * emd;
      pp = std::min(std::max(pp, MIN_PROB), MAX_PROB);
      y[j] = qnorm_acklam(pp);
    }
  }

  return p;
}

// Sum over i in [begin, end) for one replicate using incremental Halton state.
double pmvn_sum_range_incremental(ReplicateAccumulator &acc, size_t J,
                                  const std::vector<double> &lower_std,
                                  const std::vector<double> &upper_std,
                                  const std::vector<double> &C, size_t begin,
                                  size_t end, double *w_buf, double *y_buf) {

  const size_t dim = (J > 0 ? (J - 1) : 0);
  double sum = 0.0;

  if (!acc.gen_initialized || acc.G.i != begin) {
    ghalton_incremental_init(acc.G, acc.S, begin);
    acc.gen_initialized = true;
  }

  for (size_t i = begin; i < end; ++i) {
    // copy current u -> w
    for (size_t j = 0; j < dim; ++j)
      w_buf[j] = acc.G.u[j];

    sum += mvn_genz_contribution(J, lower_std, upper_std, C, w_buf, y_buf);

    if (i + 1 < end)
      ghalton_incremental_next(acc.G, acc.S);
  }

  return sum;
}

// Extend one replicate to n_target
void extend_replicate(ReplicateAccumulator &acc, size_t n_target, size_t J,
                      const std::vector<double> &lower_std,
                      const std::vector<double> &upper_std,
                      const std::vector<double> &C, double *w_buf,
                      double *y_buf) {

  if (n_target <= acc.n_done)
    return;

  acc.sum += pmvn_sum_range_incremental(acc, J, lower_std, upper_std, C,
                                        acc.n_done, n_target, w_buf, y_buf);
  acc.n_done = n_target;
}

struct ExtendReplicatesWorker : public RcppParallel::Worker {
  std::vector<ReplicateAccumulator> &reps;
  std::size_t n_target;
  std::size_t J;
  const std::vector<double> &lower_std;
  const std::vector<double> &upper_std;
  const std::vector<double> &C;
  std::vector<double> &phat;

  ExtendReplicatesWorker(std::vector<ReplicateAccumulator> &reps_,
                         std::size_t n_target_, std::size_t J_,
                         const std::vector<double> &lower_std_,
                         const std::vector<double> &upper_std_,
                         const std::vector<double> &C_,
                         std::vector<double> &phat_)
      : reps(reps_), n_target(n_target_), J(J_), lower_std(lower_std_),
        upper_std(upper_std_), C(C_), phat(phat_) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> w(J > 0 ? (J - 1) : 0);
    std::vector<double> z(J, 0.0);

    for (std::size_t rep = begin; rep < end; ++rep) {
      auto &acc = reps[rep];
      extend_replicate(acc, n_target, J, lower_std, upper_std, C, w.data(),
                       z.data());
      phat[rep] = acc.sum / static_cast<double>(n_target);
    }
  }
};

// Compute mean and SE from replicate means phat[0..R-1]
void mean_se_from_replicates(const std::vector<double> &phat, double &mean_out,
                             double &se_out) {
  const size_t R = phat.size();
  mean_sd(phat.data(), R, mean_out, se_out);
  se_out /= std::sqrt(static_cast<double>(R)); // SE of mean
}

// adaptive QMC main loop: extend replicates until error estimate meets
// tolerance or n_max reached. Returns final estimate, error, and n used.
PMVNResult pmvnorm_adaptive(size_t J, const std::vector<double> &lower_std,
                            const std::vector<double> &upper_std,
                            const std::vector<double> &C, size_t n0,
                            size_t n_max, size_t R, double abseps,
                            double releps, std::uint64_t seed, bool parallel) {

  const std::uint64_t base_seed = (seed == 0) ? seed_from_clock_u64() : seed;

  // Initialize replicate accumulators (scramble fixed per replicate)
  std::vector<ReplicateAccumulator> reps(R);
  for (size_t rep = 0; rep < R; ++rep) {
    const std::uint64_t seed_rep = splitmix64(base_seed ^ (rep + 1));
    reps[rep].S = make_scramble(J, seed_rep, primes, permTN2);
    reps[rep].sum = 0.0;
    reps[rep].n_done = 0;
    reps[rep].gen_initialized = false;
  }

  // Scratch buffers reused across all reps/iterations
  std::vector<double> w(J > 0 ? (J - 1) : 0);
  std::vector<double> y(J, 0.0);

  size_t n = n0;
  std::vector<double> phat(R);

  double pbar = std::numeric_limits<double>::quiet_NaN();
  double se = std::numeric_limits<double>::quiet_NaN();
  double error = std::numeric_limits<double>::quiet_NaN();
  const double tcrit = boost_qt(0.975, static_cast<double>(R - 1));

  if (!parallel) {
    for (;;) {
      // Extend each replicate to n
      for (size_t rep = 0; rep < R; ++rep) {
        extend_replicate(reps[rep], n, J, lower_std, upper_std, C, w.data(),
                         y.data());
        phat[rep] = reps[rep].sum / static_cast<double>(n);
      }

      mean_se_from_replicates(phat, pbar, se);
      double tol = std::max(abseps, releps * std::abs(pbar));

      // stopping rule
      error = tcrit * se;
      if (error < tol)
        break;

      // increase n for next iteration, but cap at n_max
      if (n >= n_max)
        break;
      const size_t n2 = n + n0;
      n = (n2 > n_max) ? n_max : n2;
    }
  } else {
    for (;;) {
      ExtendReplicatesWorker worker(reps, n, J, lower_std, upper_std, C, phat);
      RcppParallel::parallelFor(0, R, worker);

      mean_se_from_replicates(phat, pbar, se);
      double tol = std::max(abseps, releps * std::abs(pbar));

      // stopping rule
      error = tcrit * se;
      if (error < tol)
        break;

      // increase n for next iteration, but cap at n_max
      if (n >= n_max)
        break;
      const size_t n2 = n + n0;
      n = (n2 > n_max) ? n_max : n2;
    }
  }

  PMVNResult out;
  out.prob = pbar;
  out.method = "qmc";
  out.error = error;
  out.nsamples = n * R;
  return out;
}

PMVNCholFactor precompute_chol(const FlatMatrix &sigma) {
  size_t J = sigma.nrow;
  if (sigma.ncol != J)
    throw std::invalid_argument("sigma must be square");

  FlatMatrix sigma_copy = sigma; // copy to modify
  int rank = cholesky2(sigma_copy, J, 1e-12);

  // Handle positive semidefinite covariance matrices by adding a tiny
  // diagonal ridge until a stable full-rank factorization is obtained.
  // This keeps the original fast path unchanged for positive-definite sigma.
  if (rank < static_cast<int>(J)) {
    double max_diag = 0.0;
    for (size_t j = 0; j < J; ++j) {
      const double d = sigma(j, j);
      if (!std::isfinite(d) || d < 0.0) {
        throw std::invalid_argument(
            "sigma must have finite, non-negative diagonal");
      }
      if (d > max_diag)
        max_diag = d;
    }

    const double scale = std::max(1.0, max_diag);
    double ridge = 1e-12 * scale;
    bool ok = false;

    for (size_t attempt = 0; attempt < 12; ++attempt) {
      sigma_copy = sigma;
      for (size_t j = 0; j < J; ++j) {
        sigma_copy(j, j) += ridge;
      }

      rank = cholesky2(sigma_copy, J, 1e-12);
      if (rank == static_cast<int>(J)) {
        ok = true;
        break;
      }

      ridge *= 10.0;
    }

    if (!ok) {
      throw std::invalid_argument("sigma must be positive semidefinite (failed "
                                  "to stabilize factorization)");
    }
  }

  std::vector<double> sd(J);
  for (size_t j = 0; j < J; ++j) {
    sd[j] = std::sqrt(sigma_copy(j, j));
  }

  // inverse of D^(-1/2) * L * D^(1/2) in row-major order (lower tri part only)
  std::vector<double> C(J * (J - 1) / 2);
  size_t start = 0;
  for (size_t j = 0; j < J; ++j) {
    double denom = sd[j];
    for (size_t k = 0; k < j; ++k) {
      C[start + k] = sigma_copy(j, k) * sd[k] / denom;
    }
    start += j;
  }

  return PMVNCholFactor{J, std::move(sd), std::move(C)};
}

PMVNResult pmvnorm_with_chol(const PMVNCholFactor &chol,
                             const std::vector<double> &lower,
                             const std::vector<double> &upper,
                             const std::vector<double> &mean, size_t n0,
                             size_t n_max, size_t R, double abseps,
                             double releps, uint64_t seed, bool parallel) {

  const size_t J = chol.J;
  if (lower.size() != J || upper.size() != J || mean.size() != J)
    throw std::invalid_argument("lower/upper/mean must have length J");

  // Early exit: zero probability when any lower >= upper
  for (size_t j = 0; j < J; ++j) {
    if (lower[j] >= upper[j])
      return PMVNResult{0.0, "analytic", 0.0, 1};
  }

  // Standardise bounds using the cached standard deviations (O(J))
  std::vector<double> lower_std(J), upper_std(J);
  for (size_t j = 0; j < J; ++j) {
    lower_std[j] = (lower[j] - mean[j]) / chol.sd[j];
    upper_std[j] = (upper[j] - mean[j]) / chol.sd[j];
  }

  return pmvnorm_adaptive(J, lower_std, upper_std, chol.C, n0, n_max, R, abseps,
                          releps, seed, parallel);
}

// detect if sigma is compound symmetry with non-negative correlations
bool is_compound_symmetry(const FlatMatrix &sigma) {
  size_t J = sigma.nrow;
  if (sigma.ncol != J)
    return false;
  if (J == 0)
    throw std::invalid_argument("sigma must be non-empty");
  if (J == 1)
    return true;

  double diagonal = sigma(0, 0);
  double off_diag = sigma(0, 1);
  if (off_diag < 0.0 || off_diag >= diagonal) {
    return false;
  }

  for (size_t j = 0; j < J; ++j) {
    if (sigma(j, j) != diagonal) {
      return false;
    }
    for (size_t k = j + 1; k < J; ++k) {
      if (sigma(j, k) != off_diag) {
        return false;
      }
    }
  }

  return true;
}

bool get_nonnegative_compound_symmetry(const FlatMatrix &sigma,
                                       double &diagonal, double &off_diag) {
  size_t J = sigma.nrow;
  if (sigma.ncol != J)
    return false;
  if (J == 0)
    throw std::invalid_argument("sigma must be non-empty");

  diagonal = sigma(0, 0);
  off_diag = (J > 1 ? sigma(0, 1) : 0.0);
  if (diagonal <= 0.0 || off_diag < 0.0 || off_diag > diagonal) {
    return false;
  }

  for (size_t j = 0; j < J; ++j) {
    if (sigma(j, j) != diagonal) {
      return false;
    }
    for (size_t k = j + 1; k < J; ++k) {
      if (sigma(j, k) != off_diag) {
        return false;
      }
    }
  }

  return true;
}

// Return Pr(X >= q), where X is the sum of independent Bernoulli variables
// with success probabilities p[0], ..., p[J - 1].
double poisson_binomial_tail(const std::vector<double> &p, size_t q) {
  const size_t J = p.size();
  if (q == 0)
    return 1.0;
  if (q > J)
    return 0.0;

  std::vector<double> mass(J + 1, 0.0);
  mass[0] = 1.0;
  size_t n_done = 0;
  for (double pj : p) {
    // Update backward so mass[k - 1] still refers to the distribution before
    // adding the current Bernoulli variable.
    for (size_t k = n_done + 1; k > 0; --k) {
      mass[k] = mass[k] * (1.0 - pj) + mass[k - 1] * pj;
    }
    mass[0] *= (1.0 - pj);
    ++n_done;
  }

  double tail = 0.0;
  for (size_t k = q; k <= J; ++k) {
    tail += mass[k];
  }
  return std::min(std::max(tail, 0.0), 1.0);
}

// Return Pr(X >= q), where X follows Binomial(n, p).
double binomial_tail(size_t n, size_t q, double p) {
  if (q == 0)
    return 1.0;
  if (q > n || p <= 0.0)
    return 0.0;
  if (p >= 1.0)
    return 1.0;

  double tail = 0.0;
  const double log_p = std::log(p);
  const double log1m_p = std::log1p(-p);
  for (size_t k = q; k <= n; ++k) {
    double log_mass = std::lgamma(static_cast<double>(n) + 1.0) -
                      std::lgamma(static_cast<double>(k) + 1.0) -
                      std::lgamma(static_cast<double>(n - k) + 1.0) +
                      static_cast<double>(k) * log_p +
                      static_cast<double>(n - k) * log1m_p;
    tail += std::exp(log_mass);
  }
  return std::min(std::max(tail, 0.0), 1.0);
}

bool has_equal_entries(const std::vector<double> &x) {
  for (size_t j = 1; j < x.size(); ++j) {
    if (x[j] != x[0]) {
      return false;
    }
  }
  return true;
}

double pordmvnorm_compound_symmetry(double z, size_t r,
                                    const std::vector<double> &mean,
                                    double diagonal, double off_diag,
                                    bool lower_tail) {
  const size_t J = mean.size();
  const size_t q = lower_tail ? r : J - r + 1;
  const bool equal_means = has_equal_entries(mean);

  if (J == 1) {
    return boost_pnorm(z, mean[0], std::sqrt(diagonal), lower_tail);
  }

  if (off_diag == diagonal) {
    if (equal_means) {
      return boost_pnorm(z, mean[0], std::sqrt(off_diag), lower_tail);
    }
    std::vector<double> thresholds(J);
    for (size_t j = 0; j < J; ++j) {
      thresholds[j] = z - mean[j];
    }
    std::sort(thresholds.begin(), thresholds.end());
    return boost_pnorm(thresholds[J - r], 0.0, std::sqrt(off_diag),
                       lower_tail);
  }

  const double sigma_e = std::sqrt(diagonal - off_diag);
  std::vector<double> p(J);
  if (off_diag == 0.0) {
    if (equal_means) {
      double p0 = boost_pnorm(z, mean[0], sigma_e, lower_tail);
      return binomial_tail(J, q, p0);
    }
    for (size_t j = 0; j < J; ++j) {
      p[j] = boost_pnorm(z, mean[j], sigma_e, lower_tail);
    }
    return poisson_binomial_tail(p, q);
  }

  const double sigma_b = std::sqrt(off_diag);
  if (equal_means) {
    auto f = [&](double b) {
      double p0 = boost_pnorm(z, mean[0] + b, sigma_e, lower_tail);
      return binomial_tail(J, q, p0) * boost_dnorm(b, 0.0, sigma_b);
    };

    std::vector<double> breaks = {-8.0 * sigma_b, 0.0, 8.0 * sigma_b};
    return integrate3(f, breaks, 1e-6);
  }

  auto f = [&](double b) {
    for (size_t j = 0; j < J; ++j) {
      p[j] = boost_pnorm(z, mean[j] + b, sigma_e, lower_tail);
    }
    return poisson_binomial_tail(p, q) * boost_dnorm(b, 0.0, sigma_b);
  };

  std::vector<double> breaks = {-8.0 * sigma_b, 0.0, 8.0 * sigma_b};
  return integrate3(f, breaks, 1e-6);
}

double inclusion_exclusion_coef(size_t j, size_t q) {
  double coef = std::exp(std::lgamma(static_cast<double>(j)) -
                         std::lgamma(static_cast<double>(q)) -
                         std::lgamma(static_cast<double>(j - q + 1)));
  return ((j - q) % 2 == 0) ? coef : -coef;
}

// Advance idx to the next k-combination chosen from {0, ..., n - 1}, where
// k = idx.size() and idx is stored in increasing lexicographic order.
// Returns false when idx is already the final combination.
bool next_combination(std::vector<size_t> &idx, size_t n) {
  const size_t k = idx.size();
  for (size_t i = k; i > 0; --i) {
    const size_t pos = i - 1;
    if (idx[pos] < n - k + pos) {
      // Increment the rightmost index that can move, then reset the suffix to
      // the smallest increasing continuation.
      ++idx[pos];
      for (size_t j = pos + 1; j < k; ++j) {
        idx[j] = idx[j - 1] + 1;
      }
      return true;
    }
  }
  return false;
}

PMVNResult pordmvnormcpp(double z, size_t r, const std::vector<double> &mean,
                         const FlatMatrix &sigma, bool lower_tail, size_t n0,
                         size_t n_max, size_t R, double abseps, double releps,
                         uint64_t seed, bool parallel) {
  const size_t J = mean.size();
  if (J < 1)
    throw std::invalid_argument("J must be >= 1");
  if (r < 1 || r > J)
    throw std::invalid_argument("r must be in {1, ..., J}");
  if (sigma.nrow != J || sigma.ncol != J)
    throw std::invalid_argument("sigma must be J x J");
  if (R < 2)
    throw std::invalid_argument("R must be >= 2");
  if (n0 < 1)
    throw std::invalid_argument("n0 must be >= 1");
  if (n_max < n0)
    throw std::invalid_argument("n_max must be >= n0");
  if (abseps <= 0.0)
    throw std::invalid_argument("abseps must be positive");
  if (releps < 0.0)
    throw std::invalid_argument("releps must be non-negative");

  double diagonal = 0.0;
  double off_diag = 0.0;
  if (get_nonnegative_compound_symmetry(sigma, diagonal, off_diag)) {
    double prob = pordmvnorm_compound_symmetry(z, r, mean, diagonal, off_diag,
                                               lower_tail);
    return PMVNResult{prob, "analytic", 0.0, 1};
  }

  const size_t q = lower_tail ? r : J - r + 1;
  double prob = 0.0;
  double compensation = 0.0;
  double error = 0.0;
  size_t nsamples = 0;
  bool all_analytic = true;
  size_t subset_counter = 0;

  for (size_t j = q; j <= J; ++j) {
    std::vector<size_t> idx(j);
    std::iota(idx.begin(), idx.end(), 0);
    const double coef = inclusion_exclusion_coef(j, q);

    do {
      std::vector<double> lower(
        j, lower_tail ? -std::numeric_limits<double>::infinity() : z);
      std::vector<double> upper(
        j, lower_tail ? z : std::numeric_limits<double>::infinity());
      std::vector<double> mean_j(j);
      FlatMatrix sigma_j(j, j);
      for (size_t a = 0; a < j; ++a) {
        mean_j[a] = mean[idx[a]];
        for (size_t b = 0; b < j; ++b) {
          sigma_j(a, b) = sigma(idx[a], idx[b]);
        }
      }

      PMVNResult term =
          pmvnormcpp(lower, upper, mean_j, sigma_j, n0, n_max, R, abseps,
                     releps, seed + subset_counter, parallel);
      const double y = coef * term.prob - compensation;
      const double t = prob + y;
      compensation = (t - prob) - y;
      prob = t;
      error += std::fabs(coef) * term.error;
      nsamples += term.nsamples;
      all_analytic = all_analytic && term.method == "analytic";
      ++subset_counter;
    } while (next_combination(idx, J));
  }

  prob = std::min(std::max(prob, 0.0), 1.0);
  return PMVNResult{prob, all_analytic ? "analytic" : "qmc", error, nsamples};
}

// Main entry point: validate inputs, permute/standardize, factorize, and
// call adaptive routine.
PMVNResult pmvnormcpp(const std::vector<double> &lower,
                      const std::vector<double> &upper,
                      const std::vector<double> &mean, const FlatMatrix &sigma,
                      size_t n0, size_t n_max, size_t R, double abseps,
                      double releps, uint64_t seed, bool parallel) {

  size_t J = lower.size();
  if (J < 1)
    throw std::invalid_argument("J must be >= 1");
  if (upper.size() != J || mean.size() != J)
    throw std::invalid_argument("lower/upper/mean must have same length");
  if (sigma.nrow != J || sigma.ncol != J)
    throw std::invalid_argument("sigma must be J x J");
  if (R < 2)
    throw std::invalid_argument("R must be >= 2");
  if (n0 < 1)
    throw std::invalid_argument("n0 must be >= 1");
  if (n_max < n0)
    throw std::invalid_argument("n_max must be >= n0");
  if (abseps <= 0.0)
    throw std::invalid_argument("abseps must be positive");
  if (releps < 0.0)
    throw std::invalid_argument("releps must be non-negative");

  // handle special case of univariate MVN probability with simple formula
  if (J == 1) {
    double sd = std::sqrt(sigma(0, 0));
    double aj = (lower[0] - mean[0]) / sd;
    double bj = (upper[0] - mean[0]) / sd;
    double p = boost_pnorm(bj) - boost_pnorm(aj);
    if (p <= 0.0)
      return PMVNResult{0.0, "analytic", 0.0, 1};
    return PMVNResult{p, "analytic", 0.0, 1};
  }

  // special case of compound symmetry with positive correlations
  if (J > 1 && is_compound_symmetry(sigma)) {
    double diagonal = sigma(0, 0);
    double off_diag = sigma(0, 1);
    double sigma_e = std::sqrt(diagonal - off_diag);
    double sigma_b = std::sqrt(off_diag);
    if (off_diag == 0.0) {
      // indepedent case, just product of univariate probabilities
      double p = 1.0;
      for (size_t j = 0; j < J; ++j) {
        double aj = (lower[j] - mean[j]) / sigma_e;
        double bj = (upper[j] - mean[j]) / sigma_e;
        double mass = boost_pnorm(bj) - boost_pnorm(aj);
        if (mass <= 0.0)
          return PMVNResult{0.0, "analytic", 0.0, 1};
        p *= mass;
      }
      return PMVNResult{p, "analytic", 0.0, 1};
    }

    auto f = [&](double b) {
      double p = 1.0;
      for (size_t j = 0; j < J; ++j) {
        double aj = (lower[j] - mean[j] - b) / sigma_e;
        double bj = (upper[j] - mean[j] - b) / sigma_e;
        double mass = boost_pnorm(bj) - boost_pnorm(aj);
        if (mass <= 0.0)
          return 0.0; // handle zero prob case
        p *= mass;
      }
      return p * boost_dnorm(b, 0.0, sigma_b);
    };

    std::vector<double> breaks = {-8.0 * sigma_b, 0.0, 8.0 * sigma_b};
    double p = integrate3(f, breaks, 1e-6);
    return PMVNResult{p, "analytic", 0.0, 1};
  }

  return pmvnorm_with_chol(precompute_chol(sigma), lower, upper, mean, n0,
                           n_max, R, abseps, releps, seed, parallel);
}

// [[Rcpp::export]]
Rcpp::List pmvnormRcpp(const std::vector<double> &lower,
                       const std::vector<double> &upper,
                       const std::vector<double> &mean,
                       const Rcpp::NumericMatrix &sigma, size_t n0 = 1024,
                       size_t n_max = 16384, size_t R = 8, double abseps = 1e-4,
                       double releps = 0.0, uint64_t seed = 314159,
                       bool parallel = true) {

  auto sigma_fm = flatmatrix_from_Rmatrix(sigma);
  auto out = pmvnormcpp(lower, upper, mean, sigma_fm, n0, n_max, R, abseps,
                        releps, seed, parallel);
  return Rcpp::List::create(
      Rcpp::Named("prob") = out.prob, Rcpp::Named("method") = out.method,
      Rcpp::Named("error") = out.error, Rcpp::Named("nsamples") = out.nsamples);
}

// [[Rcpp::export]]
Rcpp::List pordmvnormRcpp(double z, size_t r, const std::vector<double> &mean,
              const Rcpp::NumericMatrix &sigma,
              bool lower_tail = true, size_t n0 = 1024,
              size_t n_max = 16384, size_t R = 8,
              double abseps = 1e-4, double releps = 0.0,
              uint64_t seed = 314159, bool parallel = true) {

  auto sigma_fm = flatmatrix_from_Rmatrix(sigma);
  auto out = pordmvnormcpp(z, r, mean, sigma_fm, lower_tail, n0, n_max, R,
               abseps, releps, seed, parallel);
  return Rcpp::List::create(
      Rcpp::Named("prob") = out.prob, Rcpp::Named("method") = out.method,
      Rcpp::Named("error") = out.error, Rcpp::Named("nsamples") = out.nsamples);
}

double qmvnormcpp(const double p, const std::vector<double> &mean,
                  const FlatMatrix &sigma, size_t n0, size_t n_max, size_t R,
                  double abseps, double releps, uint64_t seed, bool parallel) {

  if (p <= 0.0 || p >= 1.0)
    throw std::invalid_argument("p must lie in (0,1)");
  size_t J = mean.size();
  if (J < 1)
    throw std::invalid_argument("J must be >= 1");
  if (sigma.nrow != J || sigma.ncol != J)
    throw std::invalid_argument("sigma must be J x J");
  if (R < 2)
    throw std::invalid_argument("R must be >= 2");
  if (n0 < 1)
    throw std::invalid_argument("n0 must be >= 1");
  if (n_max < n0)
    throw std::invalid_argument("n_max must be >= n0");
  if (abseps <= 0.0)
    throw std::invalid_argument("abseps must be positive");
  if (releps < 0.0)
    throw std::invalid_argument("releps must be non-negative");

  double upper0 = mean_kahan(mean);
  std::vector<double> lower(J, -std::numeric_limits<double>::infinity());
  std::vector<double> upper(J, upper0);
  double x1 = mean[0] - 10.0 * std::sqrt(sigma(0, 0));
  double x2 = mean[0] + 10.0 * std::sqrt(sigma(0, 0));
  if (J == 1 || (J > 1 && is_compound_symmetry(sigma))) {
    auto f = [&](double x) {
      std::fill(upper.begin(), upper.end(), x);
      auto res = pmvnormcpp(lower, upper, mean, sigma, n0, n_max, R, abseps,
                            releps, seed, parallel);
      return res.prob - p;
    };
    return brent(f, x1, x2, 1e-6);
  } else {
    PMVNCholFactor chol = precompute_chol(sigma);
    std::vector<double> lower(J, -std::numeric_limits<double>::infinity());
    std::vector<double> upper(J);
    auto f = [&](double x) {
      std::fill(upper.begin(), upper.end(), x);
      auto res = pmvnorm_with_chol(chol, lower, upper, mean, n0, n_max, R,
                                   abseps, releps, seed, parallel);
      return res.prob - p;
    };
    return brent(f, x1, x2, 1e-6);
  }
}

// [[Rcpp::export]]
double qmvnormRcpp(const double p, const std::vector<double> &mean,
                   const Rcpp::NumericMatrix &sigma, size_t n0 = 1024,
                   size_t n_max = 16384, size_t R = 8, double abseps = 1e-4,
                   double releps = 0.0, uint64_t seed = 314159,
                   bool parallel = true) {
  auto sigma_fm = flatmatrix_from_Rmatrix(sigma);
  return qmvnormcpp(p, mean, sigma_fm, n0, n_max, R, abseps, releps, seed,
                    parallel);
}
