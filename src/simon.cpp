#include "utilities.h"
#include "dataframe_list.h"

#include <Rcpp.h>
#include <boost/math/distributions/binomial.hpp>
#include <boost/random.hpp>

#include <algorithm>     // any_of, distance, fill, min, min_element
#include <cctype>        // tolower
#include <cmath>         // fabs, isnan
#include <cstddef>       // size_t
#include <limits>        // numeric_limits
#include <numeric>       // accumulate
#include <stdexcept>     // invalid_argument
#include <string>        // string
#include <vector>        // vector
#include <unordered_map> // unordered_map


// Hash function for caching pfutile results
struct PfutileKey {
  double p;
  int n1, n2, r1, r;

  bool operator==(const PfutileKey& other) const {
    return n1 == other.n1 && n2 == other.n2 &&
      r1 == other.r1 && r == other.r &&
      std::abs(p - other.p) < 1e-10;
  }
};

namespace std {
template<>
struct hash<PfutileKey> {
  size_t operator()(const PfutileKey& k) const {
    size_t seed = 0;

    // Hash combine pattern (better distribution)
    auto hash_combine = [](size_t& seed, size_t hash) {
      seed ^= hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };

    hash_combine(seed, hash<double>()(k.p));    // Include p!
    hash_combine(seed, hash<int>()(k.n1));
    hash_combine(seed, hash<int>()(k.n2));
    hash_combine(seed, hash<int>()(k.r1));
    hash_combine(seed, hash<int>()(k.r));

    return seed;
  }
};
}

// Optimized pfutile with reusable distribution objects
class PfutileCalculator {
private:
  std::unordered_map<PfutileKey, double> cache_;

public:
  double compute(double p, int n1, int n2, int r1, int r) {
    // Check cache first
    PfutileKey key{p, n1, n2, r1, r};
    auto it = cache_.find(key);
    if (it != cache_.end()) {
      return it->second;
    }

    // Create distributions once
    boost::math::binomial_distribution<double> binom1(n1, p);
    boost::math::binomial_distribution<double> binom2(n2, p);

    double aval = boost::math::cdf(binom1, r1);

    int upper_limit = std::min(n1, r);
    for (int x = r1 + 1; x <= upper_limit; ++x) {
      double pmf_x = boost::math::pdf(binom1, x);
      double cdf_stage2 = boost::math::cdf(binom2, r - x);
      aval += pmf_x * cdf_stage2;
    }

    // Cache result
    cache_[key] = aval;
    return aval;
  }
};

// Binary search for maximum r satisfying beta constraint
int find_max_r(PfutileCalculator& calc, double pi, int n1, int n2,
               int r1, double beta) {
  int r_lower = r1;
  int r_upper = r1 + n2;
  int r_best = -1;

  // Quick check 1: If even minimum r violates, no solution exists
  if (calc.compute(pi, n1, n2, r1, r_lower) > beta) {
    return -1;
  }

  // Quick check 2: If maximum r satisfies, return it immediately
  if (calc.compute(pi, n1, n2, r1, r_upper) <= beta) {
    return r_upper;
  }

  // Binary search for the maximum r where betastar <= beta
  // Invariant: betastar is monotonically increasing in r
  while (r_lower <= r_upper) {
    int r_mid = r_lower + (r_upper - r_lower) / 2;
    double betastar = calc.compute(pi, n1, n2, r1, r_mid);

    if (betastar <= beta) {
      // Constraint satisfied at r_mid
      r_best = r_mid;
      r_lower = r_mid + 1;  // Try larger r (search right)
    } else {
      // Constraint violated at r_mid
      r_upper = r_mid - 1;  // Try smaller r (search left)
    }
  }

  return r_best;
}

DataFrameCpp simon2stagecpp(
    const double alpha,
    const double beta,
    const double piH0,
    const double pi,
    const int n_max) {

  // Input validation
  if (std::isnan(alpha)) throw std::invalid_argument("alpha must be provided");

  if (alpha < 0.00001 || alpha >= 1.0)
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");

  if (std::isnan(beta)) throw std::invalid_argument("beta must be provided");

  if (beta >= 1.0 - alpha || beta < 0.0001)
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");

  if (std::isnan(piH0)) throw std::invalid_argument("piH0 must be provided");

  if (piH0 <= 0.0 || piH0 >= 1.0)
    throw std::invalid_argument("piH0 must lie between 0 and 1");

  if (std::isnan(pi)) throw std::invalid_argument("pi must be provided");

  if (pi <= piH0 || pi >= 1.0)
    throw std::invalid_argument("pi must lie between piH0 and 1");

  double p = (piH0 + pi) / 2.0;
  double z1 = boost_qnorm(1.0 - alpha);
  double z2 = boost_qnorm(1.0 - beta);

  int n_min = static_cast<int>(std::floor(p * (1.0 - p) *
                               std::pow((z1 + z2) / (pi - piH0), 2)));
  int n_max1 = static_cast<int>(std::ceil(2.0 * n_min));

  int n_lower = static_cast<int>(std::floor(0.5 * n_min));
  int n_upper = std::min(n_max, n_max1);

  // Create reusable calculator with cache
  PfutileCalculator calc_pi, calc_piH0;

  // Storage for candidate designs
  std::vector<int> nx, n1x, r1x, rx;
  std::vector<double> en0x, pet0x, alphax, powerx;

  // Search over total sample sizes
  for (int n = n_lower; n <= n_upper; ++n) {
    double EN0_best = static_cast<double>(n);
    int n1_best = -1, r1_best = -1, r_best = -1;
    double pet0_best = 0.0, alpha_best = 0.0, power_best = 0.0;

    bool exist = false;

    // Narrow n1 search range using heuristics
    int n1_min = std::max(1, static_cast<int>(0.15 * n));
    int n1_max = std::min(n - 1, static_cast<int>(0.85 * n));

    for (int n1 = n1_min; n1 <= n1_max; ++n1) {
      int n2 = n - n1;

      // early termination: if n1 alone exceeds best EN0
      if (n1 >= EN0_best) {
        break;
      }

      // Create distribution for this n1 (reused across r1 loop)
      boost::math::binomial_distribution<double> binom_n1_piH0(n1, piH0);

      for (int r1 = 0; r1 <= n1; ++r1) {
        int r = find_max_r(calc_pi, pi, n1, n2, r1, beta);

        if (r >= r1) {
          double alphastar = 1.0 - calc_piH0.compute(piH0, n1, n2, r1, r);

          if (alphastar <= alpha) {
            exist = true;

            // Use the distribution created for this n1
            double pet0 = boost::math::cdf(binom_n1_piH0, r1);
            double en0 = n1 + (1.0 - pet0) * n2;

            if (en0 < EN0_best) {
              EN0_best = en0;
              n1_best = n1;
              r1_best = r1;
              r_best = r;
              pet0_best = pet0;
              alpha_best = alphastar;
              power_best = 1.0 - calc_pi.compute(pi, n1, n2, r1, r);
            }
          }
        }
      }
    }


    if (exist) {
      nx.push_back(n);
      n1x.push_back(n1_best);
      r1x.push_back(r1_best);
      rx.push_back(r_best);
      en0x.push_back(EN0_best);
      pet0x.push_back(pet0_best);
      alphax.push_back(alpha_best);
      powerx.push_back(power_best);
    }
  }

  if (nx.empty()) {
    throw std::runtime_error("No design found satisfying the constraints");
  }

  // [Rest of the code for convex hull and output formatting - same as before]
  auto min_it = std::min_element(en0x.begin(), en0x.end());
  int I = static_cast<int>(std::distance(en0x.begin(), min_it));

  nx.resize(I + 1);
  n1x.resize(I + 1);
  r1x.resize(I + 1);
  rx.resize(I + 1);
  en0x.resize(I + 1);
  pet0x.resize(I + 1);
  alphax.resize(I + 1);
  powerx.resize(I + 1);

  std::vector<int> u_indices;
  u_indices.push_back(0);

  int i = 0;
  while (i < I) {
    int best_j = -1;
    double min_slope = std::numeric_limits<double>::infinity();

    for (int j = i + 1; j <= I; ++j) {
      double slope = (en0x[j] - en0x[i]) / (nx[j] - nx[i]);
      if (slope < min_slope) {
        min_slope = slope;
        best_j = j;
      }
    }

    i = best_j;
    u_indices.push_back(i);
  }

  // Extract admissible designs
  subset_in_place(nx, u_indices);
  subset_in_place(n1x, u_indices);
  subset_in_place(r1x, u_indices);
  subset_in_place(rx, u_indices);
  subset_in_place(en0x, u_indices);
  subset_in_place(pet0x, u_indices);
  subset_in_place(alphax, u_indices);
  subset_in_place(powerx, u_indices);

  int m = static_cast<int>(nx.size());
  std::vector<double> w1(m), w2(m);
  std::vector<std::string> design(m);

  for (int i = 0; i < m; ++i) {
    if (i < m - 1) {
      double slope = (en0x[i + 1] - en0x[i]) / (nx[i + 1] - nx[i]);
      double w = slope / (slope - 1.0);

      if (i == 0) {
        w1[i] = w;
        w2[i] = 1.0;
        design[i] = "Minimax";
      } else {
        w1[i] = w;
        w2[i] = w1[i - 1];
        design[i] = "Admissible";
      }
    } else {
      w1[i] = 0.0;
      w2[i] = w1[i - 1];
      design[i] = "Optimal";
    }
  }

  // Build result DataFrame
  DataFrameCpp result;
  std::vector<double> piH0_vec(m, piH0);
  result.push_back(std::move(piH0_vec), "piH0");
  result.push_back(pi, "pi");
  result.push_back(alpha, "alpha");
  result.push_back(beta, "beta");
  result.push_back(std::move(nx), "n");
  result.push_back(std::move(n1x), "n1");
  result.push_back(std::move(r1x), "r1");
  result.push_back(std::move(rx), "r");
  result.push_back(std::move(en0x), "EN0");
  result.push_back(std::move(alphax), "attainedAlpha");
  result.push_back(std::move(powerx), "attainedPower");
  result.push_back(std::move(pet0x), "PET0");
  result.push_back(std::move(w1), "w_lower");
  result.push_back(std::move(w2), "w_upper");
  result.push_back(std::move(design), "design");

  return result;
}


//' @title Simon's Two-Stage Design
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
//' simon2stage(0.05, 0.2, 0.1, 0.3)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame simon2stage(
    const double alpha = NA_REAL,
    const double beta = NA_REAL,
    const double piH0 = NA_REAL,
    const double pi = NA_REAL,
    const int n_max = 110) {
  auto result = simon2stagecpp(alpha, beta, piH0, pi, n_max);
  return Rcpp::wrap(result);
}


// Helper for the analyis of Simon's Bayesian basket trials
ListCpp simonBayesAnalysiscpp(
    const int nstrata,
    const std::vector<double>& r,
    const std::vector<double>& n,
    const double lambda,
    const double gamma,
    const double phi,
    const double plo) {

  if (nstrata == INT_MIN) {
    throw std::invalid_argument("nstrata must be provided.");
  }
  if (nstrata <= 0) {
    throw std::invalid_argument("nstrata must be a positive integer.");
  }
  if (!none_na(r)) {
    throw std::invalid_argument("r must be provided.");
  }
  if (!none_na(n)) {
    throw std::invalid_argument("n must be provided.");
  }
  if (static_cast<int>(r.size()) != nstrata) {
    throw std::invalid_argument("Invalid length for r.");
  }
  if (static_cast<int>(n.size()) != nstrata) {
    throw std::invalid_argument("Invalid length for n.");
  }
  if (std::any_of(r.begin(), r.end(), [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("r must be nonnegative.");
  }

  for (std::size_t i = 0; i < r.size(); ++i) {
    if (r[i] > n[i]) {
      throw std::invalid_argument("r must be less than or equal to n.");
    }
  }

  if (std::isnan(lambda)) {
    throw std::invalid_argument("lambda must be provided.");
  }
  if (std::isnan(gamma)) {
    throw std::invalid_argument("gamma must be provided.");
  }
  if (std::isnan(phi)) {
    throw std::invalid_argument("phi must be provided.");
  }
  if (std::isnan(plo)) {
    throw std::invalid_argument("plo must be provided.");
  }
  if (lambda <= 0 || lambda >= 1) {
    throw std::invalid_argument("lambda must lie between 0 and 1.");
  }
  if (gamma <= 0 || gamma >= 1) {
    throw std::invalid_argument("gamma must lie between 0 and 1.");
  }
  if (phi <= 0 || phi >= 1) {
    throw std::invalid_argument("phi must lie between 0 and 1.");
  }
  if (plo <= 0 || plo >= 1) {
    throw std::invalid_argument("plo must lie between 0 and 1.");
  }
  if (plo >= phi) {
    throw std::invalid_argument("plo must be less than phi.");
  }

  int ncases = static_cast<int>(1u << nstrata);
  FlatMatrix incid(ncases, nstrata);
  std::vector<double> prior(ncases), like(ncases);
  for (int i = 0; i < ncases; ++i) {
    int number = ncases - i - 1;
    std::vector<signed char> cc(nstrata);
    for (int j = 0; j < nstrata; ++j) {
      cc[j] = (number >> (nstrata - 1 - j)) & 1u;
    }

    bool all_ones = std::all_of(cc.begin(), cc.end(), [](int x) { return x == 1; });
    bool all_zeros = std::all_of(cc.begin(), cc.end(), [](int x) { return x == 0; });

    if (all_ones || all_zeros) {
      prior[i] = lambda * (gamma * (cc[0] == 1) + (1 - gamma) * (cc[0] == 0)) +
        (1 - lambda) * (std::pow(gamma, nstrata) * (cc[0] == 1) +
        std::pow(1 - gamma, nstrata) * (cc[0] == 0));
    } else {
      double y = std::accumulate(cc.begin(), cc.end(), 0.0);
      prior[i] = (1 - lambda) * std::pow(gamma, y) *
        std::pow(1 - gamma, nstrata - y);
    }

    std::vector<double> x(nstrata);
    for (int j = 0; j < nstrata; ++j) {
      x[j] = phi * cc[j] + plo * (1 - cc[j]);
    }

    // Calculate sum of r*log(x) + (n-r)*log(1-x)
    double log_sum = 0.0;
    for (int j = 0; j < nstrata; ++j) {
      log_sum += r[j] * std::log(x[j]) + (n[j] - r[j]) * std::log(1 - x[j]);
    }
    like[i] = std::exp(log_sum);

    // Copy cc to row i of incid matrix
    for (int j = 0; j < nstrata; ++j) {
      incid(i, j) = cc[j];
    }
  }

  std::vector<double> post(ncases);
  for (int i = 0; i < ncases; ++i) {
    post[i] = prior[i] * like[i];
  }

  // Normalize: post = post / sum(post)
  double post_sum = std::accumulate(post.begin(), post.end(), 0.0);
  for (int i = 0; i < ncases; ++i) {
    post[i] /= post_sum;
  }

  // Calculate q_prior and q_post for each stratum
  std::vector<double> q_prior(nstrata, 0.0);
  std::vector<double> q_post(nstrata, 0.0);
  for (int j = 0; j < nstrata; ++j) {
    for (int i = 0; i < ncases; ++i) {
      q_prior[j] += incid(i, j) * prior[i];
      q_post[j] += incid(i, j) * post[i];
    }
  }

  ListCpp result;
  result.push_back(std::move(incid), "case");
  result.push_back(std::move(prior), "prior_case");
  result.push_back(std::move(q_prior), "prior_stratum");
  result.push_back(std::move(post), "post_case");
  result.push_back(std::move(q_post), "post_stratum");
  return result;
}


//' @title Analysis of Simon's Bayesian Basket Trials
//' @description Obtains the prior and posterior probabilities for
//' Simon's Bayesian basket discovery trials.
//'
//' @param nstrata The number of strata.
//' @param r The vector of number of responders across strata.
//' @param n The vector of number of subjects across strata.
//' @param lambda The prior probability that the drug activity is
//'   homogeneous across strata.
//' @param gamma The prior probability that the drug is active in a
//'   stratum.
//' @param phi The response probability for an active drug.
//' @param plo The response probability for an inactive drug.
//'
//' @return A list containing the following five components:
//'
//' * \code{case}: The matrix with each row corresponding to a combination
//'   of drug activity over strata represented by the columns.
//'
//' * \code{prior_case}: The vector of joint prior probabilities
//'   for the stratum-specific response rates.
//'
//' * \code{prior_stratum}: The vector of marginal prior probabilities
//'   for the stratum-specific response rates.
//'
//' * \code{post_case}: The vector of joint posterior probabilities
//'   for the stratum-specific response rates.
//'
//' * \code{post_stratum}: The vector of marginal posterior probabilities
//'   for the stratum-specific response rates.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' a = simonBayesAnalysis(
//'   nstrata = 10,
//'   r = c(8,0,1,1,6,2,0,0,3,3),
//'   n = c(19,10,26,8,14,7,8,5,4,14),
//'   lambda = 0.5, gamma = 0.33,
//'   phi = 0.35, plo = 0.15)
//'
//' a$post_stratum
//'
//' @export
// [[Rcpp::export]]
Rcpp::List simonBayesAnalysis(
    const int nstrata = NA_INTEGER,
    const Rcpp::NumericVector& r = NA_REAL,
    const Rcpp::NumericVector& n = NA_REAL,
    const double lambda = NA_REAL,
    const double gamma = NA_REAL,
    const double phi = NA_REAL,
    const double plo = NA_REAL) {
  auto r1 = Rcpp::as<std::vector<double>>(r);
  auto n1 = Rcpp::as<std::vector<double>>(n);
  auto result = simonBayesAnalysiscpp(nstrata, r1, n1, lambda, gamma, phi, plo);
  return Rcpp::wrap(result);
}


// Helper for the simulation of Simon's Bayesian basket trials
ListCpp simonBayesSimcpp(
    const std::vector<double>& p,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& stratumFraction,
    const double lambda,
    const double gamma,
    const double phi,
    const double plo,
    const double T,
    const int maxSubjects,
    const std::vector<int>& plannedSubjects,
    const int maxNumberOfIterations,
    const int maxNumberOfRawDatasets,
    const int seed) {

  int nstrata = static_cast<int>(stratumFraction.size());

  if (!none_na(p)) {
    throw std::invalid_argument("p must be provided.");
  }
  if (std::any_of(p.begin(), p.end(), [](double x) { return x <= 0 || x >= 1; })) {
    throw std::invalid_argument("p must lie between 0 and 1.");
  }

  if (accrualTime[0] != 0) {
    throw std::invalid_argument("accrualTime must start with 0.");
  }
  if (any_nonincreasing(accrualTime)) {
    throw std::invalid_argument("accrualTime should be increasing.");
  }
  if (!none_na(accrualIntensity)) {
    throw std::invalid_argument("accrualIntensity must be provided.");
  }
  if (accrualTime.size() != accrualIntensity.size()) {
    throw std::invalid_argument(
        "accrualTime must have the same length as accrualIntensity.");
  }
  if (std::any_of(accrualIntensity.begin(), accrualIntensity.end(),
                  [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("accrualIntensity must be non-negative.");
  }
  if (std::any_of(stratumFraction.begin(), stratumFraction.end(),
                  [](double val) { return val <= 0.0; })) {
    throw std::invalid_argument("stratumFraction must be positive.");
  }

  double sum_stratumFraction = std::accumulate(stratumFraction.begin(),
                                               stratumFraction.end(), 0.0);
  if (std::fabs(sum_stratumFraction - 1.0) > 1.0e-8) {
    throw std::invalid_argument("stratumFraction must sum to 1.");
  }

  if (std::isnan(lambda)) {
    throw std::invalid_argument("lambda must be provided.");
  }
  if (std::isnan(gamma)) {
    throw std::invalid_argument("gamma must be provided.");
  }
  if (std::isnan(phi)) {
    throw std::invalid_argument("phi must be provided.");
  }
  if (std::isnan(plo)) {
    throw std::invalid_argument("plo must be provided.");
  }
  if (std::isnan(T)) {
    throw std::invalid_argument("T must be provided.");
  }
  if (lambda <= 0 || lambda >= 1) {
    throw std::invalid_argument("lambda must lie between 0 and 1.");
  }
  if (gamma <= 0 || gamma >= 1) {
    throw std::invalid_argument("gamma must lie between 0 and 1.");
  }
  if (phi <= 0 || phi >= 1) {
    throw std::invalid_argument("phi must lie between 0 and 1.");
  }
  if (plo <= 0 || plo >= 1) {
    throw std::invalid_argument("plo must lie between 0 and 1.");
  }
  if (plo >= phi) {
    throw std::invalid_argument("plo must be less than phi.");
  }
  if (T <= 0 || T >= 1) {
    throw std::invalid_argument("T must lie between 0 and 1.");
  }
  if (maxSubjects == INT_MIN) {
    throw std::invalid_argument("maxSubjects must be provided.");
  }
  if (maxSubjects < 1) {
    throw std::invalid_argument("maxSubjects must be a positive integer.");
  }
  if (!none_na(plannedSubjects)) {
    throw std::invalid_argument("plannedSubjects must be provided.");
  }
  if (plannedSubjects[0] <= 0) {
    throw std::invalid_argument("plannedSubjects must be positive.");
  }
  if (any_nonincreasing(plannedSubjects)) {
    throw std::invalid_argument("plannedSubjects must be increasing.");
  }
  if (maxNumberOfIterations < 1) {
    throw std::invalid_argument("maxNumberOfIterations must be positive.");
  }
  if (maxNumberOfRawDatasets < 0) {
    throw std::invalid_argument("maxNumberOfRawDatasets must be non-negative.");
  }

  std::vector<unsigned char> act(p.size());
  for (std::size_t i = 0; i < p.size(); ++i) {
    act[i] = (p[i] == phi) ? 1 : 0;
  }
  int nactive = static_cast<int>(std::accumulate(act.begin(), act.end(), 0));

  std::vector<double> cumStratumFraction(nstrata);
  std::partial_sum(stratumFraction.begin(), stratumFraction.end(),
                   cumStratumFraction.begin());
  std::vector<double> post_stratum(nstrata), n(nstrata), r(nstrata);
  std::vector<unsigned char> open(nstrata), pos(nstrata), neg(nstrata);


  std::vector<double> arrivalTime(maxSubjects);
  std::vector<int> stratum(maxSubjects), y(maxSubjects);
  std::vector<int> iterationNumber(maxNumberOfIterations);
  std::vector<double> N(maxNumberOfIterations);
  std::vector<double> nact(maxNumberOfIterations), nopen(maxNumberOfIterations);
  std::vector<double> tpos(maxNumberOfIterations), fneg(maxNumberOfIterations);
  std::vector<double> fpos(maxNumberOfIterations), tneg(maxNumberOfIterations);
  std::vector<int> numberOfStrata(maxNumberOfIterations, nstrata);

  // cache for the patient-level raw data to extract
  int kMax = static_cast<int>(plannedSubjects.size());
  int nrow1 = kMax * maxNumberOfRawDatasets * maxSubjects;
  std::vector<int> iterationNumberx(nrow1);
  std::vector<int> stageNumberx(nrow1);
  std::vector<int> subjectIdx(nrow1);
  std::vector<double> arrivalTimex(nrow1);
  std::vector<int> stratumx(nrow1);
  std::vector<int> yx(nrow1);

  // cache for the summary data to extract
  int nrow2 = kMax * maxNumberOfIterations * nstrata;
  std::vector<int> iterationNumbery(nrow2);
  std::vector<int> stageNumbery(nrow2);
  std::vector<int> stratumy(nrow2);
  std::vector<unsigned char> activey(nrow2);
  std::vector<int> ny(nrow2);
  std::vector<int> ry(nrow2);
  std::vector<double> posty(nrow2);
  std::vector<unsigned char> openy(nrow2);
  std::vector<unsigned char> posy(nrow2);
  std::vector<unsigned char> negy(nrow2);

  // random number generator
  boost::random::mt19937_64 rng(seed);
  boost::random::uniform_real_distribution<double> unif(0.0, 1.0);

  int index1 = 0, index2 = 0;
  for (int iter = 0; iter < maxNumberOfIterations; ++iter) {
    // initialize the contents in each stratum
    std::fill(n.begin(), n.end(), 0.0);
    std::fill(r.begin(), r.end(), 0.0);
    std::fill(open.begin(), open.end(), 1);
    std::fill(pos.begin(), pos.end(), 0);
    std::fill(neg.begin(), neg.end(), 0);

    int k = 0;      // index of the number of subjects included in analysis
    int stage = 0;
    double enrollt = 0.0;
    for (int i = 0; i < 100000; ++i) {
      // generate accrual time
      double u = unif(rng);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt);

      // generate stratum information
      u = unif(rng);
      int j = 0;
      for (; j < nstrata; ++j) {
        if (cumStratumFraction[j] > u) break;
      }

      // if the stratum is open, generate the response for the subject
      if (open[j]) {
        arrivalTime[k] = enrollt;
        stratum[k] = j+1;
        y[k] = (unif(rng) < p[j]) ? 1 : 0;

        // update the number of subjects and responders in the stratum
        n[j] = n[j] + 1;
        r[j] = r[j] + y[k];

        ++k;

        // interim analysis
        if (std::any_of(plannedSubjects.begin(), plannedSubjects.end(),
                        [k](double val) { return val == k; })) {
          // output raw data
          if (iter < maxNumberOfRawDatasets) {
            for (int idx = 0; idx < k; ++idx) {
              iterationNumberx[index1] = iter + 1;
              stageNumberx[index1] = stage + 1;
              subjectIdx[index1] = idx + 1;
              arrivalTimex[index1] = arrivalTime[idx];
              stratumx[index1] = stratum[idx];
              yx[index1] = y[idx];

              ++index1;
            }
          }

          // calculate the posterior probabilities
          ListCpp a = simonBayesAnalysiscpp(nstrata, r, n, lambda, gamma, phi, plo);
          post_stratum = a.get<std::vector<double>>("post_stratum");

          // whether to close the stratum due to positive or negative results
          for (int l = 0; l < nstrata; ++l) {
            if (open[l]) {
              if (post_stratum[l] > T) {
                pos[l] = 1; open[l] = 0;
              } else if (post_stratum[l] < 1 - T) {
                neg[l] = 1; open[l] = 0;
              }
            }

            // output strata-level data
            iterationNumbery[index2] = iter + 1;
            stageNumbery[index2] = stage + 1;
            stratumy[index2] = l + 1;
            activey[index2] = act[l];
            ny[index2] = static_cast<int>(n[l]);
            ry[index2] = static_cast<int>(r[l]);
            posty[index2] = post_stratum[l];
            openy[index2] = open[l];
            posy[index2] = pos[l];
            negy[index2] = neg[l];

            ++index2;
          }

          ++stage;
        }

        // stop the trial if all strata are closed or max subjects reached
        bool all_closed = std::all_of(open.begin(), open.end(),
                                      [](unsigned char val) { return val == 0; });

        if (all_closed || (k == maxSubjects)) {
          iterationNumber[iter] = iter + 1;
          N[iter] = k;
          nact[iter] = nactive;

          // Calculate true positives: sum(pos & act)
          int sum_pos_act = 0;
          for (int i = 0; i < nstrata; ++i) {
            sum_pos_act += pos[i] & act[i];
          }
          tpos[iter] = sum_pos_act;

          // Calculate false negatives: sum(neg & act)
          int sum_neg_act = 0;
          for (int i = 0; i < nstrata; ++i) {
            sum_neg_act += neg[i] & act[i];
          }
          fneg[iter] = sum_neg_act;

          // Calculate false positives: sum(pos & !act)
          int sum_pos_notact = 0;
          for (int i = 0; i < nstrata; ++i) {
            sum_pos_notact += pos[i] & !act[i];
          }
          fpos[iter] = sum_pos_notact;

          // Calculate true negatives: sum(neg & !act)
          int sum_neg_notact = 0;
          for (int i = 0; i < nstrata; ++i) {
            sum_neg_notact += neg[i] & !act[i];
          }
          tneg[iter] = sum_neg_notact;

          // Calculate number of open: sum(open)
          int sum_open = std::accumulate(open.begin(), open.end(), 0);
          nopen[iter] = sum_open;

          break;
        }
      }
    }
  }

  // subset the cached data to the actual size
  if (maxNumberOfRawDatasets > 0) {
    subset_in_place(iterationNumberx, 0, index1);
    subset_in_place(stageNumberx, 0, index1);
    subset_in_place(subjectIdx, 0, index1);
    subset_in_place(arrivalTimex, 0, index1);
    subset_in_place(stratumx, 0, index1);
    subset_in_place(yx, 0, index1);
  }

  subset_in_place(iterationNumbery, 0, index2);
  subset_in_place(stageNumbery, 0, index2);
  subset_in_place(stratumy, 0, index2);
  subset_in_place(activey, 0, index2);
  subset_in_place(ny, 0, index2);
  subset_in_place(ry, 0, index2);
  subset_in_place(posty, 0, index2);
  subset_in_place(openy, 0, index2);
  subset_in_place(posy, 0, index2);
  subset_in_place(negy, 0, index2);

  // calculate the means
  double mn_nact = mean_kahan(nact);
  double mn_tpos = mean_kahan(tpos);
  double mn_fneg = mean_kahan(fneg);
  double mn_fpos = mean_kahan(fpos);
  double mn_tneg = mean_kahan(tneg);
  double mn_inde = mean_kahan(nopen);
  double mn_N = mean_kahan(N);

  DataFrameCpp rawdata;
  if (maxNumberOfRawDatasets > 0) {
    rawdata.push_back(std::move(iterationNumberx), "iterationNumber");
    rawdata.push_back(std::move(stageNumberx), "stageNumber");
    rawdata.push_back(std::move(subjectIdx), "subjectId");
    rawdata.push_back(std::move(arrivalTimex), "arrivalTime");
    rawdata.push_back(std::move(stratumx), "stratum");
    rawdata.push_back(std::move(yx), "y");
  }

  DataFrameCpp sumdata1;
  sumdata1.push_back(std::move(iterationNumbery), "iterationNumber");
  sumdata1.push_back(std::move(stageNumbery), "stageNumber");
  sumdata1.push_back(std::move(stratumy), "stratum");
  sumdata1.push_back(std::move(activey), "active");
  sumdata1.push_back(std::move(ny), "n");
  sumdata1.push_back(std::move(ry), "r");
  sumdata1.push_back(std::move(posty), "posterior");
  sumdata1.push_back(std::move(openy), "open");
  sumdata1.push_back(std::move(posy), "positive");
  sumdata1.push_back(std::move(negy), "negative");

  DataFrameCpp sumdata2;
  sumdata2.push_back(std::move(iterationNumber), "iterationNumber");
  sumdata2.push_back(std::move(numberOfStrata), "numberOfStrata");
  sumdata2.push_back(std::move(nact), "n_active_strata");
  sumdata2.push_back(std::move(tpos), "true_positive");
  sumdata2.push_back(std::move(fneg), "false_negative");
  sumdata2.push_back(std::move(fpos), "false_positive");
  sumdata2.push_back(std::move(tneg), "true_negative");
  sumdata2.push_back(std::move(nopen), "n_indet_strata");
  sumdata2.push_back(std::move(N), "numberOfSubjects");

  DataFrameCpp overview;
  overview.push_back(nstrata, "numberOfStrata");
  overview.push_back(mn_nact, "n_active_strata");
  overview.push_back(mn_tpos, "true_positive");
  overview.push_back(mn_fneg, "false_negative");
  overview.push_back(mn_fpos, "false_positive");
  overview.push_back(mn_tneg, "true_negative");
  overview.push_back(mn_inde, "n_indet_strata");
  overview.push_back(mn_N, "numberOfSubjects");

  ListCpp result;
  if (maxNumberOfRawDatasets > 0) {
    result.push_back(std::move(rawdata), "rawdata");
    result.push_back(std::move(sumdata1), "sumdata1");
    result.push_back(std::move(sumdata2), "sumdata2");
    result.push_back(std::move(overview), "overview");
  } else {
    result.push_back(std::move(sumdata1), "sumdata1");
    result.push_back(std::move(sumdata2), "sumdata2");
    result.push_back(std::move(overview), "overview");
  }
  return result;
}


//' @title Simulation of Simon's Bayesian Basket Trials
//' @description Obtains the simulated raw and summary data for Simon's
//' Bayesian basket discovery trials.
//'
//' @param p The vector of true response probabilities across strata.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_stratumFraction
//' @param lambda The prior probability that the drug activity is
//'   homogeneous across strata.
//' @param gamma The prior probability that the drug is active in a
//'   stratum.
//' @param phi The response probability for an active drug.
//' @param plo The response probability for an inactive drug.
//' @param T The threshold for a conclusive posterior probability to
//'   stop enrollment.
//' @param maxSubjects The maximum total sample size.
//' @param plannedSubjects The planned cumulative number of subjects
//'   at each stage.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasets The number of raw datasets to extract.
//' @param seed The seed to reproduce the simulation results.
//'   The seed from the environment will be used if left unspecified,
//'
//' @return A list containing the following four components:
//'
//' * \code{rawdata}: A data frame for subject-level data, containing
//'   the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage number.
//'
//'     - \code{subjectId}: The subject ID.
//'
//'     - \code{arrivalTime}: The enrollment time for the subject.
//'
//'     - \code{stratum}: The stratum for the subject.
//'
//'     - \code{y}: Whether the subject was a responder (1) or
//'       nonresponder (0).
//'
//' * \code{sumdata1}: A data frame for simulation and stratum-level
//'   summary data, containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage number.
//'
//'     - \code{stratum}: The stratum number.
//'
//'     - \code{active}: Whether the drug is active in the stratum.
//'
//'     - \code{n}: The number of subjects in the stratum.
//'
//'     - \code{r}: The number of responders in the stratum.
//'
//'     - \code{posterior}: The posterior probability that the drug is
//'       active in the stratum.
//'
//'     - \code{open}: Whether the stratum is still open for enrollment.
//'
//'     - \code{positive}: Whether the stratum has been determined to be
//'       a positive stratum.
//'
//'     - \code{negative}: Whether the stratum has been determined to be
//'       a negative stratum.
//'
//' * \code{sumdata2}: A data frame for the simulation level summary data,
//'   containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{numberOfStrata}: The total number of strata.
//'
//'     - \code{n_active_strata}: The number of active strata.
//'
//'     - \code{true_positive}: The number of true positive strata.
//'
//'     - \code{false_negative}: The number of false negative strata.
//'
//'     - \code{false_positive}: The number of false positive strata.
//'
//'     - \code{true_negative}: The number of true negative strata.
//'
//'     - \code{n_indet_strata}: The number of indeterminate strata.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//' * \code{overview}: A data frame for the summary across simulations,
//'   containing the following variables:
//'
//'     - \code{numberOfStrata}: The total number of strata.
//'
//'     - \code{n_active_strata}: The average number of active strata.
//'
//'     - \code{true_positive}: The average number of true positive strata.
//'
//'     - \code{false_negative}: The average number of false negative strata.
//'
//'     - \code{false_positive}: The average number of false positive strata.
//'
//'     - \code{true_negative}: The average number of true negative strata.
//'
//'     - \code{n_indet_strata}: The average number of indeterminate strata.
//'
//'     - \code{numberOfSubjects}: The average number of subjects.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' sim1 = simonBayesSim(
//'   p = c(0.25, 0.25, 0.05),
//'   accrualIntensity = 5,
//'   stratumFraction = c(1/3, 1/3, 1/3),
//'   lambda = 0.33, gamma = 0.5,
//'   phi = 0.25, plo = 0.05,
//'   T = 0.8, maxSubjects = 50,
//'   plannedSubjects = seq(5, 50, 5),
//'   maxNumberOfIterations = 1000,
//'   maxNumberOfRawDatasets = 1,
//'   seed = 314159)
//'
//' sim1$overview
//'
//' @export
// [[Rcpp::export]]
Rcpp::List simonBayesSim(
    const Rcpp::NumericVector& p = NA_REAL,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& stratumFraction = 1,
    const double lambda = NA_REAL,
    const double gamma = NA_REAL,
    const double phi = NA_REAL,
    const double plo = NA_REAL,
    const double T = NA_REAL,
    const int maxSubjects = NA_INTEGER,
    const Rcpp::IntegerVector& plannedSubjects = NA_INTEGER,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasets = 1,
    const int seed = 0) {
  auto p1 = Rcpp::as<std::vector<double>>(p);
  auto accrualTime1 = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualIntensity1 = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto stratumFraction1 = Rcpp::as<std::vector<double>>(stratumFraction);
  auto plannedSubjects1 = Rcpp::as<std::vector<int>>(plannedSubjects);
  auto result = simonBayesSimcpp(
    p1, accrualTime1, accrualIntensity1, stratumFraction1,
    lambda, gamma, phi, plo, T, maxSubjects, plannedSubjects1,
    maxNumberOfIterations, maxNumberOfRawDatasets, seed);
  return Rcpp::wrap(result);
}
