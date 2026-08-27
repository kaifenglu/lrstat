#include "adaptive_two_stage.h"
#include "dataframe_list.h"
#include "multiplicity.h"
#include "mvnormr.h"
#include "survival_analysis.h"
#include "thread_utils.h"
#include "utilities.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using std::size_t;

static double log_beta(double a, double b) {
  return boost::math::lgamma(a) + boost::math::lgamma(b) -
         boost::math::lgamma(a + b);
}

// P(X > Y), X ~ Beta(a,b), Y ~ Beta(c,d); requires c, d positive integers
double beta_ineq_exact(double a, double b, double c, double d) {
  double lba = log_beta(a, b), s = 0.0;
  int C = static_cast<int>(c);
  for (int i = 0; i < C; ++i) {
    s += std::exp(log_beta(a + i, b + d) - std::log(d + i) -
                  log_beta(1.0 + i, d) - lba);
  }
  return s;
}

// P(X > Y), X ~ Beta(a,b), Y ~ Beta(c,d) independent
double beta_ineq(double a, double b, double c, double d, double tol = 1e-4) {
  auto f = [a, b, c, d](double x) -> double {
    return boost::math::ibeta_derivative(a, b, x) * boost::math::ibeta(c, d, x);
  };
  std::vector<double> breaks = {0.0, a / (a + b), 1.0};
  return integrate3(f, breaks, tol);
}

// Select the optimal biological dose at the end of phase 2.
//
// x0:    number of short-term responders in the control arm.
// xe:    number of short-term responders in each of the ndose doses.
// xt:    number of subjects experiencing toxicity in each dose.
// w:     weight on the posterior mean toxicity rate in the benefit-risk
//        tradeoff; 0 selects on efficacy alone.
// n:     number of subjects per arm, with n[0] for the control arm and
//        n[i+1] for dose i, so it has length ndose + 1.
// phi_t: prespecified upper limit for the toxicity rate; 1 disables the
//        safety criterion.
// ce:    threshold for the posterior probability that a dose beats the
//        control in response rate; 0 disables the efficacy criterion.
// ct:    threshold for the posterior probability that the toxicity rate of a
//        dose is below phi_t; 0 disables the safety criterion.
// uniform_prior: true for the uniform Beta(1,1) prior (exact posterior
//        comparison via beta_ineq_exact); false for the Jeffreys Beta(0.5,
//        0.5) prior, whose non-integer posterior parameters require the
//        numerical integration in beta_ineq.
//
// Returns a 0/1 indicator vector of length ndose flagging the selected dose.
// It is all zeros when no dose satisfies both criteria.
std::vector<int> get_decision(double x0, const std::vector<double> &xe,
                              const std::vector<double> &xt, double w,
                              const std::vector<double> &n, double phi_t,
                              double ce, double ct, bool uniform_prior) {
  size_t ndose = xe.size();
  double prior = uniform_prior ? 1.0 : 0.5;

  // posterior means under a beta-binomial model with the chosen prior
  std::vector<double> ppe(ndose), ppt(ndose);
  for (size_t i = 0; i < ndose; ++i) {
    ppe[i] = (prior + xe[i]) / (2.0 * prior + n[i + 1]);
    ppt[i] = (prior + xt[i]) / (2.0 * prior + n[i + 1]);
  }

  // the set of doses satisfying the safety and efficacy criteria
  std::vector<size_t> A;
  if (ce != 0.0 || ct != 0.0) {
    for (size_t i = 0; i < ndose; ++i) {
      // P(response rate of dose i > response rate of control)
      double probe =
          uniform_prior
              ? beta_ineq_exact(prior + xe[i], prior + n[i + 1] - xe[i],
                                prior + x0, prior + n[0] - x0)
              : beta_ineq(prior + xe[i], prior + n[i + 1] - xe[i], prior + x0,
                          prior + n[0] - x0);

      // P(toxicity rate of dose i < phi_t)
      double probt =
          boost::math::ibeta(prior + xt[i], prior + n[i + 1] - xt[i], phi_t);

      if (probe > ce && probt > ct) {
        A.push_back(i);
      }
    }
  } else {
    A.resize(ndose);
    std::iota(A.begin(), A.end(), static_cast<size_t>(0));
  }

  // select the dose with the highest utility value
  std::vector<int> select(ndose, 0);
  if (!A.empty()) {
    std::vector<double> u(ndose, 0.0);
    for (size_t id : A) {
      u[id] = ppe[id] - w * ppt[id];
    }
    size_t imax =
        static_cast<size_t>(std::max_element(u.begin(), u.end()) - u.begin());
    select[imax] = 1;
  }

  return select;
}

namespace {

constexpr double ALPHA_ONE_SIDED = 0.025;

enum : size_t {
  M_CTDUNNETT = 0,
  M_CTSIMES,
  M_CTPOOLED,
  M_CER,
  M_TSSSD_K,
  M_TSSSD_UK,
  M_TSSSD_K_RANK,
  M_TSSSD_K_MODRANK,
  M_TSSSD_UK_MODRANK,
  M_TSSSD_K_CE,
  M_TSSSD_UK_CE,
  M_TSSSD_K_RANK_CE,
  M_TSSSD_K_MODRANK_CE,
  M_TSSSD_UK_MODRANK_CE,
  M_NAIVE,
  M_PH3ONLY,
  NMETHOD
};

const char *const METHOD_NAME[NMETHOD] = {
  "ctdunnett", "ctsimes",      "ctpooled",      "cer",   "tsssd.k",
  "tsssd.uk", "tsssd.k.rank", "tsssd.k.modrank", "tsssd.uk.modrank",
  "tsssd.k.ce", "tsssd.uk.ce", "tsssd.k.rank.ce", "tsssd.k.modrank.ce",
  "tsssd.uk.modrank.ce",
  "naive", "ph3only"};

const double SINGLE_ARM_BOUND = boost_qnorm(1.0 - ALPHA_ONE_SIDED);

// Final boundary after selecting the largest phase-2 Wald statistic, with no
// phase-2 efficacy stopping. The selected-arm final statistics have an
// equicorrelated multivariate normal distribution under the global null.
double final_selection_boundary(size_t m, bool corr_known,
                                double information_fraction) {
  if (m == 1)
    return SINGLE_ARM_BOUND;

  const double phase2_corr = corr_known ? 0.5 : 0.0;
  const double final_corr = 1.0 - (1.0 - phase2_corr) * information_fraction;
  std::vector<double> lower(m, -POS_INF);
  std::vector<double> upper(m);

  if (m == 2) {
    auto f = [&](double critical_value) {
      std::fill(upper.begin(), upper.end(), critical_value);
      return 1.0 - pbvnormcpp(lower, upper, final_corr) - ALPHA_ONE_SIDED;
    };
    return brent(f, 0.0, 8.0, 1e-6);
  }

  std::vector<double> mean(m, 0.0);
  FlatMatrix sigma(m, m);
  sigma.fill(final_corr);
  for (size_t i = 0; i < m; ++i)
    sigma(i, i) = 1.0;

  auto f = [&](double critical_value) {
    std::fill(upper.begin(), upper.end(), critical_value);
    return 1.0 - pmvnormcpp(lower, upper, mean, sigma).prob - ALPHA_ONE_SIDED;
  };
  return brent(f, 0.0, 8.0, 1e-6);
}

// Final boundary after selecting the rankp0-th largest phase-2 statistic.
double rank_selection_boundary(size_t m, size_t rankp0,
                               double information_fraction) {
  if (m == 1)
    return SINGLE_ARM_BOUND;

  const double phase2_corr = 0.5;
  const double final_corr =
      1.0 - (1.0 - phase2_corr) * information_fraction;
  std::vector<double> lower(m);
  std::vector<double> upper(m);
  std::vector<double> mean(m, 0.0);
  FlatMatrix sigma(m, m);
  sigma.fill(final_corr);
  for (size_t i = 0; i < m; ++i)
    sigma(i, i) = 1.0;

  auto f = [&](double critical_value) {
    double probability = 0.0;
    for (size_t i = rankp0; i <= m; ++i) {
      for (size_t j = 0; j < m; ++j) {
        lower[j] = j < i ? critical_value : -POS_INF;
        upper[j] = j < i ? POS_INF : critical_value;
      }

      double combinations = 1.0;
      for (size_t j = 1; j <= i; ++j)
        combinations *= static_cast<double>(m - i + j) / j;
      probability +=
          combinations * pmvnormcpp(lower, upper, mean, sigma).prob;
    }
    return probability - ALPHA_ONE_SIDED;
  };

  return brent(f, 0.0, 8.0, 1e-6);
}

struct RawBinaryRow {
  int iteration = 0;
  int subject = 0;
  int treatment = 0;
  int response = 0;
  int toxicity = NA_INTEGER;
};

struct RawTteRow {
  int iteration = 0;
  int subject = 0;
  int phase = 0;
  int treatment = 0;
  int response = 0;
  double arrival = 0.0;
  double survival = 0.0;
  double observed = 0.0;
  int event = 0;
};

struct BinarySummaryRow {
  int iteration = 0;
  int treatment = 0;
  int responses = 0;
  int toxicities = NA_INTEGER;
  int selected = 0;
};

struct TteSummaryRow {
  int iteration = 0;
  int selectedDose = 0;
  int phase3SampleSize = 0;
  int stage1Events = 0;
  int stage2Events = 0;
  int totalEvents = 0;
  double stage1LogRankZ = NaN;
  double stage2LogRankZ = NaN;
  double cumulativeLogRankZ = NaN;
  std::vector<int> reject;
};

// per-trial storage so that worker threads never touch shared state
struct TrialResult {
  bool completed = false;
  std::vector<int> select;
  std::vector<IntMatrix> rej; // one ngrid by M matrix per requested method
  IntMatrix events;
  std::vector<RawBinaryRow> rawBinary;
  std::vector<RawTteRow> rawTte;
  std::vector<BinarySummaryRow> binarySummary;
  std::vector<TteSummaryRow> tteSummary;
};

// one-sided log-rank statistic, positive values favoring the treated group
double logrank_zg(std::vector<double> &&time, std::vector<int> &&event,
                  std::vector<int> &&treat) {
  DataFrameCpp df;
  df.push_back(std::move(time), "time");
  df.push_back(std::move(event), "event");
  df.push_back(std::move(treat), "treat");
  DataFrameCpp out =
      lrtestcpp(df, {""}, "treat", "time", "", "event", "", false, 0.0, 0.0);
  double z = out.get<double>("logRankZ")[0];
  return std::isnan(z) ? 0.0 : -z;
}

struct SimWorker : public RcppParallel::Worker {
  const size_t M;
  const size_t n1;
  const size_t n2min;
  const size_t n2max;
  const double p0;
  const std::vector<double> &pe;
  const std::vector<double> &pt;
  const double rho;
  const std::vector<double> &the0;
  const FlatMatrix &the1;
  const double T_max;
  const double w;
  const double phi_t;
  const double ce;
  const double ct;
  const double acc_rate1;
  const double acc_rate2;
  const double T_ph2followup;
  const size_t maxRawDatasets;
  const std::vector<uint64_t> &seeds;
  const std::vector<unsigned char> &use;
  const WeightMatrix &wgtmat;
  const WeightMatrix &wgtmat1;
  const BoolMatrix &family;
  const FlatMatrix &corr;
  const StageBoundaries &sb;
  const bool uniform_prior;

  std::vector<TrialResult> *results;

  SimWorker(size_t M_, size_t n1_, size_t n2min_, size_t n2max_, double p0_,
            const std::vector<double> &pe_, const std::vector<double> &pt_,
            double rho_, const std::vector<double> &the0_,
            const FlatMatrix &the1_, double T_max_, double w_, double phi_t_,
            double ce_, double ct_, double acc_rate1_, double acc_rate2_,
            double T_ph2followup_, size_t maxRawDatasets_,
            const std::vector<uint64_t> &seeds_,
            const std::vector<unsigned char> &use_, const WeightMatrix &wgtmat_,
            const WeightMatrix &wgtmat1_, const BoolMatrix &family_,
            const FlatMatrix &corr_, const StageBoundaries &sb_,
            bool uniform_prior_, std::vector<TrialResult> *results_)
      : M(M_), n1(n1_), n2min(n2min_), n2max(n2max_), p0(p0_), pe(pe_), pt(pt_),
        rho(rho_), the0(the0_), the1(the1_), T_max(T_max_), w(w_),
        phi_t(phi_t_), ce(ce_), ct(ct_), acc_rate1(acc_rate1_),
        acc_rate2(acc_rate2_), T_ph2followup(T_ph2followup_),
        maxRawDatasets(maxRawDatasets_), seeds(seeds_),
        use(use_), wgtmat(wgtmat_), wgtmat1(wgtmat1_), family(family_),
        corr(corr_), sb(sb_), uniform_prior(uniform_prior_), results(results_) {
  }

  void operator()(std::size_t begin, std::size_t end) {
    const size_t ntot = n1 + n2max;
    const size_t narm = M + 1; // arm 0 is the control arm
    const size_t ngrid = n2max - n2min + 1;
    const size_t ntests = (static_cast<size_t>(1) << M) - 1;
    const double sqrt1mrho2 = std::sqrt(1.0 - rho * rho);

    // buffers reused across the iterations handled by this worker
    std::vector<std::vector<unsigned char>> shortv(
        narm, std::vector<unsigned char>(ntot));
    std::vector<std::vector<unsigned char>> toxv(
      narm, std::vector<unsigned char>(n1));
    std::vector<std::vector<double>> longv(narm, std::vector<double>(ntot));
    std::vector<std::vector<double>> accr(narm, std::vector<double>(ntot));
    std::vector<std::vector<double>> ytime(narm, std::vector<double>(ntot));
    std::vector<std::vector<int>> yevent(narm, std::vector<int>(ntot));
    std::vector<double> xt(M), xe1(M), narms(narm, static_cast<double>(n1));
    std::vector<double> stg1_p(M), stg1_z(M), pooled_p(ntests);

    std::vector<size_t> allInt(ntests), allHyp(M);
    std::iota(allInt.begin(), allInt.end(), static_cast<size_t>(0));
    std::iota(allHyp.begin(), allHyp.end(), static_cast<size_t>(0));

    const std::vector<size_t> noStg1Rej;

    // Precompute nominal boundaries for CE-updated TSSSD methods.
    // These depend only on phase-2 and phase-3 sample sizes (n1, n2cur).
    const bool needTsssdCe =
      use[M_TSSSD_K_CE] || use[M_TSSSD_UK_CE] ||
      use[M_TSSSD_K_RANK_CE] || use[M_TSSSD_K_MODRANK_CE] ||
      use[M_TSSSD_UK_MODRANK_CE];
    std::vector<double> tNominal(ngrid, 0.5);
    // Indexed by effective hypothesis count minus one:
    // [0] => mEff=1, [M-1] => mEff=M.
    std::vector<std::vector<double>> tsssdKnownNomByEff(
      M, std::vector<double>(ngrid, SINGLE_ARM_BOUND));
    std::vector<std::vector<double>> tsssdUnknownNomByEff(
      M, std::vector<double>(ngrid, SINGLE_ARM_BOUND));
    std::vector<std::vector<double>> tsssdKnownNomByRank(
      M, std::vector<double>(ngrid, SINGLE_ARM_BOUND));
    if (needTsssdCe) {
      for (size_t n2i = 0; n2i < ngrid; ++n2i) {
        const size_t n2cur = n2min + n2i;
        const size_t nend = n1 + n2cur;
        const double tnom = static_cast<double>(n1) / static_cast<double>(nend);
        tNominal[n2i] = std::min(std::max(tnom, 1e-6), 1.0 - 1e-6);
        for (size_t mEff = 2; mEff <= M; ++mEff) {
          const size_t effIdx = mEff - 1;
          tsssdKnownNomByEff[effIdx][n2i] =
              final_selection_boundary(mEff, true, tNominal[n2i]);
          tsssdUnknownNomByEff[effIdx][n2i] =
              final_selection_boundary(mEff, false, tNominal[n2i]);
        }
        for (size_t rankp0 = 1; rankp0 <= M; ++rankp0) {
          const size_t rankIdx = rankp0 - 1;
          tsssdKnownNomByRank[rankIdx][n2i] =
              rank_selection_boundary(M, rankp0, tNominal[n2i]);
        }
      }
    }
    for (size_t iter = begin; iter < end; ++iter) {
      try {
        boost::random::mt19937_64 rng(seeds[iter]);
        boost::random::uniform_real_distribution<double> unif(0.0, 1.0);
        boost::random::normal_distribution<double> norm(0.0, 1.0);

        TrialResult &out = (*results)[iter];
        const bool saveRaw = iter < maxRawDatasets;
        out.select.assign(M, 0);
        out.rej.assign(NMETHOD, IntMatrix());
        for (size_t m = 0; m < NMETHOD; ++m) {
          if (use[m])
            out.rej[m].resize(ngrid, M);
        }
        out.events.resize(ngrid, 3);

        // binary efficacy and toxicity endpoints from a bivariate latent
        // normal with correlation rho, and the long-term endpoint that the
        // efficacy endpoint drives
        for (size_t a = 0; a < narm; ++a) {
          double qe = boost_qnorm(a == 0 ? p0 : pe[a - 1]);
          double qt = (a == 0) ? 0.0 : boost_qnorm(pt[a - 1]);
          double ntox = 0.0;
          for (size_t i = 0; i < ntot; ++i) {
            double zE = norm(rng);
            unsigned char s = (zE <= qe) ? 1u : 0u;
            shortv[a][i] = s;

            // toxicity is only assessed for the phase 2 cohort
            if (a > 0 && i < n1) {
              double zT = rho * zE + sqrt1mrho2 * norm(rng);
              toxv[a][i] = (zT <= qt) ? 1u : 0u;
              if (toxv[a][i])
                ntox += 1.0;
            }

            double rate = (a == 0 ? the0[s] : the1(a - 1, s));
            longv[a][i] = -std::log(unif(rng)) / rate;
          }
          if (a > 0)
            xt[a - 1] = ntox;
        }

        // phase 2 accrual from a homogeneous Poisson process
        double lastPh2 = 0.0;
        for (size_t a = 0; a < narm; ++a) {
          double t = 0.0;
          for (size_t i = 0; i < n1; ++i) {
            t += -std::log(unif(rng)) / acc_rate1;
            accr[a][i] = t;
          }
          if (t > lastPh2)
            lastPh2 = t;
        }

        // phase 3 enrollment opens T_ph2followup after the last phase 2
        // enrollment across all arms
        const double ph3start = lastPh2 + T_ph2followup;
        for (size_t a = 0; a < narm; ++a) {
          double t = ph3start;
          for (size_t i = n1; i < ntot; ++i) {
            t += -std::log(unif(rng)) / acc_rate2;
            accr[a][i] = t;
          }
        }

        // phase 2 dose selection based on the short-term endpoints
        double x01 = 0.0;
        for (size_t i = 0; i < n1; ++i)
          x01 += shortv[0][i];
        for (size_t k = 0; k < M; ++k) {
          double s = 0.0;
          for (size_t i = 0; i < n1; ++i)
            s += shortv[k + 1][i];
          xe1[k] = s;
        }

        out.select =
            get_decision(x01, xe1, xt, w, narms, phi_t, ce, ct, uniform_prior);

        size_t obd = 0;
        bool selected = false;
        for (size_t k = 0; k < M; ++k) {
          if (out.select[k]) {
            obd = k;
            selected = true;
            break;
          }
        }
        out.binarySummary.reserve(narm);
        out.binarySummary.push_back(
            {static_cast<int>(iter + 1), 0, static_cast<int>(x01), NA_INTEGER,
             0});
        for (size_t k = 0; k < M; ++k) {
          out.binarySummary.push_back(
              {static_cast<int>(iter + 1), static_cast<int>(k + 1),
               static_cast<int>(xe1[k]), static_cast<int>(xt[k]),
               out.select[k]});
        }

        // administrative censoring at the end of the study
        const double anal_time = ph3start + T_max;
        for (size_t a = 0; a < narm; ++a) {
          for (size_t i = 0; i < ntot; ++i) {
            double fu = anal_time - accr[a][i];
            if (fu < 0.0)
              fu = 0.0;
            if (longv[a][i] <= fu) {
              ytime[a][i] = longv[a][i];
              yevent[a][i] = 1;
            } else {
              ytime[a][i] = fu;
              yevent[a][i] = 0;
            }
          }
        }

        if (saveRaw) {
          out.rawBinary.reserve(narm * n1);
          out.rawTte.reserve(narm * n1 + (selected ? 2 * n2max : 0));
          for (size_t a = 0; a < narm; ++a) {
            for (size_t i = 0; i < n1; ++i) {
              const int subject = static_cast<int>(a * ntot + i + 1);
              out.rawBinary.push_back(
                  {static_cast<int>(iter + 1), subject, static_cast<int>(a),
                   static_cast<int>(shortv[a][i]),
                   a == 0 ? NA_INTEGER : static_cast<int>(toxv[a][i])});
              out.rawTte.push_back(
                  {static_cast<int>(iter + 1), subject, 2, static_cast<int>(a),
                   static_cast<int>(shortv[a][i]), accr[a][i], longv[a][i],
                   ytime[a][i], yevent[a][i]});
            }
          }
          if (selected) {
            for (size_t a : std::vector<size_t>{0, obd + 1}) {
              for (size_t i = n1; i < ntot; ++i) {
                out.rawTte.push_back(
                    {static_cast<int>(iter + 1),
                     static_cast<int>(a * ntot + i + 1), 3,
                     static_cast<int>(a), static_cast<int>(shortv[a][i]),
                     accr[a][i], longv[a][i], ytime[a][i], yevent[a][i]});
              }
            }
          }
        }

        if (!selected) {
          out.tteSummary.reserve(ngrid);
          for (size_t n2i = 0; n2i < ngrid; ++n2i) {
            out.tteSummary.push_back(
                {static_cast<int>(iter + 1), 0,
                 static_cast<int>(n2min + n2i), 0, 0, 0, NaN, NaN, NaN,
                 std::vector<int>(NMETHOD, 0)});
          }
          out.completed = true;
          continue;
        }

        // Compute a one-sided log-rank statistic by pooling treatment arms
        // against control.
        //
        // trt: indices of treatment arms to pool together.
        // lo:  lower index (inclusive) of subjects to include in the analysis.
        // hi:  upper index (exclusive) of subjects to include in the analysis.
        //
        // Returns a one-sided log-rank z-statistic (positive values favor
        // treated groups). The statistic pools the selected treatment arms
        // (from trt) against the control arm (arm 0), using survival times
        // and event indicators from subjects in range [lo, hi).
        auto zg = [&](const std::vector<size_t> &trt, size_t lo, size_t hi) {
          size_t sz = (trt.size() + 1) * (hi - lo);
          std::vector<double> tv;
          tv.reserve(sz);
          std::vector<int> ev;
          ev.reserve(sz);
          std::vector<int> gv;
          gv.reserve(sz);
          for (size_t a : trt) {
            for (size_t i = lo; i < hi; ++i) {
              tv.push_back(ytime[a][i]);
              ev.push_back(yevent[a][i]);
              gv.push_back(1);
            }
          }
          for (size_t i = lo; i < hi; ++i) {
            tv.push_back(ytime[0][i]);
            ev.push_back(yevent[0][i]);
            gv.push_back(0);
          }
          return logrank_zg(std::move(tv), std::move(ev), std::move(gv));
        };

        // stage 1 one-sided p-values, one per dose
        std::vector<size_t> arms;
        for (size_t k = 0; k < M; ++k) {
          arms.assign(1, k + 1);
          stg1_z[k] = zg(arms, 0, n1);
          stg1_p[k] = boost_pnorm(stg1_z[k], 0.0, 1.0, false);
        }

        LocalPValues dunnett1, simes1, pooled1;
        if (use[M_CTDUNNETT]) {
          dunnett1 = fPCStagewiseCpp(stg1_p, wgtmat, family, corr, allInt,
                                     allHyp, wgtmat, "dunnett");
        }
        if (use[M_CTSIMES]) {
          simes1 = fPCStagewiseCpp(stg1_p, wgtmat, family, corr, allInt, allHyp,
                                   wgtmat, "simes");
        }
        if (use[M_CTPOOLED]) {
          for (size_t i = 0; i < ntests; ++i) {
            arms.clear();
            for (size_t k = 0; k < M; ++k) {
              if (wgtmat.inthyp(i, k))
                arms.push_back(k + 1);
            }
            pooled_p[i] = boost_pnorm(zg(arms, 0, n1), 0.0, 1.0, false);
          }
          pooled1 = LocalPValues{allInt, wgtmat.inthyp, pooled_p};
        }

        // only intersection hypotheses containing the selected dose can keep
        // that dose from being rejected
        std::vector<size_t> idxJ;
        for (size_t i = 0; i < ntests; ++i) {
          if (wgtmat.inthyp(i, obd))
            idxJ.push_back(i);
        }
        IntMatrix inthypJ(idxJ.size(), M);
        for (size_t r = 0; r < idxJ.size(); ++r) {
          for (size_t c = 0; c < M; ++c) {
            inthypJ(r, c) = wgtmat.inthyp(idxJ[r], c);
          }
        }
        const std::vector<size_t> stg2_elemhyp{obd};

        const bool needStg2 =
          use[M_CTDUNNETT] || use[M_CTSIMES] || use[M_CTPOOLED];

        out.tteSummary.reserve(ngrid);
        for (size_t n2i = 0; n2i < ngrid; ++n2i) {
          const size_t n2cur = n2min + n2i;
          const size_t nend = n1 + n2cur;
          arms.assign(1, obd + 1);

          int d1e = 0, d2e = 0;
          for (size_t i = 0; i < n1; ++i) {
            d1e += yevent[obd + 1][i] + yevent[0][i];
          }
          for (size_t i = n1; i < nend; ++i) {
            d2e += yevent[obd + 1][i] + yevent[0][i];
          }
          out.events(n2i, 0) = d1e;
          out.events(n2i, 1) = d2e;
          out.events(n2i, 2) = d1e + d2e;

          // stage 2 test restricted to the selected dose
          const double z2 = zg(arms, n1, nend);
          double p2 = boost_pnorm(z2, 0.0, 1.0, false);
          LocalPValues stg2;
          if (needStg2) {
            stg2 = LocalPValues{idxJ, inthypJ,
                                std::vector<double>(idxJ.size(), p2)};
          }

          // test using the pooled stage 1 and stage 2 data
          const double zgn = zg(arms, 0, nend);
          double p_cum = boost_pnorm(zgn, 0.0, 1.0, false);

          const double info_frac =
              static_cast<double>(n1) / static_cast<double>(nend);

          double t1 = (d1e + d2e > 0) ? static_cast<double>(d1e) /
                                            static_cast<double>(d1e + d2e)
                                      : 0.5;
          t1 = std::min(std::max(t1, 1e-6), 1.0 - 1e-6);

          const size_t rankp = static_cast<size_t>(
              std::count_if(stg1_p.begin(), stg1_p.end(),
                            [&](double p) { return p < stg1_p[obd]; }));
          const double z1Selected = stg1_z[obd];

          auto combtest = [&](size_t m, const LocalPValues &stg1_loc_p) {
            PCStage2Result rej =
                fPCRejCpp(stg1_loc_p, stg2, noStg1Rej, stg2_elemhyp,
                          ALPHA_ONE_SIDED, info_frac);
            for (size_t k = 0; k < M; ++k)
              out.rej[m](n2i, k) = rej.rej_elem[k];
          };

          if (use[M_CTDUNNETT])
            combtest(M_CTDUNNETT, dunnett1);
          if (use[M_CTSIMES])
            combtest(M_CTSIMES, simes1);
          if (use[M_CTPOOLED])
            combtest(M_CTPOOLED, pooled1);

          if (use[M_PH3ONLY]) {
            if (p2 < ALPHA_ONE_SIDED)
              out.rej[M_PH3ONLY](n2i, obd) = 1;
          }

          if (use[M_NAIVE]) {
            if (p_cum < ALPHA_ONE_SIDED)
              out.rej[M_NAIVE](n2i, obd) = 1;
          }

          if (use[M_CER]) {
            // stage boundaries are precomputed once outside the loop since
            // rejection at the end of phase 2 is not allowed
            CER cer = fCERCerCpp(stg1_p, wgtmat, family, corr, t1, sb.stg1_bnd,
                                 sb.stg2_bnd);
            AdjustedBoundaries nb = fCERNewBoundCpp(
                stg1_p, wgtmat, family, corr, cer.stg1_inthyp_nr_idx, cer.CER,
                stg2_elemhyp, wgtmat1, t1);
            std::vector<int> rej_cer =
                fCERRejCpp(std::vector<double>{p_cum}, cer.stg1_elemhyp_r_idx,
                           stg2_elemhyp, nb.inthyp, nb.stg2_bnd_new);
            for (size_t k = 0; k < M; ++k)
              out.rej[M_CER](n2i, k) = rej_cer[k];
          }

          const bool needFullTsssdK = use[M_TSSSD_K] ||
                (use[M_TSSSD_K_MODRANK] && rankp == 0);
          const bool needFullTsssdUk =
              use[M_TSSSD_UK] || (use[M_TSSSD_UK_MODRANK] && rankp == 0);
          if (needFullTsssdK) {
            const double knownBound = final_selection_boundary(M, true, t1);
            if (use[M_TSSSD_K] && zgn > knownBound) {
              out.rej[M_TSSSD_K](n2i, obd) = 1;
            }
            if (use[M_TSSSD_K_MODRANK] && rankp == 0 && zgn > knownBound) {
              out.rej[M_TSSSD_K_MODRANK](n2i, obd) = 1;
            }
          }

          if (use[M_TSSSD_K_RANK]) {
            const double rankBound =
                rank_selection_boundary(M, rankp + 1, t1);
            if (zgn > rankBound)
              out.rej[M_TSSSD_K_RANK](n2i, obd) = 1;
          }

          if (use[M_TSSSD_K_MODRANK] && rankp > 0) {
            const double rankBound = M - rankp == 1
                           ? SINGLE_ARM_BOUND
                           : final_selection_boundary(M - rankp, true, t1);
            if (zgn > rankBound)
              out.rej[M_TSSSD_K_MODRANK](n2i, obd) = 1;
          }

          if (needFullTsssdUk) {
            const double unknownBound = final_selection_boundary(M, false, t1);
            if (use[M_TSSSD_UK] && zgn > unknownBound) {
              out.rej[M_TSSSD_UK](n2i, obd) = 1;
            }
            if (use[M_TSSSD_UK_MODRANK] && rankp == 0 && zgn > unknownBound) {
              out.rej[M_TSSSD_UK_MODRANK](n2i, obd) = 1;
            }
          }

          if (use[M_TSSSD_UK_MODRANK] && rankp > 0) {
            const double rankBound = M - rankp == 1
                           ? SINGLE_ARM_BOUND
                           : final_selection_boundary(M - rankp, false, t1);
            if (zgn > rankBound)
              out.rej[M_TSSSD_UK_MODRANK](n2i, obd) = 1;
          }

          if (needTsssdCe) {
            auto ceAdjustedBoundary = [&](double c2Nominal) {
              double sqrtTnom = std::sqrt(tNominal[n2i]);
              double sqrt1mTnom = std::sqrt(1.0 - tNominal[n2i]);
              double b1 = (c2Nominal - z1Selected * sqrtTnom) / sqrt1mTnom;
              return b1 * std::sqrt(1.0 - t1) + z1Selected * std::sqrt(t1);
            };

            if (use[M_TSSSD_K_CE]) {
              const double c2New =
                  ceAdjustedBoundary(tsssdKnownNomByEff[M - 1][n2i]);
              if (zgn > c2New)
                out.rej[M_TSSSD_K_CE](n2i, obd) = 1;
            }

            if (use[M_TSSSD_UK_CE]) {
              const double c2New =
                  ceAdjustedBoundary(tsssdUnknownNomByEff[M - 1][n2i]);
              if (zgn > c2New)
                out.rej[M_TSSSD_UK_CE](n2i, obd) = 1;
            }

            const size_t mEff = M - rankp;
            const size_t effIdx = mEff - 1;
            if (use[M_TSSSD_K_RANK_CE]) {
              const double c2New =
                  ceAdjustedBoundary(tsssdKnownNomByRank[rankp][n2i]);
              if (zgn > c2New)
                out.rej[M_TSSSD_K_RANK_CE](n2i, obd) = 1;
            }

            if (use[M_TSSSD_K_MODRANK_CE]) {
              const double c2New =
                  ceAdjustedBoundary(tsssdKnownNomByEff[effIdx][n2i]);
              if (zgn > c2New)
                out.rej[M_TSSSD_K_MODRANK_CE](n2i, obd) = 1;
            }

            if (use[M_TSSSD_UK_MODRANK_CE]) {
              const double c2New =
                  ceAdjustedBoundary(tsssdUnknownNomByEff[effIdx][n2i]);
              if (zgn > c2New)
                out.rej[M_TSSSD_UK_MODRANK_CE](n2i, obd) = 1;
            }
          }

          std::vector<int> reject(NMETHOD, 0);
          for (size_t m = 0; m < NMETHOD; ++m) {
            if (use[m])
              reject[m] = out.rej[m](n2i, obd);
          }
          out.tteSummary.push_back(
              {static_cast<int>(iter + 1), static_cast<int>(obd + 1),
               static_cast<int>(n2cur), d1e, d2e, d1e + d2e, stg1_z[obd], z2,
               zgn, std::move(reject)});
        }

        out.completed = true;
      } catch (const std::exception &e) {
        thread_utils::push_thread_warning(
            "iteration " + std::to_string(iter + 1) + ": " + e.what());
      }
    }
  }
};

} // anonymous namespace

// Parallel entry function
ListCpp lrsim_bmTrtSel_cpp(
    const size_t M, const size_t n1, const size_t n2min, const size_t n2max,
    const double p0, const std::vector<double> &pe,
    const std::vector<double> &pt, const double rho,
    const std::vector<double> the0, const FlatMatrix &the1, const double T_max,
    const double w, const double phi_t, const double ce, const double ct,
    const bool uniform_prior, const double acc_rate1, const double acc_rate2,
    const double T_ph2followup, const std::vector<std::string> &methods,
    const size_t ntrial, const size_t maxRawDatasets, const int seed) {

  if (M < 1)
    throw std::invalid_argument("M must be at least 1");
  if (n1 <= 0)
    throw std::invalid_argument("n1 must be positive");
  if (n2min <= 0)
    throw std::invalid_argument("n2min must be positive");
  if (n2max < n2min)
    throw std::invalid_argument("n2max must be >= n2min");
  if (p0 <= 0 || p0 >= 1)
    throw std::invalid_argument("p0 must lie in (0, 1)");
  if (pe.size() != M)
    throw std::invalid_argument("pe must have length M");
  if (pt.size() != M)
    throw std::invalid_argument("pt must have length M");
  for (double p : pe) {
    if (p <= 0 || p >= 1)
      throw std::invalid_argument("pe must lie in (0, 1)");
  }
  for (double p : pt) {
    if (p < 0 || p >= 1)
      throw std::invalid_argument("pt must lie in [0, 1)");
  }
  if (rho <= -1 || rho >= 1) {
    throw std::invalid_argument("rho must lie in (-1, 1)");
  }
  if (the0.size() != 2)
    throw std::invalid_argument("the0 must have length 2");
  if (the1.nrow != M || the1.ncol != 2)
    throw std::invalid_argument("the1 must have dimensions M x 2");
  for (double t : the0) {
    if (t <= 0)
      throw std::invalid_argument("the0 must be positive");
  }
  for (size_t j = 0; j < 2; ++j) {
    for (size_t i = 0; i < M; ++i) {
      if (the1(i, j) <= 0)
        throw std::invalid_argument("the1 must be positive");
    }
  }
  if (T_max <= 0)
    throw std::invalid_argument("T_max must be positive");
  if (w < 0)
    throw std::invalid_argument("w must be nonnegative");
  if (phi_t <= 0 || phi_t > 1)
    throw std::invalid_argument("phi_t must lie in (0, 1]");
  if (ce < 0 || ce >= 1)
    throw std::invalid_argument("ce must lie in [0, 1)");
  if (ct < 0 || ct >= 1)
    throw std::invalid_argument("ct must lie in [0, 1)");
  if (acc_rate1 <= 0)
    throw std::invalid_argument("acc_rate1 must be positive");
  if (acc_rate2 <= 0)
    throw std::invalid_argument("acc_rate2 must be positive");
  if (T_ph2followup < 0)
    throw std::invalid_argument("T_ph2followup must be nonegative");
  if (ntrial <= 0)
    throw std::invalid_argument("ntrial must be positive");

  // an empty methods vector requests every method
  std::vector<unsigned char> use(NMETHOD, methods.empty() ? 1 : 0);
  for (const std::string &s : methods) {
    size_t m = 0;
    for (; m < NMETHOD; ++m) {
      if (s == METHOD_NAME[m])
        break;
    }
    if (m == NMETHOD)
      throw std::invalid_argument("unknown method: " + s +
                                  "; use lowercase method names");
    use[m] = 1;
  }

  // generate seeds for each iteration to ensure reproducibility
  std::vector<uint64_t> seeds(ntrial);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < ntrial; ++iter)
    seeds[iter] = master_rng();

  const size_t ntr = M;
  const size_t ngrid = n2max - n2min + 1;

  // location of the true OBD rather than the one selected based on data
  size_t true_id = 0;
  for (size_t k = 1; k < ntr; ++k) {
    if (pe[k] - w * pt[k] > pe[true_id] - w * pt[true_id])
      true_id = k;
  }

  // Equal weights within each intersection hypothesis.
  WeightMatrix wgtmat = fDefaultWgtmatcpp(ntr);

  // single hypothesis retained at stage 2 after treatment selection
  WeightMatrix wgtmat1 = fDefaultWgtmatcpp(1);

  BoolMatrix family(1, ntr);
  family.fill(1);

  // correlation of the log-rank statistics under equal allocation
  FlatMatrix corr(ntr, ntr);
  for (size_t i = 0; i < ntr; ++i) {
    for (size_t j = 0; j < ntr; ++j)
      corr(i, j) = (i == j) ? 1.0 : 0.5;
  }

  // rejection at the end of phase 2 is not allowed, so the stage 1 boundary
  // is 0 and the stage 2 boundary does not depend on the information
  // fraction; precompute it once before data generation for the CER method
  StageBoundaries sb =
      fCERStageBoundCpp(wgtmat, family, corr, ALPHA_ONE_SIDED, 0.0, 0.5);

  std::vector<TrialResult> results(ntrial);
  SimWorker worker(M, n1, n2min, n2max, p0, pe, pt, rho, the0, the1, T_max, w,
                   phi_t, ce, ct, acc_rate1, acc_rate2, T_ph2followup,
                   maxRawDatasets, seeds, use, wgtmat, wgtmat1, family, corr,
                   sb, uniform_prior,
                   &results);
  RcppParallel::parallelFor(0, ntrial, worker);

  std::vector<int> select_count(ntr, 0);
  std::vector<IntMatrix> rej_each(NMETHOD);
  std::vector<std::vector<int>> rej_any(NMETHOD);
  for (size_t m = 0; m < NMETHOD; ++m) {
    if (!use[m])
      continue;
    rej_each[m].resize(ngrid, ntr);
    rej_any[m].assign(ngrid, 0);
  }
  IntMatrix total_events(ngrid, 3);
  std::vector<int> sumBinIter, sumBinTreatment, sumBinResponses,
      sumBinToxicities, sumBinSelected;
  std::vector<int> sumTteIter, sumTteSelectedDose, sumTtePhase3SampleSize,
      sumTteStage1Events, sumTteStage2Events, sumTteTotalEvents;
  std::vector<double> sumTteStage1Z, sumTteStage2Z, sumTteCumulativeZ;
  std::vector<std::vector<int>> sumTteReject(NMETHOD);
  std::vector<int> rawBinIter, rawBinSubject, rawBinTreatment, rawBinResponse,
      rawBinToxicity;
  std::vector<int> rawTteIter, rawTteSubject, rawTtePhase, rawTteTreatment,
      rawTteResponse, rawTteEvent;
  std::vector<double> rawTteArrival, rawTteSurvival, rawTteObserved;

  for (size_t iter = 0; iter < ntrial; ++iter) {
    const TrialResult &out = results[iter];
    for (const BinarySummaryRow &row : out.binarySummary) {
      sumBinIter.push_back(row.iteration);
      sumBinTreatment.push_back(row.treatment);
      sumBinResponses.push_back(row.responses);
      sumBinToxicities.push_back(row.toxicities);
      sumBinSelected.push_back(row.selected);
    }
    for (const TteSummaryRow &row : out.tteSummary) {
      sumTteIter.push_back(row.iteration);
      sumTteSelectedDose.push_back(row.selectedDose);
      sumTtePhase3SampleSize.push_back(row.phase3SampleSize);
      sumTteStage1Events.push_back(row.stage1Events);
      sumTteStage2Events.push_back(row.stage2Events);
      sumTteTotalEvents.push_back(row.totalEvents);
      sumTteStage1Z.push_back(row.stage1LogRankZ);
      sumTteStage2Z.push_back(row.stage2LogRankZ);
      sumTteCumulativeZ.push_back(row.cumulativeLogRankZ);
      for (size_t m = 0; m < NMETHOD; ++m) {
        if (use[m])
          sumTteReject[m].push_back(row.reject[m]);
      }
    }
    if (maxRawDatasets > 0) {
      for (const RawBinaryRow &row : out.rawBinary) {
        rawBinIter.push_back(row.iteration);
        rawBinSubject.push_back(row.subject);
        rawBinTreatment.push_back(row.treatment);
        rawBinResponse.push_back(row.response);
        rawBinToxicity.push_back(row.toxicity);
      }
      for (const RawTteRow &row : out.rawTte) {
        rawTteIter.push_back(row.iteration);
        rawTteSubject.push_back(row.subject);
        rawTtePhase.push_back(row.phase);
        rawTteTreatment.push_back(row.treatment);
        rawTteResponse.push_back(row.response);
        rawTteArrival.push_back(row.arrival);
        rawTteSurvival.push_back(row.survival);
        rawTteObserved.push_back(row.observed);
        rawTteEvent.push_back(row.event);
      }
    }

    if (!out.completed)
      continue;

    for (size_t k = 0; k < ntr; ++k)
      select_count[k] += out.select[k];

    for (size_t m = 0; m < NMETHOD; ++m) {
      if (!use[m])
        continue;
      for (size_t n2i = 0; n2i < ngrid; ++n2i) {
        for (size_t k = 0; k < ntr; ++k) {
          rej_each[m](n2i, k) += out.rej[m](n2i, k);
          // only the selected dose enters phase 3, so at most one column is
          // nonzero and the row sum is the any-rejection indicator
          rej_any[m][n2i] += out.rej[m](n2i, k);
        }
      }
    }

    for (size_t n2i = 0; n2i < ngrid; ++n2i) {
      for (size_t j = 0; j < 3; ++j)
        total_events(n2i, j) += out.events(n2i, j);
    }
  }

  const double dntrial = static_cast<double>(ntrial);

  std::vector<double> selectProb(ntr);
  for (size_t k = 0; k < ntr; ++k) {
    selectProb[k] = select_count[k] / dntrial;
  }

  // if two or more doses have identical benefit-risk, the first one is taken
  // as the true OBD, in which case this quantity is not meaningful
  double pcs = selectProb[true_id] * 100.0;

  // average number of events in stage 1, stage 2, and combined for the
  // selected dose plus the control arm
  FlatMatrix ave_event(ngrid, 3);
  for (size_t n2i = 0; n2i < ngrid; ++n2i) {
    for (size_t j = 0; j < 3; ++j) {
      ave_event(n2i, j) = std::round(total_events(n2i, j) / dntrial);
    }
  }

  // generalized power for the true OBD, rejection probability for each dose
  // conditional on it being selected, and probability of rejecting any dose
  auto summarize = [&](const IntMatrix &rej_each,
                       const std::vector<int> &rej_any,
                       std::vector<double> &gpower, FlatMatrix &prob_each,
                       std::vector<double> &prob_any) {
    gpower.resize(ngrid);
    prob_each.resize(ngrid, ntr);
    prob_any.resize(ngrid);
    for (size_t n2i = 0; n2i < ngrid; ++n2i) {
      gpower[n2i] = rej_each(n2i, true_id) / dntrial;
      for (size_t k = 0; k < ntr; ++k) {
        prob_each(n2i, k) = (selectProb[k] > 0.0)
                                ? rej_each(n2i, k) / dntrial / selectProb[k]
                                : NaN;
      }
      prob_any[n2i] = rej_any[n2i] / dntrial;
    }
  };

  std::vector<size_t> n2(ngrid);
  for (size_t n2i = 0; n2i < ngrid; ++n2i)
    n2[n2i] = n2min + n2i;

  // one entry per requested method, in the canonical method order
  std::vector<std::string> methodNames;
  ListCpp byMethod;
  for (size_t m = 0; m < NMETHOD; ++m) {
    if (!use[m])
      continue;

    std::vector<double> gpower, prob_rej_any;
    FlatMatrix prob_rej_each;
    summarize(rej_each[m], rej_any[m], gpower, prob_rej_each, prob_rej_any);

    ListCpp res;
    res.push_back(std::move(gpower), "gpower");
    res.push_back(std::move(prob_rej_each), "prob.rej.each");
    res.push_back(std::move(prob_rej_any), "prob.rej.any");

    methodNames.push_back(METHOD_NAME[m]);
    byMethod.push_back(std::move(res), METHOD_NAME[m]);
  }

  ListCpp result;
  DataFrameCpp sumdataBIN;
  sumdataBIN.push_back(std::move(sumBinIter), "iterationNumber");
  sumdataBIN.push_back(std::move(sumBinTreatment), "treatmentGroup");
  sumdataBIN.push_back(std::move(sumBinResponses), "responses");
  sumdataBIN.push_back(std::move(sumBinToxicities), "toxicities");
  sumdataBIN.push_back(std::move(sumBinSelected), "selected");

  DataFrameCpp sumdataTTE;
  sumdataTTE.push_back(std::move(sumTteIter), "iterationNumber");
  sumdataTTE.push_back(std::move(sumTteSelectedDose), "selectedDose");
  sumdataTTE.push_back(std::move(sumTtePhase3SampleSize),
                        "phase3SampleSize");
  sumdataTTE.push_back(std::move(sumTteStage1Events), "stage1Events");
  sumdataTTE.push_back(std::move(sumTteStage2Events), "stage2Events");
  sumdataTTE.push_back(std::move(sumTteTotalEvents), "totalEvents");
  sumdataTTE.push_back(std::move(sumTteStage1Z), "stage1LogRankZ");
  sumdataTTE.push_back(std::move(sumTteStage2Z), "stage2LogRankZ");
  sumdataTTE.push_back(std::move(sumTteCumulativeZ), "cumulativeLogRankZ");
  for (size_t m = 0; m < NMETHOD; ++m) {
    if (use[m]) {
      sumdataTTE.push_back(std::move(sumTteReject[m]),
                            std::string("reject.") + METHOD_NAME[m]);
    }
  }

  result.push_back(n1, "n1");
  result.push_back(std::move(n2), "n2");
  result.push_back(ntrial, "numberOfIterations");
  result.push_back(true_id + 1, "trueOBD");
  result.push_back(std::move(selectProb), "selectionProb");
  result.push_back(pcs, "pcs");
  result.push_back(std::move(ave_event), "ave.event");
  result.push_back(std::move(methodNames), "methods");
  result.push_back(std::move(byMethod), "byMethod");
  result.push_back(std::move(sumdataBIN), "sumdataBIN");
  result.push_back(std::move(sumdataTTE), "sumdataTTE");
  if (maxRawDatasets > 0) {
    DataFrameCpp rawBinary;
    rawBinary.push_back(std::move(rawBinIter), "iterationNumber");
    rawBinary.push_back(std::move(rawBinSubject), "subjectId");
    rawBinary.push_back(std::move(rawBinTreatment), "treatmentGroup");
    rawBinary.push_back(std::move(rawBinResponse), "response");
    rawBinary.push_back(std::move(rawBinToxicity), "toxicity");

    DataFrameCpp rawTte;
    rawTte.push_back(std::move(rawTteIter), "iterationNumber");
    rawTte.push_back(std::move(rawTteSubject), "subjectId");
    rawTte.push_back(std::move(rawTtePhase), "phase");
    rawTte.push_back(std::move(rawTteTreatment), "treatmentGroup");
    rawTte.push_back(std::move(rawTteResponse), "response");
    rawTte.push_back(std::move(rawTteArrival), "arrivalTime");
    rawTte.push_back(std::move(rawTteSurvival), "survivalTime");
    rawTte.push_back(std::move(rawTteObserved), "timeUnderObservation");
    rawTte.push_back(std::move(rawTteEvent), "event");

    result.push_back(std::move(rawBinary), "rawdataBIN");
    result.push_back(std::move(rawTte), "rawdataTTE");
  }

  return result;
}
// [[Rcpp::export]]
Rcpp::List lrsim_bmTrtSel_Rcpp(
    const int phase2SampleSizePerArm = NA_INTEGER,
    const int phase3SampleSizePerArmMin = NA_INTEGER,
    const int phase3SampleSizePerArmMax = NA_INTEGER,
    const double responseProbControl = NA_REAL,
    const Rcpp::NumericVector &responseProbTreatments = NA_REAL,
    const Rcpp::NumericVector &toxicityProbTreatments = NA_REAL,
    const double corrEfficacyToxicity = 0.5,
    const Rcpp::NumericVector &hazardRateControl = NA_REAL,
    const Rcpp::NumericMatrix &hazardRateTreatments = Rcpp::NumericMatrix(),
    const double studyDurationPhase3 = NA_REAL,
    const double toxicityWeight = NA_REAL,
    const double toxicityUpperLimit = NA_REAL,
    const double efficacyThreshold = 0,
    const double safetyThreshold = 0,
    const bool useUniformPrior = true,
    const Rcpp::Nullable<Rcpp::CharacterVector> methods = R_NilValue,
    const double accrualRatePhase2 = NA_REAL,
    const double accrualRatePhase3 = NA_REAL,
    const double followupTimePhase2 = 0,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasets = 0,
    const int seed = 0) {

  if (maxNumberOfRawDatasets < 0) {
    throw std::invalid_argument(
        "maxNumberOfRawDatasets must be a non-negative integer");
  }
  if (maxNumberOfRawDatasets > maxNumberOfIterations) {
    throw std::invalid_argument(
        "maxNumberOfRawDatasets cannot exceed maxNumberOfIterations");
  }

  std::vector<double> pe(responseProbTreatments.begin(),
                         responseProbTreatments.end());
  std::vector<double> pt(toxicityProbTreatments.begin(),
                         toxicityProbTreatments.end());
  std::vector<double> the0(hazardRateControl.begin(), hazardRateControl.end());
  FlatMatrix the1 = flatmatrix_from_Rmatrix(hazardRateTreatments);

  std::vector<std::string> methodVec;
  if (methods.isNotNull()) {
    methodVec =
        Rcpp::as<std::vector<std::string>>(Rcpp::CharacterVector(methods));
  }

  auto out = lrsim_bmTrtSel_cpp(
      pe.size(), static_cast<size_t>(phase2SampleSizePerArm),
      static_cast<size_t>(phase3SampleSizePerArmMin),
      static_cast<size_t>(phase3SampleSizePerArmMax), responseProbControl, pe,
      pt, corrEfficacyToxicity, the0, the1, studyDurationPhase3, toxicityWeight,
      toxicityUpperLimit, efficacyThreshold, safetyThreshold, useUniformPrior,
      accrualRatePhase2, accrualRatePhase3, followupTimePhase2, methodVec,
      static_cast<size_t>(maxNumberOfIterations),
      static_cast<size_t>(maxNumberOfRawDatasets), seed);

  thread_utils::drain_thread_warnings_to_R();

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "lrsim_bmTrtSel";

  return result;
}
