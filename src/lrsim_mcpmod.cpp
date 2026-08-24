#include "survival_analysis.h"
#include "mvnormr.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/random.hpp>

using std::size_t;


ListCpp lrsim_mcpmod_cpp(
    const size_t M,
    const double alpha,
    const std::vector<double>& hazardRatioH0s,
    const std::vector<double>& allocations,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<std::vector<double>>& lambdas,
    const FlatMatrix& candidateHazardRatios,
    const std::vector<std::vector<double>>& gammas,
    const size_t n,
    const double followupTime,
    const bool fixedFollowup,
    const int plannedEvents,
    const double plannedTime,
    const int maxNumberOfIterations,
    const int maxNumberOfRawDatasetsPerStage,
    const int seed) {

  if (M < 2)
    throw std::invalid_argument("M must be at least 2");

  // decide planning mode
  bool useEvents;
  if (plannedEvents != INT_MIN) {
    useEvents = true;
    if (plannedEvents <= 0)
      throw std::invalid_argument("plannedEvents must be positive");
  } else if (!std::isnan(plannedTime)) {
    useEvents = false;
    if (plannedTime <= 0.0)
      throw std::invalid_argument("plannedTime must be positive");
  } else {
    throw std::invalid_argument(
        "Either plannedEvents or plannedTime must be given");
  }

  // validate other input parameters
  std::vector<double> hrH0s = expand1(hazardRatioH0s, M, "hazardRatioH0s");
  if (std::any_of(hrH0s.begin(), hrH0s.end(),
                  [](double v) { return v <= 0.0; }))
    throw std::invalid_argument(
        "All hazardRatioH0 parameters must be positive");
  std::vector<double> allocs = expand1(allocations, M + 1, "allocations");
  if (std::any_of(allocs.begin(), allocs.end(), [](double v) { return v < 1; }))
    throw std::invalid_argument("All allocation parameters must be positive");
  if (accrualTime[0] != 0.0)
    throw std::invalid_argument("accrualTime must start with 0");
  if (any_nonincreasing(accrualTime))
    throw std::invalid_argument("accrualTime should be increasing");
  if (!none_na(accrualIntensity))
    throw std::invalid_argument("accrualIntensity must be provided");
  if (accrualIntensity.size() != accrualTime.size())
    throw std::invalid_argument("Invalid length for accrualIntensity");
  for (double v : accrualIntensity) {
    if (v < 0.0)
      throw std::invalid_argument("accrualIntensity must be non-negative");
  }
  if (piecewiseSurvivalTime[0] != 0.0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");
  for (double v : stratumFraction) {
    if (v <= 0.0)
      throw std::invalid_argument("stratumFraction must be positive");
  }
  double sumf =
    std::accumulate(stratumFraction.begin(), stratumFraction.end(), 0.0);
  if (std::fabs(sumf - 1.0) > 1e-12)
    throw std::invalid_argument("stratumFraction must sum to 1");

  for (size_t m = 0; m < M + 1; ++m) {
    std::string nm = std::string("lambdas[") + std::to_string(m) + "]";
    if (!none_na(lambdas[m]))
      throw std::invalid_argument(nm + "must be provided");
    for (double v : lambdas[m]) {
      if (v < 0.0)
        throw std::invalid_argument(nm + "must be non-negative");
    }
  }
  for (size_t m = 0; m < M + 1; ++m) {
    std::string nm = std::string("gammas[") + std::to_string(m) + "]";
    if (!none_na(gammas[m]))
      throw std::invalid_argument(nm + "must be provided");
    for (double v : gammas[m]) {
      if (v < 0.0)
        throw std::invalid_argument(nm + "must be non-negative");
    }
  }

  size_t T = candidateHazardRatios.ncol; // number of candidate models
  if (candidateHazardRatios.nrow != M) {
    throw std::invalid_argument(
        "Rows of candidateHazardRatios must match active arms");
  }
  for (double x : candidateHazardRatios.data) {
    if (x <= 0.0) {
      throw std::invalid_argument("candidateHazardRatios must be postive");
    }
  }

  // log hazard ratios under H0 for each arm
  std::vector<double> betaH0(M);
  for (size_t m = 0; m < M; ++m) {
    betaH0[m] = std::log(hrH0s[m]);
  }

  // candidate log hazard ratios for each arm
  std::vector<std::vector<double>> candidateBetas(T);
  for (size_t t = 0; t < T; ++t) {
    candidateBetas[t].resize(M);
    for (size_t m = 0; m < M; ++m) {
      candidateBetas[t][m] = std::log(candidateHazardRatios(m, t));
    }
  }

  if (static_cast<int>(n) == INT_MIN)
    throw std::invalid_argument("n must be provided");
  if (n <= 0)
    throw std::invalid_argument("n must be a positive integer");
  if (fixedFollowup && std::isnan(followupTime))
    throw std::invalid_argument(
        "followupTime must be provided for fixed follow-up");
  if (fixedFollowup && followupTime <= 0.0)
    throw std::invalid_argument(
        "followupTime must be positive for fixed follow-up");
  if (maxNumberOfIterations < 1)
    throw std::invalid_argument(
        "maxNumberOfIterations must be a positive integer");
  if (maxNumberOfRawDatasetsPerStage < 0)
    throw std::invalid_argument(
        "maxNumberOfRawDatasetsPerStage must be a non-negative integer");
  if (maxNumberOfRawDatasetsPerStage > maxNumberOfIterations)
    throw std::invalid_argument(
        "maxNumberOfRawDatasetsPerStage cannot exceed maxNumberOfIterations");

  size_t maxIters = static_cast<size_t>(maxNumberOfIterations);
  size_t maxRawIters = static_cast<size_t>(maxNumberOfRawDatasetsPerStage);
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  const std::vector<double>& tau = piecewiseSurvivalTime;
  const double fu = followupTime;

  // expand stratified inputs
  FlatArray lambdasx(nintv, nstrata, M + 1);
  for (size_t m = 0; m < M + 1; ++m) {
    std::string nm = std::string("lambdas[") + std::to_string(m) + "]";
    expand_stratified_to_slice(lambdas[m], lambdasx, m, nstrata, nintv,
                               nm.c_str());
  }

  FlatArray gammasx(nintv, nstrata, M + 1);
  for (size_t m = 0; m < M + 1; ++m) {
    std::string nm = std::string("gammas[") + std::to_string(m) + "]";
    expand_stratified_to_slice(gammas[m], gammasx, m, nstrata, nintv,
                               nm.c_str());
  }

  // construct covariance matrix of log hazard ratio estimates under H0
  const double sumAlloc = std::accumulate(allocs.begin(), allocs.end(), 0.0);
  std::vector<double> p(M + 1);
  for (size_t m = 0; m < M + 1; ++m) {
    p[m] = allocs[m] / sumAlloc;
  }
  FlatMatrix S0(M, M);
  S0.fill(1.0 / p[M]);
  for (size_t m = 0; m < M; ++m) {
    S0(m, m) = 1.0 / p[m] + 1.0 / p[M];
  }

  // generate seeds for each iteration to ensure reproducibility
  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter)
    seeds[iter] = master_rng();

  // One summary (stage-level) row produced by an iteration
  struct Summary1Row {
    int iterNum = 0;
    unsigned char evNotAch = 0;
    double analysisT = 0.0;
    int trtGrp = 0;
    int accruals = 0, events = 0, dropouts = 0;
  };

  struct Summary2Row {
    int iterNum = 0;
    double analysisT = 0.0;
    std::vector<double> beta;
    FlatMatrix vbeta;
    FlatMatrix c;
    FlatMatrix vzH0;
    double q = 0.0;
    std::vector<double> z;
    int optModel = 0;
    bool reject = 0;
  };

  // One subject-level (raw) row for a particular iteration and stage
  struct RawDatasetRow {
    int iterNum = 0;
    double analysisT = 0.0;
    int subjectId = 0;
    double arrivalT = 0.0;
    int stratum = 0, trtGrp = 0;
    double survivalT = 0.0, dropoutT = 0.0, timeObs = 0.0;
    unsigned char event = 0, dropEv = 0;
  };

  // Per-iteration container written exclusively by the worker thread
  struct IterationResult {
    std::vector<Summary1Row> summary1Rows;
    std::vector<Summary2Row> summary2Rows;
    std::vector<RawDatasetRow> rawRows;
    void reserveForSummary1(size_t approxRows) {
      summary1Rows.reserve(approxRows);
    }
    void reserveForSummary2(size_t approxRows) {
      summary2Rows.reserve(approxRows);
    }
    void reserveForRaw(size_t approxRows) { rawRows.reserve(approxRows); }
  };

  // pre-size per-iteration results
  std::vector<IterationResult> results;
  results.resize(maxIters);


  // Worker that runs simulation iterations [begin, end)
  struct SimWorker : public RcppParallel::Worker {
    // inputs (const refs)
    const size_t M;
    const double alpha;
    const std::vector<double>& betaH0;
    const std::vector<double>& allocs;
    const std::vector<double>& accrualTime;
    const std::vector<double>& accrualIntensity;
    const std::vector<double>& tau;
    const std::vector<double>& stratumFraction;
    const FlatArray& lambdasx;
    const std::vector<std::vector<double>>& candidateBetas;
    const FlatArray& gammasx;
    const size_t n;
    const double fu;
    const bool fixedFollowup;
    const int plannedEvents;
    const double plannedTime;
    const size_t maxRawIters; // store raw for iter < maxRawIters
    const std::vector<uint64_t>& seeds;
    const bool useEvents;
    const size_t nintv;
    const size_t nstrata;
    const size_t T;
    const FlatMatrix& S0;

    // output pointer (pre-sized vector of IterationResult)
    std::vector<IterationResult>* results;

    SimWorker(
      size_t M_,
      double alpha_,
      const std::vector<double>& betaH0_,
      const std::vector<double>& allocs_,
      const std::vector<double>& accrualTime_,
      const std::vector<double>& accrualIntensity_,
      const std::vector<double>& tau_,
      const std::vector<double>& stratumFraction_,
      const FlatArray& lambdasx_,
      const std::vector<std::vector<double>>& candidateBetas_,
      const FlatArray& gammasx_,
      size_t n_,
      double fu_,
      bool fixedFollowup_,
      const int plannedEvents_,
      const double plannedTime_,
      size_t maxRawIters_,
      const std::vector<uint64_t>& seeds_,
      bool useEvents_,
      size_t nintv_,
      size_t nstrata_,
      size_t T_,
      const FlatMatrix& S0_,
      std::vector<IterationResult>* results_)
      : M(M_),
        alpha(alpha_),
        betaH0(betaH0_),
        allocs(allocs_),
        accrualTime(accrualTime_),
        accrualIntensity(accrualIntensity_),
        tau(tau_),
        stratumFraction(stratumFraction_),
        lambdasx(lambdasx_),
        candidateBetas(candidateBetas_),
        gammasx(gammasx_),
        n(n_),
        fu(fu_),
        fixedFollowup(fixedFollowup_),
        plannedEvents(plannedEvents_),
        plannedTime(plannedTime_),
        maxRawIters(maxRawIters_),
        seeds(seeds_),
        useEvents(useEvents_),
        nintv(nintv_),
        nstrata(nstrata_),
        T(T_),
        S0(S0_),
        results(results_)
    {}

    void operator()(std::size_t begin, std::size_t end) {

      // local buffers reused by this worker
      size_t M1 = M + 1;
      size_t M2 = M + 2;
      std::vector<int> stratum(n), trtGrp(n);
      std::vector<double> arrivalT(n), survivalT(n), dropoutT(n);
      std::vector<double> timeObs(n), totalT(n);
      std::vector<unsigned char> event(n), dropEv(n);
      FlatMatrix bs(M1, nstrata);
      IntMatrix ns(M2, nstrata);
      std::vector<int> events(M2), dropouts(M2), accruals(M2);
      std::vector<double> denom_per_stratum(nstrata);
      std::vector<double> cumF(nstrata);
      std::partial_sum(stratumFraction.begin(), stratumFraction.end(),
                       cumF.begin());

      std::vector<double> totalte;
      totalte.reserve(n);
      std::vector<size_t> sub;
      sub.reserve(n);

      std::vector<double> init(1, NaN);
      std::vector<std::string> covariates(M);
      for (size_t m = 0; m < M; ++m) {
        covariates[m] = "arm_" + std::to_string(m + 1);
      }
      std::vector<double> treated;
      treated.reserve(n);

      std::vector<double> beta(M), betaH1(M);
      FlatMatrix vbeta(M, M), ibeta(M, M);

      // optimal contrasts for candidate models
      std::vector<std::vector<double>> c(T);
      FlatMatrix cmat(M, T);
      std::vector<double> temp(T);
      std::vector<double> zero(T);
      std::vector<double> z(T);
      FlatMatrix cScH0(T, T);
      FlatMatrix vzH0(T, T);

      const double sumAlloc =
        std::accumulate(allocs.begin(), allocs.end(), 0.0);

      for (size_t iter = begin; iter < end; ++iter) {

        // deterministic per-iteration RNG
        boost::random::mt19937_64 rng_local(seeds[iter]);
        boost::random::uniform_real_distribution<double> unif(0.0, 1.0);

        // per-iteration output container
        IterationResult& out = (*results)[iter];
        out.summary1Rows.clear();
        out.summary2Rows.clear();
        out.rawRows.clear();
        if (iter < maxRawIters) out.reserveForRaw(n);
        out.reserveForSummary1(M2); // all arms and overall
        out.reserveForSummary2(1);  // all pairwise comparisons with control

        std::fill(denom_per_stratum.begin(), denom_per_stratum.end(), sumAlloc);

        // reset block randomization
        for (size_t h = 0; h < nstrata; ++h) {
          flatmatrix_set_column(bs, h, allocs);
        }

        double analysisT = 0.0;
        double enrollt = 0.0;

        // generate cohort
        for (size_t i = 0; i < n; ++i) {
          double u = unif(rng_local);
          enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt);
          arrivalT[i] = enrollt;

          u = unif(rng_local);
          size_t j = findInterval1(u, cumF);
          stratum[i] = static_cast<int>(j + 1);

          // stratified block randomization among M + 1 arms
          u = unif(rng_local);
          auto b = flatmatrix_get_column_view(bs, j);
          double running = 0.0;
          size_t k = 0;
          for (size_t m = 0; m < M1; ++m) {
            running += b[m] / denom_per_stratum[j];
            if (u < running) {
              k = m;
              break;
            }
          }
          trtGrp[i] = static_cast<int>(k + 1);
          --bs(k, j); // decrement block count for this arm and stratum
          --denom_per_stratum[j];

          if (denom_per_stratum[j] <= 0) {
            flatmatrix_set_column(bs, j, allocs);
            denom_per_stratum[j] = sumAlloc;
          }

          // get lambda and gamma for this subject's stratum and arm
          size_t offset = FlatArray::idx(0, j, k, nintv, nstrata);
          const double* lamsrc = lambdasx.data_ptr() + offset;
          auto lam = DoubleView{lamsrc, nintv};
          const double* gamsrc = gammasx.data_ptr() + offset;
          auto gam = DoubleView{gamsrc, nintv};

          // survival time
          u = unif(rng_local);
          survivalT[i] = qtpwexpcpp1(u, tau, lam);

          // dropout time
          u = unif(rng_local);
          dropoutT[i] = qtpwexpcpp1(u, tau, gam);

          // initial observed time and event indicator
          double sv = survivalT[i], dr = dropoutT[i];
          if (fixedFollowup) {
            if (sv <= dr && sv <= fu) {
              timeObs[i] = sv;
              event[i] = 1;
              dropEv[i] = 0;
            } else if (dr <= sv && dr <= fu) {
              timeObs[i] = dr;
              event[i] = 0;
              dropEv[i] = 1;
            } else {
              timeObs[i] = fu;
              event[i] = 0;
              dropEv[i] = 0;
            }
          } else {
            if (sv <= dr) {
              timeObs[i] = sv;
              event[i] = 1;
              dropEv[i] = 0;
            } else {
              timeObs[i] = dr;
              event[i] = 0;
              dropEv[i] = 1;
            }
          }
          totalT[i] = arrivalT[i] + timeObs[i];
        } // cohort generated

        // determine analysis times (overall event count for all arms combined)
        bool evNotAch = false;

        if (useEvents) {
          totalte.clear();
          int nevents = 0; // all events across all arms for this iteration
          for (size_t i = 0; i < n; ++i) {
            if (event[i]) {
              ++nevents;
              totalte.push_back(totalT[i]);
            }
          }
          if (nevents == 0) {
            thread_utils::push_thread_warning(
              std::string("No events for iteration ") +
                std::to_string(iter + 1) + " skipping this iteration.");
            // leave out.summaryRows empty to signal skipped iteration
            out.summary1Rows.clear();
            out.summary2Rows.clear();
            out.rawRows.clear();
            continue;
          }
          std::sort(totalte.begin(), totalte.end());

          if (plannedEvents >= nevents) {
            analysisT = totalte.back() + 1e-12;
          } else {
            analysisT = totalte[plannedEvents - 1] + 1e-12;
          }
          evNotAch = (nevents < plannedEvents);
        } else {
          analysisT = plannedTime;
          evNotAch = false;
        }

        double time = analysisT;

        // reset counts
        ns.fill(0);
        std::fill(events.begin(), events.end(), 0);
        std::fill(dropouts.begin(), dropouts.end(), 0);

        // censoring at analysis time and count accruals/events/dropouts
        for (size_t i = 0; i < n; ++i) {
          double ar = arrivalT[i], sv = survivalT[i], dr = dropoutT[i];

          if (ar > time) {
            timeObs[i] = 0.0;
            event[i] = 0;
            dropEv[i] = 0;
            continue;
          }

          if (fixedFollowup) {
            if (ar + sv <= time && sv <= dr && sv <= fu) {
              timeObs[i] = sv;
              event[i] = 1;
              dropEv[i] = 0;
            } else if (ar + dr <= time && dr <= sv && dr <= fu) {
              timeObs[i] = dr;
              event[i] = 0;
              dropEv[i] = 1;
            } else if (ar + fu <= time && fu <= sv && fu <= dr) {
              timeObs[i] = fu;
              event[i] = 0;
              dropEv[i] = 0;
            } else {
              timeObs[i] = time - ar;
              event[i] = 0;
              dropEv[i] = 0;
            }
          } else {
            if (ar + sv <= time && sv <= dr) {
              timeObs[i] = sv;
              event[i] = 1;
              dropEv[i] = 0;
            } else if (ar + dr <= time && dr <= sv) {
              timeObs[i] = dr;
              event[i] = 0;
              dropEv[i] = 1;
            } else {
              timeObs[i] = time - ar;
              event[i] = 0;
              dropEv[i] = 0;
            }
          }

          size_t h = static_cast<size_t>(stratum[i] - 1);
          size_t m = static_cast<size_t>(trtGrp[i] - 1);

          ++ns(m, h);
          ++ns(M1, h); // total
          if (event[i]) {
            ++events[m];
            ++events[M1]; // total
          } else if (dropEv[i]) {
            ++dropouts[m];
            ++dropouts[M1]; // total
          }
        }

        std::fill(accruals.begin(), accruals.end(), 0);
        for (size_t h = 0; h < nstrata; ++h) {
          for (size_t m = 0; m < M2; ++m) {
            accruals[m] += ns(m, h);
          }
        }

        // collect indices with positive observed time
        sub.clear();
        for (size_t i = 0; i < n; ++i) {
          if (timeObs[i] > 0.0) sub.push_back(i);
        }

        // compute stratified Cox for pairwise comparisons vs. control
        DataFrameCpp data_outcome;
        data_outcome.push_back(subset(stratum, sub), "stratum");
        data_outcome.push_back(subset(timeObs, sub), "time");
        data_outcome.push_back(subset(event, sub), "event");

        // create dummy variables for treatment arms (control is reference)
        size_t nsub = sub.size();
        treated.resize(nsub);
        for (size_t m = 0; m < M; ++m) {
          std::string nm = covariates[m];
          for (size_t i = 0; i < nsub; ++i) {
            size_t idx = sub[i];
            treated[i] = (trtGrp[idx] == static_cast<int>(m + 1) ? 1.0 : 0.0);
          }
          data_outcome.push_back(treated, nm);
        }

        ListCpp fit_outcome =
          phregcpp(data_outcome, {"stratum"}, "time", "", "event", covariates,
                   "", "", "", "efron", init, 0, 0, 0, 0, 0, alpha);

        DataFrameCpp parest = fit_outcome.get<DataFrameCpp>("parest");
        beta = parest.get<double>("beta");
        vbeta = fit_outcome.get<FlatMatrix>("vbeta");
        ibeta = invsympd(vbeta, M);

        // optimal contrasts for candidate models
        for (size_t t = 0; t < T; ++t) {
          c[t] = mat_vec_mult(ibeta, candidateBetas[t]);
        }

        for (size_t t2 = 0; t2 < T; ++t2) {
          temp = mat_vec_mult(S0, c[t2]);
          for (size_t t1 = 0; t1 < T; ++t1) {
            cScH0(t1, t2) = std::inner_product(c[t1].begin(), c[t1].end(),
                  temp.begin(), 0.0);
          }
        }

        for (size_t t2 = 0; t2 < T; ++t2) {
          for (size_t t1 = 0; t1 < T; ++t1) {
            vzH0(t1, t2) =
              cScH0(t1, t2) / std::sqrt(cScH0(t1, t1) * cScH0(t2, t2));
          }
        }

        // critical value for max test
        double q = 1.0 - qmvnormcpp(1.0 - alpha, zero, vzH0);

        // construct the z test statistics
        for (size_t m = 0; m < M; ++m) {
          betaH1[m] = beta[m] - betaH0[m];
        }

        for (size_t t = 0; t < T; ++t) {
          auto temp = mat_vec_mult(vbeta, c[t]);
          double denom =
            std::inner_product(c[t].begin(), c[t].end(), temp.begin(), 0.0);
          z[t] = std::inner_product(c[t].begin(), c[t].end(), betaH1.begin(),
                                    0.0) /
                                      std::sqrt(denom);
        }

        int optModel =
          std::distance(z.begin(), std::max_element(z.begin(), z.end()));
        double maxz = z[optModel];
        bool reject = (maxz >= q);

        // optionally append raw rows for this stage
        if (iter < maxRawIters) { // only for first maxRawIters iterations
          for (size_t i = 0; i < n; ++i) {
            // skip subjects who haven't been enrolled by analysis time
            if (arrivalT[i] > time) continue;
            RawDatasetRow rr;
            rr.iterNum = static_cast<int>(iter + 1);
            rr.analysisT = time;
            rr.subjectId = static_cast<int>(i + 1);
            rr.arrivalT = arrivalT[i];
            rr.stratum = stratum[i];
            rr.trtGrp = trtGrp[i];
            rr.survivalT = survivalT[i];
            rr.dropoutT = dropoutT[i];
            rr.timeObs = timeObs[i];
            rr.event = event[i];
            rr.dropEv = dropEv[i];
            out.rawRows.push_back(std::move(rr));
          }
        }

        // append summary row for this stage
        Summary1Row sr1;
        for (size_t m = 0; m < M2; ++m) {
          sr1.iterNum = static_cast<int>(iter + 1);
          sr1.evNotAch = evNotAch ? 1 : 0;
          sr1.analysisT = time;
          sr1.trtGrp = static_cast<int>(m + 1);
          sr1.accruals = accruals[m];
          sr1.events = events[m];
          sr1.dropouts = dropouts[m];
          out.summary1Rows.push_back(sr1);
        }

        for (size_t t = 0; t < T; ++t) {
          for (size_t m = 0; m < M; ++m) {
            cmat(m, t) = c[t][m];
          }
        }

        Summary2Row sr2;
        sr2.iterNum = static_cast<int>(iter + 1);
        sr2.analysisT = time;
        sr2.beta = beta;
        sr2.vbeta = vbeta;
        sr2.c = cmat;
        sr2.vzH0 = vzH0;
        sr2.q = q;
        sr2.z = z;
        sr2.optModel = static_cast<int>(optModel + 1);
        sr2.reject = reject;
        out.summary2Rows.push_back(sr2);
      } // per-iteration
    } // operator()
  }; // SimWorker

  // run worker in parallel
  SimWorker worker(M, alpha, betaH0, allocs, accrualTime, accrualIntensity, tau,
                   stratumFraction, lambdasx, candidateBetas, gammasx, n, fu,
                   fixedFollowup, plannedEvents, plannedTime, maxRawIters,
                   seeds, useEvents, nintv, nstrata, T, S0, &results);

  RcppParallel::parallelFor(0, maxIters, worker);

  // Flatten results
  size_t ns1r = 0, ns2r = 0, nrr = 0;
  for (size_t iter = 0; iter < maxIters; ++iter) {
    ns1r += results[iter].summary1Rows.size();
    ns2r += results[iter].summary2Rows.size();
    nrr += results[iter].rawRows.size();
  }
  if (ns1r == 0 || ns2r == 0)
    throw std::runtime_error(
        "No iterations with observed events. Unable to produce output.");

  // prepare final containers (reserve capacities)
  std::vector<int> sum1_iterNum;
  sum1_iterNum.reserve(ns1r);
  std::vector<unsigned char> sum1_evNotArch;
  sum1_evNotArch.reserve(ns1r);
  std::vector<double> sum1_analysisT;
  sum1_analysisT.reserve(ns1r);
  std::vector<int> sum1_trtGrp;
  sum1_trtGrp.reserve(ns1r);
  std::vector<int> sum1_accruals;
  sum1_accruals.reserve(ns1r);
  std::vector<int> sum1_events;
  sum1_events.reserve(ns1r);
  std::vector<int> sum1_dropouts;
  sum1_dropouts.reserve(ns1r);

  std::vector<unsigned char> sum2_reject;
  sum2_reject.reserve(ns2r);
  std::vector<ListPtr> sumdata2;
  sumdata2.reserve(ns2r);

  // raw final containers
  std::vector<int> raw_iterNum;
  raw_iterNum.reserve(nrr);
  std::vector<double> raw_analysisT;
  raw_analysisT.reserve(nrr);
  std::vector<int> raw_subjectId;
  raw_subjectId.reserve(nrr);
  std::vector<double> raw_arrivalT;
  raw_arrivalT.reserve(nrr);
  std::vector<int> raw_stratum;
  raw_stratum.reserve(nrr);
  std::vector<int> raw_trtGrp;
  raw_trtGrp.reserve(nrr);
  std::vector<double> raw_survivalT;
  raw_survivalT.reserve(nrr);
  std::vector<double> raw_dropoutT;
  raw_dropoutT.reserve(nrr);
  std::vector<double> raw_timeObs;
  raw_timeObs.reserve(nrr);
  std::vector<unsigned char> raw_event;
  raw_event.reserve(nrr);
  std::vector<unsigned char> raw_dropEv;
  raw_dropEv.reserve(nrr);

  // flatten by iteration in order (preserves iteration order)
  for (size_t iter = 0; iter < maxIters; ++iter) {
    const auto& s1rows = results[iter].summary1Rows;
    for (const auto& r : s1rows) {
      sum1_iterNum.push_back(r.iterNum);
      sum1_evNotArch.push_back(r.evNotAch);
      sum1_analysisT.push_back(r.analysisT);
      sum1_trtGrp.push_back(r.trtGrp);
      sum1_accruals.push_back(r.accruals);
      sum1_events.push_back(r.events);
      sum1_dropouts.push_back(r.dropouts);
    }

    const auto& s2rows = results[iter].summary2Rows;

    if (!s2rows.empty()) {
      ListPtr sum2_fit = std::make_shared<ListCpp>();
      for (const auto& r : s2rows) {
        sum2_fit->push_back(r.iterNum, "iterationNumber");
        sum2_fit->push_back(r.analysisT, "analysisTime");
        sum2_fit->push_back(r.beta, "beta");
        sum2_fit->push_back(r.vbeta, "vbeta");
        sum2_fit->push_back(r.c, "optimalContrasts");
        sum2_fit->push_back(r.vzH0, "vzH0");
        sum2_fit->push_back(r.q, "criticalValueMaxZ");
        sum2_fit->push_back(r.z, "z");
        sum2_fit->push_back(r.optModel, "optimalModel");
        sum2_fit->push_back(r.reject, "reject");
        sum2_reject.push_back(r.reject);
      }
      sumdata2.push_back(std::move(sum2_fit));
    }

    if (iter < maxRawIters) {
      const auto& rraw = results[iter].rawRows;
      for (const auto& rr : rraw) {
        raw_iterNum.push_back(rr.iterNum);
        raw_analysisT.push_back(rr.analysisT);
        raw_subjectId.push_back(rr.subjectId);
        raw_arrivalT.push_back(rr.arrivalT);
        raw_stratum.push_back(rr.stratum);
        raw_trtGrp.push_back(rr.trtGrp);
        raw_survivalT.push_back(rr.survivalT);
        raw_dropoutT.push_back(rr.dropoutT);
        raw_timeObs.push_back(rr.timeObs);
        raw_event.push_back(rr.event);
        raw_dropEv.push_back(rr.dropEv);
      }
    }
  }

  const size_t rowsPerIter1 = M + 2;
  const size_t niters = ns2r;

  const double* sum1_T = sum1_analysisT.data();
  const int* sum1_E = sum1_events.data();
  const int* sum1_D = sum1_dropouts.data();
  const int* sum1_A = sum1_accruals.data();
  unsigned char* reject = sum2_reject.data();

  double expTime = 0.0, expEvents = 0.0, expDropouts = 0.0, expSubjects = 0.0;
  double overallReject = 0.0;
  for (size_t iter = 0; iter < niters; ++iter) {
    const size_t i1 = iter * rowsPerIter1 + M + 1;
    expTime += sum1_T[i1];
    expEvents += sum1_E[i1];
    expDropouts += sum1_D[i1];
    expSubjects += sum1_A[i1];
    overallReject += reject[iter];
  }
  expTime /= niters;
  expEvents /= niters;
  expDropouts /= niters;
  expSubjects /= niters;
  overallReject /= niters;

  ListCpp overview;

  FlatMatrix lambdaz(nstrata * nintv, M + 1);
  for (size_t m = 0; m < M + 1; ++m) {
    std::string name = std::string("lambdas[") + std::to_string(m) + "]";
    std::vector<double> v = lambdas[m];
    if (v.size() == 1) {
      std::fill(lambdaz.data_ptr() + m * nstrata * nintv,
                lambdaz.data_ptr() + (m + 1) * nstrata * nintv, v[0]);
    } else if (v.size() == nintv) {
      for (size_t s = 0; s < nstrata; ++s) {
        std::copy(v.begin(), v.end(),
                  lambdaz.data_ptr() + m * nstrata * nintv + s * nintv);
      }
    } else if (v.size() == nstrata * nintv) {
      std::copy(v.begin(), v.end(), lambdaz.data_ptr() + m * nstrata * nintv);
    } else {
      throw std::invalid_argument(std::string("Invalid length for ") + name);
    }
  }

  FlatMatrix gammaz(nstrata * nintv, M + 1);
  for (size_t m = 0; m < M + 1; ++m) {
    std::string name = std::string("gammas[") + std::to_string(m) + "]";
    std::vector<double> v = gammas[m];
    if (v.size() == 1) {
      std::fill(gammaz.data_ptr() + m * nstrata * nintv,
                gammaz.data_ptr() + (m + 1) * nstrata * nintv, v[0]);
    } else if (v.size() == nintv) {
      for (size_t s = 0; s < nstrata; ++s) {
        std::copy(v.begin(), v.end(),
                  gammaz.data_ptr() + m * nstrata * nintv + s * nintv);
      }
    } else if (v.size() == nstrata * nintv) {
      std::copy(v.begin(), v.end(), gammaz.data_ptr() + m * nstrata * nintv);
    } else {
      throw std::invalid_argument(std::string("Invalid length for ") + name);
    }
  }

  overview.push_back(overallReject, "overallReject");
  overview.push_back(expEvents, "expectedNumberOfEvents");
  overview.push_back(expDropouts, "expectedNumberOfDropouts");
  overview.push_back(expSubjects, "expectedNumberOfSubjects");
  overview.push_back(expTime, "expectedStudyDuration");
  overview.push_back(M, "M");
  overview.push_back(alpha, "alpha");
  overview.push_back(std::move(hrH0s), "hazardRatioH0s");
  overview.push_back(std::move(allocs), "allocations");
  overview.push_back(std::move(accrualTime), "accrualTime");
  overview.push_back(std::move(accrualIntensity), "accrualIntensity");
  overview.push_back(std::move(tau), "tau");
  overview.push_back(std::move(stratumFraction), "stratumFraction");
  overview.push_back(std::move(lambdaz), "lambdas");
  overview.push_back(candidateHazardRatios, "candidateHazardRatios");
  overview.push_back(std::move(gammaz), "gammas");
  overview.push_back(n, "n");
  overview.push_back(followupTime, "followupTime");
  overview.push_back(fixedFollowup, "fixedFollowup");
  overview.push_back(plannedEvents, "plannedEvents");
  overview.push_back(plannedTime, "plannedTime");
  overview.push_back(niters, "numberOfIterations");

  DataFrameCpp sumdata1;
  sumdata1.push_back(std::move(sum1_iterNum), "iterationNumber");
  sumdata1.push_back(std::move(sum1_evNotArch), "eventsNotAchieved");
  sumdata1.push_back(std::move(sum1_analysisT), "analysisTime");
  sumdata1.push_back(std::move(sum1_trtGrp), "treatmentGroup");
  sumdata1.push_back(std::move(sum1_accruals), "accruals");
  sumdata1.push_back(std::move(sum1_events), "events");
  sumdata1.push_back(std::move(sum1_dropouts), "dropouts");

  ListCpp result;
  result.push_back(std::move(overview), "overview");
  result.push_back(std::move(sumdata1), "sumdata1");
  result.push_back(std::move(sumdata2), "sumdata2");

  if (!raw_iterNum.empty()) {
    DataFrameCpp rawdata;
    rawdata.push_back(std::move(raw_iterNum), "iterationNumber");
    rawdata.push_back(std::move(raw_analysisT), "analysisTime");
    rawdata.push_back(std::move(raw_subjectId), "subjectId");
    rawdata.push_back(std::move(raw_arrivalT), "arrivalTime");
    rawdata.push_back(std::move(raw_stratum), "stratum");
    rawdata.push_back(std::move(raw_trtGrp), "treatmentGroup");
    rawdata.push_back(std::move(raw_survivalT), "survivalTime");
    rawdata.push_back(std::move(raw_dropoutT), "dropoutTime");
    rawdata.push_back(std::move(raw_timeObs), "timeUnderObservation");
    rawdata.push_back(std::move(raw_event), "event");
    rawdata.push_back(std::move(raw_dropEv), "dropoutEvent");

    result.push_back(std::move(rawdata), "rawdata");
  }

  return result;
}


// [[Rcpp::export]]
Rcpp::List lrsim_mcpmod_Rcpp(
    const int M = 2,
    const double alpha = 0.05,
    const Rcpp::NumericVector& hazardRatioH0s = 1,
    const Rcpp::NumericVector& allocations = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::Nullable<Rcpp::List> lambdas = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericMatrix> candidateHazardRatios= R_NilValue,
    const Rcpp::Nullable<Rcpp::List> gammas = R_NilValue,
    const int n = NA_INTEGER,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = false,
    const int plannedEvents = NA_INTEGER,
    const double plannedTime = NA_REAL,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = 0) {

  FlatMatrix candHRs;
  if (candidateHazardRatios.isNotNull()) {
    Rcpp::NumericMatrix cm(candidateHazardRatios); // unwrap
    if (cm.nrow() != M) {
      throw std::invalid_argument("Invalid rows of candidateHazardRatios");
    }
    candHRs = flatmatrix_from_Rmatrix(cm);
  } else {
    throw std::invalid_argument("candidateHazardRatios must be provided");
  }

  std::vector<double> hrH0s(hazardRatioH0s.begin(), hazardRatioH0s.end());
  std::vector<double> allocs(allocations.begin(), allocations.end());
  std::vector<double> accrualT(accrualTime.begin(), accrualTime.end());
  std::vector<double> accrualInt(accrualIntensity.begin(),
                                 accrualIntensity.end());
  std::vector<double> pwSurvT(piecewiseSurvivalTime.begin(),
                              piecewiseSurvivalTime.end());
  std::vector<double> stratumFrac(stratumFraction.begin(),
                                  stratumFraction.end());

  const int arms = M + 1;
  size_t nintv = pwSurvT.size();

  // lambdas: required (fail if not provided) and must have length M+1
  if (lambdas.isNull()) {
    throw std::invalid_argument("lambdas list is required");
  }
  Rcpp::List lambdasList(lambdas);
  if (static_cast<int>(lambdasList.size()) != arms) {
    throw std::invalid_argument("lambdas list must have length M + 1");
  }
  std::vector<std::vector<double>> lambdasVec(arms);
  for (int m = 0; m < arms; ++m) {
    Rcpp::NumericVector lamVec = lambdasList[m];
    lambdasVec[m] = std::vector<double>(lamVec.begin(), lamVec.end());
  }

  // gammas: if provided, must have length M+1;
  // if missing, create M+1 zero-vectors (length = nintv)
  std::vector<std::vector<double>> gammasVec(arms);
  if (gammas.isNull()) {
    // default: M+1 zero vectors (each length nintv)
    for (int m = 0; m < arms; ++m)
      gammasVec[m] = std::vector<double>(nintv, 0.0);
  } else {
    Rcpp::List gammasList(gammas);
    if (static_cast<int>(gammasList.size()) != arms) {
      throw std::invalid_argument("gammas list must have length M + 1");
    }
    for (int m = 0; m < arms; ++m) {
      Rcpp::NumericVector gamVec = gammasList[m];
      gammasVec[m] = std::vector<double>(gamVec.begin(), gamVec.end());
    }
  }
  auto out = lrsim_mcpmod_cpp(M, alpha, hrH0s, allocs, accrualT, accrualInt,
                              pwSurvT, stratumFrac, lambdasVec, candHRs,
                              gammasVec, n, followupTime, fixedFollowup,
                              plannedEvents, plannedTime, maxNumberOfIterations,
                              maxNumberOfRawDatasetsPerStage, seed);

  thread_utils::drain_thread_warnings_to_R();

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "lrsim_mcpmod";

  return result;
}
