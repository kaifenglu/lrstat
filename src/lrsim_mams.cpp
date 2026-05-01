#include "enrollment_event.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <unordered_set>
#include <stdexcept>
#include <string>
#include <vector>

#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/random.hpp>


using std::size_t;


// Parallel entry function
ListCpp lrsim_mams_cpp(
    const size_t M,
    const size_t kMax,
    const FlatMatrix& criticalValues,
    const std::vector<double>& hazardRatioH0s,
    const std::vector<double>& allocations,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<std::vector<double>>& lambdas,
    const std::vector<std::vector<double>>& gammas,
    const size_t n,
    const double followupTime,
    const bool fixedFollowup,
    const double rho1,
    const double rho2,
    const std::vector<int>& plannedEvents,
    const std::vector<double>& plannedTime,
    const int maxNumberOfIterations,
    const int maxNumberOfRawDatasetsPerStage,
    const int seed)
{
  if (M < 1) throw std::invalid_argument("M must be at least 1");
  if (kMax < 1) throw std::invalid_argument("kMax must be at least 1");

  // decide planning mode
  bool useEvents;
  if (none_na(plannedEvents)) {
    useEvents = true;
    if (plannedEvents[0] <= 0)
      throw std::invalid_argument("plannedEvents must be positive");
    if (plannedEvents.size() != kMax)
      throw std::invalid_argument("Invalid length for plannedEvents");
    if (any_nonincreasing(plannedEvents))
      throw std::invalid_argument("plannedEvents must be increasing");
  } else if (none_na(plannedTime)) {
    useEvents = false;
    if (plannedTime[0] <= 0.0)
      throw std::invalid_argument("plannedTime must be positive");
    if (plannedTime.size() != kMax)
      throw std::invalid_argument("Invalid length for plannedTime");
    if (any_nonincreasing(plannedTime))
      throw std::invalid_argument("plannedTime must be increasing");
  } else {
    throw std::invalid_argument("Either plannedEvents or plannedTime must be given");
  }

  // validate other input parameters
  std::vector<double> hrH0s = expand1(hazardRatioH0s, M, "hazardRatioH0s");
  if (std::any_of(hrH0s.begin(), hrH0s.end(), [](double v) { return v <= 0.0; }))
    throw std::invalid_argument("All hazardRatioH0 parameters must be positive");
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
    if (v < 0.0) throw std::invalid_argument("accrualIntensity must be non-negative");
  }
  if (piecewiseSurvivalTime[0] != 0.0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");
  for (double v : stratumFraction) {
    if (v <= 0.0) throw std::invalid_argument("stratumFraction must be positive");
  }
  double sumf = std::accumulate(stratumFraction.begin(), stratumFraction.end(), 0.0);
  if (std::fabs(sumf - 1.0) > 1e-12)
    throw std::invalid_argument("stratumFraction must sum to 1");

  for (size_t m = 0; m < M + 1; ++m) {
    std::string nm = std::string("lambdas[") + std::to_string(m) + "]";
    if (!none_na(lambdas[m])) throw std::invalid_argument(nm + "must be provided");
    for (double v : lambdas[m]) {
      if (v < 0.0) throw std::invalid_argument(nm + "must be non-negative");
    }
  }
  for (size_t m = 0; m < M + 1; ++m) {
    std::string nm = std::string("gammas[") + std::to_string(m) + "]";
    if (!none_na(gammas[m])) throw std::invalid_argument(nm + "must be provided");
    for (double v : gammas[m]) {
      if (v < 0.0) throw std::invalid_argument(nm + "must be non-negative");
    }
  }

  if (static_cast<int>(n) == INT_MIN)
    throw std::invalid_argument("n must be provided");
  if (n <= 0) throw std::invalid_argument("n must be a positive integer");
  if (fixedFollowup && std::isnan(followupTime))
    throw std::invalid_argument("followupTime must be provided for fixed follow-up");
  if (fixedFollowup && followupTime <= 0.0)
    throw std::invalid_argument("followupTime must be positive for fixed follow-up");
  if (rho1 < 0.0 || rho2 < 0.0)
    throw std::invalid_argument("rho parameters must be non-negative");
  if (maxNumberOfIterations < 1)
    throw std::invalid_argument("maxNumberOfIterations must be a positive integer");
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
    expand_stratified_to_slice(lambdas[m], lambdasx, m, nstrata, nintv, nm.c_str());
  }

  FlatArray gammasx(nintv, nstrata, M + 1);
  for (size_t m = 0; m < M + 1; ++m) {
    std::string nm = std::string("gammas[") + std::to_string(m) + "]";
    expand_stratified_to_slice(gammas[m], gammasx, m, nstrata, nintv, nm.c_str());
  }

  // generate seeds for each iteration to ensure reproducibility
  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter) seeds[iter] = master_rng();

  // One summary (stage-level) row produced by an iteration
  struct StageSummary1Row {
    int iterNum = 0;
    unsigned char evNotAch = 0;
    int stageNum = 0;
    double analysisT = 0.0;
    int trtGrp = 0;
    int accruals = 0, events = 0, dropouts = 0;
  };

  struct StageSummary2Row {
    int iterNum = 0;
    int stageNum = 0;
    double analysisT = 0.0;
    int actArm = 1;
    int totAccruals = 0, totEvents = 0, totDropouts = 0;
    double uscore = 0.0, vscore = 0.0, logRank = 0.0;
  };

  // One subject-level (raw) row for a particular iteration and stage
  struct RawDatasetRow {
    int iterNum = 0;
    int stageNum = 0;
    double analysisT = 0.0;
    int subjectId = 0;
    double arrivalT = 0.0;
    int stratum = 0, trtGrp = 0;
    double survivalT = 0.0, dropoutT = 0.0, timeObs = 0.0;
    unsigned char event = 0, dropEv = 0;
  };

  // Per-iteration container written exclusively by the worker thread
  struct IterationResult {
    std::vector<StageSummary1Row> summary1Rows;
    std::vector<StageSummary2Row> summary2Rows;
    std::vector<RawDatasetRow> rawRows;
    void reserveForSummary1(size_t approxRows) { summary1Rows.reserve(approxRows); }
    void reserveForSummary2(size_t approxRows) { summary2Rows.reserve(approxRows); }
    void reserveForRaw(size_t approxRows) { rawRows.reserve(approxRows); }
  };

  // pre-size per-iteration results
  std::vector<IterationResult> results;
  results.resize(maxIters);


  // Worker that runs simulation iterations [begin, end)
  struct SimWorker : public RcppParallel::Worker {
    // inputs (const refs)
    const size_t M;
    const size_t kMax;
    const std::vector<double>& hrH0s;
    const std::vector<double>& allocs;
    const std::vector<double>& accrualTime;
    const std::vector<double>& accrualIntensity;
    const std::vector<double>& tau;
    const std::vector<double>& stratumFraction;
    const FlatArray& lambdasx;
    const FlatArray& gammasx;
    const size_t n;
    const double fu;
    const bool fixedFollowup;
    const double rho1;
    const double rho2;
    const std::vector<int>& plannedEvents;
    const std::vector<double>& plannedTime;
    const size_t maxIters;
    const size_t maxRawIters; // store raw for iter < maxRawIters
    const std::vector<uint64_t>& seeds;
    const bool useEvents;
    const size_t nintv;
    const size_t nstrata;

    // output pointer (pre-sized vector of IterationResult)
    std::vector<IterationResult>* results;

    SimWorker(
      size_t M_,
      size_t kMax_,
      const std::vector<double>& hrH0s_,
      const std::vector<double>& allocs_,
      const std::vector<double>& accrualTime_,
      const std::vector<double>& accrualIntensity_,
      const std::vector<double>& tau_,
      const std::vector<double>& stratumFraction_,
      const FlatArray& lambdasx_,
      const FlatArray& gammasx_,
      size_t n_,
      double fu_,
      bool fixedFollowup_,
      double rho1_,
      double rho2_,
      const std::vector<int>& plannedEvents_,
      const std::vector<double>& plannedTime_,
      size_t maxIters_,
      size_t maxRawIters_,
      const std::vector<uint64_t>& seeds_,
      bool useEvents_,
      size_t nintv_,
      size_t nstrata_,
      std::vector<IterationResult>* results_)
      : M(M_),
        kMax(kMax_),
        hrH0s(hrH0s_),
        allocs(allocs_),
        accrualTime(accrualTime_),
        accrualIntensity(accrualIntensity_),
        tau(tau_),
        stratumFraction(stratumFraction_),
        lambdasx(lambdasx_),
        gammasx(gammasx_),
        n(n_),
        fu(fu_),
        fixedFollowup(fixedFollowup_),
        rho1(rho1_),
        rho2(rho2_),
        plannedEvents(plannedEvents_),
        plannedTime(plannedTime_),
        maxIters(maxIters_),
        maxRawIters(maxRawIters_),
        seeds(seeds_),
        useEvents(useEvents_),
        nintv(nintv_),
        nstrata(nstrata_),
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
      FlatMatrix kms(M, nstrata);
      std::vector<int> events(M2), dropouts(M2), accruals(M2);
      std::vector<double> us(M), vs(M), zs(M);
      std::vector<double> denom_per_stratum(nstrata);
      std::vector<double> cumF(nstrata);
      std::partial_sum(stratumFraction.begin(), stratumFraction.end(), cumF.begin());

      std::vector<double> analysisT(kMax);
      std::vector<double> totalte; totalte.reserve(n);
      std::vector<size_t> sub; sub.reserve(n);

      const double sumAlloc = std::accumulate(allocs.begin(), allocs.end(), 0.0);

      for (size_t iter = begin; iter < end; ++iter) {
        // deterministic per-iteration RNG
        boost::random::mt19937_64 rng_local(seeds[iter]);
        boost::random::uniform_real_distribution<double> unif(0.0, 1.0);

        // per-iteration output container
        IterationResult& out = (*results)[iter];
        out.summary1Rows.clear();
        out.summary2Rows.clear();
        out.rawRows.clear();
        if (iter < maxRawIters) out.reserveForRaw(kMax * n);
        out.reserveForSummary1(kMax * M1); // all arms
        out.reserveForSummary2(kMax * M); // all pairwise comparisons with control

        std::fill(denom_per_stratum.begin(), denom_per_stratum.end(), sumAlloc);

        // reset block randomization
        for (size_t h = 0; h < nstrata; ++h) {
          flatmatrix_set_column(bs, h, allocs);
        }

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
              k = m; break;
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
              timeObs[i] = sv; event[i] = 1; dropEv[i] = 0;
            } else if (dr <= sv && dr <= fu) {
              timeObs[i] = dr; event[i] = 0; dropEv[i] = 1;
            } else {
              timeObs[i] = fu; event[i] = 0; dropEv[i] = 0;
            }
          } else {
            if (sv <= dr) {
              timeObs[i] = sv; event[i] = 1; dropEv[i] = 0;
            } else {
              timeObs[i] = dr; event[i] = 0; dropEv[i] = 1;
            }
          }
          totalT[i] = arrivalT[i] + timeObs[i];
        } // cohort generated

        // determine analysis times (events counted only for arms 1 and control)
        size_t nstages = kMax;
        bool evNotAch = false;

        if (useEvents) {
          totalte.clear();
          int nevents = 0; // events involving arm1 or control in this iteration
          for (size_t i = 0; i < n; ++i) {
            if (event[i] && (trtGrp[i] == 1 || trtGrp[i] == static_cast<int>(M1))) {
              ++nevents; totalte.push_back(totalT[i]);
            }
          }
          if (nevents == 0) {
            thread_utils::push_thread_warning(
              std::string("No events for iteration ") + std::to_string(iter + 1) +
                " skipping this iteration.");
            // leave out.summaryRows empty to signal skipped iteration
            out.summary1Rows.clear();
            out.summary2Rows.clear();
            out.rawRows.clear();
            continue;
          }
          std::sort(totalte.begin(), totalte.end());

          size_t j;
          for (j = 0; j < kMax; ++j) {
            if (plannedEvents[j] >= nevents) { nstages = j + 1; break; }
          }

          if (j == kMax) {
            for (size_t k = 0; k < nstages; ++k) {
              analysisT[k] = totalte[plannedEvents[k] - 1] + 1e-12;
            }
          } else {
            for (size_t k = 0; k < nstages - 1; ++k) {
              analysisT[k] = totalte[plannedEvents[k] - 1] + 1e-12;
            }
            analysisT[nstages - 1] = totalte.back() + 1e-12;
          }
          evNotAch = (nevents < plannedEvents[kMax - 1]);
        } else {
          std::copy_n(plannedTime.begin(), kMax, analysisT.begin());
          evNotAch = false;
        }

        // per-stage calculations
        for (size_t k = 0; k < nstages; ++k) {
          double time = analysisT[k];

          // reset counts
          ns.fill(0);
          std::fill(events.begin(), events.end(), 0);
          std::fill(dropouts.begin(), dropouts.end(), 0);

          // censoring at analysis time and count accruals/events/dropouts
          for (size_t i = 0; i < n; ++i) {
            double ar = arrivalT[i], sv = survivalT[i], dr = dropoutT[i];

            if (ar > time) {
              timeObs[i] = time - ar; event[i] = 0; dropEv[i] = 0; continue;
            }

            if (fixedFollowup) {
              if (ar + sv <= time && sv <= dr && sv <= fu) {
                timeObs[i] = sv; event[i] = 1; dropEv[i] = 0;
              } else if (ar + dr <= time && dr <= sv && dr <= fu) {
                timeObs[i] = dr; event[i] = 0; dropEv[i] = 1;
              } else if (ar + fu <= time && fu <= sv && fu <= dr) {
                timeObs[i] = fu; event[i] = 0; dropEv[i] = 0;
              } else {
                timeObs[i] = time - ar; event[i] = 0; dropEv[i] = 0;
              }
            } else {
              if (ar + sv <= time && sv <= dr) {
                timeObs[i] = sv; event[i] = 1; dropEv[i] = 0;
              } else if (ar + dr <= time && dr <= sv) {
                timeObs[i] = dr; event[i] = 0; dropEv[i] = 1;
              } else {
                timeObs[i] = time - ar; event[i] = 0; dropEv[i] = 0;
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

          // collect indices with positive observed time and sort them
          sub.clear();
          for (size_t i = 0; i < n; ++i) if (timeObs[i] > 0.0) sub.push_back(i);
          std::sort(sub.begin(), sub.end(), [&](size_t i, size_t j) {
            return timeObs[i] < timeObs[j];
          });

          // compute stratified log-rank for pairwise comparisons vs. control
          kms.fill(1.0);
          std::fill(us.begin(), us.end(), 0.0);
          std::fill(vs.begin(), vs.end(), 0.0);

          for (size_t i = 0; i < sub.size(); ++i) {
            size_t idx = sub[i];
            size_t h = static_cast<size_t>(stratum[idx] - 1);
            size_t g = static_cast<size_t>(trtGrp[idx] - 1);

            double n2h = static_cast<double>(ns(M, h)); // control
            for (size_t m = 0; m < M; ++m) {
              double n1h = static_cast<double>(ns(m, h)); // treatment arm m

              double n1a = n1h * hrH0s[m];
              double nt = n1h + n2h;
              double nta = n1a + n2h;

              if (event[idx]) {
                if (g == m || g == M) {
                  double wh = 1.0;
                  if (rho1 != 0.0 || rho2 != 0.0) {
                    wh = std::pow(kms(m, h), rho1) * std::pow(1.0 - kms(m, h), rho2);
                    kms(m, h) *= (1.0 - 1.0 / nt);
                  }
                  double treated = (g == m ? 1.0 : 0.0);
                  us[m] += wh * (treated - n1a / nta);
                  vs[m] += wh * wh * n1a * n2h / (nta * nta);
                }
              }
            }

            // reduce risk set
            --ns(g, h);
          } // end events loop

          for (size_t m = 0; m < M; ++m) {
            zs[m] = (vs[m] > 0.0 ? us[m] / std::sqrt(vs[m]) : 0.0);
          }

          // optionally append raw rows for this stage
          if (iter < maxRawIters) { // only for first maxRawIters iterations
            for (size_t i = 0; i < n; ++i) {
              // skip subjects who haven't been enrolled by analysis time
              if (arrivalT[i] > time) continue;
              RawDatasetRow rr;
              rr.iterNum = static_cast<int>(iter + 1);
              rr.stageNum = static_cast<int>(k + 1);
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
          StageSummary1Row sr1;
          for (size_t m = 0; m < M2; ++m) {
            sr1.iterNum = static_cast<int>(iter + 1);
            sr1.evNotAch = evNotAch ? 1 : 0;
            sr1.stageNum = static_cast<int>(k + 1);
            sr1.analysisT = time;
            sr1.trtGrp = static_cast<int>(m + 1);
            sr1.accruals = accruals[m];
            sr1.events = events[m];
            sr1.dropouts = dropouts[m];
            out.summary1Rows.push_back(sr1);
          }

          StageSummary2Row sr2;
          for (size_t m = 0; m < M; ++m) {
            sr2.iterNum = static_cast<int>(iter + 1);
            sr2.stageNum = static_cast<int>(k + 1);
            sr2.analysisT = time;
            sr2.actArm = static_cast<int>(m + 1);
            sr2.totAccruals = accruals[m] + accruals[M];
            sr2.totEvents = events[m] + events[M];
            sr2.totDropouts = dropouts[m] + dropouts[M];
            sr2.uscore = us[m];
            sr2.vscore = vs[m];
            sr2.logRank = zs[m];
            out.summary2Rows.push_back(sr2);
          }
        } // per-stage
      } // per-iteration
    } // operator()
  }; // SimWorker

  // run worker in parallel
  SimWorker worker(
      M, kMax, hrH0s, allocs, accrualTime, accrualIntensity,
      tau, stratumFraction, lambdasx, gammasx, n, fu,
      fixedFollowup, rho1, rho2, plannedEvents, plannedTime,
      maxIters, maxRawIters, seeds, useEvents, nintv, nstrata,
      &results
  );

  RcppParallel::parallelFor(0, maxIters, worker);

  // Flatten results
  size_t ns1r = 0, ns2r = 0, nrr = 0;
  for (size_t iter = 0; iter < maxIters; ++iter) {
    ns1r += results[iter].summary1Rows.size();
    ns2r += results[iter].summary2Rows.size();
    nrr += results[iter].rawRows.size();
  }
  if (ns1r == 0) throw std::runtime_error(
    "No iterations with observed events for arm 1 or common control. "
    "Unable to produce output.");

  // prepare final containers (reserve capacities)
  std::vector<int> sum1_iterNum; sum1_iterNum.reserve(ns1r);
  std::vector<unsigned char> sum1_evNotArch; sum1_evNotArch.reserve(ns1r);
  std::vector<int> sum1_stopStage; sum1_stopStage.reserve(ns1r);
  std::vector<int> sum1_stageNum; sum1_stageNum.reserve(ns1r);
  std::vector<double> sum1_analysisT; sum1_analysisT.reserve(ns1r);
  std::vector<int> sum1_trtGrp; sum1_trtGrp.reserve(ns1r);
  std::vector<int> sum1_accruals; sum1_accruals.reserve(ns1r);
  std::vector<int> sum1_events; sum1_events.reserve(ns1r);
  std::vector<int> sum1_dropouts; sum1_dropouts.reserve(ns1r);

  std::vector<int> sum2_iterNum; sum2_iterNum.reserve(ns2r);
  std::vector<int> sum2_stopStage; sum2_stopStage.reserve(ns2r);
  std::vector<int> sum2_stageNum; sum2_stageNum.reserve(ns2r);
  std::vector<double> sum2_analysisT; sum2_analysisT.reserve(ns2r);
  std::vector<int> sum2_actArm; sum2_actArm.reserve(ns2r);
  std::vector<int> sum2_totAccruals; sum2_totAccruals.reserve(ns2r);
  std::vector<int> sum2_totEvents; sum2_totEvents.reserve(ns2r);
  std::vector<int> sum2_totDropouts; sum2_totDropouts.reserve(ns2r);
  std::vector<double> sum2_uscore; sum2_uscore.reserve(ns2r);
  std::vector<double> sum2_vscore; sum2_vscore.reserve(ns2r);
  std::vector<double> sum2_logRank; sum2_logRank.reserve(ns2r);
  std::vector<unsigned char> sum2_reject; sum2_reject.reserve(ns2r);

  // raw final containers
  std::vector<int> raw_iterNum; raw_iterNum.reserve(nrr);
  std::vector<int> raw_stopStage; raw_stopStage.reserve(nrr);
  std::vector<int> raw_stageNum; raw_stageNum.reserve(nrr);
  std::vector<double> raw_analysisT; raw_analysisT.reserve(nrr);
  std::vector<int> raw_subjectId; raw_subjectId.reserve(nrr);
  std::vector<double> raw_arrivalT; raw_arrivalT.reserve(nrr);
  std::vector<int> raw_stratum; raw_stratum.reserve(nrr);
  std::vector<int> raw_trtGrp; raw_trtGrp.reserve(nrr);
  std::vector<double> raw_survivalT; raw_survivalT.reserve(nrr);
  std::vector<double> raw_dropoutT; raw_dropoutT.reserve(nrr);
  std::vector<double> raw_timeObs; raw_timeObs.reserve(nrr);
  std::vector<unsigned char> raw_event; raw_event.reserve(nrr);
  std::vector<unsigned char> raw_dropEv; raw_dropEv.reserve(nrr);

  // flatten by iteration in order (preserves iteration order)
  for (size_t iter = 0; iter < maxIters; ++iter) {
    const auto& s1rows = results[iter].summary1Rows;
    for (const auto& r : s1rows) {
      sum1_iterNum.push_back(r.iterNum);
      sum1_evNotArch.push_back(r.evNotAch);
      sum1_stageNum.push_back(r.stageNum);
      sum1_analysisT.push_back(r.analysisT);
      sum1_trtGrp.push_back(r.trtGrp);
      sum1_accruals.push_back(r.accruals);
      sum1_events.push_back(r.events);
      sum1_dropouts.push_back(r.dropouts);
      sum1_stopStage.push_back(0);
    }

    const auto& s2rows = results[iter].summary2Rows;
    for (const auto& r : s2rows) {
      sum2_iterNum.push_back(r.iterNum);
      sum2_stageNum.push_back(r.stageNum);
      sum2_analysisT.push_back(r.analysisT);
      sum2_actArm.push_back(r.actArm);
      sum2_totAccruals.push_back(r.totAccruals);
      sum2_totEvents.push_back(r.totEvents);
      sum2_totDropouts.push_back(r.totDropouts);
      sum2_uscore.push_back(r.uscore);
      sum2_vscore.push_back(r.vscore);
      sum2_logRank.push_back(r.logRank);
      sum2_stopStage.push_back(0);
      sum2_reject.push_back(0);
    }

    if (iter < maxRawIters) {
      const auto& rraw = results[iter].rawRows;
      for (const auto& rr : rraw) {
        raw_iterNum.push_back(rr.iterNum);
        raw_stageNum.push_back(rr.stageNum);
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
        raw_stopStage.push_back(0);
      }
    }
  }

  const size_t rowsPerIter1 = kMax * (M + 2);
  const size_t rowsPerIter2 = kMax * M;
  const size_t niters = sum2_iterNum.size() / rowsPerIter2;

  std::vector<double> haveStage(kMax);
  FlatMatrix rejectByArm(kMax, M + 1);
  FlatMatrix timeByArm(kMax, M + 2);
  FlatMatrix eventsByArm(kMax, M + 2);
  FlatMatrix dropoutsByArm(kMax, M + 2);
  FlatMatrix subjectsByArm(kMax, M + 2);
  std::vector<double> expTimeByArm(M + 2);
  std::vector<double> expEventsByArm(M + 2);
  std::vector<double> expDropoutsByArm(M + 2);
  std::vector<double> expSubjectsByArm(M + 2);

  int* stopr = raw_stopStage.data();
  int* stop1 = sum1_stopStage.data();
  const double* sum1_T = sum1_analysisT.data();
  const int* sum1_E = sum1_events.data();
  const int* sum1_D = sum1_dropouts.data();
  const int* sum1_A = sum1_accruals.data();
  int* stop2 = sum2_stopStage.data();
  const double* logRank = sum2_logRank.data();
  unsigned char* reject = sum2_reject.data();

  size_t rawnum = 0;
  for (size_t iter = 0; iter < niters; ++iter) {
    const size_t i1 = iter * rowsPerIter1;
    const size_t i2 = iter * rowsPerIter2;

    // find stopping stage for this iteration
    size_t stop_k = kMax - 1;
    for (size_t k = 0; k < kMax; ++k) {
      const size_t offset = i2 + k * M;

      double cut = criticalValues(k, 0); // cutoff for level-M test

      std::unordered_set<size_t> I; // set of unrejected hypotheses
      for (size_t m = 0; m < M; ++m) I.insert(m);

      // check whether there is any arm that crosses the efficacy boundary
      bool anyreject = false;
      for (auto it = I.begin(); it != I.end(); ) {
        size_t m = *it;
        if (logRank[offset + m] < -cut) {
          anyreject = true;
          reject[offset + m] = 1;
          rejectByArm(k, m) += 1;
          it = I.erase(it); // erase returns next iterator (C++11+)
        } else {
          ++it;
        }
      }

      if (anyreject) {
        stop_k = k; // stopping stage
        rejectByArm(k, M) += 1;

        // is there any other arm that crosses the relaxed boundary?
        while (anyreject) {
          anyreject = false;
          if (!I.empty()) {
            double cut = criticalValues(k, M - I.size());
            for (auto it = I.begin(); it != I.end(); ) {
              size_t m = *it;
              if (logRank[offset + m] < -cut) {
                anyreject = true;
                reject[offset + m] = 1;
                rejectByArm(k, m) += 1;
                it = I.erase(it);
              } else {
                ++it;
              }
            }
          }
        }

        break;
      }
    }

    // assign stop stage for each iteration for the summary data sets
    for (size_t k = 0; k < kMax; ++k) {
      const size_t offset1 = i1 + k * (M + 2);
      for (size_t m = 0; m < M + 2; ++m) {
        stop1[offset1 + m] = stop_k + 1;
      }

      const size_t offset2 = i2 + k * M;
      for (size_t m = 0; m < M; ++m) {
        stop2[offset2 + m] = stop_k + 1;
      }
    }

    if (iter < maxRawIters) {
      size_t n1 = 0;
      for (size_t k = 0; k < kMax; ++k) {
        n1 += sum1_A[i1 + k * (M + 2) + (M + 1)];
      }

      for (size_t i = 0; i < n1; ++i) {
        stopr[rawnum + i] = stop_k + 1;
      }
      rawnum += n1;
    }

    // tally time and number of events/dropouts/subjects
    for (size_t k = 0; k <= stop_k; ++k) {
      haveStage[k] += 1;
      const size_t offset = i1 + k * (M + 2);

      for (size_t m = 0; m < M + 2; ++m) {
        size_t idx = offset + m;
        timeByArm(k, m) += sum1_T[idx];
        eventsByArm(k, m) += sum1_E[idx];
        dropoutsByArm(k, m) += sum1_D[idx];
        subjectsByArm(k, m) += sum1_A[idx];

        if (k == stop_k) {
          expTimeByArm[m] += sum1_T[idx];
          expEventsByArm[m] += sum1_E[idx];
          expDropoutsByArm[m] += sum1_D[idx];
          expSubjectsByArm[m] += sum1_A[idx];
        }
      }
    }
  }


  // empirical cumulative rejection rates by stage and treatment
  for (size_t m = 0; m < M + 1; ++m) {
    for (size_t k = 0; k < kMax; ++k) {
      rejectByArm(k, m) /= niters;
    }
  }

  FlatMatrix cumRejectByArm(kMax, M + 1);
  for (size_t m = 0; m < M + 1; ++m) {
    cumRejectByArm(0, m) = rejectByArm(0, m);
    for (size_t k = 1; k < kMax; ++k)
      cumRejectByArm(k, m) = cumRejectByArm(k - 1, m) + rejectByArm(k, m);
  }

  double overallReject = cumRejectByArm(kMax - 1, M);

  // convert counts to proportions / averages
  for (size_t m = 0; m < M + 2; ++m) {
    expTimeByArm[m] /= niters;
    expEventsByArm[m] /= niters;
    expDropoutsByArm[m] /= niters;
    expSubjectsByArm[m] /= niters;
  }

  for (size_t m = 0; m < M + 2; ++m) {
    for (size_t k = 0; k < kMax; ++k) {
      const double denom = haveStage[k];
      timeByArm(k, m) /= denom;
      eventsByArm(k, m) /= denom;
      dropoutsByArm(k, m) /= denom;
      subjectsByArm(k, m) /= denom;
    }
  }


  ListCpp overview;
  overview.push_back(overallReject, "overallReject");
  overview.push_back(std::move(rejectByArm), "rejectPerStage");
  overview.push_back(std::move(cumRejectByArm), "cumulativeRejection");
  overview.push_back(std::move(eventsByArm), "numberOfEvents");
  overview.push_back(std::move(dropoutsByArm), "numberOfDropouts");
  overview.push_back(std::move(subjectsByArm), "numberOfSubjects");
  overview.push_back(std::move(timeByArm), "analysisTime");
  overview.push_back(std::move(expEventsByArm), "expectedNumberOfEvents");
  overview.push_back(std::move(expDropoutsByArm), "expectedNumberOfDropouts");
  overview.push_back(std::move(expSubjectsByArm), "expectedNumberOfSubjects");
  overview.push_back(std::move(expTimeByArm), "expectedStudyDuration");
  overview.push_back(criticalValues, "criticalValues");
  overview.push_back(hazardRatioH0s, "hazardRatioH0s");
  overview.push_back(useEvents, "useEvents");
  overview.push_back(niters, "numberOfIterations");
  overview.push_back(n, "n");
  overview.push_back(fixedFollowup, "fixedFollowup");
  overview.push_back(rho1, "rho1");
  overview.push_back(rho2, "rho2");
  overview.push_back(M, "M");
  overview.push_back(kMax, "kMax");

  DataFrameCpp sumdata1;
  sumdata1.push_back(std::move(sum1_iterNum), "iterationNumber");
  sumdata1.push_back(std::move(sum1_evNotArch), "eventsNotAchieved");
  sumdata1.push_back(std::move(sum1_stopStage), "stopStage");
  sumdata1.push_back(std::move(sum1_stageNum), "stageNumber");
  sumdata1.push_back(std::move(sum1_analysisT), "analysisTime");
  sumdata1.push_back(std::move(sum1_trtGrp), "treatmentGroup");
  sumdata1.push_back(std::move(sum1_accruals), "accruals");
  sumdata1.push_back(std::move(sum1_events), "events");
  sumdata1.push_back(std::move(sum1_dropouts), "dropouts");

  DataFrameCpp sumdata2;
  sumdata2.push_back(std::move(sum2_iterNum), "iterationNumber");
  sumdata2.push_back(std::move(sum2_stopStage), "stopStage");
  sumdata2.push_back(std::move(sum2_stageNum), "stageNumber");
  sumdata2.push_back(std::move(sum2_analysisT), "analysisTime");
  sumdata2.push_back(std::move(sum2_actArm), "activeArm");
  sumdata2.push_back(std::move(sum2_totAccruals), "totalAccruals");
  sumdata2.push_back(std::move(sum2_totEvents), "totalEvents");
  sumdata2.push_back(std::move(sum2_totDropouts), "totalDropouts");
  sumdata2.push_back(std::move(sum2_uscore), "uscore");
  sumdata2.push_back(std::move(sum2_vscore), "vscore");
  sumdata2.push_back(std::move(sum2_logRank), "logRankStatistic");
  sumdata2.push_back(std::move(sum2_reject), "reject");

  ListCpp result;
  result.push_back(std::move(overview), "overview");
  result.push_back(std::move(sumdata1), "sumdata1");
  result.push_back(std::move(sumdata2), "sumdata2");

  // attach raw data
  if (!raw_iterNum.empty()) {
    DataFrameCpp rawdata;
    rawdata.push_back(std::move(raw_iterNum), "iterationNumber");
    rawdata.push_back(std::move(raw_stopStage), "stopStage");
    rawdata.push_back(std::move(raw_stageNum), "stageNumber");
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
Rcpp::List lrsim_mams_Rcpp(
    const int M = 2,
    const int kMax = 1,
    const Rcpp::Nullable<Rcpp::NumericMatrix> criticalValues = R_NilValue,
    const Rcpp::NumericVector& hazardRatioH0s = 1,
    const Rcpp::NumericVector& allocations = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::Nullable<Rcpp::List> lambdas = R_NilValue,
    const Rcpp::Nullable<Rcpp::List> gammas = R_NilValue,
    const int n = NA_INTEGER,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = false,
    const double rho1 = 0,
    const double rho2 = 0,
    const Rcpp::IntegerVector& plannedEvents = NA_INTEGER,
    const Rcpp::NumericVector& plannedTime = NA_REAL,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = 0) {

  FlatMatrix critValues;
  if (criticalValues.isNotNull()) {
    Rcpp::NumericMatrix cm(criticalValues); // unwrap
    if (cm.nrow() != kMax || cm.ncol() != M) {
      throw std::invalid_argument("Invalid dimensions of criticalValues");
    }
    critValues = flatmatrix_from_Rmatrix(cm);
  } else {
    throw std::invalid_argument("criticalValues must be provided");
  }

  std::vector<double> hrH0s(hazardRatioH0s.begin(), hazardRatioH0s.end());
  std::vector<double> allocs(allocations.begin(), allocations.end());
  std::vector<double> accrualT(accrualTime.begin(), accrualTime.end());
  std::vector<double> accrualInt(accrualIntensity.begin(), accrualIntensity.end());
  std::vector<double> pwSurvT(piecewiseSurvivalTime.begin(),
                              piecewiseSurvivalTime.end());
  std::vector<double> stratumFrac(stratumFraction.begin(), stratumFraction.end());

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
    for (int m = 0; m < arms; ++m) gammasVec[m] = std::vector<double>(nintv, 0.0);
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

  std::vector<int> plannedE(plannedEvents.begin(), plannedEvents.end());
  std::vector<double> plannedT(plannedTime.begin(), plannedTime.end());

  auto out = lrsim_mams_cpp(
    M, kMax, critValues, hrH0s, allocs, accrualT, accrualInt,
    pwSurvT, stratumFrac, lambdasVec, gammasVec,
    n, followupTime, fixedFollowup, rho1, rho2, plannedE, plannedT,
    maxNumberOfIterations, maxNumberOfRawDatasetsPerStage, seed);

  thread_utils::drain_thread_warnings_to_R();

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "lrsim_mams";

  return result;
}
