// [[Rcpp::depends(RcppParallel)]]

#include "enrollment_event.h"
#include "generic_design.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/random.hpp>


using std::size_t;


// The parallel entry function
ListCpp lrsimcpp(
    const int kMax,
    const std::vector<double>& informationRates,
    const std::vector<double>& criticalValues,
    const std::vector<double>& futilityBounds,
    const double hazardRatioH0,
    const int allocation1,
    const int allocation2,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const int n,
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
  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
  size_t K = static_cast<size_t>(kMax);

  // whether to plan the analyses based on events or calendar time
  bool useEvents;
  if (none_na(plannedEvents)) {
    useEvents = true;
    if (plannedEvents[0] <= 0)
      throw std::invalid_argument("pannedEvents must be positive");
    if (plannedEvents.size() != K)
      throw std::invalid_argument("Invalid length for plannedEvents");
    if (any_nonincreasing(plannedEvents))
      throw std::invalid_argument("plannedEvents must be increasing");
  } else if (none_na(plannedTime)) {
    useEvents = false;
    if (plannedTime[0] <= 0)
      throw std::invalid_argument("plannedTime must be positive");
    if (plannedTime.size() != K)
      throw std::invalid_argument("Invalid length for plannedTime");
    if (any_nonincreasing(plannedTime))
      throw std::invalid_argument("plannedTime must be increasing");
  } else {
    throw std::invalid_argument("Either plannedEvents or plannedTime must be given");
  }

  // validate informationRates and set defaults
  std::vector<double> infoRates(K);
  if (none_na(informationRates)) {
    if (informationRates.size() != K)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[K-1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else if (useEvents) {
    double totalPlannedEvents = static_cast<double>(plannedEvents[K - 1]);
    for (size_t i = 0; i < K; ++i)
      infoRates[i] = static_cast<double>(plannedEvents[i]) / totalPlannedEvents;
  } else {
    double totalPlannedTime = plannedTime[K - 1];
    for (size_t i = 0; i < K; ++i)
      infoRates[i] = plannedTime[i] / totalPlannedTime;
  }

  // validate criticalValues and futilityBounds
  if (!none_na(criticalValues))
    throw std::invalid_argument("criticalValues must be provided");
  if (criticalValues.size() != K)
    throw std::invalid_argument("Invalid length for criticalValues");

  std::vector<double> futBounds = futilityBounds;
  if (K > 1 && !none_na(futilityBounds))
    futBounds = std::vector<double>(K - 1, -6.0);
  if (none_na(futBounds) && futBounds.size() < K - 1)
    throw std::invalid_argument("Invalid length for futilityBounds");
  if (none_na(criticalValues) && none_na(futBounds)) {
    for (size_t i = 0; i < K - 1; ++i) {
      if (futBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
  }

  // validate other parameters
  if (hazardRatioH0 <= 0.0)
    throw std::invalid_argument("hazardRatioH0 must be positive");
  if (allocation1 < 1 || allocation2 < 1)
    throw std::invalid_argument("allocations must be positive integers");
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
  if (!none_na(lambda1)) throw std::invalid_argument("lambda1 must be provided");
  if (!none_na(lambda2)) throw std::invalid_argument("lambda2 must be provided");
  for (double v : lambda1) {
    if (v < 0.0) throw std::invalid_argument("lambda1 must be non-negative");
  }
  for (double v : lambda2) {
    if (v < 0.0) throw std::invalid_argument("lambda2 must be non-negative");
  }
  for (double v : gamma1) {
    if (v < 0.0) throw std::invalid_argument("gamma1 must be non-negative");
  }
  for (double v : gamma2) {
    if (v < 0.0) throw std::invalid_argument("gamma2 must be non-negative");
  }
  if (n == INT_MIN) throw std::invalid_argument("n must be provided");
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

  size_t N = static_cast<size_t>(n);
  size_t maxIters = static_cast<size_t>(maxNumberOfIterations);
  size_t maxRawIters = static_cast<size_t>(maxNumberOfRawDatasetsPerStage);
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  const std::vector<double>& tau = piecewiseSurvivalTime;
  const double fu = followupTime;

  // expand stratified inputs
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");

  // calculate total alpha once
  std::vector<double> lb(K, -6.0), zero(K, 0.0);
  ListCpp exitprobs = exitprobcpp(criticalValues, lb, zero, infoRates);
  auto exitUpper = exitprobs.get<std::vector<double>>("exitProbUpper");
  double alpha = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);

  // generate seeds for each iteration to ensure reproducibility
  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter) seeds[iter] = master_rng();


  // One summary (stage-level) row produced by an iteration
  struct StageSummaryRow {
    int iterNum = 0;
    unsigned char evNotAch = 0;
    int stopStage = 0, stageNum = 0;
    double analysisT = 0.0;
    int accruals1 = 0, accruals2 = 0, totAccruals = 0;
    int events1 = 0, events2 = 0, totEvents = 0;
    int dropouts1 = 0, dropouts2 = 0, totDropouts = 0;
    double uscore = 0.0, vscore = 0.0, logRank = 0.0;
    unsigned char rejPerStage = 0, futPerStage = 0;
  };


  // One subject-level (raw) row for a particular iteration and stage
  struct RawDatasetRow {
    int iterNum = 0, stopStage = 0, stageNum = 0;
    double analysisT = 0.0;
    int subjectId = 0;
    double arrivalT = 0.0;
    int stratum = 0, trtGrp = 0;
    double survivalT = 0.0, dropoutT = 0.0, timeObs = 0.0;
    unsigned char event = 0, dropEv = 0;
  };

  // Per-iteration container written exclusively by the worker thread
  // responsible for that iteration
  struct IterationResult {
    std::vector<StageSummaryRow> summaryRows;
    std::vector<RawDatasetRow> rawRows; // populated only for first M iterations
    void reserveForSummary(size_t approxRows) { summaryRows.reserve(approxRows); }
    void reserveForRaw(size_t approxRows) { rawRows.reserve(approxRows); }
  };

  // pre-size per-iteration results
  std::vector<IterationResult> results;
  results.resize(maxIters);


  // Parallel worker
  struct SimWorker : public RcppParallel::Worker {
    // inputs (const refs)
    const size_t K;
    const std::vector<double>& infoRates;
    const std::vector<double>& criticalValues;
    const std::vector<double>& futBounds;
    const double hazardRatioH0;
    const int allocation1;
    const int allocation2;
    const std::vector<double>& accrualTime;
    const std::vector<double>& accrualIntensity;
    const std::vector<double>& tau;
    const std::vector<double>& stratumFraction;
    const std::vector<std::vector<double>>& lambda1x;
    const std::vector<std::vector<double>>& lambda2x;
    const std::vector<std::vector<double>>& gamma1x;
    const std::vector<std::vector<double>>& gamma2x;
    const size_t N;
    const double fu;
    const bool fixedFollowup;
    const double rho1;
    const double rho2;
    const std::vector<int>& plannedEvents;
    const std::vector<double>& plannedTime;
    const size_t maxIters;
    const size_t maxRawIters; // save raw rows only for iter < maxRawIters
    const std::vector<uint64_t>& seeds;
    const bool useEvents;
    const size_t nstrata;
    const double alpha;

    // output pointer (pre-sized vector of IterationResult with length maxIters)
    std::vector<IterationResult>* results;

    SimWorker(
      size_t K_,
      const std::vector<double>& infoRates_,
      const std::vector<double>& criticalValues_,
      const std::vector<double>& futBounds_,
      double hazardRatioH0_,
      int allocation1_,
      int allocation2_,
      const std::vector<double>& accrualTime_,
      const std::vector<double>& accrualIntensity_,
      const std::vector<double>& tau_,
      const std::vector<double>& stratumFraction_,
      const std::vector<std::vector<double>>& lambda1x_,
      const std::vector<std::vector<double>>& lambda2x_,
      const std::vector<std::vector<double>>& gamma1x_,
      const std::vector<std::vector<double>>& gamma2x_,
      size_t N_,
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
      size_t nstrata_,
      double alpha_,
      std::vector<IterationResult>* results_)
      : K(K_),
        infoRates(infoRates_),
        criticalValues(criticalValues_),
        futBounds(futBounds_),
        hazardRatioH0(hazardRatioH0_),
        allocation1(allocation1_),
        allocation2(allocation2_),
        accrualTime(accrualTime_),
        accrualIntensity(accrualIntensity_),
        tau(tau_),
        stratumFraction(stratumFraction_),
        lambda1x(lambda1x_),
        lambda2x(lambda2x_),
        gamma1x(gamma1x_),
        gamma2x(gamma2x_),
        N(N_),
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
        nstrata(nstrata_),
        alpha(alpha_),
        results(results_)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      // local buffers reused
      std::vector<int> stratum(N), trtGrp(N);
      std::vector<double> arrivalT(N), survivalT(N), dropoutT(N);
      std::vector<double> timeObs(N), totalT(N);
      std::vector<unsigned char> event(N), dropEv(N);
      std::vector<int> b1(nstrata), b2(nstrata);
      std::vector<int> n1(nstrata), n2(nstrata);
      std::vector<double> km(nstrata);
      std::vector<double> cumF(nstrata);
      std::partial_sum(stratumFraction.begin(), stratumFraction.end(), cumF.begin());
      std::vector<int> obsEvents(K);
      std::vector<double> analysisT(K), vscore(K);
      std::vector<double> lb(K, -6.0), zero(K, 0.0), I(K);
      std::vector<double> totalte; totalte.reserve(N);
      std::vector<size_t> sub; sub.reserve(N);
      std::vector<double> critValues; critValues.reserve(K);
      std::vector<double> ub; ub.reserve(K);

      for (size_t iter = begin; iter < end; ++iter) {
        // deterministic per-iteration RNG
        boost::random::mt19937_64 rng_local(seeds[iter]);
        boost::random::uniform_real_distribution<double> unif(0.0, 1.0);

        // reset per-iteration results
        IterationResult& out = (*results)[iter];
        out.summaryRows.clear();
        out.rawRows.clear();
        if (iter < maxRawIters) out.reserveForRaw(K * N);
        out.reserveForSummary(K);

        // reset block randomization
        std::fill(b1.begin(), b1.end(), allocation1);
        std::fill(b2.begin(), b2.end(), allocation2);

        double enrollt = 0.0;

        // generate cohort (subject-level)
        for (size_t i = 0; i < N; ++i) {
          double u = unif(rng_local);
          enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt);
          arrivalT[i] = enrollt;

          u = unif(rng_local);
          int j = findInterval1(u, cumF);
          stratum[i] = j + 1;

          u = unif(rng_local);
          double denom = static_cast<double>(b1[j] + b2[j]);
          double p = static_cast<double>(b1[j]) / denom;
          if (u <= p) { trtGrp[i] = 1; --b1[j]; }
          else { trtGrp[i] = 2; --b2[j]; }
          if (b1[j] + b2[j] == 0) { b1[j] = allocation1; b2[j] = allocation2; }

          u = unif(rng_local);
          if (trtGrp[i] == 1) survivalT[i] = qtpwexpcpp1(u, tau, lambda1x[j]);
          else survivalT[i] = qtpwexpcpp1(u, tau, lambda2x[j]);

          u = unif(rng_local);
          if (trtGrp[i] == 1) dropoutT[i] = qtpwexpcpp1(u, tau, gamma1x[j]);
          else dropoutT[i] = qtpwexpcpp1(u, tau, gamma1x[j]);

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

        // determine analysis times and stage count
        size_t nstages = K;
        bool evNotAch = false;

        if (useEvents) {
          totalte.clear();
          int nevents = 0;
          for (size_t i = 0; i < N; ++i) {
            if (event[i]) { ++nevents; totalte.push_back(totalT[i]); }
          }
          if (nevents == 0) {
            thread_utils::push_thread_warning(
              std::string("No events for iteration ") + std::to_string(iter + 1) +
                " skipping this iteration.");
            // keep out.summaryRows empty to signal skipped iteration
            out.summaryRows.clear();
            out.rawRows.clear();
            continue;
          }
          std::sort(totalte.begin(), totalte.end());
          size_t j;
          for (j = 0; j < K; ++j) {
            if (plannedEvents[j] >= nevents) { nstages = j + 1; break; }
          }
          if (j == K) {
            // observed >= planned: analyses at planned events
            for (size_t k = 0; k < nstages; ++k) {
              analysisT[k] = totalte[plannedEvents[k] - 1] + 1e-12;
              obsEvents[k] = plannedEvents[k];
            }
          } else {
            // last analysis uses all observed events
            for (size_t k = 0; k < nstages - 1; ++k) {
              analysisT[k] = totalte[plannedEvents[k] - 1] + 1e-12;
              obsEvents[k] = plannedEvents[k];
            }
            analysisT[nstages - 1] = totalte.back() + 1e-12;
            obsEvents[nstages - 1] = nevents;
          }
          evNotAch = (nevents < plannedEvents[K - 1]);
        } else {
          // calendar-time looks
          std::copy_n(plannedTime.begin(), K, analysisT.begin());
          evNotAch = false;
        }

        // per-stage calculations; optionally collect raw rows if iter < maxRawIters
        int nstops = 0; // number of stopping stages for this iteration
        size_t stopStage = nstages;
        for (size_t k = 0; k < nstages; ++k) {
          double time = analysisT[k];

          // reset per-stratum counts
          std::fill(n1.begin(), n1.end(), 0);
          std::fill(n2.begin(), n2.end(), 0);
          int events1 = 0, events2 = 0, dropouts1 = 0, dropouts2 = 0;

          // censoring & counts
          for (size_t i = 0; i < N; ++i) {
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
            if (trtGrp[i] == 1) {
              ++n1[h];
              if (event[i]) ++events1; else if (dropEv[i]) ++dropouts1;
            } else {
              ++n2[h];
              if (event[i]) ++events2; else if (dropEv[i]) ++dropouts2;
            }
          } // end censoring

          int accruals1 = std::accumulate(n1.begin(), n1.end(), 0);
          int accruals2 = std::accumulate(n2.begin(), n2.end(), 0);
          int totAccruals = accruals1 + accruals2;
          int totEvents = events1 + events2;
          int totDropouts = dropouts1 + dropouts2;

          // collect indices with positive observed time and sort them
          sub.clear();
          for (size_t i = 0; i < N; ++i) if (timeObs[i] > 0.0) sub.push_back(i);
          std::sort(sub.begin(), sub.end(), [&](size_t a, size_t b) {
            return timeObs[a] < timeObs[b];
          });

          // compute stratified log-rank (single TTE endpoint)
          std::fill(km.begin(), km.end(), 1.0);
          double us = 0.0, vs = 0.0;

          for (size_t i = 0; i < sub.size(); ++i) {
            size_t idx = sub[i];
            size_t h = static_cast<size_t>(stratum[idx] - 1);

            double n1h = static_cast<double>(n1[h]);
            double n2h = static_cast<double>(n2[h]);
            double n1a = n1h * hazardRatioH0;
            double nt = n1h + n2h;
            double nta = n1a + n2h;

            if (event[idx]) { // at most one event can occur at any given time
              double wh = 1.0;
              if (rho1 != 0.0 || rho2 != 0.0) {
                wh = std::pow(km[h], rho1) * std::pow(1.0 - km[h], rho2);
                km[h] *= (1.0 - 1.0 / nt);
              }
              double treated = (trtGrp[idx] == 1 ? 1.0 : 0.0);
              us += wh * (treated - n1a / nta);
              vs += wh * wh * n1a * n2h / (nta * nta);
            }

            if (trtGrp[idx] == 1) --n1[h]; else --n2[h];
          }

          double z = (vs > 0.0) ? (us / std::sqrt(vs)) : 0.0;
          vscore[k] = vs;

          // adjust critical value at final stage if needed
          critValues = criticalValues; // copy
          if (useEvents) {
            if (k == nstages - 1 && evNotAch) {
              // no change to critical values at earlier stages, or at the
              // final stage if the planned total number of events is achieved
              // otherwise assign all remaining alpha to the final stage
              ub.resize(nstages);
              if (nstages > 1)
                std::copy_n(critValues.begin(), nstages - 1, ub.begin());
              if (rho1 == 0.0 && rho2 == 0.0) { // use events for std log-rank
                std::copy_n(obsEvents.begin(), nstages, I.begin());
                auto f = [&](double aval)->double {
                  ub[nstages - 1] = aval;
                  ListCpp probs = exitprobcpp(ub, lb, zero, I);
                  auto v = probs.get<std::vector<double>>("exitProbUpper");
                  return std::accumulate(v.begin(), v.end(), 0.0) - alpha;
                };
                critValues[nstages - 1] = brent(f, 0.0, 6.0, 1e-6);
              } else { // use vscore as information for weighted log-rank
                std::copy_n(vscore.begin(), nstages, I.begin());
                auto f = [&](double aval)->double {
                  ub[nstages - 1] = aval;
                  ListCpp probs = exitprobcpp(ub, lb, zero, I);
                  auto v = probs.get<std::vector<double>>("exitProbUpper");
                  return std::accumulate(v.begin(), v.end(), 0.0) - alpha;
                };
                critValues[nstages - 1] = brent(f, 0.0, 6.0, 1e-6);
              }
            }
          }

          // make decisions using -z because we are testing for a hazard ratio < 1
          double reject = 0, futility = 0;
          if (-z > critValues[k]) reject = 1;
          else if ((k < nstages - 1 && -z < futBounds[k]) || (k == nstages - 1))
            futility = 1;

          if (reject || futility) {
            ++nstops;
            if (nstops == 1) { // the first stage at which a stop occurs
              stopStage = k + 1;
            }
          }


          // optionally append raw rows for this stage
          if (iter < maxRawIters) { // only for first maxRawIters iterations
            for (size_t i = 0; i < N; ++i) {
              // skip subjects who haven't been enrolled by analysis time
              if (arrivalT[i] > time) continue;
              RawDatasetRow rr;
              rr.iterNum = static_cast<int>(iter + 1);
              rr.stopStage = static_cast<int>(stopStage);
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


          // append summary row
          StageSummaryRow sr;
          sr.iterNum = static_cast<int>(iter + 1);
          sr.evNotAch = evNotAch ? 1 : 0;
          sr.stopStage = static_cast<int>(stopStage);
          sr.stageNum = static_cast<int>(k + 1);
          sr.analysisT = time;
          sr.accruals1 = accruals1;
          sr.accruals2 = accruals2;
          sr.totAccruals = totAccruals;
          sr.events1 = events1;
          sr.events2 = events2;
          sr.totEvents = totEvents;
          sr.dropouts1 = dropouts1;
          sr.dropouts2 = dropouts2;
          sr.totDropouts = totDropouts;
          sr.uscore = us;
          sr.vscore = vs;
          sr.logRank = z;
          sr.rejPerStage = reject;
          sr.futPerStage = futility;
          out.summaryRows.push_back(std::move(sr));
        } // per-stage
      } // per-iteration
    } // operator()
  }; // SimWorker


  // create and run worker
  SimWorker worker(
      K, informationRates, criticalValues, futBounds,
      hazardRatioH0, allocation1, allocation2,
      accrualTime, accrualIntensity, tau, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      N, fu, fixedFollowup, rho1, rho2,
      plannedEvents, plannedTime,
      maxIters, maxRawIters, seeds, useEvents, nstrata, alpha,
      &results
  );

  RcppParallel::parallelFor(0, maxIters, worker);

  // Flatten per-iteration summary rows and raw rows
  size_t nsr = 0, nrr = 0;
  for (size_t iter = 0; iter < maxIters; ++iter) {
    nsr += results[iter].summaryRows.size();
    nrr += results[iter].rawRows.size();
  }
  if (nsr == 0) throw std::runtime_error(
    "No iterations with observed events. Unable to produce output.");

  // allocate final containers
  std::vector<int> sum_iterNum; sum_iterNum.reserve(nsr);
  std::vector<unsigned char> sum_evNotArch; sum_evNotArch.reserve(nsr);
  std::vector<int> sum_stopStage; sum_stopStage.reserve(nsr);
  std::vector<int> sum_stageNum; sum_stageNum.reserve(nsr);
  std::vector<double> sum_analysisT; sum_analysisT.reserve(nsr);
  std::vector<int> sum_accruals1; sum_accruals1.reserve(nsr);
  std::vector<int> sum_accruals2; sum_accruals2.reserve(nsr);
  std::vector<int> sum_totAccruals; sum_totAccruals.reserve(nsr);
  std::vector<int> sum_events1; sum_events1.reserve(nsr);
  std::vector<int> sum_events2; sum_events2.reserve(nsr);
  std::vector<int> sum_totEvents; sum_totEvents.reserve(nsr);
  std::vector<int> sum_dropouts1; sum_dropouts1.reserve(nsr);
  std::vector<int> sum_dropouts2; sum_dropouts2.reserve(nsr);
  std::vector<int> sum_totDropouts; sum_totDropouts.reserve(nsr);
  std::vector<double> sum_uscore; sum_uscore.reserve(nsr);
  std::vector<double> sum_vscore; sum_vscore.reserve(nsr);
  std::vector<double> sum_logRank; sum_logRank.reserve(nsr);
  std::vector<unsigned char> sum_rejPerStage; sum_rejPerStage.reserve(nsr);
  std::vector<unsigned char> sum_futPerStage; sum_futPerStage.reserve(nsr);

  // final raw containers
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

  // flatten
  for (size_t iter = 0; iter < maxIters; ++iter) {
    const auto& srows = results[iter].summaryRows;
    for (const auto& r : srows) {
      sum_iterNum.push_back(r.iterNum);
      sum_evNotArch.push_back(r.evNotAch);
      sum_stopStage.push_back(r.stopStage);
      sum_stageNum.push_back(r.stageNum);
      sum_analysisT.push_back(r.analysisT);
      sum_accruals1.push_back(r.accruals1);
      sum_accruals2.push_back(r.accruals2);
      sum_totAccruals.push_back(r.totAccruals);
      sum_events1.push_back(r.events1);
      sum_events2.push_back(r.events2);
      sum_totEvents.push_back(r.totEvents);
      sum_dropouts1.push_back(r.dropouts1);
      sum_dropouts2.push_back(r.dropouts2);
      sum_totDropouts.push_back(r.totDropouts);
      sum_uscore.push_back(r.uscore);
      sum_vscore.push_back(r.vscore);
      sum_logRank.push_back(r.logRank);
      sum_rejPerStage.push_back(r.rejPerStage);
      sum_futPerStage.push_back(r.futPerStage);
    }

    if (iter < maxRawIters) {
      const auto& rraw = results[iter].rawRows;
      for (const auto& rr : rraw) {
        raw_iterNum.push_back(rr.iterNum);
        raw_stopStage.push_back(rr.stopStage);
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
      }
    }
  }

  // compute per-stage simulation summaries
  size_t index2 = sum_iterNum.size();
  double niters = 0.0;
  std::vector<double> stopPerStage(K);
  std::vector<double> rejPerStage(K), futPerStage(K);
  std::vector<double> eventsPerStage(K), dropoutsPerStage(K);
  std::vector<double> subjectsPerStage(K), analysisTimePerStage(K);
  for (size_t i = 0; i < index2; ++i) {
    int k = sum_stageNum[i] - 1;
    if (sum_stageNum[i] == sum_stopStage[i]) {
      niters += 1.0;
      stopPerStage[k] += 1.0;
      rejPerStage[k] += sum_rejPerStage[i];
      futPerStage[k] += sum_futPerStage[i];
      eventsPerStage[k] += sum_totEvents[i];
      dropoutsPerStage[k] += sum_totDropouts[i];
      subjectsPerStage[k] += sum_totAccruals[i];
      analysisTimePerStage[k] += sum_analysisT[i];
    }
  }

  for (size_t k = 0; k < K; ++k) {
    rejPerStage[k] /= niters;
    futPerStage[k] /= niters;

    if (stopPerStage[k] > 0.0) {
      eventsPerStage[k] /= stopPerStage[k];
      dropoutsPerStage[k] /= stopPerStage[k];
      subjectsPerStage[k] /= stopPerStage[k];
      analysisTimePerStage[k] /= stopPerStage[k];
    } else {
      eventsPerStage[k] = 0.0;
      dropoutsPerStage[k] = 0.0;
      subjectsPerStage[k] = 0.0;
      analysisTimePerStage[k] = 0.0;
    }
  }

  // cumulative probabilities of rejection and futility by stage
  std::vector<double> cpu(K), cpl(K);
  std::partial_sum(rejPerStage.begin(), rejPerStage.end(), cpu.begin());
  std::partial_sum(futPerStage.begin(), futPerStage.end(), cpl.begin());

  // overall probability of rejection
  double overallReject = cpu[K - 1];

  // expected number of events, dropouts, subjects, study duration at trial end
  double expNumEvents = 0.0;
  double expNumDropouts = 0.0;
  double expNumSubjects = 0.0;
  double expStudyDur = 0.0;
  for (size_t i = 0; i < index2; ++i) {
    if (sum_stageNum[i] == sum_stopStage[i]) {
      expNumEvents += sum_totEvents[i];
      expNumDropouts += sum_totDropouts[i];
      expNumSubjects += sum_totAccruals[i];
      expStudyDur += sum_analysisT[i];
    }
  }
  expNumEvents /= niters;
  expNumDropouts /= niters;
  expNumSubjects /= niters;
  expStudyDur /= niters;

  // construct final output
  ListCpp overview;
  overview.push_back(std::move(rejPerStage), "rejectPerStage");
  overview.push_back(std::move(futPerStage), "futilityPerStage");
  overview.push_back(std::move(cpu), "cumulativeRejection");
  overview.push_back(std::move(cpl), "cumulativeFutility");
  overview.push_back(std::move(eventsPerStage), "numberOfEvents");
  overview.push_back(std::move(dropoutsPerStage), "numberOfDropouts");
  overview.push_back(std::move(subjectsPerStage), "numberOfSubjects");
  overview.push_back(std::move(analysisTimePerStage), "analysisTime");
  overview.push_back(overallReject, "overallReject");
  overview.push_back(expNumEvents, "expectedNumberOfEvents");
  overview.push_back(expNumDropouts, "expectedNumberOfDropouts");
  overview.push_back(expNumSubjects, "expectedNumberOfSubjects");
  overview.push_back(expStudyDur, "expectedStudyDuration");
  overview.push_back(hazardRatioH0, "hazardRatioH0");
  overview.push_back(useEvents, "useEvents");
  overview.push_back(n, "n");
  overview.push_back(fixedFollowup, "fixedFollowup");
  overview.push_back(rho1, "rho1");
  overview.push_back(rho2, "rho2");
  overview.push_back(kMax, "kMax");

  // Build summary DataFrameCpp
  DataFrameCpp sumdata;
  sumdata.push_back(std::move(sum_iterNum), "iterationNumber");
  sumdata.push_back(std::move(sum_stopStage), "stopStage");
  sumdata.push_back(std::move(sum_evNotArch), "eventsNotAchieved");
  sumdata.push_back(std::move(sum_stageNum), "stageNumber");
  sumdata.push_back(std::move(sum_analysisT), "analysisTime");
  sumdata.push_back(std::move(sum_accruals1), "accruals1");
  sumdata.push_back(std::move(sum_accruals2), "accruals2");
  sumdata.push_back(std::move(sum_totAccruals), "totalAccruals");
  sumdata.push_back(std::move(sum_events1), "events1");
  sumdata.push_back(std::move(sum_events2), "events2");
  sumdata.push_back(std::move(sum_totEvents), "totalEvents");
  sumdata.push_back(std::move(sum_dropouts1), "dropouts1");
  sumdata.push_back(std::move(sum_dropouts2), "dropouts2");
  sumdata.push_back(std::move(sum_totDropouts), "totalDropouts");
  sumdata.push_back(std::move(sum_uscore), "uscore");
  sumdata.push_back(std::move(sum_vscore), "vscore");
  sumdata.push_back(std::move(sum_logRank), "logRankStatistic");
  sumdata.push_back(std::move(sum_rejPerStage), "rejectPerStage");
  sumdata.push_back(std::move(sum_futPerStage), "futilityPerStage");

  ListCpp result;
  result.push_back(overview, "overview");
  result.push_back(sumdata, "sumdata");

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

    result.push_back(rawdata, "rawdata");
  }

  return result;
}


//' @title Log-Rank Test Simulation
//' @description Performs simulation for two-arm group sequential
//' trials based on weighted log-rank test.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates in terms of number
//'   of events for the conventional log-rank test and in terms of
//'   the actual information for weighted log-rank tests.
//'   Fixed prior to the trial. If left unspecified, it defaults to
//'   \code{plannedEvents / plannedEvents[kMax]} when \code{plannedEvents}
//'   is provided and to \code{plannedTime / plannedTime[kMax]} otherwise.
//' @inheritParams param_criticalValues
//' @inheritParams param_futilityBounds
//' @inheritParams param_hazardRatioH0
//' @param allocation1 Number of subjects in the active treatment group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in the control group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @param n Sample size.
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @param plannedEvents The planned cumulative total number of events at
//'   each stage.
//' @param plannedTime The calendar times for the analyses. To use calendar
//'   time to plan the analyses, \code{plannedEvents} should be missing.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
//'   stage to extract.
//' @param seed The seed to reproduce the simulation results.
//'
//' @return An S3 class \code{lrsim} object with 3 components:
//'
//' * \code{overview}: A list containing the following information:
//'
//'     - \code{rejectPerStage}: The efficacy stopping probability by stage.
//'
//'     - \code{futilityPerStage}: The futility stopping probability by
//'       stage.
//'
//'     - \code{cumulativeRejection}: Cumulative efficacy stopping
//'       probability by stage.
//'
//'     - \code{cumulativeFutility}: The cumulative futility stopping
//'       probability by stage.
//'
//'     - \code{numberOfEvents}: The average number of events by stage.
//'
//'     - \code{numberOfDropouts}: The average number of dropouts by stage.
//'
//'     - \code{numberOfSubjects}: The average number of subjects by stage.
//'
//'     - \code{analysisTime}: The average analysis time by stage.
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events for
//'       the overall study.
//'
//'     - \code{expectedNumberOfDropouts}: The expected number of dropouts
//'       for the overall study.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects
//'       for the overall study.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{hazardRatioH0}: Hazard ratio under the null hypothesis for
//'       the active treatment versus control.
//'
//'     - \code{useEvents}: whether the analyses are planned
//'       based on the number of events or calendar time.
//'
//'     - \code{n}: Sample size.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//'     - \code{rho1}: The first parameter of the Fleming-Harrington family
//'       of weighted log-rank test. Defaults to 0 for conventional log-rank
//'       test.
//'
//'     - \code{rho2}: The second parameter of the Fleming-Harrington family
//'       of weighted log-rank test. Defaults to 0 for conventional log-rank
//'       test.
//'
//'     - \code{kMax}: The maximum number of stages.
//'
//' * \code{sumdata}: A data frame of summary data by iteration and stage:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{eventsNotAchieved}: Whether the final target number of events
//'       is not achieved for the iteration.
//'
//'     - \code{stopStage}: The stage at which the trial stops.
//'
//'     - \code{stageNumber}: The stage number, covering all stages even if
//'       the trial stops at an interim look.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{accruals1}: The number of subjects enrolled at the stage for
//'       the treatment group.
//'
//'     - \code{accruals2}: The number of subjects enrolled at the stage for
//'       the control group.
//'
//'     - \code{totalAccruals}: The total number of subjects enrolled at
//'       the stage.
//'
//'     - \code{events1}: The number of events at the stage for
//'       the treatment group.
//'
//'     - \code{events2}: The number of events at the stage for
//'       the control group.
//'
//'     - \code{totalEvents}: The total number of events at the stage.
//'
//'     - \code{dropouts1}: The number of dropouts at the stage for
//'       the treatment group.
//'
//'     - \code{dropouts2}: The number of dropouts at the stage for
//'       the control group.
//'
//'     - \code{totalDropouts}: The total number of dropouts at the stage.
//'
//'     - \code{uscore}: The numerator of the log-rank test statistic.
//'
//'     - \code{vscore}: The variance of the log-rank test statistic.
//'
//'     - \code{logRankStatistic}: The log-rank test Z-statistic.
//'
//'     - \code{rejectPerStage}: Whether to reject the null hypothesis
//'       at the stage.
//'
//'     - \code{futilityPerStage}: Whether to stop the trial for futility
//'       at the stage.
//'
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//'   positive integer): A data frame for subject-level data for selected
//'   replications, containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stopStage}: The stage at which the trial stops.
//'
//'     - \code{stageNumber}: The stage number, covering all stages even if
//'       the trial stops at an interim look.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{subjectId}: The subject ID.
//'
//'     - \code{arrivalTime}: The enrollment time for the subject.
//'
//'     - \code{stratum}: The stratum for the subject.
//'
//'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
//'       subject.
//'
//'     - \code{survivalTime}: The underlying survival time for the subject.
//'
//'     - \code{dropoutTime}: The underlying dropout time for the subject.
//'
//'     - \code{timeUnderObservation}: The time under observation
//'       since randomization.
//'
//'     - \code{event}: Whether the subject experienced the event.
//'
//'     - \code{dropoutEvent}: Whether the subject dropped out.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Example 1: analyses based on number of events
//'
//' sim1 = lrsim(kMax = 2, informationRates = c(0.5, 1),
//'              criticalValues = c(2.797, 1.977),
//'              accrualIntensity = 11,
//'              lambda1 = 0.018, lambda2 = 0.030,
//'              n = 132,
//'              plannedEvents = c(60, 120),
//'              maxNumberOfIterations = 1000,
//'              maxNumberOfRawDatasetsPerStage = 1,
//'              seed = 314159)
//'
//' # summary statistics
//' sim1
//'
//' # summary for each simulated data set
//' head(sim1$sumdata)
//'
//' # raw data for selected replication
//' head(sim1$rawdata)
//'
//'
//' # Example 2: analyses based on calendar time have similar power
//'
//' sim2 = lrsim(kMax = 2, informationRates = c(0.5, 1),
//'              criticalValues = c(2.797, 1.977),
//'              accrualIntensity = 11,
//'              lambda1 = 0.018, lambda2 = 0.030,
//'              n = 132,
//'              plannedTime = c(31.9, 113.2),
//'              maxNumberOfIterations = 1000,
//'              maxNumberOfRawDatasetsPerStage = 1,
//'              seed = 314159)
//'
//' # summary statistics
//' sim2
//'
//' # summary for each simulated data set
//' head(sim2$sumdata)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List lrsim(
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const Rcpp::NumericVector& futilityBounds = NA_REAL,
    const double hazardRatioH0 = 1,
    const int allocation1 = 1,
    const int allocation2 = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const int n = NA_INTEGER,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = 0,
    const double rho1 = 0,
    const double rho2 = 0,
    const Rcpp::IntegerVector& plannedEvents = NA_INTEGER,
    const Rcpp::NumericVector& plannedTime = NA_REAL,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = 0) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);
  auto plannedE = Rcpp::as<std::vector<int>>(plannedEvents);
  auto plannedT = Rcpp::as<std::vector<double>>(plannedTime);

  auto out = lrsimcpp(
    kMax, infoRates, critValues, futBounds, hazardRatioH0,
    allocation1, allocation2, accrualT, accrualInt,
    pwSurvT, stratumFrac, lam1, lam2, gam1, gam2,
    n, followupTime, fixedFollowup, rho1, rho2, plannedE, plannedT,
    maxNumberOfIterations, maxNumberOfRawDatasetsPerStage, seed);

  thread_utils::drain_thread_warnings_to_R();

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "lrsim";

  return result;
}



// Parallel entry function
ListCpp lrsim3acpp(
    const int kMax,
    const double hazardRatioH013,
    const double hazardRatioH023,
    const double hazardRatioH012,
    const int allocation1,
    const int allocation2,
    const int allocation3,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& lambda3,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const std::vector<double>& gamma3,
    const int n,
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
  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
  size_t K = static_cast<size_t>(kMax);

  // decide planning mode
  bool useEvents;
  if (none_na(plannedEvents)) {
    useEvents = true;
    if (plannedEvents[0] <= 0)
      throw std::invalid_argument("plannedEvents must be positive");
    if (plannedEvents.size() != K)
      throw std::invalid_argument("Invalid length for plannedEvents");
    if (any_nonincreasing(plannedEvents))
      throw std::invalid_argument("plannedEvents must be increasing");
  } else if (none_na(plannedTime)) {
    useEvents = false;
    if (plannedTime[0] <= 0.0)
      throw std::invalid_argument("plannedTime must be positive");
    if (plannedTime.size() != K)
      throw std::invalid_argument("Invalid length for plannedTime");
    if (any_nonincreasing(plannedTime))
      throw std::invalid_argument("plannedTime must be increasing");
  } else {
    throw std::invalid_argument("Either plannedEvents or plannedTime must be given");
  }

  // validate other input parameters
  if (hazardRatioH013 <= 0.0 || hazardRatioH023 <= 0.0 || hazardRatioH012 <= 0.0)
    throw std::invalid_argument("hazardRatioH0 parameters must be positive");
  if (allocation1 < 1 || allocation2 < 1 || allocation3 < 1)
    throw std::invalid_argument("allocations must be positive integers");
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
  if (!none_na(lambda1)) throw std::invalid_argument("lambda1 must be provided");
  if (!none_na(lambda2)) throw std::invalid_argument("lambda2 must be provided");
  if (!none_na(lambda3)) throw std::invalid_argument("lambda3 must be provided");
  for (double v : lambda1) {
    if (v < 0.0) throw std::invalid_argument("lambda1 must be non-negative");
  }
  for (double v : lambda2) {
    if (v < 0.0) throw std::invalid_argument("lambda2 must be non-negative");
  }
  for (double v : lambda3) {
    if (v < 0.0) throw std::invalid_argument("lambda3 must be non-negative");
  }
  for (double v : gamma1) {
    if (v < 0.0) throw std::invalid_argument("gamma1 must be non-negative");
  }
  for (double v : gamma2) {
    if (v < 0.0) throw std::invalid_argument("gamma2 must be non-negative");
  }
  for (double v : gamma3) {
    if (v < 0.0) throw std::invalid_argument("gamma3 must be non-negative");
  }
  if (n == INT_MIN) throw std::invalid_argument("n must be provided");
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

  size_t N = static_cast<size_t>(n);
  size_t maxIters = static_cast<size_t>(maxNumberOfIterations);
  size_t maxRawIters = static_cast<size_t>(maxNumberOfRawDatasetsPerStage);
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  const std::vector<double>& tau = piecewiseSurvivalTime;
  const double fu = followupTime;

  // expand stratified inputs
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto lambda3x = expand_stratified(lambda3, nstrata, nintv, "lambda3");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");
  auto gamma3x = expand_stratified(gamma3, nstrata, nintv, "gamma3");

  // generate seeds for each iteration to ensure reproducibility
  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter) seeds[iter] = master_rng();

  // One summary (stage-level) row produced by an iteration
  struct StageSummaryRow {
    int iterNum = 0;
    unsigned char evNotAch = 0;
    int stageNum = 0;
    double analysisT = 0.0;
    int accruals1 = 0, accruals2 = 0, accruals3 = 0, totAccruals = 0;
    int events1 = 0, events2 = 0, events3 = 0, totEvents = 0;
    int dropouts1 = 0, dropouts2 = 0, dropouts3 = 0, totDropouts = 0;
    double uscore13 = 0.0, vscore13 = 0.0, logRank13 = 0.0;
    double uscore23 = 0.0, vscore23 = 0.0, logRank23 = 0.0;
    double uscore12 = 0.0, vscore12 = 0.0, logRank12 = 0.0;
  };

  // One subject-level (raw) row for a particular iteration and stage
  struct RawDatasetRow {
    int iterNum = 0, stageNum = 0;
    double analysisT = 0.0;
    int subjectId = 0;
    double arrivalT = 0.0;
    int stratum = 0, trtGrp = 0;
    double survivalT = 0.0, dropoutT = 0.0, timeObs = 0.0;
    unsigned char event = 0, dropEv = 0;
  };

  // Per-iteration container written exclusively by the worker thread
  struct IterationResult {
    std::vector<StageSummaryRow> summaryRows;
    std::vector<RawDatasetRow> rawRows; // populated only for first M iterations
    void reserveForSummary(size_t approxRows) { summaryRows.reserve(approxRows); }
    void reserveForRaw(size_t approxRows) { rawRows.reserve(approxRows); }
  };

  // pre-size per-iteration results
  std::vector<IterationResult> results;
  results.resize(maxIters);


  // Worker that runs simulation iterations [begin, end)
  struct SimWorker : public RcppParallel::Worker {
    // inputs (const refs)
    const size_t K;
    const double hazardRatioH013;
    const double hazardRatioH023;
    const double hazardRatioH012;
    const int allocation1;
    const int allocation2;
    const int allocation3;
    const std::vector<double>& accrualTime;
    const std::vector<double>& accrualIntensity;
    const std::vector<double>& tau;
    const std::vector<double>& stratumFraction;
    const std::vector<std::vector<double>>& lambda1x;
    const std::vector<std::vector<double>>& lambda2x;
    const std::vector<std::vector<double>>& lambda3x;
    const std::vector<std::vector<double>>& gamma1x;
    const std::vector<std::vector<double>>& gamma2x;
    const std::vector<std::vector<double>>& gamma3x;
    const size_t N;
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
    const size_t nstrata;

    // output pointer (pre-sized vector of IterationResult)
    std::vector<IterationResult>* results;

    SimWorker(
      size_t K_,
      double hazardRatioH013_,
      double hazardRatioH023_,
      double hazardRatioH012_,
      int allocation1_,
      int allocation2_,
      int allocation3_,
      const std::vector<double>& accrualTime_,
      const std::vector<double>& accrualIntensity_,
      const std::vector<double>& tau_,
      const std::vector<double>& stratumFraction_,
      const std::vector<std::vector<double>>& lambda1x_,
      const std::vector<std::vector<double>>& lambda2x_,
      const std::vector<std::vector<double>>& lambda3x_,
      const std::vector<std::vector<double>>& gamma1x_,
      const std::vector<std::vector<double>>& gamma2x_,
      const std::vector<std::vector<double>>& gamma3x_,
      size_t N_,
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
      size_t nstrata_,
      std::vector<IterationResult>* results_)
      : K(K_),
        hazardRatioH013(hazardRatioH013_),
        hazardRatioH023(hazardRatioH023_),
        hazardRatioH012(hazardRatioH012_),
        allocation1(allocation1_),
        allocation2(allocation2_),
        allocation3(allocation3_),
        accrualTime(accrualTime_),
        accrualIntensity(accrualIntensity_),
        tau(tau_),
        stratumFraction(stratumFraction_),
        lambda1x(lambda1x_),
        lambda2x(lambda2x_),
        lambda3x(lambda3x_),
        gamma1x(gamma1x_),
        gamma2x(gamma2x_),
        gamma3x(gamma3x_),
        N(N_),
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
        nstrata(nstrata_),
        results(results_)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      // local buffers reused by this worker
      std::vector<int> stratum(N), trtGrp(N);
      std::vector<double> arrivalT(N), survivalT(N), dropoutT(N);
      std::vector<double> timeObs(N), totalT(N);
      std::vector<unsigned char> event(N), dropEv(N);
      std::vector<int> b1(nstrata), b2(nstrata), b3(nstrata);
      std::vector<int> n1(nstrata), n2(nstrata), n3(nstrata);
      std::vector<double> km13(nstrata), km23(nstrata), km12(nstrata);
      std::vector<double> cumF(nstrata);
      std::partial_sum(stratumFraction.begin(), stratumFraction.end(), cumF.begin());

      std::vector<double> analysisT(K);
      std::vector<double> totalte; totalte.reserve(N);
      std::vector<size_t> sub; sub.reserve(N);

      for (size_t iter = begin; iter < end; ++iter) {
        // deterministic per-iteration RNG
        boost::random::mt19937_64 rng_local(seeds[iter]);
        boost::random::uniform_real_distribution<double> unif(0.0, 1.0);

        // per-iteration output container
        IterationResult& out = (*results)[iter];
        out.summaryRows.clear();
        out.rawRows.clear();
        if (iter < maxRawIters) out.reserveForRaw(K * N);
        out.reserveForSummary(K);

        // reset block randomization
        std::fill(b1.begin(), b1.end(), allocation1);
        std::fill(b2.begin(), b2.end(), allocation2);
        std::fill(b3.begin(), b3.end(), allocation3);

        double enrollt = 0.0;

        // generate cohort
        for (size_t i = 0; i < N; ++i) {
          double u = unif(rng_local);
          enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt);
          arrivalT[i] = enrollt;

          u = unif(rng_local);
          int j = findInterval1(u, cumF);
          stratum[i] = j + 1;

          // stratified block randomization among 3 arms
          u = unif(rng_local);
          double denom = static_cast<double>(b1[j] + b2[j] + b3[j]);
          double p1 = static_cast<double>(b1[j]) / denom;
          double p2 = static_cast<double>(b1[j] + b2[j]) / denom;
          if (u <= p1) { trtGrp[i] = 1; --b1[j]; }
          else if (u <= p2) { trtGrp[i] = 2; --b2[j]; }
          else { trtGrp[i] = 3; --b3[j]; }
          if (b1[j] + b2[j] + b3[j] == 0) {
            b1[j] = allocation1; b2[j] = allocation2; b3[j] = allocation3;
          }

          // survival time
          u = unif(rng_local);
          if (trtGrp[i] == 1) survivalT[i] = qtpwexpcpp1(u, tau, lambda1x[j]);
          else if (trtGrp[i] == 2) survivalT[i] = qtpwexpcpp1(u, tau, lambda2x[j]);
          else survivalT[i] = qtpwexpcpp1(u, tau, lambda3x[j]);

          // dropout time
          u = unif(rng_local);
          if (trtGrp[i] == 1) dropoutT[i] = qtpwexpcpp1(u, tau, gamma1x[j]);
          else if (trtGrp[i] == 2) dropoutT[i] = qtpwexpcpp1(u, tau, gamma2x[j]);
          else dropoutT[i] = qtpwexpcpp1(u, tau, gamma3x[j]);

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

        // determine analysis times (events counted only for arms 1 and 3)
        size_t nstages = K;
        bool evNotAch = false;

        if (useEvents) {
          totalte.clear();
          int nevents = 0; // events involving arm1 or arm3 in this iteration
          for (size_t i = 0; i < N; ++i) {
            if (event[i] && (trtGrp[i] == 1 || trtGrp[i] == 3)) {
              ++nevents; totalte.push_back(totalT[i]);
            }
          }
          if (nevents == 0) {
            thread_utils::push_thread_warning(
              std::string("No events for iteration ") + std::to_string(iter + 1) +
                " skipping this iteration.");
            // leave out.summaryRows empty to signal skipped iteration
            out.summaryRows.clear();
            out.rawRows.clear();
            continue;
          }
          std::sort(totalte.begin(), totalte.end());

          size_t j;
          for (j = 0; j < K; ++j) {
            if (plannedEvents[j] >= nevents) { nstages = j + 1; break; }
          }

          if (j == K) {
            for (size_t k = 0; k < nstages; ++k) {
              analysisT[k] = totalte[plannedEvents[k] - 1] + 1e-12;
            }
          } else {
            for (size_t k = 0; k < nstages - 1; ++k) {
              analysisT[k] = totalte[plannedEvents[k] - 1] + 1e-12;
            }
            analysisT[nstages - 1] = totalte.back() + 1e-12;
          }
          evNotAch = (nevents < plannedEvents[K - 1]);
        } else {
          std::copy_n(plannedTime.begin(), K, analysisT.begin());
          evNotAch = false;
        }

        // per-stage calculations
        for (size_t k = 0; k < nstages; ++k) {
          double time = analysisT[k];

          // reset counts
          std::fill(n1.begin(), n1.end(), 0);
          std::fill(n2.begin(), n2.end(), 0);
          std::fill(n3.begin(), n3.end(), 0);

          int events1 = 0, events2 = 0, events3 = 0;
          int dropouts1 = 0, dropouts2 = 0, dropouts3 = 0;

          // censoring at analysis time and count accruals/events/dropouts
          for (size_t i = 0; i < N; ++i) {
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
            if (trtGrp[i] == 1) {
              ++n1[h];
              if (event[i]) ++events1; else if (dropEv[i]) ++dropouts1;
            } else if (trtGrp[i] == 2) {
              ++n2[h];
              if (event[i]) ++events2; else if (dropEv[i]) ++dropouts2;
            } else {
              ++n3[h];
              if (event[i]) ++events3; else if (dropEv[i]) ++dropouts3;
            }
          }

          int accruals1 = std::accumulate(n1.begin(), n1.end(), 0);
          int accruals2 = std::accumulate(n2.begin(), n2.end(), 0);
          int accruals3 = std::accumulate(n3.begin(), n3.end(), 0);
          int totAccruals = accruals1 + accruals2 + accruals3;
          int totEvents = events1 + events2 + events3;
          int totDropouts = dropouts1 + dropouts2 + dropouts3;

          // collect indices with positive observed time and sort them
          sub.clear();
          for (size_t i = 0; i < N; ++i) if (timeObs[i] > 0.0) sub.push_back(i);
          std::sort(sub.begin(), sub.end(), [&](size_t a, size_t b) {
            return timeObs[a] < timeObs[b];
          });

          // compute stratified log-rank for three pairwise comparisons
          std::fill(km13.begin(), km13.end(), 1.0);
          std::fill(km23.begin(), km23.end(), 1.0);
          std::fill(km12.begin(), km12.end(), 1.0);

          double us13 = 0.0, vs13 = 0.0;
          double us23 = 0.0, vs23 = 0.0;
          double us12 = 0.0, vs12 = 0.0;

          for (size_t i = 0; i < sub.size(); ++i) {
            size_t idx = sub[i];
            size_t h = static_cast<size_t>(stratum[idx] - 1);

            double n1h = static_cast<double>(n1[h]);
            double n2h = static_cast<double>(n2[h]);
            double n3h = static_cast<double>(n3[h]);

            double n13a = n1h * hazardRatioH013;
            double n23a = n2h * hazardRatioH023;
            double n12a = n1h * hazardRatioH012;

            double nt13 = n1h + n3h;
            double nt23 = n2h + n3h;
            double nt12 = n1h + n2h;

            double nt13a = n13a + n3h;
            double nt23a = n23a + n3h;
            double nt12a = n12a + n2h;

            if (event[idx]) {
              // pair 1 vs 3
              if (trtGrp[idx] == 1 || trtGrp[idx] == 3) {
                double wh = 1.0;
                if (rho1 != 0.0 || rho2 != 0.0) {
                  wh = std::pow(km13[h], rho1) * std::pow(1.0 - km13[h], rho2);
                  km13[h] *= (1.0 - 1.0 / nt13);
                }
                double treated = (trtGrp[idx] == 1 ? 1.0 : 0.0);
                us13 += wh * (treated - n13a / nt13a);
                vs13 += wh * wh * n13a * n3h / (nt13a * nt13a);
              }

              // pair 2 vs 3
              if (trtGrp[idx] == 2 || trtGrp[idx] == 3) {
                double wh = 1.0;
                if (rho1 != 0.0 || rho2 != 0.0) {
                  wh = std::pow(km23[h], rho1) * std::pow(1.0 - km23[h], rho2);
                  km23[h] *= (1.0 - 1.0 / nt23);
                }
                double treated = (trtGrp[idx] == 2 ? 1.0 : 0.0);
                us23 += wh * (treated - n23a / nt23a);
                vs23 += wh * wh * n23a * n3h / (nt23a * nt23a);
              }

              // pair 1 vs 2
              if (trtGrp[idx] == 1 || trtGrp[idx] == 2) {
                double wh = 1.0;
                if (rho1 != 0.0 || rho2 != 0.0) {
                  wh = std::pow(km12[h], rho1) * std::pow(1.0 - km12[h], rho2);
                  km12[h] *= (1.0 - 1.0 / nt12);
                }
                double treated = (trtGrp[idx] == 1 ? 1.0 : 0.0);
                us12 += wh * (treated - n12a / nt12a);
                vs12 += wh * wh * n12a * n2h / (nt12a * nt12a);
              }
            } // event[idx]

            // reduce risk set
            if (trtGrp[idx] == 1) --n1[h];
            else if (trtGrp[idx] == 2) --n2[h];
            else --n3[h];
          } // end events loop

          double z13 = (vs13 > 0.0 ? us13 / std::sqrt(vs13) : 0.0);
          double z23 = (vs23 > 0.0 ? us23 / std::sqrt(vs23) : 0.0);
          double z12 = (vs12 > 0.0 ? us12 / std::sqrt(vs12) : 0.0);

          // optionally append raw rows for this stage
          if (iter < maxRawIters) { // only for first maxRawIters iterations
            for (size_t i = 0; i < N; ++i) {
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
          StageSummaryRow sr;
          sr.iterNum = static_cast<int>(iter + 1);
          sr.evNotAch = evNotAch ? 1 : 0;
          sr.stageNum = static_cast<int>(k + 1);
          sr.analysisT = time;
          sr.accruals1 = accruals1;
          sr.accruals2 = accruals2;
          sr.accruals3 = accruals3;
          sr.totAccruals = totAccruals;
          sr.events1 = events1;
          sr.events2 = events2;
          sr.events3 = events3;
          sr.totEvents = totEvents;
          sr.dropouts1 = dropouts1;
          sr.dropouts2 = dropouts2;
          sr.dropouts3 = dropouts3;
          sr.totDropouts = totDropouts;
          sr.uscore13 = us13;
          sr.vscore13 = vs13;
          sr.logRank13 = z13;
          sr.uscore23 = us23;
          sr.vscore23 = vs23;
          sr.logRank23 = z23;
          sr.uscore12 = us12;
          sr.vscore12 = vs12;
          sr.logRank12 = z12;
          out.summaryRows.push_back(std::move(sr));
        } // per-stage
      } // per-iteration
    } // operator()
  }; // SimWorker

  // run worker in parallel
  SimWorker worker(
      K, hazardRatioH013, hazardRatioH023, hazardRatioH012,
      allocation1, allocation2, allocation3,
      accrualTime, accrualIntensity, tau, stratumFraction,
      lambda1x, lambda2x, lambda3x, gamma1x, gamma2x, gamma3x,
      N, fu, fixedFollowup, rho1, rho2,
      plannedEvents, plannedTime,
      maxIters, maxRawIters, seeds, useEvents, nstrata,
      &results
  );

  RcppParallel::parallelFor(0, maxIters, worker);

  // Flatten results
  size_t nsr = 0, nrr = 0;
  for (size_t iter = 0; iter < maxIters; ++iter) {
    nsr += results[iter].summaryRows.size();
    nrr += results[iter].rawRows.size();
  }
  if (nsr == 0) throw std::runtime_error(
    "No iterations with observed events for arm 1 or arm 3. "
    "Unable to produce output.");

  // prepare final containers (reserve capacities)
  std::vector<int> sum_iterNum; sum_iterNum.reserve(nsr);
  std::vector<unsigned char> sum_evNotArch; sum_evNotArch.reserve(nsr);
  std::vector<int> sum_stageNum; sum_stageNum.reserve(nsr);
  std::vector<double> sum_analysisT; sum_analysisT.reserve(nsr);
  std::vector<int> sum_accruals1; sum_accruals1.reserve(nsr);
  std::vector<int> sum_accruals2; sum_accruals2.reserve(nsr);
  std::vector<int> sum_accruals3; sum_accruals3.reserve(nsr);
  std::vector<int> sum_totAccruals; sum_totAccruals.reserve(nsr);
  std::vector<int> sum_events1; sum_events1.reserve(nsr);
  std::vector<int> sum_events2; sum_events2.reserve(nsr);
  std::vector<int> sum_events3; sum_events3.reserve(nsr);
  std::vector<int> sum_totEvents; sum_totEvents.reserve(nsr);
  std::vector<int> sum_dropouts1; sum_dropouts1.reserve(nsr);
  std::vector<int> sum_dropouts2; sum_dropouts2.reserve(nsr);
  std::vector<int> sum_dropouts3; sum_dropouts3.reserve(nsr);
  std::vector<int> sum_totDropouts; sum_totDropouts.reserve(nsr);
  std::vector<double> sum_uscore13; sum_uscore13.reserve(nsr);
  std::vector<double> sum_vscore13; sum_vscore13.reserve(nsr);
  std::vector<double> sum_logRank13; sum_logRank13.reserve(nsr);
  std::vector<double> sum_uscore23; sum_uscore23.reserve(nsr);
  std::vector<double> sum_vscore23; sum_vscore23.reserve(nsr);
  std::vector<double> sum_logRank23; sum_logRank23.reserve(nsr);
  std::vector<double> sum_uscore12; sum_uscore12.reserve(nsr);
  std::vector<double> sum_vscore12; sum_vscore12.reserve(nsr);
  std::vector<double> sum_logRank12; sum_logRank12.reserve(nsr);

  // raw final containers
  std::vector<int> raw_iterNum; raw_iterNum.reserve(nrr);
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
    const auto& srows = results[iter].summaryRows;
    for (const auto& r : srows) {
      sum_iterNum.push_back(r.iterNum);
      sum_evNotArch.push_back(r.evNotAch);
      sum_stageNum.push_back(r.stageNum);
      sum_analysisT.push_back(r.analysisT);
      sum_accruals1.push_back(r.accruals1);
      sum_accruals2.push_back(r.accruals2);
      sum_accruals3.push_back(r.accruals3);
      sum_totAccruals.push_back(r.totAccruals);
      sum_events1.push_back(r.events1);
      sum_events2.push_back(r.events2);
      sum_events3.push_back(r.events3);
      sum_totEvents.push_back(r.totEvents);
      sum_dropouts1.push_back(r.dropouts1);
      sum_dropouts2.push_back(r.dropouts2);
      sum_dropouts3.push_back(r.dropouts3);
      sum_totDropouts.push_back(r.totDropouts);
      sum_uscore13.push_back(r.uscore13);
      sum_vscore13.push_back(r.vscore13);
      sum_logRank13.push_back(r.logRank13);
      sum_uscore23.push_back(r.uscore23);
      sum_vscore23.push_back(r.vscore23);
      sum_logRank23.push_back(r.logRank23);
      sum_uscore12.push_back(r.uscore12);
      sum_vscore12.push_back(r.vscore12);
      sum_logRank12.push_back(r.logRank12);
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
      }
    }
  }

  // Build DataFrameCpp summary
  DataFrameCpp sumdata;
  sumdata.push_back(std::move(sum_iterNum), "iterationNumber");
  sumdata.push_back(std::move(sum_evNotArch), "eventsNotAchieved");
  sumdata.push_back(std::move(sum_stageNum), "stageNumber");
  sumdata.push_back(std::move(sum_analysisT), "analysisTime");
  sumdata.push_back(std::move(sum_accruals1), "accruals1");
  sumdata.push_back(std::move(sum_accruals2), "accruals2");
  sumdata.push_back(std::move(sum_accruals3), "accruals3");
  sumdata.push_back(std::move(sum_totAccruals), "totalAccruals");
  sumdata.push_back(std::move(sum_events1), "events1");
  sumdata.push_back(std::move(sum_events2), "events2");
  sumdata.push_back(std::move(sum_events3), "events3");
  sumdata.push_back(std::move(sum_totEvents), "totalEvents");
  sumdata.push_back(std::move(sum_dropouts1), "dropouts1");
  sumdata.push_back(std::move(sum_dropouts2), "dropouts2");
  sumdata.push_back(std::move(sum_dropouts3), "dropouts3");
  sumdata.push_back(std::move(sum_totDropouts), "totDropouts");
  sumdata.push_back(std::move(sum_uscore13), "uscore13");
  sumdata.push_back(std::move(sum_vscore13), "vscore13");
  sumdata.push_back(std::move(sum_logRank13), "logRankStatistic13");
  sumdata.push_back(std::move(sum_uscore23), "uscore23");
  sumdata.push_back(std::move(sum_vscore23), "vscore23");
  sumdata.push_back(std::move(sum_logRank23), "logRankStatistic23");
  sumdata.push_back(std::move(sum_uscore12), "uscore12");
  sumdata.push_back(std::move(sum_vscore12), "vscore12");
  sumdata.push_back(std::move(sum_logRank12), "logRankStatistic12");

  ListCpp result;
  result.push_back(sumdata, "sumdata");

  // attach raw data
  if (!raw_iterNum.empty()) {
    DataFrameCpp rawdata;
    rawdata.push_back(std::move(raw_iterNum), "iterationNumber");
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

    result.push_back(rawdata, "rawdata");
  }

  return result;
}


//' @title Log-Rank Test Simulation for Three Arms
//' @description Performs simulation for three-arm group sequential trials
//' based on weighted log-rank test. The looks are driven by the total
//' number of events in Arm A and Arm C combined. Alternatively,
//' the analyses can be planned to occur at specified calendar times.
//'
//' @inheritParams param_kMax
//' @param hazardRatioH013 Hazard ratio under the null hypothesis for arm 1
//'   versus arm 3. Defaults to 1 for superiority test.
//' @param hazardRatioH023 Hazard ratio under the null hypothesis for arm 2
//'   versus arm 3. Defaults to 1 for superiority test.
//' @param hazardRatioH012 Hazard ratio under the null hypothesis for arm 1
//'   versus arm 2. Defaults to 1 for superiority test.
//' @param allocation1 Number of subjects in Arm A in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in Arm B in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation3 Number of subjects in Arm C in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param lambda1 A vector of hazard rates for the event in each analysis
//'   time interval by stratum for arm 1.
//' @param lambda2 A vector of hazard rates for the event in each analysis
//'   time interval by stratum for arm 2.
//' @param lambda3 A vector of hazard rates for the event in each analysis
//'   time interval by stratum for arm 3.
//' @param gamma1 The hazard rate for exponential dropout. A vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for arm 1.
//' @param gamma2 The hazard rate for exponential dropout. A vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for arm 2.
//' @param gamma3 The hazard rate for exponential dropout. A vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for arm 3.
//' @param n Sample size.
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @param plannedEvents The planned cumulative total number of events at
//'   Look 1 to Look \code{kMax} for Arms A and C combined.
//' @param plannedTime The calendar times for the analyses. To use calendar
//'   time to plan the analyses, \code{plannedEvents} should be missing.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
//'   stage to extract.
//' @param seed The seed to reproduce the simulation results.
//'
//' @return A list with 2 components:
//'
//' * \code{sumdata}: A data frame of summary data by iteration and stage:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{eventsNotAchieved}: Whether the target number of events
//'       is not achieved for the iteration.
//'
//'     - \code{stageNumber}: The stage number, covering all stages even if
//'       the trial stops at an interim look.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{accruals1}: The number of subjects enrolled at the stage for
//'       the active treatment 1 group.
//'
//'     - \code{accruals2}: The number of subjects enrolled at the stage for
//'       the active treatment 2 group.
//'
//'     - \code{accruals3}: The number of subjects enrolled at the stage for
//'       the control group.
//'
//'     - \code{totalAccruals}: The total number of subjects enrolled at
//'       the stage.
//'
//'     - \code{events1}: The number of events at the stage for
//'       the active treatment 1 group.
//'
//'     - \code{events2}: The number of events at the stage for
//'       the active treatment 2 group.
//'
//'     - \code{events3}: The number of events at the stage for
//'       the control group.
//'
//'     - \code{totalEvents}: The total number of events at the stage.
//'
//'     - \code{dropouts1}: The number of dropouts at the stage for
//'       the active treatment 1 group.
//'
//'     - \code{dropouts2}: The number of dropouts at the stage for
//'       the active treatment 2 group.
//'
//'     - \code{dropouts3}: The number of dropouts at the stage for
//'       the control group.
//'
//'     - \code{totalDropouts}: The total number of dropouts at the stage.
//'
//'     - \code{uscore13}: The log-rank test score statistic comparing the
//'       active treatment 1 to the control.
//'
//'     - \code{vscore13}: The log-rank test variance statistic comparing the
//'       active treatment 1 to the control.
//'
//'     - \code{logRankStatistic13}: The log-rank test Z-statistic
//'       comparing the active treatment 1 to the control.
//'
//'     - \code{uscore23}: The log-rank test score statistic comparing the
//'       active treatment 2 to the control.
//'
//'     - \code{vscore23}: The log-rank test variance statistic comparing the
//'       active treatment 2 to the control.
//'
//'     - \code{logRankStatistic23}: The log-rank test Z-statistic
//'       comparing the active treatment 2 to the control.
//'
//'     - \code{uscore12}: The log-rank test score statistic comparing the
//'       active treatment 1 to the active treatment 2.
//'
//'     - \code{vscore12}: The log-rank test variance statistic comparing the
//'       active treatment 1 to the active treatment 2.
//'
//'     - \code{logRankStatistic12}: The log-rank test Z-statistic
//'       comparing the active treatment 1 to the active treatment 2.
//'
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//'   positive integer): A data frame for subject-level data for selected
//'   replications, containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage under consideration.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{subjectId}: The subject ID.
//'
//'     - \code{arrivalTime}: The enrollment time for the subject.
//'
//'     - \code{stratum}: The stratum for the subject.
//'
//'     - \code{treatmentGroup}: The treatment group (1, 2, or 3) for
//'       the subject.
//'
//'     - \code{survivalTime}: The underlying survival time for the subject.
//'
//'     - \code{dropoutTime}: The underlying dropout time for the subject.
//'
//'     - \code{timeUnderObservation}: The time under observation
//'       since randomization for the subject.
//'
//'     - \code{event}: Whether the subject experienced the event.
//'
//'     - \code{dropoutEvent}: Whether the subject dropped out.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' sim1 = lrsim3a(
//'   kMax = 3,
//'   allocation1 = 2,
//'   allocation2 = 2,
//'   allocation3 = 1,
//'   accrualTime = c(0, 8),
//'   accrualIntensity = c(10, 28),
//'   piecewiseSurvivalTime = 0,
//'   lambda1 = log(2)/12*0.60,
//'   lambda2 = log(2)/12*0.70,
//'   lambda3 = log(2)/12,
//'   n = 700,
//'   plannedEvents = c(186, 259, 295),
//'   maxNumberOfIterations = 1000,
//'   maxNumberOfRawDatasetsPerStage = 1,
//'   seed = 314159)
//'
//' head(sim1$sumdata)
//' head(sim1$rawdata)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List lrsim3a(
    const int kMax = 1,
    const double hazardRatioH013 = 1,
    const double hazardRatioH023 = 1,
    const double hazardRatioH012 = 1,
    const int allocation1 = 1,
    const int allocation2 = 1,
    const int allocation3 = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& lambda3 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const Rcpp::NumericVector& gamma3 = 0,
    const int n = NA_INTEGER,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = 0,
    const double rho1 = 0,
    const double rho2 = 0,
    const Rcpp::IntegerVector& plannedEvents = NA_INTEGER,
    const Rcpp::NumericVector& plannedTime = NA_REAL,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = 0) {

  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto lam3 = Rcpp::as<std::vector<double>>(lambda3);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);
  auto gam3 = Rcpp::as<std::vector<double>>(gamma3);
  auto plannedE = Rcpp::as<std::vector<int>>(plannedEvents);
  auto plannedT = Rcpp::as<std::vector<double>>(plannedTime);

  auto out = lrsim3acpp(
    kMax, hazardRatioH013, hazardRatioH023, hazardRatioH012,
    allocation1, allocation2, allocation3, accrualT, accrualInt,
    pwSurvT, stratumFrac, lam1, lam2, lam3, gam1, gam2, gam3,
    n, followupTime, fixedFollowup, rho1, rho2, plannedE, plannedT,
    maxNumberOfIterations, maxNumberOfRawDatasetsPerStage, seed);

  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(out);
}


// The parallel entry function
ListCpp lrsim2ecpp(
    const int kMax,
    const int kMaxpfs,
    const double hazardRatioH0pfs,
    const double hazardRatioH0os,
    const int allocation1,
    const int allocation2,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const double rho_pd_os,
    const std::vector<double>& lambda1pfs,
    const std::vector<double>& lambda2pfs,
    const std::vector<double>& lambda1os,
    const std::vector<double>& lambda2os,
    const std::vector<double>& gamma1pfs,
    const std::vector<double>& gamma2pfs,
    const std::vector<double>& gamma1os,
    const std::vector<double>& gamma2os,
    const int n,
    const double followupTime,
    const bool fixedFollowup,
    const double rho1,
    const double rho2,
    const std::vector<int>& plannedEvents,
    const std::vector<double>& plannedTime,
    const int maxNumberOfIterations,
    const int maxNumberOfRawDatasetsPerStage,
    const int seed) {

  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
  size_t K = static_cast<size_t>(kMax);

  int kMaxpfsx = kMaxpfs;
  if (kMaxpfsx < 0) kMaxpfsx = kMax;
  if (kMaxpfsx > kMax)
    throw std::invalid_argument("kMaxpfs must be less than or equal to kMax");
  size_t Kpfs = static_cast<size_t>(kMaxpfsx);

  // whether to plan by events or calendar time
  bool useEvents;
  if (none_na(plannedEvents)) {
    useEvents = true;
    if (plannedEvents[0] <= 0)
      throw std::invalid_argument("plannedEvents must be positive");
    if (plannedEvents.size() != K)
      throw std::invalid_argument("Invalid length for plannedEvents");
    if (Kpfs > 1) {
      for (size_t i = 1; i < Kpfs; ++i) {
        if (plannedEvents[i] <= plannedEvents[i-1])
          throw std::invalid_argument("plannedEvents for PFS must be increasing");
      }
    }
    if (K - Kpfs > 1) {
      for (size_t i = Kpfs + 1; i < plannedEvents.size(); ++i) {
        if (plannedEvents[i] <= plannedEvents[i-1])
          throw std::invalid_argument("plannedEvents for OS must be increasing");
      }
    }
  } else if (none_na(plannedTime)) {
    useEvents = false;
    if (plannedTime[0] <= 0.0)
      throw std::invalid_argument("plannedTime must be positive");
    if (plannedTime.size() != K)
      throw std::invalid_argument("Invalid length for plannedTime");
    if (any_nonincreasing(plannedTime))
      throw std::invalid_argument("plannedTime must be increasing");
  } else {
    throw std::invalid_argument("Either plannedEvents or plannedTime must be given");
  }

  // validate other parameters
  if (hazardRatioH0pfs <= 0.0)
    throw std::invalid_argument("PFS hazard ratio under H0 must be positive");
  if (hazardRatioH0os <= 0.0)
    throw std::invalid_argument("OS hazard ratio under H0 must be positive");
  if (allocation1 < 1 || allocation2 < 1)
    throw std::invalid_argument("allocations must be positive integers");
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
  if (rho_pd_os <= -1.0 || rho_pd_os >= 1.0)
    throw std::invalid_argument("rho_pd_os must lie in (-1, 1)");
  if (!none_na(lambda1pfs)) throw std::invalid_argument("lambda1pfs must be provided");
  if (!none_na(lambda2pfs)) throw std::invalid_argument("lambda2pfs must be provided");
  if (!none_na(lambda1os)) throw std::invalid_argument("lambda1os must be provided");
  if (!none_na(lambda2os)) throw std::invalid_argument("lambda2os must be provided");
  for (double v : lambda1pfs) {
    if (v < 0.0) throw std::invalid_argument("lambda1pfs must be non-negative");
  }
  for (double v : lambda2pfs) {
    if (v < 0.0) throw std::invalid_argument("lambda2pfs must be non-negative");
  }
  for (double v : lambda1os) {
    if (v < 0.0) throw std::invalid_argument("lambda1os must be non-negative");
  }
  for (double v : lambda2os) {
    if (v < 0.0) throw std::invalid_argument("lambda2os must be non-negative");
  }
  for (double v : gamma1pfs) {
    if (v < 0.0) throw std::invalid_argument("gamma1pfs must be non-negative");
  }
  for (double v : gamma2pfs) {
    if (v < 0.0) throw std::invalid_argument("gamma2pfs must be non-negative");
  }
  for (double v : gamma1os) {
    if (v < 0.0) throw std::invalid_argument("gamma1os must be non-negative");
  }
  for (double v : gamma2os) {
    if (v < 0.0) throw std::invalid_argument("gamma2os must be non-negative");
  }
  if (n == INT_MIN) throw std::invalid_argument("n must be provided");
  if (n <= 0) throw std::invalid_argument("n must be positive");
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

  size_t N = static_cast<size_t>(n);
  size_t maxIters = static_cast<size_t>(maxNumberOfIterations);
  size_t maxRawIters = static_cast<size_t>(maxNumberOfRawDatasetsPerStage);
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  size_t nintv2 = (nintv == 1 ? 10u : nintv + 10u);
  const std::vector<double>& tau = piecewiseSurvivalTime;
  const double fu = followupTime;
  const double rho_pd_os_pyth_comp = std::sqrt(1 - rho_pd_os * rho_pd_os);

  // expand stratified inputs to nested vectors (per-stratum vectors)
  auto lambda1pfsx = expand_stratified(lambda1pfs, nstrata, nintv, "lambda1pfs");
  auto lambda2pfsx = expand_stratified(lambda2pfs, nstrata, nintv, "lambda2pfs");
  auto lambda1osx  = expand_stratified(lambda1os,  nstrata, nintv, "lambda1os");
  auto lambda2osx  = expand_stratified(lambda2os,  nstrata, nintv, "lambda2os");
  auto gamma1pfsx  = expand_stratified(gamma1pfs,  nstrata, nintv, "gamma1pfs");
  auto gamma2pfsx  = expand_stratified(gamma2pfs,  nstrata, nintv, "gamma2pfs");
  auto gamma1osx   = expand_stratified(gamma1os,   nstrata, nintv, "gamma1os");
  auto gamma2osx   = expand_stratified(gamma2os,   nstrata, nintv, "gamma2os");

  // compute P(D) hazard piecewise for each stratum using hazard_pdcpp
  std::vector<std::vector<double>> tau1pdx(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> tau2pdx(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> lambda1pd(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> lambda2pd(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> gamma1pd(nstrata, std::vector<double>(nintv));
  std::vector<std::vector<double>> gamma2pd(nstrata, std::vector<double>(nintv));

  for (size_t s = 0; s < nstrata; ++s) {
    // pass per-stratum vectors to hazard_pdcpp
    ListCpp a1 = hazard_pdcpp(tau, lambda1pfsx[s], lambda1osx[s], rho_pd_os);
    ListCpp a2 = hazard_pdcpp(tau, lambda2pfsx[s], lambda2osx[s], rho_pd_os);

    tau1pdx[s] = a1.get<std::vector<double>>("piecewiseSurvivalTime");
    tau2pdx[s] = a2.get<std::vector<double>>("piecewiseSurvivalTime");
    lambda1pd[s]  = a1.get<std::vector<double>>("hazard_pd");
    lambda2pd[s]  = a2.get<std::vector<double>>("hazard_pd");

    // gamma for pd is difference pfs - os
    for (size_t t = 0; t < nintv; ++t) {
      gamma1pd[s][t] = gamma1pfsx[s][t] - gamma1osx[s][t];
      gamma2pd[s][t] = gamma2pfsx[s][t] - gamma2osx[s][t];
    }
  }

  // generate seeds for each iteration to ensure reproducibility
  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter) seeds[iter] = master_rng();


  // Per-stage summary row
  struct StageSummaryRow {
    int iterNum = 0;
    unsigned char evNotAch1 = 0, evNotAch2 = 0;
    int stageNum = 0;
    double analysisT = 0.0;
    int accruals1 = 0, accruals2 = 0, totAccruals = 0;
    int endpt = 0; // 1 == PFS, 2 == OS
    int events1 = 0, events2 = 0, totEvents = 0;
    int dropouts1 = 0, dropouts2 = 0, totDropouts = 0;
    double uscore = 0.0, vscore = 0.0, logRank = 0.0;
  };

  // Per-subject raw row
  struct RawDatasetRow {
    int iterNum = 0, stageNum = 0;
    double analysisT = 0.0;
    int subjectId = 0;
    double arrivalT = 0.0;
    int stratum = 0, trtGrp = 0;
    int endpt = 0; // 1 == PFS, 2 == OS
    double survivalT = 0.0, dropoutT = 0.0, timeObs = 0.0;
    unsigned char event = 0, dropEv = 0;
  };

  // Per-iteration container
  struct IterationResult {
    std::vector<StageSummaryRow> summaryRows;
    std::vector<RawDatasetRow> rawRows;
    void reserveForSummary(size_t approxRows) { summaryRows.reserve(approxRows); }
    void reserveForRaw(size_t approxRows) { rawRows.reserve(approxRows); }
  };

  std::vector<IterationResult> results;
  results.resize(maxIters);

  // Worker struct defined in-function to avoid symbol collisions
  struct SimWorker : public RcppParallel::Worker {
    // inputs (captured by reference)
    const size_t K;
    const size_t Kpfs;
    const double hazardRatioH0pfs;
    const double hazardRatioH0os;
    const int allocation1;
    const int allocation2;
    const std::vector<double>& accrualTime;
    const std::vector<double>& accrualIntensity;
    const std::vector<double>& tau;
    const std::vector<double>& stratumFraction;
    const double rho_pd_os;
    const std::vector<std::vector<double>>& lambda1pfsx;
    const std::vector<std::vector<double>>& lambda2pfsx;
    const std::vector<std::vector<double>>& lambda1osx;
    const std::vector<std::vector<double>>& lambda2osx;
    const std::vector<std::vector<double>>& gamma1pfsx;
    const std::vector<std::vector<double>>& gamma2pfsx;
    const std::vector<std::vector<double>>& gamma1osx;
    const std::vector<std::vector<double>>& gamma2osx;
    const std::vector<std::vector<double>>& tau1pdx;
    const std::vector<std::vector<double>>& tau2pdx;
    const std::vector<std::vector<double>>& lambda1pd;
    const std::vector<std::vector<double>>& lambda2pd;
    const std::vector<std::vector<double>>& gamma1pd;
    const std::vector<std::vector<double>>& gamma2pd;

    const size_t N;
    const double fu;
    const bool fixedFollowup;
    const double rho1;
    const double rho2;
    const std::vector<int>& plannedEvents;
    const std::vector<double>& plannedTime;
    const size_t maxIters;
    const size_t maxRawIters;
    const std::vector<uint64_t>& seeds;
    const bool useEvents;
    const size_t nstrata;
    const double rho_pd_os_pyth_comp;

    // output pointer to local results
    std::vector<IterationResult>* results;

    SimWorker(
      size_t K_,
      size_t Kpfs_,
      double hazardRatioH0pfs_,
      double hazardRatioH0os_,
      int allocation1_,
      int allocation2_,
      const std::vector<double>& accrualTime_,
      const std::vector<double>& accrualIntensity_,
      const std::vector<double>& tau_,
      const std::vector<double>& stratumFraction_,
      double rho_pd_os_,
      const std::vector<std::vector<double>>& lambda1pfsx_,
      const std::vector<std::vector<double>>& lambda2pfsx_,
      const std::vector<std::vector<double>>& lambda1osx_,
      const std::vector<std::vector<double>>& lambda2osx_,
      const std::vector<std::vector<double>>& gamma1pfsx_,
      const std::vector<std::vector<double>>& gamma2pfsx_,
      const std::vector<std::vector<double>>& gamma1osx_,
      const std::vector<std::vector<double>>& gamma2osx_,
      const std::vector<std::vector<double>>& tau1pdx_,
      const std::vector<std::vector<double>>& tau2pdx_,
      const std::vector<std::vector<double>>& lambda1pd_,
      const std::vector<std::vector<double>>& lambda2pd_,
      const std::vector<std::vector<double>>& gamma1pd_,
      const std::vector<std::vector<double>>& gamma2pd_,
      size_t N_,
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
      size_t nstrata_,
      const double rho_pd_os_pyth_comp_,
      std::vector<IterationResult>* results_)
      : K(K_),
        Kpfs(Kpfs_),
        hazardRatioH0pfs(hazardRatioH0pfs_),
        hazardRatioH0os(hazardRatioH0os_),
        allocation1(allocation1_),
        allocation2(allocation2_),
        accrualTime(accrualTime_),
        accrualIntensity(accrualIntensity_),
        tau(tau_),
        stratumFraction(stratumFraction_),
        rho_pd_os(rho_pd_os_),
        lambda1pfsx(lambda1pfsx_),
        lambda2pfsx(lambda2pfsx_),
        lambda1osx(lambda1osx_),
        lambda2osx(lambda2osx_),
        gamma1pfsx(gamma1pfsx_),
        gamma2pfsx(gamma2pfsx_),
        gamma1osx(gamma1osx_),
        gamma2osx(gamma2osx_),
        tau1pdx(tau1pdx_),
        tau2pdx(tau2pdx_),
        lambda1pd(lambda1pd_),
        lambda2pd(lambda2pd_),
        gamma1pd(gamma1pd_),
        gamma2pd(gamma2pd_),
        N(N_),
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
        nstrata(nstrata_),
        rho_pd_os_pyth_comp(rho_pd_os_pyth_comp_),
        results(results_)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      // local buffers
      std::vector<int> stratum(N), trtGrp(N);
      std::vector<double> arrivalT(N), survivalT1(N), survivalT2(N);
      std::vector<double> dropoutT1(N), dropoutT2(N);
      std::vector<double> timeObs1(N), timeObs2(N);
      std::vector<double> totalT1(N), totalT2(N);
      std::vector<unsigned char> event1(N), event2(N);
      std::vector<unsigned char> dropEv1(N), dropEv2(N);

      std::vector<int> b1(nstrata), b2(nstrata);
      std::vector<int> n1(nstrata), n2(nstrata), n1x(nstrata), n2x(nstrata);
      std::vector<double> km(nstrata);
      std::vector<double> cumF(nstrata);
      std::partial_sum(stratumFraction.begin(), stratumFraction.end(), cumF.begin());

      std::vector<double> analysisT(K);
      std::vector<double> analysisT1; analysisT1.reserve(Kpfs);
      std::vector<double> analysisT2; analysisT2.reserve(K - Kpfs);
      std::vector<double> totalte1; totalte1.reserve(N);
      std::vector<double> totalte2; totalte2.reserve(N);
      std::vector<size_t> sub; sub.reserve(N);


      for (size_t iter = begin; iter < end; ++iter) {
        // RNG for this iteration
        boost::random::mt19937_64 rng_local(seeds[iter]);
        boost::random::uniform_real_distribution<double> unif(0.0, 1.0);
        boost::random::normal_distribution<double> norm(0.0, 1.0);

        IterationResult& out = (*results)[iter];
        out.summaryRows.clear();
        out.rawRows.clear();
        if (iter < maxRawIters) out.reserveForRaw(K * N);
        out.reserveForSummary(K * 2); // up to two endpoints per stage

        // reset blocks
        std::fill(b1.begin(), b1.end(), allocation1);
        std::fill(b2.begin(), b2.end(), allocation2);

        double enrollt = 0.0;

        // generate cohort
        for (size_t i = 0; i < N; ++i) {
          double u = unif(rng_local);
          enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt);
          arrivalT[i] = enrollt;

          u = unif(rng_local);
          int j = findInterval1(u, cumF);
          stratum[i] = j + 1;

          u = unif(rng_local);
          double denom = static_cast<double>(b1[j] + b2[j]);
          double p = static_cast<double>(b1[j]) / denom;
          if (u <= p) { trtGrp[i] = 1; --b1[j]; }
          else { trtGrp[i] = 2; --b2[j]; }
          if (b1[j] + b2[j] == 0) { b1[j] = allocation1; b2[j] = allocation2; }

          // correlated normals -> uniforms
          double z1 = norm(rng_local);
          double z2 = norm(rng_local);
          double u1 = boost_pnorm(z1);
          double u2 = boost_pnorm(rho_pd_os * z1 + rho_pd_os_pyth_comp * z2);

          if (trtGrp[i] == 1) {
            survivalT1[i] = qtpwexpcpp1(u1, tau1pdx[j], lambda1pd[j]);
            survivalT2[i] = qtpwexpcpp1(u2, tau, lambda1osx[j]);
          } else {
            survivalT1[i] = qtpwexpcpp1(u1, tau2pdx[j], lambda2pd[j]);
            survivalT2[i] = qtpwexpcpp1(u2, tau, lambda2osx[j]);
          }
          // PFS includes death
          if (survivalT1[i] > survivalT2[i]) survivalT1[i] = survivalT2[i];

          // dropout times
          u1 = unif(rng_local);
          u2 = unif(rng_local);
          if (trtGrp[i] == 1) {
            dropoutT1[i] = qtpwexpcpp1(u1, tau, gamma1pd[j]);
            dropoutT2[i] = qtpwexpcpp1(u2, tau, gamma1osx[j]);
          } else {
            dropoutT1[i] = qtpwexpcpp1(u1, tau, gamma2pd[j]);
            dropoutT2[i] = qtpwexpcpp1(u2, tau, gamma2osx[j]);
          }
          if (dropoutT1[i] > dropoutT2[i]) dropoutT1[i] = dropoutT2[i];

          // initial observed times/events
          double sv1 = survivalT1[i], sv2 = survivalT2[i];
          double dr1 = dropoutT1[i], dr2 = dropoutT2[i];
          if (fixedFollowup) {
            if (sv1 <= dr1 && sv1 <= fu) {
              timeObs1[i] = sv1; event1[i] = 1; dropEv1[i] = 0;
            } else if (dr1 <= sv1 && dr1 <= fu) {
              timeObs1[i] = dr1; event1[i] = 0; dropEv1[i] = 1;
            } else {
              timeObs1[i] = fu; event1[i] = 0; dropEv1[i] = 0;
            }
            if (sv2 <= dr2 && sv2 <= fu) {
              timeObs2[i] = sv2; event2[i] = 1; dropEv2[i] = 0;
            } else if (dr2 <= sv2 && dr2 <= fu) {
              timeObs2[i] = dr2; event2[i] = 0; dropEv2[i] = 1;
            } else {
              timeObs2[i] = fu; event2[i] = 0; dropEv2[i] = 0;
            }
          } else {
            if (sv1 <= dr1) {
              timeObs1[i] = sv1; event1[i] = 1; dropEv1[i] = 0;
            } else {
              timeObs1[i] = dr1; event1[i] = 0; dropEv1[i] = 1;
            }
            if (sv2 <= dr2) {
              timeObs2[i] = sv2; event2[i] = 1; dropEv2[i] = 0;
            } else {
              timeObs2[i] = dr2; event2[i] = 0; dropEv2[i] = 1;
            }
          }
          totalT1[i] = arrivalT[i] + timeObs1[i];
          totalT2[i] = arrivalT[i] + timeObs2[i];
        } // end cohort

        // determine analysis times & stages
        size_t nstages = K;
        bool ev1NotAch = false;
        bool ev2NotAch = false;

        if (useEvents) {
          totalte1.clear(); totalte2.clear();
          int nevents1 = 0, nevents2 = 0;
          for (size_t i = 0; i < N; ++i) {
            if (event1[i]) { ++nevents1; totalte1.push_back(totalT1[i]); }
            if (event2[i]) { ++nevents2; totalte2.push_back(totalT2[i]); }
          }
          if (nevents1 == 0 || nevents2 == 0) {
            thread_utils::push_thread_warning(
              std::string("No events for iteration ") + std::to_string(iter + 1) +
                " skipping this iteration.");
            out.summaryRows.clear();
            out.rawRows.clear();
            continue;
          }
          std::sort(totalte1.begin(), totalte1.end());
          std::sort(totalte2.begin(), totalte2.end());

          // PFS looks
          analysisT1.clear();
          size_t j1 = 0;
          if (Kpfs > 0) {
            for (j1 = 0; j1 < Kpfs; ++j1) {
              if (plannedEvents[j1] >= nevents1) break;
            }

            if (j1 == Kpfs) { // total number of PFS events exceeds planned
              for (size_t k = 0; k < Kpfs; ++k) {
                analysisT1.push_back(totalte1[plannedEvents[k] - 1] + 1e-12);
              }
            } else {
              for (size_t k = 0; k < j1; ++k) {
                analysisT1.push_back(totalte1[plannedEvents[k] - 1] + 1e-12);
              }
              analysisT1.push_back(totalte1.back() + 1e-12);
            }
          }

          // OS looks -> compute analysisT2
          analysisT2.clear();
          size_t j2 = 0;
          if (K > Kpfs) {
            for (j2 = 0; j2 < (K - Kpfs); ++j2) {
              if (plannedEvents[Kpfs + j2] >= nevents2) break;
            }

            if (j2 == (K - Kpfs)) { // total number of OS events exceeds planned
              for (size_t k = 0; k < (K - Kpfs); ++k) {
                analysisT2.push_back(totalte2[plannedEvents[Kpfs + k] - 1] + 1e-12);
              }
            } else {
              for (size_t k = 0; k < j2; ++k) {
                analysisT2.push_back(totalte2[plannedEvents[Kpfs + k] - 1] + 1e-12);
              }
              analysisT2.push_back(totalte2.back() + 1e-12);
            }
          }

          // combine PFS and OS looks to determine nstages and analysisTime array
          if (Kpfs == 0) { // only OS looks
            nstages = analysisT2.size();
            std::copy_n(analysisT2.begin(), nstages, analysisT.begin());
          } else if (K == Kpfs) { // only PFS looks
            nstages = analysisT1.size();
            std::copy_n(analysisT1.begin(), nstages, analysisT.begin());
          } else { // mixed
            if (analysisT2.back() > analysisT1.back()) {
              // OS looks after last PFS look contribute.
              // NOTE: In this case, the observed number of PFS events must exceed
              // the planned number of PFS events at look Kpfs, because otherwise
              // the last PFS event would be observed at analysisT1.back().
              // However, since the last OS event occurred on or after
              // analysisT2.back() > analysisT1.back(), this is a
              // contradiction as death is part of PFS event definition.
              // It follows that analysisT1.size() == Kpfs in this case.

              // find first OS look after last PFS look
              size_t l = 0;
              for (size_t idx = 0; idx < analysisT2.size(); ++idx) {
                if (analysisT2[idx] > analysisT1.back()) { l = idx; break; }
              }
              // number of stages
              nstages = Kpfs + (analysisT2.size() - l);
              // copy PFS looks unchanged and append relevant OS looks
              // keep PFS looks [0 .. Kpfs-1], then OS looks from l onwards,
              // which are the ones after last PFS look mapped to Kpfs + l onwards
              std::copy_n(analysisT1.begin(), Kpfs, analysisT.begin());
              size_t count = analysisT2.size() - l;
              std::copy_n(analysisT2.begin() + l, count, analysisT.begin() + Kpfs);
            } else {
              // only PFS looks matter
              nstages = analysisT1.size();
              std::copy_n(analysisT1.begin(), nstages, analysisT.begin());
            }
          }

          // evNotAch: check PFS and OS targetse
          if (Kpfs > 0 && nevents1 < plannedEvents[Kpfs - 1]) ev1NotAch = true;
          if (Kpfs < K && nevents2 < plannedEvents[K - 1]) ev2NotAch = true;
        } else { // calendar time
          std::copy_n(plannedTime.begin(), K, analysisT.begin());
        }

        // per-stage calculations
        for (size_t k = 0; k < nstages; ++k) {
          double time = analysisT[k];

          std::fill(n1x.begin(), n1x.end(), 0);
          std::fill(n2x.begin(), n2x.end(), 0);

          int events1e1 = 0, events2e1 = 0, dropouts1e1 = 0, dropouts2e1 = 0;
          int events1e2 = 0, events2e2 = 0, dropouts1e2 = 0, dropouts2e2 = 0;

          // censoring & counts
          for (size_t i = 0; i < N; ++i) {
            double ar = arrivalT[i];
            double sv1 = survivalT1[i], sv2 = survivalT2[i];
            double dr1 = dropoutT1[i], dr2 = dropoutT2[i];

            if (ar > time) {
              timeObs1[i] = time - ar; event1[i] = 0; dropEv1[i] = 0;
              timeObs2[i] = time - ar; event2[i] = 0; dropEv2[i] = 0;
              continue;
            }

            // endpoint 1 censoring
            if (fixedFollowup) {
              if (ar + sv1 <= time && sv1 <= dr1 && sv1 <= fu) {
                timeObs1[i] = sv1; event1[i] = 1; dropEv1[i] = 0;
              } else if (ar + dr1 <= time && dr1 <= sv1 && dr1 <= fu) {
                timeObs1[i] = dr1; event1[i] = 0; dropEv1[i] = 1;
              } else if (ar + fu <= time && fu <= sv1 && fu <= dr1) {
                timeObs1[i] = fu; event1[i] = 0; dropEv1[i] = 0;
              } else {
                timeObs1[i] = time - ar; event1[i] = 0; dropEv1[i] = 0;
              }
            } else {
              if (ar + sv1 <= time && sv1 <= dr1) {
                timeObs1[i] = sv1; event1[i] = 1; dropEv1[i] = 0;
              } else if (ar + dr1 <= time && dr1 <= sv1) {
                timeObs1[i] = dr1; event1[i] = 0; dropEv1[i] = 1;
              } else {
                timeObs1[i] = time - ar; event1[i] = 0; dropEv1[i] = 0;
              }
            }

            // endpoint2 censoring
            if (fixedFollowup) {
              if (ar + sv2 <= time && sv2 <= dr2 && sv2 <= fu) {
                timeObs2[i] = sv2; event2[i] = 1; dropEv2[i] = 0;
              } else if (ar + dr2 <= time && dr2 <= sv2 && dr2 <= fu) {
                timeObs2[i] = dr2; event2[i] = 0; dropEv2[i] = 1;
              } else if (ar + fu <= time && fu <= sv2 && fu <= dr2) {
                timeObs2[i] = fu; event2[i] = 0; dropEv2[i] = 0;
              } else {
                timeObs2[i] = time - ar; event2[i] = 0; dropEv2[i] = 0;
              }
            } else {
              if (ar + sv2 <= time && sv2 <= dr2) {
                timeObs2[i] = sv2; event2[i] = 1; dropEv2[i] = 0;
              } else if (ar + dr2 <= time && dr2 <= sv2) {
                timeObs2[i] = dr2; event2[i] = 0; dropEv2[i] = 1;
              } else {
                timeObs2[i] = time - ar; event2[i] = 0; dropEv2[i] = 0;
              }
            }

            size_t h = static_cast<size_t>(stratum[i] - 1);
            if (trtGrp[i] == 1) { ++n1x[h];
              if (event1[i]) ++events1e1; else if (dropEv1[i]) ++dropouts1e1;
              if (event2[i]) ++events1e2; else if (dropEv2[i]) ++dropouts1e2;
            } else { ++n2x[h];
              if (event1[i]) ++events2e1; else if (dropEv1[i]) ++dropouts2e1;
              if (event2[i]) ++events2e2; else if (dropEv2[i]) ++dropouts2e2;
            }
          }

          int accruals1 = std::accumulate(n1x.begin(), n1x.end(), 0);
          int accruals2 = std::accumulate(n2x.begin(), n2x.end(), 0);
          int totAccruals = accruals1 + accruals2;

          int totEventse1 = events1e1 + events2e1;
          int totDropoutse1 = dropouts1e1 + dropouts2e1;
          int totEventse2 = events1e2 + events2e2;
          int totDropoutse2 = dropouts1e2 + dropouts2e2;

          // append raw rows
          if (iter < maxRawIters) {
            for (size_t i = 0; i < N; ++i) {
              // skip subjects who haven't been enrolled by analysis time
              if (arrivalT[i] > time) continue;
              RawDatasetRow rr1;
              rr1.iterNum = static_cast<int>(iter + 1);
              rr1.stageNum = static_cast<int>(k + 1);
              rr1.analysisT = time;
              rr1.subjectId = static_cast<int>(i + 1);
              rr1.arrivalT = arrivalT[i];
              rr1.stratum = stratum[i];
              rr1.trtGrp = trtGrp[i];
              rr1.endpt = 1;
              rr1.survivalT = survivalT1[i];
              rr1.dropoutT = dropoutT1[i];
              rr1.timeObs = timeObs1[i];
              rr1.event = event1[i];
              rr1.dropEv = dropEv1[i];
              out.rawRows.push_back(std::move(rr1));

              RawDatasetRow rr2;
              rr2.iterNum = static_cast<int>(iter + 1);
              rr2.stageNum = static_cast<int>(k + 1);
              rr2.analysisT = time;
              rr2.subjectId = static_cast<int>(i + 1);
              rr2.arrivalT = arrivalT[i];
              rr2.stratum = stratum[i];
              rr2.trtGrp = trtGrp[i];
              rr2.endpt = 2;
              rr2.survivalT = survivalT2[i];
              rr2.dropoutT = dropoutT2[i];
              rr2.timeObs = timeObs2[i];
              rr2.event = event2[i];
              rr2.dropEv = dropEv2[i];
              out.rawRows.push_back(std::move(rr2));
            }
          }

          // compute stratified log-rank for endpoints
          for (int endpt = 1; endpt <= 2; ++endpt) {
            double hazardRatioH0;
            sub.clear();
            if (endpt == 1) {
              hazardRatioH0 = hazardRatioH0pfs;
              for (size_t i = 0; i < N; ++i) {
                if (timeObs1[i] > 0.0) sub.push_back(i);
              }
              std::sort(sub.begin(), sub.end(), [&](size_t a, size_t b) {
                return timeObs1[a] < timeObs1[b];
              });
            } else {
              hazardRatioH0 = hazardRatioH0os;
              for (size_t i = 0; i < N; ++i) {
                if (timeObs2[i] > 0.0) sub.push_back(i);
              }
              std::sort(sub.begin(), sub.end(), [&](size_t a, size_t b) {
                return timeObs2[a] < timeObs2[b];
              });
            }

            n1 = n1x; n2 = n2x;
            std::fill(km.begin(), km.end(), 1.0);
            double us = 0.0, vs = 0.0;

            for (size_t i = 0; i < sub.size(); ++i) {
              size_t idx = sub[i];
              size_t h = static_cast<size_t>(stratum[idx] - 1);

              double n1h = static_cast<double>(n1[h]);
              double n2h = static_cast<double>(n2[h]);
              double n1a = n1h * hazardRatioH0;
              double nt = n1h + n2h;
              double nta = n1a + n2h;

              bool evt = (endpt == 1 ? event1[idx] : event2[idx]);
              if (evt) {
                double wh = 1.0;
                if (rho1 != 0.0 || rho2 != 0.0) {
                  wh = std::pow(km[h], rho1) * std::pow(1.0 - km[h], rho2);
                  km[h] *= (1.0 - 1.0 / nt);
                }
                double treated = (trtGrp[idx] == 1 ? 1.0 : 0.0);
                us += wh * (treated - n1a / nta);
                vs += wh * wh * n1a * n2h / (nta * nta);
              }

              if (trtGrp[idx] == 1) --n1[h]; else --n2[h];
            } // log-rank

            double z = (vs > 0.0 ? us / std::sqrt(vs) : 0.0);

            StageSummaryRow sr;
            sr.iterNum = static_cast<int>(iter + 1);
            sr.evNotAch1 = ev1NotAch ? 1 : 0;
            sr.evNotAch2 = ev2NotAch ? 1 : 0;
            sr.stageNum = static_cast<int>(k + 1);
            sr.analysisT = time;
            sr.accruals1 = accruals1;
            sr.accruals2 = accruals2;
            sr.totAccruals = totAccruals;
            sr.endpt = endpt;
            if (endpt == 1) {
              sr.events1 = events1e1;
              sr.events2 = events2e1;
              sr.totEvents = totEventse1;
              sr.dropouts1 = dropouts1e1;
              sr.dropouts2 = dropouts2e1;
              sr.totDropouts = totDropoutse1;
            } else {
              sr.events1 = events1e2;
              sr.events2 = events2e2;
              sr.totEvents = totEventse2;
              sr.dropouts1 = dropouts1e2;
              sr.dropouts2 = dropouts2e2;
              sr.totDropouts = totDropoutse2;
            }
            sr.uscore = us; sr.vscore = vs; sr.logRank = z;
            out.summaryRows.push_back(std::move(sr));
          } // endpoints loop
        } // per-stage
      } // iter
    } // operator()
  }; // SimWorker

  // construct and run worker
  SimWorker worker(
      K, Kpfs, hazardRatioH0pfs, hazardRatioH0os,
      allocation1, allocation2,
      accrualTime, accrualIntensity, tau, stratumFraction,
      rho_pd_os,
      lambda1pfsx, lambda2pfsx, lambda1osx, lambda2osx,
      gamma1pfsx, gamma2pfsx, gamma1osx, gamma2osx,
      tau1pdx, tau2pdx, lambda1pd, lambda2pd, gamma1pd, gamma2pd,
      N, fu, fixedFollowup, rho1, rho2,
      plannedEvents, plannedTime,
      maxIters, maxRawIters, seeds, useEvents, nstrata, rho_pd_os_pyth_comp,
      &results
  );

  RcppParallel::parallelFor(0, maxIters, worker);

  // Flatten results into final containers
  size_t nsr = 0, nrr = 0;
  for (size_t iter = 0; iter < maxIters; ++iter) {
    nsr += results[iter].summaryRows.size();
    nrr += results[iter].rawRows.size();
  }
  if (nsr == 0) throw std::runtime_error(
    "No iterations with observed events. Unable to produce output.");

  // prepare final containers
  std::vector<int> sum_iterNum; sum_iterNum.reserve(nsr);
  std::vector<unsigned char> sum_ev1NotAch; sum_ev1NotAch.reserve(nsr);
  std::vector<unsigned char> sum_ev2NotAch; sum_ev2NotAch.reserve(nsr);
  std::vector<int> sum_stageNum; sum_stageNum.reserve(nsr);
  std::vector<double> sum_analysisT; sum_analysisT.reserve(nsr);
  std::vector<int> sum_accruals1; sum_accruals1.reserve(nsr);
  std::vector<int> sum_accruals2; sum_accruals2.reserve(nsr);
  std::vector<int> sum_totAccruals; sum_totAccruals.reserve(nsr);
  std::vector<std::string> sum_endpt; sum_endpt.reserve(nsr);
  std::vector<int> sum_events1; sum_events1.reserve(nsr);
  std::vector<int> sum_events2; sum_events2.reserve(nsr);
  std::vector<int> sum_totEvents; sum_totEvents.reserve(nsr);
  std::vector<int> sum_dropouts1; sum_dropouts1.reserve(nsr);
  std::vector<int> sum_dropouts2; sum_dropouts2.reserve(nsr);
  std::vector<int> sum_totDropouts; sum_totDropouts.reserve(nsr);
  std::vector<double> sum_uscore; sum_uscore.reserve(nsr);
  std::vector<double> sum_vscore; sum_vscore.reserve(nsr);
  std::vector<double> sum_logRank; sum_logRank.reserve(nsr);

  // raw final containers
  std::vector<int> raw_iterNum; raw_iterNum.reserve(nrr);
  std::vector<int> raw_stageNum; raw_stageNum.reserve(nrr);
  std::vector<double> raw_analysisT; raw_analysisT.reserve(nrr);
  std::vector<int> raw_subjectId; raw_subjectId.reserve(nrr);
  std::vector<double> raw_arrivalT; raw_arrivalT.reserve(nrr);
  std::vector<int> raw_stratum; raw_stratum.reserve(nrr);
  std::vector<int> raw_trtGrp; raw_trtGrp.reserve(nrr);
  std::vector<std::string> raw_endpt; raw_endpt.reserve(nrr);
  std::vector<double> raw_survivalT; raw_survivalT.reserve(nrr);
  std::vector<double> raw_dropoutT; raw_dropoutT.reserve(nrr);
  std::vector<double> raw_timeObs; raw_timeObs.reserve(nrr);
  std::vector<unsigned char> raw_event; raw_event.reserve(nrr);
  std::vector<unsigned char> raw_dropEv; raw_dropEv.reserve(nrr);

  for (size_t iter = 0; iter < maxIters; ++iter) {
    const auto& srows = results[iter].summaryRows;
    for (const auto& r : srows) {
      sum_iterNum.push_back(r.iterNum);
      sum_ev1NotAch.push_back(r.evNotAch1);
      sum_ev2NotAch.push_back(r.evNotAch2);
      sum_stageNum.push_back(r.stageNum);
      sum_analysisT.push_back(r.analysisT);
      sum_accruals1.push_back(r.accruals1);
      sum_accruals2.push_back(r.accruals2);
      sum_totAccruals.push_back(r.totAccruals);
      sum_endpt.push_back(r.endpt == 1 ? "PFS" : "OS");
      sum_events1.push_back(r.events1);
      sum_events2.push_back(r.events2);
      sum_totEvents.push_back(r.totEvents);
      sum_dropouts1.push_back(r.dropouts1);
      sum_dropouts2.push_back(r.dropouts2);
      sum_totDropouts.push_back(r.totDropouts);
      sum_uscore.push_back(r.uscore);
      sum_vscore.push_back(r.vscore);
      sum_logRank.push_back(r.logRank);
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
        raw_endpt.push_back(rr.endpt == 1 ? "PFS" : "OS");
        raw_survivalT.push_back(rr.survivalT);
        raw_dropoutT.push_back(rr.dropoutT);
        raw_timeObs.push_back(rr.timeObs);
        raw_event.push_back(rr.event);
        raw_dropEv.push_back(rr.dropEv);
      }
    }
  }

  // Build DataFrameCpp summary
  DataFrameCpp sumdata;
  sumdata.push_back(std::move(sum_iterNum), "iterationNumber");
  sumdata.push_back(std::move(sum_ev1NotAch), "events1NotAchieved");
  sumdata.push_back(std::move(sum_ev2NotAch), "events2NotAchieved");
  sumdata.push_back(std::move(sum_stageNum), "stageNumber");
  sumdata.push_back(std::move(sum_analysisT), "analysisTime");
  sumdata.push_back(std::move(sum_accruals1), "accruals1");
  sumdata.push_back(std::move(sum_accruals2), "accruals2");
  sumdata.push_back(std::move(sum_totAccruals), "totalAccruals");
  sumdata.push_back(std::move(sum_endpt), "endpoint");
  sumdata.push_back(std::move(sum_events1), "events1");
  sumdata.push_back(std::move(sum_events2), "events2");
  sumdata.push_back(std::move(sum_totEvents), "totalEvents");
  sumdata.push_back(std::move(sum_dropouts1), "dropouts1");
  sumdata.push_back(std::move(sum_dropouts2), "dropouts2");
  sumdata.push_back(std::move(sum_totDropouts), "totalDropouts");
  sumdata.push_back(std::move(sum_uscore), "uscore");
  sumdata.push_back(std::move(sum_vscore), "vscore");
  sumdata.push_back(std::move(sum_logRank), "logRankStatistic");

  ListCpp result;
  result.push_back(sumdata, "sumdata");

  if (!raw_iterNum.empty()) {
    DataFrameCpp rawdata;
    rawdata.push_back(std::move(raw_iterNum), "iterationNumber");
    rawdata.push_back(std::move(raw_stageNum), "stageNumber");
    rawdata.push_back(std::move(raw_analysisT), "analysisTime");
    rawdata.push_back(std::move(raw_subjectId), "subjectId");
    rawdata.push_back(std::move(raw_arrivalT), "arrivalTime");
    rawdata.push_back(std::move(raw_stratum), "stratum");
    rawdata.push_back(std::move(raw_trtGrp), "treatmentGroup");
    rawdata.push_back(std::move(raw_endpt), "endpoint");
    rawdata.push_back(std::move(raw_survivalT), "survivalTime");
    rawdata.push_back(std::move(raw_dropoutT), "dropoutTime");
    rawdata.push_back(std::move(raw_timeObs), "timeUnderObservation");
    rawdata.push_back(std::move(raw_event), "event");
    rawdata.push_back(std::move(raw_dropEv), "dropoutEvent");
    result.push_back(rawdata, "rawdata");
  }

  return result;
}



//' @title Log-Rank Test Simulation for PFS and OS Endpoints
//' @description Performs simulation for two-endpoint (PFS and OS) two-arm
//' group sequential trials based on weighted log-rank test. The first
//' \code{kMaxpfs} looks are driven by the total number of PFS events in
//' two arms combined, and the subsequent looks are driven by the total
//' number of OS events in two arms combined. Alternatively,
//' the analyses can be planned to occur at specified calendar times.
//'
//' @inheritParams param_kMax
//' @param kMaxpfs Number of stages with timing determined by PFS events.
//'   Ranges from 0 (none) to \code{kMax}.
//' @param hazardRatioH0pfs Hazard ratio under the null hypothesis for the
//'   active treatment vs control for PFS. Defaults to 1 for
//'   superiority test.
//' @param hazardRatioH0os Hazard ratio under the null hypothesis for the
//'   active treatment vs control for OS. Defaults to 1 for
//'   superiority test.
//' @param allocation1 Number of subjects in the treatment group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in the control group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param rho_pd_os The correlation coefficient for the standard
//'   bivariate normal random variables used to generate time to disease
//'   progression (PD) and time to death using the inverse CDF method.
//' @param lambda1pfs A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the treatment group and PFS.
//' @param lambda2pfs A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the control group and PFS.
//' @param lambda1os A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the treatment group and OS.
//' @param lambda2os A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the control group and OS.
//' @param gamma1pfs The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the treatment group and PFS.
//' @param gamma2pfs The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the control group and PFS.
//' @param gamma1os The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the treatment group and OS.
//' @param gamma2os The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the control group and OS.
//' @param n Sample size.
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @param plannedEvents The planned cumulative total number of PFS events at
//'   Look 1 to Look \code{kMaxpfs} and the planned cumulative total number
//'   of OS events at Look \code{kMaxpfs+1} to Look \code{kMax}.
//' @param plannedTime The calendar times for the analyses. To use calendar
//'   time to plan the analyses, \code{plannedEvents} should be missing.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
//'   stage to extract.
//' @param seed The seed to reproduce the simulation results.
//'
//' @return A list with 2 components:
//'
//' * \code{sumdata}: A data frame of summary data by iteration and stage:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{eventsNotAchieved}: Whether the target number of events
//'       is not achieved for the iteration.
//'
//'     - \code{stageNumber}: The stage number, covering all stages even if
//'       the trial stops at an interim look.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{accruals1}: The number of subjects enrolled at the stage for
//'       the treatment group.
//'
//'     - \code{accruals2}: The number of subjects enrolled at the stage for
//'       the control group.
//'
//'     - \code{totalAccruals}: The total number of subjects enrolled at
//'       the stage.
//'
//'     - \code{endpoint}: The endpoint (1 for PFS or 2 for OS) under
//'       consideration.
//'
//'     - \code{events1}: The number of events at the stage for
//'       the treatment group.
//'
//'     - \code{events2}: The number of events at the stage for
//'       the control group.
//'
//'     - \code{totalEvents}: The total number of events at the stage.
//'
//'     - \code{dropouts1}: The number of dropouts at the stage for
//'       the treatment group.
//'
//'     - \code{dropouts2}: The number of dropouts at the stage for
//'       the control group.
//'
//'     - \code{totalDropouts}: The total number of dropouts at the stage.
//'
//'     - \code{uscore}: The numerator of the log-rank test statistic for
//'       the endpoint.
//'
//'     - \code{vscore}: The variance of the log-rank test statistic for
//'       the endpoint.
//'
//'     - \code{logRankStatistic}: The log-rank test Z-statistic for
//'       the endpoint.
//'
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//'   positive integer): A data frame for subject-level data for selected
//'   replications, containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage under consideration.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{subjectId}: The subject ID.
//'
//'     - \code{arrivalTime}: The enrollment time for the subject.
//'
//'     - \code{stratum}: The stratum for the subject.
//'
//'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
//'       subject.
//'
//'     - \code{endpoint}: The endpoint (1 for PFS or 2 for OS) under
//'       consideration for the row. Each subject will have two rows, one
//'       for each endpoint.
//'
//'     - \code{survivalTime}: The underlying survival time for the
//'       event endpoint for the subject.
//'
//'     - \code{dropoutTime}: The underlying dropout time for the
//'       event endpoint for the subject.
//'
//'     - \code{timeUnderObservation}: The time under observation
//'       since randomization for the event endpoint for the subject.
//'
//'     - \code{event}: Whether the subject experienced the event
//'       endpoint.
//'
//'     - \code{dropoutEvent}: Whether the subject dropped out for the
//'       endpoint.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' sim1 = lrsim2e(
//'   kMax = 3,
//'   kMaxpfs = 2,
//'   allocation1 = 2,
//'   allocation2 = 1,
//'   accrualTime = c(0, 8),
//'   accrualIntensity = c(10, 28),
//'   piecewiseSurvivalTime = 0,
//'   rho_pd_os = 0,
//'   lambda1pfs = log(2)/12*0.60,
//'   lambda2pfs = log(2)/12,
//'   lambda1os = log(2)/30*0.65,
//'   lambda2os = log(2)/30,
//'   n = 420,
//'   plannedEvents = c(186, 259, 183),
//'   maxNumberOfIterations = 1000,
//'   maxNumberOfRawDatasetsPerStage = 1,
//'   seed = 314159)
//'
//' head(sim1$sumdata)
//' head(sim1$rawdata)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List lrsim2e(
    const int kMax = 1,
    const int kMaxpfs = 1,
    const double hazardRatioH0pfs = 1,
    const double hazardRatioH0os = 1,
    const int allocation1 = 1,
    const int allocation2 = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const double rho_pd_os = 0,
    const Rcpp::NumericVector& lambda1pfs = NA_REAL,
    const Rcpp::NumericVector& lambda2pfs = NA_REAL,
    const Rcpp::NumericVector& lambda1os = NA_REAL,
    const Rcpp::NumericVector& lambda2os = NA_REAL,
    const Rcpp::NumericVector& gamma1pfs = 0,
    const Rcpp::NumericVector& gamma2pfs = 0,
    const Rcpp::NumericVector& gamma1os = 0,
    const Rcpp::NumericVector& gamma2os = 0,
    const int n = NA_INTEGER,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = 0,
    const double rho1 = 0,
    const double rho2 = 0,
    const Rcpp::IntegerVector& plannedEvents = NA_INTEGER,
    const Rcpp::NumericVector& plannedTime = NA_REAL,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = 0) {

  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1pfs = Rcpp::as<std::vector<double>>(lambda1pfs);
  auto lam2pfs = Rcpp::as<std::vector<double>>(lambda2pfs);
  auto lam1os = Rcpp::as<std::vector<double>>(lambda1os);
  auto lam2os = Rcpp::as<std::vector<double>>(lambda2os);
  auto gam1pfs = Rcpp::as<std::vector<double>>(gamma1pfs);
  auto gam2pfs = Rcpp::as<std::vector<double>>(gamma2os);
  auto gam1os = Rcpp::as<std::vector<double>>(gamma1os);
  auto gam2os = Rcpp::as<std::vector<double>>(gamma2os);
  auto plannedE = Rcpp::as<std::vector<int>>(plannedEvents);
  auto plannedT = Rcpp::as<std::vector<double>>(plannedTime);

  auto out = lrsim2ecpp(
    kMax, kMaxpfs, hazardRatioH0pfs, hazardRatioH0os,
    allocation1, allocation2, accrualT, accrualInt,
    pwSurvT, stratumFrac, rho_pd_os,
    lam1pfs, lam2pfs, lam1os, lam2os,
    gam1pfs, gam2pfs, gam1os, gam2os,
    n, followupTime, fixedFollowup, rho1, rho2,
    plannedE, plannedT, maxNumberOfIterations,
    maxNumberOfRawDatasetsPerStage, seed);

  thread_utils::drain_thread_warnings_to_R();

  return Rcpp::wrap(out);
}


ListCpp lrsim2e3acpp(
    const int kMax,
    const int kMaxpfs,
    const double hazardRatioH013pfs,
    const double hazardRatioH023pfs,
    const double hazardRatioH012pfs,
    const double hazardRatioH013os,
    const double hazardRatioH023os,
    const double hazardRatioH012os,
    const int allocation1,
    const int allocation2,
    const int allocation3,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const double rho_pd_os,
    const std::vector<double>& lambda1pfs,
    const std::vector<double>& lambda2pfs,
    const std::vector<double>& lambda3pfs,
    const std::vector<double>& lambda1os,
    const std::vector<double>& lambda2os,
    const std::vector<double>& lambda3os,
    const std::vector<double>& gamma1pfs,
    const std::vector<double>& gamma2pfs,
    const std::vector<double>& gamma3pfs,
    const std::vector<double>& gamma1os,
    const std::vector<double>& gamma2os,
    const std::vector<double>& gamma3os,
    const int n,
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
  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
  size_t K = static_cast<size_t>(kMax);

  int kMaxpfsx = kMaxpfs;
  if (kMaxpfsx < 0) kMaxpfsx = kMax;
  if (kMaxpfsx > kMax)
    throw std::invalid_argument("kMaxpfs must be less than or equal to kMax");
  size_t Kpfs = static_cast<size_t>(kMaxpfsx);

  bool useEvents;
  if (none_na(plannedEvents)) {
    useEvents = true;
    if (plannedEvents.empty() || plannedEvents[0] <= 0)
      throw std::invalid_argument("plannedEvents must be positive");
    if (plannedEvents.size() != K)
      throw std::invalid_argument("Invalid length for plannedEvents");
    if (Kpfs > 1) {
      for (size_t i = 1; i < Kpfs; ++i) {
        if (plannedEvents[i] <= plannedEvents[i-1])
          throw std::invalid_argument("plannedEvents for PFS must be increasing");
      }
    }
    if (K - Kpfs > 1) {
      for (size_t i = Kpfs + 1; i < plannedEvents.size(); ++i) {
        if (plannedEvents[i] <= plannedEvents[i-1])
          throw std::invalid_argument("plannedEvents for OS must be increasing");
      }
    }
  } else if (none_na(plannedTime)) {
    useEvents = false;
    if (plannedTime[0] <= 0.0)
      throw std::invalid_argument("plannedTime must be positive");
    if (plannedTime.size() != K)
      throw std::invalid_argument("Invalid length for plannedTime");
    if (any_nonincreasing(plannedTime))
      throw std::invalid_argument("plannedTime must be increasing");
  } else {
    throw std::invalid_argument("Either plannedEvents or plannedTime must be given");
  }

  if (hazardRatioH013pfs <= 0.0 || hazardRatioH023pfs <= 0.0 ||
      hazardRatioH012pfs <= 0.0)
    throw std::invalid_argument("PFS hazard ratios under H0 must be positive");
  if (hazardRatioH013os <= 0.0 || hazardRatioH023os <= 0.0 ||
      hazardRatioH012os <= 0.0)
    throw std::invalid_argument("OS hazard ratios under H0 must be positive");
  if (allocation1 < 1 || allocation2 < 1 || allocation3 < 1)
    throw std::invalid_argument("allocations must be positive integers");
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
  if (rho_pd_os <= -1.0 || rho_pd_os >= 1.0)
    throw std::invalid_argument("rho_pd_os must lie in (-1, 1)");

  if (!none_na(lambda1pfs)) throw std::invalid_argument("lambda1pfs must be provided");
  if (!none_na(lambda2pfs)) throw std::invalid_argument("lambda2pfs must be provided");
  if (!none_na(lambda3pfs)) throw std::invalid_argument("lambda3pfs must be provided");
  if (!none_na(lambda1os)) throw std::invalid_argument("lambda1os must be provided");
  if (!none_na(lambda2os)) throw std::invalid_argument("lambda2os must be provided");
  if (!none_na(lambda3os)) throw std::invalid_argument("lambda3os must be provided");
  for (double v : lambda1pfs) {
    if (v < 0.0) throw std::invalid_argument("lambda1pfs must be non-negative");
  }
  for (double v : lambda2pfs) {
    if (v < 0.0) throw std::invalid_argument("lambda2pfs must be non-negative");
  }
  for (double v : lambda3pfs) {
    if (v < 0.0) throw std::invalid_argument("lambda3pfs must be non-negative");
  }
  for (double v : lambda1os) {
    if (v < 0.0) throw std::invalid_argument("lambda1os must be non-negative");
  }
  for (double v : lambda2os) {
    if (v < 0.0) throw std::invalid_argument("lambda2os must be non-negative");
  }
  for (double v : lambda3os) {
    if (v < 0.0) throw std::invalid_argument("lambda3os must be non-negative");
  }
  for (double v : gamma1pfs) {
    if (v < 0.0) throw std::invalid_argument("gamma1pfs must be non-negative");
  }
  for (double v : gamma2pfs) {
    if (v < 0.0) throw std::invalid_argument("gamma2pfs must be non-negative");
  }
  for (double v : gamma3pfs) {
    if (v < 0.0) throw std::invalid_argument("gamma3pfs must be non-negative");
  }
  for (double v : gamma1os) {
    if (v < 0.0) throw std::invalid_argument("gamma1os must be non-negative");
  }
  for (double v : gamma2os) {
    if (v < 0.0) throw std::invalid_argument("gamma2os must be non-negative");
  }
  for (double v : gamma3os) {
    if (v < 0.0) throw std::invalid_argument("gamma3os must be non-negative");
  }
  if (n == INT_MIN) throw std::invalid_argument("n must be provided");
  if (n <= 0) throw std::invalid_argument("n must be positive");
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

  size_t N = static_cast<size_t>(n);
  size_t maxIters = static_cast<size_t>(maxNumberOfIterations);
  size_t maxRawIters = static_cast<size_t>(maxNumberOfRawDatasetsPerStage);
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  size_t nintv2 = (nintv == 1 ? 10u : nintv + 10u);
  const std::vector<double>& tau = piecewiseSurvivalTime;
  const double fu = followupTime;
  const double rho_pd_os_pyth_comp = std::sqrt(1 - rho_pd_os * rho_pd_os);

  // expand stratified inputs (main thread)
  auto lambda1pfsx = expand_stratified(lambda1pfs, nstrata, nintv, "lambda1pfs");
  auto lambda2pfsx = expand_stratified(lambda2pfs, nstrata, nintv, "lambda2pfs");
  auto lambda3pfsx = expand_stratified(lambda3pfs, nstrata, nintv, "lambda3pfs");
  auto lambda1osx  = expand_stratified(lambda1os,  nstrata, nintv, "lambda1os");
  auto lambda2osx  = expand_stratified(lambda2os,  nstrata, nintv, "lambda2os");
  auto lambda3osx  = expand_stratified(lambda3os,  nstrata, nintv, "lambda3os");
  auto gamma1pfsx  = expand_stratified(gamma1pfs,  nstrata, nintv, "gamma1pfs");
  auto gamma2pfsx  = expand_stratified(gamma2pfs,  nstrata, nintv, "gamma2pfs");
  auto gamma3pfsx  = expand_stratified(gamma3pfs,  nstrata, nintv, "gamma3pfs");
  auto gamma1osx   = expand_stratified(gamma1os,   nstrata, nintv, "gamma1os");
  auto gamma2osx   = expand_stratified(gamma2os,   nstrata, nintv, "gamma2os");
  auto gamma3osx   = expand_stratified(gamma3os,   nstrata, nintv, "gamma3os");

  // compute pd hazards per stratum & arm (main thread)
  std::vector<std::vector<double>> tau1pdx(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> tau2pdx(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> tau3pdx(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> lambda1pd(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> lambda2pd(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> lambda3pd(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> gamma1pd(nstrata, std::vector<double>(nintv));
  std::vector<std::vector<double>> gamma2pd(nstrata, std::vector<double>(nintv));
  std::vector<std::vector<double>> gamma3pd(nstrata, std::vector<double>(nintv));

  for (size_t s = 0; s < nstrata; ++s) {
    ListCpp a1 = hazard_pdcpp(tau, lambda1pfsx[s], lambda1osx[s], rho_pd_os);
    ListCpp a2 = hazard_pdcpp(tau, lambda2pfsx[s], lambda2osx[s], rho_pd_os);
    ListCpp a3 = hazard_pdcpp(tau, lambda3pfsx[s], lambda3osx[s], rho_pd_os);

    tau1pdx[s] = a1.get<std::vector<double>>("piecewiseSurvivalTime");
    tau2pdx[s] = a2.get<std::vector<double>>("piecewiseSurvivalTime");
    tau3pdx[s] = a3.get<std::vector<double>>("piecewiseSurvivalTime");
    lambda1pd[s] = a1.get<std::vector<double>>("hazard_pd");
    lambda2pd[s] = a2.get<std::vector<double>>("hazard_pd");
    lambda3pd[s] = a3.get<std::vector<double>>("hazard_pd");

    for (size_t t = 0; t < nintv; ++t) {
      gamma1pd[s][t] = gamma1pfsx[s][t] - gamma1osx[s][t];
      gamma2pd[s][t] = gamma2pfsx[s][t] - gamma2osx[s][t];
      gamma3pd[s][t] = gamma3pfsx[s][t] - gamma3osx[s][t];
    }
  }

  // seeds for reproducibility
  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter) seeds[iter] = master_rng();


  // One summary (stage-level) row produced by an iteration
  struct StageSummaryRow {
    int iterNum = 0;
    unsigned char evNotAch1 = 0, evNotAch2 = 0;
    int stageNum = 0;
    double analysisT = 0.0;
    int accruals1 = 0, accruals2 = 0, accruals3 = 0, totAccruals = 0;
    int endpt = 0; // 1 == PFS, 2 == OS
    int events1 = 0, events2 = 0, events3 = 0, totEvents = 0;
    int dropouts1 = 0, dropouts2 = 0, dropouts3 = 0, totDropouts = 0;
    double uscore13 = 0.0, vscore13 = 0.0, logRank13 = 0.0;
    double uscore23 = 0.0, vscore23 = 0.0, logRank23 = 0.0;
    double uscore12 = 0.0, vscore12 = 0.0, logRank12 = 0.0;
  };

  // One subject-level (raw) row for a particular iteration and stage
  struct RawDatasetRow {
    int iterNum = 0, stageNum = 0;
    double analysisT = 0.0;
    int subjectId = 0;
    double arrivalT = 0.0;
    int stratum = 0, trtGrp = 0;
    int endpt = 0; // 1 == PFS, 2 == OS
    double survivalT = 0.0, dropoutT = 0.0, timeObs = 0.0;
    unsigned char event = 0, dropEv = 0;
  };

  // Per-iteration container written exclusively by worker thread
  struct IterationResult {
    std::vector<StageSummaryRow> summaryRows;
    std::vector<RawDatasetRow> rawRows;
    void reserveForSummary(size_t approxRows) { summaryRows.reserve(approxRows); }
    void reserveForRaw(size_t approxRows) { rawRows.reserve(approxRows); }
  };

  // pre-size per-iteration results
  std::vector<IterationResult> results;
  results.resize(maxIters);

  // Worker struct defined inside function
  struct SimWorker : public RcppParallel::Worker {
    // Inputs (const refs)
    const size_t K;
    const size_t Kpfs;
    const double hazardRatioH013pfs;
    const double hazardRatioH023pfs;
    const double hazardRatioH012pfs;
    const double hazardRatioH013os;
    const double hazardRatioH023os;
    const double hazardRatioH012os;
    const int allocation1;
    const int allocation2;
    const int allocation3;
    const std::vector<double>& accrualTime;
    const std::vector<double>& accrualIntensity;
    const std::vector<double>& tau;
    const std::vector<double>& stratumFraction;
    const double rho_pd_os;
    const std::vector<std::vector<double>>& lambda1pfsx;
    const std::vector<std::vector<double>>& lambda2pfsx;
    const std::vector<std::vector<double>>& lambda3pfsx;
    const std::vector<std::vector<double>>& lambda1osx;
    const std::vector<std::vector<double>>& lambda2osx;
    const std::vector<std::vector<double>>& lambda3osx;
    const std::vector<std::vector<double>>& gamma1pfsx;
    const std::vector<std::vector<double>>& gamma2pfsx;
    const std::vector<std::vector<double>>& gamma3pfsx;
    const std::vector<std::vector<double>>& gamma1osx;
    const std::vector<std::vector<double>>& gamma2osx;
    const std::vector<std::vector<double>>& gamma3osx;
    const std::vector<std::vector<double>>& tau1pdx;
    const std::vector<std::vector<double>>& tau2pdx;
    const std::vector<std::vector<double>>& tau3pdx;
    const std::vector<std::vector<double>>& lambda1pd;
    const std::vector<std::vector<double>>& lambda2pd;
    const std::vector<std::vector<double>>& lambda3pd;
    const std::vector<std::vector<double>>& gamma1pd;
    const std::vector<std::vector<double>>& gamma2pd;
    const std::vector<std::vector<double>>& gamma3pd;

    const size_t N;
    const double fu;
    const bool fixedFollowup;
    const double rho1;
    const double rho2;
    const std::vector<int>& plannedEvents;
    const std::vector<double>& plannedTime;
    const size_t maxIters;
    const size_t maxRawIters;
    const std::vector<uint64_t>& seeds;
    const bool useEvents;
    const size_t nstrata;
    const double rho_pd_os_pyth_comp;

    // Output pointer
    std::vector<IterationResult>* results;

    SimWorker(
      size_t K_,
      size_t Kpfs_,
      double hazardRatioH013pfs_,
      double hazardRatioH023pfs_,
      double hazardRatioH012pfs_,
      double hazardRatioH013os_,
      double hazardRatioH023os_,
      double hazardRatioH012os_,
      int allocation1_,
      int allocation2_,
      int allocation3_,
      const std::vector<double>& accrualTime_,
      const std::vector<double>& accrualIntensity_,
      const std::vector<double>& tau_,
      const std::vector<double>& stratumFraction_,
      double rho_pd_os_,
      const std::vector<std::vector<double>>& lambda1pfsx_,
      const std::vector<std::vector<double>>& lambda2pfsx_,
      const std::vector<std::vector<double>>& lambda3pfsx_,
      const std::vector<std::vector<double>>& lambda1osx_,
      const std::vector<std::vector<double>>& lambda2osx_,
      const std::vector<std::vector<double>>& lambda3osx_,
      const std::vector<std::vector<double>>& gamma1pfsx_,
      const std::vector<std::vector<double>>& gamma2pfsx_,
      const std::vector<std::vector<double>>& gamma3pfsx_,
      const std::vector<std::vector<double>>& gamma1osx_,
      const std::vector<std::vector<double>>& gamma2osx_,
      const std::vector<std::vector<double>>& gamma3osx_,
      const std::vector<std::vector<double>>& tau1pdx_,
      const std::vector<std::vector<double>>& tau2pdx_,
      const std::vector<std::vector<double>>& tau3pdx_,
      const std::vector<std::vector<double>>& lambda1pd_,
      const std::vector<std::vector<double>>& lambda2pd_,
      const std::vector<std::vector<double>>& lambda3pd_,
      const std::vector<std::vector<double>>& gamma1pd_,
      const std::vector<std::vector<double>>& gamma2pd_,
      const std::vector<std::vector<double>>& gamma3pd_,
      size_t N_,
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
      size_t nstrata_,
      double rho_pd_os_pyth_comp_,
      std::vector<IterationResult>* results_)
      : K(K_),
        Kpfs(Kpfs_),
        hazardRatioH013pfs(hazardRatioH013pfs_),
        hazardRatioH023pfs(hazardRatioH023pfs_),
        hazardRatioH012pfs(hazardRatioH012pfs_),
        hazardRatioH013os(hazardRatioH013os_),
        hazardRatioH023os(hazardRatioH023os_),
        hazardRatioH012os(hazardRatioH012os_),
        allocation1(allocation1_),
        allocation2(allocation2_),
        allocation3(allocation3_),
        accrualTime(accrualTime_),
        accrualIntensity(accrualIntensity_),
        tau(tau_),
        stratumFraction(stratumFraction_),
        rho_pd_os(rho_pd_os_),
        lambda1pfsx(lambda1pfsx_),
        lambda2pfsx(lambda2pfsx_),
        lambda3pfsx(lambda3pfsx_),
        lambda1osx(lambda1osx_),
        lambda2osx(lambda2osx_),
        lambda3osx(lambda3osx_),
        gamma1pfsx(gamma1pfsx_),
        gamma2pfsx(gamma2pfsx_),
        gamma3pfsx(gamma3pfsx_),
        gamma1osx(gamma1osx_),
        gamma2osx(gamma2osx_),
        gamma3osx(gamma3osx_),
        tau1pdx(tau1pdx_),
        tau2pdx(tau2pdx_),
        tau3pdx(tau3pdx_),
        lambda1pd(lambda1pd_),
        lambda2pd(lambda2pd_),
        lambda3pd(lambda3pd_),
        gamma1pd(gamma1pd_),
        gamma2pd(gamma2pd_),
        gamma3pd(gamma3pd_),
        N(N_),
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
        nstrata(nstrata_),
        rho_pd_os_pyth_comp(rho_pd_os_pyth_comp_),
        results(results_)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      // Local per-worker buffers
      std::vector<int> stratum(N), trtGrp(N);
      std::vector<double> arrivalT(N), survivalT1(N), survivalT2(N);
      std::vector<double> dropoutT1(N), dropoutT2(N);
      std::vector<double> timeObs1(N), timeObs2(N);
      std::vector<double> totalT1(N), totalT2(N);
      std::vector<unsigned char> event1(N), event2(N);
      std::vector<unsigned char> dropEv1(N), dropEv2(N);

      std::vector<int> b1(nstrata), b2(nstrata), b3(nstrata);
      std::vector<int> n1(nstrata), n2(nstrata), n3(nstrata);
      std::vector<int> n1x(nstrata), n2x(nstrata), n3x(nstrata);
      std::vector<double> km13(nstrata), km23(nstrata), km12(nstrata);
      std::vector<double> cumF(nstrata);
      std::partial_sum(stratumFraction.begin(), stratumFraction.end(), cumF.begin());

      std::vector<double> analysisT(K);
      std::vector<double> analysisT1; analysisT1.reserve(Kpfs);
      std::vector<double> analysisT2; analysisT2.reserve(K - Kpfs);
      std::vector<double> totalte1; totalte1.reserve(N);
      std::vector<double> totalte2; totalte2.reserve(N);
      std::vector<size_t> sub; sub.reserve(N);

      for (size_t iter = begin; iter < end; ++iter) {
        // deterministic per-iteration RNG
        boost::random::mt19937_64 rng_local(seeds[iter]);
        boost::random::uniform_real_distribution<double> unif(0.0, 1.0);
        boost::random::normal_distribution<double> norm(0.0, 1.0);

        IterationResult& out = (*results)[iter];
        out.summaryRows.clear();
        out.rawRows.clear();
        if (iter < maxRawIters) out.reserveForRaw(K * N);
        out.reserveForSummary(K * 2);

        // reset blocks
        std::fill(b1.begin(), b1.end(), allocation1);
        std::fill(b2.begin(), b2.end(), allocation2);
        std::fill(b3.begin(), b3.end(), allocation3);

        double enrollt = 0.0;

        // generate cohort
        for (size_t i = 0; i < N; ++i) {
          double u = unif(rng_local);
          enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt);
          arrivalT[i] = enrollt;

          u = unif(rng_local);
          int j = findInterval1(u, cumF);
          stratum[i] = j + 1;

          // stratified block randomization among 3 arms
          u = unif(rng_local);
          double denom = static_cast<double>(b1[j] + b2[j] + b3[j]);
          double p1 = static_cast<double>(b1[j]) / denom;
          double p2 = static_cast<double>(b1[j] + b2[j]) / denom;
          if (u <= p1) { trtGrp[i] = 1; --b1[j]; }
          else if (u <= p2) { trtGrp[i] = 2; --b2[j]; }
          else { trtGrp[i] = 3; --b3[j]; }
          if (b1[j] + b2[j] + b3[j] == 0) {
            b1[j] = allocation1; b2[j] = allocation2; b3[j] = allocation3;
          }

          // correlated normals -> uniforms
          double z1 = norm(rng_local);
          double z2 = norm(rng_local);
          double u1 = boost_pnorm(z1);
          double u2 = boost_pnorm(rho_pd_os * z1 + rho_pd_os_pyth_comp * z2);

          // survival times
          if (trtGrp[i] == 1) {
            survivalT1[i] = qtpwexpcpp1(u1, tau1pdx[j], lambda1pd[j]);
            survivalT2[i] = qtpwexpcpp1(u2, tau, lambda1osx[j]);
          } else if (trtGrp[i] == 2) {
            survivalT1[i] = qtpwexpcpp1(u1, tau2pdx[j], lambda2pd[j]);
            survivalT2[i] = qtpwexpcpp1(u2, tau, lambda2osx[j]);
          } else {
            survivalT1[i] = qtpwexpcpp1(u1, tau3pdx[j], lambda3pd[j]);
            survivalT2[i] = qtpwexpcpp1(u2, tau, lambda3osx[j]);
          }
          if (survivalT1[i] > survivalT2[i]) survivalT1[i] = survivalT2[i];

          // dropout times (independent)
          u1 = unif(rng_local);
          u2 = unif(rng_local);
          if (trtGrp[i] == 1) {
            dropoutT1[i] = qtpwexpcpp1(u1, tau, gamma1pd[j]);
            dropoutT2[i] = qtpwexpcpp1(u2, tau, gamma1osx[j]);
          } else if (trtGrp[i] == 2) {
            dropoutT1[i] = qtpwexpcpp1(u1, tau, gamma2pd[j]);
            dropoutT2[i] = qtpwexpcpp1(u2, tau, gamma2osx[j]);
          } else {
            dropoutT1[i] = qtpwexpcpp1(u1, tau, gamma3pd[j]);
            dropoutT2[i] = qtpwexpcpp1(u2, tau, gamma3osx[j]);
          }
          if (dropoutT1[i] > dropoutT2[i]) dropoutT1[i] = dropoutT2[i];

          // initial observed times/events (both endpoints)
          double sv1 = survivalT1[i], sv2 = survivalT2[i];
          double dr1 = dropoutT1[i], dr2 = dropoutT2[i];
          if (fixedFollowup) {
            if (sv1 <= dr1 && sv1 <= fu) {
              timeObs1[i] = sv1; event1[i] = 1; dropEv1[i] = 0;
            } else if (dr1 <= sv1 && dr1 <= fu) {
              timeObs1[i] = dr1; event1[i] = 0; dropEv1[i] = 1;
            } else {
              timeObs1[i] = fu; event1[i] = 0; dropEv1[i] = 0;
            }
            if (sv2 <= dr2 && sv2 <= fu) {
              timeObs2[i] = sv2; event2[i] = 1; dropEv2[i] = 0;
            } else if (dr2 <= sv2 && dr2 <= fu) {
              timeObs2[i] = dr2; event2[i] = 0; dropEv2[i] = 1;
            } else {
              timeObs2[i] = fu; event2[i] = 0; dropEv2[i] = 0;
            }
          } else {
            if (sv1 <= dr1) {
              timeObs1[i] = sv1; event1[i] = 1; dropEv1[i] = 0;
            } else {
              timeObs1[i] = dr1; event1[i] = 0; dropEv1[i] = 1;
            }
            if (sv2 <= dr2) {
              timeObs2[i] = sv2; event2[i] = 1; dropEv2[i] = 0;
            } else {
              timeObs2[i] = dr2; event2[i] = 0; dropEv2[i] = 1;
            }
          }
          totalT1[i] = arrivalT[i] + timeObs1[i];
          totalT2[i] = arrivalT[i] + timeObs2[i];
        } // cohort generation

        // determine analysis times and nstages
        size_t nstages = K;
        bool ev1NotAch = false;
        bool ev2NotAch = false;

        if (useEvents) {
          totalte1.clear(); totalte2.clear();
          int nevents1 = 0, nevents2 = 0;
          for (size_t i = 0; i < N; ++i) {
            if (event1[i] && (trtGrp[i] == 1 || trtGrp[i] == 3)) {
              ++nevents1; totalte1.push_back(totalT1[i]);
            }
            if (event2[i] && (trtGrp[i] == 1 || trtGrp[i] == 3)) {
              ++nevents2; totalte2.push_back(totalT2[i]);
            }
          }
          if (nevents1 == 0 || nevents2 == 0) {
            thread_utils::push_thread_warning(
              std::string("No events for iteration ") + std::to_string(iter + 1) +
                " skipping this iteration.");
            out.summaryRows.clear();
            out.rawRows.clear();
            continue;
          }
          std::sort(totalte1.begin(), totalte1.end());
          std::sort(totalte2.begin(), totalte2.end());


          // PFS looks
          analysisT1.clear();
          size_t j1 = 0;
          if (Kpfs > 0) {
            for (j1 = 0; j1 < Kpfs; ++j1) {
              if (plannedEvents[j1] >= nevents1) break;
            }

            if (j1 == Kpfs) { // total number of PFS events exceeds planned
              for (size_t k = 0; k < Kpfs; ++k) {
                analysisT1.push_back(totalte1[plannedEvents[k] - 1] + 1e-12);
              }
            } else {
              for (size_t k = 0; k < j1; ++k) {
                analysisT1.push_back(totalte1[plannedEvents[k] - 1] + 1e-12);
              }
              analysisT1.push_back(totalte1.back() + 1e-12);
            }
          }

          // OS looks -> compute analysisT2
          analysisT2.clear();
          size_t j2 = 0;
          if (K > Kpfs) {
            for (j2 = 0; j2 < (K - Kpfs); ++j2) {
              if (plannedEvents[Kpfs + j2] >= nevents2) break;
            }

            if (j2 == (K - Kpfs)) { // total number of OS events exceeds planned
              for (size_t k = 0; k < (K - Kpfs); ++k) {
                analysisT2.push_back(totalte2[plannedEvents[Kpfs + k] - 1] + 1e-12);
              }
            } else {
              for (size_t k = 0; k < j2; ++k) {
                analysisT2.push_back(totalte2[plannedEvents[Kpfs + k] - 1] + 1e-12);
              }
              analysisT2.push_back(totalte2.back() + 1e-12);
            }
          }

          // combine PFS and OS looks to determine nstages and analysisTime array
          if (Kpfs == 0) { // only OS looks
            nstages = analysisT2.size();
            std::copy_n(analysisT2.begin(), nstages, analysisT.begin());
          } else if (K == Kpfs) { // only PFS looks
            nstages = analysisT1.size();
            std::copy_n(analysisT1.begin(), nstages, analysisT.begin());
          } else { // mixed
            if (analysisT2.back() > analysisT1.back()) {
              // OS looks after last PFS look contribute.
              // NOTE: In this case, the observed number of PFS events must exceed
              // the planned number of PFS events at look Kpfs, because otherwise
              // the last PFS event would be observed at analysisT1.back().
              // However, since the last OS event occurred on or after
              // analysisT2.back() > analysisT1.back(), this is a
              // contradiction as death is part of PFS event definition.
              // It follows that analysisT1.size() == Kpfs in this case.

              // find first OS look after last PFS look
              size_t l = 0;
              for (size_t idx = 0; idx < analysisT2.size(); ++idx) {
                if (analysisT2[idx] > analysisT1.back()) { l = idx; break; }
              }
              // number of stages
              nstages = Kpfs + (analysisT2.size() - l);
              // copy PFS looks unchanged and append relevant OS looks
              // keep PFS looks [0 .. Kpfs-1], then OS looks from l onwards,
              // which are the ones after last PFS look mapped to Kpfs + l onwards
              std::copy_n(analysisT1.begin(), Kpfs, analysisT.begin());
              size_t count = analysisT2.size() - l;
              std::copy_n(analysisT2.begin() + l, count, analysisT.begin() + Kpfs);
            } else {
              // only PFS looks matter
              nstages = analysisT1.size();
              std::copy_n(analysisT1.begin(), nstages, analysisT.begin());
            }
          }

          // evNotAch: check PFS and OS targetse
          if (Kpfs > 0 && nevents1 < plannedEvents[Kpfs - 1]) ev1NotAch = true;
          if (Kpfs < K && nevents2 < plannedEvents[K - 1]) ev2NotAch = true;
        } else { // calendar time
          std::copy_n(plannedTime.begin(), K, analysisT.begin());
        }

        // per-stage computations
        for (size_t k = 0; k < nstages; ++k) {
          double time = analysisT[k];

          std::fill(n1x.begin(), n1x.end(), 0);
          std::fill(n2x.begin(), n2x.end(), 0);
          std::fill(n3x.begin(), n3x.end(), 0);

          int events1e1 = 0, events2e1 = 0, events3e1 = 0;
          int dropouts1e1 = 0, dropouts2e1 = 0, dropouts3e1 = 0;
          int events1e2 = 0, events2e2 = 0, events3e2 = 0;
          int dropouts1e2 = 0, dropouts2e2 = 0, dropouts3e2 = 0;

          // censoring & counts
          for (size_t i = 0; i < N; ++i) {
            double ar = arrivalT[i];
            double sv1 = survivalT1[i], sv2 = survivalT2[i];
            double dr1 = dropoutT1[i], dr2 = dropoutT2[i];

            if (ar > time) {
              timeObs1[i] = time - ar; event1[i] = 0; dropEv1[i] = 0;
              timeObs2[i] = time - ar; event2[i] = 0; dropEv2[i] = 0;
              continue;
            }

            // endpoint 1 censoring
            if (fixedFollowup) {
              if (ar + sv1 <= time && sv1 <= dr1 && sv1 <= fu) {
                timeObs1[i] = sv1; event1[i] = 1; dropEv1[i] = 0;
              } else if (ar + dr1 <= time && dr1 <= sv1 && dr1 <= fu) {
                timeObs1[i] = dr1; event1[i] = 0; dropEv1[i] = 1;
              } else if (ar + fu <= time && fu <= sv1 && fu <= dr1) {
                timeObs1[i] = fu; event1[i] = 0; dropEv1[i] = 0;
              } else {
                timeObs1[i] = time - ar; event1[i] = 0; dropEv1[i] = 0;
              }
            } else {
              if (ar + sv1 <= time && sv1 <= dr1) {
                timeObs1[i] = sv1; event1[i] = 1; dropEv1[i] = 0;
              } else if (ar + dr1 <= time && dr1 <= sv1) {
                timeObs1[i] = dr1; event1[i] = 0; dropEv1[i] = 1;
              } else {
                timeObs1[i] = time - ar; event1[i] = 0; dropEv1[i] = 0;
              }
            }

            // endpoint2 censoring
            if (fixedFollowup) {
              if (ar + sv2 <= time && sv2 <= dr2 && sv2 <= fu) {
                timeObs2[i] = sv2; event2[i] = 1; dropEv2[i] = 0;
              } else if (ar + dr2 <= time && dr2 <= sv2 && dr2 <= fu) {
                timeObs2[i] = dr2; event2[i] = 0; dropEv2[i] = 1;
              } else if (ar + fu <= time && fu <= sv2 && fu <= dr2) {
                timeObs2[i] = fu; event2[i] = 0; dropEv2[i] = 0;
              } else {
                timeObs2[i] = time - ar; event2[i] = 0; dropEv2[i] = 0;
              }
            } else {
              if (ar + sv2 <= time && sv2 <= dr2) {
                timeObs2[i] = sv2; event2[i] = 1; dropEv2[i] = 0;
              } else if (ar + dr2 <= time && dr2 <= sv2) {
                timeObs2[i] = dr2; event2[i] = 0; dropEv2[i] = 1;
              } else {
                timeObs2[i] = time - ar; event2[i] = 0; dropEv2[i] = 0;
              }
            }

            size_t h = static_cast<size_t>(stratum[i] - 1);
            if (trtGrp[i] == 1) { ++n1x[h];
              if (event1[i]) ++events1e1; else if (dropEv1[i]) ++dropouts1e1;
              if (event2[i]) ++events1e2; else if (dropEv2[i]) ++dropouts1e2;
            } else if (trtGrp[i] == 2) { ++n2x[h];
              if (event1[i]) ++events2e1; else if (dropEv1[i]) ++dropouts2e1;
              if (event2[i]) ++events2e2; else if (dropEv2[i]) ++dropouts2e2;
            } else { ++n3x[h];
              if (event1[i]) ++events3e1; else if (dropEv1[i]) ++dropouts3e1;
              if (event2[i]) ++events3e2; else if (dropEv2[i]) ++dropouts3e2;
            }
          } // censoring loop

          int accruals1 = std::accumulate(n1x.begin(), n1x.end(), 0);
          int accruals2 = std::accumulate(n2x.begin(), n2x.end(), 0);
          int accruals3 = std::accumulate(n3x.begin(), n3x.end(), 0);
          int totAccruals = accruals1 + accruals2 + accruals3;

          int totEventse1 = events1e1 + events2e1 + events3e1;
          int totDropoutse1 = dropouts1e1 + dropouts2e1 + dropouts3e1;
          int totEventse2 = events1e2 + events2e2 + events3e2;
          int totDropoutse2 = dropouts1e2 + dropouts2e2 + dropouts3e2;

          // optionally append raw rows (first maxRawIters iterations only)
          if (iter < maxRawIters) {
            for (size_t i = 0; i < N; ++i) {
              // skip subjects who haven't been enrolled by analysis time
              if (arrivalT[i] > time) continue;
              RawDatasetRow rr1;
              rr1.iterNum = static_cast<int>(iter + 1);
              rr1.stageNum = static_cast<int>(k + 1);
              rr1.analysisT = time;
              rr1.subjectId = static_cast<int>(i + 1);
              rr1.arrivalT = arrivalT[i];
              rr1.stratum = stratum[i];
              rr1.trtGrp = trtGrp[i];
              rr1.endpt = 1;
              rr1.survivalT = survivalT1[i];
              rr1.dropoutT = dropoutT1[i];
              rr1.timeObs = timeObs1[i];
              rr1.event = event1[i];
              rr1.dropEv = dropEv1[i];
              out.rawRows.push_back(std::move(rr1));

              RawDatasetRow rr2;
              rr2.iterNum = static_cast<int>(iter + 1);
              rr2.stageNum = static_cast<int>(k + 1);
              rr2.analysisT = time;
              rr2.subjectId = static_cast<int>(i + 1);
              rr2.arrivalT = arrivalT[i];
              rr2.stratum = stratum[i];
              rr2.trtGrp = trtGrp[i];
              rr2.endpt = 2;
              rr2.survivalT = survivalT2[i];
              rr2.dropoutT = dropoutT2[i];
              rr2.timeObs = timeObs2[i];
              rr2.event = event2[i];
              rr2.dropEv = dropEv2[i];
              out.rawRows.push_back(std::move(rr2));
            }
          }

          // compute pairwise stratified log-rank statistics for this endpoint
          for (int endpt = 1; endpt <= 2; ++endpt) {
            double h13, h23, h12;
            sub.clear();
            if (endpt == 1) {
              h13 = hazardRatioH013pfs;
              h23 = hazardRatioH023pfs;
              h12 = hazardRatioH012pfs;
              for (size_t i = 0; i < N; ++i) {
                if (timeObs1[i] > 0.0) sub.push_back(i);
              }
              std::sort(sub.begin(), sub.end(), [&](size_t a, size_t b) {
                return timeObs1[a] < timeObs1[b];
              });
            } else {
              h13 = hazardRatioH013os;
              h23 = hazardRatioH023os;
              h12 = hazardRatioH012os;
              for (size_t i = 0; i < N; ++i) {
                if (timeObs2[i] > 0.0) sub.push_back(i);
              }
              std::sort(sub.begin(), sub.end(), [&](size_t a, size_t b) {
                return timeObs2[a] < timeObs2[b];
              });
            }

            // restore risk sets
            n1 = n1x; n2 = n2x; n3 = n3x;


            std::fill(km13.begin(), km13.end(), 1.0);
            std::fill(km23.begin(), km23.end(), 1.0);
            std::fill(km12.begin(), km12.end(), 1.0);

            double us13 = 0.0, vs13 = 0.0;
            double us23 = 0.0, vs23 = 0.0;
            double us12 = 0.0, vs12 = 0.0;

            for (size_t i = 0; i < sub.size(); ++i) {
              size_t idx = sub[i];
              size_t h = static_cast<size_t>(stratum[idx] - 1);

              double n1h = static_cast<double>(n1[h]);
              double n2h = static_cast<double>(n2[h]);
              double n3h = static_cast<double>(n3[h]);

              double n13a = n1h * h13;
              double n23a = n2h * h23;
              double n12a = n1h * h12;

              double nt13 = n1h + n3h;
              double nt23 = n2h + n3h;
              double nt12 = n1h + n2h;

              double nt13a = n13a + n3h;
              double nt23a = n23a + n3h;
              double nt12a = n12a + n2h;

              bool evt = (endpt == 1 ? event1[idx] : event2[idx]);
              if (evt) {
                // 1 vs 3
                if (trtGrp[idx] == 1 || trtGrp[idx] == 3) {
                  double wh = 1.0;
                  if (rho1 != 0.0 || rho2 != 0.0) {
                    wh = std::pow(km13[h], rho1) * std::pow(1.0 - km13[h], rho2);
                    km13[h] *= (1.0 - 1.0 / nt13);
                  }
                  double treated = (trtGrp[idx] == 1 ? 1.0 : 0.0);
                  us13 += wh * (treated - n13a / nt13a);
                  vs13 += wh * wh * n13a * n3h / (nt13a * nt13a);
                }
                // 2 vs 3
                if (trtGrp[idx] == 2 || trtGrp[idx] == 3) {
                  double wh = 1.0;
                  if (rho1 != 0.0 || rho2 != 0.0) {
                    wh = std::pow(km23[h], rho1) * std::pow(1.0 - km23[h], rho2);
                    km23[h] *= (1.0 - 1.0 / nt23);
                  }
                  double treated = (trtGrp[idx] == 2 ? 1.0 : 0.0);
                  us23 += wh * (treated - n23a / nt23a);
                  vs23 += wh * wh * n23a * n3h / (nt23a * nt23a);
                }
                // 1 vs 2
                if (trtGrp[idx] == 1 || trtGrp[idx] == 2) {
                  double wh = 1.0;
                  if (rho1 != 0.0 || rho2 != 0.0) {
                    wh = std::pow(km12[h], rho1) * std::pow(1.0 - km12[h], rho2);
                    km12[h] *= (1.0 - 1.0 / nt12);
                  }
                  double treated = (trtGrp[idx] == 1 ? 1.0 : 0.0);
                  us12 += wh * (treated - n12a / nt12a);
                  vs12 += wh * wh * n12a * n2h / (nt12a * nt12a);
                }
              }

              // reduce risk set
              if (trtGrp[idx] == 1) --n1[h];
              else if (trtGrp[idx] == 2) --n2[h];
              else --n3[h];
            } // events loop

            double z13 = (vs13 > 0.0 ? (us13 / std::sqrt(vs13)) : 0.0);
            double z23 = (vs23 > 0.0 ? (us23 / std::sqrt(vs23)) : 0.0);
            double z12 = (vs12 > 0.0 ? (us12 / std::sqrt(vs12)) : 0.0);

            // append summary row
            StageSummaryRow sr;
            sr.iterNum = static_cast<int>(iter + 1);
            sr.evNotAch1 = ev1NotAch ? 1 : 0;
            sr.evNotAch2 = ev2NotAch ? 1 : 0;
            sr.stageNum = static_cast<int>(k + 1);
            sr.analysisT = time;
            sr.accruals1 = accruals1;
            sr.accruals2 = accruals2;
            sr.accruals3 = accruals3;
            sr.totAccruals = totAccruals;
            sr.endpt = endpt;
            if (endpt == 1) {
              sr.events1 = events1e1;
              sr.events2 = events2e1;
              sr.events3 = events3e1;
              sr.totEvents = totEventse1;
              sr.dropouts1 = dropouts1e1;
              sr.dropouts2 = dropouts2e1;
              sr.dropouts3 = dropouts3e1;
              sr.totDropouts = totDropoutse1;
            } else {
              sr.events1 = events1e2;
              sr.events2 = events2e2;
              sr.events3 = events3e2;
              sr.totEvents = totEventse2;
              sr.dropouts1 = dropouts1e2;
              sr.dropouts2 = dropouts2e2;
              sr.dropouts3 = dropouts3e2;
              sr.totDropouts = totDropoutse2;
            }
            sr.uscore13 = us13; sr.vscore13 = vs13; sr.logRank13 = z13;
            sr.uscore23 = us23; sr.vscore23 = vs23; sr.logRank23 = z23;
            sr.uscore12 = us12; sr.vscore12 = vs12; sr.logRank12 = z12;
            out.summaryRows.push_back(std::move(sr));
          } // endpoints loop
        } // per-stage
      } // iter
    } // operator()
  }; // SimWorker

  // construct and run worker
  SimWorker worker(
      K, Kpfs,
      hazardRatioH013pfs, hazardRatioH023pfs, hazardRatioH012pfs,
      hazardRatioH013os, hazardRatioH023os, hazardRatioH012os,
      allocation1, allocation2, allocation3,
      accrualTime, accrualIntensity, tau, stratumFraction,
      rho_pd_os,
      lambda1pfsx, lambda2pfsx, lambda3pfsx,
      lambda1osx, lambda2osx, lambda3osx,
      gamma1pfsx, gamma2pfsx, gamma3pfsx,
      gamma1osx, gamma2osx, gamma3osx,
      tau1pdx, tau2pdx, tau3pdx,
      lambda1pd, lambda2pd, lambda3pd,
      gamma1pd, gamma2pd, gamma3pd,
      N, fu, fixedFollowup, rho1, rho2,
      plannedEvents, plannedTime,
      maxIters, maxRawIters, seeds, useEvents, nstrata, rho_pd_os_pyth_comp,
      &results
  );

  RcppParallel::parallelFor(0, maxIters, worker);

  // Flatten per-iteration results
  size_t nsr = 0, nrr = 0;
  for (size_t iter = 0; iter < maxIters; ++iter) {
    nsr += results[iter].summaryRows.size();
    nrr += results[iter].rawRows.size();
  }
  if (nsr == 0) throw std::runtime_error(
    "No iterations with observed events. Unable to produce output.");

  // Final containers
  std::vector<int> sum_iterNum; sum_iterNum.reserve(nsr);
  std::vector<unsigned char> sum_ev1NotAch; sum_ev1NotAch.reserve(nsr);
  std::vector<unsigned char> sum_ev2NotAch; sum_ev2NotAch.reserve(nsr);
  std::vector<int> sum_stageNum; sum_stageNum.reserve(nsr);
  std::vector<double> sum_analysisT; sum_analysisT.reserve(nsr);
  std::vector<int> sum_accruals1; sum_accruals1.reserve(nsr);
  std::vector<int> sum_accruals2; sum_accruals2.reserve(nsr);
  std::vector<int> sum_accruals3; sum_accruals3.reserve(nsr);
  std::vector<int> sum_totAccruals; sum_totAccruals.reserve(nsr);
  std::vector<std::string> sum_endpt; sum_endpt.reserve(nsr);
  std::vector<int> sum_events1; sum_events1.reserve(nsr);
  std::vector<int> sum_events2; sum_events2.reserve(nsr);
  std::vector<int> sum_events3; sum_events3.reserve(nsr);
  std::vector<int> sum_totEvents; sum_totEvents.reserve(nsr);
  std::vector<int> sum_dropouts1; sum_dropouts1.reserve(nsr);
  std::vector<int> sum_dropouts2; sum_dropouts2.reserve(nsr);
  std::vector<int> sum_dropouts3; sum_dropouts3.reserve(nsr);
  std::vector<int> sum_totDropouts; sum_totDropouts.reserve(nsr);
  std::vector<double> sum_uscore13; sum_uscore13.reserve(nsr);
  std::vector<double> sum_vscore13; sum_vscore13.reserve(nsr);
  std::vector<double> sum_logRank13; sum_logRank13.reserve(nsr);
  std::vector<double> sum_uscore23; sum_uscore23.reserve(nsr);
  std::vector<double> sum_vscore23; sum_vscore23.reserve(nsr);
  std::vector<double> sum_logRank23; sum_logRank23.reserve(nsr);
  std::vector<double> sum_uscore12; sum_uscore12.reserve(nsr);
  std::vector<double> sum_vscore12; sum_vscore12.reserve(nsr);
  std::vector<double> sum_logRank12; sum_logRank12.reserve(nsr);

  // raw final containers
  std::vector<int> raw_iterNum; raw_iterNum.reserve(nrr);
  std::vector<int> raw_stageNum; raw_stageNum.reserve(nrr);
  std::vector<double> raw_analysisT; raw_analysisT.reserve(nrr);
  std::vector<int> raw_subjectId; raw_subjectId.reserve(nrr);
  std::vector<double> raw_arrivalT; raw_arrivalT.reserve(nrr);
  std::vector<int> raw_stratum; raw_stratum.reserve(nrr);
  std::vector<int> raw_trtGrp; raw_trtGrp.reserve(nrr);
  std::vector<std::string> raw_endpt; raw_endpt.reserve(nrr);
  std::vector<double> raw_survivalT; raw_survivalT.reserve(nrr);
  std::vector<double> raw_dropoutT; raw_dropoutT.reserve(nrr);
  std::vector<double> raw_timeObs; raw_timeObs.reserve(nrr);
  std::vector<unsigned char> raw_event; raw_event.reserve(nrr);
  std::vector<unsigned char> raw_dropEv; raw_dropEv.reserve(nrr);

  for (size_t iter = 0; iter < maxIters; ++iter) {
    const auto& srows = results[iter].summaryRows;
    for (const auto& r : srows) {
      sum_iterNum.push_back(r.iterNum);
      sum_ev1NotAch.push_back(r.evNotAch1);
      sum_ev2NotAch.push_back(r.evNotAch2);
      sum_stageNum.push_back(r.stageNum);
      sum_analysisT.push_back(r.analysisT);
      sum_accruals1.push_back(r.accruals1);
      sum_accruals2.push_back(r.accruals2);
      sum_accruals3.push_back(r.accruals3);
      sum_totAccruals.push_back(r.totAccruals);
      sum_endpt.push_back(r.endpt == 1 ? "PFS" : "OS");
      sum_events1.push_back(r.events1);
      sum_events2.push_back(r.events2);
      sum_events3.push_back(r.events3);
      sum_totEvents.push_back(r.totEvents);
      sum_dropouts1.push_back(r.dropouts1);
      sum_dropouts2.push_back(r.dropouts2);
      sum_dropouts3.push_back(r.dropouts3);
      sum_totDropouts.push_back(r.totDropouts);
      sum_uscore13.push_back(r.uscore13);
      sum_vscore13.push_back(r.vscore13);
      sum_logRank13.push_back(r.logRank13);
      sum_uscore23.push_back(r.uscore23);
      sum_vscore23.push_back(r.vscore23);
      sum_logRank23.push_back(r.logRank23);
      sum_uscore12.push_back(r.uscore12);
      sum_vscore12.push_back(r.vscore12);
      sum_logRank12.push_back(r.logRank12);
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
        raw_endpt.push_back(rr.endpt == 1 ? "PFS" : "OS");
        raw_survivalT.push_back(rr.survivalT);
        raw_dropoutT.push_back(rr.dropoutT);
        raw_timeObs.push_back(rr.timeObs);
        raw_event.push_back(rr.event);
        raw_dropEv.push_back(rr.dropEv);
      }
    }
  }

  // Build summary DataFrameCpp
  DataFrameCpp sumdata;
  sumdata.push_back(std::move(sum_iterNum), "iterationNumber");
  sumdata.push_back(std::move(sum_ev1NotAch), "events1NotAchieved");
  sumdata.push_back(std::move(sum_ev2NotAch), "events2NotAchieved");
  sumdata.push_back(std::move(sum_stageNum), "stageNumber");
  sumdata.push_back(std::move(sum_analysisT), "analysisTime");
  sumdata.push_back(std::move(sum_accruals1), "accruals1");
  sumdata.push_back(std::move(sum_accruals2), "accruals2");
  sumdata.push_back(std::move(sum_accruals3), "accruals3");
  sumdata.push_back(std::move(sum_totAccruals), "totalAccruals");
  sumdata.push_back(std::move(sum_endpt), "endpoint");
  sumdata.push_back(std::move(sum_events1), "events1");
  sumdata.push_back(std::move(sum_events2), "events2");
  sumdata.push_back(std::move(sum_events3), "events3");
  sumdata.push_back(std::move(sum_totEvents), "totalEvents");
  sumdata.push_back(std::move(sum_dropouts1), "dropouts1");
  sumdata.push_back(std::move(sum_dropouts2), "dropouts2");
  sumdata.push_back(std::move(sum_dropouts3), "dropouts3");
  sumdata.push_back(std::move(sum_totDropouts), "totalDropouts");
  sumdata.push_back(std::move(sum_uscore13), "uscore13");
  sumdata.push_back(std::move(sum_vscore13), "vscore13");
  sumdata.push_back(std::move(sum_logRank13), "logRankStatistic13");
  sumdata.push_back(std::move(sum_uscore23), "uscore23");
  sumdata.push_back(std::move(sum_vscore23), "vscore23");
  sumdata.push_back(std::move(sum_logRank23), "logRankStatistic23");
  sumdata.push_back(std::move(sum_uscore12), "uscore12");
  sumdata.push_back(std::move(sum_vscore12), "vscore12");
  sumdata.push_back(std::move(sum_logRank12), "logRankStatistic12");

  ListCpp result;
  result.push_back(sumdata, "sumdata");

  if (!raw_iterNum.empty()) {
    DataFrameCpp rawdata;
    rawdata.push_back(std::move(raw_iterNum), "iterationNumber");
    rawdata.push_back(std::move(raw_stageNum), "stageNumber");
    rawdata.push_back(std::move(raw_analysisT), "analysisTime");
    rawdata.push_back(std::move(raw_subjectId), "subjectId");
    rawdata.push_back(std::move(raw_arrivalT), "arrivalTime");
    rawdata.push_back(std::move(raw_stratum), "stratum");
    rawdata.push_back(std::move(raw_trtGrp), "treatmentGroup");
    rawdata.push_back(std::move(raw_endpt), "endpoint");
    rawdata.push_back(std::move(raw_survivalT), "survivalTime");
    rawdata.push_back(std::move(raw_dropoutT), "dropoutTime");
    rawdata.push_back(std::move(raw_timeObs), "timeUnderObservation");
    rawdata.push_back(std::move(raw_event), "event");
    rawdata.push_back(std::move(raw_dropEv), "dropoutEvent");
    result.push_back(rawdata, "rawdata");
  }

  return result;
}


//' @title Log-Rank Test Simulation for Two Endpoints (PFS and OS) and
//' Three Arms
//' @description Performs simulation for two-endpoint (PFS and OS)
//' three-arm group sequential trials based on weighted log-rank test.
//' The first \code{kMaxpfs} looks are driven by the total number of
//' PFS events in Arm A and Arm C combined, and the subsequent looks
//' are driven by the total number of OS events in Arm A and Arm C
//' combined. Alternatively, the analyses can be planned to occur at
//' specified calendar times.
//'
//' @inheritParams param_kMax
//' @param kMaxpfs Number of stages with timing determined by PFS events.
//'   Ranges from 0 (none) to \code{kMax}.
//' @param hazardRatioH013pfs Hazard ratio under the null hypothesis for arm 1
//'   vs arm 3 for PFS. Defaults to 1 for superiority test.
//' @param hazardRatioH023pfs Hazard ratio under the null hypothesis for arm 2
//'   vs arm 3 for PFS. Defaults to 1 for superiority test.
//' @param hazardRatioH012pfs Hazard ratio under the null hypothesis for arm 1
//'   vs arm 2 for PFS. Defaults to 1 for superiority test.
//' @param hazardRatioH013os Hazard ratio under the null hypothesis for arm 1
//'   vs arm 3 for OS. Defaults to 1 for superiority test.
//' @param hazardRatioH023os Hazard ratio under the null hypothesis for arm 2
//'   vs arm 3 for OS. Defaults to 1 for superiority test.
//' @param hazardRatioH012os Hazard ratio under the null hypothesis for arm 1
//'   vs arm 2 for OS. Defaults to 1 for superiority test.
//' @param allocation1 Number of subjects in Arm A in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in Arm B in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation3 Number of subjects in Arm C in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param rho_pd_os The correlation coefficient for the standard
//'   bivariate normal random variables used to generate time to
//'   disease progression and time to death using the inverse CDF method.
//' @param lambda1pfs A vector of hazard rates for the event in each analysis
//'   time interval by stratum for arm 1 and PFS.
//' @param lambda2pfs A vector of hazard rates for the event in each analysis
//'   time interval by stratum for arm 2 and PFS.
//' @param lambda3pfs A vector of hazard rates for the event in each analysis
//'   time interval by stratum for arm 3 and PFS.
//' @param lambda1os A vector of hazard rates for the event in each analysis
//'   time interval by stratum for arm 1 and OS.
//' @param lambda2os A vector of hazard rates for the event in each analysis
//'   time interval by stratum for arm 2 and OS.
//' @param lambda3os A vector of hazard rates for the event in each analysis
//'   time interval by stratum for arm 3 and OS.
//' @param gamma1pfs The hazard rate for exponential dropout. A vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for arm 1 and PFS.
//' @param gamma2pfs The hazard rate for exponential dropout. A vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for arm 2 and PFS.
//' @param gamma3pfs The hazard rate for exponential dropout. A vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for arm 3 and PFS.
//' @param gamma1os The hazard rate for exponential dropout. A vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for arm 1 and OS.
//' @param gamma2os The hazard rate for exponential dropout. A vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for arm 2 and OS.
//' @param gamma3os The hazard rate for exponential dropout. A vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for arm 3 and OS.
//' @param n Sample size.
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @param plannedEvents The planned cumulative total number of PFS events at
//'   Look 1 to Look \code{kMaxpfs} for Arms A and C combined and the planned
//'   cumulative total number of OS events at Look \code{kMaxpfs+1} to Look
//'   \code{kMax} for Arms A and C combined.
//' @param plannedTime The calendar times for the analyses. To use calendar
//'   time to plan the analyses, \code{plannedEvents} should be missing.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
//'   stage to extract.
//' @param seed The seed to reproduce the simulation results.
//'
//' @return A list with 2 components:
//'
//' * \code{sumdata}: A data frame of summary data by iteration and stage:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{eventsNotAchieved}: Whether the target number of events
//'       is not achieved for the iteration.
//'
//'     - \code{stageNumber}: The stage number, covering all stages even if
//'       the trial stops at an interim look.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{accruals1}: The number of subjects enrolled at the stage for
//'       the active treatment 1 group.
//'
//'     - \code{accruals2}: The number of subjects enrolled at the stage for
//'       the active treatment 2 group.
//'
//'     - \code{accruals3}: The number of subjects enrolled at the stage for
//'       the control group.
//'
//'     - \code{totalAccruals}: The total number of subjects enrolled at
//'       the stage.
//'
//'     - \code{endpoint}: The endpoint (1 for PFS or 2 for OS) under
//'       consideration.
//'
//'     - \code{events1}: The number of events at the stage for
//'       the active treatment 1 group.
//'
//'     - \code{events2}: The number of events at the stage for
//'       the active treatment 2 group.
//'
//'     - \code{events3}: The number of events at the stage for
//'       the control group.
//'
//'     - \code{totalEvents}: The total number of events at the stage.
//'
//'     - \code{dropouts1}: The number of dropouts at the stage for
//'       the active treatment 1 group.
//'
//'     - \code{dropouts2}: The number of dropouts at the stage for
//'       the active treatment 2 group.
//'
//'     - \code{dropouts3}: The number of dropouts at the stage for
//'       the control group.
//'
//'     - \code{totalDropouts}: The total number of dropouts at the stage.
//'
//'     - \code{logRankStatistic13}: The log-rank test Z-statistic
//'       comparing the active treatment 1 to the control for the endpoint.
//'
//'     - \code{logRankStatistic23}: The log-rank test Z-statistic
//'       comparing the active treatment 2 to the control for the endpoint.
//'
//'     - \code{logRankStatistic12}: The log-rank test Z-statistic
//'       comparing the active treatment 1 to the active treatment 2
//'       for the endpoint.
//'
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//'   positive integer): A data frame for subject-level data for selected
//'   replications, containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage under consideration.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{subjectId}: The subject ID.
//'
//'     - \code{arrivalTime}: The enrollment time for the subject.
//'
//'     - \code{stratum}: The stratum for the subject.
//'
//'     - \code{treatmentGroup}: The treatment group (1, 2, or 3) for
//'       the subject.
//'
//'     - \code{endpoint}: The endpoint (1 for PFS or 2 for OS) under
//'       consideration.
//'
//'     - \code{survivalTime}: The underlying survival time for the
//'       event endpoint for the subject.
//'
//'     - \code{dropoutTime}: The underlying dropout time for the
//'       event endpoint for the subject.
//'
//'     - \code{timeUnderObservation}: The time under observation
//'       since randomization for the event endpoint for the subject.
//'
//'     - \code{event}: Whether the subject experienced the event
//'       endpoint.
//'
//'     - \code{dropoutEvent}: Whether the subject dropped out for
//'       the endpoint.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' sim1 = lrsim2e3a(
//'   kMax = 3,
//'   kMaxpfs = 2,
//'   allocation1 = 2,
//'   allocation2 = 2,
//'   allocation3 = 1,
//'   accrualTime = c(0, 8),
//'   accrualIntensity = c(10, 28),
//'   piecewiseSurvivalTime = 0,
//'   rho_pd_os = 0,
//'   lambda1pfs = log(2)/12*0.60,
//'   lambda2pfs = log(2)/12*0.70,
//'   lambda3pfs = log(2)/12,
//'   lambda1os = log(2)/30*0.65,
//'   lambda2os = log(2)/30*0.75,
//'   lambda3os = log(2)/30,
//'   n = 700,
//'   plannedEvents = c(186, 259, 183),
//'   maxNumberOfIterations = 1000,
//'   maxNumberOfRawDatasetsPerStage = 1,
//'   seed = 314159)
//'
//' head(sim1$sumdata)
//' head(sim1$rawdata)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List lrsim2e3a(
    const int kMax = 1,
    const int kMaxpfs = 1,
    const double hazardRatioH013pfs = 1,
    const double hazardRatioH023pfs = 1,
    const double hazardRatioH012pfs = 1,
    const double hazardRatioH013os = 1,
    const double hazardRatioH023os = 1,
    const double hazardRatioH012os = 1,
    const int allocation1 = 1,
    const int allocation2 = 1,
    const int allocation3 = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const double rho_pd_os = 0,
    const Rcpp::NumericVector& lambda1pfs = NA_REAL,
    const Rcpp::NumericVector& lambda2pfs = NA_REAL,
    const Rcpp::NumericVector& lambda3pfs = NA_REAL,
    const Rcpp::NumericVector& lambda1os = NA_REAL,
    const Rcpp::NumericVector& lambda2os = NA_REAL,
    const Rcpp::NumericVector& lambda3os = NA_REAL,
    const Rcpp::NumericVector& gamma1pfs = 0,
    const Rcpp::NumericVector& gamma2pfs = 0,
    const Rcpp::NumericVector& gamma3pfs = 0,
    const Rcpp::NumericVector& gamma1os = 0,
    const Rcpp::NumericVector& gamma2os = 0,
    const Rcpp::NumericVector& gamma3os = 0,
    const int n = NA_INTEGER,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = 0,
    const double rho1 = 0,
    const double rho2 = 0,
    const Rcpp::IntegerVector& plannedEvents = NA_INTEGER,
    const Rcpp::NumericVector& plannedTime = NA_REAL,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = 0) {

  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvTime = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1pfs = Rcpp::as<std::vector<double>>(lambda1pfs);
  auto lam2pfs = Rcpp::as<std::vector<double>>(lambda2pfs);
  auto lam3pfs = Rcpp::as<std::vector<double>>(lambda3pfs);
  auto lam1os = Rcpp::as<std::vector<double>>(lambda1os);
  auto lam2os = Rcpp::as<std::vector<double>>(lambda2os);
  auto lam3os = Rcpp::as<std::vector<double>>(lambda3os);
  auto gam1pfs = Rcpp::as<std::vector<double>>(gamma1pfs);
  auto gam2pfs = Rcpp::as<std::vector<double>>(gamma2pfs);
  auto gam3pfs = Rcpp::as<std::vector<double>>(gamma3pfs);
  auto gam1os = Rcpp::as<std::vector<double>>(gamma1os);
  auto gam2os = Rcpp::as<std::vector<double>>(gamma2os);
  auto gam3os = Rcpp::as<std::vector<double>>(gamma3os);
  auto plannedE = Rcpp::as<std::vector<int>>(plannedEvents);
  auto plannedT = Rcpp::as<std::vector<double>>(plannedTime);

  auto out = lrsim2e3acpp(
    kMax, kMaxpfs, hazardRatioH013pfs, hazardRatioH023pfs,
    hazardRatioH012pfs, hazardRatioH013os, hazardRatioH023os,
    hazardRatioH012os, allocation1, allocation2, allocation3,
    accrualT, accrualInt, pwSurvTime, stratumFrac, rho_pd_os,
    lam1pfs, lam2pfs, lam3pfs, lam1os, lam2os, lam3os,
    gam1pfs, gam2pfs, gam3pfs, gam1os, gam2os, gam3os, n,
    followupTime, fixedFollowup, rho1, rho2, plannedE, plannedT,
    maxNumberOfIterations, maxNumberOfRawDatasetsPerStage, seed);

  thread_utils::drain_thread_warnings_to_R();

  return Rcpp::wrap(out);
}


ListCpp lrsimsubcpp(
    const int kMax,
    const int kMaxitt,
    const double hazardRatioH0itt,
    const double hazardRatioH0pos,
    const double hazardRatioH0neg,
    const int allocation1,
    const int allocation2,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& p_pos,
    const std::vector<double>& lambda1itt,
    const std::vector<double>& lambda2itt,
    const std::vector<double>& lambda1pos,
    const std::vector<double>& lambda2pos,
    const std::vector<double>& gamma1itt,
    const std::vector<double>& gamma2itt,
    const std::vector<double>& gamma1pos,
    const std::vector<double>& gamma2pos,
    const int n,
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
  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
  size_t K = static_cast<size_t>(kMax);

  int kMaxittx = kMaxitt;
  if (kMaxittx < 0) kMaxittx = kMax;
  if (kMaxittx > kMax)
    throw std::invalid_argument("kMaxitt must be less than or equal to kMax");
  size_t Kitt = static_cast<size_t>(kMaxittx);

  bool useEvents;
  if (none_na(plannedEvents)) {
    useEvents = true;
    if (plannedEvents[0] <= 0)
      throw std::invalid_argument("plannedEvents must be positive");
    if (plannedEvents.size() != K)
      throw std::invalid_argument("Invalid length for plannedEvents");
    if (Kitt > 1) {
      for (size_t i = 1; i < Kitt; ++i) {
        if (plannedEvents[i] <= plannedEvents[i-1])
          throw std::invalid_argument("plannedEvents for ITT must be increasing");
      }
    }
    if (K - Kitt > 1) {
      for (size_t i = Kitt + 1; i < plannedEvents.size(); ++i) {
        if (plannedEvents[i] <= plannedEvents[i-1])
          throw std::invalid_argument(
              "plannedEvents for biomarker+ must be increasing");
      }
    }
  } else if (none_na(plannedTime)) {
    useEvents = false;
    if (plannedTime[0] <= 0.0)
      throw std::invalid_argument("plannedTime must be positive");
    if (plannedTime.size() != K)
      throw std::invalid_argument("Invalid length for plannedTime");
    if (any_nonincreasing(plannedTime))
      throw std::invalid_argument("plannedTime must be increasing");
  } else {
    throw std::invalid_argument("Either plannedEvents or plannedTime must be given");
  }

  if (hazardRatioH0itt <= 0.0)
    throw std::invalid_argument("hazardRatioH0itt must be positive");
  if (hazardRatioH0pos <= 0.0)
    throw std::invalid_argument("hazardRatioH0pos must be positive");
  if (hazardRatioH0neg <= 0.0)
    throw std::invalid_argument("hazardRatioH0neg must be positive");
  if (allocation1 < 1 || allocation2 < 1)
    throw std::invalid_argument("allocations must be positive integers");
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
  if (!none_na(p_pos)) throw std::invalid_argument("p_pos must be provided");
  for (double v : p_pos) {
    if (!(v > 0.0 && v < 1.0))
      throw std::invalid_argument("p_pos must lie between 0 and 1");
  }
  if (!none_na(lambda1itt)) throw std::invalid_argument("lambda1itt must be provided");
  if (!none_na(lambda2itt)) throw std::invalid_argument("lambda2itt must be provided");
  if (!none_na(lambda1pos)) throw std::invalid_argument("lambda1pos must be provided");
  if (!none_na(lambda2pos)) throw std::invalid_argument("lambda2pos must be provided");
  for (double v : lambda1itt) {
    if (v < 0.0) throw std::invalid_argument("lambda1itt must be non-negative");
  }
  for (double v : lambda2itt) {
    if (v < 0.0) throw std::invalid_argument("lambda2itt must be non-negative");
  }
  for (double v : lambda1pos) {
    if (v < 0.0) throw std::invalid_argument("lambda1pos must be non-negative");
  }
  for (double v : lambda2pos) {
    if (v < 0.0) throw std::invalid_argument("lambda2pos must be non-negative");
  }
  for (double v : gamma1itt) {
    if (v < 0.0) throw std::invalid_argument("gamma1itt must be non-negative");
  }
  for (double v : gamma2itt) {
    if (v < 0.0) throw std::invalid_argument("gamma2itt must be non-negative");
  }
  for (double v : gamma1pos) {
    if (v < 0.0) throw std::invalid_argument("gamma1pos must be non-negative");
  }
  for (double v : gamma2pos) {
    if (v < 0.0) throw std::invalid_argument("gamma2pos must be non-negative");
  }
  if (n == INT_MIN) throw std::invalid_argument("n must be provided");
  if (n <= 0) throw std::invalid_argument("n must be positive");
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

  size_t N = static_cast<size_t>(n);
  size_t maxIters = static_cast<size_t>(maxNumberOfIterations);
  size_t maxRawIters = static_cast<size_t>(maxNumberOfRawDatasetsPerStage);
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  size_t nintv2 = (nintv == 1 ? 10u : nintv + 10u);
  const std::vector<double>& tau = piecewiseSurvivalTime;
  const double fu = followupTime;

  std::vector<double> p_posv = expand1(p_pos, nstrata, "p_pos");

  auto lam1ittx = expand_stratified(lambda1itt, nstrata, nintv, "lambda1itt");
  auto lam2ittx = expand_stratified(lambda2itt, nstrata, nintv, "lambda2itt");
  auto lam1posx = expand_stratified(lambda1pos, nstrata, nintv, "lambda1pos");
  auto lam2posx = expand_stratified(lambda2pos, nstrata, nintv, "lambda2pos");
  auto gam1ittx  = expand_stratified(gamma1itt,  nstrata, nintv, "gamma1itt");
  auto gam2ittx  = expand_stratified(gamma2itt,  nstrata, nintv, "gamma2itt");
  auto gam1posx  = expand_stratified(gamma1pos,  nstrata, nintv, "gamma1pos");
  auto gam2posx  = expand_stratified(gamma2pos,  nstrata, nintv, "gamma2pos");

  // compute subpopulation hazards via hazard_subcpp (main thread)
  std::vector<std::vector<double>> tau1pos(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> tau2pos(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> tau1neg(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> tau2neg(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> lam1posy(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> lam2posy(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> lam1negy(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> lam2negy(nstrata, std::vector<double>(nintv2));
  std::vector<std::vector<double>> gam1posy(nstrata, std::vector<double>(nintv));
  std::vector<std::vector<double>> gam2posy(nstrata, std::vector<double>(nintv));
  std::vector<std::vector<double>> gam1negy(nstrata, std::vector<double>(nintv));
  std::vector<std::vector<double>> gam2negy(nstrata, std::vector<double>(nintv));

  for (size_t s = 0; s < nstrata; ++s) {
    ListCpp a1 = hazard_subcpp(tau, lam1ittx[s], lam1posx[s], p_posv[s]);
    ListCpp a2 = hazard_subcpp(tau, lam2ittx[s], lam2posx[s], p_posv[s]);
    ListCpp b1 = hazard_subcpp(tau, gam1ittx[s], gam1posx[s], p_posv[s]);
    ListCpp b2 = hazard_subcpp(tau, gam2ittx[s], gam2posx[s], p_posv[s]);

    tau1pos[s] = a1.get<std::vector<double>>("piecewiseSurvivalTime");
    tau2pos[s] = a2.get<std::vector<double>>("piecewiseSurvivalTime");
    tau1neg[s] = b1.get<std::vector<double>>("piecewiseSurvivalTime");
    tau2neg[s] = b2.get<std::vector<double>>("piecewiseSurvivalTime");

    lam1posy[s] = a1.get<std::vector<double>>("hazard_pos");
    lam1negy[s] = a1.get<std::vector<double>>("hazard_neg");
    lam2posy[s] = a2.get<std::vector<double>>("hazard_pos");
    lam2negy[s] = a2.get<std::vector<double>>("hazard_neg");

    gam1posy[s] = b1.get<std::vector<double>>("hazard_pos");
    gam1negy[s] = b1.get<std::vector<double>>("hazard_neg");
    gam2posy[s] = b2.get<std::vector<double>>("hazard_pos");
    gam2negy[s] = b2.get<std::vector<double>>("hazard_neg");
  }

  // prepare per-iteration seed vector
  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter) seeds[iter] = master_rng();


  // One summary (stage-level) row produced by an iteration
  struct StageSummaryRow {
    int iterNum = 0;
    unsigned char ev1NotAch = 0, ev2NotAch = 0;
    int stageNum = 0;
    double analysisT = 0.0;
    int pop = 0;
    int accruals1 = 0, accruals2 = 0, totAccruals = 0;
    int events1 = 0, events2 = 0, totEvents = 0;
    int dropouts1 = 0, dropouts2 = 0, totDropouts = 0;
    double uscore = 0.0, vscore = 0.0, logRank = 0.0;
  };

  // One subject-level (raw) row for a particular iteration and stage
  struct RawDatasetRow {
    int iterNum = 0, stageNum = 0;
    double analysisT = 0.0;
    int subjectId = 0;
    double arrivalT = 0.0;
    int stratum = 0;
    unsigned char marker = 0;
    int trtGrp = 0;
    double survivalT = 0.0, dropoutT = 0.0, timeObs = 0.0;
    unsigned char event = 0, dropEv = 0;
  };

  struct IterationResult {
    std::vector<StageSummaryRow> summaryRows;
    std::vector<RawDatasetRow> rawRows;
    void reserveForSummary(size_t approxRows) { summaryRows.reserve(approxRows); }
    void reserveForRaw(size_t approxRows) { rawRows.reserve(approxRows); }
  };

  // pre-size results
  std::vector<IterationResult> results;
  results.resize(maxIters);

  // Worker (declared inside function)
  struct SimWorker : public RcppParallel::Worker {
    const size_t K;
    const size_t Kitt;
    const double hazardRatioH0itt;
    const double hazardRatioH0pos;
    const double hazardRatioH0neg;
    const int allocation1;
    const int allocation2;
    const std::vector<double>& accrualTime;
    const std::vector<double>& accrualIntensity;
    const std::vector<double>& tau;
    const std::vector<double>& stratumFraction;
    const std::vector<double>& p_posv;
    const std::vector<std::vector<double>>& lam1ittx;
    const std::vector<std::vector<double>>& lam2ittx;
    const std::vector<std::vector<double>>& lam1posx;
    const std::vector<std::vector<double>>& lam2posx;
    const std::vector<std::vector<double>>& gam1ittx;
    const std::vector<std::vector<double>>& gam2ittx;
    const std::vector<std::vector<double>>& gam1posx;
    const std::vector<std::vector<double>>& gam2posx;
    const std::vector<std::vector<double>>& tau1pos;
    const std::vector<std::vector<double>>& tau2pos;
    const std::vector<std::vector<double>>& tau1neg;
    const std::vector<std::vector<double>>& tau2neg;
    const std::vector<std::vector<double>>& lam1posy;
    const std::vector<std::vector<double>>& lam2posy;
    const std::vector<std::vector<double>>& lam1negy;
    const std::vector<std::vector<double>>& lam2negy;
    const std::vector<std::vector<double>>& gam1posy;
    const std::vector<std::vector<double>>& gam2posy;
    const std::vector<std::vector<double>>& gam1negy;
    const std::vector<std::vector<double>>& gam2negy;

    const size_t N;
    const double fu;
    const bool fixedFollowup;
    const double rho1;
    const double rho2;
    const std::vector<int>& plannedEvents;
    const std::vector<double>& plannedTime;
    const size_t maxIters;
    const size_t maxRawIters;
    const std::vector<uint64_t>& seeds;
    const bool useEvents;
    const size_t nstrata;

    std::vector<IterationResult>* results;

    SimWorker(
      size_t K_,
      size_t Kitt_,
      double hazardRatioH0itt_,
      double hazardRatioH0pos_,
      double hazardRatioH0neg_,
      int allocation1_,
      int allocation2_,
      const std::vector<double>& accrualTime_,
      const std::vector<double>& accrualIntensity_,
      const std::vector<double>& tau_,
      const std::vector<double>& stratumFraction_,
      const std::vector<double>& p_posv_,
      const std::vector<std::vector<double>>& lam1ittx_,
      const std::vector<std::vector<double>>& lam2ittx_,
      const std::vector<std::vector<double>>& lam1posx_,
      const std::vector<std::vector<double>>& lam2posx_,
      const std::vector<std::vector<double>>& gam1ittx_,
      const std::vector<std::vector<double>>& gam2ittx_,
      const std::vector<std::vector<double>>& gam1posx_,
      const std::vector<std::vector<double>>& gam2posx_,
      const std::vector<std::vector<double>>& tau1pos_,
      const std::vector<std::vector<double>>& tau2pos_,
      const std::vector<std::vector<double>>& tau1neg_,
      const std::vector<std::vector<double>>& tau2neg_,
      const std::vector<std::vector<double>>& lam1posy_,
      const std::vector<std::vector<double>>& lam2posy_,
      const std::vector<std::vector<double>>& lam1negy_,
      const std::vector<std::vector<double>>& lam2negy_,
      const std::vector<std::vector<double>>& gam1posy_,
      const std::vector<std::vector<double>>& gam2posy_,
      const std::vector<std::vector<double>>& gam1negy_,
      const std::vector<std::vector<double>>& gam2negy_,
      size_t N_,
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
      size_t nstrata_,
      std::vector<IterationResult>* results_)
      : K(K_),
        Kitt(Kitt_),
        hazardRatioH0itt(hazardRatioH0itt_),
        hazardRatioH0pos(hazardRatioH0pos_),
        hazardRatioH0neg(hazardRatioH0neg_),
        allocation1(allocation1_),
        allocation2(allocation2_),
        accrualTime(accrualTime_),
        accrualIntensity(accrualIntensity_),
        tau(tau_),
        stratumFraction(stratumFraction_),
        p_posv(p_posv_),
        lam1ittx(lam1ittx_),
        lam2ittx(lam2ittx_),
        lam1posx(lam1posx_),
        lam2posx(lam2posx_),
        gam1ittx(gam1ittx_),
        gam2ittx(gam2ittx_),
        gam1posx(gam1posx_),
        gam2posx(gam2posx_),
        tau1pos(tau1pos_),
        tau2pos(tau2pos_),
        tau1neg(tau1neg_),
        tau2neg(tau2neg_),
        lam1posy(lam1posy_),
        lam2posy(lam2posy_),
        lam1negy(lam1negy_),
        lam2negy(lam2negy_),
        gam1posy(gam1posy_),
        gam2posy(gam2posy_),
        gam1negy(gam1negy_),
        gam2negy(gam2negy_),
        N(N_),
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
        nstrata(nstrata_),
        results(results_)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      // local buffers
      std::vector<int> stratum(N), trtGrp(N);
      std::vector<double> arrivalT(N), survivalT(N), dropoutT(N);
      std::vector<double> timeObs(N), totalT(N);
      std::vector<unsigned char> marker(N), event(N), dropEv(N);

      std::vector<int> b1(nstrata), b2(nstrata);
      std::vector<int> n1(nstrata), n2(nstrata);
      std::vector<double> km(nstrata);
      std::vector<double> cumF(nstrata);
      std::partial_sum(stratumFraction.begin(), stratumFraction.end(), cumF.begin());

      std::vector<double> analysisT(K);
      std::vector<double> analysisT1; analysisT1.reserve(Kitt);
      std::vector<double> analysisT2; analysisT2.reserve(K - Kitt);
      std::vector<double> totalte; totalte.reserve(N);
      std::vector<double> totaltepos; totaltepos.reserve(N);
      std::vector<size_t> set; set.reserve(N);
      std::vector<size_t> sub; sub.reserve(N);

      for (size_t iter = begin; iter < end; ++iter) {
        // per-iteration RNG
        boost::random::mt19937_64 rng_local(seeds[iter]);
        boost::random::uniform_real_distribution<double> unif(0.0, 1.0);

        IterationResult& out = (*results)[iter];
        out.summaryRows.clear();
        out.rawRows.clear();
        if (iter < maxRawIters) out.reserveForRaw(K * N);
        out.reserveForSummary(K * 3);

        // reset block randomization
        std::fill(b1.begin(), b1.end(), allocation1);
        std::fill(b2.begin(), b2.end(), allocation2);

        double enrollt = 0.0;

        // generate cohort
        for (size_t i = 0; i < N; ++i) {
          double u = unif(rng_local);
          enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt);
          arrivalT[i] = enrollt;

          u = unif(rng_local);
          int j = findInterval1(u, cumF);
          stratum[i] = j + 1;

          // biomarker
          u = unif(rng_local);
          marker[i] = (u <= p_posv[j] ? 1 : 0);

          // stratified 2-arm randomization
          u = unif(rng_local);
          double denom = static_cast<double>(b1[j] + b2[j]);
          double p = static_cast<double>(b1[j]) / denom;
          if (u <= p) { trtGrp[i] = 1; --b1[j]; }
          else { trtGrp[i] = 2; --b2[j]; }
          if (b1[j] + b2[j] == 0) { b1[j] = allocation1; b2[j] = allocation2; }

          // survival time
          u = unif(rng_local);
          if (marker[i]) {
            if (trtGrp[i] == 1) survivalT[i] = qtpwexpcpp1(u, tau1pos[j], lam1posy[j]);
            else survivalT[i] = qtpwexpcpp1(u, tau2pos[j], lam2posy[j]);
          } else {
            if (trtGrp[i] == 1) survivalT[i] = qtpwexpcpp1(u, tau1neg[j], lam1negy[j]);
            else survivalT[i] = qtpwexpcpp1(u, tau2neg[j], lam2negy[j]);
          }

          // dropout time
          u = unif(rng_local);
          if (marker[i]) {
            if (trtGrp[i] == 1) dropoutT[i] = qtpwexpcpp1(u, tau1pos[j], gam1posy[j]);
            else dropoutT[i] = qtpwexpcpp1(u, tau2pos[j], gam2posy[j]);
          } else {
            if (trtGrp[i] == 1) dropoutT[i] = qtpwexpcpp1(u, tau1neg[j], gam1negy[j]);
            else dropoutT[i] = qtpwexpcpp1(u, tau2neg[j], gam2negy[j]);
          }

          // observed time and event
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
        } // cohort

        // determine analysis times
        size_t nstages = K;
        bool ev1NotAch = false, ev2NotAch = false;

        if (useEvents) {
          totalte.clear(); totaltepos.clear();
          int nevents = 0, neventspos = 0;
          for (size_t i = 0; i < N; ++i) {
            if (event[i]) {
              ++nevents; totalte.push_back(totalT[i]);
            }
            if (event[i] && marker[i]) {
              ++neventspos; totaltepos.push_back(totalT[i]);
            }
          }
          if (nevents == 0 || neventspos == 0) {
            thread_utils::push_thread_warning(
              std::string("No events for iteration ") + std::to_string(iter+1) +
                " skipping this iteration.");
            out.summaryRows.clear();
            out.rawRows.clear();
            continue;
          }
          std::sort(totalte.begin(), totalte.end());
          std::sort(totaltepos.begin(), totaltepos.end());

          // ITT looks
          analysisT1.clear();
          size_t j1 = 0;
          if (Kitt > 0) {
            for (j1 = 0; j1 < Kitt; ++j1) {
              if (plannedEvents[j1] >= nevents) break;
            }
            if (j1 == Kitt) {
              for (size_t k = 0; k < Kitt; ++k) {
                analysisT1.push_back(totalte[plannedEvents[k] - 1] + 1e-12);
              }
            } else {
              for (size_t k = 0; k < j1; ++k) {
                analysisT1.push_back(totalte[plannedEvents[k] - 1] + 1e-12);
              }
              analysisT1.push_back(totalte.back() + 1e-12);
            }
          }

          // biomarker+ looks
          analysisT2.clear();
          size_t j2 = 0;
          if (K > Kitt) {
            for (j2 = 0; j2 < (K - Kitt); ++j2) {
              if (plannedEvents[Kitt + j2] >= neventspos) break;
            }
            if (j2 == (K - Kitt)) {
              for (size_t k = 0; k < (K - Kitt); ++k) {
                analysisT2.push_back(totaltepos[plannedEvents[Kitt + k] - 1] + 1e-12);
              }
            } else {
              for (size_t k = 0; k < j2; ++k) {
                analysisT2.push_back(totaltepos[plannedEvents[Kitt + k] - 1] + 1e-12);
              }
              analysisT2.push_back(totaltepos.back() + 1e-12);
            }
          }

          // combine
          if (Kitt == 0) {
            nstages = analysisT2.size();
            for (size_t k = 0; k < nstages; ++k) analysisT[k] = analysisT2[k];
          } else if (K == Kitt) {
            nstages = analysisT1.size();
            for (size_t k = 0; k < nstages; ++k) analysisT[k] = analysisT1[k];
          } else {
            if (analysisT2.back() > analysisT1.back()) {
              // NOTE: In this case, the observed number of ITT events must exceed
              // the planned number of ITT events at look Kitt, because otherwise
              // the last ITT event would be observed at analysisTime1.back().
              // However, since the last biomarker+ event occurred on or after
              // analysisTime2.back() > analysisTime1.back(), this is a
              // contradiction as biomarker+ event is part of ITT event.
              // It follows that analysisTime1.size() == Kitt in this case.

              // find first biomarker+ look after last ITT look
              size_t l = 0;
              for (size_t idx = 0; idx < analysisT2.size(); ++idx) {
                if (analysisT2[idx] > analysisT1.back()) { l = idx; break; }
              }
              nstages = Kitt + (analysisT2.size() - l);
              std::copy_n(analysisT1.begin(), Kitt, analysisT.begin());
              size_t count = analysisT2.size() - l;
              std::copy_n(analysisT2.begin() + l, count, analysisT.begin() + Kitt);
            } else {
              nstages = analysisT1.size();
              std::copy_n(analysisT1.begin(), nstages, analysisT.begin());
            }
          }

          if (Kitt > 0 && nevents < plannedEvents[Kitt - 1]) ev1NotAch = true;
          if (Kitt < K && neventspos < plannedEvents[K - 1]) ev2NotAch = true;
        } else {
          nstages = K;
          std::copy_n(plannedTime.begin(), K, analysisT.begin());
        }

        // per-stage computations
        for (size_t k = 0; k < nstages; ++k) {
          double time = analysisT[k];

          for (size_t i = 0; i < N; ++i) {
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
          } // censoring loop

          // optionally append raw rows for this stage
          if (iter < maxRawIters) {
            for (size_t i = 0; i < N; ++i) {
              // skip subjects who haven't been enrolled by analysis time
              if (arrivalT[i] > time) continue;
              RawDatasetRow rr;
              rr.iterNum = static_cast<int>(iter + 1);
              rr.stageNum = static_cast<int>(k + 1);
              rr.analysisT = time;
              rr.subjectId = static_cast<int>(i + 1);
              rr.arrivalT = arrivalT[i];
              rr.stratum = stratum[i];
              rr.marker = marker[i];
              rr.trtGrp = trtGrp[i];
              rr.survivalT = survivalT[i];
              rr.dropoutT = dropoutT[i];
              rr.timeObs = timeObs[i];
              rr.event = event[i];
              rr.dropEv = dropEv[i];
              out.rawRows.push_back(std::move(rr));
            }
          }

          // three populations: 1=ITT, 2=Biomarker+, 3=Biomarker-
          for (int pop = 1; pop <= 3; ++pop) {
            double hazardRatioH0;
            set.clear();
            if (pop == 1) {
              hazardRatioH0 = hazardRatioH0itt;
              for (size_t i = 0; i < N; ++i) set.push_back(i);
            } else if (pop == 2) {
              hazardRatioH0 = hazardRatioH0pos;
              for (size_t i = 0; i < N; ++i) if (marker[i]) set.push_back(i);
            } else {
              hazardRatioH0 = hazardRatioH0neg;
              for (size_t i = 0; i < N; ++i) if (!marker[i]) set.push_back(i);
            }

            // reset risk sets per stratum
            std::fill(n1.begin(), n1.end(), 0);
            std::fill(n2.begin(), n2.end(), 0);
            int events1 = 0, events2 = 0, dropouts1 = 0, dropouts2 = 0;

            for (size_t i = 0; i < set.size(); ++i) {
              size_t idx = set[i];
              if (arrivalT[idx] > time) continue;
              size_t h = static_cast<size_t>(stratum[idx] - 1);
              if (trtGrp[idx] == 1) {
                ++n1[h];
                if (event[idx]) ++events1; else if (dropEv[idx]) ++dropouts1;
              } else {
                ++n2[h];
                if (event[idx]) ++events2; else if (dropEv[idx]) ++dropouts2;
              }
            }

            int accruals1 = std::accumulate(n1.begin(), n1.end(), 0);
            int accruals2 = std::accumulate(n2.begin(), n2.end(), 0);
            int totAccruals = accruals1 + accruals2;
            int totEvents = events1 + events2;
            int totDropouts = dropouts1 + dropouts2;

            // build list of indices with positive observed time
            sub.clear();
            for (size_t i = 0; i < set.size(); ++i) {
              size_t idx = set[i];
              if (timeObs[idx] > 0.0) sub.push_back(idx);
            }
            // sort by observed time ascending
            std::sort(sub.begin(), sub.end(), [&](size_t a, size_t b) {
              return timeObs[a] < timeObs[b];
            });

            // stratified log-rank for this population
            std::fill(km.begin(), km.end(), 1.0);
            double us = 0.0, vs = 0.0;

            for (size_t i = 0; i < sub.size(); ++i) {
              size_t idx = sub[i];
              size_t h = static_cast<size_t>(stratum[idx] - 1);
              double n1h = static_cast<double>(n1[h]);
              double n2h = static_cast<double>(n2[h]);
              double n1a = n1h * hazardRatioH0;
              double nt = n1h + n2h;
              double nta = n1a + n2h;

              if (event[idx]) {
                double wh = 1.0;
                if (rho1 != 0.0 || rho2 != 0.0) {
                  wh = std::pow(km[h], rho1) * std::pow(1.0 - km[h], rho2);
                  km[h] *= (1.0 - 1.0 / nt);
                }
                double treated = (trtGrp[idx] == 1 ? 1.0 : 0.0);
                us += wh * (treated - n1a / nta);
                vs += wh * wh * n1a * n2h / (nta * nta);
              }

              // reduce risk set
              if (trtGrp[idx] == 1) --n1[h]; else --n2[h];
            } // events loop

            double z = (vs > 0.0 ? us / std::sqrt(vs) : 0.0);

            // append summary row
            StageSummaryRow sr;
            sr.iterNum = static_cast<int>(iter + 1);
            sr.ev1NotAch = ev1NotAch ? 1 : 0;
            sr.ev2NotAch = ev2NotAch ? 1 : 0;
            sr.stageNum = static_cast<int>(k + 1);
            sr.analysisT = time;
            sr.pop = pop;
            sr.accruals1 = accruals1;
            sr.accruals2 = accruals2;
            sr.totAccruals = totAccruals;
            sr.events1 = events1;
            sr.events2 = events2;
            sr.totEvents = totEvents;
            sr.dropouts1 = dropouts1;
            sr.dropouts2 = dropouts2;
            sr.totDropouts = totDropouts;
            sr.uscore = us;
            sr.vscore = vs;
            sr.logRank = z;
            out.summaryRows.push_back(std::move(sr));
          } // populations loop
        } // per-stage
      } // iter
    } // operator()
  }; // SimWorker

  // construct and run worker
  SimWorker worker(
      K, Kitt, hazardRatioH0itt, hazardRatioH0pos, hazardRatioH0neg,
      allocation1, allocation2,
      accrualTime, accrualIntensity, tau, stratumFraction, p_posv,
      lam1ittx, lam2ittx, lam1posx, lam2posx,
      gam1ittx, gam2ittx, gam1posx, gam2posx,
      tau1pos, tau2pos, tau1neg, tau2neg,
      lam1posy, lam2posy, lam1negy, lam2negy,
      gam1posy, gam2posy, gam1negy, gam2negy,
      N, fu, fixedFollowup, rho1, rho2,
      plannedEvents, plannedTime,
      maxIters, maxRawIters, seeds, useEvents, nstrata,
      &results
  );

  RcppParallel::parallelFor(0, maxIters, worker);

  // Flatten results
  size_t nsr = 0, nrr = 0;
  for (size_t iter = 0; iter < maxIters; ++iter) {
    nsr += results[iter].summaryRows.size();
    nrr += results[iter].rawRows.size();
  }
  if (nsr == 0) throw std::runtime_error(
    "No iterations with observed events. Unable to produce output.");

  // Prepare final containers
  std::vector<int> sum_iterNum; sum_iterNum.reserve(nsr);
  std::vector<unsigned char> sum_ev1NotAch; sum_ev1NotAch.reserve(nsr);
  std::vector<unsigned char> sum_ev2NotAch; sum_ev2NotAch.reserve(nsr);
  std::vector<int> sum_stageNum; sum_stageNum.reserve(nsr);
  std::vector<double> sum_analysisT; sum_analysisT.reserve(nsr);
  std::vector<std::string> sum_pop; sum_pop.reserve(nsr);
  std::vector<int> sum_accruals1; sum_accruals1.reserve(nsr);
  std::vector<int> sum_accruals2; sum_accruals2.reserve(nsr);
  std::vector<int> sum_totAccruals; sum_totAccruals.reserve(nsr);
  std::vector<int> sum_events1; sum_events1.reserve(nsr);
  std::vector<int> sum_events2; sum_events2.reserve(nsr);
  std::vector<int> sum_totEvents; sum_totEvents.reserve(nsr);
  std::vector<int> sum_dropouts1; sum_dropouts1.reserve(nsr);
  std::vector<int> sum_dropouts2; sum_dropouts2.reserve(nsr);
  std::vector<int> sum_totDropouts; sum_totDropouts.reserve(nsr);
  std::vector<double> sum_uscore; sum_uscore.reserve(nsr);
  std::vector<double> sum_vscore; sum_vscore.reserve(nsr);
  std::vector<double> sum_logRank; sum_logRank.reserve(nsr);

  // raw final containers
  std::vector<int> raw_iterNum; raw_iterNum.reserve(nrr);
  std::vector<int> raw_stageNum; raw_stageNum.reserve(nrr);
  std::vector<double> raw_analysisT; raw_analysisT.reserve(nrr);
  std::vector<int> raw_subjectId; raw_subjectId.reserve(nrr);
  std::vector<double> raw_arrivalT; raw_arrivalT.reserve(nrr);
  std::vector<int> raw_stratum; raw_stratum.reserve(nrr);
  std::vector<unsigned char> raw_marker; raw_marker.reserve(nrr);
  std::vector<int> raw_trtGrp; raw_trtGrp.reserve(nrr);
  std::vector<double> raw_survivalT; raw_survivalT.reserve(nrr);
  std::vector<double> raw_dropoutT; raw_dropoutT.reserve(nrr);
  std::vector<double> raw_timeObs; raw_timeObs.reserve(nrr);
  std::vector<unsigned char> raw_event; raw_event.reserve(nrr);
  std::vector<unsigned char> raw_dropEv; raw_dropEv.reserve(nrr);

  for (size_t iter = 0; iter < maxIters; ++iter) {
    const auto& srows = results[iter].summaryRows;
    for (const auto& r : srows) {
      sum_iterNum.push_back(r.iterNum);
      sum_ev1NotAch.push_back(r.ev1NotAch);
      sum_ev2NotAch.push_back(r.ev2NotAch);
      sum_stageNum.push_back(r.stageNum);
      sum_analysisT.push_back(r.analysisT);
      sum_pop.push_back(r.pop == 1 ? "ITT" : (r.pop == 2 ? "Biomarker+" :
                                                "Biomarker-"));
      sum_accruals1.push_back(r.accruals1);
      sum_accruals2.push_back(r.accruals2);
      sum_totAccruals.push_back(r.totAccruals);
      sum_events1.push_back(r.events1);
      sum_events2.push_back(r.events2);
      sum_totEvents.push_back(r.totEvents);
      sum_dropouts1.push_back(r.dropouts1);
      sum_dropouts2.push_back(r.dropouts2);
      sum_totDropouts.push_back(r.totDropouts);
      sum_uscore.push_back(r.uscore);
      sum_vscore.push_back(r.vscore);
      sum_logRank.push_back(r.logRank);
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
        raw_marker.push_back(rr.marker);
        raw_trtGrp.push_back(rr.trtGrp);
        raw_survivalT.push_back(rr.survivalT);
        raw_dropoutT.push_back(rr.dropoutT);
        raw_timeObs.push_back(rr.timeObs);
        raw_event.push_back(rr.event);
        raw_dropEv.push_back(rr.dropEv);
      }
    }
  }

  // Build DataFrameCpp summary
  DataFrameCpp sumdata;
  sumdata.push_back(std::move(sum_iterNum), "iterNumber");
  sumdata.push_back(std::move(sum_ev1NotAch), "events1NotAchieved");
  sumdata.push_back(std::move(sum_ev2NotAch), "events2NotAchieved");
  sumdata.push_back(std::move(sum_stageNum), "stageNumber");
  sumdata.push_back(std::move(sum_analysisT), "analysisTime");
  sumdata.push_back(std::move(sum_pop), "population");
  sumdata.push_back(std::move(sum_accruals1), "accruals1");
  sumdata.push_back(std::move(sum_accruals2), "accruals2");
  sumdata.push_back(std::move(sum_totAccruals), "totalAccruals");
  sumdata.push_back(std::move(sum_events1), "events1");
  sumdata.push_back(std::move(sum_events2), "events2");
  sumdata.push_back(std::move(sum_totEvents), "totalEvents");
  sumdata.push_back(std::move(sum_dropouts1), "dropouts1");
  sumdata.push_back(std::move(sum_dropouts2), "dropouts2");
  sumdata.push_back(std::move(sum_totDropouts), "totalDropouts");
  sumdata.push_back(std::move(sum_uscore), "uscore");
  sumdata.push_back(std::move(sum_vscore), "vscore");
  sumdata.push_back(std::move(sum_logRank), "logRankStatistic");

  ListCpp result;
  result.push_back(sumdata, "sumdata");

  if (!raw_iterNum.empty()) {
    DataFrameCpp rawdata;
    rawdata.push_back(std::move(raw_iterNum), "iterationNumber");
    rawdata.push_back(std::move(raw_stageNum), "stageNumber");
    rawdata.push_back(std::move(raw_analysisT), "analysisTime");
    rawdata.push_back(std::move(raw_subjectId), "subjectId");
    rawdata.push_back(std::move(raw_arrivalT), "arrivalTime");
    rawdata.push_back(std::move(raw_stratum), "stratum");
    rawdata.push_back(std::move(raw_marker), "biomarker");
    rawdata.push_back(std::move(raw_trtGrp), "treatmentGroup");
    rawdata.push_back(std::move(raw_survivalT), "survivalTime");
    rawdata.push_back(std::move(raw_dropoutT), "dropoutTime");
    rawdata.push_back(std::move(raw_timeObs), "timeUnderObservation");
    rawdata.push_back(std::move(raw_event), "event");
    rawdata.push_back(std::move(raw_dropEv), "dropoutEvent");
    result.push_back(rawdata, "rawdata");
  }

  return result;
}


//' @title Log-Rank Test Simulation for Enrichment Design
//' @description Performs simulation for two-arm group
//' sequential trials based on weighted log-rank test
//' for a biomarker enrichment design. The looks are either
//' driven by the total number of events in the ITT population
//' or the biomarker positive sub population.
//' Alternatively, the analyses can be planned to occur at
//' specified calendar times.
//'
//' @inheritParams param_kMax
//' @param kMaxitt Number of stages with timing determined by events
//'   in the ITT population. Ranges from 0 (none) to \code{kMax}.
//' @param hazardRatioH0itt Hazard ratio under the null hypothesis
//'   for the ITT population. Defaults to 1 for superiority test.
//' @param hazardRatioH0pos Hazard ratio under the null hypothesis
//'   for the biomarker positive sub population. Defaults to 1 for
//'   superiority test.
//' @param hazardRatioH0neg Hazard ratio under the null hypothesis
//'   for the biomarker negative sub population. Defaults to 1 for
//'   superiority test.
//' @param allocation1 Number of subjects in the treatment group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in the control group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param p_pos The prevalence of the biomarker positive sub population
//'   in each stratum.
//' @param lambda1itt A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the treatment group in the ITT population.
//' @param lambda2itt A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the control group in the ITT population.
//' @param lambda1pos A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the treatment group in the biomarker
//'   positive sub population.
//' @param lambda2pos A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the control group in the biomarker
//'   positive sub population.
//' @param gamma1itt The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the treatment group in the ITT population.
//' @param gamma2itt The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the control group in the ITT population.
//' @param gamma1pos The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the treatment group in the biomarker
//'   positive sub population.
//' @param gamma2pos The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the control group in the biomarker
//'   positive sub population.
//' @param n Sample size.
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @param plannedEvents The planned cumulative total number events in the
//'   ITT population at Look 1 to Look \code{kMaxitt} and the planned
//'   cumulative total number of events at Look \code{kMaxitt+1} to
//'   Look \code{kMax} in the biomarker positive sub population.
//' @param plannedTime The calendar times for the analyses. To use calendar
//'   time to plan the analyses, \code{plannedEvents} should be missing.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
//'   stage to extract.
//' @param seed The seed to reproduce the simulation results.
//'
//' @return A list with 2 components:
//'
//' * \code{sumdata}: A data frame of summary data by iteration and stage:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{eventsNotAchieved}: Whether the target number of events
//'       is not achieved for the iteration.
//'
//'     - \code{stageNumber}: The stage number, covering all stages even if
//'       the trial stops at an interim look.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{population}: The population ("ITT", "Biomarker Positive",
//'       "Biomarker Negative") under consideration.
//'
//'     - \code{accruals1}: The number of subjects enrolled at the stage for
//'       the treatment group.
//'
//'     - \code{accruals2}: The number of subjects enrolled at the stage for
//'       the control group.
//'
//'     - \code{totalAccruals}: The total number of subjects enrolled at
//'       the stage.
//'
//'     - \code{events1}: The number of events at the stage for
//'       the treatment group.
//'
//'     - \code{events2}: The number of events at the stage for
//'       the control group.
//'
//'     - \code{totalEvents}: The total number of events at the stage.
//'
//'     - \code{dropouts1}: The number of dropouts at the stage for
//'       the treatment group.
//'
//'     - \code{dropouts2}: The number of dropouts at the stage for
//'       the control group.
//'
//'     - \code{totalDropouts}: The total number of dropouts at the stage.
//'
//'     - \code{logRankStatistic}: The log-rank test Z-statistic for
//'       the population.
//'
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//'   positive integer): A data frame for subject-level data for selected
//'   replications, containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage under consideration.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{subjectId}: The subject ID.
//'
//'     - \code{arrivalTime}: The enrollment time for the subject.
//'
//'     - \code{stratum}: The stratum for the subject.
//'
//'     - \code{biomarker}: The biomarker status for the subject (1 for
//'       positive, 0 for negative).
//'
//'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
//'       subject.
//'
//'     - \code{survivalTime}: The underlying survival time for the subject.
//'
//'     - \code{dropoutTime}: The underlying dropout time for the subject.
//'
//'     - \code{timeUnderObservation}: The time under observation
//'       since randomization for the subject.
//'
//'     - \code{event}: Whether the subject experienced an event.
//'
//'     - \code{dropoutEvent}: Whether the subject dropped out.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' sim1 = lrsimsub(
//'   kMax = 2,
//'   kMaxitt = 2,
//'   allocation1 = 1,
//'   allocation2 = 1,
//'   accrualTime = seq(0,9),
//'   accrualIntensity = c(seq(10,70,10),rep(70,3)),
//'   piecewiseSurvivalTime = c(0,12,24),
//'   p_pos = 0.6,
//'   lambda1itt = c(0.00256, 0.00383, 0.00700),
//'   lambda2itt = c(0.00427, 0.00638, 0.01167),
//'   lambda1pos = c(0.00299, 0.00430, 0.01064),
//'   lambda2pos = c(0.00516, 0.00741, 0.01835),
//'   gamma1itt = -log(1-0.04)/12,
//'   gamma2itt = -log(1-0.04)/12,
//'   gamma1pos = -log(1-0.04)/12,
//'   gamma2pos = -log(1-0.04)/12,
//'   n = 500,
//'   plannedEvents = c(108,144),
//'   maxNumberOfIterations = 1000,
//'   maxNumberOfRawDatasetsPerStage = 1,
//'   seed = 314159)
//'
//' head(sim1$sumdata)
//' head(sim1$rawdata)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List lrsimsub(
    const int kMax = 1,
    const int kMaxitt = 1,
    const double hazardRatioH0itt = 1,
    const double hazardRatioH0pos = 1,
    const double hazardRatioH0neg = 1,
    const int allocation1 = 1,
    const int allocation2 = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& p_pos = NA_REAL,
    const Rcpp::NumericVector& lambda1itt = NA_REAL,
    const Rcpp::NumericVector& lambda2itt = NA_REAL,
    const Rcpp::NumericVector& lambda1pos = NA_REAL,
    const Rcpp::NumericVector& lambda2pos = NA_REAL,
    const Rcpp::NumericVector& gamma1itt = 0,
    const Rcpp::NumericVector& gamma2itt = 0,
    const Rcpp::NumericVector& gamma1pos = 0,
    const Rcpp::NumericVector& gamma2pos = 0,
    const int n = NA_INTEGER,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = 0,
    const double rho1 = 0,
    const double rho2 = 0,
    const Rcpp::IntegerVector& plannedEvents = NA_INTEGER,
    const Rcpp::NumericVector& plannedTime = NA_REAL,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = 0) {

  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvTime = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto pPos = Rcpp::as<std::vector<double>>(p_pos);
  auto lam1itt = Rcpp::as<std::vector<double>>(lambda1itt);
  auto lam2itt = Rcpp::as<std::vector<double>>(lambda2itt);
  auto lam1pos = Rcpp::as<std::vector<double>>(lambda1pos);
  auto lam2pos = Rcpp::as<std::vector<double>>(lambda2pos);
  auto gam1itt = Rcpp::as<std::vector<double>>(gamma1itt);
  auto gam2itt = Rcpp::as<std::vector<double>>(gamma2itt);
  auto gam1pos = Rcpp::as<std::vector<double>>(gamma1pos);
  auto gam2pos = Rcpp::as<std::vector<double>>(gamma2pos);
  auto plannedE = Rcpp::as<std::vector<int>>(plannedEvents);
  auto plannedT = Rcpp::as<std::vector<double>>(plannedTime);

  auto out = lrsimsubcpp(
    kMax, kMaxitt, hazardRatioH0itt, hazardRatioH0pos,
    hazardRatioH0neg, allocation1, allocation2,
    accrualT, accrualInt, pwSurvTime, stratumFrac, pPos,
    lam1itt, lam2itt, lam1pos, lam2pos,
    gam1itt, gam2itt, gam1pos, gam2pos,
    n, followupTime, fixedFollowup, rho1, rho2,
    plannedE, plannedT, maxNumberOfIterations,
    maxNumberOfRawDatasetsPerStage, seed);

  thread_utils::drain_thread_warnings_to_R();

  return Rcpp::wrap(out);
}


ListCpp binary_tte_sim_cpp(
    const int kMax1,
    const int kMax2,
    const double riskDiffH0,
    const double hazardRatioH0,
    const int allocation1,
    const int allocation2,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const double globalOddsRatio,
    const std::vector<double>& pi1,
    const std::vector<double>& pi2,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const std::vector<double>& delta1,
    const std::vector<double>& delta2,
    const double upper1,
    const double upper2,
    const int n,
    const std::vector<double>& plannedTime,
    const std::vector<int>& plannedEvents,
    const int maxNumberOfIterations,
    const int maxNumberOfRawDatasetsPerStage,
    const int seed)
{
  if (kMax1 < 1) throw std::invalid_argument("kMax1 must be a positive integer");
  if (kMax2 < 1) throw std::invalid_argument("kMax2 must be a positive integer");
  size_t K1 = static_cast<size_t>(kMax1);
  size_t K2 = static_cast<size_t>(kMax2);

  if (!(riskDiffH0 > -1.0 && riskDiffH0 < 1.0))
    throw std::invalid_argument("riskDiffH0 must lie between -1 and 1");
  if (hazardRatioH0 <= 0.0)
    throw std::invalid_argument("hazardRatioH0 must be positive");
  if (allocation1 < 1 || allocation2 < 1)
    throw std::invalid_argument("allocations must be positive integers");
  if (accrualTime[0] != 0.0)
    throw std::invalid_argument("accrualTime must start with 0");
  if (any_nonincreasing(accrualTime))
    throw std::invalid_argument("accrualTime should be increasing");
  if (!none_na(accrualIntensity))
    throw std::invalid_argument("accrualIntensity must be provided");
  if (accrualTime.size() != accrualIntensity.size())
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
  if (globalOddsRatio <= 0.0)
    throw std::invalid_argument("globalOddsRatio must be positive");
  if (!none_na(pi1)) throw std::invalid_argument("pi1 must be provided");
  if (!none_na(pi2)) throw std::invalid_argument("pi2 must be provided");
  for (double v : pi1) {
    if (!(v > 0.0 && v < 1.0))
      throw std::invalid_argument("pi1 must lie between 0 and 1");
  }
  for (double v : pi2) {
    if (!(v > 0.0 && v < 1.0))
      throw std::invalid_argument("pi2 must lie between 0 and 1");
  }
  if (!none_na(lambda1)) throw std::invalid_argument("lambda1 must be provided");
  if (!none_na(lambda2)) throw std::invalid_argument("lambda2 must be provided");
  for (double v : lambda1) {
    if (v < 0.0) throw std::invalid_argument("lambda1 must be non-negative");
  }
  for (double v : lambda2) {
    if (v < 0.0) throw std::invalid_argument("lambda2 must be non-negative");
  }
  for (double v : gamma1) {
    if (v < 0.0) throw std::invalid_argument("gamma1 must be non-negative");
  }
  for (double v : gamma2) {
    if (v < 0.0) throw std::invalid_argument("gamma2 must be non-negative");
  }
  for (double v : delta1) {
    if (v < 0.0) throw std::invalid_argument("delta1 must be non-negative");
  }
  for (double v : delta2) {
    if (v < 0.0) throw std::invalid_argument("delta2 must be non-negative");
  }
  if (upper1 <= 0.0) throw std::invalid_argument("upper1 must be positive");
  if (upper2 <= 0.0) throw std::invalid_argument("upper2 must be positive");
  if (n == INT_MIN) throw std::invalid_argument("n must be provided");
  if (n <= 0) throw std::invalid_argument("n must be positive");
  if (!none_na(plannedTime))
    throw std::invalid_argument("plannedTime must be given for endpoint 1");
  if (plannedTime[0] <= 0.0)
    throw std::invalid_argument("plannedTime must be positive");
  if (plannedTime.size() != K1)
    throw std::invalid_argument("Invalid length for plannedTime");
  if (any_nonincreasing(plannedTime))
    throw std::invalid_argument("plannedTime must be increasing");
  if (!none_na(plannedEvents))
    throw std::invalid_argument("plannedEvents must be given for endpoint 2");
  if (plannedEvents[0] <= 0)
    throw std::invalid_argument("plannedEvents must be positive");
  if (plannedEvents.size() != K2)
    throw std::invalid_argument("Invalid length for plannedEvents");
  if (any_nonincreasing(plannedEvents))
    throw std::invalid_argument("plannedEvents must be increasing");
  if (maxNumberOfIterations < 1)
    throw std::invalid_argument("maxNumberOfIterations must be a positive integer");
  if (maxNumberOfRawDatasetsPerStage < 0)
    throw std::invalid_argument(
        "maxNumberOfRawDatasetsPerStage must be a non-negative integer");
  if (maxNumberOfRawDatasetsPerStage > maxNumberOfIterations)
    throw std::invalid_argument(
        "maxNumberOfRawDatasetsPerStage cannot exceed maxNumberOfIterations");

  size_t N = static_cast<size_t>(n);
  size_t maxIters = static_cast<size_t>(maxNumberOfIterations);
  size_t maxRawIters = static_cast<size_t>(maxNumberOfRawDatasetsPerStage);
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  const std::vector<double>& tau = piecewiseSurvivalTime;

  // expand per-stratum inputs on main thread
  auto pi1v = expand1(pi1, nstrata, "pi1");
  auto pi2v = expand1(pi2, nstrata, "pi2");
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto gamma1x  = expand_stratified(gamma1,  nstrata, nintv, "gamma1");
  auto gamma2x  = expand_stratified(gamma2,  nstrata, nintv, "gamma2");
  auto delta1x  = expand_stratified(delta1,  nstrata, nintv, "delta1");
  auto delta2x  = expand_stratified(delta2,  nstrata, nintv, "delta2");

  // per-iteration result struct and raw row struct defined inside the function
  // endpoint 1 (binary)
  struct StageSummary1Row {
    int iterNum = 0, stageNum = 0;
    double analysisT = 0.0;
    int accruals1 = 0, accruals2 = 0, totAccruals = 0;
    int source1 = 0, source2 = 0, source3 = 0;
    double n1 = 0.0, n2 = 0.0, n = 0.0;
    double y1 = 0.0, y2 = 0.0, y = 0.0;
    double riskDiff = 0.0, seRiskDiff = 0.0, z = 0.0;
  };

  // endpoint 2 (TTE) fields (unused for endpoint1 rows)
  struct StageSummary2Row {
    int iterNum = 0, stageNum = 0;
    unsigned char evNotAch = 0;
    double analysisT = 0.0;
    int accruals1 = 0, accruals2 = 0, totAccruals = 0;
    int events1 = 0, events2 = 0, totEvents = 0;
    int dropouts1 = 0, dropouts2 = 0, totDropouts = 0;
    double uscore = 0.0, vscore = 0.0, logRank = 0.0;
  };


  // endpoint1 raw
  struct RawDataset1Row {
    int iterNum = 0, stageNum = 0;
    double analysisT = 0.0;
    int subjectId = 0;
    double arrivalT = 0.0;
    int stratum = 0, trtGrp = 0;
    double survivalT = 0.0, dropoutT = 0.0;
    double trtDiscT = 0.0, upper = 0.0, ptfu1T = 0.0;
    double timeObs = 0.0;
    double latentResp = 0.0;
    unsigned char responder = 0;
    int source = 0;
  };

  // endpoint2 raw
  struct RawDataset2Row {
    int iterNum = 0, stageNum = 0;
    double analysisT = 0.0;
    int subjectId = 0;
    double arrivalT = 0.0;
    int stratum = 0, trtGrp = 0;
    double survivalT = 0.0, dropoutT = 0.0, timeObs = 0.0;
    unsigned char event = 0, dropEv = 0;
  };

  struct IterationResult {
    std::vector<StageSummary1Row> summary1Rows;
    std::vector<StageSummary2Row> summary2Rows;
    std::vector<RawDataset1Row> raw1Rows;
    std::vector<RawDataset2Row> raw2Rows;
    void reserveForSummary1(size_t approxRows) { summary1Rows.reserve(approxRows); }
    void reserveForSummary2(size_t approxRows) { summary2Rows.reserve(approxRows); }
    void reserveForRaw1(size_t approxRows) { raw1Rows.reserve(approxRows); }
    void reserveForRaw2(size_t approxRows) { raw2Rows.reserve(approxRows); }
  };

  // Pre-size results
  std::vector<IterationResult> results;
  results.resize(maxIters);

  // Precompute logistic intercepts alpha per stratum for latent response sampling
  std::vector<double> alpha1v(nstrata), alpha2v(nstrata);
  for (size_t s = 0; s < nstrata; ++s) {
    alpha1v[s] = std::log(pi1v[s] / (1.0 - pi1v[s]));
    alpha2v[s] = std::log(pi2v[s] / (1.0 - pi2v[s]));
  }

  double globalOddsRatio1 = globalOddsRatio - 1.0;

  auto f = [&](double v, double u, double t)->double {
    double c1 = 1.0 + globalOddsRatio1 * (u + v);
    double c2 = 4.0 * v * globalOddsRatio * globalOddsRatio1;
    double c3 = 2.0 * c1 * globalOddsRatio1 - c2;
    double sqrtdisc = std::sqrt(c1 * c1 - c2);
    double numerator = globalOddsRatio1 - 0.5 / sqrtdisc * c3;
    double denominator = 2.0 * globalOddsRatio1;
    return numerator / denominator - t;
  };


  // seeds
  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter) seeds[iter] = master_rng();

  // Worker
  struct SimWorker : public RcppParallel::Worker {
    // inputs (const refs)
    const size_t K1;
    const size_t K2;
    const double riskDiffH0;
    const double hazardRatioH0;
    const int allocation1, allocation2;
    const std::vector<double>& accrualTime;
    const std::vector<double>& accrualIntensity;
    const std::vector<double>& tau;
    const std::vector<double>& stratumFraction;
    const double globalOddsRatio;
    const std::vector<double>& pi1v;
    const std::vector<double>& pi2v;
    const std::vector<std::vector<double>>& lambda1x;
    const std::vector<std::vector<double>>& lambda2x;
    const std::vector<std::vector<double>>& gamma1x;
    const std::vector<std::vector<double>>& gamma2x;
    const std::vector<std::vector<double>>& delta1x;
    const std::vector<std::vector<double>>& delta2x;
    const double upper1;
    const double upper2;
    const size_t N;
    const std::vector<double>& plannedTime;
    const std::vector<int>& plannedEvents;
    const size_t maxIters;
    const size_t maxRawIters;
    const std::vector<uint64_t>& seeds;
    const std::vector<double>& alpha1v;
    const std::vector<double>& alpha2v;
    const double globalOddsRatio1;
    const size_t nstrata;

    std::function<double(const double, const double, const double)> f;

    std::vector<IterationResult>* results;

    SimWorker(
      size_t K1_,
      size_t K2_,
      double riskDiffH0_,
      double hazardRatioH0_,
      int allocation1_,
      int allocation2_,
      const std::vector<double>& accrualTime_,
      const std::vector<double>& accrualIntensity_,
      const std::vector<double>& tau_,
      const std::vector<double>& stratumFraction_,
      double globalOddsRatio_,
      const std::vector<double>& pi1v_,
      const std::vector<double>& pi2v_,
      const std::vector<std::vector<double>>& lambda1x_,
      const std::vector<std::vector<double>>& lambda2x_,
      const std::vector<std::vector<double>>& gamma1x_,
      const std::vector<std::vector<double>>& gamma2x_,
      const std::vector<std::vector<double>>& delta1x_,
      const std::vector<std::vector<double>>& delta2x_,
      double upper1_,
      double upper2_,
      size_t N_,
      const std::vector<double>& plannedTime_,
      const std::vector<int>& plannedEvents_,
      size_t maxIters_,
      size_t maxRawIters_,
      const std::vector<uint64_t>& seeds_,
      const std::vector<double>& alpha1v_,
      const std::vector<double>& alpha2v_,
      double globalOddsRatio1_,
      size_t nstrata_,
      decltype(f) f_,
      std::vector<IterationResult>* results_)
      : K1(K1_),
        K2(K2_),
        riskDiffH0(riskDiffH0_),
        hazardRatioH0(hazardRatioH0_),
        allocation1(allocation1_),
        allocation2(allocation2_),
        accrualTime(accrualTime_),
        accrualIntensity(accrualIntensity_),
        tau(tau_),
        stratumFraction(stratumFraction_),
        globalOddsRatio(globalOddsRatio_),
        pi1v(pi1v_),
        pi2v(pi2v_),
        lambda1x(lambda1x_),
        lambda2x(lambda2x_),
        gamma1x(gamma1x_),
        gamma2x(gamma2x_),
        delta1x(delta1x_),
        delta2x(delta2x_),
        upper1(upper1_),
        upper2(upper2_),
        N(N_),
        plannedTime(plannedTime_),
        plannedEvents(plannedEvents_),
        maxIters(maxIters_),
        maxRawIters(maxRawIters_),
        seeds(seeds_),
        alpha1v(alpha1v_),
        alpha2v(alpha2v_),
        globalOddsRatio1(globalOddsRatio1_),
        nstrata(nstrata_),
        f(std::move(f_)),
        results(results_)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      // local buffers per worker
      std::vector<int> stratum(N), trtGrp(N);
      std::vector<double> arrivalT(N), survivalT(N), dropoutT(N);
      std::vector<double> timeObs(N), totalT(N);
      std::vector<unsigned char> event(N), dropEv(N);

      std::vector<double> latentResp(N);
      std::vector<double> trtDiscT(N), ptfu1T(N), timeObs1(N);
      std::vector<unsigned char> responder(N);
      std::vector<int> source(N);

      std::vector<int> b1(nstrata), b2(nstrata);
      std::vector<int> n1(nstrata), n2(nstrata);
      std::vector<double> cumF(nstrata);
      std::partial_sum(stratumFraction.begin(), stratumFraction.end(), cumF.begin());

      std::vector<double> n11(nstrata), n21(nstrata);
      std::vector<double> n1s(nstrata), n2s(nstrata), nss(nstrata);

      std::vector<double> analysisT1; analysisT1.reserve(K1);
      std::vector<double> analysisT2; analysisT2.reserve(K2);
      std::vector<double> totalte; totalte.reserve(N);
      std::vector<size_t> set; set.reserve(N);
      std::vector<size_t> sub; sub.reserve(N);


      for (size_t iter = begin; iter < end; ++iter) {
        // per-iteration RNG
        boost::random::mt19937_64 rng_local(seeds[iter]);
        boost::random::uniform_real_distribution<double> unif(0.0, 1.0);

        IterationResult& out = (*results)[iter];
        out.summary1Rows.clear();
        out.summary2Rows.clear();
        out.raw1Rows.clear();
        out.raw2Rows.clear();
        out.reserveForSummary1(K1);
        out.reserveForSummary2(K2);
        if (iter < maxRawIters) {
          out.reserveForRaw1(N * K1);
          out.reserveForRaw2(N * K2);
        }

        // reset block randomization
        std::fill(b1.begin(), b1.end(), allocation1);
        std::fill(b2.begin(), b2.end(), allocation2);

        double enrollt = 0.0;

        // cohort generation
        for (size_t i = 0; i < N; ++i) {
          double u = unif(rng_local);
          enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt);
          arrivalT[i] = enrollt;

          u = unif(rng_local);
          int j = findInterval1(u, cumF);
          stratum[i] = j + 1;

          // randomization
          u = unif(rng_local);
          double denom = static_cast<double>(b1[j] + b2[j]);
          double p = static_cast<double>(b1[j]) / denom;
          if (u <= p) { trtGrp[i] = 1; --b1[j]; } else { trtGrp[i] = 2; --b2[j]; }
          if (b1[j] + b2[j] == 0) { b1[j] = allocation1; b2[j] = allocation2; }

          // Plackett copula sampling for joint (binary latent and TTE u/v)
          // See equations (2.9.1), (3.3.3a) and (3.3.3b) in Nelson,
          // "An Introduction to Copulas", second edition, 2006
          // according to (3.3.3a), solving (2.9.1) for v given u and t, we have
          //   ((theta-1)-1/2*((1+(theta-1)*(u+v))^2-4*u*v*theta*(theta-1))^(-1/2)*
          //   (2*(1+(theta-1)*(u+v))*(theta-1)-4*v*theta*(theta-1))) /
          //   (2*(theta-1)) = t
          // which can be rearranged to a quadratic equation in v:
          //               a * v^2 + b * v + c = 0,
          // where
          //   a = theta + t * (1 - t) * (theta - 1)^2,
          //   b = -(theta - 2 * t * (1 - t) * (theta - 1) * (1 - (theta + 1) * u)),
          //   c = t * (1 - t) * (1 + (theta - 1) * u)^2.
          // The solution to the quadratic equation is given by
          //   v = (-b +/- sqrt(b^2 - 4*a*c)) / (2*a).
          // The feasible root is the one that satisfies equation (2.9.1).
          double u_plack = unif(rng_local);
          double t_plack = unif(rng_local);
          double v_plack = t_plack;
          if (globalOddsRatio != 1.0) {
            double s = t_plack * (1.0 - t_plack);
            double a = globalOddsRatio + s * sq(globalOddsRatio1);
            double b = -(globalOddsRatio - 2.0 * s * globalOddsRatio1 *
                         (1.0 - (globalOddsRatio + 1.0) * u_plack));
            double c = s * sq(1.0 + globalOddsRatio1 * u_plack);
            double disc = b * b - 4.0 * a * c;
            if (disc < 0.0) disc = 0.0;
            double sqrtdisc = std::sqrt(disc);
            double v1 = (-b + sqrtdisc) / (2.0 * a);
            double v2 = (-b - sqrtdisc) / (2.0 * a);
            v_plack = (std::fabs(f(v1, u_plack, t_plack)) < 1e-9) ? v1 : v2;
          }

          // latentResp and survival time: inverse CDF

          // latent response (inverse CDF sampling from logistic distribution)
          // Let X denote the latent variable following a logistic distribution with
          // location parameter loc and scale 1, then the probability of response is
          //   P(X <= 0) = exp(-loc) / (1 + exp(-loc)).
          // Given the definition of alpha1 and alpha2, we have
          //   loc = -alpha for the treatmentGroup.
          // To generate X using the inverse CDF method, let u be a uniform(0,1)
          // random variable, then we want to solve for x in
          //   u = P(X <= x) = exp((x - loc) / 1) / (1 + exp((x - loc) / 1)).
          // which leads to
          //   x = loc + log(u / (1 - u))

          // survival time (inverse CDF sampling from p.w. exponential distribution)
          if (trtGrp[i] == 1) {
            latentResp[i] = -alpha1v[j] + std::log(u_plack / (1.0 - u_plack));
            survivalT[i] = qtpwexpcpp1(v_plack, tau, lambda1x[j]);
          } else {
            latentResp[i] = -alpha2v[j] + std::log(u_plack / (1.0 - u_plack));
            survivalT[i] = qtpwexpcpp1(v_plack, tau, lambda2x[j]);
          }

          // dropout time
          u = unif(rng_local);
          if (trtGrp[i] == 1) dropoutT[i] = qtpwexpcpp1(u, tau, gamma1x[j]);
          else dropoutT[i] = qtpwexpcpp1(u, tau, gamma2x[j]);

          // treatment discontinuation and ptfu1Time for endpoint 1
          u = unif(rng_local);
          if (trtGrp[i] == 1) {
            trtDiscT[i] = qtpwexpcpp1(u, tau, delta1x[j]);
            ptfu1T[i] = std::min(trtDiscT[i], upper1);
          } else {
             trtDiscT[i] = qtpwexpcpp1(u, tau, delta2x[j]);
             ptfu1T[i] = std::min(trtDiscT[i], upper2);
          }

          // endpoint 2 observed time/event preliminary
          if (survivalT[i] <= dropoutT[i]) {
            timeObs[i] = survivalT[i]; event[i] = 1; dropEv[i] = 0;
          } else {
            timeObs[i] = dropoutT[i]; event[i] = 0; dropEv[i] = 1;
          }

          totalT[i] = arrivalT[i] + timeObs[i];
        } // cohort generation

        // Endpoint 1 calendar-time looks (plannedTime)
        analysisT1 = plannedTime;

        // compute Mantel-Haenszel for binary endpoint at each stage
        for (size_t k = 0; k < K1; ++k) {
          double time = analysisT1[k];

          // reset per-stratum counts
          std::fill(n1.begin(), n1.end(), 0);
          std::fill(n2.begin(), n2.end(), 0);

          // determine responder/source and accrual counts
          for (size_t i = 0; i < N; ++i) {
            // censoring logic at analysis time for endpoint1
            double ar = arrivalT[i];
            double sv = survivalT[i];
            double dr = dropoutT[i];
            double ptfu = ptfu1T[i];

            if (ar > time) {
              timeObs1[i] = time - ar;
              responder[i] = 255; // missing indicator
              source[i] = 0;
            } else {
              if (ar + ptfu <= time && ptfu <= sv && ptfu <= dr) {
                timeObs1[i] = ptfu;
                responder[i] = (latentResp[i] <= 0.0 ? 1 : 0);
                source[i] = 1;
              } else if (ar + sv <= time && sv <= dr && sv <= ptfu) {
                timeObs1[i] = sv;
                responder[i] = 0;
                source[i] = 2;
              } else if (ar + dr <= time && dr <= sv && dr <= ptfu) {
                timeObs1[i] = dr;
                responder[i] = 0;
                source[i] = 3;
              } else {
                timeObs1[i] = time - ar;
                responder[i] = 255;
                source[i] = 4;
              }

              size_t h = static_cast<size_t>(stratum[i] - 1);
              if (trtGrp[i] == 1) ++n1[h]; else ++n2[h];
            }
          }

          // optionally collect raw rows for endpoint1
          if (iter < maxRawIters) {
            for (size_t i = 0; i < N; ++i) {
              // skip subjects who haven't been enrolled by analysis time
              if (arrivalT[i] > time) continue;
              RawDataset1Row rr;
              rr.iterNum = static_cast<int>(iter + 1);
              rr.stageNum = static_cast<int>(k + 1);
              rr.analysisT = time;
              rr.subjectId = static_cast<int>(i + 1);
              rr.arrivalT = arrivalT[i];
              rr.stratum = stratum[i];
              rr.trtGrp = trtGrp[i];
              rr.survivalT = survivalT[i];
              rr.dropoutT = dropoutT[i];
              rr.trtDiscT = trtDiscT[i];
              rr.upper = (trtGrp[i] == 1 ? upper1 : upper2);
              rr.ptfu1T = ptfu1T[i];
              rr.timeObs = timeObs1[i];
              rr.latentResp = latentResp[i];
              rr.responder = responder[i];
              rr.source = source[i];
              out.raw1Rows.push_back(std::move(rr));
            }
          }

          int accruals1 = std::accumulate(n1.begin(), n1.end(), 0);
          int accruals2 = std::accumulate(n2.begin(), n2.end(), 0);
          int totAccruals = accruals1 + accruals2;

          // build stratified counts excluding missing (responder==255)
          std::fill(n11.begin(), n11.end(), 0.0);
          std::fill(n21.begin(), n21.end(), 0.0);
          std::fill(n1s.begin(), n1s.end(), 0.0);
          std::fill(n2s.begin(), n2s.end(), 0.0);
          std::fill(nss.begin(), nss.end(), 0.0);
          for (size_t i = 0; i < N; ++i) {
            if (responder[i] == 255) continue;
            size_t h = static_cast<size_t>(stratum[i] - 1);
            ++nss[h];
            if (trtGrp[i] == 1) {
              ++n1s[h];
              n11[h] += (responder[i] == 1 ? 1.0 : 0.0);
            }
            else {
              ++n2s[h];
              n21[h] += (responder[i] == 1 ? 1.0 : 0.0);
            }
          }

          // Mantel-Haenszel risk difference (Sato variance) computation
          double A = 0.0, B = 0.0, P = 0.0, Q = 0.0;
          for (size_t h = 0; h < nstrata; ++h) {
            if (n1s[h] <= 0.0 || n2s[h] <= 0.0 || nss[h] <= 0.0) continue;
            double dh = (n11[h] / n1s[h]) - (n21[h] / n2s[h]);
            double wh = (n1s[h] * n2s[h]) / nss[h];
            A += dh * wh;
            B += wh;
            P += (n1s[h] * n1s[h] * n21[h] - n2s[h] * n2s[h] * n11[h] +
              n1s[h] * n2s[h] * (n2s[h] - n1s[h]) * 0.5) / (nss[h] * nss[h]);
            Q += (n11[h] * (n2s[h] - n21[h]) + n21[h] * (n1s[h] - n11[h])) *
              0.5 / nss[h];
          }

          double riskDiff = A / B;
          double seRiskDiff = std::sqrt(riskDiff * P + Q) / B;
          double z = (riskDiff - riskDiffH0) / seRiskDiff;

          // collect source counts
          int source1 = 0, source2 = 0, source3 = 0;
          for (size_t i = 0; i < N; ++i) {
            if (source[i] == 1) ++source1;
            else if (source[i] == 2) ++source2;
            else if (source[i] == 3) ++source3;
          }

          // append summary row for binary endpoint
          StageSummary1Row sr;
          sr.iterNum = static_cast<int>(iter + 1);
          sr.stageNum = static_cast<int>(k + 1);
          sr.analysisT = time;
          sr.accruals1 = accruals1;
          sr.accruals2 = accruals2;
          sr.totAccruals = totAccruals;
          sr.source1 = source1;
          sr.source2 = source2;
          sr.source3 = source3;
          sr.n1 = std::accumulate(n1s.begin(), n1s.end(), 0.0);
          sr.n2 = std::accumulate(n2s.begin(), n2s.end(), 0.0);
          sr.n = std::accumulate(nss.begin(), nss.end(), 0.0);
          sr.y1 = std::accumulate(n11.begin(), n11.end(), 0.0);
          sr.y2 = std::accumulate(n21.begin(), n21.end(), 0.0);
          sr.y = sr.y1 + sr.y2;
          sr.riskDiff = riskDiff; sr.seRiskDiff = seRiskDiff; sr.z = z;
          out.summary1Rows.push_back(std::move(sr));
        } // end endpoint1 stages

        // Endpoint 2 (TTE): determine analysis times by planned and observed events
        totalte.clear();
        int nevents = 0;
        for (size_t i = 0; i < N; ++i) {
          if (event[i]) { ++nevents; totalte.push_back(totalT[i]); }
        }
        if (nevents == 0) {
          thread_utils::push_thread_warning(std::string("No events for iteration ") +
            std::to_string(iter+1) + " skipping this iteration.");
          out.summary1Rows.clear();
          out.summary2Rows.clear();
          out.raw1Rows.clear();
          out.raw2Rows.clear();
          continue;
        }
        std::sort(totalte.begin(), totalte.end());

        size_t nstages2 = K2;
        analysisT2.clear();
        size_t j = 0;
        for (j = 0; j < K2; ++j) {
          if (plannedEvents[j] >= nevents) { nstages2 = j + 1; break; }
        }
        if (j == K2) {
          for (size_t k = 0; k < nstages2; ++k) {
            analysisT2[k] = totalte[plannedEvents[k] - 1] + 1e-12;
          }
        } else {
          for (size_t k = 0; k < nstages2 - 1; ++k) {
            analysisT2[k] = totalte[plannedEvents[k] - 1] + 1e-12;
          }
          analysisT2[nstages2 - 1] = totalte.back() + 1e-12;
        }
        bool evNotAch = (nevents < plannedEvents[K2 - 1]);

        // For each TTE stage compute log-rank
        for (size_t k = 0; k < nstages2; ++k) {
          double time = analysisT2[k];
          std::fill(n1.begin(), n1.end(), 0);
          std::fill(n2.begin(), n2.end(), 0);
          int events1 = 0, events2 = 0, dropouts1 = 0, dropouts2 = 0;

          // censoring & counts
          for (size_t i = 0; i < N; ++i) {
            double ar = arrivalT[i], sv = survivalT[i], dr = dropoutT[i];
            if (ar > time) {
              timeObs[i] = time - ar; event[i] = 0; dropEv[i] = 0; continue;
            }

            if (ar + sv <= time && sv <= dr) {
              timeObs[i] = sv; event[i] = 1; dropEv[i] = 0;
            } else if (ar + dr <= time && dr <= sv) {
              timeObs[i] = dr; event[i] = 0; dropEv[i] = 1;
            } else {
              timeObs[i] = time - ar; event[i] = 0; dropEv[i] = 0;
            }

            size_t h = static_cast<size_t>(stratum[i] - 1);
            if (trtGrp[i] == 1) {
              ++n1[h];
              if (event[i]) ++events1; else if (dropEv[i]) ++dropouts1;
            } else {
              ++n2[h];
              if (event[i]) ++events2; else if (dropEv[i]) ++dropouts2;
            }
          }

          // optionally collect raw rows for endpoint2
          if (iter < maxRawIters) {
            for (size_t i = 0; i < N; ++i) {
              // skip subjects who haven't been enrolled by analysis time
              if (arrivalT[i] > time) continue;
              RawDataset2Row rr;
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
              out.raw2Rows.push_back(std::move(rr));
            }
          }

          int accruals1 = std::accumulate(n1.begin(), n1.end(), 0);
          int accruals2 = std::accumulate(n2.begin(), n2.end(), 0);
          int totAccruals = accruals1 + accruals2;
          int totEvents = events1 + events2;
          int totDropouts = dropouts1 + dropouts2;

          // sort by observed time
          sub.clear();
          for (size_t i = 0; i < N; ++i) if (timeObs[i] > 0.0) sub.push_back(i);
          std::sort(sub.begin(), sub.end(), [&](size_t a, size_t b) {
            return timeObs[a] < timeObs[b];
          });

          // compute stratified log-rank with guards
          double us = 0.0, vs = 0.0;

          for (size_t i = 0; i < sub.size(); ++i) {
            size_t idx = sub[i];
            size_t h = static_cast<size_t>(stratum[idx] - 1);
            double n1h = static_cast<double>(n1[h]);
            double n2h = static_cast<double>(n2[h]);
            double n1a = n1h * hazardRatioH0;
            double nta = n1a + n2h;
            if (event[idx]) {
              double treated = (trtGrp[idx] == 1 ? 1.0 : 0.0);
              us += (treated - n1a / nta);
              vs += (n1a * n2h) / (nta * nta);
            }
            if (trtGrp[idx] == 1) --n1[h]; else --n2[h];
          }

          double z = (vs > 0.0 ? (us / std::sqrt(vs)) : 0.0);

          // append summary row for TTE endpoint
          StageSummary2Row sr;
          sr.iterNum = static_cast<int>(iter + 1);
          sr.stageNum = static_cast<int>(k + 1);
          sr.analysisT = time;
          sr.evNotAch = evNotAch ? 1 : 0;
          sr.accruals1 = accruals1;
          sr.accruals2 = accruals2;
          sr.totAccruals = totAccruals;
          sr.events1 = events1;
          sr.events2 = events2;
          sr.totEvents = totEvents;
          sr.dropouts1 = dropouts1;
          sr.dropouts2 = dropouts2;
          sr.totDropouts = totDropouts;
          sr.uscore = us;
          sr.vscore = vs;
          sr.logRank = z;
          out.summary2Rows.push_back(std::move(sr));
        } // TTE stages
      } // iter
    } // operator()
  }; // SimWorker

  // construct and run worker
  SimWorker worker(
      K1, K2, riskDiffH0, hazardRatioH0,
      allocation1, allocation2,
      accrualTime, accrualIntensity, tau, stratumFraction,
      globalOddsRatio, pi1v, pi2v,
      lambda1x, lambda2x, gamma1x, gamma2x, delta1x, delta2x,
      upper1, upper2, N,
      plannedTime, plannedEvents,
      maxIters, maxRawIters, seeds,
      alpha1v, alpha2v, globalOddsRatio1, nstrata,
      std::function<double(const double, const double, const double)>(f),
      &results
  );

  RcppParallel::parallelFor(0, maxIters, worker);

  // Flatten results
  size_t nsr1 = 0, nsr2 = 0, nrr1 = 0, nrr2 = 0;
  for (size_t iter = 0; iter < maxIters; ++iter) {
    nsr1 += results[iter].summary1Rows.size();
    nsr2 += results[iter].summary2Rows.size();
    nrr1 += results[iter].raw1Rows.size();
    nrr2 += results[iter].raw2Rows.size();
  }
  if (nsr2 == 0) throw std::runtime_error(
    "No iterations with observed events. Unable to produce output.");

  // finalize summary and raw containers for both endpoints
  // Endpoint1 summary
  std::vector<int> sum1_iterNum; sum1_iterNum.reserve(nsr1);
  std::vector<int> sum1_stageNum; sum1_stageNum.reserve(nsr1);
  std::vector<double> sum1_analysisT; sum1_analysisT.reserve(nsr1);
  std::vector<int> sum1_accruals1; sum1_accruals1.reserve(nsr1);
  std::vector<int> sum1_accruals2; sum1_accruals2.reserve(nsr1);
  std::vector<int> sum1_totAccruals; sum1_totAccruals.reserve(nsr1);
  std::vector<int> sum1_source1; sum1_source1.reserve(nsr1);
  std::vector<int> sum1_source2; sum1_source2.reserve(nsr1);
  std::vector<int> sum1_source3; sum1_source3.reserve(nsr1);
  std::vector<double> sum1_n1; sum1_n1.reserve(nsr1);
  std::vector<double> sum1_n2; sum1_n2.reserve(nsr1);
  std::vector<double> sum1_n; sum1_n.reserve(nsr1);
  std::vector<double> sum1_y1; sum1_y1.reserve(nsr1);
  std::vector<double> sum1_y2; sum1_y2.reserve(nsr1);
  std::vector<double> sum1_y; sum1_y.reserve(nsr1);
  std::vector<double> sum1_riskDiff; sum1_riskDiff.reserve(nsr1);
  std::vector<double> sum1_seRiskDiff; sum1_seRiskDiff.reserve(nsr1);
  std::vector<double> sum1_z; sum1_z.reserve(nsr1);

  // Endpoint2 summary
  std::vector<int> sum2_iterNum; sum2_iterNum.reserve(nsr2);
  std::vector<unsigned char> sum2_evNotAch; sum2_evNotAch.reserve(nsr2);
  std::vector<int> sum2_stageNum; sum2_stageNum.reserve(nsr2);
  std::vector<double> sum2_analysisT; sum2_analysisT.reserve(nsr2);
  std::vector<int> sum2_accruals1; sum2_accruals1.reserve(nsr2);
  std::vector<int> sum2_accruals2; sum2_accruals2.reserve(nsr2);
  std::vector<int> sum2_totAccruals; sum2_totAccruals.reserve(nsr2);
  std::vector<int> sum2_events1; sum2_events1.reserve(nsr2);
  std::vector<int> sum2_events2; sum2_events2.reserve(nsr2);
  std::vector<int> sum2_totEvents; sum2_totEvents.reserve(nsr2);
  std::vector<int> sum2_dropouts1; sum2_dropouts1.reserve(nsr2);
  std::vector<int> sum2_dropouts2; sum2_dropouts2.reserve(nsr2);
  std::vector<int> sum2_totDropouts; sum2_totDropouts.reserve(nsr2);
  std::vector<double> sum2_uscore; sum2_uscore.reserve(nsr2);
  std::vector<double> sum2_vscore; sum2_vscore.reserve(nsr2);
  std::vector<double> sum2_logRank; sum2_logRank.reserve(nsr2);

  // raw containers for endpoint 1
  std::vector<int> raw1_iterNum; raw1_iterNum.reserve(nrr1);
  std::vector<int> raw1_stageNum; raw1_stageNum.reserve(nrr1);
  std::vector<double> raw1_analysisT; raw1_analysisT.reserve(nrr1);
  std::vector<int> raw1_subjectId; raw1_subjectId.reserve(nrr1);
  std::vector<double> raw1_arrivalT; raw1_arrivalT.reserve(nrr1);
  std::vector<int> raw1_stratum; raw1_stratum.reserve(nrr1);
  std::vector<int> raw1_trtGrp; raw1_trtGrp.reserve(nrr1);
  std::vector<double> raw1_survivalT; raw1_survivalT.reserve(nrr1);
  std::vector<double> raw1_dropoutT; raw1_dropoutT.reserve(nrr1);
  std::vector<double> raw1_trtDiscT; raw1_trtDiscT.reserve(nrr1);
  std::vector<double> raw1_upper; raw1_upper.reserve(nrr1);
  std::vector<double> raw1_ptfu1T; raw1_ptfu1T.reserve(nrr1);
  std::vector<double> raw1_timeObs; raw1_timeObs.reserve(nrr1);
  std::vector<double> raw1_latentResp; raw1_latentResp.reserve(nrr1);
  std::vector<unsigned char> raw1_responder; raw1_responder.reserve(nrr1);
  std::vector<int> raw1_source; raw1_source.reserve(nrr1);

  // raw containers for endpoint 2
  std::vector<int> raw2_iterNum; raw2_iterNum.reserve(nrr2);
  std::vector<int> raw2_stageNum; raw2_stageNum.reserve(nrr2);
  std::vector<double> raw2_analysisT; raw2_analysisT.reserve(nrr2);
  std::vector<int> raw2_subjectId; raw2_subjectId.reserve(nrr2);
  std::vector<double> raw2_arrivalT; raw2_arrivalT.reserve(nrr2);
  std::vector<int> raw2_stratum; raw2_stratum.reserve(nrr2);
  std::vector<int> raw2_trtGrp; raw2_trtGrp.reserve(nrr2);
  std::vector<double> raw2_survivalT; raw2_survivalT.reserve(nrr2);
  std::vector<double> raw2_dropoutT; raw2_dropoutT.reserve(nrr2);
  std::vector<double> raw2_timeObs; raw2_timeObs.reserve(nrr2);
  std::vector<unsigned char> raw2_event; raw2_event.reserve(nrr2);
  std::vector<unsigned char> raw2_dropEv; raw2_dropEv.reserve(nrr2);

  // flatten preserving iteration order
  for (size_t iter = 0; iter < maxIters; ++iter) {
    const auto& s1rows = results[iter].summary1Rows;
    for (const auto& r : s1rows) {
      sum1_iterNum.push_back(r.iterNum);
      sum1_stageNum.push_back(r.stageNum);
      sum1_analysisT.push_back(r.analysisT);
      sum1_accruals1.push_back(r.accruals1);
      sum1_accruals2.push_back(r.accruals2);
      sum1_totAccruals.push_back(r.totAccruals);
      sum1_source1.push_back(r.source1);
      sum1_source2.push_back(r.source2);
      sum1_source3.push_back(r.source3);
      sum1_n1.push_back(r.n1);
      sum1_n2.push_back(r.n2);
      sum1_n.push_back(r.n);
      sum1_y1.push_back(r.y1);
      sum1_y2.push_back(r.y2);
      sum1_y.push_back(r.y);
      sum1_riskDiff.push_back(r.riskDiff);
      sum1_seRiskDiff.push_back(r.seRiskDiff);
      sum1_z.push_back(r.z);
    }

    const auto& s2rows = results[iter].summary2Rows;
    for (const auto& r : s2rows) {
      sum2_iterNum.push_back(r.iterNum);
      sum2_evNotAch.push_back(r.evNotAch);
      sum2_stageNum.push_back(r.stageNum);
      sum2_analysisT.push_back(r.analysisT);
      sum2_accruals1.push_back(r.accruals1);
      sum2_accruals2.push_back(r.accruals2);
      sum2_totAccruals.push_back(r.totAccruals);
      sum2_events1.push_back(r.events1);
      sum2_events2.push_back(r.events2);
      sum2_totEvents.push_back(r.totEvents);
      sum2_dropouts1.push_back(r.dropouts1);
      sum2_dropouts2.push_back(r.dropouts2);
      sum2_totDropouts.push_back(r.totDropouts);
      sum2_uscore.push_back(r.uscore);
      sum2_vscore.push_back(r.vscore);
      sum2_logRank.push_back(r.logRank);
    }

    if (iter < maxRawIters) {
      const auto& r1raw = results[iter].raw1Rows;
      for (const auto& rr : r1raw) {
        raw1_iterNum.push_back(rr.iterNum);
        raw1_stageNum.push_back(rr.stageNum);
        raw1_analysisT.push_back(rr.analysisT);
        raw1_subjectId.push_back(rr.subjectId);
        raw1_arrivalT.push_back(rr.arrivalT);
        raw1_stratum.push_back(rr.stratum);
        raw1_trtGrp.push_back(rr.trtGrp);
        raw1_survivalT.push_back(rr.survivalT);
        raw1_dropoutT.push_back(rr.dropoutT);
        raw1_trtDiscT.push_back(rr.trtDiscT);
        raw1_upper.push_back(rr.upper);
        raw1_ptfu1T.push_back(rr.ptfu1T);
        raw1_timeObs.push_back(rr.timeObs);
        raw1_latentResp.push_back(rr.latentResp);
        raw1_responder.push_back(rr.responder);
        raw1_source.push_back(rr.source);
      }

      const auto& r2raw = results[iter].raw2Rows;
      for (const auto& rr : r2raw) {
        raw2_iterNum.push_back(rr.iterNum);
        raw2_stageNum.push_back(rr.stageNum);
        raw2_analysisT.push_back(rr.analysisT);
        raw2_subjectId.push_back(rr.subjectId);
        raw2_arrivalT.push_back(rr.arrivalT);
        raw2_stratum.push_back(rr.stratum);
        raw2_trtGrp.push_back(rr.trtGrp);
        raw2_survivalT.push_back(rr.survivalT);
        raw2_dropoutT.push_back(rr.dropoutT);
        raw2_timeObs.push_back(rr.timeObs);
        raw2_event.push_back(rr.event);
        raw2_dropEv.push_back(rr.dropEv);
      }
    }
  }

  // Build output DataFrames
  DataFrameCpp sumdataBIN;
  sumdataBIN.push_back(std::move(sum1_iterNum), "iterationNumber");
  sumdataBIN.push_back(std::move(sum1_stageNum), "stageNumber");
  sumdataBIN.push_back(std::move(sum1_analysisT), "analysisTime");
  sumdataBIN.push_back(std::move(sum1_accruals1), "accruals1");
  sumdataBIN.push_back(std::move(sum1_accruals2), "accruals2");
  sumdataBIN.push_back(std::move(sum1_totAccruals), "totalAccruals");
  sumdataBIN.push_back(std::move(sum1_source1), "source1");
  sumdataBIN.push_back(std::move(sum1_source2), "source2");
  sumdataBIN.push_back(std::move(sum1_source3), "source3");
  sumdataBIN.push_back(std::move(sum1_n1), "n1");
  sumdataBIN.push_back(std::move(sum1_n2), "n2");
  sumdataBIN.push_back(std::move(sum1_n), "n");
  sumdataBIN.push_back(std::move(sum1_y1), "y1");
  sumdataBIN.push_back(std::move(sum1_y2), "y2");
  sumdataBIN.push_back(std::move(sum1_y), "y");
  sumdataBIN.push_back(std::move(sum1_riskDiff), "riskDiff");
  sumdataBIN.push_back(std::move(sum1_seRiskDiff), "seRiskDiff");
  sumdataBIN.push_back(std::move(sum1_z), "mhStatistic");

  DataFrameCpp sumdataTTE;
  sumdataTTE.push_back(std::move(sum2_iterNum), "iterationNumber");
  sumdataTTE.push_back(std::move(sum2_evNotAch), "eventsNotAchieved");
  sumdataTTE.push_back(std::move(sum2_stageNum), "stageNumber");
  sumdataTTE.push_back(std::move(sum2_analysisT), "analysisTime");
  sumdataTTE.push_back(std::move(sum2_accruals1), "accruals1");
  sumdataTTE.push_back(std::move(sum2_accruals2), "accruals2");
  sumdataTTE.push_back(std::move(sum2_totAccruals), "totalAccruals");
  sumdataTTE.push_back(std::move(sum2_events1), "events1");
  sumdataTTE.push_back(std::move(sum2_events2), "events2");
  sumdataTTE.push_back(std::move(sum2_totEvents), "totalEvents");
  sumdataTTE.push_back(std::move(sum2_dropouts1), "dropouts1");
  sumdataTTE.push_back(std::move(sum2_dropouts2), "dropouts2");
  sumdataTTE.push_back(std::move(sum2_totDropouts), "totalDropouts");
  sumdataTTE.push_back(std::move(sum2_uscore), "uscore");
  sumdataTTE.push_back(std::move(sum2_vscore), "vscore");
  sumdataTTE.push_back(std::move(sum2_logRank), "logRankStatistic");

  ListCpp result;
  result.push_back(sumdataBIN, "sumdataBIN");
  result.push_back(sumdataTTE, "sumdataTTE");

  if (maxNumberOfRawDatasetsPerStage > 0) {
    DataFrameCpp rawdataBIN;
    rawdataBIN.push_back(std::move(raw1_iterNum), "iterationNumber");
    rawdataBIN.push_back(std::move(raw1_stageNum), "stageNumber");
    rawdataBIN.push_back(std::move(raw1_analysisT), "analysisTime");
    rawdataBIN.push_back(std::move(raw1_subjectId), "subjectId");
    rawdataBIN.push_back(std::move(raw1_arrivalT), "arrivalTime");
    rawdataBIN.push_back(std::move(raw1_stratum), "stratum");
    rawdataBIN.push_back(std::move(raw1_trtGrp), "treatmentGroup");
    rawdataBIN.push_back(std::move(raw1_survivalT), "survivalTime");
    rawdataBIN.push_back(std::move(raw1_dropoutT), "dropoutTime");
    rawdataBIN.push_back(std::move(raw1_trtDiscT), "trtDiscTime");
    rawdataBIN.push_back(std::move(raw1_upper), "trtDurUpperLimit");
    rawdataBIN.push_back(std::move(raw1_ptfu1T), "ptfu1Time");
    rawdataBIN.push_back(std::move(raw1_timeObs), "timeUnderObservation");
    rawdataBIN.push_back(std::move(raw1_latentResp), "latentResponse");
    rawdataBIN.push_back(std::move(raw1_responder), "responder");
    rawdataBIN.push_back(std::move(raw1_source), "source");
    result.push_back(rawdataBIN, "rawdataBIN");

    DataFrameCpp rawdataTTE;
    rawdataTTE.push_back(std::move(raw2_iterNum), "iterationNumber");
    rawdataTTE.push_back(std::move(raw2_stageNum), "stageNumber");
    rawdataTTE.push_back(std::move(raw2_analysisT), "analysisTime");
    rawdataTTE.push_back(std::move(raw2_subjectId), "subjectId");
    rawdataTTE.push_back(std::move(raw2_arrivalT), "arrivalTime");
    rawdataTTE.push_back(std::move(raw2_stratum), "stratum");
    rawdataTTE.push_back(std::move(raw2_trtGrp), "treatmentGroup");
    rawdataTTE.push_back(std::move(raw2_survivalT), "survivalTime");
    rawdataTTE.push_back(std::move(raw2_dropoutT), "dropoutTime");
    rawdataTTE.push_back(std::move(raw2_timeObs), "timeUnderObservation");
    rawdataTTE.push_back(std::move(raw2_event), "event");
    rawdataTTE.push_back(std::move(raw2_dropEv), "dropoutEvent");
    result.push_back(rawdataTTE, "rawdataTTE");
  }

  return result;
}


//' @title Simulation for a Binary and a Time-to-Event Endpoint in
//' Group Sequential Trials
//' @description Performs simulation for two-endpoint two-arm group
//' sequential trials.
//' \itemize{
//'   \item Endpoint 1: Binary endpoint, analyzed using the
//'         Mantel-Haenszel test for risk difference.
//'   \item Endpoint 2: Time-to-event endpoint, analyzed using
//'         the log-rank test for treatment effect.
//' }
//' The analysis times for the binary endpoint are based on calendar times,
//' while the time-to-event analyses are triggered by reaching the
//' pre-specified number of events. The binary endpoint is
//' assessed at the first post-treatment follow-up visit (PTFU1).
//'
//' @param kMax1 Number of stages for the binary endpoint.
//' @param kMax2 Number of stages for the time-to-event endpoint.
//' @param riskDiffH0 Risk difference under the null hypothesis for the
//'   binary endpoint.
//' @param hazardRatioH0 Hazard ratio under the null hypothesis for the
//'   time-to-event endpoint.
//' @param allocation1 Number of subjects in the treatment group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in the control group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param globalOddsRatio Global odds ratio of the Plackett copula
//'   linking the two endpoints.
//' @param pi1 Response probabilities by stratum for the treatment group
//'   for the binary endpoint.
//' @param pi2 Response probabilities by stratum for the control group
//'   for the binary endpoint.
//' @param lambda1 A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the treatment group for the time-to-event
//'   endpoint.
//' @param lambda2 A vector of hazard rates for the event in each analysis
//'   time interval by stratum for the control group for the time-to-event
//'   endpoint.
//' @param gamma1 The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the treatment group.
//' @param gamma2 The hazard rate for exponential dropout, a vector of
//'   hazard rates for piecewise exponential dropout applicable for all
//'   strata, or a vector of hazard rates for dropout in each analysis time
//'   interval by stratum for the control group.
//' @param delta1 The hazard rate for exponential treatment discontinuation,
//'   a vector of hazard rates for piecewise exponential treatment
//'   discontinuation applicable for all strata, or a vector of hazard rates
//'   for treatment discontinuation in each analysis time interval by
//'   stratum for the treatment group for the binary endpoint.
//' @param delta2 The hazard rate for exponential treatment discontinuation,
//'   a vector of hazard rates for piecewise exponential treatment
//'   discontinuation applicable for all strata, or a vector of hazard rates
//'   for treatment discontinuation in each analysis time interval by
//'   stratum for the control group for the binary endpoint.
//' @param upper1 Maximim protocol-specified treatment duration for
//'   the treatment group.
//' @param upper2 Maximum protocol-specified treatment duration for
//'   the control group.
//' @param n Sample size.
//' @param plannedTime Calendar times for the analyses of the binary
//'   endpoint.
//' @param plannedEvents Target cumulative number of events for
//'   the time-to-event analyses.
//' @param maxNumberOfIterations Number of simulation iterations to perform.
//' @param maxNumberOfRawDatasetsPerStage Number of subject-level datasets
//'   to retain per stage. Set to 0 to skip raw data saving.
//' @param seed The seed to reproduce the simulation results.
//'
//' @details We consider dual primary endpoints with endpoint 1 being a
//'   binary endpoint and endpoint 2 being a time-to-event endpoint.
//'   The analyses of endpoint 1 will be based on calendar times, while
//'   the analyses of endpoint 2 will be based on the number of events.
//'   Therefore, the analyses of the two endpoints are not at the same
//'   time points. The correlation between the two endpoints is
//'   characterized by the global odds ratio of the Plackett copula.
//'   In addition, the time-to-event endpoint will render the binary
//'   endpoint as a non-responder, and so does the dropout. In addition,
//'   the treatment discontinuation will impact the number of available
//'   subjects for analysis. The administrative censoring will exclude
//'   subjects from the analysis of the binary endpoint.
//'
//' @return A list with 4 components:
//'
//' * \code{sumdataBIN}: A data frame of summary data by iteration and stage
//'   for the binary endpoint:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage number, covering all stages even if
//'       the trial stops at an interim look.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{accruals1}: The number of subjects enrolled at the stage for
//'       the treatment group.
//'
//'     - \code{accruals2}: The number of subjects enrolled at the stage for
//'       the control group.
//'
//'     - \code{totalAccruals}: The total number of subjects enrolled at
//'       the stage.
//'
//'     - \code{source1}: The total number of subjects with response status
//'       determined by the underlying latent response variable.
//'
//'     - \code{source2}: The total number of subjects with response status
//'       (non-responder) determined by experiencing the event for the
//'       time-to-event endpoint.
//'
//'     - \code{source3}: The total number of subjects with response status
//'       (non-responder) determined by dropping out prior to the PTFU1
//'       visit.
//'
//'     - \code{n1}: The number of subjects included in the analysis of
//'       the binary endpoint for the treatment group.
//'
//'     - \code{n2}: The number of subjects included in the analysis of
//'       the binary endpoint for the control group.
//'
//'     - \code{n}: The total number of subjects included in the analysis of
//'       the binary endpoint at the stage.
//'
//'     - \code{y1}: The number of responders for the binary endpoint in
//'       the treatment group.
//'
//'     - \code{y2}: The number of responders for the binary endpoint in
//'       the control group.
//'
//'     - \code{y}: The total number of responders for the binary endpoint
//'       at the stage.
//'
//'     - \code{riskDiff}: The estimated risk difference for the binary
//'       endpoint.
//'
//'     - \code{seRiskDiff}: The standard error for risk difference based on
//'       the Sato approximation.
//'
//'     - \code{mhStatistic}: The Mantel-Haenszel test Z-statistic for
//'       the binary endpoint.
//'
//' * \code{sumdataTTE}: A data frame of summary data by iteration and stage
//'   for the time-to-event endpoint:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{eventsNotAchieved}: Whether the target number of events
//'       is not achieved for the iteration.
//'
//'     - \code{stageNumber}: The stage number, covering all stages even if
//'       the trial stops at an interim look.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{accruals1}: The number of subjects enrolled at the stage for
//'       the treatment group.
//'
//'     - \code{accruals2}: The number of subjects enrolled at the stage for
//'       the control group.
//'
//'     - \code{totalAccruals}: The total number of subjects enrolled at
//'       the stage.
//'
//'     - \code{events1}: The number of events at the stage for
//'       the treatment group.
//'
//'     - \code{events2}: The number of events at the stage for
//'       the control group.
//'
//'     - \code{totalEvents}: The total number of events at the stage.
//'
//'     - \code{dropouts1}: The number of dropouts at the stage for
//'       the treatment group.
//'
//'     - \code{dropouts2}: The number of dropouts at the stage for
//'       the control group.
//'
//'     - \code{totalDropouts}: The total number of dropouts at the stage.
//'
//'     - \code{uscore}: The numerator of the log-rank test statistic for the
//'       time-to-event endpoint.
//'
//'     - \code{vscore}: The variance of the log-rank test statistic for the
//'       time-to-event endpoint.
//'
//'     - \code{logRankStatistic}: The log-rank test Z-statistic for
//'       the time-to-event endpoint.
//'
//' * \code{rawdataBIN} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//'   positive integer): A data frame for subject-level data for the binary
//'   endpoint for selected replications, containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage under consideration.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{subjectId}: The subject ID.
//'
//'     - \code{arrivalTime}: The enrollment time for the subject.
//'
//'     - \code{stratum}: The stratum for the subject.
//'
//'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
//'       subject.
//'
//'     - \code{survivalTime}: The underlying survival time for the
//'       time-to-event endpoint for the subject.
//'
//'     - \code{dropoutTime}: The underlying dropout time for the
//'       time-to-event endpoint for the subject.
//'
//'     - \code{trtDiscTime}: The underlying treatment discontinuation time
//'       for the binary endpoint for the subject.
//'
//'     - \code{trtDurUpperLimit}: The maximum protocol-specified treatment
//'       duration for the subject based on the treatment group assignment.
//'
//'     - \code{ptfu1Time}:The underlying assessment time for the
//'       binary endpoint for the subject.
//'
//'     - \code{timeUnderObservation}: The time under observation
//'       since randomization for the binary endpoint for the subject.
//'
//'     - \code{latentResponse}: The underlying latent response variable for
//'       the binary endpoint for the subject, which determines the response
//'       status for the binary endpoint at PTFU1 visit.
//'
//'     - \code{responder}: Whether the subject is a responder for the
//'       binary endpoint.
//'
//'     - \code{source}: The source of the determination of responder
//'       status for the binary endpoint: = 1 based on the underlying
//'       latent response variable, = 2 based on the occurrence of
//'       the time-to-event endpoint before the assessment time of the
//'       binary endpoint (imputed as a non-responder), = 3 based on
//'       the dropout before the assessment time of the binary endpoint
//'       (imputed as a non-responder), = 4 excluded from analysis
//'       due to administrative censoring.
//'
//' * \code{rawdataTTE} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//'   positive integer): A data frame for subject-level data for the
//'   time-to-event endpoint for selected replications, containing the
//'   following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage under consideration.
//'
//'     - \code{analysisTime}: The time for the stage since trial start.
//'
//'     - \code{subjectId}: The subject ID.
//'
//'     - \code{arrivalTime}: The enrollment time for the subject.
//'
//'     - \code{stratum}: The stratum for the subject.
//'
//'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
//'       subject.
//'
//'     - \code{survivalTime}: The underlying survival time for the
//'       time-to-event endpoint for the subject.
//'
//'     - \code{dropoutTime}: The underlying dropout time for the
//'       time-to-event endpoint for the subject.
//'
//'     - \code{timeUnderObservation}: The time under observation
//'       since randomization for the time-to-event endpoint for the subject.
//'
//'     - \code{event}: Whether the subject experienced the event for the
//'       time-to-event endpoint.
//'
//'     - \code{dropoutEvent}: Whether the subject dropped out for the
//'       time-to-event endpoint.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' tcut = c(0, 12, 36, 48)
//' surv = c(1, 0.95, 0.82, 0.74)
//' lambda2 = (log(surv[1:3]) - log(surv[2:4]))/(tcut[2:4] - tcut[1:3])
//'
//' sim1 = binary_tte_sim(
//'   kMax1 = 1,
//'   kMax2 = 2,
//'   accrualTime = seq(0, 8),
//'   accrualIntensity = 40/9 * seq(1, 9),
//'   piecewiseSurvivalTime = c(0, 12, 36),
//'   globalOddsRatio = 1,
//'   pi1 = 0.80,
//'   pi2 = 0.65,
//'   lambda1 = 0.65*lambda2,
//'   lambda2 = lambda2,
//'   gamma1 = -log(1-0.04)/12,
//'   gamma2 = -log(1-0.04)/12,
//'   delta1 = -log(1-0.02)/12,
//'   delta2 = -log(1-0.02)/12,
//'   upper1 = 15*28/30.4,
//'   upper2 = 12*28/30.4,
//'   n = 640,
//'   plannedTime = 20 + 15*28/30.4,
//'   plannedEvents = c(130, 173),
//'   maxNumberOfIterations = 1000,
//'   maxNumberOfRawDatasetsPerStage = 1,
//'   seed = 314159)
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::List binary_tte_sim(
    const int kMax1 = 1,
    const int kMax2 = 1,
    const double riskDiffH0 = 0,
    const double hazardRatioH0 = 1,
    const int allocation1 = 1,
    const int allocation2 = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const double globalOddsRatio = 1,
    const Rcpp::NumericVector& pi1 = NA_REAL,
    const Rcpp::NumericVector& pi2 = NA_REAL,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const Rcpp::NumericVector& delta1 = 0,
    const Rcpp::NumericVector& delta2 = 0,
    const double upper1 = NA_REAL,
    const double upper2 = NA_REAL,
    const int n = NA_INTEGER,
    const Rcpp::NumericVector& plannedTime = NA_REAL,
    const Rcpp::IntegerVector& plannedEvents = NA_INTEGER,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = 0) {

  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto pi1v = Rcpp::as<std::vector<double>>(pi1);
  auto pi2v = Rcpp::as<std::vector<double>>(pi2);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);
  auto del1 = Rcpp::as<std::vector<double>>(delta1);
  auto del2 = Rcpp::as<std::vector<double>>(delta2);
  auto plannedE = Rcpp::as<std::vector<int>>(plannedEvents);
  auto plannedT = Rcpp::as<std::vector<double>>(plannedTime);

  auto out = binary_tte_sim_cpp(
    kMax1, kMax2, riskDiffH0, hazardRatioH0, allocation1, allocation2,
    accrualT, accrualInt, pwSurvT, stratumFrac, globalOddsRatio,
    pi1v, pi2v, lam1, lam2, gam1, gam2, del1, del2,
    upper1, upper2, n, plannedT, plannedE,
    maxNumberOfIterations, maxNumberOfRawDatasetsPerStage, seed);

  thread_utils::drain_thread_warnings_to_R();

  return Rcpp::wrap(out);
}
