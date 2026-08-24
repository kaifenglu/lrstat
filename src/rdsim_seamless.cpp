#include "dataframe_list.h"
#include "miettinen_nurminen.h"
#include "thread_utils.h"
#include "utilities.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/random.hpp>

using std::size_t;

// Parallel entry function
ListCpp rdsim_seamless_cpp(const size_t M, const size_t K, const size_t rankp0,
                           const std::vector<double> &criticalValues,
                           const std::vector<double> &futilityBounds,
                           const std::vector<double> &riskDiffH0s,
                           const std::vector<double> &allocations,
                           const std::vector<double> &pis,
                           const bool nullVariance, const size_t n,
                           const std::vector<int> &plannedSubjects,
                           const int maxNumberOfIterations, const int seed) {
  if (M < 1)
    throw std::invalid_argument("M must be at least 1");
  if (M < 1)
    throw std::invalid_argument("M must be at least 1");
  if (rankp0 < 1 || rankp0 > M) {
    throw std::invalid_argument("rankp0 must be an integer between 1 and M");
  }
  if (K < 1)
    throw std::invalid_argument("K must be at least 1");
  size_t kMax = K + 1;
  if (K > 0 && futilityBounds.size() < K) {
    throw std::invalid_argument("futilityBounds must have length >= K");
  }
  for (size_t k = 0; k < K; ++k) {
    if (futilityBounds[k] > criticalValues[k]) {
      throw std::invalid_argument(
          "futilityBounds must lie below criticalValues");
    }
  }

  std::vector<double> rdH0s = expand1(riskDiffH0s, M, "riskDiffH0s");
  std::vector<double> allocs = expand1(allocations, M + 1, "allocations");
  if (std::any_of(allocs.begin(), allocs.end(), [](double v) { return v < 1; }))
    throw std::invalid_argument("All allocation parameters must be positive");
  const double sumAlloc = std::accumulate(allocs.begin(), allocs.end(), 0.0);

  if (pis.size() != M + 1) {
    throw std::invalid_argument("pis must have length M + 1");
  }

  for (size_t m = 0; m < M + 1; ++m) {
    std::string nm = std::string("pis[") + std::to_string(m) + "]";
    if (std::isnan(pis[m]))
      throw std::invalid_argument(nm + "must be provided");
    if (pis[m] <= 0.0 || pis[m] >= 1.0) {
      throw std::invalid_argument(nm + "must lies in (0,1)");
    }
  }

  if (static_cast<int>(n) == INT_MIN)
    throw std::invalid_argument("n must be provided");
  if (n <= 0)
    throw std::invalid_argument("n must be a positive integer");

  // cum number of patients enrolled at each stage for arm 1 & control combined
  int maxN =
      static_cast<int>(std::round((allocs[0] + allocs[M]) / sumAlloc * n));
  if (plannedSubjects[0] <= 0)
    throw std::invalid_argument("plannedSubjects must be positive");
  if (plannedSubjects.size() != kMax)
    throw std::invalid_argument("Invalid length for plannedSubjects");
  if (any_nonincreasing(plannedSubjects))
    throw std::invalid_argument("plannedSubjects must be increasing");
  if (plannedSubjects.back() != maxN)
    throw std::invalid_argument("plannedSubjects must end with sample size for "
                                "arm 1 & control combined");

  // number of new patients enrolled in each stage for arm 1 & control combined
  std::vector<int> stagewiseSubjects(kMax);
  stagewiseSubjects[0] = plannedSubjects[0];
  for (size_t k = 1; k < kMax; ++k) {
    stagewiseSubjects[k] = plannedSubjects[k] - plannedSubjects[k - 1];
  }

  if (maxNumberOfIterations < 1)
    throw std::invalid_argument(
        "maxNumberOfIterations must be a positive integer");

  size_t maxIters = static_cast<size_t>(maxNumberOfIterations);

  // generate seeds for each iteration to ensure reproducibility
  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter)
    seeds[iter] = master_rng();

  // One summary (stage-level) row produced by an iteration
  struct StageSummary1Row {
    int iterNum = 0;
    int stageNum = 0;
    int trtGrp = 0;
    int accruals = 0, events = 0;
    double phat = 0.0;
  };

  struct StageSummary2Row {
    int iterNum = 0;
    int stageNum = 0;
    int actArm = 1;
    int totAccruals = 0, totEvents = 0;
    double riskDiff = 0.0, vriskDiff = 0.0, z = 0.0;
  };

  // Per-iteration container written exclusively by the worker thread
  struct IterationResult {
    std::vector<StageSummary1Row> summary1Rows;
    std::vector<StageSummary2Row> summary2Rows;
    void reserveForSummary1(size_t approxRows) {
      summary1Rows.reserve(approxRows);
    }
    void reserveForSummary2(size_t approxRows) {
      summary2Rows.reserve(approxRows);
    }
  };

  // pre-size per-iteration results
  std::vector<IterationResult> results;
  results.resize(maxIters);

  // Worker that runs simulation iterations [begin, end)
  struct SimWorker : public RcppParallel::Worker {
    // inputs (const refs)
    const size_t M;
    const size_t kMax;
    const std::vector<double> &rdH0s;
    const std::vector<double> &allocs;
    const std::vector<double> &pis;
    const bool nullVariance;
    const size_t n;
    const std::vector<int> &stagewiseSubjects;
    const std::vector<uint64_t> &seeds;

    // output pointer (pre-sized vector of IterationResult)
    std::vector<IterationResult> *results;

    SimWorker(size_t M_, size_t kMax_, const std::vector<double> &rdH0s_,
              const std::vector<double> &allocs_,
              const std::vector<double> &pis_, bool nullVariance_, size_t n_,
              const std::vector<int> &stagewiseSubjects_,
              const std::vector<uint64_t> &seeds_,
              std::vector<IterationResult> *results_)
        : M(M_), kMax(kMax_), rdH0s(rdH0s_), allocs(allocs_), pis(pis_),
          nullVariance(nullVariance_), n(n_),
          stagewiseSubjects(stagewiseSubjects_), seeds(seeds_),
          results(results_) {}

    void operator()(std::size_t begin, std::size_t end) {
      // local buffers reused by this worker
      size_t M1 = M + 1;
      size_t M2 = M + 2;
      std::vector<int> accruals(M2), events(M2);
      std::vector<double> phat(M2);
      std::vector<double> mr(2);

      for (size_t iter = begin; iter < end; ++iter) {
        // deterministic per-iteration RNG
        boost::random::mt19937_64 rng_local(seeds[iter]);

        // per-iteration output container
        IterationResult &out = (*results)[iter];
        out.summary1Rows.clear();
        out.summary2Rows.clear();
        out.reserveForSummary1(kMax * M2); // all arms & overall
        out.reserveForSummary2(kMax * M);  // pairwise comparisons with control

        for (size_t k = 0; k < kMax; ++k) {
          // number of new patients enrolled in stage k for arm 1 & control
          // combined per unit allocation, i.e. per "allocation share" in stage
          // k
          const double stageSize = static_cast<double>(stagewiseSubjects[k]) /
                                   (allocs[0] + allocs[M]);

          accruals[M1] = 0; // overall accruals/events for all arms combined
          events[M1] = 0;

          // simulate stage k data for all arms;
          // response counts for arm m ~ Binomial(N, pis[m])
          for (size_t m = 0; m < M1; ++m) {
            // incremental number of patients enrolled in arm m in stage k
            int N = static_cast<int>(std::round(allocs[m] * stageSize));
            boost::random::binomial_distribution<> dist(N, pis[m]);
            int x = dist(rng_local);

            // update cumulative counts for arm m; calculate phat for arm m
            if (k == 0) {
              accruals[m] = N;
              events[m] = x;
            } else {
              accruals[m] = N + accruals[m];
              events[m] = x + events[m];
            }
            phat[m] = static_cast<double>(events[m]) / accruals[m];

            // update overall counts and phat
            accruals[M1] += accruals[m];
            events[M1] += events[m];
          }
          phat[M1] = static_cast<double>(events[M1]) / accruals[M1];

          for (size_t m = 0; m < M2; ++m) {
            StageSummary1Row sr1;
            sr1.iterNum = static_cast<int>(iter + 1);
            sr1.stageNum = static_cast<int>(k + 1);
            sr1.trtGrp = static_cast<int>(m + 1);
            sr1.accruals = accruals[m];
            sr1.events = events[m];
            sr1.phat = phat[m];
            out.summary1Rows.push_back(sr1);
          }

          for (size_t m = 0; m < M; ++m) {
            StageSummary2Row sr2;
            sr2.iterNum = static_cast<int>(iter + 1);
            sr2.stageNum = static_cast<int>(k + 1);
            sr2.actArm = static_cast<int>(m + 1);
            sr2.totAccruals = accruals[m] + accruals[M];
            sr2.totEvents = events[m] + events[M];

            const double rd = phat[m] - phat[M];
            double vrd = 0.0;
            if (nullVariance) {
              const double r = allocs[m] / (allocs[m] + allocs[M]);
              // restricted maximum likelihood estimates
              mr = remlRiskDiff(r, r * pis[m], 1.0 - r, (1.0 - r) * pis[M],
                                rdH0s[m]);
              const double pT = mr[0];
              const double pC = mr[1];
              vrd =
                  pT * (1.0 - pT) / accruals[m] + pC * (1.0 - pC) / accruals[M];
            } else {
              vrd = phat[m] * (1.0 - phat[m]) / accruals[m] +
                    phat[M] * (1.0 - phat[M]) / accruals[M];
            }
            const double z = (rd - rdH0s[m]) / std::sqrt(vrd);
            sr2.riskDiff = rd;
            sr2.vriskDiff = vrd;
            sr2.z = z;
            out.summary2Rows.push_back(sr2);
          }
        } // per-stage
      } // per-iteration
    } // operator()
  }; // SimWorker

  // run worker in parallel
  SimWorker worker(M, kMax, rdH0s, allocs, pis, nullVariance, n,
                   stagewiseSubjects, seeds, &results);

  RcppParallel::parallelFor(0, maxIters, worker);

  // Flatten results
  size_t ns1r = 0, ns2r = 0;
  for (size_t iter = 0; iter < maxIters; ++iter) {
    ns1r += results[iter].summary1Rows.size();
    ns2r += results[iter].summary2Rows.size();
  }

  // prepare final containers (reserve capacities)
  std::vector<int> sum1_iterNum;
  sum1_iterNum.reserve(ns1r);
  std::vector<int> sum1_stopStage;
  sum1_stopStage.reserve(ns1r);
  std::vector<int> sum1_stageNum;
  sum1_stageNum.reserve(ns1r);
  std::vector<int> sum1_trtGrp;
  sum1_trtGrp.reserve(ns1r);
  std::vector<int> sum1_accruals;
  sum1_accruals.reserve(ns1r);
  std::vector<int> sum1_events;
  sum1_events.reserve(ns1r);
  std::vector<double> sum1_phat;
  sum1_phat.reserve(ns1r);

  std::vector<int> sum2_iterNum;
  sum2_iterNum.reserve(ns2r);
  std::vector<int> sum2_selectedArm;
  sum2_selectedArm.reserve(ns2r);
  std::vector<int> sum2_stopStage;
  sum2_stopStage.reserve(ns2r);
  std::vector<int> sum2_stageNum;
  sum2_stageNum.reserve(ns2r);
  std::vector<int> sum2_actArm;
  sum2_actArm.reserve(ns2r);
  std::vector<int> sum2_totAccruals;
  sum2_totAccruals.reserve(ns2r);
  std::vector<int> sum2_totEvents;
  sum2_totEvents.reserve(ns2r);
  std::vector<double> sum2_riskDiff;
  sum2_riskDiff.reserve(ns2r);
  std::vector<double> sum2_vriskDiff;
  sum2_vriskDiff.reserve(ns2r);
  std::vector<double> sum2_z;
  sum2_z.reserve(ns2r);
  std::vector<unsigned char> sum2_reject;
  sum2_reject.reserve(ns2r);
  std::vector<unsigned char> sum2_futility;
  sum2_futility.reserve(ns2r);

  // flatten by iteration in order (preserves iteration order)
  for (size_t iter = 0; iter < maxIters; ++iter) {
    const auto &s1rows = results[iter].summary1Rows;
    for (const auto &r : s1rows) {
      sum1_iterNum.push_back(r.iterNum);
      sum1_stopStage.push_back(0);
      sum1_stageNum.push_back(r.stageNum);
      sum1_trtGrp.push_back(r.trtGrp);
      sum1_accruals.push_back(r.accruals);
      sum1_events.push_back(r.events);
      sum1_phat.push_back(r.phat);
    }

    const auto &s2rows = results[iter].summary2Rows;
    for (const auto &r : s2rows) {
      sum2_iterNum.push_back(r.iterNum);
      sum2_selectedArm.push_back(0);
      sum2_stopStage.push_back(0);
      sum2_stageNum.push_back(r.stageNum);
      sum2_actArm.push_back(r.actArm);
      sum2_totAccruals.push_back(r.totAccruals);
      sum2_totEvents.push_back(r.totEvents);
      sum2_riskDiff.push_back(r.riskDiff);
      sum2_vriskDiff.push_back(r.vriskDiff);
      sum2_z.push_back(r.z);
      sum2_reject.push_back(0);
      sum2_futility.push_back(0);
    }
  }

  const size_t rowsPerIter1 = kMax * (M + 2);
  const size_t rowsPerIter2 = kMax * M;
  const size_t niters = sum2_iterNum.size() / rowsPerIter2;
  std::vector<double> selectionProb(M, 0.0);
  std::vector<double> selectToStage2(M, 0.0);
  double selectAnyToStage2 = 0.0;
  double expNumEvents = 0.0;
  double expNumSubjects = 0.0;

  FlatMatrix rejectByArm(kMax, M + 1);
  FlatMatrix futilityByArm(kMax, M + 1);
  FlatMatrix haveStage(kMax, M + 1);
  FlatMatrix eventsByArm(kMax, M + 1);
  FlatMatrix subjectsByArm(kMax, M + 1);

  int *stop1 = sum1_stopStage.data();
  const int *sum1_E = sum1_events.data();
  const int *sum1_A = sum1_accruals.data();
  int *stop2 = sum2_stopStage.data();
  int *selectedArmVec = sum2_selectedArm.data();
  const int *sum2_totE = sum2_totEvents.data();
  const int *sum2_totA = sum2_totAccruals.data();
  const double *zstat = sum2_z.data();
  unsigned char *reject = sum2_reject.data();
  unsigned char *futility = sum2_futility.data();

  for (size_t iter = 0; iter < niters; ++iter) {
    const size_t i1 = iter * rowsPerIter1;
    const size_t i2 = iter * rowsPerIter2;

    // choose rankp0-th arm by phase-2 z statistic (larger is better)
    std::vector<std::pair<double, size_t>> ranked(M);
    for (size_t m = 0; m < M; ++m) {
      ranked[m] = std::make_pair(zstat[i2 + m], m);
    }
    std::stable_sort(
        ranked.begin(), ranked.end(),
        [](const std::pair<double, size_t> &a,
           const std::pair<double, size_t> &b) { return a.first > b.first; });
    size_t selected_arm = ranked[rankp0 - 1].second;
    selectionProb[selected_arm] += 1.0;

    for (size_t k = 0; k < kMax; ++k) {
      for (size_t m = 0; m < M; ++m) {
        size_t idx = i2 + k * M + m;
        selectedArmVec[idx] = static_cast<int>(selected_arm + 1);
      }
    }

    // stage-1 totals (all active arms + control)
    const double events1 = sum1_E[i1 + M + 1];
    const double subjects1 = sum1_A[i1 + M + 1];

    // stage-1 totals for selected active arm + control
    const double events2 = sum2_totE[i2 + selected_arm];
    const double subjects2 = sum2_totA[i2 + selected_arm];

    // stage-1 totals for non-selected active arms
    const double devents1 = events1 - events2;
    const double dsubjects1 = subjects1 - subjects2;

    // find stopping stage for this iter (default last stage)
    size_t stop_k = kMax - 1;
    bool stoppedForFutility = false;
    bool stoppedForEfficacy = false;

    // stopping in phase 2
    const size_t phase2_offset = i2;
    const double selectedZP2 = zstat[phase2_offset + selected_arm];
    bool anyRejectP2 = (selectedZP2 > criticalValues[0]);
    bool allFutileP2 = (selectedZP2 < futilityBounds[0]);

    if (anyRejectP2) {
      stop_k = 0;
      stoppedForEfficacy = true;
      rejectByArm(0, M) += 1;
      for (size_t r = rankp0 - 1; r < M; ++r) {
        const size_t m = ranked[r].second;
        size_t idx2 = i2 + m;
        if (zstat[idx2] > criticalValues[0]) {
          reject[idx2] = 1;
          rejectByArm(0, m) += 1;
        }
      }
    } else if (allFutileP2) {
      stop_k = 0;
      stoppedForFutility = true;
      futilityByArm(0, M) += 1;
      for (size_t m = 0; m < M; ++m) {
        futilityByArm(0, m) += 1;
      }
    } else {
      // The selected arm continues to stage 2 when not stopping in phase 2.
      selectToStage2[selected_arm] += 1.0;
      selectAnyToStage2 += 1.0;

      for (size_t k = 1; k < kMax; ++k) {
        const size_t idx = i2 + k * M + selected_arm;

        if (zstat[idx] > criticalValues[k]) {
          stop_k = k;
          stoppedForEfficacy = true;
          reject[idx] = 1;
          rejectByArm(k, selected_arm) += 1;
          rejectByArm(k, M) += 1;
          break;
        }

        bool futileNow = false;
        if (k < kMax - 1) {
          futileNow = (zstat[idx] < futilityBounds[k]);
        } else { // final stage futile if selected arm cannot be rejected
          futileNow = !stoppedForEfficacy;
        }

        if (futileNow) {
          stop_k = k;
          stoppedForFutility = true;
          futilityByArm(k, selected_arm) += 1;
          futilityByArm(k, M) += 1;
          break;
        }
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

    if (stoppedForFutility) {
      const size_t offset2 = i2 + stop_k * M;
      for (size_t m = 0; m < M; ++m) {
        futility[offset2 + m] = 1;
      }
    }

    // expected values depending on stop_k
    if (stop_k == 0) {
      expNumEvents += events1;
      expNumSubjects += subjects1;
    } else {
      const size_t idx = i2 + stop_k * M + selected_arm;
      expNumEvents += devents1 + sum2_totE[idx];
      expNumSubjects += dsubjects1 + sum2_totA[idx];
    }

    // everyone reaches stage 1 by definition
    haveStage(0, selected_arm) += 1.0;
    eventsByArm(0, selected_arm) += events1;
    subjectsByArm(0, selected_arm) += subjects1;

    // subsequent stages up to stop_k
    for (size_t j = 1; j <= stop_k; ++j) {
      const size_t idx = i2 + j * M + selected_arm;
      haveStage(j, selected_arm) += 1.0;
      eventsByArm(j, selected_arm) += devents1 + sum2_totE[idx];
      subjectsByArm(j, selected_arm) += dsubjects1 + sum2_totA[idx];
    }
  }

  // compute totals for "all active arms" column
  for (size_t m = 0; m < M; ++m) {
    for (size_t k = 0; k < kMax; ++k) {
      haveStage(k, M) += haveStage(k, m);
      eventsByArm(k, M) += eventsByArm(k, m);
      subjectsByArm(k, M) += subjectsByArm(k, m);
    }
  }

  // convert counts to proportions / averages
  expNumEvents /= niters;
  expNumSubjects /= niters;
  for (size_t m = 0; m < M; ++m)
    selectionProb[m] /= niters;
  for (size_t m = 0; m < M; ++m)
    selectToStage2[m] /= niters;
  selectAnyToStage2 /= niters;

  for (size_t m = 0; m < M + 1; ++m) {
    for (size_t k = 0; k < kMax; ++k) {
      const double denom = haveStage(k, m) > 0.0 ? haveStage(k, m) : 1.0;
      eventsByArm(k, m) /= denom;
      subjectsByArm(k, m) /= denom;
      rejectByArm(k, m) /= niters;
      futilityByArm(k, m) /= niters;
    }
  }

  // cumulative rejection by stage
  FlatMatrix cumRejectByArm(kMax, M + 1);
  FlatMatrix cumFutilityByArm(kMax, M + 1);
  for (size_t m = 0; m < M + 1; ++m) {
    cumRejectByArm(0, m) = rejectByArm(0, m);
    cumFutilityByArm(0, m) = futilityByArm(0, m);
    for (size_t k = 1; k < kMax; ++k) {
      cumRejectByArm(k, m) = cumRejectByArm(k - 1, m) + rejectByArm(k, m);
      cumFutilityByArm(k, m) = cumFutilityByArm(k - 1, m) + futilityByArm(k, m);
    }
  }

  double overallReject = cumRejectByArm(kMax - 1, M);
  double overallFutility = cumFutilityByArm(kMax - 1, M);

  ListCpp overview;
  overview.push_back(std::move(selectionProb), "selectionProb");
  overview.push_back(std::move(selectToStage2), "selectToStage2");
  overview.push_back(selectAnyToStage2, "selectAnyToStage2");
  overview.push_back(std::move(rejectByArm), "rejectPerStage");
  overview.push_back(std::move(futilityByArm), "futilityPerStage");
  overview.push_back(std::move(cumRejectByArm), "cumulativeRejection");
  overview.push_back(std::move(cumFutilityByArm), "cumulativeFutility");
  overview.push_back(std::move(eventsByArm), "numberOfEvents");
  overview.push_back(std::move(subjectsByArm), "numberOfSubjects");
  overview.push_back(overallReject, "overallReject");
  overview.push_back(overallFutility, "overallFutility");
  overview.push_back(expNumEvents, "expectedNumberOfEvents");
  overview.push_back(expNumSubjects, "expectedNumberOfSubjects");
  overview.push_back(criticalValues, "criticalValues");
  overview.push_back(futilityBounds, "futilityBounds");
  overview.push_back(std::move(rdH0s), "riskDiffH0s");
  overview.push_back(nullVariance, "nullVariance");
  overview.push_back(niters, "numberOfIterations");
  overview.push_back(n, "n");
  overview.push_back(std::move(allocs), "allocations");
  overview.push_back(std::move(pis), "responseRates");
  overview.push_back(std::move(plannedSubjects), "plannedSubjects");
  overview.push_back(M, "M");
  overview.push_back(K, "K");
  overview.push_back(rankp0, "rankp0");

  DataFrameCpp sumdata1;
  sumdata1.push_back(std::move(sum1_iterNum), "iterationNumber");
  sumdata1.push_back(std::move(sum1_stopStage), "stopStage");
  sumdata1.push_back(std::move(sum1_stageNum), "stageNumber");
  sumdata1.push_back(std::move(sum1_trtGrp), "treatmentGroup");
  sumdata1.push_back(std::move(sum1_accruals), "accruals");
  sumdata1.push_back(std::move(sum1_events), "events");
  sumdata1.push_back(std::move(sum1_phat), "phat");

  DataFrameCpp sumdata2;
  sumdata2.push_back(std::move(sum2_iterNum), "iterationNumber");
  sumdata2.push_back(std::move(sum2_selectedArm), "selectedArm");
  sumdata2.push_back(std::move(sum2_stopStage), "stopStage");
  sumdata2.push_back(std::move(sum2_stageNum), "stageNumber");
  sumdata2.push_back(std::move(sum2_actArm), "activeArm");
  sumdata2.push_back(std::move(sum2_totAccruals), "totalAccruals");
  sumdata2.push_back(std::move(sum2_totEvents), "totalEvents");
  sumdata2.push_back(std::move(sum2_riskDiff), "riskDiff");
  sumdata2.push_back(std::move(sum2_vriskDiff), "vriskDiff");
  sumdata2.push_back(std::move(sum2_z), "riskDiffZ");
  sumdata2.push_back(std::move(sum2_reject), "reject");
  sumdata2.push_back(std::move(sum2_futility), "futility");

  ListCpp result;
  result.push_back(std::move(overview), "overview");
  result.push_back(std::move(sumdata1), "sumdata1");
  result.push_back(std::move(sumdata2), "sumdata2");

  return result;
}

// [[Rcpp::export]]
Rcpp::List rdsim_seamless_Rcpp(
    const int M = 2, const int K = 1, const int rankp0 = 1,
    const Rcpp::NumericVector &criticalValues = NA_REAL,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBounds = R_NilValue,
    const Rcpp::NumericVector &riskDiffH0s = 1,
    const Rcpp::NumericVector &allocations = 1,
    const Rcpp::NumericVector &pis = NA_REAL, const bool nullVariance = true,
    const int n = NA_INTEGER,
    const Rcpp::IntegerVector &plannedSubjects = NA_INTEGER,
    const int maxNumberOfIterations = 1000, const int seed = 0) {

  std::vector<double> critValues(criticalValues.begin(), criticalValues.end());

  std::vector<double> futBounds;
  if (futilityBounds.isNotNull()) {
    futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
    if (K > 0 && static_cast<int>(futBounds.size()) < K) {
      throw std::invalid_argument("futilityBounds must have length >= K");
    }
  } else {
    futBounds = std::vector<double>(std::max(0, K), -8.0);
  }

  std::vector<double> rdH0s(riskDiffH0s.begin(), riskDiffH0s.end());
  std::vector<double> allocs(allocations.begin(), allocations.end());
  std::vector<double> ps(pis.begin(), pis.end());
  std::vector<int> plannedSubjects_vec(plannedSubjects.begin(),
                                       plannedSubjects.end());

  auto out = rdsim_seamless_cpp(
      M, K, static_cast<size_t>(rankp0), critValues, futBounds, rdH0s, allocs,
      ps, nullVariance, n, plannedSubjects_vec, maxNumberOfIterations, seed);

  thread_utils::drain_thread_warnings_to_R();

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "rdsim_seamless";

  return result;
}
