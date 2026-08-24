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
#include <vector>

#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/random.hpp>

using std::size_t;

ListCpp rdsim_multiarm_cpp(const size_t M, const size_t kMax,
                           const FlatMatrix &criticalValues,
                           const std::vector<double> &futilityBounds,
                           const std::vector<double> &riskDiffH0s,
                           const std::vector<double> &allocations,
                           const std::vector<double> &pis,
                           const bool nullVariance, const size_t n,
                           const std::vector<int> &plannedSubjects,
                           const int maxNumberOfIterations, const int seed) {

  if (M < 1)
    throw std::invalid_argument("M must be at least 1");
  if (kMax < 1)
    throw std::invalid_argument("kMax must be at least 1");
  if (criticalValues.nrow != kMax || criticalValues.ncol != M) {
    throw std::invalid_argument("criticalValues must have dimension kMax x M");
  }
  if (kMax > 1 && futilityBounds.size() < kMax - 1) {
    throw std::invalid_argument("futilityBounds must have length >= kMax - 1");
  }
  for (size_t k = 0; k + 1 < kMax; ++k) {
    if (futilityBounds[k] > criticalValues(k, 0)) {
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
      throw std::invalid_argument(nm + " must be provided");
    if (pis[m] <= 0.0 || pis[m] >= 1.0) {
      throw std::invalid_argument(nm + " must lies in (0,1)");
    }
  }

  if (static_cast<int>(n) == INT_MIN)
    throw std::invalid_argument("n must be provided");
  if (n <= 0)
    throw std::invalid_argument("n must be a positive integer");

  if (plannedSubjects.size() != kMax)
    throw std::invalid_argument("Invalid length for plannedSubjects");
  if (plannedSubjects[0] <= 0)
    throw std::invalid_argument("plannedSubjects must be positive");
  if (any_nonincreasing(plannedSubjects))
    throw std::invalid_argument("plannedSubjects must be increasing");
  int maxN =
      static_cast<int>(std::round((allocs[0] + allocs[M]) / sumAlloc * n));
  if (plannedSubjects.back() != maxN) {
    throw std::invalid_argument("plannedSubjects must end with sample size for "
                                "arm 1 & control combined");
  }

  if (maxNumberOfIterations < 1)
    throw std::invalid_argument(
        "maxNumberOfIterations must be a positive integer");

  const size_t maxIters = static_cast<size_t>(maxNumberOfIterations);
  std::vector<int> stagewiseSubjects(kMax);
  stagewiseSubjects[0] = plannedSubjects[0];
  for (size_t k = 1; k < kMax; ++k) {
    stagewiseSubjects[k] = plannedSubjects[k] - plannedSubjects[k - 1];
  }

  std::vector<uint64_t> seeds(maxIters);
  boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
  for (size_t iter = 0; iter < maxIters; ++iter)
    seeds[iter] = master_rng();

  struct StageSummary1Row {
    int iterNum = 0;
    int stageNum = 0;
    int trtGrp = 0;
    int accruals = 0;
    int events = 0;
    double phat = 0.0;
  };

  struct StageSummary2Row {
    int iterNum = 0;
    int stageNum = 0;
    int actArm = 1;
    int totAccruals = 0;
    int totEvents = 0;
    double riskDiff = 0.0;
    double vriskDiff = 0.0;
    double z = 0.0;
    unsigned char reject = 0;
    unsigned char futility = 0;
  };

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

  std::vector<IterationResult> results;
  results.resize(maxIters);

  struct SimWorker : public RcppParallel::Worker {
    const size_t M;
    const size_t kMax;
    const std::vector<double> &rdH0s;
    const std::vector<double> &allocs;
    const std::vector<double> &pis;
    const bool nullVariance;
    const std::vector<int> &stagewiseSubjects;
    const std::vector<uint64_t> &seeds;
    std::vector<IterationResult> *results;

    SimWorker(size_t M_, size_t kMax_, const std::vector<double> &rdH0s_,
              const std::vector<double> &allocs_,
              const std::vector<double> &pis_, bool nullVariance_,
              const std::vector<int> &stagewiseSubjects_,
              const std::vector<uint64_t> &seeds_,
              std::vector<IterationResult> *results_)
        : M(M_), kMax(kMax_), rdH0s(rdH0s_), allocs(allocs_), pis(pis_),
          nullVariance(nullVariance_), stagewiseSubjects(stagewiseSubjects_),
          seeds(seeds_), results(results_) {}

    void operator()(std::size_t begin, std::size_t end) {
      const size_t M1 = M + 1;
      const size_t M2 = M + 2;
      std::vector<int> accruals(M2), events(M2);
      std::vector<double> phat(M2);
      std::vector<double> mr(2);

      for (size_t iter = begin; iter < end; ++iter) {
        boost::random::mt19937_64 rng_local(seeds[iter]);
        IterationResult &out = (*results)[iter];
        out.summary1Rows.clear();
        out.summary2Rows.clear();
        out.reserveForSummary1(kMax * M2);
        out.reserveForSummary2(kMax * M);

        for (size_t k = 0; k < kMax; ++k) {
          const double stageSize = static_cast<double>(stagewiseSubjects[k]) /
                                   (allocs[0] + allocs[M]);

          accruals[M1] = 0;
          events[M1] = 0;

          for (size_t m = 0; m < M1; ++m) {
            const int N = static_cast<int>(std::round(allocs[m] * stageSize));
            boost::random::binomial_distribution<> dist(N, pis[m]);
            const int x = dist(rng_local);

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

  SimWorker worker(M, kMax, rdH0s, allocs, pis, nullVariance, stagewiseSubjects,
                   seeds, &results);

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

  std::vector<double> haveStage(kMax);
  FlatMatrix rejectByArm(kMax, M + 1);
  FlatMatrix futilityByArm(kMax, M + 1);
  FlatMatrix eventsByArm(kMax, M + 2);
  FlatMatrix subjectsByArm(kMax, M + 2);
  std::vector<double> expEventsByArm(M + 2);
  std::vector<double> expSubjectsByArm(M + 2);

  size_t ntests = 1u << M; // number of possible rejection patterns among M arms
  FlatMatrix rejectBySet(kMax, ntests);
  FlatMatrix rejectByNum(kMax + 1, M + 1);

  int *stop1 = sum1_stopStage.data();
  const int *sum1_E = sum1_events.data();
  const int *sum1_A = sum1_accruals.data();
  int *stop2 = sum2_stopStage.data();
  const double *zstat = sum2_z.data();
  unsigned char *reject = sum2_reject.data();
  unsigned char *futility = sum2_futility.data();

  for (size_t iter = 0; iter < niters; ++iter) {
    const size_t i1 = iter * rowsPerIter1;
    const size_t i2 = iter * rowsPerIter2;

    // find stopping stage for this iteration
    size_t stop_k = kMax - 1;
    bool stoppedForFutility = false;

    for (size_t k = 0; k < kMax; ++k) {
      const size_t offset = i2 + k * M;

      double cut = criticalValues(k, 0); // cutoff for level-M test

      std::unordered_set<size_t> I; // set of unrejected hypotheses
      for (size_t m = 0; m < M; ++m)
        I.insert(m);

      // check whether there is any arm that crosses the efficacy boundary
      bool anyreject = false;
      for (auto it = I.begin(); it != I.end();) {
        size_t m = *it;
        if (zstat[offset + m] > cut) {
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
            double cut2 = criticalValues(k, M - I.size());
            for (auto it = I.begin(); it != I.end();) {
              size_t m = *it;
              if (zstat[offset + m] > cut2) {
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

      bool allFutile = false;
      if (k + 1 < kMax) {
        allFutile = true;
        for (size_t m = 0; m < M; ++m) {
          if (zstat[offset + m] > futilityBounds[k]) {
            allFutile = false;
            break;
          }
        }
      } else {
        allFutile = true; // final-stage futility if no active arm is rejected
      }

      if (allFutile) {
        stop_k = k;
        stoppedForFutility = true;
        futilityByArm(k, M) += 1;
        for (size_t m = 0; m < M; ++m) {
          futilityByArm(k, m) += 1;
        }
        break;
      }
    }

    // record which hypotheses are rejected at the stopping stage for this
    // iteration
    const size_t offset = i2 + stop_k * M;
    std::vector<unsigned char> cc(M);
    for (size_t m = 0; m < M; ++m) {
      cc[m] = reject[offset + m];
    }
    size_t num_reject = std::count(cc.begin(), cc.end(), 1);
    rejectByNum(stop_k, num_reject) += 1;
    for (size_t k = 0; k < kMax; ++k) {
      if (k != stop_k) {
        rejectByNum(k, 0) += 1; // no rejections at other stages
        rejectBySet(k, 0) += 1; // empty set of rejections at other stages
      }
    }

    size_t value = 0;
    for (size_t m = 0; m < M; ++m) {
      if (cc[m])
        value |= (1u << m);
    }
    rejectBySet(stop_k, value) += 1;

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

    // tally time and number of events/dropouts/subjects
    for (size_t k = 0; k <= stop_k; ++k) {
      haveStage[k] += 1;
      const size_t offset1 = i1 + k * (M + 2);

      for (size_t m = 0; m < M + 2; ++m) {
        size_t idx = offset1 + m;
        eventsByArm(k, m) += sum1_E[idx];
        subjectsByArm(k, m) += sum1_A[idx];

        if (k == stop_k) {
          expEventsByArm[m] += sum1_E[idx];
          expSubjectsByArm[m] += sum1_A[idx];
        }
      }
    }
  }

  // empirical cumulative rejection rates by stage and treatment
  for (size_t m = 0; m < M + 1; ++m) {
    for (size_t k = 0; k < kMax; ++k) {
      rejectByArm(k, m) /= niters;
      futilityByArm(k, m) /= niters;
      rejectByNum(k, m) /= niters;
    }
  }

  for (size_t s = 0; s < ntests; ++s) {
    for (size_t k = 0; k < kMax; ++k) {
      rejectBySet(k, s) /= niters;
    }
  }

  // fill the last row of rejectByNum with column sums
  for (size_t m = 1; m < M + 1; ++m) {
    double sum = 0.0;
    for (size_t k = 0; k < kMax; ++k) {
      sum += rejectByNum(k, m);
    }
    rejectByNum(kMax, m) = sum;
  }
  double sum = 0.0;
  for (size_t m = 1; m < M + 1; ++m) {
    sum += rejectByNum(kMax, m);
  }
  rejectByNum(kMax, 0) = 1.0 - sum;

  // build column sums of rejectBySet
  std::vector<double> setSums(ntests);
  for (size_t s = 1; s < ntests; ++s) {
    double sum = 0.0;
    for (size_t k = 0; k < kMax; ++k) {
      sum += rejectBySet(k, s);
    }
    setSums[s] = sum;
  }
  setSums[0] = 1.0 - std::accumulate(setSums.begin() + 1, setSums.end(), 0.0);

  std::vector<std::string> intersectHyp(ntests);
  for (size_t s = 0; s < ntests; ++s) {
    std::string str;
    for (size_t m = 0; m < M; ++m) {
      if (s & (1u << m)) {
        if (!str.empty())
          str += ",";
        str += std::to_string(m + 1);
      }
    }
    if (str.empty())
      str = "none";
    intersectHyp[s] = str;
  }

  DataFrameCpp rejectSetData;
  rejectSetData.push_back(intersectHyp, "intersectionHypotheses");
  rejectSetData.push_back(setSums, "rejectionProbability");

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

  // convert counts to proportions / averages
  for (size_t m = 0; m < M + 2; ++m) {
    expEventsByArm[m] /= niters;
    expSubjectsByArm[m] /= niters;
  }

  for (size_t m = 0; m < M + 2; ++m) {
    for (size_t k = 0; k < kMax; ++k) {
      const double denom = haveStage[k];
      eventsByArm(k, m) /= denom;
      subjectsByArm(k, m) /= denom;
    }
  }

  ListCpp overview;
  overview.push_back(overallReject, "overallReject");
  overview.push_back(overallFutility, "overallFutility");
  overview.push_back(std::move(rejectByArm), "rejectPerStage");
  overview.push_back(std::move(futilityByArm), "futilityPerStage");
  overview.push_back(std::move(cumRejectByArm), "cumulativeRejection");
  overview.push_back(std::move(cumFutilityByArm), "cumulativeFutility");
  overview.push_back(std::move(rejectByNum), "rejectByNumber");
  overview.push_back(std::move(rejectSetData), "rejectBySet");
  overview.push_back(std::move(eventsByArm), "numberOfEvents");
  overview.push_back(std::move(subjectsByArm), "numberOfSubjects");
  overview.push_back(std::move(expEventsByArm), "expectedNumberOfEvents");
  overview.push_back(std::move(expSubjectsByArm), "expectedNumberOfSubjects");
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
  overview.push_back(kMax, "kMax");

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
Rcpp::List rdsim_multiarm_Rcpp(
    const int M = 2, const int kMax = 1,
    const Rcpp::Nullable<Rcpp::NumericMatrix> criticalValues = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> futilityBounds = R_NilValue,
    const Rcpp::NumericVector &riskDiffH0s = 0,
    const Rcpp::NumericVector &allocations = 1,
    const Rcpp::NumericVector &pis = NA_REAL, const bool nullVariance = true,
    const int n = NA_INTEGER,
    const Rcpp::IntegerVector &plannedSubjects = NA_INTEGER,
    const int maxNumberOfIterations = 1000, const int seed = 0) {

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

  std::vector<double> futBounds;
  if (futilityBounds.isNotNull()) {
    futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
    if (kMax > 1 && static_cast<int>(futBounds.size()) < kMax - 1) {
      throw std::invalid_argument(
          "futilityBounds must have length >= kMax - 1");
    }
  } else {
    futBounds = std::vector<double>(std::max(0, kMax - 1), -8.0);
  }

  std::vector<double> rdH0s(riskDiffH0s.begin(), riskDiffH0s.end());
  std::vector<double> allocs(allocations.begin(), allocations.end());
  std::vector<double> ps(pis.begin(), pis.end());
  std::vector<int> plannedSubjects_vec(plannedSubjects.begin(),
                                       plannedSubjects.end());

  auto out = rdsim_multiarm_cpp(M, kMax, critValues, futBounds, rdH0s, allocs,
                                ps, nullVariance, n, plannedSubjects_vec,
                                maxNumberOfIterations, seed);

  thread_utils::drain_thread_warnings_to_R();

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "rdsim_multiarm";
  return result;
}
