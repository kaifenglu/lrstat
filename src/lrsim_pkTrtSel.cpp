#include "dataframe_list.h"
#include "thread_utils.h"
#include "utilities.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using std::size_t;

namespace {

constexpr double ALPHA_ONE_SIDED = 0.025;

enum Method : size_t {
  FOLLOWUP_WISE = 0,
  PATIENT_WISE,
  DUNNETT,
  STAGE2_ONLY,
  THREE_ARM,
  NMETHOD
};

const char *const METHOD_NAME[NMETHOD] = {"followupwise", "patientwise",
                                          "dunnett", "stage2only", "threearm"};

struct Subject {
  int arm = 0;
  int stage = 1;
  double arrival = 0.0;
  double survival = 0.0;
  double exposure = NaN;
};

struct LogrankStat {
  double score = 0.0;
  double variance = 0.0;
  int events = 0;

  double z() const {
    return variance > 0.0 ? score / std::sqrt(variance) : 0.0;
  }
};

struct TrialResult {
  bool completed = false;
  int selected = 0;
  int correct = 0;
  double meanExposureHigh = NaN;
  double meanExposureLow = NaN;
  double finalZ = NaN;
  double estimatedHazardRatio = NaN;
  double pFollowupElementary1 = NaN;
  double pFollowupIntersection1 = NaN;
  double pFollowup2 = NaN;
  double pPatientElementary1 = NaN;
  double pPatientIntersection1 = NaN;
  double pPatient2 = NaN;
  double pDunnett = NaN;
  double pStage2Only = NaN;
  double pThreeArm = NaN;
  int followupEvents1 = 0;
  int followupEvents2 = 0;
  int patientEvents1 = 0;
  int patientEvents2 = 0;
};

double piecewise_accrual_time(double cumulative_intensity) {
  constexpr double first_end = 7.0 * 4.0;
  constexpr double second_end = first_end + 15.0 * 4.0;
  if (cumulative_intensity <= first_end)
    return cumulative_intensity / 7.0;
  if (cumulative_intensity <= second_end)
    return 4.0 + (cumulative_intensity - first_end) / 15.0;
  return 8.0 + (cumulative_intensity - second_end) / 22.0;
}

double hochberg_intersection(double p1, double p2) {
  return std::min(2.0 * std::min(p1, p2), std::max(p1, p2));
}

double combination_p(double p1, double p2, double w1, double w2) {
  const double z = w1 * boost_qnorm(1.0 - p1) + w2 * boost_qnorm(1.0 - p2);
  return 1.0 - boost_pnorm(z);
}

LogrankStat logrank_stat(const std::vector<Subject> &subjects, int treatment,
                         double cutoff, int stage) {
  struct Observation {
    double time;
    int event;
    int treatment;
  };

  std::vector<Observation> observations;
  observations.reserve(subjects.size());
  for (const Subject &subject : subjects) {
    if (subject.arm != 0 && subject.arm != treatment)
      continue;
    if (stage != 0 && subject.stage != stage)
      continue;
    if (subject.arrival >= cutoff)
      continue;
    const double followup = cutoff - subject.arrival;
    observations.push_back({std::min(subject.survival, followup),
                            subject.survival <= followup ? 1 : 0,
                            subject.arm == treatment ? 1 : 0});
  }

  std::sort(observations.begin(), observations.end(),
            [](const Observation &left, const Observation &right) {
              return left.time < right.time;
            });

  int riskTreatment = 0;
  int riskControl = 0;
  for (const Observation &observation : observations) {
    if (observation.treatment)
      ++riskTreatment;
    else
      ++riskControl;
  }

  LogrankStat result;
  size_t index = 0;
  while (index < observations.size()) {
    const double time = observations[index].time;
    size_t end = index;
    int eventsTreatment = 0;
    int eventsControl = 0;
    int removeTreatment = 0;
    int removeControl = 0;
    while (end < observations.size() && observations[end].time == time) {
      if (observations[end].treatment) {
        ++removeTreatment;
        eventsTreatment += observations[end].event;
      } else {
        ++removeControl;
        eventsControl += observations[end].event;
      }
      ++end;
    }

    const int risk = riskTreatment + riskControl;
    const int events = eventsTreatment + eventsControl;
    if (events > 0 && risk > 1) {
      result.score +=
          eventsTreatment - static_cast<double>(riskTreatment) * events / risk;
      result.variance += static_cast<double>(riskTreatment) * riskControl *
                         events * (risk - events) /
                         (static_cast<double>(risk) * risk * (risk - 1.0));
      result.events += events;
    }
    riskTreatment -= removeTreatment;
    riskControl -= removeControl;
    index = end;
  }
  return result;
}

int count_events(const std::vector<Subject> &subjects, double cutoff,
                 int stage) {
  return static_cast<int>(std::count_if(
      subjects.begin(), subjects.end(), [&](const Subject &subject) {
        return (stage == 0 || subject.stage == stage) &&
               subject.arrival + subject.survival <= cutoff;
      }));
}

double dunnett_p(double selected_z) {
  const std::vector<double> lower(2, -POS_INF);
  const std::vector<double> upper(2, -selected_z);
  return 1.0 - pbvnormcpp(lower, upper, 0.5);
}

struct SimWorker : public RcppParallel::Worker {
  const std::vector<double> &medianSurvival;
  const std::vector<double> &meanExposure;
  const double exposureSD;
  const double rho;
  const int interimPatients;
  const double interimFollowup;
  const int totalPatients;
  const double finalAnalysisTime;
  const std::string &selectionRule;
  const std::vector<uint64_t> &seeds;
  std::vector<TrialResult> *results;

  SimWorker(const std::vector<double> &medianSurvival_,
            const std::vector<double> &meanExposure_, double exposureSD_,
            double rho_, int interimPatients_, double interimFollowup_,
            int totalPatients_, double finalAnalysisTime_,
            const std::string &selectionRule_,
            const std::vector<uint64_t> &seeds_,
            std::vector<TrialResult> *results_)
      : medianSurvival(medianSurvival_), meanExposure(meanExposure_),
        exposureSD(exposureSD_), rho(rho_), interimPatients(interimPatients_),
        interimFollowup(interimFollowup_), totalPatients(totalPatients_),
        finalAnalysisTime(finalAnalysisTime_), selectionRule(selectionRule_),
        seeds(seeds_), results(results_) {}

  void operator()(size_t begin, size_t end) {
    const double sqrtOneMinusRho2 = std::sqrt(1.0 - rho * rho);
    std::vector<double> logExposureMean(2), logExposureSD(2);
    for (size_t arm = 0; arm < 2; ++arm) {
      const double mean = meanExposure[arm];
      const double variance = exposureSD * exposureSD;
      logExposureMean[arm] =
          2.0 * std::log(mean) - 0.5 * std::log(mean * mean + variance);
      logExposureSD[arm] =
          std::sqrt(std::log(mean * mean + variance) - 2.0 * std::log(mean));
    }

    for (size_t iteration = begin; iteration < end; ++iteration) {
      try {
        boost::random::mt19937_64 rng(seeds[iteration]);
        boost::random::uniform_real_distribution<double> uniform(0.0, 1.0);
        boost::random::normal_distribution<double> normal(0.0, 1.0);
        TrialResult &result = (*results)[iteration];

        std::vector<double> arrivals(totalPatients);
        double cumulativeIntensity = 0.0;
        for (int index = 0; index < totalPatients; ++index) {
          cumulativeIntensity += -std::log(uniform(rng));
          arrivals[index] = piecewise_accrual_time(cumulativeIntensity);
        }
        const double interimTime =
            arrivals[interimPatients - 1] + interimFollowup;
        const int stage1Patients = interimPatients;
        const int overrunEnd = static_cast<int>(
            std::upper_bound(arrivals.begin(), arrivals.end(), interimTime) -
            arrivals.begin());

        std::vector<double> survivalZ(totalPatients);
        for (double &value : survivalZ)
          value = normal(rng);

        std::vector<int> randomization(overrunEnd);
        for (int index = 0; index < overrunEnd; ++index)
          randomization[index] = index % 3;
        std::shuffle(randomization.begin(), randomization.end(), rng);

        std::vector<Subject> subjects;
        subjects.reserve(totalPatients);
        std::vector<double> exposureSum(2, 0.0);
        std::vector<int> exposureCount(2, 0);
        for (int index = 0; index < overrunEnd; ++index) {
          const int arm = randomization[index];
          const double zSurvival = survivalZ[index];
          const double survivalUniform =
              std::min(std::max(boost_pnorm(zSurvival), 1e-12), 1.0 - 1e-12);
          const double rate =
              std::log(2.0) / medianSurvival[arm == 0 ? 2 : arm - 1];
          Subject subject;
          subject.arm = arm;
          subject.stage = index < stage1Patients ? 1 : 2;
          subject.arrival = arrivals[index];
          subject.survival = -std::log(1.0 - survivalUniform) / rate;
          if (arm > 0) {
            const double zExposure =
                rho * zSurvival + sqrtOneMinusRho2 * normal(rng);
            subject.exposure = std::exp(logExposureMean[arm - 1] +
                                        logExposureSD[arm - 1] * zExposure);
            exposureSum[arm - 1] += subject.exposure;
            ++exposureCount[arm - 1];
          }
          subjects.push_back(subject);
        }

        result.meanExposureHigh = exposureSum[0] / exposureCount[0];
        result.meanExposureLow = exposureSum[1] / exposureCount[1];
        const LogrankStat interimHigh =
            logrank_stat(subjects, 1, interimTime, 1);
        const LogrankStat interimLow =
            logrank_stat(subjects, 2, interimTime, 1);

        boost::random::mt19937_64 comparatorRng = rng;
        std::vector<Subject> threeArmSubjects = subjects;
        std::vector<int> comparatorRandomization(totalPatients - overrunEnd);
        for (size_t index = 0; index < comparatorRandomization.size(); ++index)
          comparatorRandomization[index] = static_cast<int>(index % 3);
        std::shuffle(comparatorRandomization.begin(),
                     comparatorRandomization.end(), comparatorRng);
        for (int index = overrunEnd; index < totalPatients; ++index) {
          const int arm = comparatorRandomization[index - overrunEnd];
          const double survivalUniform = std::min(
              std::max(boost_pnorm(normal(comparatorRng)), 1e-12), 1.0 - 1e-12);
          const double rate =
              std::log(2.0) / medianSurvival[arm == 0 ? 2 : arm - 1];
          threeArmSubjects.push_back({arm, 2, arrivals[index],
                                      -std::log(1.0 - survivalUniform) / rate,
                                      NaN});
        }
        const LogrankStat threeArmHigh =
            logrank_stat(threeArmSubjects, 1, finalAnalysisTime, 0);
        const LogrankStat threeArmLow =
            logrank_stat(threeArmSubjects, 2, finalAnalysisTime, 0);
        result.pThreeArm =
            dunnett_p(std::min(threeArmHigh.z(), threeArmLow.z()));

        if (selectionRule == "exposure") {
          result.selected =
              result.meanExposureHigh >= 1.5 * result.meanExposureLow ? 1 : 2;
        } else if (selectionRule == "os") {
          result.selected = interimHigh.z() < interimLow.z() ? 1 : 2;
        }

        std::vector<int> stage2Randomization(totalPatients - overrunEnd);
        for (size_t index = 0; index < stage2Randomization.size(); ++index)
          stage2Randomization[index] = static_cast<int>(index % 2);
        std::shuffle(stage2Randomization.begin(), stage2Randomization.end(),
                     rng);
        auto add_stage2 = [&](std::vector<Subject> &trialSubjects,
                              int treatment) {
          for (int index = overrunEnd; index < totalPatients; ++index) {
            const int arm =
                stage2Randomization[index - overrunEnd] == 0 ? 0 : treatment;
            const double survivalUniform = std::min(
                std::max(boost_pnorm(survivalZ[index]), 1e-12), 1.0 - 1e-12);
            const double rate =
                std::log(2.0) / medianSurvival[arm == 0 ? 2 : arm - 1];
            trialSubjects.push_back({arm, 2, arrivals[index],
                                     -std::log(1.0 - survivalUniform) / rate,
                                     NaN});
          }
        };

        if (selectionRule == "perfect") {
          std::vector<Subject> highSubjects = subjects;
          std::vector<Subject> lowSubjects = subjects;
          add_stage2(highSubjects, 1);
          add_stage2(lowSubjects, 2);
          const LogrankStat finalHigh =
              logrank_stat(highSubjects, 1, finalAnalysisTime, 0);
          const LogrankStat finalLow =
              logrank_stat(lowSubjects, 2, finalAnalysisTime, 0);
          const double logHazardRatioHigh =
              finalHigh.score / finalHigh.variance;
          const double logHazardRatioLow = finalLow.score / finalLow.variance;
          result.selected = logHazardRatioHigh < logHazardRatioLow ? 1 : 2;
          subjects = result.selected == 1 ? std::move(highSubjects)
                                          : std::move(lowSubjects);
        } else {
          add_stage2(subjects, result.selected);
        }

        const LogrankStat interimSelected =
            result.selected == 1 ? interimHigh : interimLow;
        const LogrankStat finalSelected =
            logrank_stat(subjects, result.selected, finalAnalysisTime, 0);
        const LogrankStat finalStage2Selected =
            logrank_stat(subjects, result.selected, finalAnalysisTime, 2);
        const LogrankStat finalStage1High =
            logrank_stat(subjects, 1, finalAnalysisTime, 1);
        const LogrankStat finalStage1Low =
            logrank_stat(subjects, 2, finalAnalysisTime, 1);

        const double interimPHigh = boost_pnorm(interimHigh.z());
        const double interimPLow = boost_pnorm(interimLow.z());
        result.pFollowupElementary1 = boost_pnorm(interimSelected.z());
        result.pFollowupIntersection1 =
            hochberg_intersection(interimPHigh, interimPLow);
        const double incrementVariance =
            finalSelected.variance - interimSelected.variance;
        const double incrementZ =
            incrementVariance > 0.0
                ? (finalSelected.score - interimSelected.score) /
                      std::sqrt(incrementVariance)
                : 0.0;
        result.pFollowup2 = boost_pnorm(incrementZ);

        const double patientPHigh = boost_pnorm(finalStage1High.z());
        const double patientPLow = boost_pnorm(finalStage1Low.z());
        result.pPatientElementary1 =
            result.selected == 1 ? patientPHigh : patientPLow;
        result.pPatientIntersection1 =
            hochberg_intersection(patientPHigh, patientPLow);
        result.pPatient2 = boost_pnorm(finalStage2Selected.z());
        result.pStage2Only = result.pPatient2;
        result.finalZ = finalSelected.z();
        result.pDunnett = dunnett_p(result.finalZ);
        result.estimatedHazardRatio =
            finalSelected.variance > 0.0
                ? std::exp(finalSelected.score / finalSelected.variance)
                : NaN;

        result.followupEvents1 = count_events(subjects, interimTime, 1);
        result.followupEvents2 =
            std::max(0, finalSelected.events - interimSelected.events);
        result.patientEvents1 = count_events(subjects, finalAnalysisTime, 1);
        result.patientEvents2 = finalStage2Selected.events;
        const int correctArm = medianSurvival[0] > medianSurvival[1] ? 1 : 2;
        result.correct = result.selected == correctArm;
        result.completed = true;
      } catch (const std::exception &error) {
        thread_utils::push_thread_warning(
            "iteration " + std::to_string(iteration + 1) + ": " + error.what());
      }
    }
  }
};

} // namespace

ListCpp lrsim_pkTrtSel_cpp(
    const std::vector<double> &medianSurvival,
    const std::vector<double> &meanExposure, const double exposureSD,
    const double rho, const int interimPatients, const double interimFollowup,
    const int totalPatients, const double finalAnalysisTime,
    const std::string &selectionRule, const size_t ntrial, const int seed) {
  if (medianSurvival.size() != 3)
    throw std::invalid_argument(
        "medianSurvival must contain high-dose, low-dose, and control values");
  if (std::any_of(medianSurvival.begin(), medianSurvival.end(),
                  [](double value) { return value <= 0.0; }))
    throw std::invalid_argument("medianSurvival values must be positive");
  if (meanExposure.size() != 2)
    throw std::invalid_argument(
        "meanExposure must contain high-dose and low-dose values");
  if (std::any_of(meanExposure.begin(), meanExposure.end(),
                  [](double value) { return value <= 0.0; }))
    throw std::invalid_argument("meanExposure values must be positive");
  if (exposureSD <= 0.0)
    throw std::invalid_argument("exposureSD must be positive");
  if (rho < -1.0 || rho > 1.0)
    throw std::invalid_argument("corrSurvivalExposure must lie in [-1, 1]");
  if (interimPatients < 3 || interimPatients >= totalPatients)
    throw std::invalid_argument(
        "interimPatients must be at least 3 and less than totalPatients");
  if (interimFollowup < 0.0)
    throw std::invalid_argument("interimFollowup must be nonnegative");
  if (totalPatients <= 0 || finalAnalysisTime <= 0.0 || ntrial == 0)
    throw std::invalid_argument("totalPatients, finalAnalysisTime, and "
                                "maxNumberOfIterations must be positive");
  if (selectionRule != "exposure" && selectionRule != "os" &&
      selectionRule != "perfect")
    throw std::invalid_argument(
        "selectionRule must be 'exposure', 'os', or 'perfect'");

  std::vector<uint64_t> seeds(ntrial);
  boost::random::mt19937_64 master(static_cast<uint64_t>(seed));
  for (size_t iteration = 0; iteration < ntrial; ++iteration)
    seeds[iteration] = master();

  std::vector<TrialResult> trials(ntrial);
  SimWorker worker(medianSurvival, meanExposure, exposureSD, rho,
                   interimPatients, interimFollowup, totalPatients,
                   finalAnalysisTime, selectionRule, seeds, &trials);
  RcppParallel::parallelFor(0, ntrial, worker);

  double followupEvents1 = 0.0;
  double followupEvents2 = 0.0;
  double patientEvents1 = 0.0;
  double patientEvents2 = 0.0;
  size_t completed = 0;
  for (const TrialResult &trial : trials) {
    if (!trial.completed)
      continue;
    ++completed;
    followupEvents1 += trial.followupEvents1;
    followupEvents2 += trial.followupEvents2;
    patientEvents1 += trial.patientEvents1;
    patientEvents2 += trial.patientEvents2;
  }
  if (completed == 0)
    throw std::runtime_error("no simulation iterations completed");

  auto weights = [](double events1, double events2) {
    const double total = events1 + events2;
    if (total <= 0.0)
      return std::pair<double, double>{std::sqrt(0.5), std::sqrt(0.5)};
    return std::pair<double, double>{std::sqrt(events1 / total),
                                     std::sqrt(events2 / total)};
  };
  const auto followupWeights = weights(followupEvents1, followupEvents2);
  const auto patientWeights = weights(patientEvents1, patientEvents2);

  std::vector<int> rejected(NMETHOD, 0);
  int selectedHigh = 0;
  int correctSelections = 0;
  double biasSum = 0.0;
  std::vector<int> iterationNumber, selectedRegimen, correctSelection;
  std::vector<double> meanExposureHigh, meanExposureLow, finalZ,
      estimatedHazardRatio;
  std::vector<std::vector<int>> rejectRows(NMETHOD);
  for (size_t iteration = 0; iteration < ntrial; ++iteration) {
    const TrialResult &trial = trials[iteration];
    if (!trial.completed)
      continue;
    const double followupElementary =
        combination_p(trial.pFollowupElementary1, trial.pFollowup2,
                      followupWeights.first, followupWeights.second);
    const double followupIntersection =
        combination_p(trial.pFollowupIntersection1, trial.pFollowup2,
                      followupWeights.first, followupWeights.second);
    const double patientElementary =
        combination_p(trial.pPatientElementary1, trial.pPatient2,
                      patientWeights.first, patientWeights.second);
    const double patientIntersection =
        combination_p(trial.pPatientIntersection1, trial.pPatient2,
                      patientWeights.first, patientWeights.second);
    const int reject[NMETHOD] = {followupElementary < ALPHA_ONE_SIDED &&
                                     followupIntersection < ALPHA_ONE_SIDED,
                                 patientElementary < ALPHA_ONE_SIDED &&
                                     patientIntersection < ALPHA_ONE_SIDED,
                                 trial.pDunnett < ALPHA_ONE_SIDED,
                                 trial.pStage2Only < ALPHA_ONE_SIDED,
                                 trial.pThreeArm < ALPHA_ONE_SIDED};
    for (size_t method = 0; method < NMETHOD; ++method) {
      rejected[method] += reject[method];
      rejectRows[method].push_back(reject[method]);
    }

    selectedHigh += trial.selected == 1;
    correctSelections += trial.correct;
    const double trueHazardRatio =
        medianSurvival[2] / medianSurvival[trial.selected - 1];
    biasSum += trial.estimatedHazardRatio - trueHazardRatio;
    iterationNumber.push_back(static_cast<int>(iteration + 1));
    selectedRegimen.push_back(trial.selected);
    correctSelection.push_back(trial.correct);
    meanExposureHigh.push_back(trial.meanExposureHigh);
    meanExposureLow.push_back(trial.meanExposureLow);
    finalZ.push_back(trial.finalZ);
    estimatedHazardRatio.push_back(trial.estimatedHazardRatio);
  }

  ListCpp byMethod;
  std::vector<std::string> methodNames;
  for (size_t method = 0; method < NMETHOD; ++method) {
    ListCpp summary;
    summary.push_back(static_cast<double>(rejected[method]) / completed,
                      "rejectionProbability");
    methodNames.push_back(METHOD_NAME[method]);
    byMethod.push_back(std::move(summary), METHOD_NAME[method]);
  }

  DataFrameCpp sumdata;
  sumdata.push_back(std::move(iterationNumber), "iterationNumber");
  sumdata.push_back(std::move(selectedRegimen), "selectedRegimen");
  sumdata.push_back(std::move(correctSelection), "correctSelection");
  sumdata.push_back(std::move(meanExposureHigh), "meanExposureHigh");
  sumdata.push_back(std::move(meanExposureLow), "meanExposureLow");
  sumdata.push_back(std::move(finalZ), "finalLogRankZ");
  sumdata.push_back(std::move(estimatedHazardRatio), "estimatedHazardRatio");
  for (size_t method = 0; method < NMETHOD; ++method)
    sumdata.push_back(std::move(rejectRows[method]),
                      std::string("reject.") + METHOD_NAME[method]);

  ListCpp result;
  result.push_back(interimPatients, "interimPatients");
  result.push_back(interimPatients, "stage1Patients");
  result.push_back(totalPatients, "totalPatients");
  result.push_back(finalAnalysisTime, "finalAnalysisTime");
  result.push_back(selectionRule, "selectionRule");
  result.push_back(completed, "numberOfIterations");
  result.push_back(static_cast<double>(selectedHigh) / completed,
                   "probabilitySelectHigh");
  result.push_back(static_cast<double>(correctSelections) / completed,
                   "probabilityCorrectSelection");
  result.push_back(biasSum / completed, "selectionBias");
  result.push_back(
      std::vector<double>{followupWeights.first, followupWeights.second},
      "followupWiseWeights");
  result.push_back(
      std::vector<double>{patientWeights.first, patientWeights.second},
      "patientWiseWeights");
  result.push_back(std::move(methodNames), "methods");
  result.push_back(std::move(byMethod), "byMethod");
  result.push_back(std::move(sumdata), "sumdata");
  return result;
}

// [[Rcpp::export]]
Rcpp::List
lrsim_pkTrtSel_Rcpp(const Rcpp::NumericVector &medianSurvival,
                    const Rcpp::NumericVector &meanExposure,
                    const double exposureSD, const double corrSurvivalExposure,
                    const int interimPatients, const double interimFollowup,
                    const int totalPatients, const double finalAnalysisTime,
                    const std::string &selectionRule,
                    const int maxNumberOfIterations, const int seed) {
  if (maxNumberOfIterations <= 0)
    throw std::invalid_argument("maxNumberOfIterations must be positive");
  auto output = lrsim_pkTrtSel_cpp(
      Rcpp::as<std::vector<double>>(medianSurvival),
      Rcpp::as<std::vector<double>>(meanExposure), exposureSD,
      corrSurvivalExposure, interimPatients, interimFollowup, totalPatients,
      finalAnalysisTime, selectionRule,
      static_cast<size_t>(maxNumberOfIterations), seed);
  thread_utils::drain_thread_warnings_to_R();
  Rcpp::List result = Rcpp::wrap(output);
  result.attr("class") = "lrsim_pkTrtSel";
  return result;
}