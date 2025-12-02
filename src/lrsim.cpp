#include "utilities.h"

using namespace Rcpp;


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
//'   The seed from the environment will be used if left unspecified.
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
//'     - \code{stopStage}: The stage at which the trial stops.
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
List lrsim(const int kMax = 1,
           const NumericVector& informationRates = NA_REAL,
           const NumericVector& criticalValues = NA_REAL,
           const NumericVector& futilityBounds = NA_REAL,
           const double hazardRatioH0 = 1,
           const int allocation1 = 1,
           const int allocation2 = 1,
           const NumericVector& accrualTime = 0,
           const NumericVector& accrualIntensity = NA_REAL,
           const NumericVector& piecewiseSurvivalTime = 0,
           const NumericVector& stratumFraction = 1,
           const NumericVector& lambda1 = NA_REAL,
           const NumericVector& lambda2 = NA_REAL,
           const NumericVector& gamma1 = 0,
           const NumericVector& gamma2 = 0,
           const int n = NA_INTEGER,
           const double followupTime = NA_REAL,
           const bool fixedFollowup = 0,
           const double rho1 = 0,
           const double rho2 = 0,
           const IntegerVector& plannedEvents = NA_INTEGER,
           const NumericVector& plannedTime = NA_REAL,
           const int maxNumberOfIterations = 1000,
           const int maxNumberOfRawDatasetsPerStage = 0,
           const int seed = NA_INTEGER) {

  int nstrata = static_cast<int>(stratumFraction.size());
  int nints = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nints;
  NumericVector lambda1x(nsi), lambda2x(nsi);
  NumericVector gamma1x(nsi), gamma2x(nsi);

  bool useEvents, eventsNotAchieved;
  NumericVector informationRates1 = clone(informationRates);
  NumericVector futilityBounds1 = clone(futilityBounds);


  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  // whether to plan the analyses based on events or calendar time
  if (is_false(any(is_na(plannedEvents)))) {
    useEvents = 1;
    if (plannedEvents[0] <= 0) {
      stop("Elements of plannedEvents must be positive");
    }

    if (plannedEvents.size() != kMax) {
      stop("Invalid length for plannedEvents");
    }

    if (kMax > 1 && is_true(any(diff(plannedEvents) <= 0))) {
      stop("Elements of plannedEvents must be increasing");
    }
  } else if (is_false(any(is_na(plannedTime)))) {
    useEvents = 0;
    if (plannedTime[0] <= 0) {
      stop("Elements of plannedTime must be positive");
    }

    if (plannedTime.size() != kMax) {
      stop("Invalid length for plannedTime");
    }

    if (kMax > 1 && is_true(any(diff(plannedTime) <= 0))) {
      stop("Elements of plannedTime must be increasing");
    }
  } else {
    stop("Either plannedEvents or plannedTime must be given");
  }


  // set default informationRates
  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else if (useEvents) {
    informationRates1 = NumericVector(plannedEvents)/
      (plannedEvents[kMax-1]+0.0);
  } else {
    informationRates1 = plannedTime/plannedTime[kMax-1];
  }


  if (is_true(any(is_na(criticalValues)))) {
    stop("criticalValues must be provided");
  }

  if (criticalValues.size() != kMax) {
    stop("Invalid length for criticalValues");
  }


  if (kMax > 1 && is_true(any(is_na(futilityBounds)))) {
    futilityBounds1 = rep(-6.0, kMax-1);
  }

  if (is_false(any(is_na(futilityBounds1)))) {
    if (futilityBounds1.size() < kMax-1) {
      stop("Invalid length for futilityBounds");
    }
  }

  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds1)))) {
    for (int i=0; i<kMax-1; ++i) {
      if (futilityBounds1[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }
  }

  if (hazardRatioH0 <= 0) {
    stop("hazardRatioH0 must be positive");
  }

  if (allocation1 < 1) {
    stop("allocation1 must be a positive integer");
  }

  if (allocation2 < 1) {
    stop("allocation2 must be a positive integer");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nints > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }

  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }


  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nsi);
  } else if (lambda1.size() == nints) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nints) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }

  if (gamma1.size() == 1) {
    gamma1x = rep(gamma1, nsi);
  } else if (gamma1.size() == nints) {
    gamma1x = rep(gamma1, nstrata);
  } else if (gamma1.size() == nsi) {
    gamma1x = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() == 1) {
    gamma2x = rep(gamma2, nsi);
  } else if (gamma2.size() == nints) {
    gamma2x = rep(gamma2, nstrata);
  } else if (gamma2.size() == nsi) {
    gamma2x = gamma2;
  } else {
    stop("Invalid length for gamma2");
  }

  if (n == NA_INTEGER) {
    stop("n must be provided");
  }

  if (n <= 0) {
    stop("n must be positive");
  }

  if (fixedFollowup) {
    if (R_isnancpp(followupTime)) {
      stop("followupTime must be provided for fixed follow-up");
    }

    if (followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }
  }

  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }

  if (rho2 < 0) {
    stop("rho2 must be non-negative");
  }

  if (maxNumberOfIterations < 1) {
    stop("maxNumberOfIterations must be a positive integer");
  }

  if (maxNumberOfRawDatasetsPerStage < 0) {
    stop("maxNumberOfRawDatasetsPerStage must be a non-negative integer");
  }


  // declare variables
  int nevents, nstages, stopStage;
  int index1=0, index2=0;

  double u, enrollt, time, uscore1, vscore1;

  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n);

  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n),
  timeUnderObservation(n), totalTime(n), totalt(n);

  LogicalVector event(n), dropoutEvent(n);


  // stratum information
  IntegerVector b1(nstrata), b2(nstrata), n1(nstrata), n2(nstrata);
  NumericVector nt(nstrata), n1a(nstrata), nta(nstrata);
  NumericVector km(nstrata), w(nstrata);
  NumericVector cumStratumFraction = cumsum(stratumFraction);


  // within-stratum hazard rates
  NumericVector lam1(nints), lam2(nints);
  NumericVector gam1(nints), gam2(nints);

  // stage-wise information
  IntegerVector accruals1(kMax), accruals2(kMax), totalAccruals(kMax),
  events1(kMax), events2(kMax), totalEvents(kMax),
  dropouts1(kMax), dropouts2(kMax), totalDropouts(kMax),
  niter(kMax), obsEvents(kMax);

  NumericVector analysisTime(kMax), uscore(kMax), vscore(kMax),
  lrstat(kMax), adjCriticalValues(kMax);

  LogicalVector rejectPerStage(kMax), futilityPerStage(kMax);


  // cache for the patient-level raw data to extract
  int nrow1 = std::min(n*kMax*maxNumberOfRawDatasetsPerStage,
                       n*maxNumberOfIterations);

  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stopStagex(nrow1);
  NumericVector analysisTimex(nrow1);
  IntegerVector subjectIdx(nrow1);
  NumericVector arrivalTimex(nrow1);
  IntegerVector stratumx(nrow1);
  IntegerVector treatmentGroupx(nrow1);
  NumericVector survivalTimex(nrow1);
  NumericVector dropoutTimex(nrow1);
  NumericVector timeUnderObservationx(nrow1);
  LogicalVector eventx(nrow1);
  LogicalVector dropoutEventx(nrow1);

  // cache for the simulation-level summary data to extract
  int nrow2 = kMax*maxNumberOfIterations;

  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector stopStagey(nrow2);
  LogicalVector eventsNotAchievedy(nrow2);
  IntegerVector stageNumbery(nrow2);
  NumericVector analysisTimey(nrow2);
  IntegerVector accruals1y(nrow2);
  IntegerVector accruals2y(nrow2);
  IntegerVector totalAccrualsy(nrow2);
  IntegerVector events1y(nrow2);
  IntegerVector events2y(nrow2);
  IntegerVector totalEventsy(nrow2);
  IntegerVector dropouts1y(nrow2);
  IntegerVector dropouts2y(nrow2);
  IntegerVector totalDropoutsy(nrow2);
  NumericVector uscorey(nrow2);
  NumericVector vscorey(nrow2);
  NumericVector logRankStatisticy(nrow2);
  LogicalVector rejectPerStagey(nrow2);
  LogicalVector futilityPerStagey(nrow2);


  // total alpha to adjust the critical value at the final stage
  NumericVector lb(kMax, -6.0);
  NumericVector theta(kMax);
  List p1 = exitprobcpp(criticalValues, lb, theta, informationRates1);
  double alpha = sum(NumericVector(p1[0]));


  // set up random seed
  if (seed != NA_INTEGER) {
    set_seed(seed);
  }


  for (int iter=0; iter<maxNumberOfIterations; ++iter) {
    int nstops = 0;

    b1.fill(allocation1);
    b2.fill(allocation2);

    enrollt = 0;
    for (int i=0; i<n; ++i) {

      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);
      arrivalTime[i] = enrollt;

      // generate stratum information
      u = R::runif(0,1);
      int j;
      for (j=0; j<nstrata; ++j) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }

      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+0.0)) {
        treatmentGroup[i] = 1;
        --b1[j];
      } else {
        treatmentGroup[i] = 2;
        --b2[j];
      }

      // start a new block after depleting the current block
      if (b1[j]+b2[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
      }

      // stratum-specific hazard rates for event and dropout
      Range jj = Range(j*nints, (j+1)*nints-1);

      lam1 = lambda1x[jj];
      lam2 = lambda2x[jj];
      gam1 = gamma1x[jj];
      gam2 = gamma2x[jj];

      // generate survival time
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        survivalTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, lam1, 0, 1, 0);
      } else {
        survivalTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, lam2, 0, 1, 0);
      }

      // generate dropout time
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, gam1, 0, 1, 0);
      } else {
        dropoutTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, gam2, 0, 1, 0);
      }

      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime[i] <= dropoutTime[i] &&
            survivalTime[i] <= followupTime) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else if (dropoutTime[i] <= survivalTime[i] &&
          dropoutTime[i] <= followupTime) {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        } else {
          timeUnderObservation[i] = followupTime;
          event[i] = 0;
          dropoutEvent[i] = 0;
        }
      } else { // variable follow-up design
        if (survivalTime[i] <= dropoutTime[i]) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        }
      }

      totalTime[i] = arrivalTime[i] + timeUnderObservation[i];

    }


    // find the analysis time for each stage
    if (useEvents) {
      nevents = sum(event);
      totalt = stl_sort(totalTime[event]);
      nstages = kMax;

      int j;
      for (j=0; j<kMax; ++j) {
        if (plannedEvents[j] >= nevents) {
          nstages = j+1;
          break;
        }
      }


      if (j==kMax) { // total number of events exceeds planned
        for (int k=0; k<nstages; ++k) {
          analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
          obsEvents[k] = plannedEvents[k];
        }
      } else {
        for (int k=0; k<nstages; ++k) {
          if (k < nstages-1) {
            analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
            obsEvents[k] = plannedEvents[k];
          } else {
            analysisTime[k] = totalt[nevents-1] + 1e-12;
            obsEvents[k] = nevents;
          }
        }
      }

      // observed total number of events less than planned
      eventsNotAchieved = (nevents < plannedEvents[kMax-1]);
    } else {
      nstages = kMax;
      analysisTime = clone(plannedTime);
      eventsNotAchieved = 0;
    }


    // construct the log-rank test statistic at each stage
    stopStage = nstages;
    for (int k=0; k<nstages; ++k) {
      time = analysisTime[k];

      n1.fill(0);  // number of subjects in each stratum by treatment
      n2.fill(0);
      events1[k] = 0;
      events2[k] = 0;
      dropouts1[k] = 0;
      dropouts2[k] = 0;

      // censor at analysis time
      for (int i=0; i<n; ++i) {
        int h = stratum[i]-1;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation[i] = time - arrivalTime[i];
          event[i] = 0;
          dropoutEvent[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            ++n1[h];
          } else if (treatmentGroup[i]==2) {
            ++n2[h];
          }

          if (fixedFollowup) {
            // the first three cases correspond to arrivalTime[i] +
            // min(survivalTime[i], dropoutTime[i], followupTime) <= time
            if (arrivalTime[i] + survivalTime[i] <= time &&
                survivalTime[i] <= dropoutTime[i] &&
                survivalTime[i] <= followupTime) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] <= time &&
              dropoutTime[i] <= survivalTime[i] &&
              dropoutTime[i] <= followupTime) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else if (arrivalTime[i] + followupTime <= time &&
              followupTime <= survivalTime[i] &&
              followupTime <= dropoutTime[i]) {
              timeUnderObservation[i] = followupTime;
              event[i] = 0;
              dropoutEvent[i] = 0;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime[i] <= time &&
                survivalTime[i] <= dropoutTime[i]) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] <= time &&
              dropoutTime[i] <= survivalTime[i]) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          }

          if (treatmentGroup[i]==1 && event[i]) ++events1[k];
          if (treatmentGroup[i]==2 && event[i]) ++events2[k];
          if (treatmentGroup[i]==1 && dropoutEvent[i]) ++dropouts1[k];
          if (treatmentGroup[i]==2 && dropoutEvent[i]) ++dropouts2[k];
        }
      }

      // number of accrued patients and total number of events
      accruals1[k] = sum(n1);
      accruals2[k] = sum(n2);
      totalAccruals[k] = accruals1[k] + accruals2[k];

      totalEvents[k] = events1[k] + events2[k];
      totalDropouts[k] = dropouts1[k] + dropouts2[k];

      // order the data by time under observation
      NumericVector timeUnderObservationSorted =
        stl_sort(timeUnderObservation);
      IntegerVector sortedIndex = match(timeUnderObservationSorted,
                                        timeUnderObservation);
      sortedIndex = sortedIndex - 1;
      IntegerVector stratumSorted = stratum[sortedIndex];
      IntegerVector treatmentGroupSorted = treatmentGroup[sortedIndex];
      LogicalVector eventSorted = event[sortedIndex];

      LogicalVector subSorted = (timeUnderObservationSorted > 0);
      stratumSorted = stratumSorted[subSorted];
      treatmentGroupSorted = treatmentGroupSorted[subSorted];
      eventSorted = eventSorted[subSorted];
      int nsubSorted = static_cast<int>(eventSorted.size());

      // calculate the stratified log-rank test
      uscore1 = 0;
      vscore1 = 0;
      km.fill(1);  // km(t-) estimate by stratum
      for (int i=0; i<nsubSorted; ++i) {
        int h = stratumSorted[i] - 1;
        n1a[h] = n1[h]*hazardRatioH0;
        nt[h] = n1[h] + n2[h];
        nta[h] = n1a[h] + n2[h];

        if (eventSorted[i]) { // at most 1 event can occur at a given time
          w[h] = pow(km[h], rho1)*pow(1-km[h], rho2);
          uscore1 += w[h]*((treatmentGroupSorted[i]==1)-n1a[h]/nta[h]);
          vscore1 += w[h]*w[h]*n1a[h]*n2[h]/(nta[h]*nta[h]);
          km[h] *= (1-1/nt[h]); // update km estimate
        }

        // reduce the risk set
        if (treatmentGroupSorted[i]==1) {
          --n1[h];
        } else {
          --n2[h];
        }
      }

      uscore[k] = uscore1;
      vscore[k] = vscore1;

      // log-rank z statistic
      lrstat[k] = uscore1/sqrt(vscore1);

      if (useEvents) {
        // adjust the critical value at the final stage if the planned total
        // number of events is not achieved
        if (k < nstages-1 || !eventsNotAchieved) {
          // no change to the critical
          // values at earlier stages, or at the final stage if the planned
          // total number of events is achieved (the number of stages is also
          // the same as planned in this case)
          adjCriticalValues[k] = criticalValues[k];
        } else { // assign all remaining alpha to the final stage
          if (rho1 == 0 && rho2 == 0) { // conventional log-rank test
            auto f = [criticalValues, alpha, &obsEvents,
                      &nstages](double aval)->double {
                        NumericVector u(nstages);
                        for (int i=0; i<nstages-1; ++i) {
                          u[i] = criticalValues[i];
                        }
                        u[nstages-1] = aval;
                        NumericVector l = rep(-6.0, nstages);
                        NumericVector theta = rep(0.0, nstages);
                        NumericVector I = NumericVector(obsEvents)[
                        Range(0,nstages-1)];
                        List p2 = exitprobcpp(u, l, theta, I);
                        return sum(NumericVector(p2[0])) - alpha;
                      };

            adjCriticalValues[nstages-1] = brent(f, 0, 6, 1e-6);
          } else { // weighted log-rank test
            auto f = [criticalValues, alpha, &vscore,
                      &nstages](double aval)->double {
                        NumericVector u(nstages);
                        for (int i=0; i<nstages-1; ++i) {
                          u[i] = criticalValues[i];
                        }
                        u[nstages-1] = aval;
                        NumericVector l = rep(-6.0, nstages);
                        NumericVector theta = rep(0.0, nstages);
                        NumericVector I = vscore[Range(0,nstages-1)];
                        List p2 = exitprobcpp(u, l, theta, I);
                        return sum(NumericVector(p2[0])) - alpha;
                      };

            adjCriticalValues[nstages-1] = brent(f, 0, 6, 1e-6);
          }
        }

      } else {
        adjCriticalValues[k] = criticalValues[k];
      }


      // compare to the critical values to make decisions
      rejectPerStage[k] = 0;
      futilityPerStage[k] = 0;
      if (-lrstat[k] > adjCriticalValues[k]) {
        rejectPerStage[k] = 1;
      } else if ((k < nstages-1 && -lrstat[k] < futilityBounds1[k])
                   || (k == nstages-1))  {
        futilityPerStage[k] = 1;
      }


      if (rejectPerStage[k]==1 || futilityPerStage[k]==1) {
        ++nstops;

        if (nstops == 1) { // extract at most one raw data set per iteration

          // add raw data to output
          if (niter[k] < maxNumberOfRawDatasetsPerStage) {
            for (int i=0; i<n; ++i) {
              iterationNumberx[index1] = iter+1;
              stopStagex[index1] = k+1;
              analysisTimex[index1] = time;
              subjectIdx[index1] = i+1;
              arrivalTimex[index1] = arrivalTime[i];
              stratumx[index1] = stratum[i];
              treatmentGroupx[index1] = treatmentGroup[i];
              survivalTimex[index1] = survivalTime[i];
              dropoutTimex[index1] = dropoutTime[i];
              timeUnderObservationx[index1] = timeUnderObservation[i];
              eventx[index1] = event[i];
              dropoutEventx[index1] = dropoutEvent[i];
              ++index1;
            }

            // update the number of stage k dataset to extract
            ++niter[k];
          }

          stopStage = k+1;

        }

      }

    }

    // add summary data to output
    for (int k=0; k<nstages; ++k) {
      iterationNumbery[index2] = iter+1;
      stopStagey[index2] = stopStage;
      eventsNotAchievedy[index2] = eventsNotAchieved;
      stageNumbery[index2] = k+1;
      analysisTimey[index2] = analysisTime[k];
      accruals1y[index2] = accruals1[k];
      accruals2y[index2] = accruals2[k];
      totalAccrualsy[index2] = totalAccruals[k];
      events1y[index2] = events1[k];
      events2y[index2] = events2[k];
      totalEventsy[index2] = totalEvents[k];
      dropouts1y[index2] = dropouts1[k];
      dropouts2y[index2] = dropouts2[k];
      totalDropoutsy[index2] = totalDropouts[k];
      uscorey[index2] = uscore[k];
      vscorey[index2] = vscore[k];
      logRankStatisticy[index2] = lrstat[k];
      rejectPerStagey[index2] = rejectPerStage[k];
      futilityPerStagey[index2] = futilityPerStage[k];
      ++index2;
    }


  }

  // only keep nonmissing records
  LogicalVector sub2 = !is_na(iterationNumbery);
  iterationNumbery = iterationNumbery[sub2];
  stopStagey = stopStagey[sub2];
  eventsNotAchievedy = eventsNotAchievedy[sub2];
  stageNumbery = stageNumbery[sub2];
  analysisTimey = analysisTimey[sub2];
  accruals1y = accruals1y[sub2];
  accruals2y = accruals2y[sub2];
  totalAccrualsy = totalAccrualsy[sub2];
  events1y = events1y[sub2];
  events2y = events2y[sub2];
  totalEventsy = totalEventsy[sub2];
  dropouts1y = dropouts1y[sub2];
  dropouts2y = dropouts2y[sub2];
  totalDropoutsy = totalDropoutsy[sub2];
  uscorey = uscorey[sub2];
  vscorey = vscorey[sub2];
  logRankStatisticy = logRankStatisticy[sub2];
  rejectPerStagey = rejectPerStagey[sub2];
  futilityPerStagey = futilityPerStagey[sub2];



  // simulation results on power and expected sample size

  NumericVector pRejectPerStage(kMax), pFutilityPerStage(kMax),
  nEventsPerStage(kMax), nDropoutsPerStage(kMax), nSubjectsPerStage(kMax),
  analysisTimePerStage(kMax);


  // number of observations in the summary dataset
  int nrow3 = static_cast<int>(stageNumbery.size());

  for (int i=0; i<nrow3; ++i) {
    int k = stageNumbery[i] - 1;
    if (stageNumbery[i] == stopStagey[i]) {
      pRejectPerStage[k] += rejectPerStagey[i];
      pFutilityPerStage[k] += futilityPerStagey[i];
    }

    nEventsPerStage[k] += totalEventsy[i];
    nDropoutsPerStage[k] += totalDropoutsy[i];
    nSubjectsPerStage[k] += totalAccrualsy[i];
    analysisTimePerStage[k] += analysisTimey[i];
  }


  for (int k=0; k<kMax; ++k) {
    pRejectPerStage[k] /= maxNumberOfIterations;
    pFutilityPerStage[k] /= maxNumberOfIterations;
    nEventsPerStage[k] /= maxNumberOfIterations;
    nDropoutsPerStage[k] /= maxNumberOfIterations;
    nSubjectsPerStage[k] /= maxNumberOfIterations;
    analysisTimePerStage[k] /= maxNumberOfIterations;
  }

  NumericVector cpu = cumsum(pRejectPerStage);
  NumericVector cpl = cumsum(pFutilityPerStage);

  double pOverallReject = sum(pRejectPerStage);

  double expectedNumberOfEvents=0, expectedNumberOfDropouts=0,
    expectedNumberOfSubjects=0, expectedStudyDuration=0;

  for (int i=0; i<nrow3; ++i) {
    if (stageNumbery[i] == stopStagey[i]) {
      expectedNumberOfEvents += totalEventsy[i];
      expectedNumberOfDropouts += totalDropoutsy[i];
      expectedNumberOfSubjects += totalAccrualsy[i];
      expectedStudyDuration += analysisTimey[i];
    }
  }

  expectedNumberOfEvents /= maxNumberOfIterations;
  expectedNumberOfDropouts /= maxNumberOfIterations;
  expectedNumberOfSubjects /= maxNumberOfIterations;
  expectedStudyDuration /= maxNumberOfIterations;

  List overview = List::create(
    _["rejectPerStage"] = pRejectPerStage,
    _["futilityPerStage"] = pFutilityPerStage,
    _["cumulativeRejection"] = cpu,
    _["cumulativeFutility"] = cpl,
    _["numberOfEvents"] = nEventsPerStage,
    _["numberOfDropouts"] = nDropoutsPerStage,
    _["numberOfSubjects"] = nSubjectsPerStage,
    _["analysisTime"] = analysisTimePerStage,
    _["overallReject"] = pOverallReject,
    _["expectedNumberOfEvents"] = expectedNumberOfEvents,
    _["expectedNumberOfDropouts"] = expectedNumberOfDropouts,
    _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
    _["expectedStudyDuration"] = expectedStudyDuration,
    _["hazardRatioH0"] = hazardRatioH0,
    _["useEvents"] = useEvents,
    _["n"] = n,
    _["fixedFollowup"] = fixedFollowup,
    _["rho1"] = rho1,
    _["rho2"] = rho2,
    _["kMax"] = kMax);



  // simulation datasets
  DataFrame sumdata = DataFrame::create(
    _["iterationNumber"] = iterationNumbery,
    _["stopStage"] = stopStagey,
    _["eventsNotAchieved"] = eventsNotAchievedy,
    _["stageNumber"] = stageNumbery,
    _["analysisTime"] = analysisTimey,
    _["accruals1"] = accruals1y,
    _["accruals2"] = accruals2y,
    _["totalAccruals"] = totalAccrualsy,
    _["events1"] = events1y,
    _["events2"] = events2y,
    _["totalEvents"] = totalEventsy,
    _["dropouts1"] = dropouts1y,
    _["dropouts2"] = dropouts2y,
    _["totalDropouts"] = totalDropoutsy,
    _["uscore"] = uscorey,
    _["vscore"] = vscorey,
    _["logRankStatistic"] = logRankStatisticy,
    _["rejectPerStage"] = rejectPerStagey,
    _["futilityPerStage"] = futilityPerStagey);


  List result;

  if (maxNumberOfRawDatasetsPerStage > 0) {
    LogicalVector sub1 = !is_na(iterationNumberx);
    iterationNumberx = iterationNumberx[sub1];
    stopStagex = stopStagex[sub1];
    analysisTimex = analysisTimex[sub1];
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    survivalTimex = survivalTimex[sub1];
    dropoutTimex = dropoutTimex[sub1];
    timeUnderObservationx = timeUnderObservationx[sub1];
    eventx = eventx[sub1];
    dropoutEventx = dropoutEventx[sub1];

    DataFrame rawdata = DataFrame::create(
      _["iterationNumber"] = iterationNumberx,
      _["stopStage"] = stopStagex,
      _["analysisTime"] = analysisTimex,
      _["subjectId"] = subjectIdx,
      _["arrivalTime"] = arrivalTimex,
      _["stratum"] = stratumx,
      _["treatmentGroup"] = treatmentGroupx,
      _["survivalTime"] = survivalTimex,
      _["dropoutTime"] = dropoutTimex,
      _["timeUnderObservation"] = timeUnderObservationx,
      _["event"] = eventx,
      _["dropoutEvent"] = dropoutEventx);

    result = List::create(_["overview"] = overview,
                          _["sumdata"] = sumdata,
                          _["rawdata"] = rawdata);
  } else {
    result = List::create(_["overview"] = overview,
                          _["sumdata"] = sumdata);
  }

  result.attr("class") = "lrsim";


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
//'   The seed from the environment will be used if left unspecified.
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
//'     - \code{logRankStatistic13}: The log-rank test Z-statistic
//'       comparing the active treatment 1 to the control.
//'
//'     - \code{logRankStatistic23}: The log-rank test Z-statistic
//'       comparing the active treatment 2 to the control.
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
List lrsim3a(const int kMax = 1,
             const double hazardRatioH013 = 1,
             const double hazardRatioH023 = 1,
             const double hazardRatioH012 = 1,
             const int allocation1 = 1,
             const int allocation2 = 1,
             const int allocation3 = 1,
             const NumericVector& accrualTime = 0,
             const NumericVector& accrualIntensity = NA_REAL,
             const NumericVector& piecewiseSurvivalTime = 0,
             const NumericVector& stratumFraction = 1,
             const NumericVector& lambda1 = NA_REAL,
             const NumericVector& lambda2 = NA_REAL,
             const NumericVector& lambda3 = NA_REAL,
             const NumericVector& gamma1 = 0,
             const NumericVector& gamma2 = 0,
             const NumericVector& gamma3 = 0,
             const int n = NA_INTEGER,
             const double followupTime = NA_REAL,
             const bool fixedFollowup = 0,
             const double rho1 = 0,
             const double rho2 = 0,
             const IntegerVector& plannedEvents = NA_INTEGER,
             const NumericVector& plannedTime = NA_REAL,
             const int maxNumberOfIterations = 1000,
             const int maxNumberOfRawDatasetsPerStage = 0,
             const int seed = NA_INTEGER) {

  // check input parameters
  int nstrata = static_cast<int>(stratumFraction.size());
  int nints = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nints;

  NumericVector lambda1x(nsi), lambda2x(nsi), lambda3x(nsi);
  NumericVector gamma1x(nsi), gamma2x(nsi), gamma3x(nsi);

  bool useEvents, eventsNotAchieved;


  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  // whether to plan the analyses based on events or calendar time
  if (is_false(any(is_na(plannedEvents)))) {
    useEvents = 1;
    if (plannedEvents[0] <= 0) {
      stop("Elements of plannedEvents must be positive");
    }

    if (plannedEvents.size() != kMax) {
      stop("Invalid length for plannedEvents");
    }

    if (kMax > 1 && is_true(any(diff(plannedEvents) <= 0))) {
      stop("Elements of plannedEvents must be increasing");
    }
  } else if (is_false(any(is_na(plannedTime)))) {
    useEvents = 0;
    if (plannedTime[0] <= 0) {
      stop("Elements of plannedTime must be positive");
    }

    if (plannedTime.size() != kMax) {
      stop("Invalid length for plannedTime");
    }

    if (kMax > 1 && is_true(any(diff(plannedTime) <= 0))) {
      stop("Elements of plannedTime must be increasing");
    }
  } else {
    stop("Either plannedEvents or plannedTime must be given");
  }

  if (hazardRatioH013 <= 0) {
    stop("hazardRatioH013 must be positive");
  }

  if (hazardRatioH023 <= 0) {
    stop("hazardRatioH023 must be positive");
  }

  if (hazardRatioH012 <= 0) {
    stop("hazardRatioH012 must be positive");
  }

  if (allocation1 < 1) {
    stop("allocation1 must be a positive integer");
  }

  if (allocation2 < 1) {
    stop("allocation2 must be a positive integer");
  }

  if (allocation3 < 1) {
    stop("allocation3 must be a positive integer");
  }


  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }


  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nints > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }


  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }

  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(is_na(lambda3)))) {
    stop("lambda3 must be provided");
  }

  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }

  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }

  if (is_true(any(lambda3 < 0))) {
    stop("lambda3 must be non-negative");
  }


  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (is_true(any(gamma3 < 0))) {
    stop("gamma3 must be non-negative");
  }



  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nsi);
  } else if (lambda1.size() == nints) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }


  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nints) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }


  if (lambda3.size() == 1) {
    lambda3x = rep(lambda3, nsi);
  } else if (lambda3.size() == nints) {
    lambda3x = rep(lambda3, nstrata);
  } else if (lambda3.size() == nsi) {
    lambda3x = lambda3;
  } else {
    stop("Invalid length for lambda3");
  }



  if (gamma1.size() == 1) {
    gamma1x = rep(gamma1, nsi);
  } else if (gamma1.size() == nints) {
    gamma1x = rep(gamma1, nstrata);
  } else if (gamma1.size() == nsi) {
    gamma1x = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }


  if (gamma2.size() == 1) {
    gamma2x = rep(gamma2, nsi);
  } else if (gamma2.size() == nints) {
    gamma2x = rep(gamma2, nstrata);
  } else if (gamma2.size() == nsi) {
    gamma2x = gamma2;
  } else {
    stop("Invalid length for gamma2");
  }


  if (gamma3.size() == 1) {
    gamma3x = rep(gamma3, nsi);
  } else if (gamma3.size() == nints) {
    gamma3x = rep(gamma3, nstrata);
  } else if (gamma3.size() == nsi) {
    gamma3x = gamma3;
  } else {
    stop("Invalid length for gamma3");
  }


  if (n == NA_INTEGER) {
    stop("n must be provided");
  }

  if (n <= 0) {
    stop("n must be positive");
  }

  if (fixedFollowup) {
    if (R_isnancpp(followupTime)) {
      stop("followupTime must be provided for fixed follow-up");
    }

    if (followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }
  }

  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }

  if (rho2 < 0) {
    stop("rho2 must be non-negative");
  }

  if (maxNumberOfIterations < 1) {
    stop("maxNumberOfIterations must be a positive integer");
  }

  if (maxNumberOfRawDatasetsPerStage < 0) {
    stop("maxNumberOfRawDatasetsPerStage must be a non-negative integer");
  }


  // declare variables
  int nevents, nstages;
  int accruals1, accruals2, accruals3, totalAccruals;
  int events1, events2, events3, totalEvents;
  int dropouts1, dropouts2, dropouts3, totalDropouts;
  int index1=0, index2=0;

  double enrollt, u, time;
  double uscore13, uscore23, uscore12;
  double vscore13, vscore23, vscore12;


  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n);

  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n),
  timeUnderObservation(n), totalTime(n), totalt(n);

  LogicalVector event(n), dropoutEvent(n), eventac(n);


  // stratum information
  IntegerVector b1(nstrata), b2(nstrata), b3(nstrata);
  IntegerVector n1(nstrata), n2(nstrata), n3(nstrata);
  NumericVector nt13(nstrata), nt23(nstrata), nt12(nstrata);
  NumericVector n13a(nstrata), n23a(nstrata), n12a(nstrata);
  NumericVector nt13a(nstrata), nt23a(nstrata), nt12a(nstrata);

  NumericVector km13(nstrata), km23(nstrata), km12(nstrata);
  NumericVector w13(nstrata), w23(nstrata), w12(nstrata);
  NumericVector cumStratumFraction = cumsum(stratumFraction);

  // within-stratum hazard rates
  NumericVector lam1(nints), lam2(nints), lam3(nints);
  NumericVector gam1(nints), gam2(nints), gam3(nints);


  // stage-wise information
  IntegerVector niter(kMax);
  NumericVector analysisTime(kMax);


  // cache for the patient-level raw data to extract
  int nrow1 = n*kMax*maxNumberOfRawDatasetsPerStage;

  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stageNumberx(nrow1);
  NumericVector analysisTimex(nrow1);
  IntegerVector subjectIdx(nrow1);
  NumericVector arrivalTimex(nrow1);
  IntegerVector stratumx(nrow1);
  IntegerVector treatmentGroupx(nrow1);
  NumericVector survivalTimex(nrow1);
  NumericVector dropoutTimex(nrow1);
  NumericVector timeUnderObservationx(nrow1);
  LogicalVector eventx(nrow1);
  LogicalVector dropoutEventx(nrow1);

  // cache for the simulation-level summary data to extract
  int nrow2 = kMax*maxNumberOfIterations*2;

  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  LogicalVector eventsNotAchievedy(nrow2);
  IntegerVector stageNumbery(nrow2);
  NumericVector analysisTimey(nrow2);
  IntegerVector accruals1y(nrow2);
  IntegerVector accruals2y(nrow2);
  IntegerVector accruals3y(nrow2);
  IntegerVector totalAccrualsy(nrow2);
  IntegerVector events1y(nrow2);
  IntegerVector events2y(nrow2);
  IntegerVector events3y(nrow2);
  IntegerVector totalEventsy(nrow2);
  IntegerVector dropouts1y(nrow2);
  IntegerVector dropouts2y(nrow2);
  IntegerVector dropouts3y(nrow2);
  IntegerVector totalDropoutsy(nrow2);
  NumericVector logRankStatistic13y(nrow2);
  NumericVector logRankStatistic23y(nrow2);
  NumericVector logRankStatistic12y(nrow2);


  // set up random seed
  if (seed != NA_INTEGER) {
    set_seed(seed);
  }


  for (int iter=0; iter<maxNumberOfIterations; ++iter) {

    b1.fill(allocation1);
    b2.fill(allocation2);
    b3.fill(allocation3);

    enrollt = 0;
    for (int i=0; i<n; ++i) {

      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);
      arrivalTime[i] = enrollt;

      // generate stratum information
      u = R::runif(0,1);
      int j;
      for (j=0; j<nstrata; ++j) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }

      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+b3[j]+0.0)) {
        treatmentGroup[i] = 1;
        --b1[j];
      } else if (u <= (b1[j]+b2[j])/(b1[j]+b2[j]+b3[j]+0.0)) {
        treatmentGroup[i] = 2;
        --b2[j];
      } else {
        treatmentGroup[i] = 3;
        --b3[j];
      }

      // start a new block after depleting the current block
      if (b1[j]+b2[j]+b3[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
        b3[j] = allocation3;
      }

      // stratum-specific hazard rates for event and dropout
      Range jj = Range(j*nints, (j+1)*nints-1);

      lam1 = lambda1x[jj];
      lam2 = lambda2x[jj];
      lam3 = lambda3x[jj];

      gam1 = gamma1x[jj];
      gam2 = gamma2x[jj];
      gam3 = gamma3x[jj];


      // generate survival time
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        survivalTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, lam1, 0, 1, 0);
      } else if (treatmentGroup[i]==2) {
        survivalTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, lam2, 0, 1, 0);
      } else {
        survivalTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, lam3, 0, 1, 0);
      }


      // generate dropout time
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, gam1, 0, 1, 0);
      } else if (treatmentGroup[i]==2) {
        dropoutTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, gam2, 0, 1, 0);
      } else {
        dropoutTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, gam3, 0, 1, 0);
      }


      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime[i] <= dropoutTime[i] &&
            survivalTime[i] <= followupTime) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else if (dropoutTime[i] <= survivalTime[i] &&
          dropoutTime[i] <= followupTime) {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        } else {
          timeUnderObservation[i] = followupTime;
          event[i] = 0;
          dropoutEvent[i] = 0;
        }
      } else { // variable follow-up design
        if (survivalTime[i] <= dropoutTime[i]) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        }
      }

      totalTime[i] = arrivalTime[i] + timeUnderObservation[i];

    }


    // find the analysis time for each stage based on Arm A vs. Arm C
    if (useEvents) {
      eventac = event & ((treatmentGroup==1) | (treatmentGroup==3));
      nevents = sum(eventac);
      totalt = stl_sort(totalTime[eventac]);
      nstages = kMax;

      int j;
      for (j=0; j<kMax; ++j) {
        if (plannedEvents[j] >= nevents) {
          nstages = j+1;
          break;
        }
      }

      if (j==kMax) { // total number of events exceeds planned
        for (int k=0; k<nstages; ++k) {
          analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
        }
      } else {
        for (int k=0; k<nstages; ++k) {
          if (k < nstages-1) {
            analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
          } else {
            analysisTime[k] = totalt[nevents-1] + 1e-12;
          }
        }
      }

      // observed total number of events less than planned
      eventsNotAchieved = (nevents < plannedEvents[kMax-1]);
    } else {
      nstages = kMax;
      analysisTime = clone(plannedTime);
      eventsNotAchieved = 0;
    }


    // construct the log-rank test statistic at each stage
    for (int k=0; k<nstages; ++k) {
      time = analysisTime[k];

      n1.fill(0);  // number of subjects in each stratum by treatment
      n2.fill(0);
      n3.fill(0);

      events1 = 0;
      events2 = 0;
      events3 = 0;

      dropouts1 = 0;
      dropouts2 = 0;
      dropouts3 = 0;

      // censor at analysis time
      for (int i=0; i<n; ++i) {
        int h = stratum[i]-1;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation[i] = time - arrivalTime[i];
          event[i] = 0;
          dropoutEvent[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            ++n1[h];
          } else if (treatmentGroup[i]==2) {
            ++n2[h];
          } else {
            ++n3[h];
          }

          if (fixedFollowup) {
            // the first three cases correspond to arrivalTime[i] +
            // min(survivalTime[i], dropoutTime[i], followupTime) <= time
            if (arrivalTime[i] + survivalTime[i] <= time &&
                survivalTime[i] <= dropoutTime[i] &&
                survivalTime[i] <= followupTime) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] <= time &&
              dropoutTime[i] <= survivalTime[i] &&
              dropoutTime[i] <= followupTime) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else if (arrivalTime[i] + followupTime <= time &&
              followupTime <= survivalTime[i] &&
              followupTime <= dropoutTime[i]) {
              timeUnderObservation[i] = followupTime;
              event[i] = 0;
              dropoutEvent[i] = 0;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime[i] <= time &&
                survivalTime[i] <= dropoutTime[i]) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] <= time &&
              dropoutTime[i] <= survivalTime[i]) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          }

          if (treatmentGroup[i]==1 && event[i]) ++events1;
          if (treatmentGroup[i]==2 && event[i]) ++events2;
          if (treatmentGroup[i]==3 && event[i]) ++events3;
          if (treatmentGroup[i]==1 && dropoutEvent[i]) ++dropouts1;
          if (treatmentGroup[i]==2 && dropoutEvent[i]) ++dropouts2;
          if (treatmentGroup[i]==3 && dropoutEvent[i]) ++dropouts3;
        }
      }


      // add raw data to output
      if (niter[k] < maxNumberOfRawDatasetsPerStage) {
        for (int i=0; i<n; ++i) {
          iterationNumberx[index1] = iter+1;
          stageNumberx[index1] = k+1;
          analysisTimex[index1] = time;
          subjectIdx[index1] = i+1;
          arrivalTimex[index1] = arrivalTime[i];
          stratumx[index1] = stratum[i];
          treatmentGroupx[index1] = treatmentGroup[i];
          survivalTimex[index1] = survivalTime[i];
          dropoutTimex[index1] = dropoutTime[i];
          timeUnderObservationx[index1] = timeUnderObservation[i];
          eventx[index1] = event[i];
          dropoutEventx[index1] = dropoutEvent[i];

          ++index1;
        }

        // update the number of stage k dataset to extract
        ++niter[k];
      }


      // number of accrued patients and total number of events
      accruals1 = sum(n1);
      accruals2 = sum(n2);
      accruals3 = sum(n3);
      totalAccruals = accruals1 + accruals2 + accruals3;

      totalEvents = events1 + events2 + events3;
      totalDropouts = dropouts1 + dropouts2 + dropouts3;


      // order the data by time under observation
      NumericVector timeUnderObservationSorted =
        stl_sort(timeUnderObservation);
      IntegerVector sortedIndex = match(timeUnderObservationSorted,
                                        timeUnderObservation);
      sortedIndex = sortedIndex - 1;
      IntegerVector stratumSorted = stratum[sortedIndex];
      IntegerVector treatmentGroupSorted = treatmentGroup[sortedIndex];
      LogicalVector eventSorted = event[sortedIndex];

      LogicalVector subSorted = (timeUnderObservationSorted > 0);
      stratumSorted = stratumSorted[subSorted];
      treatmentGroupSorted = treatmentGroupSorted[subSorted];
      eventSorted = eventSorted[subSorted];
      int nsubSorted = static_cast<int>(eventSorted.size());

      // calculate the stratified log-rank test
      uscore13 = 0;
      vscore13 = 0;
      uscore23 = 0;
      vscore23 = 0;
      uscore12 = 0;
      vscore12 = 0;
      km13.fill(1);
      km23.fill(1);
      km12.fill(1);
      for (int i=0; i<nsubSorted; ++i) {
        int h = stratumSorted[i] - 1;
        nt13[h] = n1[h] + n3[h];
        nt23[h] = n2[h] + n3[h];
        nt12[h] = n1[h] + n2[h];

        n13a[h] = n1[h]*hazardRatioH013;
        n23a[h] = n2[h]*hazardRatioH023;
        n12a[h] = n1[h]*hazardRatioH012;

        nt13a[h] = n13a[h] + n3[h];
        nt23a[h] = n23a[h] + n3[h];
        nt12a[h] = n12a[h] + n2[h];


        if (eventSorted[i] && (treatmentGroupSorted[i]==1 ||
            treatmentGroupSorted[i]==3)) {
          w13[h] = pow(km13[h], rho1)*pow(1-km13[h], rho2);
          uscore13 += w13[h]*((treatmentGroupSorted[i]==1)
                                - n13a[h]/nt13a[h]);
          vscore13 += w13[h]*w13[h]*n13a[h]*n3[h]/(nt13a[h]*nt13a[h]);
          km13[h] *= (1-1/nt13[h]); // update km estimate
        }

        if (eventSorted[i] && (treatmentGroupSorted[i]==2 ||
            treatmentGroupSorted[i]==3)) {
          w23[h] = pow(km23[h], rho1)*pow(1-km23[h], rho2);
          uscore23 += w23[h]*((treatmentGroupSorted[i]==2)
                                - n23a[h]/nt23a[h]);
          vscore23 += w23[h]*w23[h]*n23a[h]*n3[h]/(nt23a[h]*nt23a[h]);
          km23[h] *= (1-1/nt23[h]); // update km estimate
        }

        if (eventSorted[i] && (treatmentGroupSorted[i]==1 ||
            treatmentGroupSorted[i]==2)) {
          w12[h] = pow(km12[h], rho1)*pow(1-km12[h], rho2);
          uscore12 += w12[h]*((treatmentGroupSorted[i]==1)
                                - n12a[h]/nt12a[h]);
          vscore12 += w12[h]*w12[h]*n12a[h]*n2[h]/(nt12a[h]*nt12a[h]);
          km12[h] *= (1-1/nt12[h]); // update km estimate
        }

        // reduce the risk set
        if (treatmentGroupSorted[i]==1) {
          --n1[h];
        } else if (treatmentGroupSorted[i]==2) {
          --n2[h];
        } else {
          --n3[h];
        }
      }



      // add summary data to output
      iterationNumbery[index2] = iter+1;
      eventsNotAchievedy[index2] = eventsNotAchieved;
      stageNumbery[index2] = k+1;
      analysisTimey[index2] = time;
      accruals1y[index2] = accruals1;
      accruals2y[index2] = accruals2;
      accruals3y[index2] = accruals3;
      totalAccrualsy[index2] = totalAccruals;

      events1y[index2] = events1;
      events2y[index2] = events2;
      events3y[index2] = events3;
      totalEventsy[index2] = totalEvents;
      dropouts1y[index2] = dropouts1;
      dropouts2y[index2] = dropouts2;
      dropouts3y[index2] = dropouts3;
      totalDropoutsy[index2] = totalDropouts;

      logRankStatistic13y[index2] = uscore13/sqrt(vscore13);
      logRankStatistic23y[index2] = uscore23/sqrt(vscore23);
      logRankStatistic12y[index2] = uscore12/sqrt(vscore12);

      ++index2;


    } // end of stage

  } // end of iteration


  // only keep nonmissing records
  LogicalVector sub2 = !is_na(iterationNumbery);
  iterationNumbery = iterationNumbery[sub2];
  eventsNotAchievedy = eventsNotAchievedy[sub2];
  stageNumbery = stageNumbery[sub2];
  analysisTimey = analysisTimey[sub2];
  accruals1y = accruals1y[sub2];
  accruals2y = accruals2y[sub2];
  accruals3y = accruals3y[sub2];
  totalAccrualsy = totalAccrualsy[sub2];
  events1y = events1y[sub2];
  events2y = events2y[sub2];
  events3y = events3y[sub2];
  totalEventsy = totalEventsy[sub2];
  dropouts1y = dropouts1y[sub2];
  dropouts2y = dropouts2y[sub2];
  dropouts3y = dropouts3y[sub2];
  totalDropoutsy = totalDropoutsy[sub2];
  logRankStatistic13y = logRankStatistic13y[sub2];
  logRankStatistic23y = logRankStatistic23y[sub2];
  logRankStatistic12y = logRankStatistic12y[sub2];

  DataFrame sumdata = DataFrame::create(
    _["iterationNumber"] = iterationNumbery,
    _["eventsNotAchieved"] = eventsNotAchievedy,
    _["stageNumber"] = stageNumbery,
    _["analysisTime"] = analysisTimey,
    _["accruals1"] = accruals1y,
    _["accruals2"] = accruals2y,
    _["accruals3"] = accruals3y,
    _["totalAccruals"] = totalAccrualsy,
    _["events1"] = events1y,
    _["events2"] = events2y,
    _["events3"] = events3y,
    _["totalEvents"] = totalEventsy,
    _["dropouts1"] = dropouts1y,
    _["dropouts2"] = dropouts2y,
    _["dropouts3"] = dropouts3y,
    _["totalDropouts"] = totalDropoutsy,
    _["logRankStatistic13"] = logRankStatistic13y,
    _["logRankStatistic23"] = logRankStatistic23y,
    _["logRankStatistic12"] = logRankStatistic12y);


  List result;

  if (maxNumberOfRawDatasetsPerStage > 0) {
    LogicalVector sub1 = !is_na(iterationNumberx);
    iterationNumberx = iterationNumberx[sub1];
    stageNumberx = stageNumberx[sub1];
    analysisTimex = analysisTimex[sub1];
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    survivalTimex = survivalTimex[sub1];
    dropoutTimex = dropoutTimex[sub1];
    timeUnderObservationx = timeUnderObservationx[sub1];
    eventx = eventx[sub1];
    dropoutEventx = dropoutEventx[sub1];

    DataFrame rawdata = DataFrame::create(
      _["iterationNumber"] = iterationNumberx,
      _["stageNumber"] = stageNumberx,
      _["analysisTime"] = analysisTimex,
      _["subjectId"] = subjectIdx,
      _["arrivalTime"] = arrivalTimex,
      _["stratum"] = stratumx,
      _["treatmentGroup"] = treatmentGroupx,
      _["survivalTime"] = survivalTimex,
      _["dropoutTime"] = dropoutTimex,
      _["timeUnderObservation"] = timeUnderObservationx,
      _["event"] = eventx,
      _["dropoutEvent"] = dropoutEventx);

    result = List::create(_["sumdata"] = sumdata,
                          _["rawdata"] = rawdata);
  } else {
    result = List::create(_["sumdata"] = sumdata);
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
//'   The seed from the environment will be used if left unspecified.
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
//'     - \code{survivalTime1}: The underlying survival time for
//'       event endpoint 1 (PFS) for the subject.
//'
//'     - \code{dropoutTime1}: The underlying dropout time for
//'       event endpoint 1 (PFS) for the subject.
//'
//'     - \code{timeUnderObservation1}: The time under observation
//'       since randomization for event endpoint 1 (PFS) for the subject.
//'
//'     - \code{event1}: Whether the subject experienced event
//'       endpoint 1 (PFS).
//'
//'     - \code{dropoutEvent1}: Whether the subject dropped out for
//'       endpoint 1 (PFS).
//'
//'     - \code{survivalTime2}: The underlying survival time for
//'       event endpoint 2 (OS) for the subject.
//'
//'     - \code{dropoutTime2}: The underlying dropout time for
//'       event endpoint 2 (OS) for the subject.
//'
//'     - \code{timeUnderObservation2}: The time under observation
//'       since randomization for event endpoint 2 (OS) for the subject.
//'
//'     - \code{event2}: Whether the subject experienced event
//'       endpoint 2 (OS).
//'
//'     - \code{dropoutEvent2}: Whether the subject dropped out for
//'       endpoint 2 (OS).
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
List lrsim2e(const int kMax = 1,
             const int kMaxpfs = 1,
             const double hazardRatioH0pfs = 1,
             const double hazardRatioH0os = 1,
             const int allocation1 = 1,
             const int allocation2 = 1,
             const NumericVector& accrualTime = 0,
             const NumericVector& accrualIntensity = NA_REAL,
             const NumericVector& piecewiseSurvivalTime = 0,
             const NumericVector& stratumFraction = 1,
             const double rho_pd_os = 0,
             const NumericVector& lambda1pfs = NA_REAL,
             const NumericVector& lambda2pfs = NA_REAL,
             const NumericVector& lambda1os = NA_REAL,
             const NumericVector& lambda2os = NA_REAL,
             const NumericVector& gamma1pfs = 0,
             const NumericVector& gamma2pfs = 0,
             const NumericVector& gamma1os = 0,
             const NumericVector& gamma2os = 0,
             const int n = NA_INTEGER,
             const double followupTime = NA_REAL,
             const bool fixedFollowup = 0,
             const double rho1 = 0,
             const double rho2 = 0,
             const IntegerVector& plannedEvents = NA_INTEGER,
             const NumericVector& plannedTime = NA_REAL,
             const int maxNumberOfIterations = 1000,
             const int maxNumberOfRawDatasetsPerStage = 0,
             const int seed = NA_INTEGER) {

  // check input parameters
  int kMaxpfsx = kMaxpfs;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nints = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nints;
  int nints2 = nints == 1 ? 10 : nints + 10;
  int nsi2 = nstrata*nints2;

  NumericVector lambda1pfsx(nsi), lambda2pfsx(nsi);
  NumericVector lambda1osx(nsi), lambda2osx(nsi);
  NumericVector gamma1pfsx(nsi), gamma2pfsx(nsi);
  NumericVector gamma1osx(nsi), gamma2osx(nsi);
  NumericVector lambda1pd(nsi2), lambda2pd(nsi2);
  NumericVector gamma1pd(nsi), gamma2pd(nsi);

  NumericVector tau1pd(nsi2), tau2pd(nsi2);

  bool useEvents, eventsNotAchieved;


  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (kMaxpfs < 0) {
    kMaxpfsx = kMax;
  }

  if (kMaxpfsx > kMax) {
    stop("kMaxpfs must be less than or equal to kMax");
  }

  // whether to plan the analyses based on events or calendar time
  if (is_false(any(is_na(plannedEvents)))) {
    useEvents = 1;

    if (plannedEvents[0] <= 0) {
      stop("Elements of plannedEvents must be positive");
    }

    if (plannedEvents.size() != kMax) {
      stop("Invalid length for plannedEvents");
    }

    if (kMaxpfsx > 1) {
      IntegerVector plannedEvents1 = plannedEvents[Range(0,kMaxpfsx-1)];
      if (is_true(any(diff(plannedEvents1) <= 0))) {
        stop("plannedEvents for PFS must be increasing");
      }
    }

    if (kMax - kMaxpfsx > 1) {
      IntegerVector plannedEvents2 = plannedEvents[Range(kMaxpfsx, kMax-1)];
      if (is_true(any(diff(plannedEvents2) <= 0))) {
        stop("plannedEvents for OS must be increasing");
      }
    }
  } else if (is_false(any(is_na(plannedTime)))) {
    useEvents = 0;
    if (plannedTime[0] <= 0) {
      stop("Elements of plannedTime must be positive");
    }

    if (plannedTime.size() != kMax) {
      stop("Invalid length for plannedTime");
    }

    if (kMax > 1 && is_true(any(diff(plannedTime) <= 0))) {
      stop("Elements of plannedTime must be increasing");
    }
  } else {
    stop("Either plannedEvents or plannedTime must be given");
  }


  if (hazardRatioH0pfs <= 0) {
    stop("hazardRatioH0pfs must be positive");
  }

  if (hazardRatioH0os <= 0) {
    stop("hazardRatioH0os must be positive");
  }


  if (allocation1 < 1) {
    stop("allocation1 must be a positive integer");
  }

  if (allocation2 < 1) {
    stop("allocation2 must be a positive integer");
  }


  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }


  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nints > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }


  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }


  if (rho_pd_os <= -1 || rho_pd_os >= 1) {
    stop("rho_pd_os must lie in (-1, 1)");
  }


  if (is_true(any(is_na(lambda1pfs)))) {
    stop("lambda1pfs must be provided");
  }

  if (is_true(any(is_na(lambda2pfs)))) {
    stop("lambda2pfs must be provided");
  }

  if (is_true(any(is_na(lambda1os)))) {
    stop("lambda1os must be provided");
  }

  if (is_true(any(is_na(lambda2os)))) {
    stop("lambda2os must be provided");
  }

  if (is_true(any(lambda1pfs < 0))) {
    stop("lambda1pfs must be non-negative");
  }

  if (is_true(any(lambda2pfs < 0))) {
    stop("lambda2pfs must be non-negative");
  }

  if (is_true(any(lambda1os < 0))) {
    stop("lambda1os must be non-negative");
  }

  if (is_true(any(lambda2os < 0))) {
    stop("lambda2os must be non-negative");
  }


  if (is_true(any(gamma1pfs < 0))) {
    stop("gamma1pfs must be non-negative");
  }

  if (is_true(any(gamma2pfs < 0))) {
    stop("gamma2pfs must be non-negative");
  }

  if (is_true(any(gamma1os < 0))) {
    stop("gamma1os must be non-negative");
  }

  if (is_true(any(gamma2os < 0))) {
    stop("gamma2os must be non-negative");
  }


  if (lambda1pfs.size() == 1) {
    lambda1pfsx = rep(lambda1pfs, nsi);
  } else if (lambda1pfs.size() == nints) {
    lambda1pfsx = rep(lambda1pfs, nstrata);
  } else if (lambda1pfs.size() == nsi) {
    lambda1pfsx = lambda1pfs;
  } else {
    stop("Invalid length for lambda1pfs");
  }

  if (lambda2pfs.size() == 1) {
    lambda2pfsx = rep(lambda2pfs, nsi);
  } else if (lambda2pfs.size() == nints) {
    lambda2pfsx = rep(lambda2pfs, nstrata);
  } else if (lambda2pfs.size() == nsi) {
    lambda2pfsx = lambda2pfs;
  } else {
    stop("Invalid length for lambda2pfs");
  }



  if (lambda1os.size() == 1) {
    lambda1osx = rep(lambda1os, nsi);
  } else if (lambda1os.size() == nints) {
    lambda1osx = rep(lambda1os, nstrata);
  } else if (lambda1os.size() == nsi) {
    lambda1osx = lambda1os;
  } else {
    stop("Invalid length for lambda1os");
  }

  if (lambda2os.size() == 1) {
    lambda2osx = rep(lambda2os, nsi);
  } else if (lambda2os.size() == nints) {
    lambda2osx = rep(lambda2os, nstrata);
  } else if (lambda2os.size() == nsi) {
    lambda2osx = lambda2os;
  } else {
    stop("Invalid length for lambda2os");
  }



  if (gamma1pfs.size() == 1) {
    gamma1pfsx = rep(gamma1pfs, nsi);
  } else if (gamma1pfs.size() == nints) {
    gamma1pfsx = rep(gamma1pfs, nstrata);
  } else if (gamma1pfs.size() == nsi) {
    gamma1pfsx = gamma1pfs;
  } else {
    stop("Invalid length for gamma1pfs");
  }

  if (gamma2pfs.size() == 1) {
    gamma2pfsx = rep(gamma2pfs, nsi);
  } else if (gamma2pfs.size() == nints) {
    gamma2pfsx = rep(gamma2pfs, nstrata);
  } else if (gamma2pfs.size() == nsi) {
    gamma2pfsx = gamma2pfs;
  } else {
    stop("Invalid length for gamma2pfs");
  }


  if (gamma1os.size() == 1) {
    gamma1osx = rep(gamma1os, nsi);
  } else if (gamma1os.size() == nints) {
    gamma1osx = rep(gamma1os, nstrata);
  } else if (gamma1os.size() == nsi) {
    gamma1osx = gamma1os;
  } else {
    stop("Invalid length for gamma1os");
  }

  if (gamma2os.size() == 1) {
    gamma2osx = rep(gamma2os, nsi);
  } else if (gamma2os.size() == nints) {
    gamma2osx = rep(gamma2os, nstrata);
  } else if (gamma2os.size() == nsi) {
    gamma2osx = gamma2os;
  } else {
    stop("Invalid length for gamma2os");
  }

  if (is_true(any(lambda1pfsx <= lambda1osx))) {
    stop("lambda1pfs must be greater than lambda1os");
  }

  if (is_true(any(lambda2pfsx <= lambda2osx))) {
    stop("lambda2pfs must be greater than lambda2os");
  }

  if (is_true(any(gamma1pfsx < gamma1osx))) {
    stop("gamma1pfs must be greater than or equal to gamma1os");
  }

  if (is_true(any(gamma2pfsx < gamma2osx))) {
    stop("gamma2pfs must be greater than or equal to gamma2os");
  }

  for (int j=0; j<nstrata; ++j) {
    Range jj = Range(j*nints, (j+1)*nints-1);

    NumericVector lam1pfs = lambda1pfsx[jj];
    NumericVector lam1os = lambda1osx[jj];

    List a1 = hazard_pdcpp(piecewiseSurvivalTime,
                           lam1pfs, lam1os, rho_pd_os);

    NumericVector time1pd = a1["piecewiseSurvivalTime"];
    NumericVector lam1pd = a1["hazard_pd"];

    NumericVector lam2pfs = lambda2pfsx[jj];
    NumericVector lam2os = lambda2osx[jj];

    List a2 = hazard_pdcpp(piecewiseSurvivalTime,
                           lam2pfs, lam2os, rho_pd_os);

    NumericVector time2pd = a2["piecewiseSurvivalTime"];
    NumericVector lam2pd = a2["hazard_pd"];

    for (int k=0; k<nints2; ++k) {
      tau1pd[j*nints2 + k] = time1pd[k];
      tau2pd[j*nints2 + k] = time2pd[k];
      lambda1pd[j*nints + k] = lam1pd[k];
      lambda2pd[j*nints + k] = lam2pd[k];
    }
  }

  gamma1pd = gamma1pfsx - gamma1osx;
  gamma2pd = gamma2pfsx - gamma2osx;


  if (n == NA_INTEGER) {
    stop("n must be provided");
  }

  if (n <= 0) {
    stop("n must be positive");
  }

  if (fixedFollowup) {
    if (R_isnancpp(followupTime)) {
      stop("followupTime must be provided for fixed follow-up");
    }

    if (followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }
  }

  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }

  if (rho2 < 0) {
    stop("rho2 must be non-negative");
  }

  if (maxNumberOfIterations < 1) {
    stop("maxNumberOfIterations must be a positive integer");
  }

  if (maxNumberOfRawDatasetsPerStage < 0) {
    stop("maxNumberOfRawDatasetsPerStage must be a non-negative integer");
  }


  // declare variables
  int nevents1, nevents2, nstages;
  int accruals1, accruals2, totalAccruals;
  int events1e1, events2e1, totalEventse1;
  int events1e2, events2e2, totalEventse2;
  int dropouts1e1, dropouts2e1, totalDropoutse1;
  int dropouts1e2, dropouts2e2, totalDropoutse2;
  int index1=0, index2=0;

  double enrollt, u, u1, u2, time, uscore, vscore;


  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n);

  NumericVector arrivalTime(n), survivalTime1(n), survivalTime2(n),
  dropoutTime1(n), dropoutTime2(n), timeUnderObservation1(n),
  timeUnderObservation2(n), totalTime1(n), totalTime2(n),
  totalt1(n), totalt2(n);

  LogicalVector event1(n), event2(n), dropoutEvent1(n), dropoutEvent2(n),
  event1ac(n), event2ac(n);


  // stratum information
  IntegerVector b1(nstrata), b2(nstrata);
  IntegerVector n1(nstrata), n2(nstrata);

  // original copy of n1 and n2 when looping over the endpoints
  IntegerVector n1x(nstrata), n2x(nstrata);

  // hazardRatioH0 adjusted n1 and nt for calculating the log-rank statistic
  NumericVector nt(nstrata), n1a(nstrata), nta(nstrata);

  NumericVector km(nstrata), w(nstrata);
  NumericVector cumStratumFraction = cumsum(stratumFraction);

  // within-stratum hazard rates
  NumericVector time1pd(nints2), time2pd(nints2);
  NumericVector lam1pd(nints2), lam2pd(nints2);
  NumericVector lam1os(nints), lam2os(nints);
  NumericVector gam1pd(nints), gam2pd(nints);
  NumericVector gam1os(nints), gam2os(nints);


  // stage-wise information
  IntegerVector niter(kMax);
  NumericVector analysisTime(kMax);


  // cache for the patient-level raw data to extract
  int nrow1 = n*kMax*maxNumberOfRawDatasetsPerStage;

  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stageNumberx(nrow1);
  NumericVector analysisTimex(nrow1);
  IntegerVector subjectIdx(nrow1);
  NumericVector arrivalTimex(nrow1);
  IntegerVector stratumx(nrow1);
  IntegerVector treatmentGroupx(nrow1);
  NumericVector survivalTime1x(nrow1);
  NumericVector survivalTime2x(nrow1);
  NumericVector dropoutTime1x(nrow1);
  NumericVector dropoutTime2x(nrow1);
  NumericVector timeUnderObservation1x(nrow1);
  NumericVector timeUnderObservation2x(nrow1);
  LogicalVector event1x(nrow1);
  LogicalVector event2x(nrow1);
  LogicalVector dropoutEvent1x(nrow1);
  LogicalVector dropoutEvent2x(nrow1);


  // cache for the simulation-level summary data to extract
  int nrow2 = kMax*maxNumberOfIterations*2;

  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  LogicalVector eventsNotAchievedy(nrow2);
  IntegerVector stageNumbery(nrow2);
  NumericVector analysisTimey(nrow2);
  IntegerVector accruals1y(nrow2);
  IntegerVector accruals2y(nrow2);
  IntegerVector totalAccrualsy(nrow2);
  IntegerVector endpointy(nrow2);
  IntegerVector events1y(nrow2);
  IntegerVector events2y(nrow2);
  IntegerVector totalEventsy(nrow2);
  IntegerVector dropouts1y(nrow2);
  IntegerVector dropouts2y(nrow2);
  IntegerVector totalDropoutsy(nrow2);
  NumericVector uscorey(nrow2);
  NumericVector vscorey(nrow2);
  NumericVector logRankStatisticy(nrow2);


  // set up random seed
  if (seed != NA_INTEGER) {
    set_seed(seed);
  }


  // simulation
  for (int iter=0; iter<maxNumberOfIterations; ++iter) {

    b1.fill(allocation1);
    b2.fill(allocation2);

    enrollt = 0;
    for (int i=0; i<n; ++i) {

      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);
      arrivalTime[i] = enrollt;

      // generate stratum information
      u = R::runif(0,1);
      int j;
      for (j=0; j<nstrata; ++j) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }

      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+0.0)) {
        treatmentGroup[i] = 1;
        --b1[j];
      } else {
        treatmentGroup[i] = 2;
        --b2[j];
      }

      // start a new block after depleting the current block
      if (b1[j]+b2[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
      }

      // stratum-specific hazard rates for event and dropout
      Range jj = Range(j*nints, (j+1)*nints-1);
      Range jj2 = Range(j*nints2, (j+1)*nints2-1);

      time1pd = tau1pd[jj2];
      time2pd = tau2pd[jj2];

      lam1pd = lambda1pd[jj2];
      lam2pd = lambda2pd[jj2];

      lam1os = lambda1osx[jj];
      lam2os = lambda2osx[jj];

      gam1pd = gamma1pd[jj];
      gam2pd = gamma2pd[jj];

      gam1os = gamma1osx[jj];
      gam2os = gamma2osx[jj];

      // standard bivariate normal with correlation rho_pd_os
      u1 = R::rnorm(0,1);
      u2 = R::rnorm(rho_pd_os*u1, sqrt(1-rho_pd_os*rho_pd_os));

      // transform to uniform
      u1 = R::pnorm(u1, 0, 1, 1, 0);
      u2 = R::pnorm(u2, 0, 1, 1, 0);

      // generate survival times
      if (treatmentGroup[i]==1) {
        survivalTime1[i] = qtpwexpcpp1(u1, time1pd, lam1pd, 0, 1, 0);
        survivalTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,lam1os,0,1,0);
      } else {
        survivalTime1[i] = qtpwexpcpp1(u1, time2pd, lam2pd, 0, 1, 0);
        survivalTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,lam2os,0,1,0);
      }
      // PFS includes death
      survivalTime1[i] = std::min(survivalTime1[i], survivalTime2[i]);


      // generate dropout times
      u1 = R::runif(0,1);
      u2 = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime1[i] = qtpwexpcpp1(u1,piecewiseSurvivalTime,gam1pd,0,1,0);
        dropoutTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,gam1os,0,1,0);
      } else {
        dropoutTime1[i] = qtpwexpcpp1(u1,piecewiseSurvivalTime,gam2pd,0,1,0);
        dropoutTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,gam2os,0,1,0);
      }
      // whatever censors OS will also censor PFS
      dropoutTime1[i] = std::min(dropoutTime1[i], dropoutTime2[i]);


      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime1[i] <= dropoutTime1[i] &&
            survivalTime1[i] <= followupTime) {
          timeUnderObservation1[i] = survivalTime1[i];
          event1[i] = 1;
          dropoutEvent1[i] = 0;
        } else if (dropoutTime1[i] <= survivalTime1[i] &&
          dropoutTime1[i] <= followupTime) {
          timeUnderObservation1[i] = dropoutTime1[i];
          event1[i] = 0;
          dropoutEvent1[i] = 1;
        } else {
          timeUnderObservation1[i] = followupTime;
          event1[i] = 0;
          dropoutEvent1[i] = 0;
        }
      } else {
        if (survivalTime1[i] <= dropoutTime1[i]) {
          timeUnderObservation1[i] = survivalTime1[i];
          event1[i] = 1;
          dropoutEvent1[i] = 0;
        } else {
          timeUnderObservation1[i] = dropoutTime1[i];
          event1[i] = 0;
          dropoutEvent1[i] = 1;
        }
      }

      totalTime1[i] = arrivalTime[i] + timeUnderObservation1[i];


      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime2[i] <= dropoutTime2[i] &&
            survivalTime2[i] <= followupTime) {
          timeUnderObservation2[i] = survivalTime2[i];
          event2[i] = 1;
          dropoutEvent2[i] = 0;
        } else if (dropoutTime2[i] <= survivalTime2[i] &&
          dropoutTime2[i] <= followupTime) {
          timeUnderObservation2[i] = dropoutTime2[i];
          event2[i] = 0;
          dropoutEvent2[i] = 1;
        } else {
          timeUnderObservation2[i] = followupTime;
          event2[i] = 0;
          dropoutEvent2[i] = 0;
        }
      } else {
        if (survivalTime2[i] <= dropoutTime2[i]) {
          timeUnderObservation2[i] = survivalTime2[i];
          event2[i] = 1;
          dropoutEvent2[i] = 0;
        } else {
          timeUnderObservation2[i] = dropoutTime2[i];
          event2[i] = 0;
          dropoutEvent2[i] = 1;
        }
      }

      totalTime2[i] = arrivalTime[i] + timeUnderObservation2[i];

    }


    // find the analysis time for each stage
    if (useEvents) {
      nevents1 = sum(event1);
      nevents2 = sum(event2);
      totalt1 = stl_sort(totalTime1[event1]);
      totalt2 = stl_sort(totalTime2[event2]);

      int j1 = kMaxpfsx, j2 = kMax - kMaxpfsx;

      // PFS looks
      if (kMaxpfsx > 0) {
        for (j1=0; j1<kMaxpfsx; ++j1) {
          if (plannedEvents[j1] >= nevents1) {
            break;
          }
        }

        if (j1==kMaxpfsx) { // total number of PFS events exceeds planned
          for (int k=0; k<kMaxpfsx; ++k) {
            analysisTime[k] = totalt1[plannedEvents[k]-1] + 1e-12;
          }
        } else {
          for (int k=0; k<=j1; ++k) {
            if (k < j1) {
              analysisTime[k] = totalt1[plannedEvents[k]-1] + 1e-12;
            } else {
              analysisTime[k] = totalt1[nevents1-1] + 1e-12;
            }
          }
        }
      }

      // OS looks
      NumericVector analysisTime2(kMax - kMaxpfsx);

      if (kMax > kMaxpfsx) {
        for (j2=0; j2<kMax-kMaxpfsx; ++j2) {
          if (plannedEvents[kMaxpfsx+j2] >= nevents2) {
            break;
          }
        }

        if (j2==kMax-kMaxpfsx) { // total number of OS events exceeds planned
          for (int k=0; k<kMax-kMaxpfsx; ++k) {
            analysisTime2[k] = totalt2[plannedEvents[kMaxpfsx+k]-1] + 1e-12;
          }
        } else {
          for (int k=0; k<=j2; ++k) {
            if (k < j2) {
              analysisTime2[k] = totalt2[plannedEvents[kMaxpfsx+k]-1] + 1e-12;
            } else {
              analysisTime2[k] = totalt2[nevents2-1] + 1e-12;
            }
          }
        }
      }

      // determine the number of looks and timing of the looks
      if (kMaxpfsx == 0) { // all looks based on OS events
        if (j2 == kMax - kMaxpfsx) {
          nstages = kMax - kMaxpfsx;
        } else {
          nstages = j2 + 1;
        }

        for (int k=0; k<nstages; ++k) {
          analysisTime[k] = analysisTime2[k];
        }
      } else if (kMax == kMaxpfsx) { // all looks based on PFS events
        if (j1 == kMaxpfsx) {
          nstages = kMaxpfsx;
        } else {
          nstages = j1 + 1;
        }
      } else {
        if (analysisTime2[kMax-kMaxpfsx-1] > analysisTime[kMaxpfsx-1]) {
          // only OS looks that occur after the last PFS look contribute
          int l = static_cast<int>(which_max(analysisTime2 >
                                               analysisTime[kMaxpfsx-1]));
          nstages = kMax-l;
          for (int k=kMaxpfsx; k<kMax-l; ++k) {
            analysisTime[k] = analysisTime2[k-kMaxpfsx+l];
          }
        } else {
          if (j1 == kMaxpfsx) {
            nstages = kMaxpfsx;
          } else {
            nstages = j1 + 1;
          }
        }
      }

      // whether the target PFS and OS events are achieved
      if (kMaxpfsx > 0 && nevents1 < plannedEvents[kMaxpfsx-1]) {
        eventsNotAchieved = 1;
      } else if (kMaxpfsx < kMax && nevents2 < plannedEvents[kMax-1]) {
        eventsNotAchieved = 1;
      } else {
        eventsNotAchieved = 0;
      }
    } else { // looks based on calendar time
      nstages = kMax;
      analysisTime = clone(plannedTime);
      eventsNotAchieved = 0;
    }


    // construct the log-rank test statistic at each stage
    for (int k=0; k<nstages; ++k) {
      time = analysisTime[k];

      n1x.fill(0);  // number of subjects in each stratum by treatment
      n2x.fill(0);

      events1e1 = 0;
      events2e1 = 0;

      dropouts1e1 = 0;
      dropouts2e1 = 0;

      events1e2 = 0;
      events2e2 = 0;

      dropouts1e2 = 0;
      dropouts2e2 = 0;

      // censor at analysis time
      for (int i=0; i<n; ++i) {
        int h = stratum[i]-1;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation1[i] = time - arrivalTime[i];
          event1[i] = 0;
          dropoutEvent1[i] = 0;

          timeUnderObservation2[i] = time - arrivalTime[i];
          event2[i] = 0;
          dropoutEvent2[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            ++n1x[h];
          } else {
            ++n2x[h];
          }


          // censored time for endpoint 1
          if (fixedFollowup) {
            if (arrivalTime[i] + survivalTime1[i] <= time &&
                survivalTime1[i] <= dropoutTime1[i] &&
                survivalTime1[i] <= followupTime) {
              timeUnderObservation1[i] = survivalTime1[i];
              event1[i] = 1;
              dropoutEvent1[i] = 0;
            } else if (arrivalTime[i] + dropoutTime1[i] <= time &&
              dropoutTime1[i] <= survivalTime1[i] &&
              dropoutTime1[i] <= followupTime) {
              timeUnderObservation1[i] = dropoutTime1[i];
              event1[i] = 0;
              dropoutEvent1[i] = 1;
            } else if (arrivalTime[i] + followupTime <= time &&
              followupTime <= survivalTime1[i] &&
              followupTime <= dropoutTime1[i]) {
              timeUnderObservation1[i] = followupTime;
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            } else {
              timeUnderObservation1[i] = time - arrivalTime[i];
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime1[i] <= time &&
                survivalTime1[i] <= dropoutTime1[i]) {
              timeUnderObservation1[i] = survivalTime1[i];
              event1[i] = 1;
              dropoutEvent1[i] = 0;
            } else if (arrivalTime[i] + dropoutTime1[i] <= time &&
              dropoutTime1[i] <= survivalTime1[i]) {
              timeUnderObservation1[i] = dropoutTime1[i];
              event1[i] = 0;
              dropoutEvent1[i] = 1;
            } else {
              timeUnderObservation1[i] = time - arrivalTime[i];
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            }
          }

          if (treatmentGroup[i]==1 && event1[i]) ++events1e1;
          if (treatmentGroup[i]==2 && event1[i]) ++events2e1;
          if (treatmentGroup[i]==1 && dropoutEvent1[i]) ++dropouts1e1;
          if (treatmentGroup[i]==2 && dropoutEvent1[i]) ++dropouts2e1;


          // censored time for endpoint 2
          if (fixedFollowup) {
            if (arrivalTime[i] + survivalTime2[i] <= time &&
                survivalTime2[i] <= dropoutTime2[i] &&
                survivalTime2[i] <= followupTime) {
              timeUnderObservation2[i] = survivalTime2[i];
              event2[i] = 1;
              dropoutEvent2[i] = 0;
            } else if (arrivalTime[i] + dropoutTime2[i] <= time &&
              dropoutTime2[i] <= survivalTime2[i] &&
              dropoutTime2[i] <= followupTime) {
              timeUnderObservation2[i] = dropoutTime2[i];
              event2[i] = 0;
              dropoutEvent2[i] = 1;
            } else if (arrivalTime[i] + followupTime <= time &&
              followupTime <= survivalTime2[i] &&
              followupTime <= dropoutTime2[i]) {
              timeUnderObservation2[i] = followupTime;
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            } else {
              timeUnderObservation2[i] = time - arrivalTime[i];
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime2[i] <= time &&
                survivalTime2[i] <= dropoutTime2[i]) {
              timeUnderObservation2[i] = survivalTime2[i];
              event2[i] = 1;
              dropoutEvent2[i] = 0;
            } else if (arrivalTime[i] + dropoutTime2[i] <= time &&
              dropoutTime2[i] <= survivalTime2[i]) {
              timeUnderObservation2[i] = dropoutTime2[i];
              event2[i] = 0;
              dropoutEvent2[i] = 1;
            } else {
              timeUnderObservation2[i] = time - arrivalTime[i];
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            }
          }

          if (treatmentGroup[i]==1 && event2[i]) ++events1e2;
          if (treatmentGroup[i]==2 && event2[i]) ++events2e2;
          if (treatmentGroup[i]==1 && dropoutEvent2[i]) ++dropouts1e2;
          if (treatmentGroup[i]==2 && dropoutEvent2[i]) ++dropouts2e2;
        }
      }


      // add raw data to output
      if (niter[k] < maxNumberOfRawDatasetsPerStage) {
        for (int i=0; i<n; ++i) {
          iterationNumberx[index1] = iter+1;
          stageNumberx[index1] = k+1;
          analysisTimex[index1] = time;
          subjectIdx[index1] = i+1;
          arrivalTimex[index1] = arrivalTime[i];
          stratumx[index1] = stratum[i];
          treatmentGroupx[index1] = treatmentGroup[i];

          survivalTime1x[index1] = survivalTime1[i];
          dropoutTime1x[index1] = dropoutTime1[i];
          timeUnderObservation1x[index1] = timeUnderObservation1[i];
          event1x[index1] = event1[i];
          dropoutEvent1x[index1] = dropoutEvent1[i];

          survivalTime2x[index1] = survivalTime2[i];
          dropoutTime2x[index1] = dropoutTime2[i];
          timeUnderObservation2x[index1] = timeUnderObservation2[i];
          event2x[index1] = event2[i];
          dropoutEvent2x[index1] = dropoutEvent2[i];

          ++index1;
        }

        // update the number of stage k dataset to extract
        ++niter[k];
      }


      // number of accrued patients and total number of events
      accruals1 = sum(n1x);
      accruals2 = sum(n2x);
      totalAccruals = accruals1 + accruals2;

      totalEventse1 = events1e1 + events2e1;
      totalDropoutse1 = dropouts1e1 + dropouts2e1;

      totalEventse2 = events1e2 + events2e2;
      totalDropoutse2 = dropouts1e2 + dropouts2e2;


      for (int endpoint=1; endpoint<=2; ++endpoint) {
        n1 = clone(n1x);
        n2 = clone(n2x);

        double hazardRatioH0;
        if (endpoint == 1) {
          hazardRatioH0 = hazardRatioH0pfs;
        } else {
          hazardRatioH0 = hazardRatioH0os;
        }

        // order the data by time under observation
        NumericVector timeUnderObservationSorted;
        IntegerVector sortedIndex;
        LogicalVector eventSorted;
        if (endpoint == 1) {
          timeUnderObservationSorted = stl_sort(timeUnderObservation1);
          sortedIndex = match(timeUnderObservationSorted,
                              timeUnderObservation1);
          sortedIndex = sortedIndex - 1;
          eventSorted = event1[sortedIndex];
        } else {
          timeUnderObservationSorted = stl_sort(timeUnderObservation2);
          sortedIndex = match(timeUnderObservationSorted,
                              timeUnderObservation2);
          sortedIndex = sortedIndex - 1;
          eventSorted = event2[sortedIndex];
        }

        IntegerVector stratumSorted = stratum[sortedIndex];
        IntegerVector treatmentGroupSorted = treatmentGroup[sortedIndex];

        LogicalVector subSorted = (timeUnderObservationSorted > 0);
        stratumSorted = stratumSorted[subSorted];
        treatmentGroupSorted = treatmentGroupSorted[subSorted];
        eventSorted = eventSorted[subSorted];
        int nsubSorted = static_cast<int>(eventSorted.size());

        // calculate the stratified log-rank test
        uscore = 0;
        vscore = 0;
        km.fill(1);
        for (int i=0; i<nsubSorted; ++i) {
          int h = stratumSorted[i] - 1;
          nt[h] = n1[h] + n2[h];

          n1a[h] = n1[h]*hazardRatioH0;
          nta[h] = n1a[h] + n2[h];

          if (eventSorted[i]) {
            w[h] = pow(km[h], rho1)*pow(1-km[h], rho2);
            uscore += w[h]*((treatmentGroupSorted[i]==1) - n1a[h]/nta[h]);
            vscore += w[h]*w[h]*n1a[h]*n2[h]/(nta[h]*nta[h]);
            km[h] *= (1-1/nt[h]); // update km estimate
          }

          // reduce the risk set
          if (treatmentGroupSorted[i]==1) {
            --n1[h];
          } else {
            --n2[h];
          }
        }



        // add summary data to output
        iterationNumbery[index2] = iter+1;
        eventsNotAchievedy[index2] = eventsNotAchieved;
        stageNumbery[index2] = k+1;
        analysisTimey[index2] = time;
        accruals1y[index2] = accruals1;
        accruals2y[index2] = accruals2;
        totalAccrualsy[index2] = totalAccruals;
        endpointy[index2] = endpoint;

        if (endpoint == 1) {
          events1y[index2] = events1e1;
          events2y[index2] = events2e1;
          totalEventsy[index2] = totalEventse1;
          dropouts1y[index2] = dropouts1e1;
          dropouts2y[index2] = dropouts2e1;
          totalDropoutsy[index2] = totalDropoutse1;
        } else {
          events1y[index2] = events1e2;
          events2y[index2] = events2e2;
          totalEventsy[index2] = totalEventse2;
          dropouts1y[index2] = dropouts1e2;
          dropouts2y[index2] = dropouts2e2;
          totalDropoutsy[index2] = totalDropoutse2;
        }

        uscorey[index2] = uscore;
        vscorey[index2] = vscore;
        logRankStatisticy[index2] = uscore/sqrt(vscore);
        ++index2;

      } // end of endpoint

    } // end of stage

  } // end of iteration


  // simulation summary data set
  LogicalVector sub2 = !is_na(iterationNumbery);
  iterationNumbery = iterationNumbery[sub2];
  eventsNotAchievedy = eventsNotAchievedy[sub2];
  stageNumbery = stageNumbery[sub2];
  analysisTimey = analysisTimey[sub2];
  accruals1y = accruals1y[sub2];
  accruals2y = accruals2y[sub2];
  totalAccrualsy = totalAccrualsy[sub2];
  endpointy = endpointy[sub2];
  events1y = events1y[sub2];
  events2y = events2y[sub2];
  totalEventsy = totalEventsy[sub2];
  dropouts1y = dropouts1y[sub2];
  dropouts2y = dropouts2y[sub2];
  totalDropoutsy = totalDropoutsy[sub2];
  uscorey = uscorey[sub2];
  vscorey = vscorey[sub2];
  logRankStatisticy = logRankStatisticy[sub2];

  DataFrame sumdata = DataFrame::create(
    _["iterationNumber"] = iterationNumbery,
    _["eventsNotAchieved"] = eventsNotAchievedy,
    _["stageNumber"] = stageNumbery,
    _["analysisTime"] = analysisTimey,
    _["accruals1"] = accruals1y,
    _["accruals2"] = accruals2y,
    _["totalAccruals"] = totalAccrualsy,
    _["endpoint"] = endpointy,
    _["events1"] = events1y,
    _["events2"] = events2y,
    _["totalEvents"] = totalEventsy,
    _["dropouts1"] = dropouts1y,
    _["dropouts2"] = dropouts2y,
    _["totalDropouts"] = totalDropoutsy,
    _["uscore"] = uscorey,
    _["vscore"] = vscorey,
    _["logRankStatistic"] = logRankStatisticy);


  List result;

  if (maxNumberOfRawDatasetsPerStage > 0) {
    LogicalVector sub1 = !is_na(iterationNumberx);
    iterationNumberx = iterationNumberx[sub1];
    stageNumberx = stageNumberx[sub1];
    analysisTimex = analysisTimex[sub1];
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    survivalTime1x = survivalTime1x[sub1];
    dropoutTime1x = dropoutTime1x[sub1];
    timeUnderObservation1x = timeUnderObservation1x[sub1];
    event1x = event1x[sub1];
    dropoutEvent1x = dropoutEvent1x[sub1];
    survivalTime2x = survivalTime2x[sub1];
    dropoutTime2x = dropoutTime2x[sub1];
    timeUnderObservation2x = timeUnderObservation2x[sub1];
    event2x = event2x[sub1];
    dropoutEvent2x = dropoutEvent2x[sub1];

    DataFrame rawdata = DataFrame::create(
      _["iterationNumber"] = iterationNumberx,
      _["stageNumber"] = stageNumberx,
      _["analysisTime"] = analysisTimex,
      _["subjectId"] = subjectIdx,
      _["arrivalTime"] = arrivalTimex,
      _["stratum"] = stratumx,
      _["treatmentGroup"] = treatmentGroupx,
      _["survivalTime1"] = survivalTime1x,
      _["dropoutTime1"] = dropoutTime1x,
      _["timeUnderObservation1"] = timeUnderObservation1x,
      _["event1"] = event1x,
      _["dropoutEvent1"] = dropoutEvent1x,
      _["survivalTime2"] = survivalTime2x,
      _["dropoutTime2"] = dropoutTime2x,
      _["timeUnderObservation2"] = timeUnderObservation2x,
      _["event2"] = event2x,
      _["dropoutEvent2"] = dropoutEvent2x);

    result = List::create(_["sumdata"] = sumdata,
                          _["rawdata"] = rawdata);
  } else {
    result = List::create(_["sumdata"] = sumdata);
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
//'   The seed from the environment will be used if left unspecified.
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
//'     - \code{survivalTime1}: The underlying survival time for
//'       event endpoint 1 (PFS) for the subject.
//'
//'     - \code{dropoutTime1}: The underlying dropout time for
//'       event endpoint 1 (PFS) for the subject.
//'
//'     - \code{timeUnderObservation1}: The time under observation
//'       since randomization for event endpoint 1 (PFS) for the subject.
//'
//'     - \code{event1}: Whether the subject experienced event
//'       endpoint 1 (PFS) .
//'
//'     - \code{dropoutEvent1}: Whether the subject dropped out for
//'       endpoint 1 (PFS).
//'
//'     - \code{survivalTime2}: The underlying survival time for
//'       event endpoint 2 (OS) for the subject.
//'
//'     - \code{dropoutTime2}: The underlying dropout time for
//'       event endpoint 2 (OS) for the subject.
//'
//'     - \code{timeUnderObservation2}: The time under observation
//'       since randomization for event endpoint 2 (OS) for the subject.
//'
//'     - \code{event2}: Whether the subject experienced event
//'       endpoint 2 (OS).
//'
//'     - \code{dropoutEvent2}: Whether the subject dropped out for
//'       endpoint 2 (OS).
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
//'   maxNumberOfIterations = 500,
//'   maxNumberOfRawDatasetsPerStage = 1,
//'   seed = 314159)
//'
//' head(sim1$sumdata)
//' head(sim1$rawdata)
//'
//' @export
// [[Rcpp::export]]
List lrsim2e3a(const int kMax = 1,
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
               const NumericVector& accrualTime = 0,
               const NumericVector& accrualIntensity = NA_REAL,
               const NumericVector& piecewiseSurvivalTime = 0,
               const NumericVector& stratumFraction = 1,
               const double rho_pd_os = 0,
               const NumericVector& lambda1pfs = NA_REAL,
               const NumericVector& lambda2pfs = NA_REAL,
               const NumericVector& lambda3pfs = NA_REAL,
               const NumericVector& lambda1os = NA_REAL,
               const NumericVector& lambda2os = NA_REAL,
               const NumericVector& lambda3os = NA_REAL,
               const NumericVector& gamma1pfs = 0,
               const NumericVector& gamma2pfs = 0,
               const NumericVector& gamma3pfs = 0,
               const NumericVector& gamma1os = 0,
               const NumericVector& gamma2os = 0,
               const NumericVector& gamma3os = 0,
               const int n = NA_INTEGER,
               const double followupTime = NA_REAL,
               const bool fixedFollowup = 0,
               const double rho1 = 0,
               const double rho2 = 0,
               const IntegerVector& plannedEvents = NA_INTEGER,
               const NumericVector& plannedTime = NA_REAL,
               const int maxNumberOfIterations = 1000,
               const int maxNumberOfRawDatasetsPerStage = 0,
               const int seed = NA_INTEGER) {

  // check input parameters
  int kMaxpfsx = kMaxpfs;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nints = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nints;
  int nints2 = nints == 1 ? 10 : nints + 10;
  int nsi2 = nstrata*nints2;

  NumericVector lambda1pfsx(nsi), lambda2pfsx(nsi), lambda3pfsx(nsi);
  NumericVector lambda1osx(nsi), lambda2osx(nsi), lambda3osx(nsi);
  NumericVector gamma1pfsx(nsi), gamma2pfsx(nsi), gamma3pfsx(nsi);
  NumericVector gamma1osx(nsi), gamma2osx(nsi), gamma3osx(nsi);
  NumericVector lambda1pd(nsi2), lambda2pd(nsi2), lambda3pd(nsi2);
  NumericVector gamma1pd(nsi), gamma2pd(nsi), gamma3pd(nsi);

  NumericVector tau1pd(nsi2), tau2pd(nsi2), tau3pd(nsi2);

  bool useEvents, eventsNotAchieved;


  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (kMaxpfs < 0) {
    kMaxpfsx = kMax;
  }

  if (kMaxpfsx > kMax) {
    stop("kMaxpfs must be less than or equal to kMax");
  }


  // whether to plan the analyses based on events or calendar time
  if (is_false(any(is_na(plannedEvents)))) {
    useEvents = 1;

    if (plannedEvents[0] <= 0) {
      stop("Elements of plannedEvents must be positive");
    }

    if (plannedEvents.size() != kMax) {
      stop("Invalid length for plannedEvents");
    }

    if (kMaxpfsx > 1) {
      IntegerVector plannedEvents1 = plannedEvents[Range(0,kMaxpfsx-1)];
      if (is_true(any(diff(plannedEvents1) <= 0))) {
        stop("plannedEvents for PFS must be increasing");
      }
    }

    if (kMax - kMaxpfsx > 1) {
      IntegerVector plannedEvents2 = plannedEvents[Range(kMaxpfsx, kMax-1)];
      if (is_true(any(diff(plannedEvents2) <= 0))) {
        stop("plannedEvents for OS must be increasing");
      }
    }
  } else if (is_false(any(is_na(plannedTime)))) {
    useEvents = 0;
    if (plannedTime[0] <= 0) {
      stop("Elements of plannedTime must be positive");
    }

    if (plannedTime.size() != kMax) {
      stop("Invalid length for plannedTime");
    }

    if (kMax > 1 && is_true(any(diff(plannedTime) <= 0))) {
      stop("Elements of plannedTime must be increasing");
    }
  } else {
    stop("Either plannedEvents or plannedTime must be given");
  }

  if (hazardRatioH013pfs <= 0) {
    stop("hazardRatioH013pfs must be positive");
  }

  if (hazardRatioH023pfs <= 0) {
    stop("hazardRatioH023pfs must be positive");
  }

  if (hazardRatioH012pfs <= 0) {
    stop("hazardRatioH012pfs must be positive");
  }

  if (hazardRatioH013os <= 0) {
    stop("hazardRatioH013os must be positive");
  }

  if (hazardRatioH023os <= 0) {
    stop("hazardRatioH023os must be positive");
  }

  if (hazardRatioH012os <= 0) {
    stop("hazardRatioH012os must be positive");
  }

  if (allocation1 < 1) {
    stop("allocation1 must be a positive integer");
  }

  if (allocation2 < 1) {
    stop("allocation2 must be a positive integer");
  }

  if (allocation3 < 1) {
    stop("allocation3 must be a positive integer");
  }


  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }


  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nints > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }


  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }


  if (rho_pd_os <= -1 || rho_pd_os >= 1) {
    stop("rho_pd_os must lie in (-1, 1)");
  }


  if (is_true(any(is_na(lambda1pfs)))) {
    stop("lambda1pfs must be provided");
  }

  if (is_true(any(is_na(lambda2pfs)))) {
    stop("lambda2pfs must be provided");
  }

  if (is_true(any(is_na(lambda3pfs)))) {
    stop("lambda3pfs must be provided");
  }

  if (is_true(any(is_na(lambda1os)))) {
    stop("lambda1os must be provided");
  }

  if (is_true(any(is_na(lambda2os)))) {
    stop("lambda2os must be provided");
  }

  if (is_true(any(is_na(lambda3os)))) {
    stop("lambda3os must be provided");
  }

  if (is_true(any(lambda1pfs < 0))) {
    stop("lambda1pfs must be non-negative");
  }

  if (is_true(any(lambda2pfs < 0))) {
    stop("lambda2pfs must be non-negative");
  }

  if (is_true(any(lambda3pfs < 0))) {
    stop("lambda3pfs must be non-negative");
  }

  if (is_true(any(lambda1os < 0))) {
    stop("lambda1os must be non-negative");
  }

  if (is_true(any(lambda2os < 0))) {
    stop("lambda2os must be non-negative");
  }

  if (is_true(any(lambda3os < 0))) {
    stop("lambda3os must be non-negative");
  }


  if (is_true(any(gamma1pfs < 0))) {
    stop("gamma1pfs must be non-negative");
  }

  if (is_true(any(gamma2pfs < 0))) {
    stop("gamma2pfs must be non-negative");
  }

  if (is_true(any(gamma3pfs < 0))) {
    stop("gamma3pfs must be non-negative");
  }

  if (is_true(any(gamma1os < 0))) {
    stop("gamma1os must be non-negative");
  }

  if (is_true(any(gamma2os < 0))) {
    stop("gamma2os must be non-negative");
  }

  if (is_true(any(gamma3os < 0))) {
    stop("gamma3os must be non-negative");
  }


  if (lambda1pfs.size() == 1) {
    lambda1pfsx = rep(lambda1pfs, nsi);
  } else if (lambda1pfs.size() == nints) {
    lambda1pfsx = rep(lambda1pfs, nstrata);
  } else if (lambda1pfs.size() == nsi) {
    lambda1pfsx = lambda1pfs;
  } else {
    stop("Invalid length for lambda1pfs");
  }

  if (lambda2pfs.size() == 1) {
    lambda2pfsx = rep(lambda2pfs, nsi);
  } else if (lambda2pfs.size() == nints) {
    lambda2pfsx = rep(lambda2pfs, nstrata);
  } else if (lambda2pfs.size() == nsi) {
    lambda2pfsx = lambda2pfs;
  } else {
    stop("Invalid length for lambda2pfs");
  }

  if (lambda3pfs.size() == 1) {
    lambda3pfsx = rep(lambda3pfs, nsi);
  } else if (lambda3pfs.size() == nints) {
    lambda3pfsx = rep(lambda3pfs, nstrata);
  } else if (lambda3pfs.size() == nsi) {
    lambda3pfsx = lambda3pfs;
  } else {
    stop("Invalid length for lambda3pfs");
  }


  if (lambda1os.size() == 1) {
    lambda1osx = rep(lambda1os, nsi);
  } else if (lambda1os.size() == nints) {
    lambda1osx = rep(lambda1os, nstrata);
  } else if (lambda1os.size() == nsi) {
    lambda1osx = lambda1os;
  } else {
    stop("Invalid length for lambda1os");
  }

  if (lambda2os.size() == 1) {
    lambda2osx = rep(lambda2os, nsi);
  } else if (lambda2os.size() == nints) {
    lambda2osx = rep(lambda2os, nstrata);
  } else if (lambda2os.size() == nsi) {
    lambda2osx = lambda2os;
  } else {
    stop("Invalid length for lambda2os");
  }

  if (lambda3os.size() == 1) {
    lambda3osx = rep(lambda3os, nsi);
  } else if (lambda3os.size() == nints) {
    lambda3osx = rep(lambda3os, nstrata);
  } else if (lambda3os.size() == nsi) {
    lambda3osx = lambda3os;
  } else {
    stop("Invalid length for lambda3os");
  }


  if (gamma1pfs.size() == 1) {
    gamma1pfsx = rep(gamma1pfs, nsi);
  } else if (gamma1pfs.size() == nints) {
    gamma1pfsx = rep(gamma1pfs, nstrata);
  } else if (gamma1pfs.size() == nsi) {
    gamma1pfsx = gamma1pfs;
  } else {
    stop("Invalid length for gamma1pfs");
  }

  if (gamma2pfs.size() == 1) {
    gamma2pfsx = rep(gamma2pfs, nsi);
  } else if (gamma2pfs.size() == nints) {
    gamma2pfsx = rep(gamma2pfs, nstrata);
  } else if (gamma2pfs.size() == nsi) {
    gamma2pfsx = gamma2pfs;
  } else {
    stop("Invalid length for gamma2pfs");
  }

  if (gamma3pfs.size() == 1) {
    gamma3pfsx = rep(gamma3pfs, nsi);
  } else if (gamma3pfs.size() == nints) {
    gamma3pfsx = rep(gamma3pfs, nstrata);
  } else if (gamma3pfs.size() == nsi) {
    gamma3pfsx = gamma3pfs;
  } else {
    stop("Invalid length for gamma3pfs");
  }


  if (gamma1os.size() == 1) {
    gamma1osx = rep(gamma1os, nsi);
  } else if (gamma1os.size() == nints) {
    gamma1osx = rep(gamma1os, nstrata);
  } else if (gamma1os.size() == nsi) {
    gamma1osx = gamma1os;
  } else {
    stop("Invalid length for gamma1os");
  }

  if (gamma2os.size() == 1) {
    gamma2osx = rep(gamma2os, nsi);
  } else if (gamma2os.size() == nints) {
    gamma2osx = rep(gamma2os, nstrata);
  } else if (gamma2os.size() == nsi) {
    gamma2osx = gamma2os;
  } else {
    stop("Invalid length for gamma2os");
  }

  if (gamma3os.size() == 1) {
    gamma3osx = rep(gamma3os, nsi);
  } else if (gamma3os.size() == nints) {
    gamma3osx = rep(gamma3os, nstrata);
  } else if (gamma3os.size() == nsi) {
    gamma3osx = gamma3os;
  } else {
    stop("Invalid length for gamma3os");
  }

  if (is_true(any(lambda1pfsx <= lambda1osx))) {
    stop("lambda1pfs must be greater than lambda1os");
  }

  if (is_true(any(lambda2pfsx <= lambda2osx))) {
    stop("lambda2pfs must be greater than lambda2os");
  }

  if (is_true(any(lambda3pfsx <= lambda3osx))) {
    stop("lambda3pfs must be greater than lambda3os");
  }

  if (is_true(any(gamma1pfsx < gamma1osx))) {
    stop("gamma1pfs must be greater than or equal to gamma1os");
  }

  if (is_true(any(gamma2pfsx < gamma2osx))) {
    stop("gamma2pfs must be greater than or equal to gamma2os");
  }

  if (is_true(any(gamma3pfsx < gamma3osx))) {
    stop("gamma3pfs must be greater than or equal to gamma3os");
  }

  for (int j=0; j<nstrata; ++j) {
    Range jj = Range(j*nints, (j+1)*nints-1);

    NumericVector lam1pfs = lambda1pfsx[jj];
    NumericVector lam1os = lambda1osx[jj];

    List a1 = hazard_pdcpp(piecewiseSurvivalTime,
                           lam1pfs, lam1os, rho_pd_os);

    NumericVector time1pd = a1["piecewiseSurvivalTime"];
    NumericVector lam1pd = a1["hazard_pd"];

    NumericVector lam2pfs = lambda2pfsx[jj];
    NumericVector lam2os = lambda2osx[jj];

    List a2 = hazard_pdcpp(piecewiseSurvivalTime,
                           lam2pfs, lam2os, rho_pd_os);

    NumericVector time2pd = a2["piecewiseSurvivalTime"];
    NumericVector lam2pd = a2["hazard_pd"];

    NumericVector lam3pfs = lambda3pfsx[jj];
    NumericVector lam3os = lambda3osx[jj];

    List a3 = hazard_pdcpp(piecewiseSurvivalTime,
                           lam3pfs, lam3os, rho_pd_os);

    NumericVector time3pd = a3["piecewiseSurvivalTime"];
    NumericVector lam3pd = a3["hazard_pd"];

    for (int k=0; k<nints2; ++k) {
      tau1pd[j*nints2 + k] = time1pd[k];
      tau2pd[j*nints2 + k] = time2pd[k];
      tau3pd[j*nints2 + k] = time3pd[k];
      lambda1pd[j*nints + k] = lam1pd[k];
      lambda2pd[j*nints + k] = lam2pd[k];
      lambda3pd[j*nints + k] = lam3pd[k];
    }
  }

  gamma1pd = gamma1pfsx - gamma1osx;
  gamma2pd = gamma2pfsx - gamma2osx;
  gamma3pd = gamma3pfsx - gamma3osx;


  if (n == NA_INTEGER) {
    stop("n must be provided");
  }

  if (n <= 0) {
    stop("n must be positive");
  }

  if (fixedFollowup) {
    if (R_isnancpp(followupTime)) {
      stop("followupTime must be provided for fixed follow-up");
    }

    if (followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }
  }

  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }

  if (rho2 < 0) {
    stop("rho2 must be non-negative");
  }


  if (maxNumberOfIterations < 1) {
    stop("maxNumberOfIterations must be a positive integer");
  }

  if (maxNumberOfRawDatasetsPerStage < 0) {
    stop("maxNumberOfRawDatasetsPerStage must be a non-negative integer");
  }


  // declare variables
  int nevents1, nevents2, nstages;
  int accruals1, accruals2, accruals3, totalAccruals;
  int events1e1, events2e1, events3e1, totalEventse1;
  int events1e2, events2e2, events3e2, totalEventse2;
  int dropouts1e1, dropouts2e1, dropouts3e1, totalDropoutse1;
  int dropouts1e2, dropouts2e2, dropouts3e2, totalDropoutse2;
  int index1=0, index2=0;

  double enrollt, u, u1, u2, time;
  double uscore13, uscore23, uscore12;
  double vscore13, vscore23, vscore12;


  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n);

  NumericVector arrivalTime(n), survivalTime1(n), survivalTime2(n),
  dropoutTime1(n), dropoutTime2(n), timeUnderObservation1(n),
  timeUnderObservation2(n), totalTime1(n), totalTime2(n),
  totalt1(n), totalt2(n);

  LogicalVector event1(n), event2(n), dropoutEvent1(n), dropoutEvent2(n),
  event1ac(n), event2ac(n);


  // stratum information
  IntegerVector b1(nstrata), b2(nstrata), b3(nstrata);
  IntegerVector n1(nstrata), n2(nstrata), n3(nstrata);
  NumericVector nt13(nstrata), nt23(nstrata), nt12(nstrata);

  // original copy of n1, n2, and n3 when looping over the endpoints
  IntegerVector n1x(nstrata), n2x(nstrata), n3x(nstrata);

  // hazardRatioH0 adjusted at risk for calculating the log-rank statistic
  NumericVector n13a(nstrata), n23a(nstrata), n12a(nstrata);
  NumericVector nt13a(nstrata), nt23a(nstrata), nt12a(nstrata);

  NumericVector km13(nstrata), km23(nstrata), km12(nstrata);
  NumericVector w13(nstrata), w23(nstrata), w12(nstrata);
  NumericVector cumStratumFraction = cumsum(stratumFraction);

  // within-stratum hazard rates
  NumericVector time1pd(nints2), time2pd(nints2), time3pd(nints2);
  NumericVector lam1pd(nints2), lam2pd(nints2), lam3pd(nints2);
  NumericVector lam1os(nints), lam2os(nints), lam3os(nints);
  NumericVector gam1pd(nints), gam2pd(nints), gam3pd(nints);
  NumericVector gam1os(nints), gam2os(nints), gam3os(nints);


  // stage-wise information
  IntegerVector niter(kMax);
  NumericVector analysisTime(kMax);


  // cache for the patient-level raw data to extract
  int nrow1 = n*kMax*maxNumberOfRawDatasetsPerStage;

  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stageNumberx(nrow1);
  NumericVector analysisTimex(nrow1);
  IntegerVector subjectIdx(nrow1);
  NumericVector arrivalTimex(nrow1);
  IntegerVector stratumx(nrow1);
  IntegerVector treatmentGroupx(nrow1);
  NumericVector survivalTime1x(nrow1);
  NumericVector survivalTime2x(nrow1);
  NumericVector dropoutTime1x(nrow1);
  NumericVector dropoutTime2x(nrow1);
  NumericVector timeUnderObservation1x(nrow1);
  NumericVector timeUnderObservation2x(nrow1);
  LogicalVector event1x(nrow1);
  LogicalVector event2x(nrow1);
  LogicalVector dropoutEvent1x(nrow1);
  LogicalVector dropoutEvent2x(nrow1);


  // cache for the simulation-level summary data to extract
  int nrow2 = kMax*maxNumberOfIterations*2;

  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  LogicalVector eventsNotAchievedy(nrow2);
  IntegerVector stageNumbery(nrow2);
  NumericVector analysisTimey(nrow2);
  IntegerVector accruals1y(nrow2);
  IntegerVector accruals2y(nrow2);
  IntegerVector accruals3y(nrow2);
  IntegerVector totalAccrualsy(nrow2);
  IntegerVector endpointy(nrow2);
  IntegerVector events1y(nrow2);
  IntegerVector events2y(nrow2);
  IntegerVector events3y(nrow2);
  IntegerVector totalEventsy(nrow2);
  IntegerVector dropouts1y(nrow2);
  IntegerVector dropouts2y(nrow2);
  IntegerVector dropouts3y(nrow2);
  IntegerVector totalDropoutsy(nrow2);
  NumericVector logRankStatistic13y(nrow2);
  NumericVector logRankStatistic23y(nrow2);
  NumericVector logRankStatistic12y(nrow2);


  // set up random seed
  if (seed != NA_INTEGER) {
    set_seed(seed);
  }


  // simulation
  for (int iter=0; iter<maxNumberOfIterations; ++iter) {

    b1.fill(allocation1);
    b2.fill(allocation2);
    b3.fill(allocation3);

    enrollt = 0;
    for (int i=0; i<n; ++i) {

      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);
      arrivalTime[i] = enrollt;

      // generate stratum information
      u = R::runif(0,1);
      int j;
      for (j=0; j<nstrata; ++j) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }

      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+b3[j]+0.0)) {
        treatmentGroup[i] = 1;
        --b1[j];
      } else if (u <= (b1[j]+b2[j])/(b1[j]+b2[j]+b3[j]+0.0)) {
        treatmentGroup[i] = 2;
        --b2[j];
      } else {
        treatmentGroup[i] = 3;
        --b3[j];
      }

      // start a new block after depleting the current block
      if (b1[j]+b2[j]+b3[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
        b3[j] = allocation3;
      }

      // stratum-specific hazard rates for event and dropout
      Range jj = Range(j*nints, (j+1)*nints-1);
      Range jj2 = Range(j*nints2, (j+1)*nints2-1);

      time1pd = tau1pd[jj2];
      time2pd = tau2pd[jj2];
      time3pd = tau3pd[jj2];

      lam1pd = lambda1pd[jj2];
      lam2pd = lambda2pd[jj2];
      lam3pd = lambda3pd[jj2];

      lam1os = lambda1osx[jj];
      lam2os = lambda2osx[jj];
      lam3os = lambda3osx[jj];

      gam1pd = gamma1pd[jj];
      gam2pd = gamma2pd[jj];
      gam3pd = gamma3pd[jj];

      gam1os = gamma1osx[jj];
      gam2os = gamma2osx[jj];
      gam3os = gamma3osx[jj];

      // standard bivariate normal with correlation rho_pd_os
      u1 = R::rnorm(0,1);
      u2 = R::rnorm(rho_pd_os*u1, sqrt(1-rho_pd_os*rho_pd_os));

      // transform to uniform
      u1 = R::pnorm(u1, 0, 1, 1, 0);
      u2 = R::pnorm(u2, 0, 1, 1, 0);

      // generate survival times
      if (treatmentGroup[i]==1) {
        survivalTime1[i] = qtpwexpcpp1(u1, time1pd, lam1pd, 0, 1, 0);
        survivalTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,lam1os,0,1,0);
      } else if (treatmentGroup[i]==2) {
        survivalTime1[i] = qtpwexpcpp1(u1, time2pd, lam2pd, 0, 1, 0);
        survivalTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,lam2os,0,1,0);
      } else {
        survivalTime1[i] = qtpwexpcpp1(u1, time3pd, lam3pd, 0, 1, 0);
        survivalTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,lam3os,0,1,0);
      }
      // PFS includes death
      survivalTime1[i] = std::min(survivalTime1[i], survivalTime2[i]);


      // generate dropout times
      u1 = R::runif(0,1);
      u2 = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime1[i] = qtpwexpcpp1(u1,piecewiseSurvivalTime,gam1pd,0,1,0);
        dropoutTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,gam1os,0,1,0);
      } else if (treatmentGroup[i]==2) {
        dropoutTime1[i] = qtpwexpcpp1(u1,piecewiseSurvivalTime,gam2pd,0,1,0);
        dropoutTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,gam2os,0,1,0);
      } else {
        dropoutTime1[i] = qtpwexpcpp1(u1,piecewiseSurvivalTime,gam3pd,0,1,0);
        dropoutTime2[i] = qtpwexpcpp1(u2,piecewiseSurvivalTime,gam3os,0,1,0);
      }
      // whatever censors OS will also censor PFS
      dropoutTime1[i] = std::min(dropoutTime1[i], dropoutTime2[i]);


      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime1[i] <= dropoutTime1[i] &&
            survivalTime1[i] <= followupTime) {
          timeUnderObservation1[i] = survivalTime1[i];
          event1[i] = 1;
          dropoutEvent1[i] = 0;
        } else if (dropoutTime1[i] <= survivalTime1[i] &&
          dropoutTime1[i] <= followupTime) {
          timeUnderObservation1[i] = dropoutTime1[i];
          event1[i] = 0;
          dropoutEvent1[i] = 1;
        } else {
          timeUnderObservation1[i] = followupTime;
          event1[i] = 0;
          dropoutEvent1[i] = 0;
        }
      } else {
        if (survivalTime1[i] <= dropoutTime1[i]) {
          timeUnderObservation1[i] = survivalTime1[i];
          event1[i] = 1;
          dropoutEvent1[i] = 0;
        } else {
          timeUnderObservation1[i] = dropoutTime1[i];
          event1[i] = 0;
          dropoutEvent1[i] = 1;
        }
      }

      totalTime1[i] = arrivalTime[i] + timeUnderObservation1[i];


      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime2[i] <= dropoutTime2[i] &&
            survivalTime2[i] <= followupTime) {
          timeUnderObservation2[i] = survivalTime2[i];
          event2[i] = 1;
          dropoutEvent2[i] = 0;
        } else if (dropoutTime2[i] <= survivalTime2[i] &&
          dropoutTime2[i] <= followupTime) {
          timeUnderObservation2[i] = dropoutTime2[i];
          event2[i] = 0;
          dropoutEvent2[i] = 1;
        } else {
          timeUnderObservation2[i] = followupTime;
          event2[i] = 0;
          dropoutEvent2[i] = 0;
        }
      } else {
        if (survivalTime2[i] <= dropoutTime2[i]) {
          timeUnderObservation2[i] = survivalTime2[i];
          event2[i] = 1;
          dropoutEvent2[i] = 0;
        } else {
          timeUnderObservation2[i] = dropoutTime2[i];
          event2[i] = 0;
          dropoutEvent2[i] = 1;
        }
      }

      totalTime2[i] = arrivalTime[i] + timeUnderObservation2[i];

    }




    // find the analysis time for each stage
    if (useEvents) {
      event1ac = event1 & ((treatmentGroup==1) | (treatmentGroup==3));
      nevents1 = sum(event1ac);
      totalt1 = stl_sort(totalTime1[event1ac]);

      event2ac = event2 & ((treatmentGroup==1) | (treatmentGroup==3));
      nevents2 = sum(event2ac);
      totalt2 = stl_sort(totalTime2[event2ac]);

      int j1 = kMaxpfsx, j2 = kMax - kMaxpfsx;

      // PFS looks
      if (kMaxpfsx > 0) {
        for (j1=0; j1<kMaxpfsx; ++j1) {
          if (plannedEvents[j1] >= nevents1) {
            break;
          }
        }

        if (j1==kMaxpfsx) { // total number of PFS events exceeds planned
          for (int k=0; k<kMaxpfsx; ++k) {
            analysisTime[k] = totalt1[plannedEvents[k]-1] + 1e-12;
          }
        } else {
          for (int k=0; k<=j1; ++k) {
            if (k < j1) {
              analysisTime[k] = totalt1[plannedEvents[k]-1] + 1e-12;
            } else {
              analysisTime[k] = totalt1[nevents1-1] + 1e-12;
            }
          }
        }
      }

      // OS looks
      NumericVector analysisTime2(kMax - kMaxpfsx);

      if (kMax > kMaxpfsx) {
        for (j2=0; j2<kMax-kMaxpfsx; ++j2) {
          if (plannedEvents[kMaxpfsx+j2] >= nevents2) {
            break;
          }
        }

        if (j2==kMax-kMaxpfsx) { // total number of OS events exceeds planned
          for (int k=0; k<kMax-kMaxpfsx; ++k) {
            analysisTime2[k] = totalt2[plannedEvents[kMaxpfsx+k]-1] + 1e-12;
          }
        } else {
          for (int k=0; k<=j2; ++k) {
            if (k < j2) {
              analysisTime2[k] = totalt2[plannedEvents[kMaxpfsx+k]-1] + 1e-12;
            } else {
              analysisTime2[k] = totalt2[nevents2-1] + 1e-12;
            }
          }
        }
      }

      // determine the number of looks and timing of the looks
      if (kMaxpfsx == 0) { // all looks based on OS events
        if (j2 == kMax - kMaxpfsx) {
          nstages = kMax - kMaxpfsx;
        } else {
          nstages = j2 + 1;
        }

        for (int k=0; k<nstages; ++k) {
          analysisTime[k] = analysisTime2[k];
        }
      } else if (kMax == kMaxpfsx) { // all looks based on PFS events
        if (j1 == kMaxpfsx) {
          nstages = kMaxpfsx;
        } else {
          nstages = j1 + 1;
        }
      } else {
        if (analysisTime2[kMax-kMaxpfsx-1] > analysisTime[kMaxpfsx-1]) {
          // only OS looks that occur after the last PFS look contribute
          int l = static_cast<int>(which_max(analysisTime2 >
                                               analysisTime[kMaxpfsx-1]));
          nstages = kMax-l;
          for (int k=kMaxpfsx; k<kMax-l; ++k) {
            analysisTime[k] = analysisTime2[k-kMaxpfsx+l];
          }
        } else {
          if (j1 == kMaxpfsx) {
            nstages = kMaxpfsx;
          } else {
            nstages = j1 + 1;
          }
        }
      }

      // whether the target PFS and OS events are achieved
      if (kMaxpfsx > 0 && nevents1 < plannedEvents[kMaxpfsx-1]) {
        eventsNotAchieved = 1;
      } else if (kMaxpfsx < kMax && nevents2 < plannedEvents[kMax-1]) {
        eventsNotAchieved = 1;
      } else {
        eventsNotAchieved = 0;
      }
    } else { // looks based on calendar time
      nstages = kMax;
      analysisTime = clone(plannedTime);
      eventsNotAchieved = 0;
    }


    // construct the log-rank test statistic at each stage
    for (int k=0; k<nstages; ++k) {
      time = analysisTime[k];

      n1x.fill(0);  // number of subjects in each stratum by treatment
      n2x.fill(0);
      n3x.fill(0);

      events1e1 = 0;
      events2e1 = 0;
      events3e1 = 0;

      dropouts1e1 = 0;
      dropouts2e1 = 0;
      dropouts3e1 = 0;

      events1e2 = 0;
      events2e2 = 0;
      events3e2 = 0;

      dropouts1e2 = 0;
      dropouts2e2 = 0;
      dropouts3e2 = 0;

      // censor at analysis time
      for (int i=0; i<n; ++i) {
        int h = stratum[i]-1;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation1[i] = time - arrivalTime[i];
          event1[i] = 0;
          dropoutEvent1[i] = 0;

          timeUnderObservation2[i] = time - arrivalTime[i];
          event2[i] = 0;
          dropoutEvent2[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            ++n1x[h];
          } else if (treatmentGroup[i]==2) {
            ++n2x[h];
          } else {
            ++n3x[h];
          }


          // censored time for endpoint 1
          if (fixedFollowup) {
            if (arrivalTime[i] + survivalTime1[i] <= time &&
                survivalTime1[i] <= dropoutTime1[i] &&
                survivalTime1[i] <= followupTime) {
              timeUnderObservation1[i] = survivalTime1[i];
              event1[i] = 1;
              dropoutEvent1[i] = 0;
            } else if (arrivalTime[i] + dropoutTime1[i] <= time &&
              dropoutTime1[i] <= survivalTime1[i] &&
              dropoutTime1[i] <= followupTime) {
              timeUnderObservation1[i] = dropoutTime1[i];
              event1[i] = 0;
              dropoutEvent1[i] = 1;
            } else if (arrivalTime[i] + followupTime <= time &&
              followupTime <= survivalTime1[i] &&
              followupTime <= dropoutTime1[i]) {
              timeUnderObservation1[i] = followupTime;
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            } else {
              timeUnderObservation1[i] = time - arrivalTime[i];
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime1[i] <= time &&
                survivalTime1[i] <= dropoutTime1[i]) {
              timeUnderObservation1[i] = survivalTime1[i];
              event1[i] = 1;
              dropoutEvent1[i] = 0;
            } else if (arrivalTime[i] + dropoutTime1[i] <= time &&
              dropoutTime1[i] <= survivalTime1[i]) {
              timeUnderObservation1[i] = dropoutTime1[i];
              event1[i] = 0;
              dropoutEvent1[i] = 1;
            } else {
              timeUnderObservation1[i] = time - arrivalTime[i];
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            }
          }

          if (treatmentGroup[i]==1 && event1[i]) ++events1e1;
          if (treatmentGroup[i]==2 && event1[i]) ++events2e1;
          if (treatmentGroup[i]==3 && event1[i]) ++events3e1;
          if (treatmentGroup[i]==1 && dropoutEvent1[i]) ++dropouts1e1;
          if (treatmentGroup[i]==2 && dropoutEvent1[i]) ++dropouts2e1;
          if (treatmentGroup[i]==3 && dropoutEvent1[i]) ++dropouts3e1;


          // censored time for endpoint 2
          if (fixedFollowup) {
            if (arrivalTime[i] + survivalTime2[i] <= time &&
                survivalTime2[i] <= dropoutTime2[i] &&
                survivalTime2[i] <= followupTime) {
              timeUnderObservation2[i] = survivalTime2[i];
              event2[i] = 1;
              dropoutEvent2[i] = 0;
            } else if (arrivalTime[i] + dropoutTime2[i] <= time &&
              dropoutTime2[i] <= survivalTime2[i] &&
              dropoutTime2[i] <= followupTime) {
              timeUnderObservation2[i] = dropoutTime2[i];
              event2[i] = 0;
              dropoutEvent2[i] = 1;
            } else if (arrivalTime[i] + followupTime <= time &&
              followupTime <= survivalTime2[i] &&
              followupTime <= dropoutTime2[i]) {
              timeUnderObservation2[i] = followupTime;
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            } else {
              timeUnderObservation2[i] = time - arrivalTime[i];
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime2[i] <= time &&
                survivalTime2[i] <= dropoutTime2[i]) {
              timeUnderObservation2[i] = survivalTime2[i];
              event2[i] = 1;
              dropoutEvent2[i] = 0;
            } else if (arrivalTime[i] + dropoutTime2[i] <= time &&
              dropoutTime2[i] <= survivalTime2[i]) {
              timeUnderObservation2[i] = dropoutTime2[i];
              event2[i] = 0;
              dropoutEvent2[i] = 1;
            } else {
              timeUnderObservation2[i] = time - arrivalTime[i];
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            }
          }

          if (treatmentGroup[i]==1 && event2[i]) ++events1e2;
          if (treatmentGroup[i]==2 && event2[i]) ++events2e2;
          if (treatmentGroup[i]==3 && event2[i]) ++events3e2;
          if (treatmentGroup[i]==1 && dropoutEvent2[i]) ++dropouts1e2;
          if (treatmentGroup[i]==2 && dropoutEvent2[i]) ++dropouts2e2;
          if (treatmentGroup[i]==3 && dropoutEvent2[i]) ++dropouts3e2;
        }
      }


      // add raw data to output
      if (niter[k] < maxNumberOfRawDatasetsPerStage) {
        for (int i=0; i<n; ++i) {
          iterationNumberx[index1] = iter+1;
          stageNumberx[index1] = k+1;
          analysisTimex[index1] = time;
          subjectIdx[index1] = i+1;
          arrivalTimex[index1] = arrivalTime[i];
          stratumx[index1] = stratum[i];
          treatmentGroupx[index1] = treatmentGroup[i];

          survivalTime1x[index1] = survivalTime1[i];
          dropoutTime1x[index1] = dropoutTime1[i];
          timeUnderObservation1x[index1] = timeUnderObservation1[i];
          event1x[index1] = event1[i];
          dropoutEvent1x[index1] = dropoutEvent1[i];

          survivalTime2x[index1] = survivalTime2[i];
          dropoutTime2x[index1] = dropoutTime2[i];
          timeUnderObservation2x[index1] = timeUnderObservation2[i];
          event2x[index1] = event2[i];
          dropoutEvent2x[index1] = dropoutEvent2[i];

          ++index1;
        }

        // update the number of stage k dataset to extract
        ++niter[k];
      }


      // number of accrued patients and total number of events
      accruals1 = sum(n1x);
      accruals2 = sum(n2x);
      accruals3 = sum(n3x);
      totalAccruals = accruals1 + accruals2 + accruals3;

      totalEventse1 = events1e1 + events2e1 + events3e1;
      totalDropoutse1 = dropouts1e1 + dropouts2e1 + dropouts3e1;

      totalEventse2 = events1e2 + events2e2 + events3e2;
      totalDropoutse2 = dropouts1e2 + dropouts2e2 + dropouts3e2;


      for (int endpoint=1; endpoint<=2; ++endpoint) {
        n1 = clone(n1x);
        n2 = clone(n2x);
        n3 = clone(n3x);

        double hazardRatioH013, hazardRatioH023, hazardRatioH012;
        if (endpoint == 1) {
          hazardRatioH013 = hazardRatioH013pfs;
          hazardRatioH023 = hazardRatioH023pfs;
          hazardRatioH012 = hazardRatioH012pfs;
        } else {
          hazardRatioH013 = hazardRatioH013os;
          hazardRatioH023 = hazardRatioH023os;
          hazardRatioH012 = hazardRatioH012os;
        }


        // order the data by time under observation
        NumericVector timeUnderObservationSorted;
        IntegerVector sortedIndex;
        LogicalVector eventSorted;
        if (endpoint == 1) {
          timeUnderObservationSorted = stl_sort(timeUnderObservation1);
          sortedIndex = match(timeUnderObservationSorted,
                              timeUnderObservation1);
          sortedIndex = sortedIndex - 1;
          eventSorted = event1[sortedIndex];
        } else {
          timeUnderObservationSorted = stl_sort(timeUnderObservation2);
          sortedIndex = match(timeUnderObservationSorted,
                              timeUnderObservation2);
          sortedIndex = sortedIndex - 1;
          eventSorted = event2[sortedIndex];
        }

        IntegerVector stratumSorted = stratum[sortedIndex];
        IntegerVector treatmentGroupSorted = treatmentGroup[sortedIndex];

        LogicalVector subSorted = (timeUnderObservationSorted > 0);
        stratumSorted = stratumSorted[subSorted];
        treatmentGroupSorted = treatmentGroupSorted[subSorted];
        eventSorted = eventSorted[subSorted];
        int nsubSorted = static_cast<int>(eventSorted.size());

        // calculate the stratified log-rank test
        uscore13 = 0;
        vscore13 = 0;
        uscore23 = 0;
        vscore23 = 0;
        uscore12 = 0;
        vscore12 = 0;
        km13.fill(1);
        km23.fill(1);
        km12.fill(1);
        for (int i=0; i<nsubSorted; ++i) {
          int h = stratumSorted[i] - 1;
          nt13[h] = n1[h] + n3[h];
          nt23[h] = n2[h] + n3[h];
          nt12[h] = n1[h] + n2[h];

          n13a[h] = n1[h]*hazardRatioH013;
          n23a[h] = n2[h]*hazardRatioH023;
          n12a[h] = n1[h]*hazardRatioH012;

          nt13a[h] = n13a[h] + n3[h];
          nt23a[h] = n23a[h] + n3[h];
          nt12a[h] = n12a[h] + n2[h];


          if (eventSorted[i] && (treatmentGroupSorted[i]==1 ||
              treatmentGroupSorted[i]==3)) {
            w13[h] = pow(km13[h], rho1)*pow(1-km13[h], rho2);
            uscore13 += w13[h]*((treatmentGroupSorted[i]==1)
                                  - n13a[h]/nt13a[h]);
            vscore13 += w13[h]*w13[h]*n13a[h]*n3[h]/(nt13a[h]*nt13a[h]);
            km13[h] *= (1-1/nt13[h]); // update km estimate
          }

          if (eventSorted[i] && (treatmentGroupSorted[i]==2 ||
              treatmentGroupSorted[i]==3)) {
            w23[h] = pow(km23[h], rho1)*pow(1-km23[h], rho2);
            uscore23 += w23[h]*((treatmentGroupSorted[i]==2)
                                  - n23a[h]/nt23a[h]);
            vscore23 += w23[h]*w23[h]*n23a[h]*n3[h]/(nt23a[h]*nt23a[h]);
            km23[h] *= (1-1/nt23[h]); // update km estimate
          }

          if (eventSorted[i] && (treatmentGroupSorted[i]==1 ||
              treatmentGroupSorted[i]==2)) {
            w12[h] = pow(km12[h], rho1)*pow(1-km12[h], rho2);
            uscore12 += w12[h]*((treatmentGroupSorted[i]==1)
                                  - n12a[h]/nt12a[h]);
            vscore12 += w12[h]*w12[h]*n12a[h]*n2[h]/(nt12a[h]*nt12a[h]);
            km12[h] *= (1-1/nt12[h]); // update km estimate
          }

          // reduce the risk set
          if (treatmentGroupSorted[i]==1) {
            --n1[h];
          } else if (treatmentGroupSorted[i]==2) {
            --n2[h];
          } else {
            --n3[h];
          }
        }



        // add summary data to output
        iterationNumbery[index2] = iter+1;
        eventsNotAchievedy[index2] = eventsNotAchieved;
        stageNumbery[index2] = k+1;
        analysisTimey[index2] = time;
        accruals1y[index2] = accruals1;
        accruals2y[index2] = accruals2;
        accruals3y[index2] = accruals3;
        totalAccrualsy[index2] = totalAccruals;
        endpointy[index2] = endpoint;

        if (endpoint == 1) {
          events1y[index2] = events1e1;
          events2y[index2] = events2e1;
          events3y[index2] = events3e1;
          totalEventsy[index2] = totalEventse1;
          dropouts1y[index2] = dropouts1e1;
          dropouts2y[index2] = dropouts2e1;
          dropouts3y[index2] = dropouts3e1;
          totalDropoutsy[index2] = totalDropoutse1;
        } else {
          events1y[index2] = events1e2;
          events2y[index2] = events2e2;
          events3y[index2] = events3e2;
          totalEventsy[index2] = totalEventse2;
          dropouts1y[index2] = dropouts1e2;
          dropouts2y[index2] = dropouts2e2;
          dropouts3y[index2] = dropouts3e2;
          totalDropoutsy[index2] = totalDropoutse2;
        }

        logRankStatistic13y[index2] = uscore13/sqrt(vscore13);
        logRankStatistic23y[index2] = uscore23/sqrt(vscore23);
        logRankStatistic12y[index2] = uscore12/sqrt(vscore12);

        ++index2;

      } // end of endpoint

    } // end of stage

  } // end of iteration


  // simulation summary data set
  LogicalVector sub2 = !is_na(iterationNumbery);
  iterationNumbery = iterationNumbery[sub2];
  eventsNotAchievedy = eventsNotAchievedy[sub2];
  stageNumbery = stageNumbery[sub2];
  analysisTimey = analysisTimey[sub2];
  accruals1y = accruals1y[sub2];
  accruals2y = accruals2y[sub2];
  accruals3y = accruals3y[sub2];
  totalAccrualsy = totalAccrualsy[sub2];
  endpointy = endpointy[sub2];
  events1y = events1y[sub2];
  events2y = events2y[sub2];
  events3y = events3y[sub2];
  totalEventsy = totalEventsy[sub2];
  dropouts1y = dropouts1y[sub2];
  dropouts2y = dropouts2y[sub2];
  dropouts3y = dropouts3y[sub2];
  totalDropoutsy = totalDropoutsy[sub2];
  logRankStatistic13y = logRankStatistic13y[sub2];
  logRankStatistic23y = logRankStatistic23y[sub2];
  logRankStatistic12y = logRankStatistic12y[sub2];

  DataFrame sumdata = DataFrame::create(
    _["iterationNumber"] = iterationNumbery,
    _["eventsNotAchieved"] = eventsNotAchievedy,
    _["stageNumber"] = stageNumbery,
    _["analysisTime"] = analysisTimey,
    _["accruals1"] = accruals1y,
    _["accruals2"] = accruals2y,
    _["accruals3"] = accruals3y,
    _["totalAccruals"] = totalAccrualsy,
    _["endpoint"] = endpointy,
    _["events1"] = events1y,
    _["events2"] = events2y,
    _["events3"] = events3y,
    _["totalEvents"] = totalEventsy,
    _["dropouts1"] = dropouts1y,
    _["dropouts2"] = dropouts2y,
    _["dropouts3"] = dropouts3y,
    _["totalDropouts"] = totalDropoutsy,
    _["logRankStatistic13"] = logRankStatistic13y,
    _["logRankStatistic23"] = logRankStatistic23y,
    _["logRankStatistic12"] = logRankStatistic12y);


  List result;

  if (maxNumberOfRawDatasetsPerStage > 0) {
    LogicalVector sub1 = !is_na(iterationNumberx);
    iterationNumberx = iterationNumberx[sub1];
    stageNumberx = stageNumberx[sub1];
    analysisTimex = analysisTimex[sub1];
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    survivalTime1x = survivalTime1x[sub1];
    dropoutTime1x = dropoutTime1x[sub1];
    timeUnderObservation1x = timeUnderObservation1x[sub1];
    event1x = event1x[sub1];
    dropoutEvent1x = dropoutEvent1x[sub1];
    survivalTime2x = survivalTime2x[sub1];
    dropoutTime2x = dropoutTime2x[sub1];
    timeUnderObservation2x = timeUnderObservation2x[sub1];
    event2x = event2x[sub1];
    dropoutEvent2x = dropoutEvent2x[sub1];

    DataFrame rawdata = DataFrame::create(
      _["iterationNumber"] = iterationNumberx,
      _["stageNumber"] = stageNumberx,
      _["analysisTime"] = analysisTimex,
      _["subjectId"] = subjectIdx,
      _["arrivalTime"] = arrivalTimex,
      _["stratum"] = stratumx,
      _["treatmentGroup"] = treatmentGroupx,
      _["survivalTime1"] = survivalTime1x,
      _["dropoutTime1"] = dropoutTime1x,
      _["timeUnderObservation1"] = timeUnderObservation1x,
      _["event1"] = event1x,
      _["dropoutEvent1"] = dropoutEvent1x,
      _["survivalTime2"] = survivalTime2x,
      _["dropoutTime2"] = dropoutTime2x,
      _["timeUnderObservation2"] = timeUnderObservation2x,
      _["event2"] = event2x,
      _["dropoutEvent2"] = dropoutEvent2x);

    result = List::create(_["sumdata"] = sumdata,
                          _["rawdata"] = rawdata);
  } else {
    result = List::create(_["sumdata"] = sumdata);
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
//' @param seed Random seed for reproducibility. If not specified,
//'   the current R environment seed is used.
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
//'     - \code{mnStatistic}: The Mantel-Haenszel test Z-statistic for
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
//'     - \code{ptfu1Time}:The underlying assessment time for the
//'       binary endpoint for the subject.
//'
//'     - \code{timeUnderObservation}: The time under observation
//'       since randomization for the binary endpoint for the subject.
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
//'   accrualTime = 0:8,
//'   accrualIntensity = c(((1:8) - 0.5)/8, 1)*40,
//'   piecewiseSurvivalTime = c(0,12,36),
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
List binary_tte_sim(
    const int kMax1 = 1,
    const int kMax2 = 1,
    const double riskDiffH0 = 0,
    const double hazardRatioH0 = 1,
    const int allocation1 = 1,
    const int allocation2 = 1,
    const NumericVector& accrualTime = 0,
    const NumericVector& accrualIntensity = NA_REAL,
    const NumericVector& piecewiseSurvivalTime = 0,
    const NumericVector& stratumFraction = 1,
    const double globalOddsRatio = 1,
    const NumericVector& pi1 = NA_REAL,
    const NumericVector& pi2 = NA_REAL,
    const NumericVector& lambda1 = NA_REAL,
    const NumericVector& lambda2 = NA_REAL,
    const NumericVector& gamma1 = 0,
    const NumericVector& gamma2 = 0,
    const NumericVector& delta1 = 0,
    const NumericVector& delta2 = 0,
    const double upper1 = NA_REAL,
    const double upper2 = NA_REAL,
    const int n = NA_INTEGER,
    const NumericVector& plannedTime = NA_REAL,
    const IntegerVector& plannedEvents = NA_INTEGER,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = NA_INTEGER) {

  // check input parameters
  int nstrata = static_cast<int>(stratumFraction.size());
  int nints = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nints;

  NumericVector pi1x(nstrata), pi2x(nstrata);
  NumericVector lambda1x(nsi), lambda2x(nsi);
  NumericVector gamma1x(nsi), gamma2x(nsi);
  NumericVector delta1x(nsi), delta2x(nsi);

  bool eventsNotAchieved;

  if (kMax1 < 1) {
    stop("kMax1 must be a positive integer");
  }

  if (kMax2 < 1) {
    stop("kMax2 must be a positive integer");
  }


  if ((riskDiffH0 <= -1) || (riskDiffH0 >= 1)) {
    stop("riskDiffH0 must lie between -1 and 1");
  }

  if (hazardRatioH0 <= 0) {
    stop("hazardRatioH0 must be positive");
  }


  if (allocation1 < 1) {
    stop("allocation1 must be a positive integer");
  }

  if (allocation2 < 1) {
    stop("allocation2 must be a positive integer");
  }


  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }


  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nints > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }


  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }


  if (globalOddsRatio <= 0) {
    stop("globalOddsRatio must be positive");
  }


  if (is_true(any(is_na(pi1)))) {
    stop("pi1 must be provided");
  }

  if (is_true(any(is_na(pi2)))) {
    stop("pi2 must be provided");
  }

  if (is_true(any((pi1 <= 0) | (pi1 >= 1)))) {
    stop("pi1 must lie between 0 and 1");
  }

  if (is_true(any((pi2 <= 0) | (pi2 >= 1)))) {
    stop("pi2 must lie between 0 and 1");
  }

  if (pi1.size() == 1) {
    pi1x = rep(pi1, nstrata);
  } else if (pi1.size() == nstrata) {
    pi1x = pi1;
  } else {
    stop("Invalid length for pi1");
  }

  if (pi2.size() == 1) {
    pi2x = rep(pi2, nstrata);
  } else if (pi2.size() == nstrata) {
    pi2x = pi2;
  } else {
    stop("Invalid length for pi2");
  }


  if (is_true(any(is_na(lambda1)))) {
    stop("lambda1 must be provided");
  }

  if (is_true(any(is_na(lambda2)))) {
    stop("lambda2 must be provided");
  }

  if (is_true(any(lambda1 < 0))) {
    stop("lambda1 must be non-negative");
  }

  if (is_true(any(lambda2 < 0))) {
    stop("lambda2 must be non-negative");
  }

  if (is_true(any(gamma1 < 0))) {
    stop("gamma1 must be non-negative");
  }

  if (is_true(any(gamma2 < 0))) {
    stop("gamma2 must be non-negative");
  }

  if (is_true(any(delta1 < 0))) {
    stop("delta1 must be non-negative");
  }

  if (is_true(any(delta2 < 0))) {
    stop("delta2 must be non-negative");
  }

  if (lambda1.size() == 1) {
    lambda1x = rep(lambda1, nsi);
  } else if (lambda1.size() == nints) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }

  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nints) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }

  if (gamma1.size() == 1) {
    gamma1x = rep(gamma1, nsi);
  } else if (gamma1.size() == nints) {
    gamma1x = rep(gamma1, nstrata);
  } else if (gamma1.size() == nsi) {
    gamma1x = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }

  if (gamma2.size() == 1) {
    gamma2x = rep(gamma2, nsi);
  } else if (gamma2.size() == nints) {
    gamma2x = rep(gamma2, nstrata);
  } else if (gamma2.size() == nsi) {
    gamma2x = gamma2;
  } else {
    stop("Invalid length for gamma2");
  }

  if (delta1.size() == 1) {
    delta1x = rep(delta1, nsi);
  } else if (delta1.size() == nints) {
    delta1x = rep(delta1, nstrata);
  } else if (delta1.size() == nsi) {
    delta1x = delta1;
  } else {
    stop("Invalid length for delta1");
  }

  if (delta2.size() == 1) {
    delta2x = rep(delta2, nsi);
  } else if (delta2.size() == nints) {
    delta2x = rep(delta2, nstrata);
  } else if (delta2.size() == nsi) {
    delta2x = delta2;
  } else {
    stop("Invalid length for delta2");
  }


  if (upper1 <= 0) {
    stop("upper1 must be positive");
  }

  if (upper2 <= 0) {
    stop("upper2 must be positive");
  }


  if (n == NA_INTEGER) {
    stop("n must be provided");
  }

  if (n <= 0) {
    stop("n must be positive");
  }


  if (is_true(any(is_na(plannedTime)))) {
    stop("plannedTime must be given for endpoint 1");
  }

  if (plannedTime[0] <= 0) {
    stop("Elements of plannedTime must be positive");
  }

  if (plannedTime.size() != kMax1) {
    stop("Invalid length for plannedTime");
  }

  if (kMax1 > 1 && is_true(any(diff(plannedTime) <= 0))) {
    stop("plannedTime must be increasing");
  }


  if (is_true(any(is_na(plannedEvents)))) {
    stop("plannedEvents must be given for endpoint 2");
  }

  if (plannedEvents[0] <= 0) {
    stop("Elements of plannedEvents must be positive");
  }

  if (plannedEvents.size() != kMax2) {
    stop("Invalid length for plannedEvents");
  }

  if (kMax2 > 1 && is_true(any(diff(plannedEvents) <= 0))) {
    stop("plannedEvents must be increasing");
  }


  if (maxNumberOfIterations < 1) {
    stop("maxNumberOfIterations must be a positive integer");
  }

  if (maxNumberOfRawDatasetsPerStage < 0) {
    stop("maxNumberOfRawDatasetsPerStage must be a non-negative integer");
  }


  // declare variables
  int nevents, nstages1, nstages2;
  int accruals1, accruals2, totalAccruals;
  int events1, events2, totalEvents;
  int dropouts1, dropouts2, totalDropouts;
  int index1x=0, index2x=0, index1y=0, index2y=0;

  double enrollt, u, u1, u2, time, uscore, vscore;


  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n);
  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n);
  NumericVector timeUnderObservation(n), totalTime(n), totalt(n);
  LogicalVector event(n), dropoutEvent(n);

  NumericVector latentResponse(n);
  NumericVector trtDiscTime(n), ptfu1Time(n);
  NumericVector timeUnderObservation1(n);
  LogicalVector responder(n);
  IntegerVector source(n);

  // stratum information
  IntegerVector b1(nstrata), b2(nstrata);
  IntegerVector n1(nstrata), n2(nstrata);

  // original copy of n1 and n2 when looping over the endpoints
  IntegerVector n1x(nstrata), n2x(nstrata);

  // hazardRatioH0 adjusted n1 and nt for calculating the log-rank statistic
  NumericVector nt(nstrata), n1a(nstrata), nta(nstrata);

  NumericVector cumStratumFraction = cumsum(stratumFraction);

  // within-stratum response rates and hazard rates
  double alpha0, alpha1; // parameters on the logit scale for endpoint 1
  NumericVector lam1(nints), lam2(nints);
  NumericVector gam1(nints), gam2(nints);
  NumericVector del1(nints), del2(nints);


  // stage-wise information
  NumericVector analysisTime1(kMax1), analysisTime2(kMax2);
  IntegerVector niter1(kMax1), niter2(kMax2);


  // cache for the patient-level raw data for endpoint 1 (uMRD) to extract
  int nrow1x = n*kMax1*maxNumberOfRawDatasetsPerStage;

  IntegerVector iterationNumber1x = IntegerVector(nrow1x, NA_INTEGER);
  IntegerVector stageNumber1x(nrow1x);
  NumericVector analysisTime1x(nrow1x);
  IntegerVector subjectId1x(nrow1x);
  NumericVector arrivalTime1x(nrow1x);
  IntegerVector stratum1x(nrow1x);
  IntegerVector treatmentGroup1x(nrow1x);
  NumericVector survivalTime1x(nrow1x);
  NumericVector dropoutTime1x(nrow1x);
  NumericVector ptfu1Timex(nrow1x);
  NumericVector timeUnderObservation1x(nrow1x);
  LogicalVector responderx(nrow1x);
  IntegerVector sourcex(nrow1x);


  // cache for the patient-level raw data for endpoint 2 (PFS) to extract
  int nrow2x = n*kMax2*maxNumberOfRawDatasetsPerStage;

  IntegerVector iterationNumber2x = IntegerVector(nrow2x, NA_INTEGER);
  IntegerVector stageNumber2x(nrow2x);
  NumericVector analysisTime2x(nrow2x);
  IntegerVector subjectId2x(nrow2x);
  NumericVector arrivalTime2x(nrow2x);
  IntegerVector stratum2x(nrow2x);
  IntegerVector treatmentGroup2x(nrow2x);
  NumericVector survivalTime2x(nrow2x);
  NumericVector dropoutTime2x(nrow2x);
  NumericVector timeUnderObservation2x(nrow2x);
  LogicalVector eventx(nrow2x);
  LogicalVector dropoutEventx(nrow2x);


  // cache for the simulation-level summary data for endpoint 1 to extract
  int nrow1y = kMax1*maxNumberOfIterations;

  IntegerVector iterationNumber1y = IntegerVector(nrow1y, NA_INTEGER);
  IntegerVector stageNumber1y(nrow1y);
  NumericVector analysisTime1y(nrow1y);
  IntegerVector accruals11y(nrow1y);
  IntegerVector accruals21y(nrow1y);
  IntegerVector totalAccruals1y(nrow1y);
  IntegerVector source1y(nrow1y);
  IntegerVector source2y(nrow1y);
  IntegerVector source3y(nrow1y);
  NumericVector n1y(nrow1y);
  NumericVector n2y(nrow1y);
  NumericVector ny(nrow1y);
  NumericVector y1y(nrow1y);
  NumericVector y2y(nrow1y);
  NumericVector yy(nrow1y);
  NumericVector riskDiffy(nrow1y);
  NumericVector seRiskDiffy(nrow1y);
  NumericVector mhStatisticy(nrow1y);


  int nrow2y = kMax2*maxNumberOfIterations;

  IntegerVector iterationNumber2y = IntegerVector(nrow2y, NA_INTEGER);
  LogicalVector eventsNotAchievedy(nrow2y);
  IntegerVector stageNumber2y(nrow2y);
  NumericVector analysisTime2y(nrow2y);
  IntegerVector accruals1y(nrow2y);
  IntegerVector accruals2y(nrow2y);
  IntegerVector totalAccrualsy(nrow2y);
  IntegerVector events1y(nrow2y);
  IntegerVector events2y(nrow2y);
  IntegerVector totalEventsy(nrow2y);
  IntegerVector dropouts1y(nrow2y);
  IntegerVector dropouts2y(nrow2y);
  IntegerVector totalDropoutsy(nrow2y);
  NumericVector logRankStatisticy(nrow2y);


  // set up random seed
  if (seed != NA_INTEGER) {
    set_seed(seed);
  }

  double theta = globalOddsRatio;
  NumericVector alpha0x = log(pi2x/(1-pi2x));
  NumericVector alpha1x = log(pi1x/(1-pi1x)) - alpha0x;
  double S, b, c, d;

  for (int iter=0; iter<maxNumberOfIterations; ++iter) {

    b1.fill(allocation1);
    b2.fill(allocation2);

    enrollt = 0;
    for (int i=0; i<n; ++i) {

      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);
      arrivalTime[i] = enrollt;

      // generate stratum information
      u = R::runif(0,1);
      int j;
      for (j=0; j<nstrata; ++j) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }

      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+0.0)) {
        treatmentGroup[i] = 1;
        --b1[j];
      } else {
        treatmentGroup[i] = 2;
        --b2[j];
      }

      // start a new block after depleting the current block
      if (b1[j]+b2[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
      }

      // stratum-specific response rates and hazard rates
      alpha0 = alpha0x[j];
      alpha1 = alpha1x[j];

      Range jj = Range(j*nints, (j+1)*nints-1);

      lam1 = lambda1x[jj];
      lam2 = lambda2x[jj];

      gam1 = gamma1x[jj];
      gam2 = gamma2x[jj];

      del1 = delta1x[jj];
      del2 = delta2x[jj];

      // Plackett copula
      u1 = R::runif(0,1);
      u2 = R::runif(0,1);
      S = u2*(1-u2);
      b = theta + S*(theta-1)*(theta-1);
      c = 2*S*(1 + (theta*theta - 1)*u1) + (1 - 2*S)*theta;
      d = sqrt(theta*(theta + 4*S*(theta-1)*(theta-1)*u1*(1-u1)));
      u2 = (c - (1 - 2*u2)*d)/(2*b);

      // back transform
      double location = -alpha0 - alpha1*(treatmentGroup[i] == 1);
      latentResponse[i] = R::qlogis(u1, location, 1, 1, 0);

      // generate survival times
      NumericVector lam = lam2 + (lam1 - lam2)*(treatmentGroup[i] == 1);
      survivalTime[i] = qtpwexpcpp1(u2, piecewiseSurvivalTime, lam, 0, 1, 0);

      // generate dropout times
      u = R::runif(0,1);
      NumericVector gam = gam2 + (gam1 - gam2)*(treatmentGroup[i] == 1);
      dropoutTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, gam, 0, 1, 0);

      // generate treatment discontinuation times
      u = R::runif(0,1);
      NumericVector del = del2 + (del1 - del2)*(treatmentGroup[i] == 1);
      trtDiscTime[i] = qtpwexpcpp1(u, piecewiseSurvivalTime, del, 0, 1, 0);

      double upper = upper2 + (upper1 - upper2)*(treatmentGroup[i] == 1);
      ptfu1Time[i] = std::min(trtDiscTime[i], upper);

      // initial observed time and event indicator for endpoint 2 (TTE)
      if (survivalTime[i] <= dropoutTime[i]) {
        timeUnderObservation[i] = survivalTime[i];
        event[i] = 1;
        dropoutEvent[i] = 0;
      } else {
        timeUnderObservation[i] = dropoutTime[i];
        event[i] = 0;
        dropoutEvent[i] = 1;
      }

      totalTime[i] = arrivalTime[i] + timeUnderObservation[i];
    }


    // find the analysis time for each stage for endpoint 1
    nstages1 = kMax1;
    analysisTime1 = clone(plannedTime);


    // construct the Mantel-Haenszel statistic for endpoint 1 at each stage
    for (int k=0; k<nstages1; ++k) {
      time = analysisTime1[k];

      n1x.fill(0);  // number of subjects in each stratum by treatment
      n2x.fill(0);

      for (int i=0; i<n; ++i) {
        int h = stratum[i]-1;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation1[i] = time - arrivalTime[i];
          responder[i] = NA_LOGICAL;
          source[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            ++n1x[h];
          } else {
            ++n2x[h];
          }

          // censor at analysis time
          if (arrivalTime[i] + ptfu1Time[i] <= time &&
              ptfu1Time[i] <= survivalTime[i] &&
              ptfu1Time[i] <= dropoutTime[i]) {
            timeUnderObservation1[i] = ptfu1Time[i];
            responder[i] = (latentResponse[i] <= 0);
            source[i] = 1;
          } else if (arrivalTime[i] + survivalTime[i] <= time &&
            survivalTime[i] <= dropoutTime[i] &&
            survivalTime[i] <= ptfu1Time[i]) {
            timeUnderObservation1[i] = survivalTime[i];
            responder[i] = 0;
            source[i] = 2;
          } else if (arrivalTime[i] + dropoutTime[i] <= time &&
            dropoutTime[i] <= survivalTime[i] &&
            dropoutTime[i] <= ptfu1Time[i]) {
            timeUnderObservation1[i] = dropoutTime[i];
            responder[i] = 0;
            source[i] = 3;
          } else {
            timeUnderObservation1[i] = time - arrivalTime[i];
            responder[i] = NA_LOGICAL;
            source[i] = 4;
          }
        }
      }

      // add raw data to output
      if (niter1[k] < maxNumberOfRawDatasetsPerStage) {
        for (int i=0; i<n; ++i) {
          iterationNumber1x[index1x] = iter+1;
          stageNumber1x[index1x] = k+1;
          analysisTime1x[index1x] = time;
          subjectId1x[index1x] = i+1;
          arrivalTime1x[index1x] = arrivalTime[i];
          stratum1x[index1x] = stratum[i];
          treatmentGroup1x[index1x] = treatmentGroup[i];
          survivalTime1x[index1x] = survivalTime[i];
          dropoutTime1x[index1x] = dropoutTime[i];
          ptfu1Timex[index1x] = ptfu1Time[i];
          timeUnderObservation1x[index1x] = timeUnderObservation1[i];
          responderx[index1x] = responder[i];
          sourcex[index1x] = source[i];

          ++index1x;
        }

        // update the number of stage k dataset to extract
        ++niter1[k];
      }


      accruals1 = sum(n1x);
      accruals2 = sum(n2x);
      totalAccruals = accruals1 + accruals2;

      // exclude subjects administratively censored
      LogicalVector subSorted = !is_na(responder);
      IntegerVector stratumSorted = stratum[subSorted];
      IntegerVector treatmentGroupSorted = treatmentGroup[subSorted];
      LogicalVector eventSorted = responder[subSorted];
      int nsubSorted = static_cast<int>(eventSorted.size());

      // obtain the Mantel-Haenszel statistic for stratified risk difference
      NumericVector n11(nstrata), n21(nstrata), n1s(nstrata), n2s(nstrata);
      for (int i=0; i<nsubSorted; ++i) {
        int h = stratumSorted[i] - 1;
        if (treatmentGroupSorted[i] == 1) {
          ++n1s[h];
          n11[h] += eventSorted[i];
        } else {
          ++n2s[h];
          n21[h] += eventSorted[i];
        }
      }
      NumericVector nss = n1s + n2s;

      double A = 0, B = 0, P = 0, Q = 0;
      for (int h=0; h<nstrata; ++h) {
        double dh = n11[h]/n1s[h] - n21[h]/n2s[h];
        double wh = n1s[h]*n2s[h]/nss[h];
        A += dh*wh;
        B += wh;
        P += (n1s[h]*n1s[h]*n21[h] - n2s[h]*n2s[h]*n11[h] +
          n1s[h]*n2s[h]*(n2s[h] - n1s[h])*0.5)/(nss[h]*nss[h]);
        Q += (n11[h]*(n2s[h]-n21[h]) + n21[h]*(n1s[h]-n11[h]))/(2*nss[h]);
      }

      double riskDiff = A/B;
      double seRiskDiff = sqrt(riskDiff*P + Q)/B;

      // add summary data to output
      iterationNumber1y[index1y] = iter+1;
      stageNumber1y[index1y] = k+1;
      analysisTime1y[index1y] = analysisTime1[k];
      accruals11y[index1y] = accruals1;
      accruals21y[index1y] = accruals2;
      totalAccruals1y[index1y] = totalAccruals;
      source1y[index1y] = sum(source == 1);
      source2y[index1y] = sum(source == 2);
      source3y[index1y] = sum(source == 3);
      n1y[index1y] = sum(n1s);
      n2y[index1y] = sum(n2s);
      ny[index1y] = sum(nss);
      y1y[index1y] = sum(n11);
      y2y[index1y] = sum(n21);
      yy[index1y] = sum(n11+n21);
      riskDiffy[index1y] = riskDiff;
      seRiskDiffy[index1y] = seRiskDiff;
      mhStatisticy[index1y] = (riskDiff - riskDiffH0)/seRiskDiff;
      ++index1y;
    } // end of stage for endpoint 1


    // find the analysis time for each stage for endpoint 2
    nevents = sum(event);
    totalt = stl_sort(totalTime[event]);
    nstages2 = kMax2;

    int j;
    for (j=0; j<kMax2; ++j) {
      if (plannedEvents[j] >= nevents) {
        nstages2 = j+1;
        break;
      }
    }

    if (j==kMax2) { // total number of events exceeds planned
      for (int k=0; k<nstages2; ++k) {
        analysisTime2[k] = totalt[plannedEvents[k]-1] + 1e-12;
      }
    } else {
      for (int k=0; k<nstages2; ++k) {
        if (k < nstages2-1) {
          analysisTime2[k] = totalt[plannedEvents[k]-1] + 1e-12;
        } else {
          analysisTime2[k] = totalt[nevents-1] + 1e-12;
        }
      }
    }

    // observed total number of events less than planned
    eventsNotAchieved = (nevents < plannedEvents[kMax2-1]);


    // construct the log-rank test statistic for endpoint 2 at each stage
    for (int k=0; k<nstages2; ++k) {
      time = analysisTime2[k];

      n1x.fill(0);  // number of subjects in each stratum by treatment
      n2x.fill(0);

      events1 = 0;
      events2 = 0;

      dropouts1 = 0;
      dropouts2 = 0;

      // censor at analysis time
      for (int i=0; i<n; ++i) {
        int h = stratum[i]-1;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation[i] = time - arrivalTime[i];
          event[i] = 0;
          dropoutEvent[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            ++n1x[h];
          } else {
            ++n2x[h];
          }

          // censored time for endpoint 2
          if (arrivalTime[i] + survivalTime[i] <= time &&
              survivalTime[i] <= dropoutTime[i]) {
            timeUnderObservation[i] = survivalTime[i];
            event[i] = 1;
            dropoutEvent[i] = 0;
          } else if (arrivalTime[i] + dropoutTime[i] <= time &&
            dropoutTime[i] <= survivalTime[i]) {
            timeUnderObservation[i] = dropoutTime[i];
            event[i] = 0;
            dropoutEvent[i] = 1;
          } else {
            timeUnderObservation[i] = time - arrivalTime[i];
            event[i] = 0;
            dropoutEvent[i] = 0;
          }

          if (treatmentGroup[i]==1 && event[i]) ++events1;
          if (treatmentGroup[i]==2 && event[i]) ++events2;
          if (treatmentGroup[i]==1 && dropoutEvent[i]) ++dropouts1;
          if (treatmentGroup[i]==2 && dropoutEvent[i]) ++dropouts2;
        }
      }


      // add raw data to output
      if (niter2[k] < maxNumberOfRawDatasetsPerStage) {
        for (int i=0; i<n; ++i) {
          iterationNumber2x[index2x] = iter+1;
          stageNumber2x[index2x] = k+1;
          analysisTime2x[index2x] = time;
          subjectId2x[index2x] = i+1;
          arrivalTime2x[index2x] = arrivalTime[i];
          stratum2x[index2x] = stratum[i];
          treatmentGroup2x[index2x] = treatmentGroup[i];
          survivalTime2x[index2x] = survivalTime[i];
          dropoutTime2x[index2x] = dropoutTime[i];
          timeUnderObservation2x[index2x] = timeUnderObservation[i];
          eventx[index2x] = event[i];
          dropoutEventx[index2x] = dropoutEvent[i];

          ++index2x;
        }

        // update the number of stage k dataset to extract
        ++niter2[k];
      }


      // number of accrued patients and total number of events
      accruals1 = sum(n1x);
      accruals2 = sum(n2x);
      totalAccruals = accruals1 + accruals2;
      totalEvents = events1 + events2;
      totalDropouts = dropouts1 + dropouts2;

      // reinitiate n1 and n2 for stage k
      n1 = clone(n1x);
      n2 = clone(n2x);

      // order the data by time under observation
      NumericVector timeUnderObservationSorted =
        stl_sort(timeUnderObservation);
      IntegerVector sortedIndex = match(timeUnderObservationSorted,
                                        timeUnderObservation);
      sortedIndex = sortedIndex - 1;
      IntegerVector stratumSorted = stratum[sortedIndex];
      IntegerVector treatmentGroupSorted = treatmentGroup[sortedIndex];
      LogicalVector eventSorted = event[sortedIndex];

      LogicalVector subSorted = (timeUnderObservationSorted > 0);
      eventSorted = eventSorted[subSorted];
      stratumSorted = stratumSorted[subSorted];
      treatmentGroupSorted = treatmentGroupSorted[subSorted];
      int nsubSorted = static_cast<int>(eventSorted.size());

      // calculate the stratified log-rank test
      uscore = 0;
      vscore = 0;
      for (int i=0; i<nsubSorted; ++i) {
        int h = stratumSorted[i] - 1;
        nt[h] = n1[h] + n2[h];

        n1a[h] = n1[h]*hazardRatioH0;
        nta[h] = n1a[h] + n2[h];

        if (eventSorted[i]) {
          uscore += (treatmentGroupSorted[i]==1) - n1a[h]/nta[h];
          vscore += n1a[h]*n2[h]/(nta[h]*nta[h]);
        }

        // reduce the risk set
        if (treatmentGroupSorted[i]==1) {
          --n1[h];
        } else {
          --n2[h];
        }
      }

      // add summary data to output
      iterationNumber2y[index2y] = iter+1;
      eventsNotAchievedy[index2y] = eventsNotAchieved;
      stageNumber2y[index2y] = k+1;
      analysisTime2y[index2y] = time;
      accruals1y[index2y] = accruals1;
      accruals2y[index2y] = accruals2;
      totalAccrualsy[index2y] = totalAccruals;
      events1y[index2y] = events1;
      events2y[index2y] = events2;
      totalEventsy[index2y] = totalEvents;
      dropouts1y[index2y] = dropouts1;
      dropouts2y[index2y] = dropouts2;
      totalDropoutsy[index2y] = totalDropouts;
      logRankStatisticy[index2y] = uscore/sqrt(vscore);
      ++index2y;

    } // end of stage for endpoint 2

  } // end of iteration



  // simulation summary data set for endpoint 1
  LogicalVector sub1 = !is_na(iterationNumber1y);
  iterationNumber1y = iterationNumber1y[sub1];
  stageNumber1y = stageNumber1y[sub1];
  analysisTime1y = analysisTime1y[sub1];
  accruals11y = accruals11y[sub1];
  accruals21y = accruals21y[sub1];
  totalAccruals1y = totalAccruals1y[sub1];
  source1y = source1y[sub1];
  source2y = source2y[sub1];
  source3y = source3y[sub1];
  n1y = n1y[sub1];
  n2y = n2y[sub1];
  ny = ny[sub1];
  y1y = y1y[sub1];
  y2y = y2y[sub1];
  yy = yy[sub1];
  riskDiffy = riskDiffy[sub1];
  seRiskDiffy = seRiskDiffy[sub1];
  mhStatisticy = mhStatisticy[sub1];

  DataFrame sumdataBIN = DataFrame::create(
    _["iterationNumber"] = iterationNumber1y,
    _["stageNumber"] = stageNumber1y,
    _["analysisTime"] = analysisTime1y,
    _["accruals1"] = accruals11y,
    _["accruals2"] = accruals21y,
    _["totalAccruals"] = totalAccruals1y,
    _["source1"] = source1y,
    _["source2"] = source2y,
    _["source3"] = source3y,
    _["n1"] = n1y,
    _["n2"] = n2y,
    _["n"] = ny,
    _["y1"] = y1y,
    _["y2"] = y2y,
    _["y"] = yy,
    _["riskDiff"] = riskDiffy,
    _["seRiskDiff"] = seRiskDiffy,
    _["mhStatistic"] = mhStatisticy);


  // simulation summary data set for endpoint 2
  LogicalVector sub2 = !is_na(iterationNumber2y);
  iterationNumber2y = iterationNumber2y[sub2];
  eventsNotAchievedy = eventsNotAchievedy[sub2];
  stageNumber2y = stageNumber2y[sub2];
  analysisTime2y = analysisTime2y[sub2];
  accruals1y = accruals1y[sub2];
  accruals2y = accruals2y[sub2];
  totalAccrualsy = totalAccrualsy[sub2];
  events1y = events1y[sub2];
  events2y = events2y[sub2];
  totalEventsy = totalEventsy[sub2];
  dropouts1y = dropouts1y[sub2];
  dropouts2y = dropouts2y[sub2];
  totalDropoutsy = totalDropoutsy[sub2];
  logRankStatisticy = logRankStatisticy[sub2];

  DataFrame sumdataTTE = DataFrame::create(
    _["iterationNumber"] = iterationNumber2y,
    _["eventsNotAchieved"] = eventsNotAchievedy,
    _["stageNumber"] = stageNumber2y,
    _["analysisTime"] = analysisTime2y,
    _["accruals1"] = accruals1y,
    _["accruals2"] = accruals2y,
    _["totalAccruals"] = totalAccrualsy,
    _["events1"] = events1y,
    _["events2"] = events2y,
    _["totalEvents"] = totalEventsy,
    _["dropouts1"] = dropouts1y,
    _["dropouts2"] = dropouts2y,
    _["totalDropouts"] = totalDropoutsy,
    _["logRankStatistic"] = logRankStatisticy);


  List result;

  if (maxNumberOfRawDatasetsPerStage > 0) {
    LogicalVector sub1 = !is_na(iterationNumber1x);
    iterationNumber1x = iterationNumber1x[sub1];
    stageNumber1x = stageNumber1x[sub1];
    analysisTime1x = analysisTime1x[sub1];
    subjectId1x = subjectId1x[sub1];
    arrivalTime1x = arrivalTime1x[sub1];
    stratum1x = stratum1x[sub1];
    treatmentGroup1x = treatmentGroup1x[sub1];
    survivalTime1x = survivalTime1x[sub1];
    dropoutTime1x = dropoutTime1x[sub1];
    ptfu1Timex = ptfu1Timex[sub1];
    timeUnderObservation1x = timeUnderObservation1x[sub1];
    responderx = responderx[sub1];
    sourcex = sourcex[sub1];

    DataFrame rawdataBIN = DataFrame::create(
      _["iterationNumber"] = iterationNumber1x,
      _["stageNumber"] = stageNumber1x,
      _["analysisTime"] = analysisTime1x,
      _["subjectId"] = subjectId1x,
      _["arrivalTime"] = arrivalTime1x,
      _["stratum"] = stratum1x,
      _["treatmentGroup"] = treatmentGroup1x,
      _["survivalTime"] = survivalTime1x,
      _["dropoutTime"] = dropoutTime1x,
      _["ptfu1Time"] = ptfu1Timex,
      _["timeUnderObservation"] = timeUnderObservation1x,
      _["responder"] = responderx,
      _["source"] = sourcex);


    LogicalVector sub2 = !is_na(iterationNumber2x);
    iterationNumber2x = iterationNumber2x[sub2];
    stageNumber2x = stageNumber2x[sub2];
    analysisTime2x = analysisTime2x[sub2];
    subjectId2x = subjectId2x[sub2];
    arrivalTime2x = arrivalTime2x[sub2];
    stratum2x = stratum2x[sub2];
    treatmentGroup2x = treatmentGroup2x[sub2];
    survivalTime2x = survivalTime2x[sub2];
    dropoutTime2x = dropoutTime2x[sub2];
    timeUnderObservation2x = timeUnderObservation2x[sub2];
    eventx = eventx[sub2];
    dropoutEventx = dropoutEventx[sub2];

    DataFrame rawdataTTE = DataFrame::create(
      _["iterationNumber"] = iterationNumber2x,
      _["stageNumber"] = stageNumber2x,
      _["analysisTime"] = analysisTime2x,
      _["subjectId"] = subjectId2x,
      _["arrivalTime"] = arrivalTime2x,
      _["stratum"] = stratum2x,
      _["treatmentGroup"] = treatmentGroup2x,
      _["survivalTime"] = survivalTime2x,
      _["dropoutTime"] = dropoutTime2x,
      _["timeUnderObservation"] = timeUnderObservation2x,
      _["event"] = eventx,
      _["dropoutEvent"] = dropoutEventx);

    result = List::create(_["sumdataBIN"] = sumdataBIN,
                          _["sumdataTTE"] = sumdataTTE,
                          _["rawdataBIN"] = rawdataBIN,
                          _["rawdataTTE"] = rawdataTTE);
  } else {
    result = List::create(_["sumdataBIN"] = sumdataBIN,
                          _["sumdataTTE"] = sumdataTTE);
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
//'   The seed from the environment will be used if left unspecified.
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
List lrsimsub(const int kMax = 1,
              const int kMaxitt = 1,
              const double hazardRatioH0itt = 1,
              const double hazardRatioH0pos = 1,
              const double hazardRatioH0neg = 1,
              const int allocation1 = 1,
              const int allocation2 = 1,
              const NumericVector& accrualTime = 0,
              const NumericVector& accrualIntensity = NA_REAL,
              const NumericVector& piecewiseSurvivalTime = 0,
              const NumericVector& stratumFraction = 1,
              const NumericVector& p_pos = NA_REAL,
              const NumericVector& lambda1itt = NA_REAL,
              const NumericVector& lambda2itt = NA_REAL,
              const NumericVector& lambda1pos = NA_REAL,
              const NumericVector& lambda2pos = NA_REAL,
              const NumericVector& gamma1itt = 0,
              const NumericVector& gamma2itt = 0,
              const NumericVector& gamma1pos = 0,
              const NumericVector& gamma2pos = 0,
              const int n = NA_INTEGER,
              const double followupTime = NA_REAL,
              const bool fixedFollowup = 0,
              const double rho1 = 0,
              const double rho2 = 0,
              const IntegerVector& plannedEvents = NA_INTEGER,
              const NumericVector& plannedTime = NA_REAL,
              const int maxNumberOfIterations = 1000,
              const int maxNumberOfRawDatasetsPerStage = 0,
              const int seed = NA_INTEGER) {

  // check input parameters
  int kMaxittx = kMaxitt;

  int nstrata = static_cast<int>(stratumFraction.size());
  int nints = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nints;
  int nints2 = nints == 1 ? 10 : nints + 10;
  int nsi2 = nstrata*nints2;

  NumericVector lambda1ittx(nsi), lambda2ittx(nsi);
  NumericVector lambda1posx(nsi), lambda2posx(nsi);
  NumericVector lambda1posy(nsi2), lambda2posy(nsi2);
  NumericVector lambda1negy(nsi2), lambda2negy(nsi2);
  NumericVector gamma1ittx(nsi), gamma2ittx(nsi);
  NumericVector gamma1posx(nsi), gamma2posx(nsi);
  NumericVector gamma1posy(nsi2), gamma2posy(nsi2);
  NumericVector gamma1negy(nsi2), gamma2negy(nsi2);

  NumericVector tau1pos(nsi2), tau2pos(nsi2);
  NumericVector tau1neg(nsi2), tau2neg(nsi2);

  NumericVector p_posx(nstrata);

  bool useEvents, eventsNotAchieved;


  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }

  if (kMaxitt < 0) {
    kMaxittx = kMax;
  }

  if (kMaxittx > kMax) {
    stop("kMaxitt must be less than or equal to kMax");
  }

  // whether to plan the analyses based on events or calendar time
  if (is_false(any(is_na(plannedEvents)))) {
    useEvents = 1;

    if (plannedEvents[0] <= 0) {
      stop("Elements of plannedEvents must be positive");
    }

    if (plannedEvents.size() != kMax) {
      stop("Invalid length for plannedEvents");
    }

    if (kMaxittx > 1) {
      IntegerVector plannedEvents1 = plannedEvents[Range(0,kMaxittx-1)];
      if (is_true(any(diff(plannedEvents1) <= 0))) {
        stop("plannedEvents for ITT must be increasing");
      }
    }

    if (kMax - kMaxittx > 1) {
      IntegerVector plannedEvents2 = plannedEvents[Range(kMaxittx, kMax-1)];
      if (is_true(any(diff(plannedEvents2) <= 0))) {
        stop("plannedEvents for biomarker+ sub population must be increasing");
      }
    }
  } else if (is_false(any(is_na(plannedTime)))) {
    useEvents = 0;
    if (plannedTime[0] <= 0) {
      stop("Elements of plannedTime must be positive");
    }

    if (plannedTime.size() != kMax) {
      stop("Invalid length for plannedTime");
    }

    if (kMax > 1 && is_true(any(diff(plannedTime) <= 0))) {
      stop("Elements of plannedTime must be increasing");
    }
  } else {
    stop("Either plannedEvents or plannedTime must be given");
  }


  if (hazardRatioH0itt <= 0) {
    stop("hazardRatioH0itt must be positive");
  }

  if (hazardRatioH0pos <= 0) {
    stop("hazardRatioH0pos must be positive");
  }

  if (hazardRatioH0neg <= 0) {
    stop("hazardRatioH0neg must be positive");
  }


  if (allocation1 < 1) {
    stop("allocation1 must be a positive integer");
  }

  if (allocation2 < 1) {
    stop("allocation2 must be a positive integer");
  }


  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }


  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }

  if (nints > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }


  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }


  if (is_true(any(is_na(p_pos)))) {
    stop("p_pos must be provided");
  }

  if (p_pos.size() == 1) {
    p_posx = rep(p_pos, nstrata);
  } else if (p_pos.size() == nstrata) {
    p_posx = p_pos;
  } else {
    stop("Invalid length for p_pos");
  }

  if (is_true(any((p_posx <= 0) | (p_posx >= 1)))) {
    stop("p_pos must lie between 0 and 1");
  }


  if (is_true(any(is_na(lambda1itt)))) {
    stop("lambda1itt must be provided");
  }

  if (is_true(any(is_na(lambda2itt)))) {
    stop("lambda2itt must be provided");
  }

  if (is_true(any(is_na(lambda1pos)))) {
    stop("lambda1pos must be provided");
  }

  if (is_true(any(is_na(lambda2pos)))) {
    stop("lambda2pos must be provided");
  }

  if (is_true(any(lambda1itt < 0))) {
    stop("lambda1itt must be non-negative");
  }

  if (is_true(any(lambda2itt < 0))) {
    stop("lambda2itt must be non-negative");
  }

  if (is_true(any(lambda1pos < 0))) {
    stop("lambda1pos must be non-negative");
  }

  if (is_true(any(lambda2pos < 0))) {
    stop("lambda2pos must be non-negative");
  }


  if (is_true(any(gamma1itt < 0))) {
    stop("gamma1itt must be non-negative");
  }

  if (is_true(any(gamma2itt < 0))) {
    stop("gamma2itt must be non-negative");
  }

  if (is_true(any(gamma1pos < 0))) {
    stop("gamma1pos must be non-negative");
  }

  if (is_true(any(gamma2pos < 0))) {
    stop("gamma2pos must be non-negative");
  }


  if (lambda1itt.size() == 1) {
    lambda1ittx = rep(lambda1itt, nsi);
  } else if (lambda1itt.size() == nints) {
    lambda1ittx = rep(lambda1itt, nstrata);
  } else if (lambda1itt.size() == nsi) {
    lambda1ittx = lambda1itt;
  } else {
    stop("Invalid length for lambda1itt");
  }

  if (lambda2itt.size() == 1) {
    lambda2ittx = rep(lambda2itt, nsi);
  } else if (lambda2itt.size() == nints) {
    lambda2ittx = rep(lambda2itt, nstrata);
  } else if (lambda2itt.size() == nsi) {
    lambda2ittx = lambda2itt;
  } else {
    stop("Invalid length for lambda2itt");
  }


  if (lambda1pos.size() == 1) {
    lambda1posx = rep(lambda1pos, nsi);
  } else if (lambda1pos.size() == nints) {
    lambda1posx = rep(lambda1pos, nstrata);
  } else if (lambda1pos.size() == nsi) {
    lambda1posx = lambda1pos;
  } else {
    stop("Invalid length for lambda1pos");
  }

  if (lambda2pos.size() == 1) {
    lambda2posx = rep(lambda2pos, nsi);
  } else if (lambda2pos.size() == nints) {
    lambda2posx = rep(lambda2pos, nstrata);
  } else if (lambda2pos.size() == nsi) {
    lambda2posx = lambda2pos;
  } else {
    stop("Invalid length for lambda2pos");
  }


  if (gamma1itt.size() == 1) {
    gamma1ittx = rep(gamma1itt, nsi);
  } else if (gamma1itt.size() == nints) {
    gamma1ittx = rep(gamma1itt, nstrata);
  } else if (gamma1itt.size() == nsi) {
    gamma1ittx = gamma1itt;
  } else {
    stop("Invalid length for gamma1itt");
  }

  if (gamma2itt.size() == 1) {
    gamma2ittx = rep(gamma2itt, nsi);
  } else if (gamma2itt.size() == nints) {
    gamma2ittx = rep(gamma2itt, nstrata);
  } else if (gamma2itt.size() == nsi) {
    gamma2ittx = gamma2itt;
  } else {
    stop("Invalid length for gamma2itt");
  }


  if (gamma1pos.size() == 1) {
    gamma1posx = rep(gamma1pos, nsi);
  } else if (gamma1pos.size() == nints) {
    gamma1posx = rep(gamma1pos, nstrata);
  } else if (gamma1pos.size() == nsi) {
    gamma1posx = gamma1pos;
  } else {
    stop("Invalid length for gamma1pos");
  }

  if (gamma2pos.size() == 1) {
    gamma2posx = rep(gamma2pos, nsi);
  } else if (gamma2pos.size() == nints) {
    gamma2posx = rep(gamma2pos, nstrata);
  } else if (gamma2pos.size() == nsi) {
    gamma2posx = gamma2pos;
  } else {
    stop("Invalid length for gamma2pos");
  }


  for (int j=0; j<nstrata; ++j) {
    Range jj = Range(j*nints, (j+1)*nints-1);

    double ppos = p_posx[j];

    NumericVector lam1itt = lambda1ittx[jj];
    NumericVector lam1pos = lambda1posx[jj];

    List a1 = hazard_subcpp(piecewiseSurvivalTime,
                            lam1itt, lam1pos, ppos);

    NumericVector time1pos = a1["piecewiseSurvivalTime"];
    NumericVector lam1posy = a1["hazard_pos"];
    NumericVector lam1negy = a1["hazard_neg"];

    NumericVector lam2itt = lambda2ittx[jj];
    NumericVector lam2pos = lambda2posx[jj];

    List a2 = hazard_subcpp(piecewiseSurvivalTime,
                            lam2itt, lam2pos, ppos);

    NumericVector time2pos = a2["piecewiseSurvivalTime"];
    NumericVector lam2posy = a2["hazard_pos"];
    NumericVector lam2negy = a2["hazard_neg"];

    NumericVector gam1itt = gamma1ittx[jj];
    NumericVector gam1pos = gamma1posx[jj];

    List b1 = hazard_subcpp(piecewiseSurvivalTime,
                            gam1itt, gam1pos, ppos);

    NumericVector time1neg = b1["piecewiseSurvivalTime"];
    NumericVector gam1posy = b1["hazard_pos"];
    NumericVector gam1negy = b1["hazard_neg"];

    NumericVector gam2itt = gamma2ittx[jj];
    NumericVector gam2pos = gamma2posx[jj];

    List b2 = hazard_subcpp(piecewiseSurvivalTime,
                            gam2itt, gam2pos, ppos);

    NumericVector time2neg = b2["piecewiseSurvivalTime"];
    NumericVector gam2posy = b2["hazard_pos"];
    NumericVector gam2negy = b2["hazard_neg"];

    for (int k=0; k<nints2; ++k) {
      tau1pos[j*nints2 + k] = time1pos[k];
      tau2pos[j*nints2 + k] = time2pos[k];
      tau1neg[j*nints2 + k] = time1neg[k];
      tau2neg[j*nints2 + k] = time2neg[k];
      lambda1posy[j*nints2 + k] = lam1posy[k];
      lambda2posy[j*nints2 + k] = lam2posy[k];
      lambda1negy[j*nints2 + k] = lam1negy[k];
      lambda2negy[j*nints2 + k] = lam2negy[k];
      gamma1posy[j*nints2 + k] = gam1posy[k];
      gamma2posy[j*nints2 + k] = gam2posy[k];
      gamma1negy[j*nints2 + k] = gam1negy[k];
      gamma2negy[j*nints2 + k] = gam2negy[k];
    }
  }


  if (n == NA_INTEGER) {
    stop("n must be provided");
  }

  if (n <= 0) {
    stop("n must be positive");
  }

  if (fixedFollowup) {
    if (R_isnancpp(followupTime)) {
      stop("followupTime must be provided for fixed follow-up");
    }

    if (followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }
  }

  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }

  if (rho2 < 0) {
    stop("rho2 must be non-negative");
  }

  if (maxNumberOfIterations < 1) {
    stop("maxNumberOfIterations must be a positive integer");
  }

  if (maxNumberOfRawDatasetsPerStage < 0) {
    stop("maxNumberOfRawDatasetsPerStage must be a non-negative integer");
  }


  // declare variables
  int nevents, neventspos, nstages;
  int accruals1, accruals2, totalAccruals;
  int events1, events2, totalEvents;
  int dropouts1, dropouts2, totalDropouts;
  int index1=0, index2=0;

  double enrollt, u, time, uscore, vscore, hazardRatioH0;


  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n);

  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n),
  timeUnderObservation(n), totalTime(n), totalt(n), totaltpos(n);

  LogicalVector biomarker(n), event(n), dropoutEvent(n), eventpos(n);

  // stratum information
  IntegerVector b1(nstrata), b2(nstrata);
  IntegerVector n1(nstrata), n2(nstrata);

  // hazardRatioH0 adjusted n1 and nt for calculating the log-rank statistic
  NumericVector nt(nstrata), n1a(nstrata), nta(nstrata);

  NumericVector km(nstrata), w(nstrata);
  NumericVector cumStratumFraction = cumsum(stratumFraction);


  // within-stratum hazard rates
  NumericVector time1(nints2), time2(nints2);
  NumericVector lam1(nints2), lam2(nints2);
  NumericVector gam1(nints2), gam2(nints2);

  // stage-wise information
  IntegerVector niter(kMax);
  NumericVector analysisTime(kMax);


  // cache for the patient-level raw data to extract
  int nrow1 = n*kMax*maxNumberOfRawDatasetsPerStage;

  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stageNumberx(nrow1);
  NumericVector analysisTimex(nrow1);
  IntegerVector subjectIdx(nrow1);
  NumericVector arrivalTimex(nrow1);
  IntegerVector stratumx(nrow1);
  LogicalVector biomarkerx(nrow1);
  IntegerVector treatmentGroupx(nrow1);
  NumericVector survivalTimex(nrow1);
  NumericVector dropoutTimex(nrow1);
  NumericVector timeUnderObservationx(nrow1);
  LogicalVector eventx(nrow1);
  LogicalVector dropoutEventx(nrow1);


  // cache for the simulation-level summary data to extract
  int nrow2 = kMax*maxNumberOfIterations*3;

  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  LogicalVector eventsNotAchievedy(nrow2);
  IntegerVector stageNumbery(nrow2);
  NumericVector analysisTimey(nrow2);
  StringVector populationy(nrow2);
  IntegerVector accruals1y(nrow2);
  IntegerVector accruals2y(nrow2);
  IntegerVector totalAccrualsy(nrow2);
  IntegerVector events1y(nrow2);
  IntegerVector events2y(nrow2);
  IntegerVector totalEventsy(nrow2);
  IntegerVector dropouts1y(nrow2);
  IntegerVector dropouts2y(nrow2);
  IntegerVector totalDropoutsy(nrow2);
  NumericVector logRankStatisticy(nrow2);


  // set up random seed
  if (seed != NA_INTEGER) {
    set_seed(seed);
  }


  // simulation
  for (int iter=0; iter<maxNumberOfIterations; ++iter) {

    b1.fill(allocation1);
    b2.fill(allocation2);

    enrollt = 0;
    for (int i=0; i<n; ++i) {

      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);
      arrivalTime[i] = enrollt;

      // generate stratum information
      u = R::runif(0,1);
      int j;
      for (j=0; j<nstrata; ++j) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }

      // generate biomarker status
      u = R::runif(0,1);
      if (u <= p_pos[j]) {
        biomarker[i] = 1; // positive
      } else {
        biomarker[i] = 0; // negative
      }

      // generate treatment group
      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+0.0)) {
        treatmentGroup[i] = 1;
        --b1[j];
      } else {
        treatmentGroup[i] = 2;
        --b2[j];
      }

      // start a new block after depleting the current block
      if (b1[j]+b2[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
      }

      // stratum-specific hazard rates for event and dropout
      Range jj2 = Range(j*nints2, (j+1)*nints2-1);

      if (biomarker[i]) { // biomarker positive
        time1 = tau1pos[jj2];
        time2 = tau2pos[jj2];
        lam1 = lambda1posy[jj2];
        lam2 = lambda2posy[jj2];
        gam1 = gamma1posy[jj2];
        gam2 = gamma2posy[jj2];
      } else { // biomarker negative
        time1 = tau1neg[jj2];
        time2 = tau2neg[jj2];
        lam1 = lambda1negy[jj2];
        lam2 = lambda2negy[jj2];
        gam1 = gamma1negy[jj2];
        gam2 = gamma2negy[jj2];
      }


      // generate survival times
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        survivalTime[i] = qtpwexpcpp1(u, time1, lam1, 0, 1, 0);
      } else {
        survivalTime[i] = qtpwexpcpp1(u, time2, lam2, 0, 1, 0);
      }

      // generate dropout times
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime[i] = qtpwexpcpp1(u, time1, gam1, 0, 1, 0);
      } else {
        dropoutTime[i] = qtpwexpcpp1(u, time2, gam2, 0, 1, 0);
      }

      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime[i] <= dropoutTime[i] &&
            survivalTime[i] <= followupTime) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else if (dropoutTime[i] <= survivalTime[i] &&
          dropoutTime[i] <= followupTime) {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        } else {
          timeUnderObservation[i] = followupTime;
          event[i] = 0;
          dropoutEvent[i] = 0;
        }
      } else {
        if (survivalTime[i] <= dropoutTime[i]) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        }
      }

      totalTime[i] = arrivalTime[i] + timeUnderObservation[i];
    }

    // find the analysis time for each stage
    if (useEvents) {
      eventpos = event & biomarker;
      nevents = sum(event);
      neventspos = sum(eventpos);
      totalt = stl_sort(totalTime[event]);
      totaltpos = stl_sort(totalTime[eventpos]);

      int j1 = kMaxittx, j2 = kMax - kMaxittx;

      // ITT looks
      if (kMaxittx > 0) {
        for (j1=0; j1<kMaxittx; ++j1) {
          if (plannedEvents[j1] >= nevents) {
            break;
          }
        }

        if (j1==kMaxittx) { // total number of events exceeds planned
          for (int k=0; k<kMaxittx; ++k) { // analyses occur at planned events
            analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
          }
        } else {
          for (int k=0; k<=j1; ++k) {
            if (k < j1) {
              analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
            } else { // the last look may have events <= planned
              analysisTime[k] = totalt[nevents-1] + 1e-12;
            }
          }
        }
      }

      // biomarker positive looks
      NumericVector analysisTime2(kMax - kMaxittx);

      if (kMax > kMaxittx) {
        for (j2=0; j2<kMax-kMaxittx; ++j2) {
          if (plannedEvents[kMaxittx+j2] >= neventspos) {
            break;
          }
        }

        if (j2==kMax-kMaxittx) { // total number of events exceeds planned
          for (int k=0; k<kMax-kMaxittx; ++k) {
            analysisTime2[k] = totaltpos[plannedEvents[kMaxittx+k]-1] + 1e-12;
          }
        } else {
          for (int k=0; k<=j2; ++k) {
            if (k < j2) {
              analysisTime2[k] = totaltpos[plannedEvents[kMaxittx+k]-1]+1e-12;
            } else {
              analysisTime2[k] = totaltpos[neventspos-1] + 1e-12;
            }
          }
        }
      }

      // determine the number of looks and timing of the looks
      if (kMaxittx == 0) { // all looks based on biomarker positive events
        if (j2 == kMax - kMaxittx) {
          nstages = kMax - kMaxittx;
        } else {
          nstages = j2 + 1;
        }

        for (int k=0; k<nstages; ++k) {
          analysisTime[k] = analysisTime2[k];
        }
      } else if (kMax == kMaxittx) { // all looks based on ITT events
        if (j1 == kMaxittx) {
          nstages = kMaxittx;
        } else {
          nstages = j1 + 1;
        }
      } else {
        if (analysisTime2[kMax-kMaxittx-1] > analysisTime[kMaxittx-1]) {
          // only biomarker positive looks that occur after the last ITT
          // look contribute
          int l = static_cast<int>(which_max(analysisTime2 >
                                               analysisTime[kMaxittx-1]));
          nstages = kMax-l;
          for (int k=kMaxittx; k<kMax-l; ++k) {
            analysisTime[k] = analysisTime2[k-kMaxittx+l];
          }
        } else { // no biomarker positive looks after the last ITT look
          if (j1 == kMaxittx) {
            nstages = kMaxittx;
          } else {
            nstages = j1 + 1;
          }
        }
      }

      // whether the target PFS and OS events are achieved
      if (kMaxittx > 0 && nevents < plannedEvents[kMaxittx-1]) {
        eventsNotAchieved = 1;
      } else if (kMaxittx < kMax && neventspos < plannedEvents[kMax-1]) {
        eventsNotAchieved = 1;
      } else {
        eventsNotAchieved = 0;
      }
    } else { // looks based on calendar time
      nstages = kMax;
      analysisTime = clone(plannedTime);
      eventsNotAchieved = 0;
    }


    // construct the log-rank test statistic at each stage
    for (int k=0; k<nstages; ++k) {
      time = analysisTime[k];

      for (int pop=1; pop<=3; ++pop) {
        // subset subjects for the population of interest
        LogicalVector sub(n);
        if (pop == 1) { // ITT
          sub = LogicalVector(n, 1);
          hazardRatioH0 = hazardRatioH0itt;
        } else if (pop == 2) { // biomarker positive
          sub = biomarker;
          hazardRatioH0 = hazardRatioH0pos;
        } else { // biomarker negative
          sub = !biomarker;
          hazardRatioH0 = hazardRatioH0neg;
        }

        int nsub = sum(sub);
        IntegerVector stratumSub = stratum[sub];
        IntegerVector treatmentGroupSub = treatmentGroup[sub];
        NumericVector arrivalTimeSub = arrivalTime[sub];
        NumericVector survivalTimeSub = survivalTime[sub];
        NumericVector dropoutTimeSub = dropoutTime[sub];
        NumericVector timeUnderObservationSub = timeUnderObservation[sub];
        LogicalVector eventSub = event[sub];
        LogicalVector dropoutEventSub = dropoutEvent[sub];

        n1.fill(0);
        n2.fill(0);
        events1 = 0;
        events2 = 0;
        dropouts1 = 0;
        dropouts2 = 0;

        // censor at analysis time
        for (int i=0; i<nsub; ++i) {
          int h = stratumSub[i]-1;
          if (arrivalTimeSub[i] > time) { // patients not yet enrolled
            timeUnderObservationSub[i] = time - arrivalTimeSub[i];
            eventSub[i] = 0;
            dropoutEventSub[i] = 0;
          } else {
            if (treatmentGroupSub[i]==1) {
              ++n1[h];
            } else {
              ++n2[h];
            }

            if (fixedFollowup) {
              // the first three cases correspond to arrivalTime[i] +
              // min(survivalTime[i], dropoutTime[i], followupTime) <= time
              if (arrivalTimeSub[i] + survivalTimeSub[i] <= time &&
                  survivalTimeSub[i] <= dropoutTimeSub[i] &&
                  survivalTimeSub[i] <= followupTime) {
                timeUnderObservationSub[i] = survivalTimeSub[i];
                eventSub[i] = 1;
                dropoutEventSub[i] = 0;
              } else if (arrivalTimeSub[i] + dropoutTimeSub[i] <= time &&
                dropoutTimeSub[i] <= survivalTimeSub[i] &&
                dropoutTimeSub[i] <= followupTime) {
                timeUnderObservationSub[i] = dropoutTimeSub[i];
                eventSub[i] = 0;
                dropoutEventSub[i] = 1;
              } else if (arrivalTimeSub[i] + followupTime <= time &&
                followupTime <= survivalTimeSub[i] &&
                followupTime <= dropoutTimeSub[i]) {
                timeUnderObservationSub[i] = followupTime;
                eventSub[i] = 0;
                dropoutEventSub[i] = 0;
              } else {
                timeUnderObservationSub[i] = time - arrivalTimeSub[i];
                eventSub[i] = 0;
                dropoutEventSub[i] = 0;
              }
            } else {
              if (arrivalTimeSub[i] + survivalTimeSub[i] <= time &&
                  survivalTimeSub[i] <= dropoutTimeSub[i]) {
                timeUnderObservationSub[i] = survivalTimeSub[i];
                eventSub[i] = 1;
                dropoutEventSub[i] = 0;
              } else if (arrivalTimeSub[i] + dropoutTimeSub[i] <= time &&
                dropoutTimeSub[i] <= survivalTimeSub[i]) {
                timeUnderObservationSub[i] = dropoutTimeSub[i];
                eventSub[i] = 0;
                dropoutEventSub[i] = 1;
              } else {
                timeUnderObservationSub[i] = time - arrivalTimeSub[i];
                eventSub[i] = 0;
                dropoutEventSub[i] = 0;
              }
            }

            if (treatmentGroupSub[i]==1 && eventSub[i]) ++events1;
            if (treatmentGroupSub[i]==2 && eventSub[i]) ++events2;
            if (treatmentGroupSub[i]==1 && dropoutEventSub[i]) ++dropouts1;
            if (treatmentGroupSub[i]==2 && dropoutEventSub[i]) ++dropouts2;
          }
        }


        // add raw data to output
        if (pop == 1) {
          if (niter[k] < maxNumberOfRawDatasetsPerStage) {
            for (int i=0; i<n; ++i) {
              iterationNumberx[index1] = iter+1;
              stageNumberx[index1] = k+1;
              analysisTimex[index1] = time;
              subjectIdx[index1] = i+1;
              arrivalTimex[index1] = arrivalTime[i];
              stratumx[index1] = stratum[i];
              biomarkerx[index1] = biomarker[i];
              treatmentGroupx[index1] = treatmentGroup[i];
              survivalTimex[index1] = survivalTime[i];
              dropoutTimex[index1] = dropoutTime[i];
              timeUnderObservationx[index1] = timeUnderObservationSub[i];
              eventx[index1] = eventSub[i];
              dropoutEventx[index1] = dropoutEventSub[i];

              ++index1;
            }

            // update the number of stage k dataset to extract
            ++niter[k];
          }
        }

        // number of accrued patients and total number of events
        accruals1 = sum(n1);
        accruals2 = sum(n2);
        totalAccruals = accruals1 + accruals2;

        totalEvents = events1 + events2;
        totalDropouts = dropouts1 + dropouts2;

        // order the data by time under observation
        NumericVector timeUnderObservationSorted =
          stl_sort(timeUnderObservationSub);
        IntegerVector sortedIndex = match(timeUnderObservationSorted,
                                          timeUnderObservationSub);
        sortedIndex = sortedIndex - 1;
        IntegerVector stratumSorted = stratumSub[sortedIndex];
        IntegerVector treatmentGroupSorted = treatmentGroupSub[sortedIndex];
        LogicalVector eventSorted = eventSub[sortedIndex];

        LogicalVector subSorted = (timeUnderObservationSorted > 0);
        stratumSorted = stratumSorted[subSorted];
        treatmentGroupSorted = treatmentGroupSorted[subSorted];
        eventSorted = eventSorted[subSorted];
        int nsubSorted = static_cast<int>(eventSorted.size());

        // calculate the stratified log-rank test
        uscore = 0;
        vscore = 0;
        km.fill(1);
        for (int i=0; i<nsubSorted; ++i) {
          int h = stratumSorted[i] - 1;
          n1a[h] = n1[h]*hazardRatioH0;
          nt[h] = n1[h] + n2[h];
          nta[h] = n1a[h] + n2[h];

          if (eventSorted[i]) {
            w[h] = pow(km[h], rho1)*pow(1-km[h], rho2);
            uscore += w[h]*((treatmentGroupSorted[i]==1) - n1a[h]/nta[h]);
            vscore += w[h]*w[h]*n1a[h]*n2[h]/(nta[h]*nta[h]);
            km[h] *= (1-1/nt[h]); // update km estimate
          }

          // reduce the risk set
          if (treatmentGroupSorted[i]==1) {
            --n1[h];
          } else {
            --n2[h];
          }
        }

        // add summary data to output
        iterationNumbery[index2] = iter+1;
        eventsNotAchievedy[index2] = eventsNotAchieved;
        stageNumbery[index2] = k+1;
        analysisTimey[index2] = time;
        populationy[index2] = pop == 1 ? "ITT" :
          (pop == 2 ? "Biomarker Positive" : "Biomarker Negative");
        accruals1y[index2] = accruals1;
        accruals2y[index2] = accruals2;
        totalAccrualsy[index2] = totalAccruals;
        events1y[index2] = events1;
        events2y[index2] = events2;
        totalEventsy[index2] = totalEvents;
        dropouts1y[index2] = dropouts1;
        dropouts2y[index2] = dropouts2;
        totalDropoutsy[index2] = totalDropouts;
        logRankStatisticy[index2] = uscore/sqrt(vscore);
        ++index2;
      } // end of population of interest
    } // end of stage
  } // end of iteration



  // only keep nonmissing records
  LogicalVector sub2 = !is_na(iterationNumbery);
  iterationNumbery = iterationNumbery[sub2];
  eventsNotAchievedy = eventsNotAchievedy[sub2];
  stageNumbery = stageNumbery[sub2];
  analysisTimey = analysisTimey[sub2];
  populationy = populationy[sub2];
  accruals1y = accruals1y[sub2];
  accruals2y = accruals2y[sub2];
  totalAccrualsy = totalAccrualsy[sub2];
  events1y = events1y[sub2];
  events2y = events2y[sub2];
  totalEventsy = totalEventsy[sub2];
  dropouts1y = dropouts1y[sub2];
  dropouts2y = dropouts2y[sub2];
  totalDropoutsy = totalDropoutsy[sub2];
  logRankStatisticy = logRankStatisticy[sub2];

  DataFrame sumdata = DataFrame::create(
    _["iterationNumber"] = iterationNumbery,
    _["eventsNotAchieved"] = eventsNotAchievedy,
    _["stageNumber"] = stageNumbery,
    _["analysisTime"] = analysisTimey,
    _["population"] = populationy,
    _["accruals1"] = accruals1y,
    _["accruals2"] = accruals2y,
    _["totalAccruals"] = totalAccrualsy,
    _["events1"] = events1y,
    _["events2"] = events2y,
    _["totalEvents"] = totalEventsy,
    _["dropouts1"] = dropouts1y,
    _["dropouts2"] = dropouts2y,
    _["totalDropouts"] = totalDropoutsy,
    _["logRankStatistic"] = logRankStatisticy);


  List result;

  if (maxNumberOfRawDatasetsPerStage > 0) {
    LogicalVector sub1 = !is_na(iterationNumberx);
    iterationNumberx = iterationNumberx[sub1];
    stageNumberx = stageNumberx[sub1];
    analysisTimex = analysisTimex[sub1];
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    biomarkerx = biomarkerx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    survivalTimex = survivalTimex[sub1];
    dropoutTimex = dropoutTimex[sub1];
    timeUnderObservationx = timeUnderObservationx[sub1];
    eventx = eventx[sub1];
    dropoutEventx = dropoutEventx[sub1];

    DataFrame rawdata = DataFrame::create(
      _["iterationNumber"] = iterationNumberx,
      _["stageNumber"] = stageNumberx,
      _["analysisTime"] = analysisTimex,
      _["subjectId"] = subjectIdx,
      _["arrivalTime"] = arrivalTimex,
      _["stratum"] = stratumx,
      _["biomarker"] = biomarkerx,
      _["treatmentGroup"] = treatmentGroupx,
      _["survivalTime"] = survivalTimex,
      _["dropoutTime"] = dropoutTimex,
      _["timeUnderObservation"] = timeUnderObservationx,
      _["event"] = eventx,
      _["dropoutEvent"] = dropoutEventx);

    result = List::create(_["sumdata"] = sumdata,
                          _["rawdata"] = rawdata);
  } else {
    result = List::create(_["sumdata"] = sumdata);
  }

  return result;
}
