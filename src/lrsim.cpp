#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;


//' @title Log-rank test simulation
//' @description Performs simulation for two-arm group sequential
//' trials based on weighted log-rank test.
//'
//' @inheritParams param_kMax
//' @param informationTime Information time in terms of variance of
//'   weighted log-rank test score statistic under the null hypothesis.
//'   Same as \code{informationRates} in terms of number of events for
//'   the conventional log-rank test. Use \code{caltime} and \code{lrstat}
//'   to derive the information time for weighted log-rank tests.
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
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @param plannedEvents The planned cumulative total number of events at each
//'   stage.
//' @param plannedTime The calendar times for the analyses. To use calendar  
//'   time to plan the analyses, \code{plannedEvents} should be missing.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per stage
//'   to extract. Defaults to 1.
//' @param seed The seed to reproduce the simulation results.
//'   The computer clock will be used if left unspecified,
//'
//' @return An S3 class \code{lrsim} object with 3 components:
//'
//' * \code{overview}: A list containing the following information:
//' 
//'   - \code{rejectPerStage}: The efficacy stopping probability by stage.
//'   
//'   - \code{futilityPerStage}: The futility stopping probability by stage.
//'   
//'   - \code{cumulativeRejection}: Cumulative efficacy stopping 
//'   probability by stage.
//'   
//'   - \code{cumulativeFutility}: The cumulative futility stopping 
//'   probability by stage.
//'   
//'   - \code{numberOfEvents}: The average number of events by stage.
//'   
//'   - \code{numberOfDropouts}: The average number of dropouts by stage.
//'   
//'   - \code{numberOfSubjects}: The average number of subjects by stage.
//'   
//'   - \code{analysisTime}: The average analysis time by stage.
//'   
//'   - \code{overallReject}: The overall rejection probability. 
//'   
//'   - \code{expectedNumberOfEvents}: The expected number of events for 
//'   the overall study.
//'   
//'   - \code{expectedNumberOfDropouts}: The expected number of dropouts for
//'   the overall study.
//'   
//'   - \code{expectedNumberOfSubjects}: The expected number of subjects for 
//'   the overall study. 
//'   
//'   - \code{expectedStudyDuration}: The expected study duration.
//'   
//'   - \code{hazardRatioH0}: Hazard ratio under the null hypothesis for 
//'   the active treatment versus control.
//'   
//'   - \code{useEvents}: whether the analyses are planned 
//'   based on the number of events or calendar time. 
//'   
//'   - \code{accrualDuration}: Duration of the enrollment period.
//'   
//'   - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'   
//'   - \code{rho1}: The first parameter of the Fleming-Harrington family 
//'   of weighted log-rank test. Defaults to 0 for conventional log-rank test.
//'   
//'   - \code{rho2}: The second parameter of the Fleming-Harrington family 
//'   of weighted log-rank test. Defaults to 0 for conventional log-rank test.
//'   
//'   - \code{kMax}: The maximum number of stages.
//'
//' * \code{sumdata}: A data frame of summary data by iteration and stage: 
//' 
//'   - \code{iterationNumber}: The iteration number.
//'   
//'   - \code{stopStage}: The stage at which the trial stops.
//'   
//'   - \code{eventsNotAchieved}: Whether the target number of events 
//'   is not achieved for the iteration.
//'   
//'   - \code{stageNumber}: The stage number, covering all stages even if 
//'   the trial stops at an interim look.
//'   
//'   - \code{analysisTime}: The time for the stage since trial start.
//'   
//'   - \code{accruals1}: The number of subjects enrolled at the stage for
//'   the treatment group.
//'   
//'   - \code{accruals2}: The number of subjects enrolled at the stage for
//'   the control group.
//'   
//'   - \code{totalAccruals}: The total number of subjects enrolled at 
//'   the stage.
//'   
//'   - \code{events1}: The number of events at the stage for 
//'   the treatment group.
//'   
//'   - \code{events2}: The number of events at the stage for 
//'   the control group.
//'   
//'   - \code{totalEvents}: The total number of events at the stage.
//'   
//'   - \code{dropouts1}: The number of dropouts at the stage for 
//'   the treatment group. 
//'   
//'   - \code{dropouts2}: The number of dropouts at the stage for 
//'   the control group.
//'   
//'   - \code{totalDropouts}: The total number of dropouts at the stage.
//'   
//'   - \code{uscore}: The numerator of the log-rank test statistic. 
//'   
//'   - \code{vscore}: The variance of the log-rank test statistic.
//'   
//'   - \code{logRankStatistic}: The log-rank test Z-statistic.
//'   
//'   - \code{rejectPerStage}: Whether to reject the null hypothesis 
//'   at the stage.
//'   
//'   - \code{futilityPerStage}: Whether to stop the trial for futility  
//'   at the stage.
//'
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//' positive integer): A data frame for subject-level data for selected
//' replications, containing the following variables:
//' 
//'   - \code{iterationNumber}: The iteration number.
//'   
//'   - \code{stopStage}: The stage at which the trial stops.
//'   
//'   - \code{analysisTime}: The time for the stage since trial start.
//'   
//'   - \code{subjectId}: The subject ID.
//'   
//'   - \code{arrivalTime}: The enrollment time for the subject.
//'   
//'   - \code{stratum}: The stratum for the subject.
//'   
//'   - \code{treatmentGroup}: The treatment group (1 or 2) for the subject.
//'   
//'   - \code{survivalTime}: The underlying survival time for the subject.
//'   
//'   - \code{dropoutTime}: The underlying dropout time for the subject.
//'   
//'   - \code{timeUnderObservation}: The time under observation since    
//'   since randomization.
//'   
//'   - \code{event}: Whether the subject experienced the event.
//'   
//'   - \code{dropoutEvent}: Whether the subject dropped out.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Example 1: analyses based on number of events
//' 
//' sim1 = lrsim(kMax = 2, informationTime = c(0.5, 1),
//'              criticalValues = c(2.797, 1.977),
//'              accrualIntensity = 11,
//'              lambda1 = 0.018, lambda2 = 0.030,
//'              accrualDuration = 12,
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
//' sim2 = lrsim(kMax = 2, informationTime = c(0.5, 1),
//'              criticalValues = c(2.797, 1.977),
//'              accrualIntensity = 11,
//'              lambda1 = 0.018, lambda2 = 0.030,
//'              accrualDuration = 12,
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
List lrsim(const int kMax = NA_INTEGER,
           const NumericVector& informationTime = NA_REAL,
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
           const double accrualDuration = NA_REAL,
           const double followupTime = NA_REAL,
           const bool fixedFollowup = 0,
           const double rho1 = 0,
           const double rho2 = 0,
           const IntegerVector& plannedEvents = NA_INTEGER,
           const NumericVector& plannedTime = NA_REAL,
           const int maxNumberOfIterations = 1000,
           const int maxNumberOfRawDatasetsPerStage = 0,
           int seed = NA_INTEGER) {
  
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  NumericVector lambda1x(nsi), lambda2x(nsi);
  NumericVector gamma1x(nsi), gamma2x(nsi);
  
  bool useEvents, eventsNotAchieved;
  NumericVector informationTime1 = clone(informationTime);
  NumericVector futilityBounds1 = clone(futilityBounds);
  
  
  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }
  
  // whether to plan the analyses based on events or calender time
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
  
  
  // set default informationTime
  if (is_false(any(is_na(informationTime)))) {
    if (informationTime.size() != kMax) {
      stop("Invalid length for informationTime");
    } else if (informationTime[0] <= 0) {
      stop("Elements of informationTime must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationTime) <= 0))) {
      stop("Elements of informationTime must be increasing");
    } else if (informationTime[kMax-1] != 1) {
      stop("informationTime must end with 1");
    }
  } else if (useEvents) {
    informationTime1 = as<NumericVector>(plannedEvents)/
      (plannedEvents[kMax-1]+0.0);
  } else {
    informationTime1 = plannedTime/plannedTime[kMax-1];
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
    for (int i=0; i<kMax-1; i++) {
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
  
  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }
  
  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }
  
  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  
  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }
  
  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }
  
  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
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
  } else if (lambda1.size() == nintervals) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }
  
  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nintervals) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }
  
  if (gamma1.size() == 1) {
    gamma1x = rep(gamma1, nsi);
  } else if (gamma1.size() == nintervals) {
    gamma1x = rep(gamma1, nstrata);
  } else if (gamma1.size() == nsi) {
    gamma1x = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }
  
  if (gamma2.size() == 1) {
    gamma2x = rep(gamma2, nsi);
  } else if (gamma2.size() == nintervals) {
    gamma2x = rep(gamma2, nstrata);
  } else if (gamma2.size() == nsi) {
    gamma2x = gamma2;
  } else {
    stop("Invalid length for gamma2");
  }
  
  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  }
  
  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
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
  int i, iter, j, k, h, nevents, nstages, stopStage, nsub;
  int index1=0, index2=0;
  
  double u, enrollt, time, uscore1, vscore1;
  
  
  // maximum number of subjects to enroll
  int m = accrualTime.size();
  double s = 0;
  for (i=0; i<m; i++) {
    if (i<m-1 && accrualTime[i+1] < accrualDuration) {
      s += accrualIntensity[i]*(accrualTime[i+1] - accrualTime[i]);
    } else {
      s += accrualIntensity[i]*(accrualDuration - accrualTime[i]);
      break;
    }
  }
  int n = floor(s + 0.5);
  
  
  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n), sortedIndex(n), 
  stratumSorted(n), treatmentGroupSorted(n);
  
  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n),
  timeUnderObservation(n), totalTime(n), totalt(n), 
  observationTime(n), timeUnderObservationSorted(n);
  
  LogicalVector event(n), dropoutEvent(n), sub(n), eventSorted(n);
  
  
  // stratum information 
  IntegerVector b1(nstrata), b2(nstrata), n1(nstrata), n1x(nstrata), 
  n2(nstrata), nt(nstrata), ntx(nstrata); 
  
  NumericVector km(nstrata), w(nstrata);
  NumericVector cumStratumFraction = cumsum(stratumFraction);

  
  // within-stratum hazard rates 
  NumericVector lam1(nintervals), lam2(nintervals); 
  NumericVector gam1(nintervals), gam2(nintervals);
  
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
  IntegerVector stopStagex = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector subjectIdx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector arrivalTimex = NumericVector(nrow1, NA_REAL);
  IntegerVector stratumx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector treatmentGroupx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector survivalTimex = NumericVector(nrow1, NA_REAL);
  NumericVector dropoutTimex = NumericVector(nrow1, NA_REAL);
  NumericVector observationTimex = NumericVector(nrow1, NA_REAL);
  NumericVector timeUnderObservationx = NumericVector(nrow1, NA_REAL);
  LogicalVector eventx = LogicalVector(nrow1, NA_LOGICAL);
  LogicalVector dropoutEventx = LogicalVector(nrow1, NA_LOGICAL);
  
  // cache for the simulation-level summary data to extract
  int nrow2 = kMax*maxNumberOfIterations;
  
  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector stopStagey = IntegerVector(nrow2, NA_INTEGER);
  LogicalVector eventsNotAchievedy = LogicalVector(nrow2, NA_LOGICAL);
  IntegerVector stageNumbery = IntegerVector(nrow2, NA_INTEGER);
  NumericVector analysisTimey = NumericVector(nrow2, NA_REAL);
  IntegerVector accruals1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector accruals2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalAccrualsy = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalEventsy = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalDropoutsy = IntegerVector(nrow2, NA_INTEGER);
  NumericVector uscorey = NumericVector(nrow2, NA_REAL);
  NumericVector vscorey = NumericVector(nrow2, NA_REAL);
  NumericVector logRankStatisticy = NumericVector(nrow2, NA_REAL);
  LogicalVector rejectPerStagey = LogicalVector(nrow2, NA_LOGICAL);
  LogicalVector futilityPerStagey = LogicalVector(nrow2, NA_LOGICAL);
  
  
  // total alpha to adjust the critical value at the final stage
  NumericVector lb = rep(-6.0, kMax);
  NumericVector theta = rep(0.0, kMax);
  List p1 = exitprobcpp(criticalValues, lb, theta, informationTime1);
  double alpha = sum(NumericVector(p1[0]));
  
  
  // set up random seed
  if (seed==NA_INTEGER) {
    set_seed(std::time(0));
  } else {
    set_seed(seed);
  }
  
  
  for (iter=0; iter<maxNumberOfIterations; iter++) {
    int nstops = 0;
    
    b1.fill(allocation1);
    b2.fill(allocation2);
    
    enrollt = 0;
    for (i=0; i<n; i++) {
      
      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp(u, accrualTime, accrualIntensity, enrollt);
      arrivalTime[i] = enrollt;
      
      // generate stratum information
      u = R::runif(0,1);
      for (j=0; j<nstrata; j++) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }
      
      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+0.0)) {
        treatmentGroup[i] = 1;
        b1[j]--;
      } else {
        treatmentGroup[i] = 2;
        b2[j]--;
      }
      
      // start a new block after depleting the current block
      if (b1[j]+b2[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
      }

      // stratum-specific hazard rates for event and dropout      
      Range jj = Range(j*nintervals, (j+1)*nintervals-1);
      
      lam1 = lambda1x[jj];
      lam2 = lambda2x[jj];
      gam1 = gamma1x[jj];
      gam2 = gamma2x[jj];
      
      // generate survival time
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        survivalTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, lam1, 0);
      } else {
        survivalTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, lam2, 0);
      }
      
      // generate dropout time
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, gam1, 0);
      } else {
        dropoutTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, gam2, 0);
      }
      
      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime[i] < dropoutTime[i] &&
            survivalTime[i] < followupTime) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else if (dropoutTime[i] < survivalTime[i] &&
          dropoutTime[i] < followupTime) {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        } else {
          timeUnderObservation[i] = followupTime;
          event[i] = 0;
          dropoutEvent[i] = 0;
        }
      } else { // variable follow-up design
        if (survivalTime[i] < dropoutTime[i]) {
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
      
      for (j=0; j<kMax; j++) {
        if (plannedEvents[j] >= nevents) {
          nstages = j+1;
          break;
        }
      }
      
      
      if (j==kMax) { // total number of events exceeds planned
        for (k=0; k<nstages; k++) {
          analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
          obsEvents[k] = plannedEvents[k];
        }
      } else {
        for (k=0; k<nstages; k++) {
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
    for (k=0; k<nstages; k++) {
      time = analysisTime[k];
      
      n1.fill(0);  // number of subjects in each stratum by treatment
      n2.fill(0);
      events1[k] = 0;
      events2[k] = 0;
      dropouts1[k] = 0;
      dropouts2[k] = 0;
      
      // censor at analysis time
      for (i=0; i<n; i++) {
        h = stratum[i]-1;
        observationTime[i] = time;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation[i] = time - arrivalTime[i];
          event[i] = 0;
          dropoutEvent[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            n1[h]++;
          } else if (treatmentGroup[i]==2) {
            n2[h]++;
          }
          
          if (fixedFollowup) {
            // the first three cases correspond to arrivalTime[i] + 
            // min(survivalTime[i], dropoutTime[i], followupTime) < time
            if (arrivalTime[i] + survivalTime[i] < time &&
                survivalTime[i] < dropoutTime[i] &&
                survivalTime[i] < followupTime) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] < time &&
              dropoutTime[i] < survivalTime[i] &&
              dropoutTime[i] < followupTime) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else if (arrivalTime[i] + followupTime < time &&
              followupTime < survivalTime[i] &&
              followupTime < dropoutTime[i]) {
              timeUnderObservation[i] = followupTime;
              event[i] = 0;
              dropoutEvent[i] = 0;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime[i] < time &&
                survivalTime[i] < dropoutTime[i]) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] < time &&
              dropoutTime[i] < survivalTime[i]) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          }
          
          if (treatmentGroup[i]==1 && event[i]) events1[k]++;
          if (treatmentGroup[i]==2 && event[i]) events2[k]++;
          if (treatmentGroup[i]==1 && dropoutEvent[i]) dropouts1[k]++;
          if (treatmentGroup[i]==2 && dropoutEvent[i]) dropouts2[k]++;
        }
      }
      
      // number of accrued patients and total number of events
      accruals1[k] = sum(n1);
      accruals2[k] = sum(n2);
      totalAccruals[k] = accruals1[k] + accruals2[k];
      
      totalEvents[k] = events1[k] + events2[k];
      totalDropouts[k] = dropouts1[k] + dropouts2[k];
      
      // order the data by time under observation
      timeUnderObservationSorted = stl_sort(timeUnderObservation);
      sortedIndex = match(timeUnderObservationSorted, timeUnderObservation);
      sortedIndex = sortedIndex - 1;
      eventSorted = event[sortedIndex];
      stratumSorted = stratum[sortedIndex];
      treatmentGroupSorted = treatmentGroup[sortedIndex];
      sub = (timeUnderObservationSorted > 0);
      eventSorted = eventSorted[sub];
      stratumSorted = stratumSorted[sub];
      treatmentGroupSorted = treatmentGroupSorted[sub];
      nsub = eventSorted.size();
      
      // calculate the stratified log-rank test
      uscore1 = 0;
      vscore1 = 0;
      km.fill(1);  // km(t-) estimate by stratum
      for (i=0; i<nsub; i++) {
        h = stratumSorted[i] - 1;
        n1x[h] = n1[h]*hazardRatioH0;
        nt[h] = n1[h] + n2[h];
        ntx[h] = n1x[h] + n2[h];
        
        if (eventSorted[i]) { // at most 1 event can occur at a given time
          w[h] = pow(km[h], rho1)*pow(1-km[h], rho2);
          uscore1 += w[h]*((treatmentGroupSorted[i]==1)-n1x[h]/(ntx[h]+0.0));
          vscore1 += w[h]*w[h]*n1x[h]*n2[h]/(ntx[h]*ntx[h]+0.0);
          km[h] *= (1-1/(nt[h]+0.0)); // update km estimate
        }
        
        // reduce the risk set
        if (treatmentGroupSorted[i]==1) {
          n1[h]--;
        } else {
          n2[h]--;
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
                        for (int i=0; i<nstages-1; i++) {
                          u[i] = criticalValues[i];
                        }
                        u[nstages-1] = aval;
                        NumericVector l = rep(-6.0, nstages);
                        NumericVector theta = rep(0.0, nstages);
                        NumericVector I = as<NumericVector>(obsEvents)[
                        Range(0,nstages-1)];
                        List p2 = exitprobcpp(u, l, theta, I);
                        return sum(NumericVector(p2[0])) - alpha;
                      };
            
            adjCriticalValues[nstages-1] = brent(f, 0, 6, 1e-6);
          } else { // weighted log-rank test
            auto f = [criticalValues, alpha, &vscore, 
                      &nstages](double aval)->double {
                        NumericVector u(nstages);
                        for (int i=0; i<nstages-1; i++) {
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
        nstops++;
        
        if (nstops == 1) { // extract at most one raw data set per iteration
          
          // add raw data to output
          if (niter[k] < maxNumberOfRawDatasetsPerStage) {
            for (i=0; i<n; i++) {
              iterationNumberx[index1] = iter+1;
              stopStagex[index1] = k+1;
              subjectIdx[index1] = i+1;
              arrivalTimex[index1] = arrivalTime[i];
              stratumx[index1] = stratum[i];
              treatmentGroupx[index1] = treatmentGroup[i];
              survivalTimex[index1] = survivalTime[i];
              dropoutTimex[index1] = dropoutTime[i];
              observationTimex[index1] = observationTime[i];
              timeUnderObservationx[index1] = timeUnderObservation[i];
              eventx[index1] = event[i];
              dropoutEventx[index1] = dropoutEvent[i];
              index1++;
            }
            
            // update the number of stage k dataset to extract
            niter[k]++;
          }
          
          stopStage = k+1;
          
        }
        
      }
      
    }
    
    // add summary data to output
    for (k=0; k<nstages; k++) {
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
      index2++;
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
  int nrow3 = stageNumbery.size();
  
  for (i=0; i<nrow3; i++) {
    k = stageNumbery[i] - 1;
    if (stageNumbery[i] == stopStagey[i]) {
      pRejectPerStage[k] += rejectPerStagey[i];
      pFutilityPerStage[k] += futilityPerStagey[i];
    }
    
    nEventsPerStage[k] += totalEventsy[i];
    nDropoutsPerStage[k] += totalDropoutsy[i];
    nSubjectsPerStage[k] += totalAccrualsy[i];
    analysisTimePerStage[k] += analysisTimey[i];
  }
  
  
  for (k=0; k<kMax; k++) {
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
  
  for (i=0; i<nrow3; i++) {
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
    _["accrualDuration"] = accrualDuration,
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
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    survivalTimex = survivalTimex[sub1];
    dropoutTimex = dropoutTimex[sub1];
    observationTimex = observationTimex[sub1];
    timeUnderObservationx = timeUnderObservationx[sub1];
    eventx = eventx[sub1];
    dropoutEventx = dropoutEventx[sub1];
    
    DataFrame rawdata = DataFrame::create(
      _["iterationNumber"] = iterationNumberx,
      _["stopStage"] = stopStagex,
      _["analysisTime"] = observationTimex,
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






//' @title Log-rank test simulation for three arms
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
//' @inheritParams param_accrualDuration
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
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per stage
//'   to extract. Defaults to 1.
//' @param seed The seed to reproduce the simulation results.
//'   The computer clock will be used if left unspecified,
//'
//' @return A list with 2 components:
//'
//' * \code{sumdata}: A data frame of summary data by iteration and stage: 
//' 
//'   - \code{iterationNumber}: The iteration number.
//'   
//'   - \code{eventsNotAchieved}: Whether the target number of events 
//'   is not achieved for the iteration.
//'   
//'   - \code{stageNumber}: The stage number, covering all stages even if 
//'   the trial stops at an interim look.
//'   
//'   - \code{analysisTime}: The time for the stage since trial start.
//'   
//'   - \code{accruals1}: The number of subjects enrolled at the stage for
//'   the active treatment 1 group.
//'   
//'   - \code{accruals2}: The number of subjects enrolled at the stage for
//'   the active treatment 2 group.
//'   
//'   - \code{accruals3}: The number of subjects enrolled at the stage for
//'   the control group.
//'   
//'   - \code{totalAccruals}: The total number of subjects enrolled at 
//'   the stage.
//'   
//'   - \code{events1}: The number of events at the stage for 
//'   the active treatment 1 group.
//'   
//'   - \code{events2}: The number of events at the stage for 
//'   the active treatment 2 group.
//'   
//'   - \code{events3}: The number of events at the stage for 
//'   the control group.
//'   
//'   - \code{totalEvents}: The total number of events at the stage.
//'   
//'   - \code{dropouts1}: The number of dropouts at the stage for 
//'   the active treatment 1 group. 
//'   
//'   - \code{dropouts2}: The number of dropouts at the stage for 
//'   the active treatment 2 group.
//'   
//'   - \code{dropouts3}: The number of dropouts at the stage for 
//'   the control group.
//'   
//'   - \code{totalDropouts}: The total number of dropouts at the stage.
//'   
//'   - \code{logRankStatistic13}: The log-rank test Z-statistic 
//'   comparing the active treatment 1 to the control.
//'   
//'   - \code{logRankStatistic23}: The log-rank test Z-statistic 
//'   comparing the active treatment 2 to the control.
//'   
//'   - \code{logRankStatistic12}: The log-rank test Z-statistic 
//'   comparing the active treatment 1 to the active treatment 2.
//'   
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//' positive integer): A data frame for subject-level data for selected
//' replications, containing the following variables:
//' 
//'   - \code{iterationNumber}: The iteration number.
//'   
//'   - \code{stageNumber}: The stage under consideration.
//'   
//'   - \code{analysisTime}: The time for the stage since trial start.
//'   
//'   - \code{subjectId}: The subject ID.
//'   
//'   - \code{arrivalTime}: The enrollment time for the subject.
//'   
//'   - \code{stratum}: The stratum for the subject.
//'   
//'   - \code{treatmentGroup}: The treatment group (1, 2, or 3) for 
//'   the subject.
//'   
//'   - \code{survivalTime}: The underlying survival time for the subject.
//'   
//'   - \code{dropoutTime}: The underlying dropout time for the subject.
//'   
//'   - \code{timeUnderObservation}: The time under observation since 
//'   since randomization for the subject.
//'   
//'   - \code{event}: Whether the subject experienced the event.
//'   
//'   - \code{dropoutEvent}: Whether the subject dropped out. 
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
//'   accrualDuration = 30.143,
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
List lrsim3a(const int kMax = NA_INTEGER,
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
             const double accrualDuration = NA_REAL,
             const double followupTime = NA_REAL,
             const bool fixedFollowup = 0,
             const double rho1 = 0,
             const double rho2 = 0,
             const IntegerVector& plannedEvents = NA_INTEGER,
             const NumericVector& plannedTime = NA_REAL,
             const int maxNumberOfIterations = 1000,
             const int maxNumberOfRawDatasetsPerStage = 0,
             int seed = NA_INTEGER) {
  
  // check input parameters
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  
  NumericVector lambda1x(nsi), lambda2x(nsi), lambda3x(nsi);
  NumericVector gamma1x(nsi), gamma2x(nsi), gamma3x(nsi);
  
  bool useEvents, eventsNotAchieved;
  
  
  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }
  
  // whether to plan the analyses based on events or calender time
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
  
  
  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }
  
  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }
  
  
  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  
  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }
  
  
  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }
  
  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
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
  } else if (lambda1.size() == nintervals) {
    lambda1x = rep(lambda1, nstrata);
  } else if (lambda1.size() == nsi) {
    lambda1x = lambda1;
  } else {
    stop("Invalid length for lambda1");
  }
  
  
  if (lambda2.size() == 1) {
    lambda2x = rep(lambda2, nsi);
  } else if (lambda2.size() == nintervals) {
    lambda2x = rep(lambda2, nstrata);
  } else if (lambda2.size() == nsi) {
    lambda2x = lambda2;
  } else {
    stop("Invalid length for lambda2");
  }
  
  
  if (lambda3.size() == 1) {
    lambda3x = rep(lambda3, nsi);
  } else if (lambda3.size() == nintervals) {
    lambda3x = rep(lambda3, nstrata);
  } else if (lambda3.size() == nsi) {
    lambda3x = lambda3;
  } else {
    stop("Invalid length for lambda3");
  }
  
  
  
  if (gamma1.size() == 1) {
    gamma1x = rep(gamma1, nsi);
  } else if (gamma1.size() == nintervals) {
    gamma1x = rep(gamma1, nstrata);
  } else if (gamma1.size() == nsi) {
    gamma1x = gamma1;
  } else {
    stop("Invalid length for gamma1");
  }
  
  
  if (gamma2.size() == 1) {
    gamma2x = rep(gamma2, nsi);
  } else if (gamma2.size() == nintervals) {
    gamma2x = rep(gamma2, nstrata);
  } else if (gamma2.size() == nsi) {
    gamma2x = gamma2;
  } else {
    stop("Invalid length for gamma2");
  }
  
  
  if (gamma3.size() == 1) {
    gamma3x = rep(gamma3, nsi);
  } else if (gamma3.size() == nintervals) {
    gamma3x = rep(gamma3, nstrata);
  } else if (gamma3.size() == nsi) {
    gamma3x = gamma3;
  } else {
    stop("Invalid length for gamma3");
  }
  
  
  
  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  }
  
  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
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
  int i, iter, j, k, h, nevents, nstages, nsub;
  int accruals1, accruals2, accruals3, totalAccruals; 
  int events1, events2, events3, totalEvents;
  int dropouts1, dropouts2, dropouts3, totalDropouts;
  int index1=0, index2=0;
  
  double enrollt, u, time;
  double uscore13, uscore23, uscore12;
  double vscore13, vscore23, vscore12;
  
  
  // maximum number of subjects to enroll
  int m = accrualTime.size();
  double s = 0;
  for (i=0; i<m; i++) {
    if (i<m-1 && accrualTime[i+1] < accrualDuration) {
      s += accrualIntensity[i]*(accrualTime[i+1] - accrualTime[i]);
    } else {
      s += accrualIntensity[i]*(accrualDuration - accrualTime[i]);
      break;
    }
  }
  int n = floor(s + 0.5);
  
  
  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n), sortedIndex(n), 
  stratumSorted(n), treatmentGroupSorted(n);
  
  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n), 
  timeUnderObservation(n), totalTime(n), totalt(n), 
  observationTime(n), timeUnderObservationSorted(n);
  
  LogicalVector event(n), dropoutEvent(n), eventac(n), sub(n), eventSorted(n);
  
  
  // stratum information 
  IntegerVector b1(nstrata), b2(nstrata), b3(nstrata);
  IntegerVector n1(nstrata), n2(nstrata), n3(nstrata);
  IntegerVector nt13(nstrata), nt23(nstrata), nt12(nstrata);
  IntegerVector n13a(nstrata), n23a(nstrata), n12a(nstrata);
  IntegerVector nt13a(nstrata), nt23a(nstrata), nt12a(nstrata);
  
  NumericVector km13(nstrata), km23(nstrata), km12(nstrata);
  NumericVector w13(nstrata), w23(nstrata), w12(nstrata);
  NumericVector cumStratumFraction = cumsum(stratumFraction);
  
  // within-stratum hazard rates
  NumericVector lam1(nintervals), lam2(nintervals), lam3(nintervals); 
  NumericVector gam1(nintervals), gam2(nintervals), gam3(nintervals);
  
  
  // stage-wise information
  IntegerVector niter(kMax);
  NumericVector analysisTime(kMax);
  
  
  // cache for the patient-level raw data to extract
  int nrow1 = n*kMax*maxNumberOfRawDatasetsPerStage;
  
  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stageNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector subjectIdx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector arrivalTimex = NumericVector(nrow1, NA_REAL);
  IntegerVector stratumx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector treatmentGroupx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector survivalTimex = NumericVector(nrow1, NA_REAL);
  NumericVector dropoutTimex = NumericVector(nrow1, NA_REAL);
  NumericVector observationTimex = NumericVector(nrow1, NA_REAL);
  NumericVector timeUnderObservationx = NumericVector(nrow1, NA_REAL);
  LogicalVector eventx = LogicalVector(nrow1, NA_LOGICAL);
  LogicalVector dropoutEventx = LogicalVector(nrow1, NA_LOGICAL);
  
  // cache for the simulation-level summary data to extract
  int nrow2 = kMax*maxNumberOfIterations*2;
  
  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  LogicalVector eventsNotAchievedy = LogicalVector(nrow2, NA_LOGICAL);
  IntegerVector stageNumbery = IntegerVector(nrow2, NA_INTEGER);
  NumericVector analysisTimey = NumericVector(nrow2, NA_REAL);
  IntegerVector accruals1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector accruals2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector accruals3y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalAccrualsy = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events3y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalEventsy = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts3y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalDropoutsy = IntegerVector(nrow2, NA_INTEGER);
  NumericVector logRankStatistic13y = NumericVector(nrow2, NA_REAL);
  NumericVector logRankStatistic23y = NumericVector(nrow2, NA_REAL);
  NumericVector logRankStatistic12y = NumericVector(nrow2, NA_REAL);
  
  
  // set up random seed
  if (seed==NA_INTEGER) {
    set_seed(std::time(0));
  } else {
    set_seed(seed);
  }
  

  for (iter=0; iter<maxNumberOfIterations; iter++) {
    
    b1.fill(allocation1);
    b2.fill(allocation2);
    b3.fill(allocation3);
    
    enrollt = 0;
    for (i=0; i<n; i++) {
      
      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp(u, accrualTime, accrualIntensity, enrollt);
      arrivalTime[i] = enrollt;
      
      // generate stratum information
      u = R::runif(0,1);
      for (j=0; j<nstrata; j++) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }
      
      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+b3[j]+0.0)) {
        treatmentGroup[i] = 1;
        b1[j]--;
      } else if (u <= (b1[j]+b2[j])/(b1[j]+b2[j]+b3[j]+0.0)) {
        treatmentGroup[i] = 2;
        b2[j]--;
      } else {
        treatmentGroup[i] = 3;
        b3[j]--;
      }
      
      // start a new block after depleting the current block
      if (b1[j]+b2[j]+b3[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
        b3[j] = allocation3;
      }
      
      // stratum-specific hazard rates for event and dropout
      Range jj = Range(j*nintervals, (j+1)*nintervals-1);
      
      lam1 = lambda1x[jj];
      lam2 = lambda2x[jj];
      lam3 = lambda3x[jj];
      
      gam1 = gamma1x[jj];
      gam2 = gamma2x[jj];
      gam3 = gamma3x[jj];
      
      
      // generate survival time
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        survivalTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, lam1, 0);
      } else if (treatmentGroup[i]==2) {
        survivalTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, lam2, 0);
      } else {
        survivalTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, lam3, 0);
      }
      
      
      // generate dropout time
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, gam1, 0);
      } else if (treatmentGroup[i]==2) {
        dropoutTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, gam2, 0);
      } else {
        dropoutTime[i] = qtpwexpcpp(u, piecewiseSurvivalTime, gam3, 0);
      }
      
      
      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime[i] < dropoutTime[i] &&
            survivalTime[i] < followupTime) {
          timeUnderObservation[i] = survivalTime[i];
          event[i] = 1;
          dropoutEvent[i] = 0;
        } else if (dropoutTime[i] < survivalTime[i] &&
          dropoutTime[i] < followupTime) {
          timeUnderObservation[i] = dropoutTime[i];
          event[i] = 0;
          dropoutEvent[i] = 1;
        } else {
          timeUnderObservation[i] = followupTime;
          event[i] = 0;
          dropoutEvent[i] = 0;
        }
      } else { // variable follow-up design
        if (survivalTime[i] < dropoutTime[i]) {
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
      
      for (j=0; j<kMax; j++) {
        if (plannedEvents[j] >= nevents) {
          nstages = j+1;
          break;
        }
      }
      
      if (j==kMax) { // total number of events exceeds planned
        for (k=0; k<nstages; k++) {
          analysisTime[k] = totalt[plannedEvents[k]-1] + 1e-12;
        }
      } else {
        for (k=0; k<nstages; k++) {
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
    for (k=0; k<nstages; k++) {
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
      for (i=0; i<n; i++) {
        h = stratum[i]-1;
        observationTime[i] = time;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation[i] = time - arrivalTime[i];
          event[i] = 0;
          dropoutEvent[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            n1[h]++;
          } else if (treatmentGroup[i]==2) {
            n2[h]++;
          } else {
            n3[h]++;
          }
          
          if (fixedFollowup) {
            // the first three cases correspond to arrivalTime[i] + 
            // min(survivalTime[i], dropoutTime[i], followupTime) < time
            if (arrivalTime[i] + survivalTime[i] < time &&
                survivalTime[i] < dropoutTime[i] &&
                survivalTime[i] < followupTime) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] < time &&
              dropoutTime[i] < survivalTime[i] &&
              dropoutTime[i] < followupTime) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else if (arrivalTime[i] + followupTime < time &&
              followupTime < survivalTime[i] &&
              followupTime < dropoutTime[i]) {
              timeUnderObservation[i] = followupTime;
              event[i] = 0;
              dropoutEvent[i] = 0;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime[i] < time &&
                survivalTime[i] < dropoutTime[i]) {
              timeUnderObservation[i] = survivalTime[i];
              event[i] = 1;
              dropoutEvent[i] = 0;
            } else if (arrivalTime[i] + dropoutTime[i] < time &&
              dropoutTime[i] < survivalTime[i]) {
              timeUnderObservation[i] = dropoutTime[i];
              event[i] = 0;
              dropoutEvent[i] = 1;
            } else {
              timeUnderObservation[i] = time - arrivalTime[i];
              event[i] = 0;
              dropoutEvent[i] = 0;
            }
          }
          
          if (treatmentGroup[i]==1 && event[i]) events1++;
          if (treatmentGroup[i]==2 && event[i]) events2++;
          if (treatmentGroup[i]==3 && event[i]) events3++;
          if (treatmentGroup[i]==1 && dropoutEvent[i]) dropouts1++;
          if (treatmentGroup[i]==2 && dropoutEvent[i]) dropouts2++;
          if (treatmentGroup[i]==3 && dropoutEvent[i]) dropouts3++;
        }
      }
      
      
      // add raw data to output
      if (niter[k] < maxNumberOfRawDatasetsPerStage) {
        for (i=0; i<n; i++) {
          iterationNumberx[index1] = iter+1;
          stageNumberx[index1] = k+1;
          subjectIdx[index1] = i+1;
          arrivalTimex[index1] = arrivalTime[i];
          stratumx[index1] = stratum[i];
          treatmentGroupx[index1] = treatmentGroup[i];
          observationTimex[index1] = observationTime[i];
          
          survivalTimex[index1] = survivalTime[i];
          dropoutTimex[index1] = dropoutTime[i];
          timeUnderObservationx[index1] = timeUnderObservation[i];
          eventx[index1] = event[i];
          dropoutEventx[index1] = dropoutEvent[i];
          
          index1++;
        }
        
        // update the number of stage k dataset to extract
        niter[k]++;
      }
      
      
      // number of accrued patients and total number of events
      accruals1 = sum(n1);
      accruals2 = sum(n2);
      accruals3 = sum(n3);
      totalAccruals = accruals1 + accruals2 + accruals3;
      
      totalEvents = events1 + events2 + events3;
      totalDropouts = dropouts1 + dropouts2 + dropouts3;
      
      
      // order the data by time under observation
      timeUnderObservationSorted = stl_sort(timeUnderObservation);
      sortedIndex = match(timeUnderObservationSorted, timeUnderObservation);
      sortedIndex = sortedIndex - 1;
      eventSorted = event[sortedIndex];
      stratumSorted = stratum[sortedIndex];
      treatmentGroupSorted = treatmentGroup[sortedIndex];
      sub = (timeUnderObservationSorted > 0);
      eventSorted = eventSorted[sub];
      stratumSorted = stratumSorted[sub];
      treatmentGroupSorted = treatmentGroupSorted[sub];
      nsub = eventSorted.size();
      
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
      for (i=0; i<nsub; i++) {
        h = stratumSorted[i] - 1;
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
                                - n13a[h]/(nt13a[h]+0.0));
          vscore13 += w13[h]*w13[h]*n13a[h]*n3[h]/(nt13a[h]*nt13a[h]+0.0);
          km13[h] *= (1-1/(nt13[h]+0.0)); // update km estimate
        }
        
        if (eventSorted[i] && (treatmentGroupSorted[i]==2 || 
            treatmentGroupSorted[i]==3)) { 
          w23[h] = pow(km23[h], rho1)*pow(1-km23[h], rho2);
          uscore23 += w23[h]*((treatmentGroupSorted[i]==2) 
                                - n23a[h]/(nt23a[h]+0.0));
          vscore23 += w23[h]*w23[h]*n23a[h]*n3[h]/(nt23a[h]*nt23a[h]+0.0);
          km23[h] *= (1-1/(nt23[h]+0.0)); // update km estimate
        }
        
        if (eventSorted[i] && (treatmentGroupSorted[i]==1 || 
            treatmentGroupSorted[i]==2)) { 
          w12[h] = pow(km12[h], rho1)*pow(1-km12[h], rho2);
          uscore12 += w12[h]*((treatmentGroupSorted[i]==1) 
                                - n12a[h]/(nt12a[h]+0.0));
          vscore12 += w12[h]*w12[h]*n12a[h]*n2[h]/(nt12a[h]*nt12a[h]+0.0);
          km12[h] *= (1-1/(nt12[h]+0.0)); // update km estimate
        }
        
        // reduce the risk set
        if (treatmentGroupSorted[i]==1) {
          n1[h]--;
        } else if (treatmentGroupSorted[i]==2) {
          n2[h]--;
        } else {
          n3[h]--;
        }
      }
      
      
      
      // add summary data to output  
      iterationNumbery[index2] = iter+1;
      eventsNotAchievedy[index2] = eventsNotAchieved;
      stageNumbery[index2] = k+1;
      analysisTimey[index2] = analysisTime[k];
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
      
      index2++;
      
      
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
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    observationTimex = observationTimex[sub1];
    survivalTimex = survivalTimex[sub1];
    dropoutTimex = dropoutTimex[sub1];
    timeUnderObservationx = timeUnderObservationx[sub1];
    eventx = eventx[sub1];
    dropoutEventx = dropoutEventx[sub1];
    
    DataFrame rawdata = DataFrame::create(
      _["iterationNumber"] = iterationNumberx,
      _["stageNumber"] = stageNumberx,
      _["analysisTime"] = observationTimex,
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





//' @title Log-rank test simulation for two endpoints and two arms
//' @description Performs simulation for two-endpoint two-arm group 
//' sequential trials based on weighted log-rank test. The first \code{kMaxe1} 
//' looks are driven by the total number of PFS events in two arms 
//' combined, and the subsequent looks are driven by the total 
//' number of OS events in two arms combined. Alternatively, 
//' the analyses can be planned to occur at specified calendar times.
//' 
//' @inheritParams param_kMax
//' @param kMaxe1 Number of stages with timing determined by PFS events. 
//'   Ranges from 0 (none) to \code{kMax}. 
//' @param hazardRatioH0e1 Hazard ratio under the null hypothesis for the 
//'   active treatment vs control for endpoint 1 (PFS). Defaults to 1 for 
//'   superiority test.
//' @param hazardRatioH0e2 Hazard ratio under the null hypothesis for the 
//'   active treatment vs control for endpoint 2 (OS). Defaults to 1 for 
//'   superiority test.
//' @param allocation1 Number of subjects in the treatment group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in the control group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param rho The correlation coefficient for the standard bivariate normal
//'   random variables used to generate time to disease progression and time 
//'   to death using the inverse CDF method.
//' @param lambda1e1 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for the treatment group and endpoint 1 (PFS).
//' @param lambda2e1 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for the control group and endpoint 1 (PFS).
//' @param lambda1e2 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for the treatment group and endpoint 2 (OS).
//' @param lambda2e2 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for the control group and endpoint 2 (OS).
//' @param gamma1e1 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for the treatment group and endpoint 1 (PFS).
//' @param gamma2e1 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for the control group and endpoint 1 (PFS).
//' @param gamma1e2 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for the treatment group and endpoint 2 (OS).
//' @param gamma2e2 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for the control group and endpoint 2 (OS).
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @param plannedEvents The planned cumulative total number of PFS events at 
//'   Look 1 to Look kMaxe1 and the planned cumulative total number of OS 
//'   events at Look kMaxe1+1 to Look kMax.
//' @param plannedTime The calendar times for the analyses. To use calendar  
//'   time to plan the analyses, \code{plannedEvents} should be missing.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per stage
//'   to extract. Defaults to 1.
//' @param seed The seed to reproduce the simulation results.
//'   The computer clock will be used if left unspecified,
//'
//' @return A list with 2 components:
//'
//' * \code{sumdata}: A data frame of summary data by iteration and stage: 
//' 
//'   - \code{iterationNumber}: The iteration number.
//'   
//'   - \code{eventsNotAchieved}: Whether the target number of events 
//'   is not achieved for the iteration.
//'   
//'   - \code{stageNumber}: The stage number, covering all stages even if 
//'   the trial stops at an interim look.
//'   
//'   - \code{analysisTime}: The time for the stage since trial start.
//'   
//'   - \code{accruals1}: The number of subjects enrolled at the stage for
//'   the treatment group.
//'   
//'   - \code{accruals2}: The number of subjects enrolled at the stage for
//'   the control group.
//'   
//'   - \code{totalAccruals}: The total number of subjects enrolled at 
//'   the stage.
//'   
//'   - \code{endpoint}: The endpoint (1 or 2) under consideration.
//'   
//'   - \code{events1}: The number of events at the stage for 
//'   the treatment group.
//'   
//'   - \code{events2}: The number of events at the stage for 
//'   the control group.
//'   
//'   - \code{totalEvents}: The total number of events at the stage.
//'   
//'   - \code{dropouts1}: The number of dropouts at the stage for 
//'   the treatment group. 
//'   
//'   - \code{dropouts2}: The number of dropouts at the stage for 
//'   the control group.
//'   
//'   - \code{totalDropouts}: The total number of dropouts at the stage.
//'   
//'   - \code{logRankStatistic}: The log-rank test Z-statistic for 
//'   the endpoint.
//'
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//' positive integer): A data frame for subject-level data for selected
//' replications, containing the following variables:
//' 
//'   - \code{iterationNumber}: The iteration number.
//'   
//'   - \code{stageNumber}: The stage under consideration.
//'   
//'   - \code{analysisTime}: The time for the stage since trial start.
//'   
//'   - \code{subjectId}: The subject ID.
//'   
//'   - \code{arrivalTime}: The enrollment time for the subject.
//'   
//'   - \code{stratum}: The stratum for the subject.
//'   
//'   - \code{treatmentGroup}: The treatment group (1 or 2) for the subject.
//'   
//'   - \code{survivalTime1}: The underlying survival time for 
//'   event endpoint 1 for the subject.
//'   
//'   - \code{dropoutTime1}: The underlying dropout time for 
//'   event endpoint 1 for the subject.
//'   
//'   - \code{timeUnderObservation1}: The time under observation since 
//'   since randomization for event endpoint 1 for the subject.
//'   
//'   - \code{event1}: Whether the subject experienced event endpoint 1.
//'   
//'   - \code{dropoutEvent1}: Whether the subject dropped out for 
//'   endpoint 1.
//'   
//'   - \code{survivalTime2}: The underlying survival time for 
//'   event endpoint 2 for the subject.
//'   
//'   - \code{dropoutTime2}: The underlying dropout time for 
//'   event endpoint 2 for the subject.
//'   
//'   - \code{timeUnderObservation2}: The time under observation since 
//'   since randomization for event endpoint 2 for the subject.
//'   
//'   - \code{event2}: Whether the subject experienced event endpoint 2.
//'   
//'   - \code{dropoutEvent2}: Whether the subject dropped out for 
//'   endpoint 2.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//' 
//' @examples
//' 
//' sim1 = lrsim2e(
//'   kMax = 3, 
//'   kMaxe1 = 2,
//'   allocation1 = 2, 
//'   allocation2 = 1, 
//'   accrualTime = c(0, 8), 
//'   accrualIntensity = c(10, 28),
//'   piecewiseSurvivalTime = 0, 
//'   rho = 0,
//'   lambda1e1 = log(2)/12*0.60,
//'   lambda2e1 = log(2)/12,
//'   lambda1e2 = log(2)/30*0.65, 
//'   lambda2e2 = log(2)/30,
//'   accrualDuration = 20.143,
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
List lrsim2e(const int kMax = NA_INTEGER,
             const int kMaxe1 = NA_INTEGER,
             const double hazardRatioH0e1 = 1,
             const double hazardRatioH0e2 = 1,
             const int allocation1 = 1,
             const int allocation2 = 1,
             const NumericVector& accrualTime = 0,
             const NumericVector& accrualIntensity = NA_REAL,
             const NumericVector& piecewiseSurvivalTime = 0,
             const NumericVector& stratumFraction = 1,
             const double rho = 0, 
             const NumericVector& lambda1e1 = NA_REAL,
             const NumericVector& lambda2e1 = NA_REAL,
             const NumericVector& lambda1e2 = NA_REAL,
             const NumericVector& lambda2e2 = NA_REAL,
             const NumericVector& gamma1e1 = 0,
             const NumericVector& gamma2e1 = 0,
             const NumericVector& gamma1e2 = 0,
             const NumericVector& gamma2e2 = 0,
             const double accrualDuration = NA_REAL,
             const double followupTime = NA_REAL,
             const bool fixedFollowup = 0,
             const double rho1 = 0,
             const double rho2 = 0,
             const IntegerVector& plannedEvents = NA_INTEGER,
             const NumericVector& plannedTime = NA_REAL,
             const int maxNumberOfIterations = 1000,
             const int maxNumberOfRawDatasetsPerStage = 0,
             int seed = NA_INTEGER) {
  
  // check input parameters
  int kMaxe1x = kMaxe1;
  
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  
  NumericVector lambda1e1x(nsi), lambda2e1x(nsi);
  NumericVector lambda1e2x(nsi), lambda2e2x(nsi);
  NumericVector gamma1e1x(nsi), gamma2e1x(nsi);
  NumericVector gamma1e2x(nsi), gamma2e2x(nsi);
  NumericVector lambda1e1d(nsi), lambda2e1d(nsi);
  NumericVector gamma1e1d(nsi), gamma2e1d(nsi);
  
  bool useEvents, eventsNotAchieved;
  
  
  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }
  
  if (kMaxe1 < 0) {
    kMaxe1x = kMax;
  }
  
  if (kMaxe1x > kMax) {
    stop("kMaxe1 must be less than or equal to kMax");
  }
  
  // whether to plan the analyses based on events or calender time
  if (is_false(any(is_na(plannedEvents)))) {
    useEvents = 1;
    
    if (plannedEvents[0] <= 0) {
      stop("Elements of plannedEvents must be positive");
    }
    
    if (plannedEvents.size() != kMax) {
      stop("Invalid length for plannedEvents");
    }
    
    if (kMaxe1x > 1) {
      IntegerVector plannedEvents1 = plannedEvents[Range(0,kMaxe1x-1)];
      if (is_true(any(diff(plannedEvents1) <= 0))) {
        stop("plannedEvents for endpoint 1 must be increasing");
      }
    }
    
    if (kMax - kMaxe1x > 1) {
      IntegerVector plannedEvents2 = plannedEvents[Range(kMaxe1x, kMax-1)];
      if (is_true(any(diff(plannedEvents2) <= 0))) {
        stop("plannedEvents for endpoint 2 must be increasing");
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
  
  
  if (hazardRatioH0e1 <= 0) {
    stop("hazardRatioH0e1 must be positive");
  }
  
  if (hazardRatioH0e2 <= 0) {
    stop("hazardRatioH0e2 must be positive");
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
  
  
  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }
  
  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }
  
  
  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  
  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }
  
  
  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }
  
  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }
  
  
  if (rho <= -1 || rho >= 1) {
    stop("rho must lie in (-1, 1)");
  }
  
  
  if (is_true(any(lambda1e1 < 0))) {
    stop("lambda1e1 must be non-negative");
  }
  
  if (is_true(any(lambda2e1 < 0))) {
    stop("lambda2e1 must be non-negative");
  }
  
  if (is_true(any(lambda1e2 < 0))) {
    stop("lambda1e2 must be non-negative");
  }
  
  if (is_true(any(lambda2e2 < 0))) {
    stop("lambda2e2 must be non-negative");
  }
  
  
  if (is_true(any(gamma1e1 < 0))) {
    stop("gamma1e1 must be non-negative");
  }
  
  if (is_true(any(gamma2e1 < 0))) {
    stop("gamma2e1 must be non-negative");
  }
  
  if (is_true(any(gamma1e2 < 0))) {
    stop("gamma1e2 must be non-negative");
  }
  
  if (is_true(any(gamma2e2 < 0))) {
    stop("gamma2e2 must be non-negative");
  }
  
  
  if (lambda1e1.size() == 1) {
    lambda1e1x = rep(lambda1e1, nsi);
  } else if (lambda1e1.size() == nintervals) {
    lambda1e1x = rep(lambda1e1, nstrata);
  } else if (lambda1e1.size() == nsi) {
    lambda1e1x = lambda1e1;
  } else {
    stop("Invalid length for lambda1e1");
  }
  
  if (lambda2e1.size() == 1) {
    lambda2e1x = rep(lambda2e1, nsi);
  } else if (lambda2e1.size() == nintervals) {
    lambda2e1x = rep(lambda2e1, nstrata);
  } else if (lambda2e1.size() == nsi) {
    lambda2e1x = lambda2e1;
  } else {
    stop("Invalid length for lambda2e1");
  }
  
  
  
  if (lambda1e2.size() == 1) {
    lambda1e2x = rep(lambda1e2, nsi);
  } else if (lambda1e2.size() == nintervals) {
    lambda1e2x = rep(lambda1e2, nstrata);
  } else if (lambda1e2.size() == nsi) {
    lambda1e2x = lambda1e2;
  } else {
    stop("Invalid length for lambda1e2");
  }
  
  if (lambda2e2.size() == 1) {
    lambda2e2x = rep(lambda2e2, nsi);
  } else if (lambda2e2.size() == nintervals) {
    lambda2e2x = rep(lambda2e2, nstrata);
  } else if (lambda2e2.size() == nsi) {
    lambda2e2x = lambda2e2;
  } else {
    stop("Invalid length for lambda2e2");
  }
  
  
  
  if (gamma1e1.size() == 1) {
    gamma1e1x = rep(gamma1e1, nsi);
  } else if (gamma1e1.size() == nintervals) {
    gamma1e1x = rep(gamma1e1, nstrata);
  } else if (gamma1e1.size() == nsi) {
    gamma1e1x = gamma1e1;
  } else {
    stop("Invalid length for gamma1e1");
  }
  
  if (gamma2e1.size() == 1) {
    gamma2e1x = rep(gamma2e1, nsi);
  } else if (gamma2e1.size() == nintervals) {
    gamma2e1x = rep(gamma2e1, nstrata);
  } else if (gamma2e1.size() == nsi) {
    gamma2e1x = gamma2e1;
  } else {
    stop("Invalid length for gamma2e1");
  }
  
  
  if (gamma1e2.size() == 1) {
    gamma1e2x = rep(gamma1e2, nsi);
  } else if (gamma1e2.size() == nintervals) {
    gamma1e2x = rep(gamma1e2, nstrata);
  } else if (gamma1e2.size() == nsi) {
    gamma1e2x = gamma1e2;
  } else {
    stop("Invalid length for gamma1e2");
  }
  
  if (gamma2e2.size() == 1) {
    gamma2e2x = rep(gamma2e2, nsi);
  } else if (gamma2e2.size() == nintervals) {
    gamma2e2x = rep(gamma2e2, nstrata);
  } else if (gamma2e2.size() == nsi) {
    gamma2e2x = gamma2e2;
  } else {
    stop("Invalid length for gamma2e2");
  }
  
  
  
  lambda1e1d = lambda1e1x - lambda1e2x;
  lambda2e1d = lambda2e1x - lambda2e2x;
  gamma1e1d = gamma1e1x - gamma1e2x;
  gamma2e1d = gamma2e1x - gamma2e2x;
  
  if (is_true(any(lambda1e1d < 0))) {
    stop("lambda1e1 must be greater than or equal to lambda1e2");
  }
  
  if (is_true(any(lambda2e1d < 0))) {
    stop("lambda2e1 must be greater than or equal to lambda2e2");
  }
  
  
  if (is_true(any(gamma1e1d < 0))) {
    stop("gamma1e1 must be greater than or equal to gamma1e2");
  }
  
  if (is_true(any(gamma2e1d < 0))) {
    stop("gamma2e1 must be greater than or equal to gamma2e2");
  }
  
  
  
  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  }
  
  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
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
  int i, iter, j, k, h, nevents1, nevents2, nstages, nsub;
  int accruals1, accruals2, totalAccruals;
  int events1e1, events2e1, totalEventse1;
  int events1e2, events2e2, totalEventse2;
  int dropouts1e1, dropouts2e1, totalDropoutse1;
  int dropouts1e2, dropouts2e2, totalDropoutse2;
  int index1=0, index2=0;
  
  double enrollt, u, u1, u2, time, uscore, vscore;
  
  
  // maximum number of subjects to enroll
  int m = accrualTime.size();
  double s = 0;
  for (i=0; i<m; i++) {
    if (i<m-1 && accrualTime[i+1] < accrualDuration) {
      s += accrualIntensity[i]*(accrualTime[i+1] - accrualTime[i]);
    } else {
      s += accrualIntensity[i]*(accrualDuration - accrualTime[i]);
      break;
    }
  }
  int n = floor(s + 0.5);
  
  
  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n), sortedIndex(n), 
  stratumSorted(n), treatmentGroupSorted(n);
  
  NumericVector arrivalTime(n), survivalTime1(n), survivalTime2(n), 
  dropoutTime1(n), dropoutTime2(n), timeUnderObservation1(n), 
  timeUnderObservation2(n), totalTime1(n), totalTime2(n), 
  totalt1(n), totalt2(n), observationTime(n), timeUnderObservationSorted(n);
  
  LogicalVector event1(n), event2(n), dropoutEvent1(n), dropoutEvent2(n), 
  event1ac(n), event2ac(n), sub(n), eventSorted(n);
  
  
  // stratum information 
  IntegerVector b1(nstrata), b2(nstrata); 
  IntegerVector n1(nstrata), n2(nstrata), nt(nstrata);
  
  // original copy of n1 and n2 when looping over the endpoints
  IntegerVector n1x(nstrata), n2x(nstrata);  
  
  // hazardRatioH0 adjusted n1 and nt for calculating the log-rank statistic
  IntegerVector n1a(nstrata), nta(nstrata);

  NumericVector km(nstrata), w(nstrata);
  NumericVector cumStratumFraction = cumsum(stratumFraction);
  
  
  // within-stratum hazard rates
  NumericVector lam1e1(nintervals), lam2e1(nintervals); 
  NumericVector lam1e2(nintervals), lam2e2(nintervals); 
  NumericVector gam1e1(nintervals), gam2e1(nintervals); 
  NumericVector gam1e2(nintervals), gam2e2(nintervals);
  
  
  // stage-wise information
  IntegerVector niter(kMax);
  NumericVector analysisTime(kMax);
  
  
  // cache for the patient-level raw data to extract
  int nrow1 = n*kMax*maxNumberOfRawDatasetsPerStage;
  
  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stageNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector subjectIdx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector arrivalTimex = NumericVector(nrow1, NA_REAL);
  IntegerVector stratumx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector treatmentGroupx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector survivalTime1x = NumericVector(nrow1, NA_REAL);
  NumericVector survivalTime2x = NumericVector(nrow1, NA_REAL);
  NumericVector dropoutTime1x = NumericVector(nrow1, NA_REAL);
  NumericVector dropoutTime2x = NumericVector(nrow1, NA_REAL);
  NumericVector observationTimex = NumericVector(nrow1, NA_REAL);
  NumericVector timeUnderObservation1x = NumericVector(nrow1, NA_REAL);
  NumericVector timeUnderObservation2x = NumericVector(nrow1, NA_REAL);
  LogicalVector event1x = LogicalVector(nrow1, NA_LOGICAL);
  LogicalVector event2x = LogicalVector(nrow1, NA_LOGICAL);
  LogicalVector dropoutEvent1x = LogicalVector(nrow1, NA_LOGICAL);
  LogicalVector dropoutEvent2x = LogicalVector(nrow1, NA_LOGICAL);
  
  
  // cache for the simulation-level summary data to extract
  int nrow2 = kMax*maxNumberOfIterations*2;
  
  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  LogicalVector eventsNotAchievedy = LogicalVector(nrow2, NA_LOGICAL);
  IntegerVector stageNumbery = IntegerVector(nrow2, NA_INTEGER);
  NumericVector analysisTimey = NumericVector(nrow2, NA_REAL);
  IntegerVector accruals1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector accruals2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalAccrualsy = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector endpointy = IntegerVector(nrow2, NA_INTEGER);  
  IntegerVector events1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalEventsy = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalDropoutsy = IntegerVector(nrow2, NA_INTEGER);
  NumericVector logRankStatisticy = NumericVector(nrow2, NA_REAL);
  
  
  // set up random seed
  if (seed==NA_INTEGER) {
    set_seed(std::time(0));
  } else {
    set_seed(seed);
  }
  
  
  // simulation 
  for (iter=0; iter<maxNumberOfIterations; iter++) {
    
    b1.fill(allocation1);
    b2.fill(allocation2);
    
    enrollt = 0;
    for (i=0; i<n; i++) {
      
      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp(u, accrualTime, accrualIntensity, enrollt);
      arrivalTime[i] = enrollt;
      
      // generate stratum information
      u = R::runif(0,1);
      for (j=0; j<nstrata; j++) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }
      
      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+0.0)) {
        treatmentGroup[i] = 1;
        b1[j]--;
      } else {
        treatmentGroup[i] = 2;
        b2[j]--;
      }
      
      // start a new block after depleting the current block
      if (b1[j]+b2[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
      }
      
      // stratum-specific hazard rates for event and dropout
      Range jj = Range(j*nintervals, (j+1)*nintervals-1);
      
      lam1e1 = lambda1e1d[jj];
      lam2e1 = lambda2e1d[jj];
      
      lam1e2 = lambda1e2x[jj];
      lam2e2 = lambda2e2x[jj];
      
      gam1e1 = gamma1e1d[jj];
      gam2e1 = gamma2e1d[jj];
      
      gam1e2 = gamma1e2x[jj];
      gam2e2 = gamma2e2x[jj];
      
      // standard bivariate normal with correlation rho
      u1 = R::rnorm(0,1);
      u2 = R::rnorm(rho*u1, sqrt(1-rho*rho));
      
      // transform to uniform
      u1 = R::pnorm(u1, 0, 1, 1, 0);
      u2 = R::pnorm(u2, 0, 1, 1, 0);
      
      // generate survival times
      if (treatmentGroup[i]==1) {
        survivalTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, lam1e1, 0);
        survivalTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, lam1e2, 0);
      } else {
        survivalTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, lam2e1, 0);
        survivalTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, lam2e2, 0);
      }
      // PFS includes death
      survivalTime1[i] = std::min(survivalTime1[i], survivalTime2[i]);
      
      
      // generate dropout times
      u1 = R::runif(0,1);
      u2 = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, gam1e1, 0);
        dropoutTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, gam1e2, 0);
      } else {
        dropoutTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, gam2e1, 0);
        dropoutTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, gam2e2, 0);
      }
      // whatever censors OS will also censors PFS
      dropoutTime1[i] = std::min(dropoutTime1[i], dropoutTime2[i]);
      
      
      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime1[i] < dropoutTime1[i] &&
            survivalTime1[i] < followupTime) {
          timeUnderObservation1[i] = survivalTime1[i];
          event1[i] = 1;
          dropoutEvent1[i] = 0;
        } else if (dropoutTime1[i] < survivalTime1[i] &&
          dropoutTime1[i] < followupTime) {
          timeUnderObservation1[i] = dropoutTime1[i];
          event1[i] = 0;
          dropoutEvent1[i] = 1;
        } else {
          timeUnderObservation1[i] = followupTime;
          event1[i] = 0;
          dropoutEvent1[i] = 0;
        }
      } else {
        if (survivalTime1[i] < dropoutTime1[i]) {
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
        if (survivalTime2[i] < dropoutTime2[i] &&
            survivalTime2[i] < followupTime) {
          timeUnderObservation2[i] = survivalTime2[i];
          event2[i] = 1;
          dropoutEvent2[i] = 0;
        } else if (dropoutTime2[i] < survivalTime2[i] &&
          dropoutTime2[i] < followupTime) {
          timeUnderObservation2[i] = dropoutTime2[i];
          event2[i] = 0;
          dropoutEvent2[i] = 1;
        } else {
          timeUnderObservation2[i] = followupTime;
          event2[i] = 0;
          dropoutEvent2[i] = 0;
        }
      } else {
        if (survivalTime2[i] < dropoutTime2[i]) {
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
      
      int j1 = kMaxe1x, j2 = kMax - kMaxe1x;
      
      // PFS looks
      if (kMaxe1x > 0) {
        for (j1=0; j1<kMaxe1x; j1++) {
          if (plannedEvents[j1] >= nevents1) {
            break;
          }
        }
        
        if (j1==kMaxe1x) { // total number of PFS events exceeds planned
          for (k=0; k<kMaxe1x; k++) {
            analysisTime[k] = totalt1[plannedEvents[k]-1] + 1e-12;
          }
        } else {
          for (k=0; k<=j1; k++) {
            if (k < j1) {
              analysisTime[k] = totalt1[plannedEvents[k]-1] + 1e-12;
            } else {
              analysisTime[k] = totalt1[nevents1-1] + 1e-12;
            }
          }
        }
      }
      
      // OS looks
      NumericVector analysisTime2(kMax - kMaxe1x);
      
      if (kMax > kMaxe1x) {
        for (j2=0; j2<kMax-kMaxe1x; j2++) {
          if (plannedEvents[kMaxe1x+j2] >= nevents2) {
            break;
          }
        }
        
        if (j2==kMax-kMaxe1x) { // total number of OS events exceeds planned
          for (k=0; k<kMax-kMaxe1x; k++) {
            analysisTime2[k] = totalt2[plannedEvents[kMaxe1x+k]-1] + 1e-12;
          }
        } else {
          for (k=0; k<=j2; k++) {
            if (k < j2) {
              analysisTime2[k] = totalt2[plannedEvents[kMaxe1x+k]-1] + 1e-12;
            } else {
              analysisTime2[k] = totalt2[nevents2-1] + 1e-12;
            }
          }
        }
      }
      
      // determine the number of looks and timing of the looks
      if (kMaxe1x == 0) { // all looks based on OS events
        if (j2 == kMax - kMaxe1x) {
          nstages = kMax - kMaxe1x;
        } else {
          nstages = j2 + 1;
        }
        
        for (k=0; k<nstages; k++) {
          analysisTime[k] = analysisTime2[k];
        }
      } else if (kMax == kMaxe1x) { // all looks based on PFS events
        if (j1 == kMaxe1x) {
          nstages = kMaxe1x;
        } else {
          nstages = j1 + 1;
        }
      } else {
        if (analysisTime2[kMax-kMaxe1x-1] > analysisTime[kMaxe1x-1]) {
          // only OS looks that occur after the last PFS look contribute
          int l = which_min(analysisTime2 > analysisTime[kMaxe1x-1]);
          nstages = kMax-l;
          for (k=kMaxe1x; k<kMax-l; k++) {
            analysisTime[k] = analysisTime2[k-kMaxe1x+l];
          }
        } else {
          if (j1 == kMaxe1x) {
            nstages = kMaxe1x;
          } else {
            nstages = j1 + 1;
          }
        }
      }
      
      // whether the target PFS and OS events are achieved
      if (kMaxe1x > 0 && nevents1 < plannedEvents[kMaxe1x-1]) {
        eventsNotAchieved = 1;
      } else if (nevents2 < plannedEvents[kMax-1]) {
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
    for (k=0; k<nstages; k++) {
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
      for (i=0; i<n; i++) {
        h = stratum[i]-1;
        observationTime[i] = time;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation1[i] = time - arrivalTime[i];
          event1[i] = 0;
          dropoutEvent1[i] = 0;
          
          timeUnderObservation2[i] = time - arrivalTime[i];
          event2[i] = 0;
          dropoutEvent2[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            n1x[h]++;
          } else {
            n2x[h]++;
          }
          
          
          // censored time for endpoint 1
          if (fixedFollowup) {
            if (arrivalTime[i] + survivalTime1[i] < time &&
                survivalTime1[i] < dropoutTime1[i] &&
                survivalTime1[i] < followupTime) {
              timeUnderObservation1[i] = survivalTime1[i];
              event1[i] = 1;
              dropoutEvent1[i] = 0;
            } else if (arrivalTime[i] + dropoutTime1[i] < time &&
              dropoutTime1[i] < survivalTime1[i] &&
              dropoutTime1[i] < followupTime) {
              timeUnderObservation1[i] = dropoutTime1[i];
              event1[i] = 0;
              dropoutEvent1[i] = 1;
            } else if (arrivalTime[i] + followupTime < time &&
              followupTime < survivalTime1[i] &&
              followupTime < dropoutTime1[i]) {
              timeUnderObservation1[i] = followupTime;
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            } else {
              timeUnderObservation1[i] = time - arrivalTime[i];
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime1[i] < time &&
                survivalTime1[i] < dropoutTime1[i]) {
              timeUnderObservation1[i] = survivalTime1[i];
              event1[i] = 1;
              dropoutEvent1[i] = 0;
            } else if (arrivalTime[i] + dropoutTime1[i] < time &&
              dropoutTime1[i] < survivalTime1[i]) {
              timeUnderObservation1[i] = dropoutTime1[i];
              event1[i] = 0;
              dropoutEvent1[i] = 1;
            } else {
              timeUnderObservation1[i] = time - arrivalTime[i];
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            }
          }
          
          if (treatmentGroup[i]==1 && event1[i]) events1e1++;
          if (treatmentGroup[i]==2 && event1[i]) events2e1++;
          if (treatmentGroup[i]==1 && dropoutEvent1[i]) dropouts1e1++;
          if (treatmentGroup[i]==2 && dropoutEvent1[i]) dropouts2e1++;
          
          
          // censored time for endpoint 2
          if (fixedFollowup) {
            if (arrivalTime[i] + survivalTime2[i] < time &&
                survivalTime2[i] < dropoutTime2[i] &&
                survivalTime2[i] < followupTime) {
              timeUnderObservation2[i] = survivalTime2[i];
              event2[i] = 1;
              dropoutEvent2[i] = 0;
            } else if (arrivalTime[i] + dropoutTime2[i] < time &&
              dropoutTime2[i] < survivalTime2[i] &&
              dropoutTime2[i] < followupTime) {
              timeUnderObservation2[i] = dropoutTime2[i];
              event2[i] = 0;
              dropoutEvent2[i] = 1;
            } else if (arrivalTime[i] + followupTime < time &&
              followupTime < survivalTime2[i] &&
              followupTime < dropoutTime2[i]) {
              timeUnderObservation2[i] = followupTime;
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            } else {
              timeUnderObservation2[i] = time - arrivalTime[i];
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime2[i] < time &&
                survivalTime2[i] < dropoutTime2[i]) {
              timeUnderObservation2[i] = survivalTime2[i];
              event2[i] = 1;
              dropoutEvent2[i] = 0;
            } else if (arrivalTime[i] + dropoutTime2[i] < time &&
              dropoutTime2[i] < survivalTime2[i]) {
              timeUnderObservation2[i] = dropoutTime2[i];
              event2[i] = 0;
              dropoutEvent2[i] = 1;
            } else {
              timeUnderObservation2[i] = time - arrivalTime[i];
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            }
          }
          
          if (treatmentGroup[i]==1 && event2[i]) events1e2++;
          if (treatmentGroup[i]==2 && event2[i]) events2e2++;
          if (treatmentGroup[i]==1 && dropoutEvent2[i]) dropouts1e2++;
          if (treatmentGroup[i]==2 && dropoutEvent2[i]) dropouts2e2++;
        }
      }
      
      
      // add raw data to output
      if (niter[k] < maxNumberOfRawDatasetsPerStage) {
        for (i=0; i<n; i++) {
          iterationNumberx[index1] = iter+1;
          stageNumberx[index1] = k+1;
          subjectIdx[index1] = i+1;
          arrivalTimex[index1] = arrivalTime[i];
          stratumx[index1] = stratum[i];
          treatmentGroupx[index1] = treatmentGroup[i];
          observationTimex[index1] = observationTime[i];
          
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
          
          index1++;
        }
        
        // update the number of stage k dataset to extract
        niter[k]++;
      }
      
      
      // number of accrued patients and total number of events
      accruals1 = sum(n1x);
      accruals2 = sum(n2x);
      totalAccruals = accruals1 + accruals2;
      
      totalEventse1 = events1e1 + events2e1;
      totalDropoutse1 = dropouts1e1 + dropouts2e1;
      
      totalEventse2 = events1e2 + events2e2;
      totalDropoutse2 = dropouts1e2 + dropouts2e2;
      
      
      for (int endpoint=1; endpoint<=2; endpoint++) {
        n1 = clone(n1x);
        n2 = clone(n2x);
        
        double hazardRatioH0;
        if (endpoint == 1) {
          hazardRatioH0 = hazardRatioH0e1;
        } else {
          hazardRatioH0 = hazardRatioH0e2;
        }
        
        // order the data by time under observation
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
        
        stratumSorted = stratum[sortedIndex];
        treatmentGroupSorted = treatmentGroup[sortedIndex];
        sub = (timeUnderObservationSorted > 0);
        eventSorted = eventSorted[sub];
        stratumSorted = stratumSorted[sub];
        treatmentGroupSorted = treatmentGroupSorted[sub];
        nsub = eventSorted.size();
        
        // calculate the stratified log-rank test
        uscore = 0;
        vscore = 0;
        km.fill(1);
        for (i=0; i<nsub; i++) {
          h = stratumSorted[i] - 1;
          nt[h] = n1[h] + n2[h];
          
          n1a[h] = n1[h]*hazardRatioH0;
          nta[h] = n1a[h] + n2[h];
          
          if (eventSorted[i]) { 
            w[h] = pow(km[h], rho1)*pow(1-km[h], rho2);
            uscore += w[h]*((treatmentGroupSorted[i]==1) 
                              - n1a[h]/(nta[h]+0.0));
            vscore += w[h]*w[h]*n1a[h]*n2[h]/(nta[h]*nta[h]+0.0);
            km[h] *= (1-1/(nt[h]+0.0)); // update km estimate
          }
          
          // reduce the risk set
          if (treatmentGroupSorted[i]==1) {
            n1[h]--;
          } else {
            n2[h]--;
          }
        }
        
        
        
        // add summary data to output  
        iterationNumbery[index2] = iter+1;
        eventsNotAchievedy[index2] = eventsNotAchieved;
        stageNumbery[index2] = k+1;
        analysisTimey[index2] = analysisTime[k];
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
        
        logRankStatisticy[index2] = uscore/sqrt(vscore);
        index2++;
        
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
    _["logRankStatistic"] = logRankStatisticy);
  
  
  List result;
  
  if (maxNumberOfRawDatasetsPerStage > 0) {
    LogicalVector sub1 = !is_na(iterationNumberx);
    iterationNumberx = iterationNumberx[sub1];
    stageNumberx = stageNumberx[sub1];
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    observationTimex = observationTimex[sub1];
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
      _["analysisTime"] = observationTimex,
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



//' @title Log-rank test simulation for two endpoints and three arms
//' @description Performs simulation for two-endpoint three-arm group 
//' sequential trials based on weighted log-rank test. The first \code{kMaxe1} 
//' looks are driven by the total number of PFS events in Arm A 
//' and Arm C combined, and the subsequent looks are driven by the total 
//' number of OS events in Arm A and Arm C combined. Alternatively, 
//' the analyses can be planned to occur at specified calendar times. 
//' 
//' @inheritParams param_kMax
//' @param kMaxe1 Number of stages with timing determined by PFS events. 
//'   Ranges from 0 (none) to \code{kMax}. 
//' @param hazardRatioH013e1 Hazard ratio under the null hypothesis for arm 1 
//'   vs arm 3 for endpoint 1 (PFS). Defaults to 1 for superiority test.
//' @param hazardRatioH023e1 Hazard ratio under the null hypothesis for arm 2 
//'   vs arm 3 for endpoint 1 (PFS). Defaults to 1 for superiority test.
//' @param hazardRatioH012e1 Hazard ratio under the null hypothesis for arm 1 
//'   vs arm 2 for endpoint 1 (PFS). Defaults to 1 for superiority test.
//' @param hazardRatioH013e2 Hazard ratio under the null hypothesis for arm 1 
//'   vs arm 3 for endpoint 2 (OS). Defaults to 1 for superiority test.
//' @param hazardRatioH023e2 Hazard ratio under the null hypothesis for arm 2 
//'   vs arm 3 for endpoint 2 (OS). Defaults to 1 for superiority test.
//' @param hazardRatioH012e2 Hazard ratio under the null hypothesis for arm 1 
//'   vs arm 2 for endpoint 2 (OS). Defaults to 1 for superiority test.
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
//' @param rho The correlation coefficient for the standard bivariate normal
//'   random variables used to generate time to disease progression and time 
//'   to death using the inverse CDF method.
//' @param lambda1e1 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for arm 1 and endpoint 1 (PFS).
//' @param lambda2e1 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for arm 2 and endpoint 1 (PFS).
//' @param lambda3e1 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for arm 3 and endpoint 1 (PFS).
//' @param lambda1e2 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for arm 1 and endpoint 2 (OS).
//' @param lambda2e2 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for arm 2 and endpoint 2 (OS).
//' @param lambda3e2 A vector of hazard rates for the event in each analysis  
//'   time interval by stratum for arm 3 and endpoint 2 (OS).
//' @param gamma1e1 The hazard rate for exponential dropout. A vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 1 and endpoint 1 (PFS).
//' @param gamma2e1 The hazard rate for exponential dropout. A vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 2 and endpoint 1 (PFS).
//' @param gamma3e1 The hazard rate for exponential dropout. A vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 3 and endpoint 1 (PFS).
//' @param gamma1e2 The hazard rate for exponential dropout. A vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 1 and endpoint 2 (OS).
//' @param gamma2e2 The hazard rate for exponential dropout. A vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 2 and endpoint 2 (OS).
//' @param gamma3e2 The hazard rate for exponential dropout. A vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 3 and endpoint 2 (OS).
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @param plannedEvents The planned cumulative total number of PFS events at 
//'   Look 1 to Look kMaxe1 for Arms A and C combined and the planned 
//'   cumulative total number of OS events at Look kMaxe1+1 to Look kMax 
//'   for Arms A and C combined.
//' @param plannedTime The calendar times for the analyses. To use calendar  
//'   time to plan the analyses, \code{plannedEvents} should be missing.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per stage
//'   to extract. Defaults to 1.
//' @param seed The seed to reproduce the simulation results.
//'   The computer clock will be used if left unspecified,
//'
//' @return A list with 2 components:
//'
//' * \code{sumdata}: A data frame of summary data by iteration and stage: 
//' 
//'   - \code{iterationNumber}: The iteration number.
//'   
//'   - \code{eventsNotAchieved}: Whether the target number of events 
//'   is not achieved for the iteration.
//'   
//'   - \code{stageNumber}: The stage number, covering all stages even if 
//'   the trial stops at an interim look.
//'   
//'   - \code{analysisTime}: The time for the stage since trial start.
//'   
//'   - \code{accruals1}: The number of subjects enrolled at the stage for
//'   the active treatment 1 group.
//'   
//'   - \code{accruals2}: The number of subjects enrolled at the stage for
//'   the active treatment 2 group.
//'   
//'   - \code{accruals3}: The number of subjects enrolled at the stage for
//'   the control group.
//'   
//'   - \code{totalAccruals}: The total number of subjects enrolled at 
//'   the stage.
//'   
//'   - \code{endpoint}: The endpoint (1 or 2) under consideration.
//'   
//'   - \code{events1}: The number of events at the stage for 
//'   the active treatment 1 group.
//'   
//'   - \code{events2}: The number of events at the stage for 
//'   the active treatment 2 group.
//'   
//'   - \code{events3}: The number of events at the stage for 
//'   the control group.
//'   
//'   - \code{totalEvents}: The total number of events at the stage.
//'   
//'   - \code{dropouts1}: The number of dropouts at the stage for 
//'   the active treatment 1 group. 
//'   
//'   - \code{dropouts2}: The number of dropouts at the stage for 
//'   the active treatment 2 group.
//'   
//'   - \code{dropouts3}: The number of dropouts at the stage for 
//'   the control group.
//'   
//'   - \code{totalDropouts}: The total number of dropouts at the stage.
//'   
//'   - \code{logRankStatistic13}: The log-rank test Z-statistic 
//'   comparing the active treatment 1 to the control for the endpoint.
//'   
//'   - \code{logRankStatistic23}: The log-rank test Z-statistic 
//'   comparing the active treatment 2 to the control for the endpoint.
//'   
//'   - \code{logRankStatistic12}: The log-rank test Z-statistic 
//'   comparing the active treatment 1 to the active treatment 2 
//'   for the endpoint.
//'   
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//' positive integer): A data frame for subject-level data for selected
//' replications, containing the following variables:
//' 
//'   - \code{iterationNumber}: The iteration number.
//'   
//'   - \code{stageNumber}: The stage under consideration.
//'   
//'   - \code{analysisTime}: The time for the stage since trial start.
//'   
//'   - \code{subjectId}: The subject ID.
//'   
//'   - \code{arrivalTime}: The enrollment time for the subject.
//'   
//'   - \code{stratum}: The stratum for the subject.
//'   
//'   - \code{treatmentGroup}: The treatment group (1, 2, or 3) for 
//'   the subject.
//'   
//'   - \code{survivalTime1}: The underlying survival time for 
//'   event endpoint 1 for the subject.
//'   
//'   - \code{dropoutTime1}: The underlying dropout time for 
//'   event endpoint 1 for the subject.
//'   
//'   - \code{timeUnderObservation1}: The time under observation since 
//'   since randomization for event endpoint 1 for the subject.
//'   
//'   - \code{event1}: Whether the subject experienced event endpoint 1.
//'   
//'   - \code{dropoutEvent1}: Whether the subject dropped out for 
//'   endpoint 1.
//'   
//'   - \code{survivalTime2}: The underlying survival time for 
//'   event endpoint 2 for the subject.
//'   
//'   - \code{dropoutTime2}: The underlying dropout time for 
//'   event endpoint 2 for the subject.
//'   
//'   - \code{timeUnderObservation2}: The time under observation since 
//'   since randomization for event endpoint 2 for the subject.
//'   
//'   - \code{event2}: Whether the subject experienced event endpoint 2.
//'   
//'   - \code{dropoutEvent2}: Whether the subject dropped out for 
//'   endpoint 2.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//' 
//' @examples
//' 
//' sim1 = lrsim2e3a(
//'   kMax = 3, 
//'   kMaxe1 = 2,
//'   allocation1 = 2, 
//'   allocation2 = 2, 
//'   allocation3 = 1,
//'   accrualTime = c(0, 8), 
//'   accrualIntensity = c(10, 28),
//'   piecewiseSurvivalTime = 0, 
//'   rho = 0,
//'   lambda1e1 = log(2)/12*0.60,
//'   lambda2e1 = log(2)/12*0.70,
//'   lambda3e1 = log(2)/12,
//'   lambda1e2 = log(2)/30*0.65, 
//'   lambda2e2 = log(2)/30*0.75,
//'   lambda3e2 = log(2)/30,
//'   accrualDuration = 30.143,
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
List lrsim2e3a(const int kMax = NA_INTEGER,
               const int kMaxe1 = NA_INTEGER,
               const double hazardRatioH013e1 = 1,
               const double hazardRatioH023e1 = 1,
               const double hazardRatioH012e1 = 1,
               const double hazardRatioH013e2 = 1,
               const double hazardRatioH023e2 = 1,
               const double hazardRatioH012e2 = 1,
               const int allocation1 = 1,
               const int allocation2 = 1,
               const int allocation3 = 1,
               const NumericVector& accrualTime = 0,
               const NumericVector& accrualIntensity = NA_REAL,
               const NumericVector& piecewiseSurvivalTime = 0,
               const NumericVector& stratumFraction = 1,
               const double rho = 0, 
               const NumericVector& lambda1e1 = NA_REAL,
               const NumericVector& lambda2e1 = NA_REAL,
               const NumericVector& lambda3e1 = NA_REAL,
               const NumericVector& lambda1e2 = NA_REAL,
               const NumericVector& lambda2e2 = NA_REAL,
               const NumericVector& lambda3e2 = NA_REAL,
               const NumericVector& gamma1e1 = 0,
               const NumericVector& gamma2e1 = 0,
               const NumericVector& gamma3e1 = 0,
               const NumericVector& gamma1e2 = 0,
               const NumericVector& gamma2e2 = 0,
               const NumericVector& gamma3e2 = 0,
               const double accrualDuration = NA_REAL,
               const double followupTime = NA_REAL,
               const bool fixedFollowup = 0,
               const double rho1 = 0,
               const double rho2 = 0,
               const IntegerVector& plannedEvents = NA_INTEGER,
               const NumericVector& plannedTime = NA_REAL,
               const int maxNumberOfIterations = 1000,
               const int maxNumberOfRawDatasetsPerStage = 0,
               int seed = NA_INTEGER) {
  
  // check input parameters
  int kMaxe1x = kMaxe1;
  
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  
  NumericVector lambda1e1x(nsi), lambda2e1x(nsi), lambda3e1x(nsi); 
  NumericVector lambda1e2x(nsi), lambda2e2x(nsi), lambda3e2x(nsi); 
  NumericVector gamma1e1x(nsi), gamma2e1x(nsi), gamma3e1x(nsi);
  NumericVector gamma1e2x(nsi), gamma2e2x(nsi), gamma3e2x(nsi);
  NumericVector lambda1e1d(nsi), lambda2e1d(nsi), lambda3e1d(nsi); 
  NumericVector gamma1e1d(nsi), gamma2e1d(nsi), gamma3e1d(nsi);
  
  bool useEvents, eventsNotAchieved;
  
  
  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }
  
  if (kMaxe1 < 0) {
    kMaxe1x = kMax;
  }
  
  if (kMaxe1x > kMax) {
    stop("kMaxe1 must be less than or equal to kMax");
  }
  
  
  // whether to plan the analyses based on events or calender time
  if (is_false(any(is_na(plannedEvents)))) {
    useEvents = 1;
    
    if (plannedEvents[0] <= 0) {
      stop("Elements of plannedEvents must be positive");
    }
    
    if (plannedEvents.size() != kMax) {
      stop("Invalid length for plannedEvents");
    }
    
    if (kMaxe1x > 1) {
      IntegerVector plannedEvents1 = plannedEvents[Range(0,kMaxe1x-1)];
      if (is_true(any(diff(plannedEvents1) <= 0))) {
        stop("plannedEvents for endpoint 1 must be increasing");
      }
    }
    
    if (kMax - kMaxe1x > 1) {
      IntegerVector plannedEvents2 = plannedEvents[Range(kMaxe1x, kMax-1)];
      if (is_true(any(diff(plannedEvents2) <= 0))) {
        stop("plannedEvents for endpoint 2 must be increasing");
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
  
  if (hazardRatioH013e1 <= 0) {
    stop("hazardRatioH013e1 must be positive");
  }
  
  if (hazardRatioH023e1 <= 0) {
    stop("hazardRatioH023e1 must be positive");
  }
  
  if (hazardRatioH012e1 <= 0) {
    stop("hazardRatioH012e1 must be positive");
  }
  
  if (hazardRatioH013e2 <= 0) {
    stop("hazardRatioH013e2 must be positive");
  }
  
  if (hazardRatioH023e2 <= 0) {
    stop("hazardRatioH023e2 must be positive");
  }
  
  if (hazardRatioH012e2 <= 0) {
    stop("hazardRatioH012e2 must be positive");
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
  
  
  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }
  
  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }
  
  
  if (piecewiseSurvivalTime[0] != 0) {
    stop("piecewiseSurvivalTime must start with 0");
  }
  
  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
    stop("piecewiseSurvivalTime should be increasing");
  }
  
  
  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }
  
  if (sum(stratumFraction) != 1) {
    stop("stratumFraction must sum to 1");
  }
  
  
  if (rho <= -1 || rho >= 1) {
    stop("rho must lie in (-1, 1)");
  }
  
  
  if (is_true(any(lambda1e1 < 0))) {
    stop("lambda1e1 must be non-negative");
  }
  
  if (is_true(any(lambda2e1 < 0))) {
    stop("lambda2e1 must be non-negative");
  }
  
  if (is_true(any(lambda3e1 < 0))) {
    stop("lambda3e1 must be non-negative");
  }
  
  if (is_true(any(lambda1e2 < 0))) {
    stop("lambda1e2 must be non-negative");
  }
  
  if (is_true(any(lambda2e2 < 0))) {
    stop("lambda2e2 must be non-negative");
  }
  
  if (is_true(any(lambda3e2 < 0))) {
    stop("lambda3e2 must be non-negative");
  }
  
  
  if (is_true(any(gamma1e1 < 0))) {
    stop("gamma1e1 must be non-negative");
  }
  
  if (is_true(any(gamma2e1 < 0))) {
    stop("gamma2e1 must be non-negative");
  }
  
  if (is_true(any(gamma3e1 < 0))) {
    stop("gamma3e1 must be non-negative");
  }
  
  if (is_true(any(gamma1e2 < 0))) {
    stop("gamma1e2 must be non-negative");
  }
  
  if (is_true(any(gamma2e2 < 0))) {
    stop("gamma2e2 must be non-negative");
  }
  
  if (is_true(any(gamma3e2 < 0))) {
    stop("gamma3e2 must be non-negative");
  }  
  
  
  if (lambda1e1.size() == 1) {
    lambda1e1x = rep(lambda1e1, nsi);
  } else if (lambda1e1.size() == nintervals) {
    lambda1e1x = rep(lambda1e1, nstrata);
  } else if (lambda1e1.size() == nsi) {
    lambda1e1x = lambda1e1;
  } else {
    stop("Invalid length for lambda1e1");
  }
  
  if (lambda2e1.size() == 1) {
    lambda2e1x = rep(lambda2e1, nsi);
  } else if (lambda2e1.size() == nintervals) {
    lambda2e1x = rep(lambda2e1, nstrata);
  } else if (lambda2e1.size() == nsi) {
    lambda2e1x = lambda2e1;
  } else {
    stop("Invalid length for lambda2e1");
  }
  
  if (lambda3e1.size() == 1) {
    lambda3e1x = rep(lambda3e1, nsi);
  } else if (lambda3e1.size() == nintervals) {
    lambda3e1x = rep(lambda3e1, nstrata);
  } else if (lambda3e1.size() == nsi) {
    lambda3e1x = lambda3e1;
  } else {
    stop("Invalid length for lambda3e1");
  }
  
  
  if (lambda1e2.size() == 1) {
    lambda1e2x = rep(lambda1e2, nsi);
  } else if (lambda1e2.size() == nintervals) {
    lambda1e2x = rep(lambda1e2, nstrata);
  } else if (lambda1e2.size() == nsi) {
    lambda1e2x = lambda1e2;
  } else {
    stop("Invalid length for lambda1e2");
  }
  
  if (lambda2e2.size() == 1) {
    lambda2e2x = rep(lambda2e2, nsi);
  } else if (lambda2e2.size() == nintervals) {
    lambda2e2x = rep(lambda2e2, nstrata);
  } else if (lambda2e2.size() == nsi) {
    lambda2e2x = lambda2e2;
  } else {
    stop("Invalid length for lambda2e2");
  }
  
  if (lambda3e2.size() == 1) {
    lambda3e2x = rep(lambda3e2, nsi);
  } else if (lambda3e2.size() == nintervals) {
    lambda3e2x = rep(lambda3e2, nstrata);
  } else if (lambda3e2.size() == nsi) {
    lambda3e2x = lambda3e2;
  } else {
    stop("Invalid length for lambda3e2");
  }
  
  
  if (gamma1e1.size() == 1) {
    gamma1e1x = rep(gamma1e1, nsi);
  } else if (gamma1e1.size() == nintervals) {
    gamma1e1x = rep(gamma1e1, nstrata);
  } else if (gamma1e1.size() == nsi) {
    gamma1e1x = gamma1e1;
  } else {
    stop("Invalid length for gamma1e1");
  }
  
  if (gamma2e1.size() == 1) {
    gamma2e1x = rep(gamma2e1, nsi);
  } else if (gamma2e1.size() == nintervals) {
    gamma2e1x = rep(gamma2e1, nstrata);
  } else if (gamma2e1.size() == nsi) {
    gamma2e1x = gamma2e1;
  } else {
    stop("Invalid length for gamma2e1");
  }
  
  if (gamma3e1.size() == 1) {
    gamma3e1x = rep(gamma3e1, nsi);
  } else if (gamma3e1.size() == nintervals) {
    gamma3e1x = rep(gamma3e1, nstrata);
  } else if (gamma3e1.size() == nsi) {
    gamma3e1x = gamma3e1;
  } else {
    stop("Invalid length for gamma3e1");
  }
  
  
  if (gamma1e2.size() == 1) {
    gamma1e2x = rep(gamma1e2, nsi);
  } else if (gamma1e2.size() == nintervals) {
    gamma1e2x = rep(gamma1e2, nstrata);
  } else if (gamma1e2.size() == nsi) {
    gamma1e2x = gamma1e2;
  } else {
    stop("Invalid length for gamma1e2");
  }
  
  if (gamma2e2.size() == 1) {
    gamma2e2x = rep(gamma2e2, nsi);
  } else if (gamma2e2.size() == nintervals) {
    gamma2e2x = rep(gamma2e2, nstrata);
  } else if (gamma2e2.size() == nsi) {
    gamma2e2x = gamma2e2;
  } else {
    stop("Invalid length for gamma2e2");
  }
  
  if (gamma3e2.size() == 1) {
    gamma3e2x = rep(gamma3e2, nsi);
  } else if (gamma3e2.size() == nintervals) {
    gamma3e2x = rep(gamma3e2, nstrata);
  } else if (gamma3e2.size() == nsi) {
    gamma3e2x = gamma3e2;
  } else {
    stop("Invalid length for gamma3e2");
  }
  
  
  lambda1e1d = lambda1e1x - lambda1e2x;
  lambda2e1d = lambda2e1x - lambda2e2x;
  lambda3e1d = lambda3e1x - lambda3e2x;
  gamma1e1d = gamma1e1x - gamma1e2x;
  gamma2e1d = gamma2e1x - gamma2e2x;
  gamma3e1d = gamma3e1x - gamma3e2x;
  
  if (is_true(any(lambda1e1d < 0))) {
    stop("lambda1e1 must be greater than or equal to lambda1e2");
  }
  
  if (is_true(any(lambda2e1d < 0))) {
    stop("lambda2e1 must be greater than or equal to lambda2e2");
  }
  
  if (is_true(any(lambda3e1d < 0))) {
    stop("lambda3e1 must be greater than or equal to lambda3e2");
  }
  
  if (is_true(any(gamma1e1d < 0))) {
    stop("gamma1e1 must be greater than or equal to gamma1e2");
  }
  
  if (is_true(any(gamma2e1d < 0))) {
    stop("gamma2e1 must be greater than or equal to gamma2e2");
  }
  
  if (is_true(any(gamma3e1d < 0))) {
    stop("gamma3e1 must be greater than or equal to gamma3e2");
  }
  
  
  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  }
  
  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
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
  int i, iter, j, k, h, nevents1, nevents2, nstages, nsub; 
  int accruals1, accruals2, accruals3, totalAccruals; 
  int events1e1, events2e1, events3e1, totalEventse1; 
  int events1e2, events2e2, events3e2, totalEventse2; 
  int dropouts1e1, dropouts2e1, dropouts3e1, totalDropoutse1;
  int dropouts1e2, dropouts2e2, dropouts3e2, totalDropoutse2;
  int index1=0, index2=0;
  
  double enrollt, u, u1, u2, time;
  double uscore13, uscore23, uscore12; 
  double vscore13, vscore23, vscore12;
  
  
  // maximum number of subjects to enroll
  int m = accrualTime.size();
  double s = 0;
  for (i=0; i<m; i++) {
    if (i<m-1 && accrualTime[i+1] < accrualDuration) {
      s += accrualIntensity[i]*(accrualTime[i+1] - accrualTime[i]);
    } else {
      s += accrualIntensity[i]*(accrualDuration - accrualTime[i]);
      break;
    }
  }
  int n = floor(s + 0.5);
  
  
  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n), sortedIndex(n), 
  stratumSorted(n), treatmentGroupSorted(n);
  
  NumericVector arrivalTime(n), survivalTime1(n), survivalTime2(n), 
  dropoutTime1(n), dropoutTime2(n), timeUnderObservation1(n), 
  timeUnderObservation2(n), totalTime1(n), totalTime2(n), 
  totalt1(n), totalt2(n), observationTime(n), timeUnderObservationSorted(n);
  
  LogicalVector event1(n), event2(n), dropoutEvent1(n), dropoutEvent2(n), 
  event1ac(n), event2ac(n), sub(n), eventSorted(n);
  
  
  // stratum information 
  IntegerVector b1(nstrata), b2(nstrata), b3(nstrata); 
  IntegerVector n1(nstrata), n2(nstrata), n3(nstrata); 
  IntegerVector nt13(nstrata), nt23(nstrata), nt12(nstrata);

  // original copy of n1 and n2 when looping over the endpoints
  IntegerVector n1x(nstrata), n2x(nstrata), n3x(nstrata); 

  // hazardRatioH0 adjusted at risk for calculating the log-rank statistic
  IntegerVector n13a(nstrata), n23a(nstrata), n12a(nstrata);
  IntegerVector nt13a(nstrata), nt23a(nstrata), nt12a(nstrata);
  
  NumericVector km13(nstrata), km23(nstrata), km12(nstrata);
  NumericVector w13(nstrata), w23(nstrata), w12(nstrata);
  NumericVector cumStratumFraction = cumsum(stratumFraction);
  
  // within-stratum hazard rates
  NumericVector lam1e1(nintervals), lam2e1(nintervals), lam3e1(nintervals); 
  NumericVector lam1e2(nintervals), lam2e2(nintervals), lam3e2(nintervals); 
  NumericVector gam1e1(nintervals), gam2e1(nintervals), gam3e1(nintervals);
  NumericVector gam1e2(nintervals), gam2e2(nintervals), gam3e2(nintervals);
  
  
  // stage-wise information
  IntegerVector niter(kMax);
  NumericVector analysisTime(kMax);
  
  
  // cache for the patient-level raw data to extract
  int nrow1 = n*kMax*maxNumberOfRawDatasetsPerStage;
  
  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stageNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector subjectIdx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector arrivalTimex = NumericVector(nrow1, NA_REAL);
  IntegerVector stratumx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector treatmentGroupx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector survivalTime1x = NumericVector(nrow1, NA_REAL);
  NumericVector survivalTime2x = NumericVector(nrow1, NA_REAL);
  NumericVector dropoutTime1x = NumericVector(nrow1, NA_REAL);
  NumericVector dropoutTime2x = NumericVector(nrow1, NA_REAL);
  NumericVector observationTimex = NumericVector(nrow1, NA_REAL);
  NumericVector timeUnderObservation1x = NumericVector(nrow1, NA_REAL);
  NumericVector timeUnderObservation2x = NumericVector(nrow1, NA_REAL);
  LogicalVector event1x = LogicalVector(nrow1, NA_LOGICAL);
  LogicalVector event2x = LogicalVector(nrow1, NA_LOGICAL);
  LogicalVector dropoutEvent1x = LogicalVector(nrow1, NA_LOGICAL);
  LogicalVector dropoutEvent2x = LogicalVector(nrow1, NA_LOGICAL);
  
  
  // cache for the simulation-level summary data to extract
  int nrow2 = kMax*maxNumberOfIterations*2;
  
  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  LogicalVector eventsNotAchievedy = LogicalVector(nrow2, NA_LOGICAL);
  IntegerVector stageNumbery = IntegerVector(nrow2, NA_INTEGER);
  NumericVector analysisTimey = NumericVector(nrow2, NA_REAL);
  IntegerVector accruals1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector accruals2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector accruals3y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalAccrualsy = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector endpointy = IntegerVector(nrow2, NA_INTEGER);  
  IntegerVector events1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector events3y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalEventsy = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts1y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts2y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector dropouts3y = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector totalDropoutsy = IntegerVector(nrow2, NA_INTEGER);
  NumericVector logRankStatistic13y = NumericVector(nrow2, NA_REAL);
  NumericVector logRankStatistic23y = NumericVector(nrow2, NA_REAL);
  NumericVector logRankStatistic12y = NumericVector(nrow2, NA_REAL);
  
  
  // set up random seed
  if (seed==NA_INTEGER) {
    set_seed(std::time(0));
  } else {
    set_seed(seed);
  }
  
  
  // simulation 
  for (iter=0; iter<maxNumberOfIterations; iter++) {
    
    b1.fill(allocation1);
    b2.fill(allocation2);
    b3.fill(allocation3);
    
    enrollt = 0;
    for (i=0; i<n; i++) {
      
      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp(u, accrualTime, accrualIntensity, enrollt);
      arrivalTime[i] = enrollt;
      
      // generate stratum information
      u = R::runif(0,1);
      for (j=0; j<nstrata; j++) {
        if (cumStratumFraction[j] > u) {
          stratum[i] = j+1;
          break;
        }
      }
      
      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1[j]/(b1[j]+b2[j]+b3[j]+0.0)) {
        treatmentGroup[i] = 1;
        b1[j]--;
      } else if (u <= (b1[j]+b2[j])/(b1[j]+b2[j]+b3[j]+0.0)) {
        treatmentGroup[i] = 2;
        b2[j]--;
      } else {
        treatmentGroup[i] = 3;
        b3[j]--;
      }
      
      // start a new block after depleting the current block
      if (b1[j]+b2[j]+b3[j]==0) {
        b1[j] = allocation1;
        b2[j] = allocation2;
        b3[j] = allocation3;
      }
      
      // stratum-specific hazard rates for event and dropout
      Range jj = Range(j*nintervals, (j+1)*nintervals-1);
      
      lam1e1 = lambda1e1d[jj];
      lam2e1 = lambda2e1d[jj];
      lam3e1 = lambda3e1d[jj];
      
      lam1e2 = lambda1e2x[jj];
      lam2e2 = lambda2e2x[jj];
      lam3e2 = lambda3e2x[jj];
      
      gam1e1 = gamma1e1d[jj];
      gam2e1 = gamma2e1d[jj];
      gam3e1 = gamma3e1d[jj];
      
      gam1e2 = gamma1e2x[jj];
      gam2e2 = gamma2e2x[jj];
      gam3e2 = gamma3e2x[jj];
      
      // standard bivariate normal with correlation rho
      u1 = R::rnorm(0,1);
      u2 = R::rnorm(rho*u1, sqrt(1-rho*rho));
      
      // transform to uniform
      u1 = R::pnorm(u1, 0, 1, 1, 0);
      u2 = R::pnorm(u2, 0, 1, 1, 0);
      
      // generate survival times
      if (treatmentGroup[i]==1) {
        survivalTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, lam1e1, 0);
        survivalTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, lam1e2, 0);
      } else if (treatmentGroup[i]==2) {
        survivalTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, lam2e1, 0);
        survivalTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, lam2e2, 0);
      } else {
        survivalTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, lam3e1, 0);
        survivalTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, lam3e2, 0);
      }
      // PFS includes death
      survivalTime1[i] = std::min(survivalTime1[i], survivalTime2[i]);
      
      
      // generate dropout times
      u1 = R::runif(0,1);
      u2 = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, gam1e1, 0);
        dropoutTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, gam1e2, 0);
      } else if (treatmentGroup[i]==2) {
        dropoutTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, gam2e1, 0);
        dropoutTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, gam2e2, 0);
      } else {
        dropoutTime1[i] = qtpwexpcpp(u1, piecewiseSurvivalTime, gam3e1, 0);
        dropoutTime2[i] = qtpwexpcpp(u2, piecewiseSurvivalTime, gam3e2, 0);
      }
      // whatever censors OS will also censors PFS
      dropoutTime1[i] = std::min(dropoutTime1[i], dropoutTime2[i]);
      
      
      // initial observed time and event indicator
      if (fixedFollowup) { // fixed follow-up design
        if (survivalTime1[i] < dropoutTime1[i] &&
            survivalTime1[i] < followupTime) {
          timeUnderObservation1[i] = survivalTime1[i];
          event1[i] = 1;
          dropoutEvent1[i] = 0;
        } else if (dropoutTime1[i] < survivalTime1[i] &&
          dropoutTime1[i] < followupTime) {
          timeUnderObservation1[i] = dropoutTime1[i];
          event1[i] = 0;
          dropoutEvent1[i] = 1;
        } else {
          timeUnderObservation1[i] = followupTime;
          event1[i] = 0;
          dropoutEvent1[i] = 0;
        }
      } else {
        if (survivalTime1[i] < dropoutTime1[i]) {
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
        if (survivalTime2[i] < dropoutTime2[i] &&
            survivalTime2[i] < followupTime) {
          timeUnderObservation2[i] = survivalTime2[i];
          event2[i] = 1;
          dropoutEvent2[i] = 0;
        } else if (dropoutTime2[i] < survivalTime2[i] &&
          dropoutTime2[i] < followupTime) {
          timeUnderObservation2[i] = dropoutTime2[i];
          event2[i] = 0;
          dropoutEvent2[i] = 1;
        } else {
          timeUnderObservation2[i] = followupTime;
          event2[i] = 0;
          dropoutEvent2[i] = 0;
        }
      } else {
        if (survivalTime2[i] < dropoutTime2[i]) {
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
      
      int j1 = kMaxe1x, j2 = kMax - kMaxe1x;
      
      // PFS looks
      if (kMaxe1x > 0) {
        for (j1=0; j1<kMaxe1x; j1++) {
          if (plannedEvents[j1] >= nevents1) {
            break;
          }
        }
        
        if (j1==kMaxe1x) { // total number of PFS events exceeds planned
          for (k=0; k<kMaxe1x; k++) {
            analysisTime[k] = totalt1[plannedEvents[k]-1] + 1e-12;
          }
        } else {
          for (k=0; k<=j1; k++) {
            if (k < j1) {
              analysisTime[k] = totalt1[plannedEvents[k]-1] + 1e-12;
            } else {
              analysisTime[k] = totalt1[nevents1-1] + 1e-12;
            }
          }
        }
      }
      
      // OS looks
      NumericVector analysisTime2(kMax - kMaxe1x);
      
      if (kMax > kMaxe1x) {
        for (j2=0; j2<kMax-kMaxe1x; j2++) {
          if (plannedEvents[kMaxe1x+j2] >= nevents2) {
            break;
          }
        }
        
        if (j2==kMax-kMaxe1x) { // total number of OS events exceeds planned
          for (k=0; k<kMax-kMaxe1x; k++) {
            analysisTime2[k] = totalt2[plannedEvents[kMaxe1x+k]-1] + 1e-12;
          }
        } else {
          for (k=0; k<=j2; k++) {
            if (k < j2) {
              analysisTime2[k] = totalt2[plannedEvents[kMaxe1x+k]-1] + 1e-12;
            } else {
              analysisTime2[k] = totalt2[nevents2-1] + 1e-12;
            }
          }
        }
      }
      
      // determine the number of looks and timing of the looks
      if (kMaxe1x == 0) { // all looks based on OS events
        if (j2 == kMax - kMaxe1x) {
          nstages = kMax - kMaxe1x;
        } else {
          nstages = j2 + 1;
        }
        
        for (k=0; k<nstages; k++) {
          analysisTime[k] = analysisTime2[k];
        }
      } else if (kMax == kMaxe1x) { // all looks based on PFS events
        if (j1 == kMaxe1x) {
          nstages = kMaxe1x;
        } else {
          nstages = j1 + 1;
        }
      } else {
        if (analysisTime2[kMax-kMaxe1x-1] > analysisTime[kMaxe1x-1]) {
          // only OS looks that occur after the last PFS look contribute
          int l = which_min(analysisTime2 > analysisTime[kMaxe1x-1]);
          nstages = kMax-l;
          for (k=kMaxe1x; k<kMax-l; k++) {
            analysisTime[k] = analysisTime2[k-kMaxe1x+l];
          }
        } else {
          if (j1 == kMaxe1x) {
            nstages = kMaxe1x;
          } else {
            nstages = j1 + 1;
          }
        }
      }
      
      // whether the target PFS and OS events are achieved
      if (kMaxe1x > 0 && nevents1 < plannedEvents[kMaxe1x-1]) {
        eventsNotAchieved = 1;
      } else if (nevents2 < plannedEvents[kMax-1]) {
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
    for (k=0; k<nstages; k++) {
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
      for (i=0; i<n; i++) {
        h = stratum[i]-1;
        observationTime[i] = time;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation1[i] = time - arrivalTime[i];
          event1[i] = 0;
          dropoutEvent1[i] = 0;
          
          timeUnderObservation2[i] = time - arrivalTime[i];
          event2[i] = 0;
          dropoutEvent2[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            n1x[h]++;
          } else if (treatmentGroup[i]==2) {
            n2x[h]++;
          } else {
            n3x[h]++;
          }
          
          
          // censored time for endpoint 1
          if (fixedFollowup) {
            if (arrivalTime[i] + survivalTime1[i] < time &&
                survivalTime1[i] < dropoutTime1[i] &&
                survivalTime1[i] < followupTime) {
              timeUnderObservation1[i] = survivalTime1[i];
              event1[i] = 1;
              dropoutEvent1[i] = 0;
            } else if (arrivalTime[i] + dropoutTime1[i] < time &&
              dropoutTime1[i] < survivalTime1[i] &&
              dropoutTime1[i] < followupTime) {
              timeUnderObservation1[i] = dropoutTime1[i];
              event1[i] = 0;
              dropoutEvent1[i] = 1;
            } else if (arrivalTime[i] + followupTime < time &&
              followupTime < survivalTime1[i] &&
              followupTime < dropoutTime1[i]) {
              timeUnderObservation1[i] = followupTime;
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            } else {
              timeUnderObservation1[i] = time - arrivalTime[i];
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime1[i] < time &&
                survivalTime1[i] < dropoutTime1[i]) {
              timeUnderObservation1[i] = survivalTime1[i];
              event1[i] = 1;
              dropoutEvent1[i] = 0;
            } else if (arrivalTime[i] + dropoutTime1[i] < time &&
              dropoutTime1[i] < survivalTime1[i]) {
              timeUnderObservation1[i] = dropoutTime1[i];
              event1[i] = 0;
              dropoutEvent1[i] = 1;
            } else {
              timeUnderObservation1[i] = time - arrivalTime[i];
              event1[i] = 0;
              dropoutEvent1[i] = 0;
            }
          }
          
          if (treatmentGroup[i]==1 && event1[i]) events1e1++;
          if (treatmentGroup[i]==2 && event1[i]) events2e1++;
          if (treatmentGroup[i]==3 && event1[i]) events3e1++;
          if (treatmentGroup[i]==1 && dropoutEvent1[i]) dropouts1e1++;
          if (treatmentGroup[i]==2 && dropoutEvent1[i]) dropouts2e1++;
          if (treatmentGroup[i]==3 && dropoutEvent1[i]) dropouts3e1++;
          
          
          // censored time for endpoint 2
          if (fixedFollowup) {
            if (arrivalTime[i] + survivalTime2[i] < time &&
                survivalTime2[i] < dropoutTime2[i] &&
                survivalTime2[i] < followupTime) {
              timeUnderObservation2[i] = survivalTime2[i];
              event2[i] = 1;
              dropoutEvent2[i] = 0;
            } else if (arrivalTime[i] + dropoutTime2[i] < time &&
              dropoutTime2[i] < survivalTime2[i] &&
              dropoutTime2[i] < followupTime) {
              timeUnderObservation2[i] = dropoutTime2[i];
              event2[i] = 0;
              dropoutEvent2[i] = 1;
            } else if (arrivalTime[i] + followupTime < time &&
              followupTime < survivalTime2[i] &&
              followupTime < dropoutTime2[i]) {
              timeUnderObservation2[i] = followupTime;
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            } else {
              timeUnderObservation2[i] = time - arrivalTime[i];
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            }
          } else {
            if (arrivalTime[i] + survivalTime2[i] < time &&
                survivalTime2[i] < dropoutTime2[i]) {
              timeUnderObservation2[i] = survivalTime2[i];
              event2[i] = 1;
              dropoutEvent2[i] = 0;
            } else if (arrivalTime[i] + dropoutTime2[i] < time &&
              dropoutTime2[i] < survivalTime2[i]) {
              timeUnderObservation2[i] = dropoutTime2[i];
              event2[i] = 0;
              dropoutEvent2[i] = 1;
            } else {
              timeUnderObservation2[i] = time - arrivalTime[i];
              event2[i] = 0;
              dropoutEvent2[i] = 0;
            }
          }
          
          if (treatmentGroup[i]==1 && event2[i]) events1e2++;
          if (treatmentGroup[i]==2 && event2[i]) events2e2++;
          if (treatmentGroup[i]==3 && event2[i]) events3e2++;
          if (treatmentGroup[i]==1 && dropoutEvent2[i]) dropouts1e2++;
          if (treatmentGroup[i]==2 && dropoutEvent2[i]) dropouts2e2++;
          if (treatmentGroup[i]==3 && dropoutEvent2[i]) dropouts3e2++;
        }
      }
      
      
      // add raw data to output
      if (niter[k] < maxNumberOfRawDatasetsPerStage) {
        for (i=0; i<n; i++) {
          iterationNumberx[index1] = iter+1;
          stageNumberx[index1] = k+1;
          subjectIdx[index1] = i+1;
          arrivalTimex[index1] = arrivalTime[i];
          stratumx[index1] = stratum[i];
          treatmentGroupx[index1] = treatmentGroup[i];
          observationTimex[index1] = observationTime[i];
          
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
          
          index1++;
        }
        
        // update the number of stage k dataset to extract
        niter[k]++;
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
      
      
      for (int endpoint=1; endpoint<=2; endpoint++) {
        n1 = clone(n1x);
        n2 = clone(n2x);
        n3 = clone(n3x);
        
        double hazardRatioH013, hazardRatioH023, hazardRatioH012;
        if (endpoint == 1) {
          hazardRatioH013 = hazardRatioH013e1;
          hazardRatioH023 = hazardRatioH023e1;
          hazardRatioH012 = hazardRatioH012e1;
        } else {
          hazardRatioH013 = hazardRatioH013e2;
          hazardRatioH023 = hazardRatioH023e2;
          hazardRatioH012 = hazardRatioH012e2;
        }
        
        
        // order the data by time under observation
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
        
        stratumSorted = stratum[sortedIndex];
        treatmentGroupSorted = treatmentGroup[sortedIndex];
        sub = (timeUnderObservationSorted > 0);
        eventSorted = eventSorted[sub];
        stratumSorted = stratumSorted[sub];
        treatmentGroupSorted = treatmentGroupSorted[sub];
        nsub = eventSorted.size();
        
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
        for (i=0; i<nsub; i++) {
          h = stratumSorted[i] - 1;
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
                                  - n13a[h]/(nt13a[h]+0.0));
            vscore13 += w13[h]*w13[h]*n13a[h]*n3[h]/(nt13a[h]*nt13a[h]+0.0);
            km13[h] *= (1-1/(nt13[h]+0.0)); // update km estimate
          }
          
          if (eventSorted[i] && (treatmentGroupSorted[i]==2 || 
              treatmentGroupSorted[i]==3)) { 
            w23[h] = pow(km23[h], rho1)*pow(1-km23[h], rho2);
            uscore23 += w23[h]*((treatmentGroupSorted[i]==2) 
                                  - n23a[h]/(nt23a[h]+0.0));
            vscore23 += w23[h]*w23[h]*n23a[h]*n3[h]/(nt23a[h]*nt23a[h]+0.0);
            km23[h] *= (1-1/(nt23[h]+0.0)); // update km estimate
          }
          
          if (eventSorted[i] && (treatmentGroupSorted[i]==1 || 
              treatmentGroupSorted[i]==2)) { 
            w12[h] = pow(km12[h], rho1)*pow(1-km12[h], rho2);
            uscore12 += w12[h]*((treatmentGroupSorted[i]==1) 
                                  - n12a[h]/(nt12a[h]+0.0));
            vscore12 += w12[h]*w12[h]*n12a[h]*n2[h]/(nt12a[h]*nt12a[h]+0.0);
            km12[h] *= (1-1/(nt12[h]+0.0)); // update km estimate
          }
          
          // reduce the risk set
          if (treatmentGroupSorted[i]==1) {
            n1[h]--;
          } else if (treatmentGroupSorted[i]==2) {
            n2[h]--;
          } else {
            n3[h]--;
          }
        }
        
        
        
        // add summary data to output  
        iterationNumbery[index2] = iter+1;
        eventsNotAchievedy[index2] = eventsNotAchieved;
        stageNumbery[index2] = k+1;
        analysisTimey[index2] = analysisTime[k];
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
        
        index2++;
        
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
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    treatmentGroupx = treatmentGroupx[sub1];
    observationTimex = observationTimex[sub1];
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
      _["analysisTime"] = observationTimex,
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




