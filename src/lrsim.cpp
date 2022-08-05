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
//'   Same as informationRates in terms of number of events for
//'   the conventional log-rank test. Use \code{caltime} and \code{lrstat}
//'   to derive the information time for weighted log-rank tests.
//'   Fixed prior to the trial. Defaults to
//'   \code{plannedEvents / sum(plannedEvents)} if left unspecified.
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
//' @param plannedTime The planned analysis time for each stage needed for  
//'   analyses planned at calendar times, in which case, plannedEvents should 
//'   be missing.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per stage
//'   to extract. Defaults to 1.
//' @param seed The seed to reproduce the simulation results.
//'   The computer clock will be used if left unspecified,
//'
//' @return A list of S3 class \code{lrsim} with 3 components:
//'
//' * \code{overview} is a list of containing incremental and cumulative
//' efficacy and futility stopping probabilities by stage, expected number
//' of events, number of dropouts, number of subjects, and analysis time
//' by stage, overall rejection probability, overall expected number of
//' events, number of dropouts, number of subjects, and study duration,
//' the hazard ratio under H0, and whether the analyses are planned 
//' based on the number of events or calendar time.
//'
//' * \code{sumdata} is a data frame of summary data by stage for each
//' iteration, containing at which stage the trial stops, whether the target
//' number of events is achieved, the analysis time, number of accrued
//' subjects overall and by treatment group, number of events overall and by
//' treatment group, number of dropouts overall and by treatment group,
//' numerator and variance of weighted log-rank score statistic, log-rank
//' test Z-statistic, and whether the trial stops for efficacy or futility
//' at the stage.
//'
//'
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//' positive integer) is a data frame for subject-level data for selected
//' replications, containing the subject number, arrival time, stratum,
//' treatment group, survival time, dropout time, observation time when
//' the trial stops, time under observation, and event and dropout
//' indicators.
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
  
  NumericVector informationTime1 = clone(informationTime);
  NumericVector futilityBounds1 = clone(futilityBounds);
  
  int i, iter, j, j1, j2, k, h, nsub;
  int nstrata = stratumFraction.size();
  int nintervals = piecewiseSurvivalTime.size();
  int nsi = nstrata*nintervals;
  
  NumericVector lambda1x(nsi), lambda2x(nsi), gamma1x(nsi), gamma2x(nsi);
  
  double u, enrollt;
  bool useEvents;
  
  // b1 and b2 are the available slots for the two treatments in a block
  IntegerVector b1(nstrata);
  IntegerVector b2(nstrata);
  
  if (R_isnancpp(kMax)) {
    stop("kMax must be provided");
  }
  
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
    
    if (is_true(any(diff(plannedEvents) <= 0))) {
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
    
    if (is_true(any(diff(plannedTime) <= 0))) {
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
  
  
  if (kMax > 1) {
    if (is_true(any(is_na(futilityBounds)))) {
      futilityBounds1 = rep(-6.0, kMax-1);
    }
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
  
  
  IntegerVector stratum(n), treatmentGroup(n);
  
  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n),
  observationTime(n), timeUnderObservation(n), totalTime(n), totalt(n);
  
  LogicalVector event(n), dropoutEvent(n);
  
  NumericVector cumStratumFraction = cumsum(stratumFraction);
  
  NumericVector lam1(nintervals), lam2(nintervals);
  NumericVector gam1(nintervals), gam2(nintervals);
  
  int nevents, nstages;
  NumericVector analysisTime(kMax);
  
  IntegerVector accruals1(kMax), accruals2(kMax), totalAccruals(kMax),
  events1(kMax), events2(kMax), totalEvents(kMax),
  dropouts1(kMax), dropouts2(kMax), totalDropouts(kMax);
  
  NumericVector timeUnderObservationSorted(n);
  IntegerVector sortedIndex(n), stratumSorted(n), treatmentGroupSorted(n);
  LogicalVector eventSorted(n);
  
  double uscore1, vscore1;
  NumericVector km1(nstrata), w1(nstrata);
  
  NumericVector uscore(kMax), vscore(kMax), lrstat(kMax);
  
  LogicalVector rejectPerStage(kMax), futilityPerStage(kMax);
  int stopStage;
  
  LogicalVector sub(n);
  
  IntegerVector n1(nstrata), n1x(nstrata), n2(nstrata),
  nt(nstrata), ntx(nstrata);
  
  // cache for the number of raw data sets per stage to extract
  IntegerVector niter(kMax);
  
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
  
  
  
  // observed number of events by stage and adjusted critical values
  IntegerVector obsEvents(kMax);
  NumericVector adjCriticalValues(kMax);
  
  NumericVector lb = rep(-6.0, kMax);
  NumericVector theta = rep(0.0, kMax);
  List p1 = exitprob(criticalValues, lb, theta, informationTime1);
  double alpha = sum(NumericVector(p1[0]));
  
  int index1=0, index2=0;
  double time;
  bool eventsNotAchieved;
  
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
      enrollt = qtpwexp(u, accrualTime, accrualIntensity, enrollt);
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
      
      // generate survival time
      j1 = j*nintervals;
      j2 = j1 + nintervals - 1;
      
      lam1 = lambda1x[Range(j1,j2)];
      lam2 = lambda2x[Range(j1,j2)];
      gam1 = gamma1x[Range(j1,j2)];
      gam2 = gamma2x[Range(j1,j2)];
      
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        survivalTime[i] = qtpwexp(u, piecewiseSurvivalTime, lam1, 0);
      } else {
        survivalTime[i] = qtpwexp(u, piecewiseSurvivalTime, lam2, 0);
      }
      
      // generate dropout time
      u = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime[i] = qtpwexp(u, piecewiseSurvivalTime, gam1, 0);
      } else {
        dropoutTime[i] = qtpwexp(u, piecewiseSurvivalTime, gam2, 0);
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
      eventsNotAchieved = (obsEvents[nstages-1] < plannedEvents[kMax-1]);
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
              timeUnderObservation[i] = std::min(time - arrivalTime[i],
                                                 followupTime);
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
      km1.fill(1);  // km(t-) estimate by stratum
      for (i=0; i<nsub; i++) {
        h = stratumSorted[i] - 1;
        n1x[h] = n1[h]*hazardRatioH0;
        nt[h] = n1[h] + n2[h];
        ntx[h] = n1x[h] + n2[h];
        
        if (eventSorted[i]) { // at most 1 event can occur at a given time
          w1[h] = pow(km1[h], rho1)*pow(1-km1[h], rho2);
          uscore1 += w1[h]*((treatmentGroupSorted[i]==1)-n1x[h]/(ntx[h]+0.0));
          vscore1 += w1[h]*w1[h]*n1x[h]*n2[h]/(ntx[h]*ntx[h]+0.0);
          km1[h] *= (1-1/(nt[h]+0.0)); // update km estimate
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
      lrstat[k] = uscore[k]/sqrt(vscore[k]);
      
      if (useEvents) {
        // adjust the critical value at the final stage if the planned total
        // number of events is not achieved
        if (k < nstages-1 || !eventsNotAchieved) { // no change to the critical
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
                        List p2 = exitprob(u, l, theta, I);
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
                        List p2 = exitprob(u, l, theta, I);
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
        
        if (nstops == 1) {
          
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
    _["rejectPerStage"]=pRejectPerStage,
    _["futilityPerStage"]=pFutilityPerStage,
    _["cumulativeRejection"]=cpu,
    _["cumulativeFutility"]=cpl,
    _["numberOfEvents"]=nEventsPerStage,
    _["numberOfDropouts"]=nDropoutsPerStage,
    _["numberOfSubjects"]=nSubjectsPerStage,
    _["analysisTime"]=analysisTimePerStage,
    _["overallReject"]=pOverallReject,
    _["expectedNumberOfEvents"]=expectedNumberOfEvents,
    _["expectedNumberOfDropouts"]=expectedNumberOfDropouts,
    _["expectedNumberOfSubjects"]=expectedNumberOfSubjects,
    _["expectedStudyDuration"]=expectedStudyDuration,
    _["hazardRatioH0"]=hazardRatioH0,
    _["useEvents"]=useEvents);
  
  
  
  // simulation datasets
  DataFrame sumdata = DataFrame::create(
    _["iterationNumber"]=iterationNumbery,
    _["stageNumber"]=stageNumbery,
    _["stopStage"]=stopStagey,
    _["eventsNotAchieved"]=eventsNotAchievedy,
    _["analysisTime"]=analysisTimey,
    _["accruals1"]=accruals1y,
    _["accruals2"]=accruals2y,
    _["totalAccruals"]=totalAccrualsy,
    _["events1"]=events1y,
    _["events2"]=events2y,
    _["totalEvents"]=totalEventsy,
    _["dropouts1"]=dropouts1y,
    _["dropouts2"]=dropouts2y,
    _["totalDropouts"]=totalDropoutsy,
    _["uscore"]=uscorey,
    _["vscore"]=vscorey,
    _["logRankStatistic"]=logRankStatisticy,
    _["rejectPerStage"]=rejectPerStagey,
    _["futilityPerStage"]=futilityPerStagey);
  
  
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
      _["iterationNumber"]=iterationNumberx,
      _["stopStage"]=stopStagex,
      _["subjectId"]=subjectIdx,
      _["arrivalTime"]=arrivalTimex,
      _["stratum"]=stratumx,
      _["treatmentGroup"]=treatmentGroupx,
      _["survivalTime"]=survivalTimex,
      _["dropoutTime"]=dropoutTimex,
      _["observationTime"]=observationTimex,
      _["timeUnderObservation"]=timeUnderObservationx,
      _["event"]=eventx,
      _["dropoutEvent"]=dropoutEventx);
    
    result = List::create(_["overview"]=overview,
                          _["sumdata"]=sumdata,
                          _["rawdata"]=rawdata);
  } else {
    result = List::create(_["overview"]=overview,
                          _["sumdata"]=sumdata);
  }
  
  result.attr("class") = "lrsim";
  
  
  return result;
}



