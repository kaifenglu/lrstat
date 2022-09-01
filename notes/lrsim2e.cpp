#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;


//' @title Log-rank test simulation for three arms and two endpoints
//' @description Performs simulation for three-arm two-endpoint group 
//' sequential trials based on weighted log-rank test. The first kMaxe1 
//' looks are driven by the total number of PFS events in Arm A 
//' and Arm C combined, and the subsequent looks are driven by the total 
//' number of OS events in Arm A and Arm C combined. 
//' 
//' @inheritParams param_kMax
//' @param kMaxe1 Number of stages for PFS comparison of Arm A vs. Arm C.
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
//' @param gamma1e1 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 1 and endpoint 1 (PFS).
//' @param gamma2e1 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 2 and endpoint 1 (PFS).
//' @param gamma3e1 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 3 and endpoint 1 (PFS).
//' @param gamma1e2 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 1 and endpoint 2 (OS).
//' @param gamma2e2 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 2 and endpoint 2 (OS).
//' @param gamma3e2 The hazard rate for exponential dropout, a vector of 
//'   hazard rates for piecewise exponential dropout applicable for all 
//'   strata, or a vector of hazard rates for dropout in each analysis time 
//'   interval by stratum for arm 3 and endpoint 2 (OS).
//' @inheritParams param_accrualDuration
//' @param plannedEvents The planned cumulative total number of PFS events at 
//'   Look 1 to Look kMaxe1 for Arms A and C combined and the planned 
//'   cumulative total number of OS events at Look kMaxe1+1 to Look kMax 
//'   for Arms A and C combined.
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
//' @return A list with 2 components:
//'
//' * \code{sumdata} is a data frame of summary data by stage for each
//' iteration, containing the analysis time, number of accrued subjects 
//' overall and by treatment group, and number of events overall and 
//' by treatment group, number of dropouts overall and by treatment group, 
//' and log-rank test statistic for each comparison by endpoint.
//'
//' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
//' positive integer) is a data frame for subject-level data for selected
//' replications, containing the stage number, subject number, arrival time, 
//' stratum, treatment group, observation time, and survival time, 
//' dropout time, time under observation, event and dropout indicators 
//' for each endpoint.
//' 
//' 
//' @export
// [[Rcpp::export]]
List lrsim2e(const int kMax = NA_INTEGER,
             const int kMaxe1 = NA_INTEGER,
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
  
  NumericVector lambda1e1x(nsi), lambda2e1x(nsi), lambda3e1x(nsi), 
  lambda1e2x(nsi), lambda2e2x(nsi), lambda3e2x(nsi), 
  gamma1e1x(nsi), gamma2e1x(nsi), gamma3e1x(nsi), 
  gamma1e2x(nsi), gamma2e2x(nsi), gamma3e2x(nsi),
  lambda1e1d(nsi), lambda2e1d(nsi), lambda3e1d(nsi), 
  gamma1e1d(nsi), gamma2e1d(nsi), gamma3e1d(nsi);

  bool useEvents;
  
  if (R_isnancpp(kMax)) {
    stop("kMax must be provided");
  }
  
  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }
  
  if (R_isnancpp(kMaxe1)){
    kMaxe1x = kMax;
  }
  
  if (kMaxe1x < 1) {
    stop("kMaxe1 must be a positive integer");
  }
  
  if (kMaxe1x > kMax) {
    stop("kMaxe1 must be less than or equal to kMax");
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
  
  // whether to plan the analyses based on events or calender time
  if (is_false(any(is_na(plannedEvents)))) {
    useEvents = 1;
    
    if (plannedEvents.size() != kMax) {
      stop("Invalid length for plannedEvents");
    }
    
    if (plannedEvents[0] <= 0) {
      stop("Elements of plannedEvents must be positive");
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
    
    if (is_true(any(diff(plannedTime) <= 0))) {
      stop("Elements of plannedTime must be increasing");
    }
  } else {
    stop("Either plannedEvents or plannedTime must be given");
  }
  
  
  if (maxNumberOfIterations < 1) {
    stop("maxNumberOfIterations must be a positive integer");
  }
  
  if (maxNumberOfRawDatasetsPerStage < 0) {
    stop("maxNumberOfRawDatasetsPerStage must be a non-negative integer");
  }
  
  
  // declare variables 
  int i, m, n, iter, j, k, nstages, h, index1=0, index2=0, nsub, 
  accruals1, accruals2, accruals3, totalAccruals, 
  events1e1, events2e1, events3e1, totalEventse1, 
  events1e2, events2e2, events3e2, totalEventse2, 
  dropouts1e1, dropouts2e1, dropouts3e1, totalDropoutse1, 
  dropouts1e2, dropouts2e2, dropouts3e2, totalDropoutse2;
  
  double s, enrollt, u, u1, u2, time, uscore13, uscore23, uscore12, 
  vscore13, vscore23, vscore12;
  
  
  // stratum information 
  IntegerVector b1(nstrata), b2(nstrata), b3(nstrata), 
  n1(nstrata), n2(nstrata), n3(nstrata), 
  n1x(nstrata), n2x(nstrata), n3x(nstrata), 
  nt13(nstrata), nt23(nstrata), nt12(nstrata);
  
  NumericVector cumStratumFraction = cumsum(stratumFraction);

  
  // maximum number of subjects to enroll
  m = accrualTime.size();
  s = 0;
  for (i=0; i<m; i++) {
    if (i<m-1 && accrualTime[i+1] < accrualDuration) {
      s += accrualIntensity[i]*(accrualTime[i+1] - accrualTime[i]);
    } else {
      s += accrualIntensity[i]*(accrualDuration - accrualTime[i]);
      break;
    }
  }
  n = floor(s + 0.5);
  
  
  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n), sortedIndex(n), 
  stratumSorted(n), treatmentGroupSorted(n);
  
  NumericVector arrivalTime(n), survivalTime1(n), survivalTime2(n), 
  dropoutTime1(n), dropoutTime2(n), timeUnderObservation1(n), 
  timeUnderObservation2(n), totalTime1(n), totalTime2(n), 
  totalt1(n), totalt2(n), observationTime(n), timeUnderObservationSorted(n);
  
  LogicalVector event1(n), event2(n), dropoutEvent1(n), dropoutEvent2(n), 
  event1ac(n), event2ac(n), sub(n), eventSorted(n);
  

  // within-stratum hazard rates
  IntegerVector jj(nintervals);  
  
  NumericVector lam1e1(nintervals), lam2e1(nintervals), lam3e1(nintervals), 
  lam1e2(nintervals), lam2e2(nintervals), lam3e2(nintervals), 
  gam1e1(nintervals), gam2e1(nintervals), gam3e1(nintervals),
  gam1e2(nintervals), gam2e2(nintervals), gam3e2(nintervals);
  
  
  // stage-wise information
  IntegerVector niter(kMax);
  NumericVector analysisTime(kMax);
 
  
  // cache for the number of raw data sets per stage to extract
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
  
  
  // cache for the summary data sets to extract
  int nrow2 = kMax*maxNumberOfIterations*2;
  
  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
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
      
      // generate survival time
      jj = Range(j*nintervals, (j+1)*nintervals-1);

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
        survivalTime1[i] = qtpwexp(u1, piecewiseSurvivalTime, lam1e1, 0);
        survivalTime2[i] = qtpwexp(u2, piecewiseSurvivalTime, lam1e2, 0);
      } else if (treatmentGroup[i]==2) {
        survivalTime1[i] = qtpwexp(u1, piecewiseSurvivalTime, lam2e1, 0);
        survivalTime2[i] = qtpwexp(u2, piecewiseSurvivalTime, lam2e2, 0);
      } else {
        survivalTime1[i] = qtpwexp(u1, piecewiseSurvivalTime, lam3e1, 0);
        survivalTime2[i] = qtpwexp(u2, piecewiseSurvivalTime, lam3e2, 0);
      }
      // PFS includes death
      survivalTime1[i] = std::min(survivalTime1[i], survivalTime2[i]);
      
      
      // generate dropout times
      u1 = R::runif(0,1);
      u2 = R::runif(0,1);
      if (treatmentGroup[i]==1) {
        dropoutTime1[i] = qtpwexp(u1, piecewiseSurvivalTime, gam1e1, 0);
        dropoutTime2[i] = qtpwexp(u2, piecewiseSurvivalTime, gam1e2, 0);
      } else if (treatmentGroup[i]==2) {
        dropoutTime1[i] = qtpwexp(u1, piecewiseSurvivalTime, gam2e1, 0);
        dropoutTime2[i] = qtpwexp(u2, piecewiseSurvivalTime, gam2e2, 0);
      } else {
        dropoutTime1[i] = qtpwexp(u1, piecewiseSurvivalTime, gam3e1, 0);
        dropoutTime2[i] = qtpwexp(u2, piecewiseSurvivalTime, gam3e2, 0);
      }
      // whatever censors OS will also censors PFS
      dropoutTime1[i] = std::min(dropoutTime1[i], dropoutTime2[i]);
      
      
      // initial observed time and event indicator
      if (survivalTime1[i] < dropoutTime1[i]) {
        timeUnderObservation1[i] = survivalTime1[i];
        event1[i] = 1;
        dropoutEvent1[i] = 0;
      } else {
        timeUnderObservation1[i] = dropoutTime1[i];
        event1[i] = 0;
        dropoutEvent1[i] = 1;
      }
      totalTime1[i] = arrivalTime[i] + timeUnderObservation1[i];
      
      if (survivalTime2[i] < dropoutTime2[i]) {
        timeUnderObservation2[i] = survivalTime2[i];
        event2[i] = 1;
        dropoutEvent2[i] = 0;
      } else {
        timeUnderObservation2[i] = dropoutTime2[i];
        event2[i] = 0;
        dropoutEvent2[i] = 1;
      }
      totalTime2[i] = arrivalTime[i] + timeUnderObservation2[i];
      
    }
    
    
    // find the analysis time for each stage based on Arm A vs. Arm C
    event1ac = event1 & ((treatmentGroup==1) | (treatmentGroup==3));
    totalt1 = stl_sort(totalTime1[event1ac]);
    
    event2ac = event2 & ((treatmentGroup==1) | (treatmentGroup==3));
    totalt2 = stl_sort(totalTime2[event2ac]);
    

    // find the analysis time for each stage
    if (useEvents) {
      // PFS looks
      for (k=0; k<kMaxe1x; k++) {
        analysisTime[k] = totalt1[plannedEvents[k]-1] + 1e-12;
      }
      
      if (kMax > kMaxe1x) {
        // planned OS looks after the PFS looks
        NumericVector analysisTime2(kMax - kMaxe1x);
        for (k=kMaxe1x; k<kMax; k++) {
          analysisTime2[k-kMaxe1x] = totalt2[plannedEvents[k]-1] + 1e-12;
        }
        
        // whether the planned OS looks occur after the PFS looks
        if (analysisTime2[kMax-kMaxe1x-1] > analysisTime[kMaxe1x-1]) {
          int l = which_min(analysisTime2 > analysisTime[kMaxe1x-1]);
          for (k=kMaxe1x; k<kMax-l; k++) {
            analysisTime[k] = analysisTime2[k-kMaxe1x+l];
          }
          nstages = kMax-l;
        } else {
          nstages = kMaxe1x;
        }
      } else {
        nstages = kMax;
      }
    } else {
      nstages = kMax;
      analysisTime = clone(plannedTime);
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
          
          if (treatmentGroup[i]==1 && event1[i]) events1e1++;
          if (treatmentGroup[i]==2 && event1[i]) events2e1++;
          if (treatmentGroup[i]==3 && event1[i]) events3e1++;
          if (treatmentGroup[i]==1 && dropoutEvent1[i]) dropouts1e1++;
          if (treatmentGroup[i]==2 && dropoutEvent1[i]) dropouts2e1++;
          if (treatmentGroup[i]==3 && dropoutEvent1[i]) dropouts3e1++;
          
          
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
        for (i=0; i<nsub; i++) {
          h = stratumSorted[i] - 1;
          nt13[h] = n1[h] + n3[h];
          nt23[h] = n2[h] + n3[h];
          nt12[h] = n1[h] + n2[h];
          
          if (eventSorted[i] && (treatmentGroupSorted[i]==1 || 
              treatmentGroupSorted[i]==3)) { 
            uscore13 += (treatmentGroupSorted[i]==1) - n1[h]/(nt13[h]+0.0);
            vscore13 += n1[h]*n3[h]/(nt13[h]*nt13[h]+0.0);
          }
          
          if (eventSorted[i] && (treatmentGroupSorted[i]==2 || 
              treatmentGroupSorted[i]==3)) { 
            uscore23 += (treatmentGroupSorted[i]==2) - n2[h]/(nt23[h]+0.0);
            vscore23 += n2[h]*n3[h]/(nt23[h]*nt23[h]+0.0);
          }
          
          if (eventSorted[i] && (treatmentGroupSorted[i]==1 || 
              treatmentGroupSorted[i]==2)) { 
            uscore12 += (treatmentGroupSorted[i]==1) - n1[h]/(nt12[h]+0.0);
            vscore12 += n1[h]*n2[h]/(nt12[h]*nt12[h]+0.0);
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
      _["subjectId"] = subjectIdx,
      _["arrivalTime"] = arrivalTimex,
      _["stratum"] = stratumx,
      _["treatmentGroup"] = treatmentGroupx,
      _["observationTime"] = observationTimex,
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
    
    result = List::create(_["sumdata"]=sumdata,
                          _["rawdata"]=rawdata);
  } else {
    result = List::create(_["sumdata"]=sumdata);
  }
  
  return result;
}



