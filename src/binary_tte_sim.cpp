#include "utilities.h"
using namespace Rcpp;


//' @title Simulation for a binary endpoint and a time-to-event endpoint
//' @description Performs simulation for two-endpoint two-arm group
//' sequential trials. The first endpoint is a binary endpoint and
//' the Mantel-Haenszel test is used to test risk difference.
//' The second endpoint is a time-to-event endpoint and the log-rank
//' test is used to test the treatment difference. The analysis times
//' of the first endpoint are determined by the specified calendar times,
//' while the analysis times for the second endpoint is based on the
//' planned number of events at each look. The binary endpoint is
//' assessed at the first post-treatment follow-up visit.
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
//' @param globalOddsRatio The global odds ratio of the Plackett copula.
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
//' @param upper1 The protocol-specified treatment duration for the treatment
//'   group.
//' @param upper2 The protocol-specified treatment duration for the control
//'   group.
//' @inheritParams param_accrualDuration
//' @param plannedTime The calendar times for the analyses of the binary
//'   endpoint.
//' @param plannedEvents The planned cumulative total number of events for
//'   the time-to-event endpoint.
//' @param maxNumberOfIterations The number of simulation iterations.
//' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
//'   stage to extract.
//' @param seed The seed to reproduce the simulation results.
//'   The seed from the environment will be used if left unspecified,
//'
//' @details We consider dual primary endpoints with endpoint 1 being a
//'   binary endpoint and endpoint 2 being a time-to-event endpoint.
//'   The analyses of endpoint 1 will be based on calendar times, while
//'   the analyses of endpoint 2 will be based on the number of events.
//'   Therefor the analyses of the two endpoints are not at the same
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
//'   accrualDuration = 20,
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
    const double accrualDuration = NA_REAL,
    const NumericVector& plannedTime = NA_REAL,
    const IntegerVector& plannedEvents = NA_INTEGER,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasetsPerStage = 0,
    const int seed = NA_INTEGER) {

  // check input parameters
  int nstrata = static_cast<int>(stratumFraction.size());
  int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
  int nsi = nstrata*nintervals;

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

  if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
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

  if (delta1.size() == 1) {
    delta1x = rep(delta1, nsi);
  } else if (delta1.size() == nintervals) {
    delta1x = rep(delta1, nstrata);
  } else if (delta1.size() == nsi) {
    delta1x = delta1;
  } else {
    stop("Invalid length for delta1");
  }

  if (delta2.size() == 1) {
    delta2x = rep(delta2, nsi);
  } else if (delta2.size() == nintervals) {
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


  if (R_isnancpp(accrualDuration)) {
    stop("accrualDuration must be provided");
  }

  if (accrualDuration <= 0) {
    stop("accrualDuration must be positive");
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

  if (plannedTime[kMax1-1] != accrualDuration + std::max(upper1, upper2)) {
    stop("plannedTime must end at the maximum follow-up time for endpoint 1");
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
  int h, i, iter, j, k, nevents, nstages1, nstages2, nsub;
  int accruals1, accruals2, totalAccruals;
  int events1, events2, totalEvents;
  int dropouts1, dropouts2, totalDropouts;
  int index1x=0, index2x=0, index1y=0, index2y=0;

  double enrollt, u, u1, u2, time, uscore, vscore;


  // maximum number of subjects to enroll
  int m = static_cast<int>(accrualTime.size());
  double s = 0;
  for (i=0; i<m; i++) {
    if (i<m-1 && accrualTime[i+1] < accrualDuration) {
      s += accrualIntensity[i]*(accrualTime[i+1] - accrualTime[i]);
    } else {
      s += accrualIntensity[i]*(accrualDuration - accrualTime[i]);
      break;
    }
  }
  int n = static_cast<int>(floor(s + 0.5));


  // subject-level raw data set for one simulation
  IntegerVector stratum(n), treatmentGroup(n);
  IntegerVector order(n), stratumSorted(n), treatmentGroupSorted(n);

  NumericVector arrivalTime(n), survivalTime(n), dropoutTime(n);
  NumericVector timeUnderObservation(n), totalTime(n), totalt(n);
  NumericVector observationTime(n), timeUnderObservationSorted(n);

  LogicalVector event(n), dropoutEvent(n), sub(n), eventSorted(n);

  NumericVector latentResponse(n);
  NumericVector trtDiscTime(n), ptfu1Time(n);
  NumericVector observationTime1(n), timeUnderObservation1(n);
  LogicalVector responder(n);
  IntegerVector source(n);

  // stratum information
  IntegerVector b1(nstrata), b2(nstrata);
  IntegerVector n1(nstrata), n2(nstrata), nt(nstrata);

  // original copy of n1 and n2 when looping over the endpoints
  IntegerVector n1x(nstrata), n2x(nstrata);

  // hazardRatioH0 adjusted n1 and nt for calculating the log-rank statistic
  NumericVector n1a(nstrata), nta(nstrata);

  NumericVector cumStratumFraction = cumsum(stratumFraction);

  // within-stratum response rates and hazard rates
  double alpha0, alpha1; // parameters on the logit scale for endpoint 1
  NumericVector lam1(nintervals), lam2(nintervals);
  NumericVector gam1(nintervals), gam2(nintervals);
  NumericVector del1(nintervals), del2(nintervals);


  // stage-wise information
  NumericVector analysisTime1(kMax1), analysisTime2(kMax2);
  IntegerVector niter1(kMax1), niter2(kMax2);


  // cache for the patient-level raw data for endpoint 1 (uMRD) to extract
  int nrow1x = n*kMax1*maxNumberOfRawDatasetsPerStage;

  IntegerVector iterationNumber1x = IntegerVector(nrow1x, NA_INTEGER);
  IntegerVector stageNumber1x(nrow1x);
  IntegerVector subjectId1x(nrow1x);
  NumericVector arrivalTime1x(nrow1x);
  IntegerVector stratum1x(nrow1x);
  IntegerVector treatmentGroup1x(nrow1x);
  NumericVector survivalTime1x(nrow1x);
  NumericVector dropoutTime1x(nrow1x);
  NumericVector ptfu1Timex(nrow1x);
  NumericVector observationTime1x(nrow1x);
  NumericVector timeUnderObservation1x(nrow1x);
  LogicalVector responderx(nrow1x);
  IntegerVector sourcex(nrow1x);


  // cache for the patient-level raw data for endpoint 2 (PFS) to extract
  int nrow2x = n*kMax2*maxNumberOfRawDatasetsPerStage;

  IntegerVector iterationNumber2x = IntegerVector(nrow2x, NA_INTEGER);
  IntegerVector stageNumber2x(nrow2x);
  IntegerVector subjectId2x(nrow2x);
  NumericVector arrivalTime2x(nrow2x);
  IntegerVector stratum2x(nrow2x);
  IntegerVector treatmentGroup2x(nrow2x);
  NumericVector survivalTime2x(nrow2x);
  NumericVector dropoutTime2x(nrow2x);
  NumericVector observationTime2x(nrow2x);
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

  for (iter=0; iter<maxNumberOfIterations; iter++) {

    b1.fill(allocation1);
    b2.fill(allocation2);

    enrollt = 0;
    for (i=0; i<n; i++) {

      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);
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

      // stratum-specific response rates and hazard rates
      alpha0 = alpha0x[j];
      alpha1 = alpha1x[j];

      Range jj = Range(j*nintervals, (j+1)*nintervals-1);

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
      if (survivalTime[i] < dropoutTime[i]) {
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
    for (k=0; k<nstages1; k++) {
      time = analysisTime1[k];

      n1x.fill(0);  // number of subjects in each stratum by treatment
      n2x.fill(0);

      for (i=0; i<n; i++) {
        h = stratum[i]-1;
        observationTime1[i] = time;
        if (arrivalTime[i] > time) { // patients not yet enrolled
          timeUnderObservation1[i] = time - arrivalTime[i];
          responder[i] = NA_LOGICAL;
          source[i] = 0;
        } else {
          if (treatmentGroup[i]==1) {
            n1x[h]++;
          } else {
            n2x[h]++;
          }

          // censor at analysis time
          if (arrivalTime[i] + ptfu1Time[i] < time &&
              ptfu1Time[i] < survivalTime[i] &&
              ptfu1Time[i] < dropoutTime[i]) {
            timeUnderObservation1[i] = ptfu1Time[i];
            responder[i] = (latentResponse[i] <= 0);
            source[i] = 1;
          } else if (arrivalTime[i] + survivalTime[i] < time &&
            survivalTime[i] < dropoutTime[i] &&
            survivalTime[i] < ptfu1Time[i]) {
            timeUnderObservation1[i] = survivalTime[i];
            responder[i] = 0;
            source[i] = 2;
          } else if (arrivalTime[i] + dropoutTime[i] < time &&
            dropoutTime[i] < survivalTime[i] &&
            dropoutTime[i] < ptfu1Time[i]) {
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
        for (i=0; i<n; i++) {
          iterationNumber1x[index1x] = iter+1;
          stageNumber1x[index1x] = k+1;
          subjectId1x[index1x] = i+1;
          arrivalTime1x[index1x] = arrivalTime[i];
          stratum1x[index1x] = stratum[i];
          treatmentGroup1x[index1x] = treatmentGroup[i];
          survivalTime1x[index1x] = survivalTime[i];
          dropoutTime1x[index1x] = dropoutTime[i];
          ptfu1Timex[index1x] = ptfu1Time[i];
          observationTime1x[index1x] = observationTime1[i];
          timeUnderObservation1x[index1x] = timeUnderObservation1[i];
          responderx[index1x] = responder[i];
          sourcex[index1x] = source[i];

          index1x++;
        }

        // update the number of stage k dataset to extract
        niter1[k]++;
      }


      accruals1 = sum(n1x);
      accruals2 = sum(n2x);
      totalAccruals = accruals1 + accruals2;




      sub = !is_na(responder); // exclude subjects administratively censored
      eventSorted = responder[sub];
      stratumSorted = stratum[sub];
      treatmentGroupSorted = treatmentGroup[sub];
      nsub = static_cast<int>(eventSorted.size());

      // obtain the Mantel-Haenszel statistic for stratified risk difference
      NumericVector n11(nstrata), n21(nstrata), n1s(nstrata), n2s(nstrata);
      for (i=0; i<nsub; i++) {
        h = stratumSorted[i] - 1;
        if (treatmentGroupSorted[i] == 1) {
          n1s[h]++;
          n11[h] += eventSorted[i];
        } else {
          n2s[h]++;
          n21[h] += eventSorted[i];
        }
      }
      NumericVector nss = n1s + n2s;

      double A = 0, B = 0, P = 0, Q = 0;
      for (h=0; h<nstrata; h++) {
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
      index1y++;
    } // end of stage for endpoint 1


    // find the analysis time for each stage for endpoint 2
    nevents = sum(event);
    totalt = stl_sort(totalTime[event]);
    nstages2 = kMax2;

    for (j=0; j<kMax2; j++) {
      if (plannedEvents[j] >= nevents) {
        nstages2 = j+1;
        break;
      }
    }

    if (j==kMax2) { // total number of events exceeds planned
      for (k=0; k<nstages2; k++) {
        analysisTime2[k] = totalt[plannedEvents[k]-1] + 1e-12;
      }
    } else {
      for (k=0; k<nstages2; k++) {
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
    for (k=0; k<nstages2; k++) {
      time = analysisTime2[k];

      n1x.fill(0);  // number of subjects in each stratum by treatment
      n2x.fill(0);

      events1 = 0;
      events2 = 0;

      dropouts1 = 0;
      dropouts2 = 0;

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
            n1x[h]++;
          } else {
            n2x[h]++;
          }

          // censored time for endpoint 2
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

          if (treatmentGroup[i]==1 && event[i]) events1++;
          if (treatmentGroup[i]==2 && event[i]) events2++;
          if (treatmentGroup[i]==1 && dropoutEvent[i]) dropouts1++;
          if (treatmentGroup[i]==2 && dropoutEvent[i]) dropouts2++;
        }
      }


      // add raw data to output
      if (niter2[k] < maxNumberOfRawDatasetsPerStage) {
        for (i=0; i<n; i++) {
          iterationNumber2x[index2x] = iter+1;
          stageNumber2x[index2x] = k+1;
          subjectId2x[index2x] = i+1;
          arrivalTime2x[index2x] = arrivalTime[i];
          stratum2x[index2x] = stratum[i];
          treatmentGroup2x[index2x] = treatmentGroup[i];
          survivalTime2x[index2x] = survivalTime[i];
          dropoutTime2x[index2x] = dropoutTime[i];
          observationTime2x[index2x] = observationTime[i];
          timeUnderObservation2x[index2x] = timeUnderObservation[i];
          eventx[index2x] = event[i];
          dropoutEventx[index2x] = dropoutEvent[i];

          index2x++;
        }

        // update the number of stage k dataset to extract
        niter2[k]++;
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
      order = seq(0, n-1);
      std::sort(order.begin(), order.end(), [&](int i, int j) {
        return timeUnderObservation[i] < timeUnderObservation[j];
      });

      timeUnderObservationSorted = timeUnderObservation[order];
      eventSorted = event[order];
      stratumSorted = stratum[order];
      treatmentGroupSorted = treatmentGroup[order];

      sub = (timeUnderObservationSorted > 0);
      eventSorted = eventSorted[sub];
      stratumSorted = stratumSorted[sub];
      treatmentGroupSorted = treatmentGroupSorted[sub];
      nsub = static_cast<int>(eventSorted.size());

      // calculate the stratified log-rank test
      uscore = 0;
      vscore = 0;
      for (i=0; i<nsub; i++) {
        h = stratumSorted[i] - 1;
        nt[h] = n1[h] + n2[h];

        n1a[h] = n1[h]*hazardRatioH0;
        nta[h] = n1a[h] + n2[h];

        if (eventSorted[i]) {
          uscore += (treatmentGroupSorted[i]==1) - n1a[h]/nta[h];
          vscore += n1a[h]*n2[h]/(nta[h]*nta[h]);
        }

        // reduce the risk set
        if (treatmentGroupSorted[i]==1) {
          n1[h]--;
        } else {
          n2[h]--;
        }
      }

      // add summary data to output
      iterationNumber2y[index2y] = iter+1;
      eventsNotAchievedy[index2y] = eventsNotAchieved;
      stageNumber2y[index2y] = k+1;
      analysisTime2y[index2y] = analysisTime2[k];
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
      index2y++;

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
    subjectId1x = subjectId1x[sub1];
    arrivalTime1x = arrivalTime1x[sub1];
    stratum1x = stratum1x[sub1];
    treatmentGroup1x = treatmentGroup1x[sub1];
    survivalTime1x = survivalTime1x[sub1];
    dropoutTime1x = dropoutTime1x[sub1];
    ptfu1Timex = ptfu1Timex[sub1];
    observationTime1x = observationTime1x[sub1];
    timeUnderObservation1x = timeUnderObservation1x[sub1];
    responderx = responderx[sub1];
    sourcex = sourcex[sub1];

    DataFrame rawdataBIN = DataFrame::create(
      _["iterationNumber"] = iterationNumber1x,
      _["stageNumber"] = stageNumber1x,
      _["analysisTime"] = observationTime1x,
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
    subjectId2x = subjectId2x[sub2];
    arrivalTime2x = arrivalTime2x[sub2];
    stratum2x = stratum2x[sub2];
    treatmentGroup2x = treatmentGroup2x[sub2];
    survivalTime2x = survivalTime2x[sub2];
    dropoutTime2x = dropoutTime2x[sub2];
    observationTime2x = observationTime2x[sub2];
    timeUnderObservation2x = timeUnderObservation2x[sub2];
    eventx = eventx[sub2];
    dropoutEventx = dropoutEventx[sub2];

    DataFrame rawdataTTE = DataFrame::create(
      _["iterationNumber"] = iterationNumber2x,
      _["stageNumber"] = stageNumber2x,
      _["analysisTime"] = observationTime2x,
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

