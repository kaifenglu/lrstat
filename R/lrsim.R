#' @title Log-Rank Test Simulation
#' @description Performs simulation for two-arm group sequential
#' trials based on weighted log-rank test.
#'
#' @inheritParams param_kMax
#' @param informationRates The information rates in terms of number
#'   of events for the conventional log-rank test and in terms of
#'   the actual information for weighted log-rank tests.
#'   Fixed prior to the trial. If left unspecified, it defaults to
#'   \code{plannedEvents / plannedEvents[kMax]} when \code{plannedEvents}
#'   is provided and to \code{plannedTime / plannedTime[kMax]} otherwise.
#' @inheritParams param_criticalValues
#' @inheritParams param_futilityBounds
#' @inheritParams param_hazardRatioH0
#' @param allocation1 Number of subjects in the active treatment group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param allocation2 Number of subjects in the control group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_stratumFraction
#' @inheritParams param_lambda1_stratified
#' @inheritParams param_lambda2_stratified
#' @inheritParams param_gamma1_stratified
#' @inheritParams param_gamma2_stratified
#' @param n Sample size.
#' @inheritParams param_followupTime
#' @inheritParams param_fixedFollowup
#' @inheritParams param_rho1
#' @inheritParams param_rho2
#' @param plannedEvents The planned cumulative total number of events at
#'   each stage.
#' @param plannedTime The calendar times for the analyses. To use calendar
#'   time to plan the analyses, \code{plannedEvents} should be missing.
#' @param maxNumberOfIterations The number of simulation iterations.
#'   Defaults to 1000.
#' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
#'   stage to extract.
#' @param seed The seed to reproduce the simulation results.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return An S3 class \code{lrsim} object with 3 components:
#'
#' * \code{overview}: A list containing the following information:
#'
#'     - \code{rejectPerStage}: The efficacy stopping probability by stage.
#'
#'     - \code{futilityPerStage}: The futility stopping probability by
#'       stage.
#'
#'     - \code{cumulativeRejection}: Cumulative efficacy stopping
#'       probability by stage.
#'
#'     - \code{cumulativeFutility}: The cumulative futility stopping
#'       probability by stage.
#'
#'     - \code{numberOfEvents}: The average number of events by stage.
#'
#'     - \code{numberOfDropouts}: The average number of dropouts by stage.
#'
#'     - \code{numberOfSubjects}: The average number of subjects by stage.
#'
#'     - \code{analysisTime}: The average analysis time by stage.
#'
#'     - \code{overallReject}: The overall rejection probability.
#'
#'     - \code{expectedNumberOfEvents}: The expected number of events for
#'       the overall study.
#'
#'     - \code{expectedNumberOfDropouts}: The expected number of dropouts
#'       for the overall study.
#'
#'     - \code{expectedNumberOfSubjects}: The expected number of subjects
#'       for the overall study.
#'
#'     - \code{expectedStudyDuration}: The expected study duration.
#'
#'     - \code{hazardRatioH0}: Hazard ratio under the null hypothesis for
#'       the active treatment versus control.
#'
#'     - \code{useEvents}: whether the analyses are planned
#'       based on the number of events or calendar time.
#'
#'     - \code{n}: Sample size.
#'
#'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
#'
#'     - \code{rho1}: The first parameter of the Fleming-Harrington family
#'       of weighted log-rank test. Defaults to 0 for conventional log-rank
#'       test.
#'
#'     - \code{rho2}: The second parameter of the Fleming-Harrington family
#'       of weighted log-rank test. Defaults to 0 for conventional log-rank
#'       test.
#'
#'     - \code{kMax}: The maximum number of stages.
#'
#' * \code{sumdata}: A data frame of summary data by iteration and stage:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{eventsNotAchieved}: Whether the final target number of events
#'       is not achieved for the iteration.
#'
#'     - \code{stopStage}: The stage at which the trial stops.
#'
#'     - \code{stageNumber}: The stage number, covering all stages even if
#'       the trial stops at an interim look.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{accruals1}: The number of subjects enrolled at the stage for
#'       the treatment group.
#'
#'     - \code{accruals2}: The number of subjects enrolled at the stage for
#'       the control group.
#'
#'     - \code{totalAccruals}: The total number of subjects enrolled at
#'       the stage.
#'
#'     - \code{events1}: The number of events at the stage for
#'       the treatment group.
#'
#'     - \code{events2}: The number of events at the stage for
#'       the control group.
#'
#'     - \code{totalEvents}: The total number of events at the stage.
#'
#'     - \code{dropouts1}: The number of dropouts at the stage for
#'       the treatment group.
#'
#'     - \code{dropouts2}: The number of dropouts at the stage for
#'       the control group.
#'
#'     - \code{totalDropouts}: The total number of dropouts at the stage.
#'
#'     - \code{uscore}: The numerator of the log-rank test statistic.
#'
#'     - \code{vscore}: The variance of the log-rank test statistic.
#'
#'     - \code{logRankStatistic}: The log-rank test Z-statistic.
#'
#'     - \code{rejectPerStage}: Whether to reject the null hypothesis
#'       at the stage.
#'
#'     - \code{futilityPerStage}: Whether to stop the trial for futility
#'       at the stage.
#'
#' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
#'   positive integer): A data frame for subject-level data for selected
#'   replications, containing the following variables:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{stopStage}: The stage at which the trial stops.
#'
#'     - \code{stageNumber}: The stage number, covering all stages even if
#'       the trial stops at an interim look.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{subjectId}: The subject ID.
#'
#'     - \code{arrivalTime}: The enrollment time for the subject.
#'
#'     - \code{stratum}: The stratum for the subject.
#'
#'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
#'       subject.
#'
#'     - \code{survivalTime}: The underlying survival time for the subject.
#'
#'     - \code{dropoutTime}: The underlying dropout time for the subject.
#'
#'     - \code{timeUnderObservation}: The time under observation
#'       since randomization.
#'
#'     - \code{event}: Whether the subject experienced the event.
#'
#'     - \code{dropoutEvent}: Whether the subject dropped out.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' # Example 1: analyses based on number of events
#'
#' sim1 = lrsim(
#'   kMax = 2, informationRates = c(0.5, 1),
#'   criticalValues = c(2.797, 1.977),
#'   accrualIntensity = 11,
#'   lambda1 = 0.018, lambda2 = 0.030,
#'   n = 132,
#'   plannedEvents = c(60, 120),
#'   maxNumberOfIterations = 1000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 1)
#'
#' # summary statistics
#' sim1
#'
#' # summary for each simulated data set
#' head(sim1$sumdata)
#'
#' # raw data for selected replication
#' head(sim1$rawdata)
#'
#'
#' # Example 2: analyses based on calendar time have similar power
#'
#' sim2 = lrsim(
#'   kMax = 2, informationRates = c(0.5, 1),
#'   criticalValues = c(2.797, 1.977),
#'   accrualIntensity = 11,
#'   lambda1 = 0.018, lambda2 = 0.030,
#'   n = 132,
#'   plannedTime = c(31.9, 113.2),
#'   maxNumberOfIterations = 1000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 1)
#'
#' # summary statistics
#' sim2
#'
#' # summary for each simulated data set
#' head(sim2$sumdata)
#'
#' @export
lrsim <- function(
    kMax = 1,
    informationRates = NA,
    criticalValues = NA,
    futilityBounds = NA,
    hazardRatioH0 = 1,
    allocation1 = 1,
    allocation2 = 1,
    accrualTime = 0,
    accrualIntensity = NA,
    piecewiseSurvivalTime = 0,
    stratumFraction = 1,
    lambda1 = NA,
    lambda2 = NA,
    gamma1 = 0,
    gamma2 = 0,
    n = NA,
    followupTime = NA,
    fixedFollowup = FALSE,
    rho1 = 0,
    rho2 = 0,
    plannedEvents = NA,
    plannedTime = NA,
    maxNumberOfIterations = 1000,
    maxNumberOfRawDatasetsPerStage = 0,
    seed = 0,
    nthreads = 0) {

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  lrsimRcpp(kMax, informationRates, criticalValues, futilityBounds,
            hazardRatioH0, allocation1, allocation2,
            accrualTime, accrualIntensity,
            piecewiseSurvivalTime, stratumFraction,
            lambda1, lambda2, gamma1, gamma2,
            n, followupTime, fixedFollowup, rho1, rho2,
            plannedEvents, plannedTime, maxNumberOfIterations,
            maxNumberOfRawDatasetsPerStage, seed)
}


#' @title Log-Rank Test Simulation for Three Arms
#' @description Performs simulation for three-arm group sequential trials
#' based on weighted log-rank test. The looks are driven by the total
#' number of events in Arm A and Arm C combined. Alternatively,
#' the analyses can be planned to occur at specified calendar times.
#'
#' @inheritParams param_kMax
#' @param hazardRatioH013 Hazard ratio under the null hypothesis for arm 1
#'   versus arm 3. Defaults to 1 for superiority test.
#' @param hazardRatioH023 Hazard ratio under the null hypothesis for arm 2
#'   versus arm 3. Defaults to 1 for superiority test.
#' @param hazardRatioH012 Hazard ratio under the null hypothesis for arm 1
#'   versus arm 2. Defaults to 1 for superiority test.
#' @param allocation1 Number of subjects in Arm A in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param allocation2 Number of subjects in Arm B in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param allocation3 Number of subjects in Arm C in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_stratumFraction
#' @param lambda1 A vector of hazard rates for the event in each analysis
#'   time interval by stratum for arm 1.
#' @param lambda2 A vector of hazard rates for the event in each analysis
#'   time interval by stratum for arm 2.
#' @param lambda3 A vector of hazard rates for the event in each analysis
#'   time interval by stratum for arm 3.
#' @param gamma1 The hazard rate for exponential dropout. A vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for arm 1.
#' @param gamma2 The hazard rate for exponential dropout. A vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for arm 2.
#' @param gamma3 The hazard rate for exponential dropout. A vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for arm 3.
#' @param n Sample size.
#' @inheritParams param_followupTime
#' @inheritParams param_fixedFollowup
#' @inheritParams param_rho1
#' @inheritParams param_rho2
#' @param plannedEvents The planned cumulative total number of events at
#'   Look 1 to Look \code{kMax} for Arms A and C combined.
#' @param plannedTime The calendar times for the analyses. To use calendar
#'   time to plan the analyses, \code{plannedEvents} should be missing.
#' @param maxNumberOfIterations The number of simulation iterations.
#'   Defaults to 1000.
#' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
#'   stage to extract.
#' @param seed The seed to reproduce the simulation results.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return A list with 2 components:
#'
#' * \code{sumdata}: A data frame of summary data by iteration and stage:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{eventsNotAchieved}: Whether the target number of events
#'       is not achieved for the iteration.
#'
#'     - \code{stageNumber}: The stage number, covering all stages even if
#'       the trial stops at an interim look.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{accruals1}: The number of subjects enrolled at the stage for
#'       the active treatment 1 group.
#'
#'     - \code{accruals2}: The number of subjects enrolled at the stage for
#'       the active treatment 2 group.
#'
#'     - \code{accruals3}: The number of subjects enrolled at the stage for
#'       the control group.
#'
#'     - \code{totalAccruals}: The total number of subjects enrolled at
#'       the stage.
#'
#'     - \code{events1}: The number of events at the stage for
#'       the active treatment 1 group.
#'
#'     - \code{events2}: The number of events at the stage for
#'       the active treatment 2 group.
#'
#'     - \code{events3}: The number of events at the stage for
#'       the control group.
#'
#'     - \code{totalEvents}: The total number of events at the stage.
#'
#'     - \code{dropouts1}: The number of dropouts at the stage for
#'       the active treatment 1 group.
#'
#'     - \code{dropouts2}: The number of dropouts at the stage for
#'       the active treatment 2 group.
#'
#'     - \code{dropouts3}: The number of dropouts at the stage for
#'       the control group.
#'
#'     - \code{totalDropouts}: The total number of dropouts at the stage.
#'
#'     - \code{uscore13}: The log-rank test score statistic comparing the
#'       active treatment 1 to the control.
#'
#'     - \code{vscore13}: The log-rank test variance statistic comparing the
#'       active treatment 1 to the control.
#'
#'     - \code{logRankStatistic13}: The log-rank test Z-statistic
#'       comparing the active treatment 1 to the control.
#'
#'     - \code{uscore23}: The log-rank test score statistic comparing the
#'       active treatment 2 to the control.
#'
#'     - \code{vscore23}: The log-rank test variance statistic comparing the
#'       active treatment 2 to the control.
#'
#'     - \code{logRankStatistic23}: The log-rank test Z-statistic
#'       comparing the active treatment 2 to the control.
#'
#'     - \code{uscore12}: The log-rank test score statistic comparing the
#'       active treatment 1 to the active treatment 2.
#'
#'     - \code{vscore12}: The log-rank test variance statistic comparing the
#'       active treatment 1 to the active treatment 2.
#'
#'     - \code{logRankStatistic12}: The log-rank test Z-statistic
#'       comparing the active treatment 1 to the active treatment 2.
#'
#' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
#'   positive integer): A data frame for subject-level data for selected
#'   replications, containing the following variables:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{stageNumber}: The stage under consideration.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{subjectId}: The subject ID.
#'
#'     - \code{arrivalTime}: The enrollment time for the subject.
#'
#'     - \code{stratum}: The stratum for the subject.
#'
#'     - \code{treatmentGroup}: The treatment group (1, 2, or 3) for
#'       the subject.
#'
#'     - \code{survivalTime}: The underlying survival time for the subject.
#'
#'     - \code{dropoutTime}: The underlying dropout time for the subject.
#'
#'     - \code{timeUnderObservation}: The time under observation
#'       since randomization for the subject.
#'
#'     - \code{event}: Whether the subject experienced the event.
#'
#'     - \code{dropoutEvent}: Whether the subject dropped out.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' sim1 = lrsim3a(
#'   kMax = 3,
#'   allocation1 = 2,
#'   allocation2 = 2,
#'   allocation3 = 1,
#'   accrualTime = c(0, 8),
#'   accrualIntensity = c(10, 28),
#'   piecewiseSurvivalTime = 0,
#'   lambda1 = log(2)/12*0.60,
#'   lambda2 = log(2)/12*0.70,
#'   lambda3 = log(2)/12,
#'   n = 700,
#'   plannedEvents = c(186, 259, 295),
#'   maxNumberOfIterations = 1000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 1)
#'
#' head(sim1$sumdata)
#' head(sim1$rawdata)
#'
#' @export
lrsim3a <- function(
  kMax = 1,
  hazardRatioH013 = 1,
  hazardRatioH023 = 1,
  hazardRatioH012 = 1,
  allocation1 = 1,
  allocation2 = 1,
  allocation3 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  lambda1 = NA,
  lambda2 = NA,
  lambda3 = NA,
  gamma1 = 0,
  gamma2 = 0,
  gamma3 = 0,
  n = NA,
  followupTime = NA,
  fixedFollowup = FALSE,
  rho1 = 0,
  rho2 = 0,
  plannedEvents = NA,
  plannedTime = NA,
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 0,
  seed = 0,
  nthreads = 0) {

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  lrsim3aRcpp(kMax, hazardRatioH013, hazardRatioH023, hazardRatioH012,
              allocation1, allocation2, allocation3,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction,
              lambda1, lambda2, lambda3, gamma1, gamma2, gamma3,
              n, followupTime, fixedFollowup, rho1, rho2,
              plannedEvents, plannedTime, maxNumberOfIterations,
              maxNumberOfRawDatasetsPerStage, seed)
}


#' @title Log-Rank Test Simulation for PFS and OS Endpoints
#' @description Performs simulation for two-endpoint (PFS and OS) two-arm
#' group sequential trials based on weighted log-rank test. The first
#' \code{kMaxpfs} looks are driven by the total number of PFS events in
#' two arms combined, and the subsequent looks are driven by the total
#' number of OS events in two arms combined. Alternatively,
#' the analyses can be planned to occur at specified calendar times.
#'
#' @inheritParams param_kMax
#' @param kMaxpfs Number of stages with timing determined by PFS events.
#'   Ranges from 0 (none) to \code{kMax}.
#' @param hazardRatioH0pfs Hazard ratio under the null hypothesis for the
#'   active treatment vs control for PFS. Defaults to 1 for
#'   superiority test.
#' @param hazardRatioH0os Hazard ratio under the null hypothesis for the
#'   active treatment vs control for OS. Defaults to 1 for
#'   superiority test.
#' @param allocation1 Number of subjects in the treatment group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param allocation2 Number of subjects in the control group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_stratumFraction
#' @param rho_pd_os The correlation coefficient for the standard
#'   bivariate normal random variables used to generate time to disease
#'   progression (PD) and time to death using the inverse CDF method.
#' @param lambda1pfs A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the treatment group and PFS.
#' @param lambda2pfs A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the control group and PFS.
#' @param lambda1os A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the treatment group and OS.
#' @param lambda2os A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the control group and OS.
#' @param gamma1pfs The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the treatment group and PFS.
#' @param gamma2pfs The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the control group and PFS.
#' @param gamma1os The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the treatment group and OS.
#' @param gamma2os The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the control group and OS.
#' @param n Sample size.
#' @inheritParams param_followupTime
#' @inheritParams param_fixedFollowup
#' @inheritParams param_rho1
#' @inheritParams param_rho2
#' @param plannedEvents The planned cumulative total number of PFS events at
#'   Look 1 to Look \code{kMaxpfs} and the planned cumulative total number
#'   of OS events at Look \code{kMaxpfs+1} to Look \code{kMax}.
#' @param plannedTime The calendar times for the analyses. To use calendar
#'   time to plan the analyses, \code{plannedEvents} should be missing.
#' @param maxNumberOfIterations The number of simulation iterations.
#'   Defaults to 1000.
#' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
#'   stage to extract.
#' @param seed The seed to reproduce the simulation results.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return A list with 2 components:
#'
#' * \code{sumdata}: A data frame of summary data by iteration and stage:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{eventsNotAchieved}: Whether the target number of events
#'       is not achieved for the iteration.
#'
#'     - \code{stageNumber}: The stage number, covering all stages even if
#'       the trial stops at an interim look.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{accruals1}: The number of subjects enrolled at the stage for
#'       the treatment group.
#'
#'     - \code{accruals2}: The number of subjects enrolled at the stage for
#'       the control group.
#'
#'     - \code{totalAccruals}: The total number of subjects enrolled at
#'       the stage.
#'
#'     - \code{endpoint}: The endpoint (1 for PFS or 2 for OS) under
#'       consideration.
#'
#'     - \code{events1}: The number of events at the stage for
#'       the treatment group.
#'
#'     - \code{events2}: The number of events at the stage for
#'       the control group.
#'
#'     - \code{totalEvents}: The total number of events at the stage.
#'
#'     - \code{dropouts1}: The number of dropouts at the stage for
#'       the treatment group.
#'
#'     - \code{dropouts2}: The number of dropouts at the stage for
#'       the control group.
#'
#'     - \code{totalDropouts}: The total number of dropouts at the stage.
#'
#'     - \code{uscore}: The numerator of the log-rank test statistic for
#'       the endpoint.
#'
#'     - \code{vscore}: The variance of the log-rank test statistic for
#'       the endpoint.
#'
#'     - \code{logRankStatistic}: The log-rank test Z-statistic for
#'       the endpoint.
#'
#' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
#'   positive integer): A data frame for subject-level data for selected
#'   replications, containing the following variables:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{stageNumber}: The stage under consideration.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{subjectId}: The subject ID.
#'
#'     - \code{arrivalTime}: The enrollment time for the subject.
#'
#'     - \code{stratum}: The stratum for the subject.
#'
#'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
#'       subject.
#'
#'     - \code{endpoint}: The endpoint (1 for PFS or 2 for OS) under
#'       consideration for the row. Each subject will have two rows, one
#'       for each endpoint.
#'
#'     - \code{survivalTime}: The underlying survival time for the
#'       event endpoint for the subject.
#'
#'     - \code{dropoutTime}: The underlying dropout time for the
#'       event endpoint for the subject.
#'
#'     - \code{timeUnderObservation}: The time under observation
#'       since randomization for the event endpoint for the subject.
#'
#'     - \code{event}: Whether the subject experienced the event
#'       endpoint.
#'
#'     - \code{dropoutEvent}: Whether the subject dropped out for the
#'       endpoint.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' sim1 = lrsim2e(
#'   kMax = 3,
#'   kMaxpfs = 2,
#'   allocation1 = 2,
#'   allocation2 = 1,
#'   accrualTime = c(0, 8),
#'   accrualIntensity = c(10, 28),
#'   piecewiseSurvivalTime = 0,
#'   rho_pd_os = 0,
#'   lambda1pfs = log(2)/12*0.60,
#'   lambda2pfs = log(2)/12,
#'   lambda1os = log(2)/30*0.65,
#'   lambda2os = log(2)/30,
#'   n = 420,
#'   plannedEvents = c(186, 259, 183),
#'   maxNumberOfIterations = 1000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 1)
#'
#' head(sim1$sumdata)
#' head(sim1$rawdata)
#'
#' @export
lrsim2e <- function(
  kMax = 1,
  kMaxpfs = 1,
  hazardRatioH0pfs = 1,
  hazardRatioH0os = 1,
  allocation1 = 1,
  allocation2 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  rho_pd_os = 0,
  lambda1pfs = NA,
  lambda2pfs = NA,
  lambda1os = NA,
  lambda2os = NA,
  gamma1pfs = 0,
  gamma2pfs = 0,
  gamma1os = 0,
  gamma2os = 0,
  n = NA,
  followupTime = NA,
  fixedFollowup = FALSE,
  rho1 = 0,
  rho2 = 0,
  plannedEvents = NA,
  plannedTime = NA,
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 0,
  seed = 0,
  nthreads = 0) {

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  lrsim2eRcpp(kMax, kMaxpfs, hazardRatioH0pfs, hazardRatioH0os,
              allocation1, allocation2,
              accrualTime, accrualIntensity,
              piecewiseSurvivalTime, stratumFraction, rho_pd_os,
              lambda1pfs, lambda2pfs, lambda1os, lambda2os,
              gamma1pfs, gamma2pfs, gamma1os, gamma2os,
              n, followupTime, fixedFollowup, rho1, rho2,
              plannedEvents, plannedTime, maxNumberOfIterations,
              maxNumberOfRawDatasetsPerStage, seed)
}


#' @title Log-Rank Test Simulation for Two Endpoints (PFS and OS) and
#' Three Arms
#' @description Performs simulation for two-endpoint (PFS and OS)
#' three-arm group sequential trials based on weighted log-rank test.
#' The first \code{kMaxpfs} looks are driven by the total number of
#' PFS events in Arm A and Arm C combined, and the subsequent looks
#' are driven by the total number of OS events in Arm A and Arm C
#' combined. Alternatively, the analyses can be planned to occur at
#' specified calendar times.
#'
#' @inheritParams param_kMax
#' @param kMaxpfs Number of stages with timing determined by PFS events.
#'   Ranges from 0 (none) to \code{kMax}.
#' @param hazardRatioH013pfs Hazard ratio under the null hypothesis for arm 1
#'   vs arm 3 for PFS. Defaults to 1 for superiority test.
#' @param hazardRatioH023pfs Hazard ratio under the null hypothesis for arm 2
#'   vs arm 3 for PFS. Defaults to 1 for superiority test.
#' @param hazardRatioH012pfs Hazard ratio under the null hypothesis for arm 1
#'   vs arm 2 for PFS. Defaults to 1 for superiority test.
#' @param hazardRatioH013os Hazard ratio under the null hypothesis for arm 1
#'   vs arm 3 for OS. Defaults to 1 for superiority test.
#' @param hazardRatioH023os Hazard ratio under the null hypothesis for arm 2
#'   vs arm 3 for OS. Defaults to 1 for superiority test.
#' @param hazardRatioH012os Hazard ratio under the null hypothesis for arm 1
#'   vs arm 2 for OS. Defaults to 1 for superiority test.
#' @param allocation1 Number of subjects in Arm A in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param allocation2 Number of subjects in Arm B in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param allocation3 Number of subjects in Arm C in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_stratumFraction
#' @param rho_pd_os The correlation coefficient for the standard
#'   bivariate normal random variables used to generate time to
#'   disease progression and time to death using the inverse CDF method.
#' @param lambda1pfs A vector of hazard rates for the event in each analysis
#'   time interval by stratum for arm 1 and PFS.
#' @param lambda2pfs A vector of hazard rates for the event in each analysis
#'   time interval by stratum for arm 2 and PFS.
#' @param lambda3pfs A vector of hazard rates for the event in each analysis
#'   time interval by stratum for arm 3 and PFS.
#' @param lambda1os A vector of hazard rates for the event in each analysis
#'   time interval by stratum for arm 1 and OS.
#' @param lambda2os A vector of hazard rates for the event in each analysis
#'   time interval by stratum for arm 2 and OS.
#' @param lambda3os A vector of hazard rates for the event in each analysis
#'   time interval by stratum for arm 3 and OS.
#' @param gamma1pfs The hazard rate for exponential dropout. A vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for arm 1 and PFS.
#' @param gamma2pfs The hazard rate for exponential dropout. A vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for arm 2 and PFS.
#' @param gamma3pfs The hazard rate for exponential dropout. A vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for arm 3 and PFS.
#' @param gamma1os The hazard rate for exponential dropout. A vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for arm 1 and OS.
#' @param gamma2os The hazard rate for exponential dropout. A vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for arm 2 and OS.
#' @param gamma3os The hazard rate for exponential dropout. A vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for arm 3 and OS.
#' @param n Sample size.
#' @inheritParams param_followupTime
#' @inheritParams param_fixedFollowup
#' @inheritParams param_rho1
#' @inheritParams param_rho2
#' @param plannedEvents The planned cumulative total number of PFS events at
#'   Look 1 to Look \code{kMaxpfs} for Arms A and C combined and the planned
#'   cumulative total number of OS events at Look \code{kMaxpfs+1} to Look
#'   \code{kMax} for Arms A and C combined.
#' @param plannedTime The calendar times for the analyses. To use calendar
#'   time to plan the analyses, \code{plannedEvents} should be missing.
#' @param maxNumberOfIterations The number of simulation iterations.
#'   Defaults to 1000.
#' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
#'   stage to extract.
#' @param seed The seed to reproduce the simulation results.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return A list with 2 components:
#'
#' * \code{sumdata}: A data frame of summary data by iteration and stage:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{eventsNotAchieved}: Whether the target number of events
#'       is not achieved for the iteration.
#'
#'     - \code{stageNumber}: The stage number, covering all stages even if
#'       the trial stops at an interim look.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{accruals1}: The number of subjects enrolled at the stage for
#'       the active treatment 1 group.
#'
#'     - \code{accruals2}: The number of subjects enrolled at the stage for
#'       the active treatment 2 group.
#'
#'     - \code{accruals3}: The number of subjects enrolled at the stage for
#'       the control group.
#'
#'     - \code{totalAccruals}: The total number of subjects enrolled at
#'       the stage.
#'
#'     - \code{endpoint}: The endpoint (1 for PFS or 2 for OS) under
#'       consideration.
#'
#'     - \code{events1}: The number of events at the stage for
#'       the active treatment 1 group.
#'
#'     - \code{events2}: The number of events at the stage for
#'       the active treatment 2 group.
#'
#'     - \code{events3}: The number of events at the stage for
#'       the control group.
#'
#'     - \code{totalEvents}: The total number of events at the stage.
#'
#'     - \code{dropouts1}: The number of dropouts at the stage for
#'       the active treatment 1 group.
#'
#'     - \code{dropouts2}: The number of dropouts at the stage for
#'       the active treatment 2 group.
#'
#'     - \code{dropouts3}: The number of dropouts at the stage for
#'       the control group.
#'
#'     - \code{totalDropouts}: The total number of dropouts at the stage.
#'
#'     - \code{logRankStatistic13}: The log-rank test Z-statistic
#'       comparing the active treatment 1 to the control for the endpoint.
#'
#'     - \code{logRankStatistic23}: The log-rank test Z-statistic
#'       comparing the active treatment 2 to the control for the endpoint.
#'
#'     - \code{logRankStatistic12}: The log-rank test Z-statistic
#'       comparing the active treatment 1 to the active treatment 2
#'       for the endpoint.
#'
#' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
#'   positive integer): A data frame for subject-level data for selected
#'   replications, containing the following variables:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{stageNumber}: The stage under consideration.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{subjectId}: The subject ID.
#'
#'     - \code{arrivalTime}: The enrollment time for the subject.
#'
#'     - \code{stratum}: The stratum for the subject.
#'
#'     - \code{treatmentGroup}: The treatment group (1, 2, or 3) for
#'       the subject.
#'
#'     - \code{endpoint}: The endpoint (1 for PFS or 2 for OS) under
#'       consideration.
#'
#'     - \code{survivalTime}: The underlying survival time for the
#'       event endpoint for the subject.
#'
#'     - \code{dropoutTime}: The underlying dropout time for the
#'       event endpoint for the subject.
#'
#'     - \code{timeUnderObservation}: The time under observation
#'       since randomization for the event endpoint for the subject.
#'
#'     - \code{event}: Whether the subject experienced the event
#'       endpoint.
#'
#'     - \code{dropoutEvent}: Whether the subject dropped out for
#'       the endpoint.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' sim1 = lrsim2e3a(
#'   kMax = 3,
#'   kMaxpfs = 2,
#'   allocation1 = 2,
#'   allocation2 = 2,
#'   allocation3 = 1,
#'   accrualTime = c(0, 8),
#'   accrualIntensity = c(10, 28),
#'   piecewiseSurvivalTime = 0,
#'   rho_pd_os = 0,
#'   lambda1pfs = log(2)/12*0.60,
#'   lambda2pfs = log(2)/12*0.70,
#'   lambda3pfs = log(2)/12,
#'   lambda1os = log(2)/30*0.65,
#'   lambda2os = log(2)/30*0.75,
#'   lambda3os = log(2)/30,
#'   n = 700,
#'   plannedEvents = c(186, 259, 183),
#'   maxNumberOfIterations = 1000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 1)
#'
#' head(sim1$sumdata)
#' head(sim1$rawdata)
#'
#' @export
lrsim2e3a <- function(
  kMax = 1,
  kMaxpfs = 1,
  hazardRatioH013pfs = 1,
  hazardRatioH023pfs = 1,
  hazardRatioH012pfs = 1,
  hazardRatioH013os = 1,
  hazardRatioH023os = 1,
  hazardRatioH012os = 1,
  allocation1 = 1,
  allocation2 = 1,
  allocation3 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  rho_pd_os = 0,
  lambda1pfs = NA,
  lambda2pfs = NA,
  lambda3pfs = NA,
  lambda1os = NA,
  lambda2os = NA,
  lambda3os = NA,
  gamma1pfs = 0,
  gamma2pfs = 0,
  gamma3pfs = 0,
  gamma1os = 0,
  gamma2os = 0,
  gamma3os = 0,
  n = NA,
  followupTime = NA,
  fixedFollowup = FALSE,
  rho1 = 0,
  rho2 = 0,
  plannedEvents = NA,
  plannedTime = NA,
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 0,
  seed = 0,
  nthreads = 0) {

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  lrsim2e3aRcpp(kMax, kMaxpfs,
                hazardRatioH013pfs, hazardRatioH023pfs, hazardRatioH012pfs,
                hazardRatioH013os, hazardRatioH023os, hazardRatioH012os,
                allocation1, allocation2, allocation3,
                accrualTime, accrualIntensity,
                piecewiseSurvivalTime, stratumFraction, rho_pd_os,
                lambda1pfs, lambda2pfs, lambda3pfs,
                lambda1os, lambda2os, lambda3os,
                gamma1pfs, gamma2pfs, gamma3pfs,
                gamma1os, gamma2os, gamma3os,
                n, followupTime, fixedFollowup, rho1, rho2,
                plannedEvents, plannedTime, maxNumberOfIterations,
                maxNumberOfRawDatasetsPerStage, seed)
}



#' @title Log-Rank Test Simulation for Enrichment Design
#' @description Performs simulation for two-arm group
#' sequential trials based on weighted log-rank test
#' for a biomarker enrichment design. The looks are either
#' driven by the total number of events in the ITT population
#' or the biomarker positive sub population.
#' Alternatively, the analyses can be planned to occur at
#' specified calendar times.
#'
#' @inheritParams param_kMax
#' @param kMaxitt Number of stages with timing determined by events
#'   in the ITT population. Ranges from 0 (none) to \code{kMax}.
#' @param hazardRatioH0itt Hazard ratio under the null hypothesis
#'   for the ITT population. Defaults to 1 for superiority test.
#' @param hazardRatioH0pos Hazard ratio under the null hypothesis
#'   for the biomarker positive sub population. Defaults to 1 for
#'   superiority test.
#' @param hazardRatioH0neg Hazard ratio under the null hypothesis
#'   for the biomarker negative sub population. Defaults to 1 for
#'   superiority test.
#' @param allocation1 Number of subjects in the treatment group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param allocation2 Number of subjects in the control group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_stratumFraction
#' @param p_pos The prevalence of the biomarker positive sub population
#'   in each stratum.
#' @param lambda1itt A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the treatment group in the ITT population.
#' @param lambda2itt A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the control group in the ITT population.
#' @param lambda1pos A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the treatment group in the biomarker
#'   positive sub population.
#' @param lambda2pos A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the control group in the biomarker
#'   positive sub population.
#' @param gamma1itt The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the treatment group in the ITT population.
#' @param gamma2itt The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the control group in the ITT population.
#' @param gamma1pos The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the treatment group in the biomarker
#'   positive sub population.
#' @param gamma2pos The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the control group in the biomarker
#'   positive sub population.
#' @param n Sample size.
#' @inheritParams param_followupTime
#' @inheritParams param_fixedFollowup
#' @inheritParams param_rho1
#' @inheritParams param_rho2
#' @param plannedEvents The planned cumulative total number events in the
#'   ITT population at Look 1 to Look \code{kMaxitt} and the planned
#'   cumulative total number of events at Look \code{kMaxitt+1} to
#'   Look \code{kMax} in the biomarker positive sub population.
#' @param plannedTime The calendar times for the analyses. To use calendar
#'   time to plan the analyses, \code{plannedEvents} should be missing.
#' @param maxNumberOfIterations The number of simulation iterations.
#'   Defaults to 1000.
#' @param maxNumberOfRawDatasetsPerStage The number of raw datasets per
#'   stage to extract.
#' @param seed The seed to reproduce the simulation results.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return A list with 2 components:
#'
#' * \code{sumdata}: A data frame of summary data by iteration and stage:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{eventsNotAchieved}: Whether the target number of events
#'       is not achieved for the iteration.
#'
#'     - \code{stageNumber}: The stage number, covering all stages even if
#'       the trial stops at an interim look.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{population}: The population ("ITT", "Biomarker Positive",
#'       "Biomarker Negative") under consideration.
#'
#'     - \code{accruals1}: The number of subjects enrolled at the stage for
#'       the treatment group.
#'
#'     - \code{accruals2}: The number of subjects enrolled at the stage for
#'       the control group.
#'
#'     - \code{totalAccruals}: The total number of subjects enrolled at
#'       the stage.
#'
#'     - \code{events1}: The number of events at the stage for
#'       the treatment group.
#'
#'     - \code{events2}: The number of events at the stage for
#'       the control group.
#'
#'     - \code{totalEvents}: The total number of events at the stage.
#'
#'     - \code{dropouts1}: The number of dropouts at the stage for
#'       the treatment group.
#'
#'     - \code{dropouts2}: The number of dropouts at the stage for
#'       the control group.
#'
#'     - \code{totalDropouts}: The total number of dropouts at the stage.
#'
#'     - \code{logRankStatistic}: The log-rank test Z-statistic for
#'       the population.
#'
#' * \code{rawdata} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
#'   positive integer): A data frame for subject-level data for selected
#'   replications, containing the following variables:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{stageNumber}: The stage under consideration.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{subjectId}: The subject ID.
#'
#'     - \code{arrivalTime}: The enrollment time for the subject.
#'
#'     - \code{stratum}: The stratum for the subject.
#'
#'     - \code{biomarker}: The biomarker status for the subject (1 for
#'       positive, 0 for negative).
#'
#'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
#'       subject.
#'
#'     - \code{survivalTime}: The underlying survival time for the subject.
#'
#'     - \code{dropoutTime}: The underlying dropout time for the subject.
#'
#'     - \code{timeUnderObservation}: The time under observation
#'       since randomization for the subject.
#'
#'     - \code{event}: Whether the subject experienced an event.
#'
#'     - \code{dropoutEvent}: Whether the subject dropped out.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' sim1 = lrsimsub(
#'   kMax = 2,
#'   kMaxitt = 2,
#'   allocation1 = 1,
#'   allocation2 = 1,
#'   accrualTime = seq(0,9),
#'   accrualIntensity = c(seq(10,70,10),rep(70,3)),
#'   piecewiseSurvivalTime = c(0,12,24),
#'   p_pos = 0.6,
#'   lambda1itt = c(0.00256, 0.00383, 0.00700),
#'   lambda2itt = c(0.00427, 0.00638, 0.01167),
#'   lambda1pos = c(0.00299, 0.00430, 0.01064),
#'   lambda2pos = c(0.00516, 0.00741, 0.01835),
#'   gamma1itt = -log(1-0.04)/12,
#'   gamma2itt = -log(1-0.04)/12,
#'   gamma1pos = -log(1-0.04)/12,
#'   gamma2pos = -log(1-0.04)/12,
#'   n = 500,
#'   plannedEvents = c(108,144),
#'   maxNumberOfIterations = 1000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 1)
#'
#' head(sim1$sumdata)
#' head(sim1$rawdata)
#'
#' @export
lrsimsub <- function(
  kMax = 1,
  kMaxitt = 1,
  hazardRatioH0itt = 1,
  hazardRatioH0pos = 1,
  hazardRatioH0neg = 1,
  allocation1 = 1,
  allocation2 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  p_pos = NA,
  lambda1itt = NA,
  lambda2itt = NA,
  lambda1pos = NA,
  lambda2pos = NA,
  gamma1itt = 0,
  gamma2itt = 0,
  gamma1pos = 0,
  gamma2pos = 0,
  n = NA,
  followupTime = NA,
  fixedFollowup = FALSE,
  rho1 = 0,
  rho2 = 0,
  plannedEvents = NA,
  plannedTime = NA,
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 0,
  seed = 0,
  nthreads = 0) {

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  lrsimsubRcpp(kMax, kMaxitt,
               hazardRatioH0itt, hazardRatioH0pos, hazardRatioH0neg,
               allocation1, allocation2,
               accrualTime, accrualIntensity,
               piecewiseSurvivalTime, stratumFraction, p_pos,
               lambda1itt, lambda2itt, lambda1pos, lambda2pos,
               gamma1itt, gamma2itt, gamma1pos, gamma2pos,
               n, followupTime, fixedFollowup, rho1, rho2,
               plannedEvents, plannedTime, maxNumberOfIterations,
               maxNumberOfRawDatasetsPerStage, seed)
}


#' @title Simulation for a Binary and a Time-to-Event Endpoint in
#' Group Sequential Trials
#' @description Performs simulation for two-endpoint two-arm group
#' sequential trials.
#' \itemize{
#'   \item Endpoint 1: Binary endpoint, analyzed using the
#'         Mantel-Haenszel test for risk difference.
#'   \item Endpoint 2: Time-to-event endpoint, analyzed using
#'         the log-rank test for treatment effect.
#' }
#' The analysis times for the binary endpoint are based on calendar times,
#' while the time-to-event analyses are triggered by reaching the
#' pre-specified number of events. The binary endpoint is
#' assessed at the first post-treatment follow-up visit (PTFU1).
#'
#' @param kMax1 Number of stages for the binary endpoint.
#' @param kMax2 Number of stages for the time-to-event endpoint.
#' @param riskDiffH0 Risk difference under the null hypothesis for the
#'   binary endpoint.
#' @param hazardRatioH0 Hazard ratio under the null hypothesis for the
#'   time-to-event endpoint.
#' @param allocation1 Number of subjects in the treatment group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param allocation2 Number of subjects in the control group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_stratumFraction
#' @param globalOddsRatio Global odds ratio of the Plackett copula
#'   linking the two endpoints.
#' @param pi1 Response probabilities by stratum for the treatment group
#'   for the binary endpoint.
#' @param pi2 Response probabilities by stratum for the control group
#'   for the binary endpoint.
#' @param lambda1 A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the treatment group for the time-to-event
#'   endpoint.
#' @param lambda2 A vector of hazard rates for the event in each analysis
#'   time interval by stratum for the control group for the time-to-event
#'   endpoint.
#' @param gamma1 The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the treatment group.
#' @param gamma2 The hazard rate for exponential dropout, a vector of
#'   hazard rates for piecewise exponential dropout applicable for all
#'   strata, or a vector of hazard rates for dropout in each analysis time
#'   interval by stratum for the control group.
#' @param delta1 The hazard rate for exponential treatment discontinuation,
#'   a vector of hazard rates for piecewise exponential treatment
#'   discontinuation applicable for all strata, or a vector of hazard rates
#'   for treatment discontinuation in each analysis time interval by
#'   stratum for the treatment group for the binary endpoint.
#' @param delta2 The hazard rate for exponential treatment discontinuation,
#'   a vector of hazard rates for piecewise exponential treatment
#'   discontinuation applicable for all strata, or a vector of hazard rates
#'   for treatment discontinuation in each analysis time interval by
#'   stratum for the control group for the binary endpoint.
#' @param upper1 Maximim protocol-specified treatment duration for
#'   the treatment group.
#' @param upper2 Maximum protocol-specified treatment duration for
#'   the control group.
#' @param n Sample size.
#' @param plannedTime Calendar times for the analyses of the binary
#'   endpoint.
#' @param plannedEvents Target cumulative number of events for
#'   the time-to-event analyses.
#' @param maxNumberOfIterations Number of simulation iterations to perform.
#' @param maxNumberOfRawDatasetsPerStage Number of subject-level datasets
#'   to retain per stage. Set to 0 to skip raw data saving.
#' @param seed The seed to reproduce the simulation results.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @details We consider dual primary endpoints with endpoint 1 being a
#'   binary endpoint and endpoint 2 being a time-to-event endpoint.
#'   The analyses of endpoint 1 will be based on calendar times, while
#'   the analyses of endpoint 2 will be based on the number of events.
#'   Therefore, the analyses of the two endpoints are not at the same
#'   time points. The correlation between the two endpoints is
#'   characterized by the global odds ratio of the Plackett copula.
#'   In addition, the time-to-event endpoint will render the binary
#'   endpoint as a non-responder, and so does the dropout. In addition,
#'   the treatment discontinuation will impact the number of available
#'   subjects for analysis. The administrative censoring will exclude
#'   subjects from the analysis of the binary endpoint.
#'
#' @return A list with 4 components:
#'
#' * \code{sumdataBIN}: A data frame of summary data by iteration and stage
#'   for the binary endpoint:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{stageNumber}: The stage number, covering all stages even if
#'       the trial stops at an interim look.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{accruals1}: The number of subjects enrolled at the stage for
#'       the treatment group.
#'
#'     - \code{accruals2}: The number of subjects enrolled at the stage for
#'       the control group.
#'
#'     - \code{totalAccruals}: The total number of subjects enrolled at
#'       the stage.
#'
#'     - \code{source1}: The total number of subjects with response status
#'       determined by the underlying latent response variable.
#'
#'     - \code{source2}: The total number of subjects with response status
#'       (non-responder) determined by experiencing the event for the
#'       time-to-event endpoint.
#'
#'     - \code{source3}: The total number of subjects with response status
#'       (non-responder) determined by dropping out prior to the PTFU1
#'       visit.
#'
#'     - \code{n1}: The number of subjects included in the analysis of
#'       the binary endpoint for the treatment group.
#'
#'     - \code{n2}: The number of subjects included in the analysis of
#'       the binary endpoint for the control group.
#'
#'     - \code{n}: The total number of subjects included in the analysis of
#'       the binary endpoint at the stage.
#'
#'     - \code{y1}: The number of responders for the binary endpoint in
#'       the treatment group.
#'
#'     - \code{y2}: The number of responders for the binary endpoint in
#'       the control group.
#'
#'     - \code{y}: The total number of responders for the binary endpoint
#'       at the stage.
#'
#'     - \code{riskDiff}: The estimated risk difference for the binary
#'       endpoint.
#'
#'     - \code{seRiskDiff}: The standard error for risk difference based on
#'       the Sato approximation.
#'
#'     - \code{mhStatistic}: The Mantel-Haenszel test Z-statistic for
#'       the binary endpoint.
#'
#' * \code{sumdataTTE}: A data frame of summary data by iteration and stage
#'   for the time-to-event endpoint:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{eventsNotAchieved}: Whether the target number of events
#'       is not achieved for the iteration.
#'
#'     - \code{stageNumber}: The stage number, covering all stages even if
#'       the trial stops at an interim look.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{accruals1}: The number of subjects enrolled at the stage for
#'       the treatment group.
#'
#'     - \code{accruals2}: The number of subjects enrolled at the stage for
#'       the control group.
#'
#'     - \code{totalAccruals}: The total number of subjects enrolled at
#'       the stage.
#'
#'     - \code{events1}: The number of events at the stage for
#'       the treatment group.
#'
#'     - \code{events2}: The number of events at the stage for
#'       the control group.
#'
#'     - \code{totalEvents}: The total number of events at the stage.
#'
#'     - \code{dropouts1}: The number of dropouts at the stage for
#'       the treatment group.
#'
#'     - \code{dropouts2}: The number of dropouts at the stage for
#'       the control group.
#'
#'     - \code{totalDropouts}: The total number of dropouts at the stage.
#'
#'     - \code{uscore}: The numerator of the log-rank test statistic for the
#'       time-to-event endpoint.
#'
#'     - \code{vscore}: The variance of the log-rank test statistic for the
#'       time-to-event endpoint.
#'
#'     - \code{logRankStatistic}: The log-rank test Z-statistic for
#'       the time-to-event endpoint.
#'
#' * \code{rawdataBIN} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
#'   positive integer): A data frame for subject-level data for the binary
#'   endpoint for selected replications, containing the following variables:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{stageNumber}: The stage under consideration.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{subjectId}: The subject ID.
#'
#'     - \code{arrivalTime}: The enrollment time for the subject.
#'
#'     - \code{stratum}: The stratum for the subject.
#'
#'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
#'       subject.
#'
#'     - \code{survivalTime}: The underlying survival time for the
#'       time-to-event endpoint for the subject.
#'
#'     - \code{dropoutTime}: The underlying dropout time for the
#'       time-to-event endpoint for the subject.
#'
#'     - \code{trtDiscTime}: The underlying treatment discontinuation time
#'       for the binary endpoint for the subject.
#'
#'     - \code{trtDurUpperLimit}: The maximum protocol-specified treatment
#'       duration for the subject based on the treatment group assignment.
#'
#'     - \code{ptfu1Time}:The underlying assessment time for the
#'       binary endpoint for the subject.
#'
#'     - \code{timeUnderObservation}: The time under observation
#'       since randomization for the binary endpoint for the subject.
#'
#'     - \code{latentResponse}: The underlying latent response variable for
#'       the binary endpoint for the subject, which determines the response
#'       status for the binary endpoint at PTFU1 visit.
#'
#'     - \code{responder}: Whether the subject is a responder for the
#'       binary endpoint.
#'
#'     - \code{source}: The source of the determination of responder
#'       status for the binary endpoint: = 1 based on the underlying
#'       latent response variable, = 2 based on the occurrence of
#'       the time-to-event endpoint before the assessment time of the
#'       binary endpoint (imputed as a non-responder), = 3 based on
#'       the dropout before the assessment time of the binary endpoint
#'       (imputed as a non-responder), = 4 excluded from analysis
#'       due to administrative censoring.
#'
#' * \code{rawdataTTE} (exists if \code{maxNumberOfRawDatasetsPerStage} is a
#'   positive integer): A data frame for subject-level data for the
#'   time-to-event endpoint for selected replications, containing the
#'   following variables:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{stageNumber}: The stage under consideration.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{subjectId}: The subject ID.
#'
#'     - \code{arrivalTime}: The enrollment time for the subject.
#'
#'     - \code{stratum}: The stratum for the subject.
#'
#'     - \code{treatmentGroup}: The treatment group (1 or 2) for the
#'       subject.
#'
#'     - \code{survivalTime}: The underlying survival time for the
#'       time-to-event endpoint for the subject.
#'
#'     - \code{dropoutTime}: The underlying dropout time for the
#'       time-to-event endpoint for the subject.
#'
#'     - \code{timeUnderObservation}: The time under observation
#'       since randomization for the time-to-event endpoint for the subject.
#'
#'     - \code{event}: Whether the subject experienced the event for the
#'       time-to-event endpoint.
#'
#'     - \code{dropoutEvent}: Whether the subject dropped out for the
#'       time-to-event endpoint.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' tcut = c(0, 12, 36, 48)
#' surv = c(1, 0.95, 0.82, 0.74)
#' lambda2 = (log(surv[1:3]) - log(surv[2:4]))/(tcut[2:4] - tcut[1:3])
#'
#' sim1 = binary_tte_sim(
#'   kMax1 = 1,
#'   kMax2 = 2,
#'   accrualTime = seq(0, 8),
#'   accrualIntensity = 40/9 * seq(1, 9),
#'   piecewiseSurvivalTime = c(0, 12, 36),
#'   globalOddsRatio = 1,
#'   pi1 = 0.80,
#'   pi2 = 0.65,
#'   lambda1 = 0.65*lambda2,
#'   lambda2 = lambda2,
#'   gamma1 = -log(1-0.04)/12,
#'   gamma2 = -log(1-0.04)/12,
#'   delta1 = -log(1-0.02)/12,
#'   delta2 = -log(1-0.02)/12,
#'   upper1 = 15*28/30.4,
#'   upper2 = 12*28/30.4,
#'   n = 640,
#'   plannedTime = 20 + 15*28/30.4,
#'   plannedEvents = c(130, 173),
#'   maxNumberOfIterations = 1000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 1)
#'
#' @export
binary_tte_sim <- function(
  kMax1 = 1,
  kMax2 = 1,
  riskDiffH0 = 0,
  hazardRatioH0 = 1,
  allocation1 = 1,
  allocation2 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  globalOddsRatio = NA,
  pi1 = NA,
  pi2 = NA,
  lambda1 = NA,
  lambda2 = NA,
  gamma1 = NA,
  gamma2 = NA,
  delta1 = NA,
  delta2 = NA,
  upper1 = NA,
  upper2 = NA,
  n = NA,
  plannedTime = NA,
  plannedEvents = NA,
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 0,
  seed = 0,
  nthreads = 0) {

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  binary_tte_simRcpp(kMax1, kMax2, riskDiffH0, hazardRatioH0,
                     allocation1, allocation2,
                     accrualTime, accrualIntensity,
                     piecewiseSurvivalTime, stratumFraction,
                     globalOddsRatio, pi1, pi2,
                     lambda1, lambda2, gamma1, gamma2,
                     delta1, delta2, upper1, upper2,
                     n, plannedTime, plannedEvents, maxNumberOfIterations,
                     maxNumberOfRawDatasetsPerStage, seed)
}
