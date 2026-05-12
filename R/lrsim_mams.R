#' @title Log-Rank Test Simulation for Multi-Arm Multi-Stage Design
#' @description Simulate a multi-arm multi-stage design using a weighted
#' log-rank test. Analyses can be triggered either by the cumulative number
#' of events (combined for an active arm and the common control) or by
#' pre-specified calendar times.
#'
#' @param M Number of active treatment arms.
#' @param kMax Number of sequential looks.
#' @param criticalValues The matrix of by-level upper boundaries on the
#'   max z-test statistic scale for efficacy stopping.
#'   The first column is for level \eqn{M}, the second column is for
#'   level \eqn{M - 1}, and so on, with the last column for level 1.
#'   Decision rule:
#'   - At Look 1, compute the Wald statistic for each active arm versus the
#'     common control. If the maximum of these statistics exceeds the Look 1
#'     critical value, stop for efficacy and check whether
#'     there is any other active arm which can be rejected using
#'     a relaxed boundary under the closed testing principle.
#'   - If the Look 1 stopping rule is not met, continue to the next look;
#'     if the maximum of Wald statistics at this look exceeds the
#'     corresponding level \eqn{M} critical value, stop for efficacy;
#'     otherwise continue.
#'   - If no critical value is exceeded by Look \code{kMax}, the procedure
#'     ends without rejection.
#' @param futilityBounds Numeric vector of length \code{kMax - 1} giving the
#'   futility boundaries on the max-Z scale for the first \code{kMax - 1}
#'   analyses. At an interim look, the study stops for futility if all active
#'   treatment arms cross the futility boundary. At the final look, the study
#'   is counted as stopping for futility if none of the active treatment arms
#'   can be rejected. If omitted, no interim futility stopping is applied.
#' @param hazardRatioH0s Numeric vector of length \eqn{M}. Hazard ratios
#'   under \eqn{H_0} for each active arm versus the common control. Defaults
#'   to 1 for superiority tests.
#' @param allocations Integer or integer vector of length \eqn{M + 1}.
#'   Number of subjects per arm within a randomization block. A single value
#'   implies equal allocation; defaults to 1.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_stratumFraction
#' @param lambdas List of length \eqn{M} (one element per arm). Each element
#'   is a scalar or a numeric vector of event hazard rates for the
#'   corresponding arm, given by analysis interval and stratum as required
#'   by the simulation.
#' @param gammas List of length \eqn{M} (one element per arm). Each element
#'   is a scalar or a numeric vector of dropout hazard rates for the
#'   corresponding arm, by analysis interval and stratum.
#' @param n Planned total sample size across all active arms and control.
#' @inheritParams param_followupTime
#' @inheritParams param_fixedFollowup
#' @inheritParams param_rho1
#' @inheritParams param_rho2
#' @param plannedEvents Numeric vector of length \eqn{K + 1} giving the
#'   planned cumulative number of events to trigger Look 1 through Look
#'   \eqn{K + 1}. Each entry refers to the combined events for the first
#'   active arm and the common control. Use \code{plannedEvents} to schedule
#'   event-driven looks.
#' @param plannedTime Numeric vector of calendar times for the analyses.
#'   If \code{plannedTime} is supplied, analyses are scheduled by calendar
#'   time and \code{plannedEvents} should be left missing.
#' @param maxNumberOfIterations Number of Monte Carlo replications.
#'   Defaults to 1000.
#' @param maxNumberOfRawDatasetsPerStage Number of subject-level
#'   raw datasets to retain per stage (for selected replications).
#' @param seed Random seed for reproducibility.
#' @param nthreads Number of threads for parallel simulation.
#'   Use 0 to accept the default RcppParallel behavior.
#'
#' @return An S3 object of class \code{"lrsim_seamless"} with these components:
#'
#' * \code{overview}: A list summarizing trial-level results and settings:
#'     - \code{overallReject}: Overall probability of rejecting the null
#'       by trial end.
#'     - \code{overallFutility}: Overall probability of stopping for futility
#'       by trial end.
#'     - \code{rejectPerStage}: Probability of rejecting the null for each
#'       active arm at each stage.
#'     - \code{futilityPerStage}: Probability of futility stopping for each arm
#'       at each stage.
#'     - \code{cumulativeRejection}: Cumulative probability of rejection
#'       for each active arm by stage.
#'     - \code{cumulativeFutility}: Cumulative probability of futility
#'       stopping for each active arm by stage.
#'     - \code{numberOfEvents}: Cumulative event counts by stage and arm.
#'     - \code{numberOfDropouts}: Cumulative dropouts by stage and arm.
#'     - \code{numberOfSubjects}: Cumulative enrollments by stage and arm.
#'     - \code{analysisTime}: Average calendar time for each stage by arm
#'       among replications that reached that stage.
#'     - \code{expectedNumberOfEvents}: Expected cumulative events at trial end.
#'     - \code{expectedNumberOfDropouts}: Expected cumulative dropouts at trial end.
#'     - \code{expectedNumberOfSubjects}: Expected cumulative enrollments
#'       at trial end.
#'     - \code{expectedStudyDuration}: Expected study duration.
#'     - \code{criticalValues}: The input matrix of by-level critical boundaries.
#'     - \code{futilityBounds}: The input vector of futility boundaries.
#'     - \code{hazardRatioH0s}: The input hazard ratios under \eqn{H_0}.
#'     - \code{useEvents}: Logical indicating whether analyses were event-driven.
#'     - \code{numberOfIterations}: Number of simulation iterations performed.
#'     - \code{n}: Planned total sample size.
#'     - \code{fixedFollowup}: Logical indicating whether fixed follow-up was used.
#'     - \code{rho1}, \code{rho2}: Fleming–Harrington weighting parameters used.
#'     - \code{M}: Number of active arms in Phase 2.
#'     - \code{K}: Number of sequential looks in Phase 3.
#'
#' * \code{sumdata1}: Data frame summarizing each iteration, stage, and
#'   treatment group:
#'     - \code{iterationNumber}, \code{eventsNotAchieved},
#'       \code{stopStage}, \code{stageNumber},
#'       \code{analysisTime}, \code{treatmentGroup}, \code{accruals},
#'       \code{events}, \code{dropouts}.
#'     - For each stage the final row summarizes the overall study
#'       (all arms combined).
#'
#' * \code{summdata2}: Data frame summarizing log-rank statistics by iteration,
#'   stage, and active arm:
#'     - \code{iterationNumber}, \code{stopStage}, \code{stageNumber},
#'       \code{analysisTime}, \code{activeArm},
#'       \code{totalAccruals}, \code{totalEvents}, \code{totalDropouts},
#'       \code{uscore}, \code{vscore}, \code{logRankStatistic},
#'       \code{reject}, \code{futility}.
#'     - For each active arm, total accruals, events, and dropouts refer to
#'       the combined counts for that arm and the common control at that stage.
#'
#' * \code{rawdata} (present when \code{maxNumberOfRawDatasetsPerStage} > 0):
#'   Subject-level data for selected replications with variables:
#'     - \code{iterationNumber}, \code{stopStage}, \code{stageNumber},
#'       \code{analysisTime}, \code{subjectId}, \code{arrivalTime},
#'       \code{stratum}, \code{treatmentGroup},
#'       \code{survivalTime}, \code{dropoutTime}, \code{timeUnderObservation},
#'       \code{event}, \code{dropoutEvent}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' (sim1 <- lrsim_mams(
#'   M = 2,
#'   kMax = 3,
#'   criticalValues = matrix(c(3.880, 2.747, 2.275,
#'                             3.710, 2.511, 1.993), 3, 2),
#'   futilityBounds = c(0.074, 1.207),
#'   accrualTime = c(0, 8),
#'   accrualIntensity = c(10, 28),
#'   piecewiseSurvivalTime = 0,
#'   lambdas = list(log(2)/12*0.5, log(2)/12*0.75, log(2)/12),
#'   n = 700,
#'   plannedEvents = c(36, 72, 108),
#'   maxNumberOfIterations = 10000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 0))
#'
#' @export
lrsim_mams <- function(
    M = 2,
    kMax = 1,
    criticalValues = NULL,
    futilityBounds = NULL,
    hazardRatioH0s = 1,
    allocations = 1,
    accrualTime = 0,
    accrualIntensity = NA,
    piecewiseSurvivalTime = 0,
    stratumFraction = 1,
    lambdas = NULL,
    gammas = NULL,
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

  lrsim_mams_Rcpp(
    M, kMax, criticalValues, futilityBounds, hazardRatioH0s,
    allocations, accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction, lambdas, gammas,
    n, followupTime, fixedFollowup, rho1, rho2,
    plannedEvents, plannedTime, maxNumberOfIterations,
    maxNumberOfRawDatasetsPerStage, seed)
}
